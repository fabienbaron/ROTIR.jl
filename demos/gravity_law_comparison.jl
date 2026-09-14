#!/usr/bin/env julia
# =======================================================================================
# Which gravity-darkening law does the data prefer?
# =======================================================================================
# A rapid rotator's surface brightness is set by two things that are hard to tell apart:
# GRAVITY DARKENING, which cools the equator because the effective gravity is lower there, and
# LIMB DARKENING, which dims the edge of the disc because the line of sight leaves the
# photosphere at a shallower angle. Both take flux out of the limb, and an interferometer
# measures their sum.
#
# That matters because the gravity-darkening law is not settled. von Zeipel (1924) follows from
# assuming the star is barotropic, which is incompatible with radiative equilibrium in a
# rotating star, and for a fast rotator it predicts too much pole-to-equator contrast.
# Espinosa Lara & Rieutord (2011, A&A 533, A43) replace barotropy with "the radiative flux is
# anti-parallel to the local effective gravity" — accurate to better than half a degree even
# near break-up — and get a warmer equator at the same exponent:
#
#     frac_escapevel   T_eq/T_pole vZ   T_eq/T_pole ELR
#          0.50             0.958            0.961
#          0.90             0.788            0.831
#          0.95             0.719            0.789
#
# A FIT FORCED TO USE THE WRONG LAW HAS TO BUY THE DIFFERENCE SOMEWHERE, and limb darkening is
# what is for sale. Asking von Zeipel for a fast rotator's equatorial flux drives the limb
# darkening down, and a coefficient fitted that way can come out NEGATIVE — which is a
# statement about the law it was fitted under, not about the star.
#
# Both laws carry the same free exponent β here, so they have the SAME parameter count and the
# difference in log-evidence is a Bayes factor between them directly.
#
#   julia --project=demos demos/gravity_law_comparison.jl
#   METHOD=nautilus NEFF=2000 NSIDE=4 julia --project=demos demos/gravity_law_comparison.jl
#
# TWO MODES, because the two questions need different samplers. The default runs NUTS, which
# needs only Zygote and takes minutes: it gives the posterior and, with it, the β /
# limb-darkening CORRELATION — the number that says whether the two are trading against each
# other. It cannot answer "which law", because a Hamiltonian sampler computes no evidence.
#
# `METHOD=nautilus` (or `ultranest`) runs a NESTED sampler, which does compute one, and then
# the two runs differ by a Bayes factor. Neither sampler is a dependency of the demos
# environment — add `Nautilus` to `demos/Project.toml` for the pure-Julia one.

# `AdvancedHMC` and `LogDensityProblems` are what activate `ROTIRHMCExt`, which is where
# `_fit_hmc` lives; Zygote supplies the gradient it integrates.
using ROTIR, Zygote, AdvancedHMC, LogDensityProblems, Printf, Statistics, LinearAlgebra

const NSIDE  = parse(Int, get(ENV, "NSIDE", "3"))
const NEFF   = parse(Int, get(ENV, "NEFF", "1000"))
const METHOD = Symbol(get(ENV, "METHOD", "hmc"))

if METHOD !== :hmc
    # Loaded here rather than at the top so the default mode needs nothing beyond Zygote.
    try
        METHOD === :nautilus ? (@eval using Nautilus) : (@eval using PythonCall)
    catch
        error("""
            METHOD=$(METHOD) needs $(METHOD === :nautilus ? "Nautilus.jl" : "PythonCall + UltraNest"),
            which is not in the demos environment. Add it there, or run the default
            METHOD=hmc, which needs only Zygote and still shows the beta / limb-darkening
            correlation (it cannot compute an evidence, so it cannot rank the two laws).""")
    end
end

# lambda And, the multi-epoch set the rest of the demos use. It is a slow rotator, so this
# script is a WORKED EXAMPLE of the comparison rather than a claim about that star: expect the
# two laws to be indistinguishable here, which is itself the right answer at this rotation.
DATA = joinpath(@__DIR__, "data")
files = sort([joinpath(DATA, f) for f in readdir(DATA) if occursin("lam_And", f)])
data_all = readoifits_multiepochs(files; T = Float64)
data = data_all[1, :]
tepochs = Float64[d.mean_mjd for d in data]; tepochs .-= tepochs[1]
tess = tessellation_healpix(NSIDE; T = Float64)
@printf("%d epochs, HEALPix %d (%d tessels), method %s\n",
        length(data), NSIDE, tess.npix, METHOD)

# THE SAME FREE SET UNDER BOTH LAWS, and `beta` is in it. Espinosa Lara & Rieutord derive their
# law with the exponent pinned at 1/4, which assumes a grey radiative atmosphere; leaving it
# pinned would give their law one parameter fewer than von Zeipel and make the evidence ratio a
# comparison of two different things. To recover their published law exactly, drop "beta" from
# this list and set θ0[5] = 0.25.
free   = ["rpole", "frac_escapevel", "inclination", "position_angle", "beta", "ld1"]
θnames = parametric_param_names(; tpole_free = true)
θ0     = [1.37, 0.50, 78.0, 24.0, 0.25, 0.23, 0.0, 4800.0]
idx    = parametric_free_indices(free; tpole_free = true)
lo, hi = default_parametric_bounds(; tpole_free = true)
lo = copy(lo); hi = copy(hi)
lo[1], hi[1] = 0.5, 3.0            # a nested sampler needs finite, and informative, bounds
lo[2], hi[2] = 0.0, 0.995
lo[5], hi[5] = 0.0, 0.5
lo[6], hi[6] = -0.5, 1.0           # WIDE ENOUGH TO GO NEGATIVE, which is the point

results = Dict{Symbol,Any}()
for law in (:vonzeipel, :elr)
    spec = gravity_law_spec(law)
    base = default_star_params(2; ldtype = 3, rotation_period = 54.8, tpole = 4800.0,
                               gravity_law = law)
    @printf("\n=== %s — %s\n", spec.label, spec.reference)

    if METHOD === :hmc
        r = ROTIR._fit_hmc(data, tess, tepochs, base; θ0 = θ0, free = free,
                           lb = lo, ub = hi, tpole_free = true,
                           n_samples = 400, n_adapt = 300, verb = true)
        results[law] = r
        # `_fit_hmc` returns the free entries in `parametric_free_indices` order, which is
        # sorted by position in θ rather than by the order `free` was written in.
        nm = θnames[idx]
        for (j, n) in enumerate(nm)
            @printf("  %-16s %10.4f  (%.4f .. %.4f)\n", n, r.median[j], r.q16[j], r.q84[j])
        end
        # The correlation the negative-limb-darkening question turns on.
        jb = findfirst(==("beta"), nm); jl = findfirst(==("ld1"), nm)
        if jb !== nothing && jl !== nothing
            C = cor(r.samples[:, jb], r.samples[:, jl])
            @printf("  corr(beta, ld1) = %+.3f\n", C)
        end
    else
        # Nested sampling wants a χ², not a log-posterior, and only the free coordinates.
        logπ = build_parametric_logπ(data, tess, tepochs, base; tpole_free = true)
        chi2_of = function (z)
            θ = copy(θ0); θ[idx] .= z
            return -2 * logπ(θ)
        end
        r = ROTIR._fit_nested(METHOD, chi2_of, θnames[idx], lo[idx], hi[idx];
                              verb = true, n_eff = NEFF, n_live = max(NEFF ÷ 4, 200),
                              min_num_live_points = max(NEFF ÷ 4, 200))
        results[law] = r
        for (j, n) in enumerate(θnames[idx])
            @printf("  %-16s %10.4f  (%.4f .. %.4f)\n", n, r.median[j], r.q16[j], r.q84[j])
        end
        @printf("  log Z = %.3f +/- %.3f\n", r.logz, r.logzerr)
    end
end

# The comparison. On Jeffreys' scale a log-evidence difference above about 5 is decisive and
# below 1 is "the data cannot tell" — in which case the slow-rotation law is the honest choice,
# since nothing in the data is paying for the extra physics.
if METHOD !== :hmc
    dlz = results[:elr].logz - results[:vonzeipel].logz
    @printf("\nlog Z(ELR) - log Z(von Zeipel) = %+.3f\n", dlz)
    println(abs(dlz) < 1 ? "  the data do not distinguish the two laws" :
            dlz > 0      ? "  the data prefer Espinosa Lara & Rieutord" :
                           "  the data prefer von Zeipel")
end
jl = findfirst(==("ld1"), θnames[idx])
if jl !== nothing
    println("\nlimb darkening under each law:")
    for law in (:vonzeipel, :elr)
        @printf("  ld1, %-22s = %+.4f\n", string(law), results[law].median[jl])
    end
    println("  A coefficient that is negative under one law and positive under the other is " *
            "the\n  gravity-darkening law showing up as limb darkening, not the star.")
end
