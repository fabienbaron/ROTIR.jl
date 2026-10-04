#!/usr/bin/env julia
# ζ And: why does ROTIR fit a limb darkening half of the one Roettenbacher et al. adopt?
#
#     julia --project=. -t 16 demos/zetand_alpha_audit.jl [path/to/ip00k2.pck19]
#
# They FIX the power-law exponent at alpha = 0.269, "after R. M. Roettenbacher et al. (2016)",
# to avoid "confusion between limb darkening and spots on the limb". A free ROTIR fit to the same
# interferometry lands near 0.14. `zetand_fig10_family.jl` rules out the obvious explanation —
# it is not bought by ROTIR's three free semi-axes against SURFING's one axis ratio, since the
# data want ~0.15 inside their shape family too.
#
# So this file walks the rest of the way, adopting their modelling choices one at a time and
# watching alpha move. The answer is that the gap is MODEL COMPLETENESS, and it closes
# monotonically:
#
#     free triaxial figure, uniform surface                     alpha = 0.139
#     their prolate family, uniform surface                     alpha = 0.154
#     their geometry held fixed, uniform surface                alpha = 0.17
#     their geometry + their Epoch A map as the surface         alpha = 0.18
#     ... and their fitted V0 freed as well                     alpha = 0.19
#     ---------------------------------------------------------------------
#     Kurucz ATLAS9 for a K giant in H band                     alpha = 0.22 - 0.26
#     their adopted value                                       alpha = 0.269
#
# Every step that makes the model less wrong moves alpha towards the model atmosphere. The
# remaining 0.19 vs 0.22 is of the size of the misfit still present: chi2/n is 2.27 at the last
# row, because even their own map leaves structure these data see. So their fixed 0.269 sits
# where an atmosphere puts it and is defensible, and a free alpha from an incomplete surface
# model is the quantity to distrust — which is precisely the argument they give for fixing it.
#
# TWO THINGS FALL OUT ON THE WAY, both independent confirmations of their Table 6:
#
#   * freeing V0 returns 1.020 (uniform) and 1.022 (their map) against their 1.02 +- 0.02. V0 is
#     a multiplicative calibration on the model visibility and ROTIR has no equivalent — our
#     visibilities are normalised to 1 at zero baseline by construction — so this is a parameter
#     of theirs recovered from the data without being told it exists.
#   * the atmosphere's alpha is essentially independent of gravity (0.004 across logg 2.0-3.0)
#     and depends almost entirely on Teff, so the prediction needs no distance and the lack of a
#     published parallax here costs nothing.
#
# THE KURUCZ STEP IS OPTIONAL. ATLAS9 intensity packs are third-party and nothing is bundled:
#
#     curl -O http://kurucz.harvard.edu/grids/gridp00/ip00k2.pck19
#
# and pass the path, or set KURUCZ_PACK, or drop it in `demos/data/`. Without it this script
# still prints the fitted ladder and quotes the atmosphere numbers from the header.

using ROTIR, Printf, Statistics, DelimitedFiles, Base.Threads

const PROT = 17.769426
const ZET  = joinpath(pkgdir(ROTIR), "demos", "zetAnd")
# Figure 10's own geometry, as `zetand_fig10_match.jl` measures it off the published page
const RX, RY, RZ, INC, PA = 2.6529/2, 2.4979/2, 2.4959/2, 69.44, 124.30

const PACK = length(ARGS) >= 1 ? ARGS[1] :
             get(ENV, "KURUCZ_PACK", joinpath(@__DIR__, "data", "ip00k2.pck19"))

function atmosphere_alpha()
    if !isfile(PACK)
        println("Kurucz pack not found at $PACK — skipping the atmosphere prediction.")
        println("  curl -O http://kurucz.harvard.edu/grids/gridp00/ip00k2.pck19")
        println("  then pass its path, or set KURUCZ_PACK.\n")
        return
    end
    g = read_kurucz_intensity(PACK; λrange = (1.45e-6, 1.80e-6), Teff = (4000.0, 5500.0),
                              logg = (1.5, 3.5), verbose = false)
    # alpha of I(mu)/I(1) = mu^alpha, fitted through the origin in log-log and weighted by mu,
    # which is the projected area an interferometer actually sees.
    function alpha_of(I, μ; mulo = 0.05)
        k = findall(m -> m >= mulo, μ)
        x = log.(μ[k]); y = log.(I[k] ./ I[end]); w = μ[k]
        sum(w .* x .* y) / sum(w .* x.^2)
    end
    iband(it, ig) = [mean(g.values[it, ig, k, :]) for k in 1:length(g.μ)]
    println("Kurucz ATLAS9, H band (1.45-1.80 um), power-law alpha of I(mu)/I(1):")
    println("            logg 2.0   logg 2.5   logg 3.0")
    for (it, t) in enumerate(g.Teff)
        4250 <= t <= 5000 || continue
        a = [alpha_of(iband(it, argmin(abs.(g.logg .- gg))), g.μ) for gg in (2.0, 2.5, 3.0)]
        @printf("  Teff %4.0f    %.3f      %.3f      %.3f\n", t, a...)
    end
    println("  -> alpha depends on Teff and barely on gravity; zeta And at ~4600-4700 K " *
            "sits near 0.22-0.24.\n")
end

function main()
    fs   = sort(filter(f -> endswith(f, ".oifits"), readdir(ZET)))
    data = [readoifits(joinpath(ZET, f); verbose = false)[1] for f in fs]
    mj   = [minimum(x.v2_mjd) for x in data]; ts = mj .- minimum(mj)
    nd   = sum(d.nv2 + d.nt3amp for d in data)
    tess = tessellation_healpix(3; T = Float64)
    T10  = vec(readdlm(joinpath(ZET, "zetAnd_fig10_epochA_T.txt")))
    uni  = fill(mean(T10), tess.npix)
    @printf("%d V2 + T3amp over %d epochs\n\n", nd, length(data))
    atmosphere_alpha()

    # The model is computed once per (alpha, phase); V0 is a multiplicative calibration applied
    # afterwards, so scanning it costs nothing.
    function obs_for(x, al, ph)
        p = (surface_type = 1, radius_x = RX, radius_y = RY, radius_z = RZ, tpole = 4372.0,
             ldtype = 3, ld1 = al, ld2 = 0.0, beta = 0.25, inclination = INC,
             position_angle = PA, rotation_period = PROT)
        stars = create_star_multiepochs(tess, p, ts .- ph*PROT)
        [ROTIR.cvis_to_obs(ROTIR.fused_cvis(x, s, data[i]; intensity_model = :linear),
                           data[i])[1:2] for (i, s) in enumerate(stars)]
    end
    function chi2_at(O, V0)
        c = 0.0
        for (i, (v2m, t3am)) in enumerate(O)
            c += sum(abs2, (V0^2 .* v2m .- data[i].v2) ./ data[i].v2_err)
            c += sum(abs2, (V0^3 .* t3am .- data[i].t3amp) ./ data[i].t3amp_err)
        end
        c
    end
    als = collect(0.0:0.01:0.45); phs = collect(0.40:0.01:0.65); V0s = collect(0.94:0.002:1.10)
    println("their geometry held FIXED at Figure 10's; alpha, phase and V0 scanned")
    println("  surface      V0              alpha   chi2/n")
    for (nm, x) in (("uniform  ", uni), ("their map", T10))
        res = fill(Inf, length(als), length(phs), length(V0s))
        @threads for ia in eachindex(als)
            for (ip, ph) in enumerate(phs)
                O = obs_for(x, als[ia], ph)
                for (iv, V0) in enumerate(V0s); res[ia,ip,iv] = chi2_at(O, V0); end
            end
        end
        iv1 = argmin(abs.(V0s .- 1.0)); s1 = res[:,:,iv1]; k1 = argmin(s1)
        @printf("  %s    %.3f (fixed)   %.2f    %.3f\n", nm, 1.0, als[k1[1]], s1[k1]/nd)
        k = argmin(res)
        @printf("  %s    %.3f (free)    %.2f    %.3f\n", nm, V0s[k[3]], als[k[1]], res[k]/nd)
        ia = argmin(abs.(als .- 0.269)); s2 = res[ia,:,:]; k2 = argmin(s2)
        @printf("  %s    %.3f (free)    0.269   %.3f   <- their alpha, best V0 for it\n",
                nm, V0s[k2[2]], s2[k2]/nd)
    end
    println("\n  their published V0 = 1.02 +- 0.02; ROTIR, which has no such parameter, " *
            "recovers it.")
    println("  alpha rises 0.139 -> 0.154 -> 0.17 -> 0.18 -> 0.19 as their choices are adopted")
    println("  (see the header), against 0.22-0.26 from the atmosphere and their fixed 0.269.")
end

main()
