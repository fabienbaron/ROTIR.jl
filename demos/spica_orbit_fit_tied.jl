# Spica (alpha Vir) — orbit fit through ROTIR's generic path
#
#   julia --project=demos demos/spica_orbit_fit_tied.jl
#
#   APSIDAL=1 julia --project=demos demos/spica_orbit_fit_tied.jl   # fit apsidal motion too
#   METHOD=neldermead ...                                           # quick smoke test
#   NLIVE=800 ...                                                   # tighter posterior
#
# ---------------------------------------------------------------------------------------
# Scope
# ---------------------------------------------------------------------------------------
# Two limb-darkened disks on an ECCENTRIC orbit. Both stars are resolved, so unlike beta Lyr
# there is no accretion disc and no obvious tie to impose — the demonstration here is that
# the same `fit_orbit` entry point covers a quite different system, with `e` and `omega`
# free rather than pinned at zero.
#
# DELIBERATELY NOT MODELLED: mutual irradiation. Each star heats the inner face of the
# other, which on this system is a ~3 % temperature effect on the secondary's day side
# (see docs/src/guides/reflection.md and demos/spica_irradiation_detect.jl). A uniform or
# limb-darkened DISK cannot represent a one-sided brightness gradient at all, so that signal
# is absorbed into whatever the fitted diameters and flux ratio can accommodate. This script
# fits the ORBIT; for the irradiation question use the tessellated forward model.
#
# ---------------------------------------------------------------------------------------
# Apsidal motion
# ---------------------------------------------------------------------------------------
# Spica's periapsis advances at ~2.6 deg/yr (U ~ 105-139 yr), and the CHARA data span
# 2007-2015 — so omega drifts by roughly 21 deg ACROSS THE DATASET. Ignoring it does not
# merely inflate chi2; it biases omega and T0, because the fit is forced to find a single
# angle that compromises between the early and late epochs.
#
# `domega` is off by default in `orbit_fit_spec` (it is unmeasurable over a short baseline
# and would be a free parameter doing nothing). Here it is worth freeing: set APSIDAL=1.
# The starting value is the published U = 139 yr.

using ROTIR, Printf, Statistics, DelimitedFiles

const HERE    = @__DIR__
const OIFITS  = joinpath(HERE, "data", "2007_2012_2015.Spica.oifits")
const METHOD  = Symbol(get(ENV, "METHOD", "ultranest"))
const NLIVE   = parse(Int, get(ENV, "NLIVE", "400"))
const APSIDAL = get(ENV, "APSIDAL", "0") == "1"
const OUT     = joinpath(HERE, "results")
isfile(OIFITS) || error("Spica OIFITS not found at $OIFITS")

# Split into nightly epochs on 0.5 d MJD gaps — the same rule as
# demos/spica_binary_roche.jl and demos/spica_irradiation_detect.jl. Nightly matters here:
# the period is 4.01 d, so a night is ~6 % of an orbit and lumping a whole campaign into one
# block would smear the separation that carries the astrometry.
raw   = readoifits(OIFITS; T = Float64)[1, 1]
mjds  = sort(raw.v2_mjd)
jumps = findall(diff(mjds) .> 0.5)
starts, stops = mjds[[1; jumps .+ 1]], mjds[[jumps; length(mjds)]]
@printf("Spica: %d V², %d T3amp, %d T3φ over %.1f yr\n",
        raw.nv2, raw.nt3amp, raw.nt3phi, (maximum(mjds) - minimum(mjds))/365.25)
data = [filter_data(raw, set_data_filter(raw; mjd_range = [starts[i] - 0.01, stops[i] + 0.01]))
        for i in eachindex(starts)]
data = [d for d in data if d.nv2 > 0]
@printf("       %d nightly epochs, %d points\n\n",
        length(data), sum(d.nv2 + d.nt3amp + d.nt3phi for d in data))

# Aufdenberg et al. 2015 / demos/spica_params.jl
const R_SUN_M = 6.957e8
const PC_M    = 3.0856775814913673e16
const RAD_MAS = 180 * 3600 * 1000 / pi
D_PC = 77.0
# Angular DIAMETER in mas. `OrbitComponent` diameters are diameters; spica_params.jl works
# in polar RADII, so the factor 2 is not optional.
ang(r_rsun) = 2 * r_rsun * R_SUN_M / (D_PC * PC_M) * RAD_MAS
phys(θ_diam_mas) = θ_diam_mas / 2 / RAD_MAS * (D_PC * PC_M) / R_SUN_M
elements = (a = 1.54, i = 116.0, Omega = 309.938, omega = 255.0, e = 0.123,
            P = 4.0145, T0 = 2454189.40, dP = 0.0,
            domega = APSIDAL ? 360.0/(139.0*365.25) : 0.0)

# Hestroffer power law in H, matching the tessellated demos (LD is band-dependent).
star1 = LimbDarkenedDisk(diameter = ang(7.40), law = :power, ld1 = 0.15)
star2 = LimbDarkenedDisk(diameter = ang(3.74), law = :power, ld1 = 0.15)

free = Symbol[:a, :i, :Omega, :omega, :e, :T0, :f, :c1_diameter, :c2_diameter]
APSIDAL && push!(free, :domega)

@printf("free: %s\n", join(string.(free), ", "))
@printf("apsidal motion: %s\n\n", APSIDAL ? "FITTED" : "fixed at 0 (see header)")

res = fit_orbit(data, star1, star2; elements = elements, flux_ratio = 0.30,
                free = free, method = METHOD, min_num_live_points = NLIVE,
                verbose = true,
                # Omega over the FULL circle, overriding the spec default of [0,180). That
                # default exists because (Omega, T0) -> (Omega+180, T0+P/2) is an exact
                # degeneracy AT e = 0, which is the beta Lyr case. Spica has e = 0.123, and
                # closure phases besides, so the two are distinguishable and folding the
                # published 309.938 into 129.938 would fit a mirrored orbit.
                # P is pinned by spectroscopy.
                bounds = Dict(:a => (1.2, 1.9), :i => (100.0, 135.0),
                              :Omega => (0.0, 360.0),
                              :e => (0.0, 0.30), :omega => (200.0, 310.0),
                              :c1_diameter => (0.5, 1.4), :c2_diameter => (0.2, 0.9),
                              :f => (0.05, 0.9)))

@printf("\n%-16s %14s\n", "parameter", "fitted")
for (n, v) in zip(res.names, res.params)
    @printf("%-16s %14.6f%s\n", n, v, n in free ? "" : "  (fixed)")
end
@printf("\nchi2/n = %.4f\n", res.chi2_red)
@printf("R1 = %.3f Rsun, R2 = %.3f Rsun  (published 7.40, 3.74)\n",
        phys(res.params[findfirst(==(:c1_diameter), res.names)]),
        phys(res.params[findfirst(==(:c2_diameter), res.names)]))

mkpath(OUT)
tag = APSIDAL ? "apsidal" : "fixedw"
writedlm(joinpath(OUT, "spica_orbitfit_$(tag)_best.txt"), hcat(string.(res.names), res.params))
haskey(res, :posterior) &&
    writedlm(joinpath(OUT, "spica_orbitfit_$(tag)_posterior.txt"), res.posterior)
@info "wrote results to $OUT"
