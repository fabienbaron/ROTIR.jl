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
# METHOD defaults to :ultranest, which needs PythonCall to load ROTIRUltraNestExt.
using PythonCall

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

# H-band flux ratio from the published radii and effective temperatures, through a Planck
# ratio at 1.65 um: (R2/R1)^2 * B(T2)/B(T1) = 0.196. The previous starting guess of 0.30 was
# 50 % high — which matters, because it is also where the sampler starts looking.
const F_H = (3.74/7.40)^2 * (exp(14387.77/1.65/25300) - 1) / (exp(14387.77/1.65/20585) - 1)

# LIMB DARKENING: fitted for the PRIMARY ONLY, and the asymmetry is measured, not assumed.
#
# Away from eclipse the only thing carrying information about a star's surface-brightness
# profile is the SHAPE of its visibility curve, which is what the limb-darkening exponent
# controls. Pinning it at a tabulated value does not remove a parameter, it hides a
# degeneracy: a more limb-darkened star of larger diameter gives nearly the same V² as a
# less darkened smaller one, so the tabulated value's error passes into the radius while the
# quoted uncertainty knows nothing about it. Same LD-radius degeneracy the radial
# regularizers in docs/src/api/radial_regularizers.md exist to break for single-star imaging.
#
# BUT THE TWO STARS ARE NOT COMPARABLY RESOLVED. At the longest CHARA baseline here
# (338.5 m, lambda = 1.609 um) the resolution is 0.490 mas, and:
#
#   primary    0.902 mas = 1.84 resolution elements, V² falls to 0.069  -> shape measurable
#   secondary  0.452 mas = 0.92 resolution elements, V² only to 0.578   -> barely resolved
#
# Neither reaches its first null (1.33x and 2.65x beyond B_max). Changing the exponent from
# 0.15 to 0.35 moves V² by 0.015 for the primary against a median sigma of 0.078 — about
# 0.19 sigma per point — but the secondary's equivalent signal is diluted by its 14.5 % flux
# share to 0.023 sigma per point. Over ~2500 points with 8 correlated spectral channels
# (effective N ~ 300) that is roughly 3 sigma for the primary and 0.5 sigma for the
# secondary. So `c2_ld1` stays fixed: freeing it adds prior volume and a strong correlation
# with `c2_diameter` in exchange for no information.
free = Symbol[:a, :i, :Omega, :omega, :e, :T0, :f,
              :c1_diameter, :c1_ld1, :c2_diameter]
APSIDAL && push!(free, :domega)

@printf("free: %s\n", join(string.(free), ", "))
@printf("apsidal motion: %s\n\n", APSIDAL ? "FITTED" : "fixed at 0 (see header)")

res = fit_orbit(data, star1, star2; elements = elements, flux_ratio = F_H,
                free = free, method = METHOD, min_num_live_points = NLIVE,
                verbose = true,
                # ---------------------------------------------------------------------
                # Boxes from PUBLISHED uncertainties, not from habit. Nested sampling pays
                # for prior volume in wall clock, and an unjustified box is a slow fit, not
                # a cautious one. Against the previous set these cut the volume by ~2.3e3.
                #
                #  a       1.25-1.86 mas   a sini = 23.14 +/- 0.88 Rsun [Robinette &
                #                          Aufdenberg 2015] with parallax 13.06 +/- 0.70
                #                          [van Leeuwen 2007] gives 1.555 +/- 0.102; ~3 sigma
                #  i       105-128 deg     116 [Aufdenberg+2015]. No published sigma, so this
                #                          stays generous — but MOST sees a GRAZING eclipse
                #                          (Desmet+2009), which pins i far better than this
                #                          box; tighten if you adopt that constraint.
                #  Omega   290-330 deg     309.938. Was the FULL circle, which is 9x the
                #                          volume for no reason: the node is not in doubt,
                #                          only its precision.
                #  omega   225-285 deg     255 at T0. Note it DRIFTS ~21 deg across the data
                #                          through apsidal motion, so this is the value at
                #                          the reference epoch, not a constant.
                #  e       0.02-0.16       CONTESTED: 0.065 +/- 0.015 [Wages & Aufdenberg],
                #                          0.116 +/- 0.004 [Robinette & Aufdenberg],
                #                          0.118-0.119 [Aufdenberg+2015]. This spans all of
                #                          them at >3 sigma. The lower edge is 0.02 rather
                #                          than 0 deliberately: at e = 0 omega is undefined,
                #                          so including it hands the sampler a degenerate
                #                          direction. If the posterior piles up at 0.02, the
                #                          data prefer a circular orbit — that is a result,
                #                          not a bound to relax.
                #  T0      +/-P/2          left at the default. Omega is now restricted to
                #                          one branch, so a full period is exactly one
                #                          fundamental domain of (Omega,T0) -> (Omega+180,
                #                          T0+P/2).
                #  c1/c2   +/-4 sigma      R1 = 7.47 +/- 0.54, R2 = 3.74 +/- 0.53 Rsun
                #                          [Tkachenko+2016 Tab 3] at d = 77 pc, as DIAMETERS
                #                          0.902 +/- 0.065 and 0.452 +/- 0.064 mas.
                #  f       0.06-0.45       0.196 +/- 0.063 from the same radii and
                #                          Teff = 25300 +/- 500 / 20585 +/- 850 K through a
                #                          PLANCK ratio at 1.65 um. Rayleigh-Jeans would give
                #                          0.204, 4 % high — hc/lambda k = 8720 K is not
                #                          negligible against 20-25 kK.
                #  domega  5e-3 - 1.1e-2   only with APSIDAL=1. Covers U = 139 +/- 6 yr
                #                          (7.09e-3 deg/d) and U = 105 +/- 2 yr (9.39e-3),
                #                          the two competing determinations. Was
                #                          (-0.05, 0.05) from the library default: 17x wider
                #                          and admitting retrograde apsidal motion.
                bounds = Dict(:a => (1.25, 1.86), :i => (105.0, 128.0),
                              :Omega => (290.0, 330.0), :omega => (225.0, 285.0),
                              :e => (0.02, 0.16),
                              :c1_diameter => (0.64, 1.17), :c2_diameter => (0.20, 0.71),
                              # Hestroffer exponent I(mu) ~ mu^ld1. Tabulated H-band value
                              # for both is 0.15 (spica_params.jl); the box spans uniform
                              # (0, no darkening) to well past any hot-star H-band
                              # prediction. The lower edge is 0 rather than negative: mu^ld1
                              # with ld1 < 0 is limb BRIGHTENING, which does not happen for
                              # a B star in the near-IR.
                              # Hestroffer exponent I(mu) ~ mu^ld1, primary only (see the
                              # resolution argument above). Lower edge 0 rather than
                              # negative: mu^ld1 with ld1 < 0 is limb BRIGHTENING, which a
                              # B star does not do in the near-IR.
                              :c1_ld1 => (0.0, 0.6),
                              :f => (0.06, 0.45),
                              :domega => (5.0e-3, 1.1e-2)))

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
