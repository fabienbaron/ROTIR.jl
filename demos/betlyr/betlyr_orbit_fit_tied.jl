# beta Lyrae — orbit fit through ROTIR's generic path, with the disc PA TIED to the node
#
#   julia --project=demos demos/betlyr/betlyr_orbit_fit_tied.jl
#
#   DATASET=2007 julia --project=demos demos/betlyr/betlyr_orbit_fit_tied.jl  # Zhao nights only
#   DATASET=2013 ...                                                        # MIRC-X nights only
#   NLIVE=800 ...                                                           # tighter run
#   TIED=0    ...                                                           # free pa, for comparison
#   FITP=0    ...                                                           # hold P at the ephemeris
#   METHOD=neldermead ...                                                   # quick smoke test
#
# ---------------------------------------------------------------------------------------
# What is different from betlyr_orbit_fit_ultranest.jl
# ---------------------------------------------------------------------------------------
# That driver fits the bespoke model in `betlyr_model.jl`, in which all ten parameters are
# independent. Here the same system is expressed with `orbit_fit_spec`, which allows a
# parameter to be TIED to another rather than fitted:
#
#     ties = Dict(:c2_pa => "-Omega")
#
# A disc in the orbital plane has its projected major axis along the line of nodes, so its
# position angle is not an independent quantity. Fitting it freely both spends a parameter
# and permits a disc oriented inconsistently with the orbit fitted alongside it. On the
# 2006-2007 data the freely-fitted angle lands 1.9 deg from the tied prediction, so this
# costs almost nothing in fit quality and removes a degree of freedom.
#
# SIGN. `EllipticalGaussian` takes `pa` as the angle whose axis is scaled by `ratio < 1` —
# the MINOR axis — measured from +RA toward +Dec, while Omega is the node measured North
# through East. Major axis on the node line is therefore `-Omega`, not `Omega`. Check this
# before copying the tie to another system: a wrong tie cannot be spotted from the chi2.
#
# NOT A REPLACEMENT for the bespoke driver. This model has no occultation term, and the
# components overlap on the sky at 3 of 7 epochs by half-FWHM (6 of 7 using the Gaussian's
# 2-sigma extent), with the disc in front in October and behind in July. Whatever this fit
# reports, that systematic is still present — see the caveats in docs/src/guides/orbits.md.

using ROTIR, Printf, Statistics, DelimitedFiles
# METHOD defaults to :ultranest, which needs PythonCall to load ROTIRUltraNestExt.
using PythonCall
import OITOOLS: oifits_prep

const HERE    = @__DIR__
const OLDDIR  = joinpath(HERE, "older_data")   # Zhao et al.'s 2006-2007 MIRC nights
const TIED    = get(ENV, "TIED", "1") == "1"
const METHOD  = Symbol(get(ENV, "METHOD", "ultranest"))
const NLIVE   = parse(Int, get(ENV, "NLIVE", "400"))
const DATASET = get(ENV, "DATASET", "all")
const OUT     = joinpath(HERE, "results")
DATASET in ("all", "2007", "2013") || error("DATASET must be all|2007|2013")

# ---------------------------------------------------------------------------------------
# Calibration systematics, per MIRC epoch
# ---------------------------------------------------------------------------------------
# Monnier et al. (2012, ApJ 761, 3 — the Vega paper), §2-3, as implemented in
# ~/SOFTWARE/SplatOI.jl (`mirc_floors`). The 2007 data predate MIRC's photometric-channel
# upgrade and carry ~3x worse calibration, which is PHYSICAL rather than bookkeeping:
#
#   epoch                        mult V²   mult T3amp   add V²    CP floor
#   :mirc4  2007, pre-upgrade      20 %       30 %      2e-4        1 deg
#   :mirc6  2012+, 6 telescopes    6.6 %      10 %      2e-4        1 deg
#
# THIS IS WHY THE TWO DATASETS CANNOT SIMPLY BE CONCATENATED. The median relative V²
# uncertainty already in the files is 20.6 % for 2006/2007 against 13.6 % for 2013, so
# without epoch-appropriate floors the 2013 nights dominate the likelihood by ~9x in χ² for
# a model of identical fidelity — the fit would quietly become a 2013-only fit with the
# older nights along for the ride. Applying the floors at load puts both on one error model.
#
# NOT APPLIED (also absent from SplatOI): Monnier reduces the statistical weight of V²,
# T3amp and CP by the number of spectral channels (8 for MIRC) because they are strongly
# correlated across channels, and down-weights V²/T3amp by a further factor 2 since T3amp is
# derivable from V². Both change the χ² NORMALISATION, so treat absolute χ²ᵣ here as
# uncalibrated — differences between fits on the same data remain meaningful.
mirc_floors(epoch::Symbol; t3phi = 1.0) =
    (min_t3phi_err_add = t3phi,
     min_v2_err_rel    = epoch === :mirc4 ? 0.20 : 0.066,
     min_v2_err_add    = 2e-4,
     min_t3amp_err_rel = epoch === :mirc4 ? 0.30 : 0.10,
     quad = true)

function load_set(dir, epoch)
    fs = sort!(filter(f -> endswith(f, ".oifits"), readdir(dir, join = true)))
    out = OIdata{Float64}[]
    for f in fs
        d = readoifits(f; T = Float64, verbose = false)[1, 1]
        isempty(d.v2_mjd) && continue
        med0 = median(d.v2_err ./ max.(abs.(d.v2), 1e-12))
        oifits_prep(d; mirc_floors(epoch)...)      # mutates in place
        med1 = median(d.v2_err ./ max.(abs.(d.v2), 1e-12))
        @printf("  %-58s  σ(V²)/V² %5.1f%% -> %5.1f%%\n", basename(f), 100med0, 100med1)
        push!(out, d)
    end
    return out
end

data = OIdata{Float64}[]
println("loading with Monnier+2012 calibration floors:")
if DATASET in ("all", "2007")
    println(" 2006-2007 MIRC (:mirc4, pre-upgrade)")
    append!(data, load_set(OLDDIR, :mirc4))
end
if DATASET in ("all", "2013")
    println(" 2013 MIRC-X (:mirc6)")
    append!(data, load_set(HERE, :mirc6))
end
isempty(data) && error("no usable OIFITS for DATASET=$DATASET")
sort!(data, by = d -> minimum(d.v2_mjd))
span = (maximum(maximum(d.v2_mjd) for d in data) -
        minimum(minimum(d.v2_mjd) for d in data)) / 365.25
@printf("\n%d epochs, %d data points, %.2f yr baseline (%.0f orbits)\n",
        length(data), sum(d.nv2 + d.nt3amp + d.nt3phi for d in data), span, span*365.25/12.9414)

# ---------------------------------------------------------------------------------------
# The ephemeris, and why the baseline decides what is fittable
# ---------------------------------------------------------------------------------------
# P and T0 are Z08's, NOT the values the bespoke run recovered. That fit moved P to 12.9616
# while it was free; starting from its output would measure `a`, `i` and `Omega` against an
# orbit the paper does not describe.
#
# BETA LYR'S PERIOD IS NOT CONSTANT. Mass transfer at ~2e-5 Msun/yr lengthens the orbit at
# Pdot ~ 18.93 s/yr = 6.0e-7 d/d. What that costs depends entirely on the baseline:
#
#   DATASET=2007   269 d,   21 orbits  — quadratic term moves phase by ~0.5 deg. Negligible;
#                                        P and dP both stay fixed.
#   DATASET=all    6.7 yr, 189 orbits  — the quadratic term moves phase by ~0.14 orbits,
#                                        i.e. ~50 deg (betlyr_params.jl). Ignoring dP here
#                                        does not inflate chi2 politely, it forces T0 and P
#                                        to absorb a drift they cannot represent.
#
# But MATTERING IS NOT THE SAME AS BEING MEASURABLE. The quadratic ephemeris comes from a
# century of eclipse timings: Ak+2007's C = 3.87265e-6 +/- 0.00369e-6 d/E^2 gives
# dP = 5.9977e-7 +/- 5.7e-10 d/d, i.e. dP is known to 0.10 %. Six years of visibilities
# cannot improve on that by any margin, so dP is APPLIED, not fitted — leaving it free over
# any honest box is 3500 sigma of prior volume bought for nothing.
#
# Same for P. Both are opt-in consistency checks (FITP=1 / FITDP=1): the useful question is
# not "what are they" but "do the other elements move when they are freed".
const DP_PUB = 2 * 3.87265e-6 / 12.913779        # d/d = 5.9977e-7, Ak+2007 via [M18] Eq.(1)
const FITP  = get(ENV, "FITP",  "0") == "1"
const FITDP = get(ENV, "FITDP", "0") == "1"

elements = (a = 0.8958, i = 91.244, Omega = 73.529, omega = 90.0, e = 0.0,
            P = 12.9414, T0 = 2454283.0430, dP = DP_PUB)
donor = UniformDisk(diameter = 0.590)                        # B6-8 II
disc  = EllipticalGaussian(fwhm = 0.5825, ratio = 0.6423, pa = 108.35)

free = Symbol[:a, :i, :Omega, :T0, :f, :c1_diameter, :c2_fwhm, :c2_ratio]
ties = TIED ? Dict(:c2_pa => "-Omega") : Dict{Symbol,String}()
TIED || push!(free, :c2_pa)
FITP  && push!(free, :P)
FITDP && push!(free, :dP)

@printf("\n%s: %d free parameters%s\n", TIED ? "TIED" : "FREE",
        length(free), TIED ? " + c2_pa tied to -Omega" : "")
@printf("P  = %.4f d  %s   (Z08 ephemeris)\n", elements.P, FITP ? "FREE" : "fixed")
@printf("dP = %.3e d/d %s   (%.2f s/yr, Ak+2007)\n", elements.dP,
        FITDP ? "FREE" : "fixed", elements.dP*365.25*86400)
@printf("quadratic term over this baseline: %.3f orbits (%.0f deg)\n",
        0.5*elements.dP*(span*365.25)^2/elements.P,
        360*0.5*elements.dP*(span*365.25)^2/elements.P)

res = fit_orbit(data, donor, disc; elements = elements, flux_ratio = 0.8179,
                free = free, ties = ties, method = METHOD,
                min_num_live_points = NLIVE, verbose = true,
                # `a` is held near the published value so the tie is what is being
                # tested. The P prior is scaled to the physics rather than to what the data
                # might see: the observable is accumulated phase drift n*dP/P, and over 21
                # orbits even 1e-3 d moves phase by well under a degree, so anything wider
                # asserts an uncertainty the eclipse ephemeris contradicts.
                # ---------------------------------------------------------------------
                # Every box below is set from a PUBLISHED uncertainty, not from habit.
                # Nested sampling has to compress the prior down to the posterior, so an
                # unjustified box is paid for in wall clock. Against the defaults these cut
                # the prior volume by ~1.8e6, dominated by the three that were never
                # overridden at all (f and the two component sizes).
                #
                #  param        box            justification
                #  a            +/-0.08 mas    0.865 +/- 0.048 [Z08 Tab.3], ~1.7 sigma;
                #                              the combined fit sits at 0.918
                #  i            88-97 deg      93.5 +/- 1.0 [M18 §6]. NOTE the old (85,95)
                #                              EXCLUDED the +3 sigma edge at 96.5 — a prior
                #                              overruling the literature. Now covers it and
                #                              is still narrower.
                #  Omega        68-80 deg      73.7 +/- 1.0 [M18 §6], ~6 sigma
                #  T0           +/-0.5 d       propagated ephemeris sigma is 0.075 d at
                #                              E=3561 (T0 0.015, E*sP, E^2*sC in quadrature),
                #                              so this is ~6.6 sigma. The default was +/-P/2,
                #                              i.e. 86 sigma. Safe to shrink because Omega is
                #                              restricted to one branch, so the
                #                              (Omega,T0) -> (Omega+180,T0+P/2) mirror is
                #                              already outside the box.
                #  c1_diameter  0.3-1.0 mas    donor UD 0.57; it fills its Roche lobe, so it
                #                              cannot be much larger
                #  c2_fwhm      0.3-1.5 mas    [Z08 §3] 1.04 +/- 0.11 major axis
                #  c2_ratio     0.35-0.90      0.63/1.04 = 0.606, sigma ~0.08
                #  f            0.3-2.0        0.806 published. The default was (0,100) —
                #                              59x the volume, and it allowed a disc
                #                              outshining the donor 100-fold.
                #  P            +/-0.001 d     only used when FITP=1
                bounds = Dict(:a => (0.82, 0.98), :i => (88.0, 97.0),
                              :Omega => (68.0, 80.0),
                              :T0 => (elements.T0 - 0.5, elements.T0 + 0.5),
                              :c1_diameter => (0.30, 1.00),
                              :c2_fwhm => (0.30, 1.50), :c2_ratio => (0.35, 0.90),
                              :f => (0.30, 2.00),
                              :P  => (12.9404, 12.9424),
                              # +/-2e-8 = +/-35 sigma on a quantity known to 0.10 %.
                              # Only reached with FITDP=1, and if the posterior fills
                              # this box the data are not constraining dP at all —
                              # which is the expected answer, not a failure.
                              :dP => (DP_PUB - 2e-8, DP_PUB + 2e-8)))

@printf("\n%-16s %14s\n", "parameter", "fitted")
for (n, v) in zip(res.names, res.params)
    tag = n in keys(ties) ? "  (tied)" : (n in free ? "" : "  (fixed)")
    @printf("%-16s %14.6f%s\n", n, v, tag)
end
@printf("\nchi2/n = %.4f   (V2 %.2f, T3amp %.2f, T3phi %.2f)\n",
        res.chi2_red, res.chi2_split[1], res.chi2_split[2], res.chi2_split[3])

# The tie is only meaningful if the convention is right: print the major axis in the same
# frame as Omega so the two can be compared directly.
jpa = findfirst(==(:c2_pa), res.names)
@printf("disc major axis astro PA = %.2f deg   line of nodes Omega = %.2f deg\n",
        mod(-res.params[jpa], 180), mod(res.params[3], 180))

mkpath(OUT)
# A `neldermead` run writes a best fit but NO posterior, so letting it share the output name
# leaves a directory in which best.txt and posterior.txt come from different fits — and
# silently replaces a converged result with a smoke test. Tag the method.
tag = "$(DATASET)_" * (TIED ? "tied" : "free") * (FITP ? "_fitP" : "") *
      (FITDP ? "_fitdP" : "") * (METHOD === :ultranest ? "" : "_$(METHOD)")
writedlm(joinpath(OUT, "orbitfit_$(tag)_best.txt"), hcat(string.(res.names), res.params))
haskey(res, :posterior) &&
    writedlm(joinpath(OUT, "orbitfit_$(tag)_posterior.txt"), res.posterior)
@info "wrote results to $OUT"
