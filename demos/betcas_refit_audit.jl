#!/usr/bin/env julia
# =======================================================================================
# beta Cas refit: what this dataset can and cannot settle
# =======================================================================================
# The intended headline was "predicted limb darkening makes beta come out consistent with the
# radiative 0.25". That is NOT what the data say, and this script is the honest record of why.
# It exists because the negative result is more useful than the number would have been.
#
# The information-content question — does removing ld1 break the beta degeneracy — IS answered
# cleanly, at fixed parameters, by demos/betcas_fisher_predicted_ld.jl. That answer stands.
# What does not follow from it is a converged beta from this one epoch.
#
# WHAT IS ESTABLISHED
#
#   * The rapid-rotator model is strongly REQUIRED. chi2r falls 76.8 (uniform disc, 1 param)
#     -> 39.6 (limb-darkened disc, 2 params) -> 3.44 (rotator, 6 params). An 11.5x improvement
#     over the best circular model, so the oblateness and gravity darkening are doing real
#     work and this is not six parameters fitting noise. (All three over the SAME observables
#     and denominator — see `chi2r_of`; an earlier version of this script did not, and
#     overstated the factor as 16x.)
#   * Incidentally, the best-fit power-law exponent for the circular disc is alpha = 0.20,
#     against the alpha ~ 0.16 the Korg H-band grid PREDICTS. A loose but free consistency
#     check on the atmosphere grid.
#   * With FITTED limb darkening and a beta = 0.25 start, the fit lands essentially on the
#     published solution: inc = 20.0 deg, beta = 0.148 under ELR (published 0.146), and a
#     derived mass of 2.49 Msun. So the machinery reproduces the known answer.
#   * With PREDICTED limb darkening the optimiser LEAVES that basin, to inc = 38.6 deg,
#     beta = 0.350 and a derived 6.7 Msun. chi2r is lower there (3.435 vs 3.563) but the
#     solution is physically excluded. This is a problem for the method ON THIS DATASET and it
#     is not papered over here: either predicted LD genuinely prefers the flatter geometry, or
#     something else in the model is absorbing the difference. It is not grid clamping — the
#     surface spans Teff 5888-7206 K and logg 3.29-3.75, inside the grid's 5400-7600 / 2.9-4.1.
#     Deciding between those needs more data, not a tighter prior.
#   * One configuration (predicted LD + von Zeipel) fails in OptimPack's line search with an
#     AssertionError. The `fit_parametric` docstring already warns that these objectives stop
#     on x/f tests rather than on the gradient; this is the same fragility.
#
# WHAT THIS DATASET CANNOT SETTLE
#
#   * chi2r = 3.44 at best. Either the errors are underestimated by ~1.85x or the model is
#     missing something real (a spot, a companion, calibration systematics). Formal sigmas and
#     formal Delta-chi2 between models both inherit that.
#   * The chi2 surface has ~18 distinct local minima over 27 starts, with chi2r from 3.44 to
#     6.07 (plus five starts that diverge past 6000). Only 4 of 27 reach the best basin.
#   * The deepest basin is PHYSICALLY EXCLUDED: inc ~ 38-48 deg gives vsini 130-158 km/s
#     against beta Cas's spectroscopic ~70, and a derived mass of 6.7-8.8 Msun against ~1.9.
#   * The published-style basin (inc ~ 20 deg, fev ~ 0.85, beta ~ 0.12) IS a genuine local
#     minimum, at chi2r = 3.5713 — only 0.13 worse. Which basin the optimiser reaches depends
#     on the STARTING beta: from beta = 0.25 it stays at inc = 20.0; from beta = 0.146 it
#     walks out to inc = 38.
#   * A vsini prior does not rescue it. Moving to the published basin costs Delta-chi2/2 ~ 81
#     in log-likelihood while a 70 +- 8 km/s prior can save at most ~40. Tightening the prior
#     until it wins would be choosing the answer, not measuring it.
#
# WHAT WOULD BE NEEDED: more epochs (one epoch of 540 V2 + 720 T3phi does not constrain six
# parameters through a degenerate valley), an explanation for chi2r = 3.44 before any formal
# error is quoted, and the vsini and mass constraints applied with rescaled errors.
#
# The information-content question is answered separately and cleanly, at fixed parameters,
# by demos/betcas_fisher_predicted_ld.jl — which is why a Fisher analysis was the right way to
# decide whether predicted LD helps, and a single fit was not.
#
#   julia --project=demos demos/betcas_refit_audit.jl
# =======================================================================================

using ROTIR, Zygote, Printf, LinearAlgebra, Statistics

const D = readoifits_multiepochs([joinpath(@__DIR__, "data",
    "MEDIAN5.MIRCX_L2.2025Oct30.HD_432.MIRCX_IDL.bet_Cas.AVG10m.oifits")]; T = Float64)[1, 1]
const DATA = [D]
const TESS = tessellation_healpix(4, T = Float64)
const TE   = [0.0]
const BAND = band_of(D)
const GRID = joinpath(@__DIR__, "data", "betcas_H_korg.fits")
const DPC, PROT = 16.8, 1/1.12

# MUST match what `fit_parametric` minimises, or the circular baselines below are not
# comparable to the rotator's chi2r. It optimises V2 + T3amp + T3phi (OI_DEFAULT_WEIGHTS =
# [1,1,1,0,0,0,0]) and divides by `npts - nfree`. An earlier version of this script used
# V2 + T3phi over 1260 points while quoting the rotator's figure over 1980 — apples to
# oranges, and it overstated the rotator's advantage as 16x when it is 11.5x.
const NPTS = D.nv2 + D.nt3amp + D.nt3phi
chi2r_of(cv, npar) = ROTIR.cvis_chi2(cv, D) / max(NPTS - npar, 1)

println("="^92)
println("[1] Is a rapid rotator required at all?  (circular baselines)")
println("="^92)
# In functions, not top-level loops: Julia's soft scope makes an assignment inside a
# top-level `for` local to it, so the accumulator would never be updated.
function scan_ud()
    best = (Inf, 0.0)
    for diam in 1.70:0.002:2.60
        c = chi2r_of(visibility_ud([diam], D.uv), 1); c < best[1] && (best = (c, diam))
    end
    best
end
function scan_ld()
    best = (Inf, 0.0, 0.0)
    for diam in 1.70:0.002:2.60, α in 0.0:0.02:0.8
        c = chi2r_of(visibility_ldpow([diam, α], D.uv), 2); c < best[1] && (best = (c, diam, α))
    end
    best
end
bud = scan_ud(); bld = scan_ld()
@printf("  uniform disc        (1 param): chi2r = %8.3f   diam = %.3f mas\n", bud[1], bud[2])
@printf("  limb-darkened disc  (2 param): chi2r = %8.3f   diam = %.3f mas, alpha = %.2f\n",
        bld[1], bld[2], bld[3])
println("  (the Korg H-band grid predicts alpha ~ 0.16 — a loose but free check on the grid)")
@printf("  all three chi2r over the SAME %d points (V2+T3amp+T3phi), chi2/(npts-npar)\n", NPTS)

base(law, ldt) = default_star_params(2; rpole = 0.849, d = DPC, frac_escapevel = 0.92,
    rotation_period = PROT, tpole = 7208.0, inclination = 19.9, position_angle = -7.09,
    beta = 0.25, gravity_law = law, ldtype = ldt, ld1 = 0.21, ld2 = 0.0)

println("\n" * "="^92)
println("[2] How many minima does the rotator's chi2 surface have?")
println("="^92)
bp = base(1, 3); free = ["rpole","omega","inc","PA","beta","ld1"]
function multistart(bp, free)
    sols = Tuple{Float64,Vector{Float64}}[]
    for rp in (0.80, 0.849, 0.95), fv in (0.5, 0.75, 0.92), ic in (20.0, 45.0, 70.0)
        try
            θ, c, _ = fit_parametric(DATA, TESS, TE, bp; θ0 = [rp,fv,ic,-7.09,0.25,0.21,0.0],
                                     free = free, intensity_model = :planck, band = BAND,
                                     maxiter = 300)
            push!(sols, (c, θ))
        catch e; end
    end
    sols
end
sols = multistart(bp, free)
cs = first.(sols)
@printf("  %d starts; chi2r min %.4f, median %.4f, max %.4g\n",
        length(sols), minimum(cs), median(cs), maximum(cs))
@printf("  %d of %d land within 0.1%% of the best\n",
        count(<(minimum(cs)*1.001), cs), length(sols))
println("  distinct chi2r (1e-4): ", sort(unique(round.(cs, digits = 4))))
bi = argmin(cs); θb = sols[bi][2]
qb = derived_quantities(merge(bp, (rpole=θb[1], frac_escapevel=θb[2], inclination=θb[3])))
@printf("  deepest basin: rpole %.4f fev %.4f inc %.2f beta %.4f -> M = %.2f Msun, vsini = %.0f km/s\n",
        θb[1], θb[2], θb[3], θb[5], qb.mass, qb.vsini)
println("  beta Cas is ~1.9 Msun with vsini ~ 70 km/s, so the DEEPEST basin is excluded.")

println("\n" * "="^92)
println("[3] Predicted vs fitted limb darkening: what happens to beta")
println("="^92)
if !isfile(GRID)
    @warn "betcas_H_korg.fits missing — run demos/betcas_fisher_predicted_ld.jl first"
else
    prov = TabulatedProvider(load_intensity_grid(GRID))
    @printf("  %-30s %9s %9s %9s %8s\n", "configuration", "chi2r", "beta", "inc", "M/Msun")
    function ld_comparison(prov, free)
    for (nm, law, pred) in (("fitted power-law LD, vZ",1,false), ("fitted power-law LD, ELR",2,false),
                            ("predicted (atmosphere), vZ",1,true), ("predicted (atmosphere), ELR",2,true))
        bpx = base(law, pred ? 0 : 3)
        L = pred ? parametric_layout(ldtype=0, distance_free=true) : parametric_layout()
        fr = pred ? ["rpole","omega","inc","PA","beta"] : free
        θ0 = pred ? [0.849,0.92,19.9,-7.09,0.25,DPC] : [0.849,0.92,19.9,-7.09,0.25,0.21,0.0]
        try
            θ, c, _ = fit_parametric(DATA, TESS, TE, bpx; θ0=θ0, free=fr, layout=L,
                          provider = pred ? prov : nothing,
                          intensity_model = pred ? :linear : :planck, band=BAND, maxiter=500)
            q = derived_quantities(layout_merge(L, bpx, θ))
            @printf("  %-30s %9.4f %9.4f %9.2f %8.2f\n", nm, c, θ[θindex(L,"beta")], θ[3], q.mass)
        catch e
            @printf("  %-30s  line search failed (%s)\n", nm, typeof(e))
        end
    end
    end
    ld_comparison(prov, free)
    println("\n  Fitted LD reproduces the published solution (inc 20.0, beta 0.148 vs 0.146,")
    println("  M 2.5 Msun). Predicted LD LEAVES that basin for inc 38.6 / 6.7 Msun at a lower")
    println("  chi2r — physically excluded, and not explained by grid clamping. The two rows")
    println("  are therefore NOT an apples-to-apples beta comparison, and no conclusion about")
    println("  beta should be drawn from them.")
end
println("\n" * "="^92)
println("Bottom line: the model is required, the parameters are not determined by one epoch,")
println("and no formal error from this fit should be quoted until chi2r = 3.44 is explained.")
println("="^92)
