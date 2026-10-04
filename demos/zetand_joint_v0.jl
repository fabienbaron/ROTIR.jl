#!/usr/bin/env julia
# ζ And: the joint model + imaging fit, with the visibility calibration V0 fitted alongside.
#
#     julia --project=. -t 16 demos/zetand_joint_v0.jl
#     ZETAND_NOUTER=3 ZETAND_XITER=60 julia --project=. -t 16 demos/zetand_joint_v0.jl   # quick
#
# WHAT THIS IS FOR. Every ζ And result so far has fitted a UNIFORM ellipsoid: χ²/n bottoms out
# at 3.6 because such a model cannot represent a spotted star, and everything derived from it
# inherits that misfit — which is why `zetand_alpha_audit.jl` cannot quote an error on α and
# why `zetand_triaxial.jl` has to measure a spot-induced bias on the axis ratio rather than
# ignore it. The joint fit gives the spots a surface to live on and lets the shape move with
# them, which is what the residual χ² has been asking for.
#
# WHY V0 BELONGS HERE. Roettenbacher et al. fit a "visibility at the origin", V0 = 1.02 ± 0.02;
# ROTIR has no equivalent, because `compute_polyflux_and_cvis!` divides by the total flux and
# the model visibility is therefore exactly 1 at zero baseline by construction. Anything that
# scales the MEASURED visibilities uniformly — calibrator-diameter error, coherence loss — has
# nowhere to go and is absorbed into the geometry and the limb darkening instead. Freeing it
# on a uniform model already recovers 1.020–1.022 (`zetand_alpha_audit.jl`); here it is fitted
# jointly, where the map can no longer stand in for it.
#
# HOW V0 REACHES BOTH STEPS, because the alternation would otherwise descend two different
# criteria and converge quietly on neither:
#
#   θ step    `shape_chi2_fg!` carries it in the MODEL, `cvis = F·V0/flux`, and differentiates
#             it: ∂χ²/∂V0 = Σₖ Re(adj_cvis[k]·Fₖ/flux). Validated against finite differences
#             to 4e-10.
#   map step  `image_reconstruct_oi` has no V0, so the DATA are rescaled instead, which is
#             EXACTLY the same χ²:
#                 ((V0²|V|² − v2)/σ)²  ≡  ((|V|² − v2/V0²)/(σ/V0²))²
#             and likewise T3amp against V0³, with T3φ untouched since a real positive scale
#             cannot rotate a phase. Verified equal to 1.6e-9 at V0 = 0.94, 1.00 and 1.07.
#
# THE LIMB DARKENING IS FIXED AT THE ATMOSPHERE'S VALUE, α = 0.22, not fitted. Kurucz ATLAS9
# gives 0.221 at 4750 K / logg 2.5 in H band (`zetand_alpha_audit.jl`), and that file shows a
# free α on an incomplete surface model comes out low — 0.139 on a uniform triaxial figure,
# climbing to 0.19 as the model is completed. With a free map there is nothing to stop α and
# the map trading limb flux, so predicting it is the point rather than a convenience.
#
# T3φ IS WEIGHTED IN, unlike every uniform-ellipsoid fit in this repo. A centrosymmetric model
# predicts closure phases of 0° or 180° only, so including them there drags the fit; here the
# asymmetry is exactly what the map is for.
#
# THE REGULARISATION WEIGHT IS NOT INHERITED, IT IS MEASURED. `joint_reconstruct_oi`'s default
# of 10 was calibrated on RW Cep, and its own comment says the balance against χ² scales with
# the number of data points. Used here unchanged it gives χ²/n = 1.05 and a map spanning
# 0–21539 K on a 4400 K star: the map absorbing noise, and a χ² — and a V0 — that mean nothing.
# So this file runs the ladder first. MEASURED (nouter 4, 14218 data):
#
#   weight   χ²/n   on their errors   map range (K)   sd/mean      V0
#        1   1.09        0.40            0 – 17925     0.946     1.0249
#      100   1.20        0.44            0 – 14134     0.455     1.0272
#     1000   1.46        0.54         1463 –  5798     0.122     1.0332
#    10000   1.73        0.64         3711 –  4699     0.0367    1.0420
#    30000   2.26        0.83         3724 –  4591     0.0285    1.0428
#   100000   4.12        1.51         3981 –  4489     0.0157    1.0372
#
# against their published Epoch A map at 3742–4589 K, sd/mean 0.0352. Three criteria that do
# not share an assumption agree on 1e4–3e4: it reproduces their map's contrast (sd 160.5 vs
# 154.1), it is the knee of the χ² curve, and it brackets the misfit they themselves accepted
# — they inflate their errors by 1.65 to reach χ² ~ 1, i.e. χ²/n ≈ 2.72 on the PUBLISHED
# errors, which is where 3e4 sits. 1e4 is adopted, and the ladder is printed so the choice is
# visible rather than asserted.
#
# WHAT SURVIVES THE CHOICE, and it is the point of the exercise: **V0 moves only from 1.025 to
# 1.043 across four orders of magnitude of regularisation**. A calibration scale should be
# insensitive to the image prior, and it is — so the V0 result does not rest on the weight.

using ROTIR, Printf, Statistics, DelimitedFiles, Random

const PROT = 17.769426
const ZET  = joinpath(pkgdir(ROTIR), "demos", "zetAnd")
const NSIDE_EXP = 3            # 768 tessels, "matching the CHARA Array's spatial resolution"
const ALPHA = 0.22             # Kurucz ATLAS9, H band, 4750 K / logg 2.5
# THEIR phase zero, so the reconstruction lands in THEIR longitude system and the comparison
# needs no hand-tuned shift. T0 is their time of nodal passage from the ellipsoidal light
# curve; the half period is because their ephemeris refers to the OTHER node. Measured once
# (chi2/n 19.34 at T0, 2.11 at T0 + P/2, and a full scan optimal at exactly +0.500 P) and used
# as a definition since. See docs/src/guides/conventions.md.
const THEIR_T0 = 49992.281
const OUR_T0   = THEIR_T0 + PROT/2

getenvi(k, d) = parse(Int, get(ENV, k, string(d)))
getenvf(k, d) = parse(Float64, get(ENV, k, string(d)))

function main()
    files = sort(filter(f -> endswith(f, ".oifits"), readdir(ZET)))
    data  = [readoifits(joinpath(ZET, f); verbose = false)[1] for f in files]
    # MJD passed straight through: `t0` in `base` below carries the phase zero, so there is no
    # arbitrary "first epoch is zero" origin and no fitted offset anywhere in this file.
    ts    = [minimum(d.v2_mjd) for d in data]
    nd    = sum(d.nv2 + d.nt3amp + d.nt3phi for d in data)
    tess  = tessellation_healpix(NSIDE_EXP; T = Float64)
    @printf("%d epochs, %d data (V² + T3amp + T3φ), %d tessels\n\n", length(data), nd, tess.npix)

    base = (surface_type = 1, radius_x = 1.3275, radius_y = 1.2663, radius_z = 1.2120,
            tpole = 4372.0, ldtype = 3, ld1 = ALPHA, ld2 = 0.0, beta = 0.25,
            inclination = 64.7, position_angle = 122.8, rotation_period = PROT, t0 = OUR_T0)

    # Starting point: the uniform-ellipsoid fit of `zetand_triaxial.jl`, and a UNIFORM map —
    # deliberately not their published one, so the result does not inherit it.
    θ0 = [base.radius_x, base.radius_y, base.radius_z, base.inclination, base.position_angle]
    x0 = fill(Float64(base.tpole), tess.npix)
    lo = [1.0, 1.0, 1.0, 30.0,  90.0, 0.90]
    hi = [1.8, 1.8, 1.8, 89.0, 170.0, 1.10]

    nouter = getenvi("ZETAND_NOUTER", 6)
    xiter  = getenvi("ZETAND_XITER", 200)
    titer  = getenvi("ZETAND_TITER", 40)
    regw   = getenvf("ZETAND_REGW", 1e4)        # from the ladder in the header
    ladder = get(ENV, "ZETAND_LADDER", "0") == "1"

    T10 = vec(readdlm(joinpath(ZET, "zetAnd_fig10_epochA_T.txt")))
    @printf("their published Epoch A map: %.0f-%.0f K, sd %.1f, sd/mean %.4f\n",
            minimum(T10), maximum(T10), std(T10), std(T10)/mean(T10))
    @printf("limb darkening FIXED at alpha = %.2f (Kurucz ATLAS9, H band)\n", ALPHA)

    if ladder
        println("\n=== regularisation ladder (ZETAND_LADDER=1) ===")
        @printf("%10s %9s %9s %9s %9s %9s %8s %8s\n", "reg_weight", "chi2/n", "on theirs",
                "map min", "map max", "map sd", "sd/mean", "V0")
        for w in (1.0, 1e2, 1e3, 1e4, 3e4, 1e5)
            xm, θ = joint_reconstruct_oi(copy(x0), vcat(θ0, 1.0), data, tess, base, ts;
                        maxiter_xmap = xiter, maxiter_θ = titer, nouter = 4, reg_weight = w,
                        regularizer = "sobel", θ_lower = lo, θ_upper = hi,
                        fit_v0 = true, verbose = false)
            c = shape_chi2_fg!(zeros(6), zeros(tess.npix), copy(xm), θ, data, tess, base, ts;
                               verbose = false, fit_v0 = true)
            @printf("%10.0f %9.4f %9.4f %9.0f %9.0f %9.1f %8.4f %8.4f\n", w, c/nd,
                    c/nd/1.65^2, minimum(xm), maximum(xm), std(xm), std(xm)/mean(xm), θ[6])
        end
    end

    results = Dict{String,Any}()
    for (tag, fitv0) in (("V0 fixed at 1", false), ("V0 fitted", true))
        println("\n", "="^78, "\n  ", tag, "  (reg_weight = ", regw, ")\n", "="^78)
        θstart = fitv0 ? vcat(θ0, 1.0) : copy(θ0)
        el = @elapsed xmap, θ = joint_reconstruct_oi(copy(x0), θstart, data, tess, base, ts;
                                 maxiter_xmap = xiter, maxiter_θ = titer, nouter = nouter,
                                 reg_weight = regw, regularizer = "sobel",
                                 θ_lower = fitv0 ? lo : lo[1:5],
                                 θ_upper = fitv0 ? hi : hi[1:5],
                                 fit_v0 = fitv0, verbose = false)
        g = zeros(length(θ)); gx = zeros(tess.npix)
        c = shape_chi2_fg!(g, gx, copy(xmap), θ, data, tess, base, ts;
                           verbose = false, fit_v0 = fitv0)
        results[tag] = (; xmap, θ, chi2 = c, el)
        @printf("\n  rx = %.4f  ry = %.4f  rz = %.4f mas\n", θ[1], θ[2], θ[3])
        @printf("  inclination = %.2f deg   PA = %.2f deg\n", θ[4], θ[5])
        fitv0 && @printf("  V0 = %.4f            (Roettenbacher et al.: 1.02 +/- 0.02)\n", θ[6])
        @printf("  major/minor = %.4f     (their Table 6: 1.063 +/- 0.013)\n",
                maximum(θ[1:3]) / minimum(θ[1:3]))
        @printf("  chi2/n = %.4f on published errors, %.4f on their inflated ones\n",
                c/nd, c/nd/1.65^2)
        @printf("  map: %.0f - %.0f K, mean %.0f, sd %.1f, sd/mean %.4f (theirs %.4f)\n",
                minimum(xmap), maximum(xmap), mean(xmap), std(xmap), std(xmap)/mean(xmap),
                std(T10)/mean(T10))
        writedlm(joinpath(ZET, "zetAnd_joint_$(fitv0 ? "v0" : "fixed").txt"), xmap)
    end

    a, b = results["V0 fixed at 1"], results["V0 fitted"]
    println("\n", "="^78)
    @printf("freeing V0: χ²/n %.4f → %.4f  (Δχ² = %.1f for one parameter over %d data)\n",
            a.chi2/nd, b.chi2/nd, a.chi2 - b.chi2, nd)
    @printf("the shape moves by %.2f %% in the major axis and %.2f deg in inclination\n",
            100*abs(maximum(b.θ[1:3])/maximum(a.θ[1:3]) - 1), abs(b.θ[4] - a.θ[4]))
    println("\nmaps written to demos/zetAnd/zetAnd_joint_{fixed,v0}.txt")
    println("CAVEAT: χ²/n is on the PUBLISHED errors. Roettenbacher et al. inflate theirs by")
    println("1.65 to reach χ² ~ 1, so divide by 2.72 before comparing with their numbers.")
end

main()
