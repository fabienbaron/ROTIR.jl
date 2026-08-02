#!/usr/bin/env julia
# Is the bootstrap σ honest? — calibration against simulated truth.
#
#   julia --project=demos -t auto demos/bootstrap_calibration.jl
#
# Cross-checking ROTIR against another code only shows the two agree; it cannot show
# either is right. The only real test is to generate data from ROTIR's *own* forward model
# at known parameters — using ρ Cas's actual uv coverage and error bars — and ask whether
# the quoted uncertainty matches the scatter the parameters actually have.
#
# Two regimes, because they answer different questions:
#
#   A. independent Gaussian noise, exactly as the error bars claim.
#      Here the analytic/Monte-Carlo answer is correct by construction, and a calibrated
#      bootstrap must reproduce it. Ratio ≈ 1 expected.
#
#   B. the same, plus a correlated multiplicative error per (baseline, exposure) block —
#      a calibration systematic, the thing real interferometric data actually suffers from
#      and the reason the block bootstrap exists. A parametric Monte Carlo that trusts the
#      error bars stays blind to it; the block bootstrap should widen. Ratio ≈ 1 expected
#      for the block bootstrap, ≪ 1 for the naive alternative.
#
# Cost: (NTRUTH + 2·NBOOT) fits. At n=3/Float64 a warm-started fit is ~10-30 s of core
# time, so the defaults below are an hour or two on many cores. Reduce NTRUTH first.

using ROTIR, Zygote, Printf, Statistics, Random
include(joinpath(@__DIR__, "posterior_utils.jl"))

const NTRUTH   = parse(Int, get(ENV, "NTRUTH", "60"))   # noise realisations = true scatter
const NBOOT    = parse(Int, get(ENV, "NBOOT",  "200"))  # bootstrap replicates, ONE realisation
const CAL_ERR  = 0.03    # regime B: 3% rms correlated error per block
const SEED     = 20260801

# ── Setup ───────────────────────────────────────────────────────────────────
data0   = readoifits_multiepochs(["./demos/data/rho_Cas_example.oifits"]; T = Float64)[1, :]
tepochs = [0.0]
tessels = tessellation_healpix(3, T = Float64)

base = (surface_type = 2, rpole = 1.25, tpole = 4000.0,
        ldtype = 3, ld1 = 1.75, ld2 = 0.0,
        inclination = 0.0, position_angle = 0.0, rotation_period = 1e6,
        beta = 0.08, frac_escapevel = 0.0, B_rot = 0.0)

θ_true = [1.2528, 0.0, 0.0, 0.0, 0.08, 1.7454, 0.0]     # rpole ω inc PA β ld1 ld2
free   = ["rpole", "ld1"]
ifree  = parametric_free_indices(free)

# ── Noiseless model observables on the real uv coverage ─────────────────────
truth_params = merge(base, (rpole = θ_true[1], ld1 = θ_true[6]))
stars = create_star_multiepochs(tessels, truth_params, tepochs, T = Float64)
setup_oi!(data0, stars)
xmap = parametric_temperature_map(truth_params, stars[1])

function noiseless(data, stars, xmap)
    out = deepcopy(data)
    for e in eachindex(out)
        cvis = poly_to_cvis(xmap, stars[e])
        out[e].v2 .= cvis_to_v2(cvis, out[e].indx_v2)
        _, t3a, t3p = cvis_to_t3(cvis, out[e].indx_t3_1, out[e].indx_t3_2,
                                 out[e].indx_t3_3; T = Float64)
        out[e].t3amp .= t3a
        out[e].t3phi .= t3p
    end
    return out
end

model_data = noiseless(data0, stars, xmap)
blocks     = epoch_blocks(model_data)
@printf("%d points, %d :config blocks, truth: diameter %.4f mas, ld1 %.4f\n",
        sum(d -> d.nv2 + d.nt3amp + d.nt3phi, model_data), sum(length, blocks),
        2θ_true[1], θ_true[6])

# ── Noise models ────────────────────────────────────────────────────────────
# A: independent, exactly as the error bars claim (OITOOLS' perturb_data).
realise_A(rng) = [perturb_data(d; rng = rng) for d in model_data]

# B: the same, plus one correlated multiplicative factor per block. V² and T3amp scale
# with the factor, T3phi is left alone — a transmission/calibration error, not a phase one.
function realise_B(rng)
    d = realise_A(rng)
    for e in eachindex(d)
        b = blocks[e]
        for k in 1:length(b)
            f = 1 + CAL_ERR * randn(rng)
            for i in b.idx_v2[k];  d[e].v2[i]    *= f;   end
            for i in b.idx_t3[k];  d[e].t3amp[i] *= f;   end
        end
    end
    return d
end

fit1(d, θ0) = fit_parametric(d, tessels, tepochs, base;
                             θ0 = θ0, free = free, maxiter = 400)[1][ifree]

# ── True scatter: refit many independent realisations ───────────────────────
function true_scatter(realise, label)
    S = fill(NaN, NTRUTH, length(ifree))
    tick = progress_counter(NTRUTH; label = "truth fits:")
    Threads.@threads for i in 1:NTRUTH
        try
            S[i, :] = fit1(realise(Xoshiro(SEED + 1000i)), θ_true)
        catch e
            e isa InterruptException && rethrow(e)
        end
        tick()
    end
    ok = [all(isfinite, @view S[i, :]) for i in 1:NTRUTH]
    σ  = [0.5 * (quantile(S[ok, j], 0.84) - quantile(S[ok, j], 0.16)) for j in eachindex(ifree)]
    @printf("  %s truth: %d/%d fits, σ = %s\n", label, count(ok), NTRUTH,
            string(round.(σ, sigdigits = 3)))
    return σ
end

# ── Run both regimes ────────────────────────────────────────────────────────
for (label, realise) in (("A independent noise ", realise_A),
                         ("B + 3% correlated   ", realise_B))
    println("\n=== regime $label ===")
    σ_true = true_scatter(realise, label)

    one = realise(Xoshiro(SEED))              # a single dataset, as an observer would have
    θ̂   = fit_parametric(one, tessels, tepochs, base; θ0 = θ_true, free = free,
                         maxiter = 400)[1]

    b = bootstrap_parametric(one, tessels, tepochs, base; θ0 = θ̂, free = free,
                             nboot = NBOOT, seed = SEED, sigma_clipping = 4.5,
                             maxiter = 400, verb = false)
    # the naive alternative: resample the *values* from the error bars, not the blocks
    mc = bootstrap_parametric(rdata -> fit1(rdata, θ̂), one; nboot = NBOOT, seed = SEED,
                              x_opt = θ̂[ifree], list_free_params = free,
                              mode = :weights, granularity = :point, verb = false)

    for (j, p) in enumerate(free)
        @printf("  %-6s  true σ %.5g | block bootstrap %.5g (ratio %.2f) | i.i.d. %.5g (ratio %.2f)\n",
                p, σ_true[j], b.sigma[j], b.sigma[j] / σ_true[j],
                mc.sigma[j], mc.sigma[j] / σ_true[j])
    end
end

println("\nRatio 1.00 = calibrated. In regime B the block bootstrap should stay near 1")
println("while the i.i.d. one falls well below it — that gap is the correlated error")
println("budget, and it is exactly what a posterior computed from the quoted error bars")
println("also misses.")
