#!/usr/bin/env julia
# Block-bootstrap uncertainties on the ρ Cas parametric fit.
#
#   julia --project=demos -t auto demos/rho_cas_bootstrap.jl
#   NBOOT=500 julia --project=demos -t auto demos/rho_cas_bootstrap.jl
#
# On a fresh machine, add the registries first (OptimPackNextGen and friends are not in
# General — see README.md):
#   pkg"registry add General"
#   pkg"registry add https://github.com/emmt/EmmtRegistry"
#
# Resamples the data in blocks of (MJD, telescope configuration) and refits, so the error
# bars do not assume the quoted OIFITS errors are correct or uncorrelated. This is the one
# method in demos/rho_cas_compare.jl that questions the error bars themselves; the three
# samplers all take them at face value.
#
# It is also the one method that cannot leave the basin it starts in — every replicate is
# warm-started at the full-data θ̂, deliberately, so the spread measures the data rather
# than the optimiser. Run demos/rho_cas_basins.jl first to know which basin that is.

using ROTIR, Zygote, Printf, Statistics
include(joinpath(@__DIR__, "rho_cas_model.jl"))
include(joinpath(@__DIR__, "posterior_utils.jl"))

const NBOOT = parse(Int, get(ENV, "NBOOT", "200"))

describe_model()
blocks = epoch_blocks(DATA)
@printf("%d :config blocks over %d epoch(s)\n", sum(length, blocks), length(DATA))

# ── Full-data fit ───────────────────────────────────────────────────────────
(θ̂, chi2r, info), t_fit = timed(() ->
    fit_parametric(DATA, TESSELS, TEPOCHS, BASE; θ0 = THETA0, free = FREE_NAMES,
                   maxiter = 400))
@printf("\nfull-data fit (%.1f s, %d evaluations): χ²ᵣ = %.4f\n", t_fit, info.evaluations, chi2r)
for (j, l) in enumerate(LABELS)
    @printf("  %-8s %.5f\n", l, θ̂[IFREE][j])
end

# ── Bootstrap ───────────────────────────────────────────────────────────────
# No σ-clipping by default. The percentile σ here is ~1e-3 mas, so a 4.5σ window is
# ~0.003 mas wide and the three-pass clip throws away a third of the replicates — it is
# rejecting the distribution's own tail, not outliers. Set SIGMA_CLIP to re-enable it.
const SIGMA_CLIP = haskey(ENV, "SIGMA_CLIP") ? parse(Float64, ENV["SIGMA_CLIP"]) : nothing

b, wall = timed(() -> bootstrap_parametric(DATA, TESSELS, TEPOCHS, BASE;
    θ0 = θ̂, free = FREE_NAMES, refit_full = false, nboot = NBOOT, seed = 42,
    lb = [0.5, 0.0, 0.0, -180.0, 0.0, 0.0, -1.0], ub = [4.0, 0.99, 180.0, 180.0, 1.0, 2.0, 1.0],
    sigma_clipping = SIGMA_CLIP, maxiter = 400, verb = true))

@printf("\nwall time %.1f s   %d/%d replicates kept\n", wall, count(b.mask), NBOOT)
let d = [abs(b.samples[i, 1] - b.median[1]) / max(b.sigma[1], eps()) for i in 1:NBOOT if b.mask[i]]
    @printf("replicate spread: %d beyond 4.5σ, %d beyond 10σ (of %d) — heavy tails here are\n",
            count(>(4.5), d), count(>(10), d), length(d))
    println("the multimodality, so clipping them narrows the error bar by fiat.")
end
summarise(b.samples[b.mask, :], LABELS; title = "Block bootstrap:")
@printf("correlation:\n%s\n", string(round.(b.correlation, digits = 3)))
save_posterior("bootstrap", b.samples[b.mask, :], LABELS;
               wall_seconds = wall, chi2r = chi2r, nboot = NBOOT,
               nblocks = sum(length, blocks))
corner_plot(b.samples[b.mask, :], LABELS, "bootstrap_corner.png";
            truths = θ̂[IFREE], title = "ρ Cas — block bootstrap")

# ── Is the data correlated? ─────────────────────────────────────────────────
# An i.i.d. (:point) bootstrap ignores the correlation between the wavelength channels and
# baselines of an exposure. If it returns visibly smaller uncertainties than the :config
# block bootstrap, that gap *is* the correlated part of the error budget.
println("\nSecond bootstrap (:point granularity, another $NBOOT fits) — this takes as long")
println("as the first one:")
bp, _ = timed(() -> bootstrap_parametric(DATA, TESSELS, TEPOCHS, BASE;
    θ0 = θ̂, free = FREE_NAMES, refit_full = false, nboot = NBOOT, seed = 42,
    granularity = :point, sigma_clipping = SIGMA_CLIP, maxiter = 400, verb = true))
for (j, l) in enumerate(LABELS)
    @printf("σ(%s):  :config blocks %.6g   :point %.6g   (ratio %.2f)\n",
            l, b.sigma[j], bp.sigma[j], b.sigma[j] / bp.sigma[j])
end

if chi2r > 2
    @printf("\nNOTE χ²ᵣ = %.2f: this model does not describe the star. The bootstrap gives\n", chi2r)
    println("honest uncertainties *of this model's parameters*, not of the star's radius.")
end
