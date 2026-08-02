#!/usr/bin/env julia
# Block-bootstrap uncertainties on a parametric fit of ρ Cas (CHARA).
#
#   julia --project=demos -t auto demos/rho_cas_bootstrap.jl
#
# Fits the angular size and limb darkening of a resolved single star, then resamples the
# data in blocks of (MJD, telescope configuration) and refits, so the error bars do not
# assume the quoted OIFITS errors are correct or uncorrelated.
#
# Companion: demos/rapid_rotator_betCas_nuts.jl samples the posterior of the same model
# with NUTS. The posterior trusts the error bars; the bootstrap does not. Run both.

using ROTIR, Zygote, Printf, Statistics

const NBOOT = 100          # ≥1000 for publication-grade percentiles; 100 runs in minutes

# ── Data ────────────────────────────────────────────────────────────────────
data    = readoifits_multiepochs(["./demos/data/rho_Cas_example.oifits"]; T = Float64)[1, :]
tepochs = [0.0]            # single epoch: no rotational phase information
tessels = tessellation_healpix(3, T = Float64)   # 768 tessels

blocks = epoch_blocks(data)
@printf("%d epoch(s), %d data points, %d :config blocks\n",
        length(data), sum(d -> d.nv2 + d.nt3amp + d.nt3phi, data), sum(length, blocks))

# ── Model ───────────────────────────────────────────────────────────────────
# Rapid-rotator surface with ω frozen at 0 is exactly a limb-darkened sphere: the radius
# reduces to rpole everywhere, and with a uniform von Zeipel map the inclination, position
# angle, β and tpole are inert. Fitting them would only add flat directions.
base = (surface_type = 2, rpole = 1.25, tpole = 4000.0,
        ldtype = 3, ld1 = 1.75, ld2 = 0.0,          # Hestroffer power law I ∝ μ^ld1
        inclination = 0.0, position_angle = 0.0, rotation_period = 1e6,
        beta = 0.08, frac_escapevel = 0.0, B_rot = 0.0)

θ0   = [1.25, 0.0, 0.0, 0.0, 0.08, 1.75, 0.0]       # rpole ω inc PA β ld1 ld2
free = ["rpole", "ld1"]

# ── Full-data fit ───────────────────────────────────────────────────────────
# The χ² of a resolved star is multimodal in diameter (this one is resolved past its first
# null: closure phases have an rms of 88°). Start close to the expected size — a far start
# converges to a different basin entirely.
θ̂, chi2r, info = fit_parametric(data, tessels, tepochs, base;
                                θ0 = θ0, free = free, maxiter = 400)
@printf("\nFull-data fit: diameter = %.4f mas   ld1 = %.4f   χ²ᵣ = %.3f  (%d evaluations)\n",
        2θ̂[1], θ̂[6], chi2r, info.evaluations)

# ── Bootstrap ───────────────────────────────────────────────────────────────
b = bootstrap_parametric(data, tessels, tepochs, base;
                         θ0 = θ̂, free = free, nboot = NBOOT, seed = 42,
                         sigma_clipping = 4.5, maxiter = 400, verb = true)

@printf("\n%d/%d replicates kept\n", count(b.mask), NBOOT)
@printf("diameter = %.4f +%.4f -%.4f mas\n",
        2b.median[1], 2b.sigma_plus[1], 2b.sigma_minus[1])
@printf("ld1      = %.4f +%.4f -%.4f\n", b.median[2], b.sigma_plus[2], b.sigma_minus[2])
@printf("correlation(diameter, ld1) = %+.3f\n", b.correlation[1, 2])

# ── Is the data correlated? ─────────────────────────────────────────────────
# An i.i.d. (:point) bootstrap ignores the correlation between the wavelength channels and
# baselines of an exposure. If it returns visibly smaller uncertainties than the :config
# block bootstrap, that gap *is* the correlated part of the error budget.
bp = bootstrap_parametric(data, tessels, tepochs, base;
                          θ0 = θ̂, free = free, nboot = NBOOT, seed = 42,
                          granularity = :point, sigma_clipping = 4.5,
                          maxiter = 400, verb = false)
@printf("\nσ(diameter):  :config blocks %.5f mas   vs   :point %.5f mas   (ratio %.2f)\n",
        2b.sigma[1], 2bp.sigma[1], b.sigma[1] / bp.sigma[1])

# ── Caveat worth printing next to any number above ──────────────────────────
if chi2r > 2
    @printf("\nNOTE χ²ᵣ = %.2f: this model does not describe the star (ρ Cas has surface\n", chi2r)
    println("structure a smooth limb-darkened disk cannot reproduce). The bootstrap gives")
    println("honest uncertainties *of this model's parameters*, not of the star's radius.")
end
