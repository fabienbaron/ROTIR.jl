#!/usr/bin/env julia
# Validation for parametric_gradient.jl
#   (1) hand-coded leaf derivatives (von Zeipel, LD) vs finite differences  [no Zygote]
#   (2) forward consistency: build_parametric_logπ vs spheroid_parametric_f
#   (3) end-to-end Zygote.gradient vs central finite differences
#   (4) tpole null-space gate (linear intensity) and Planck breaking it
# Runs in Float64 for FD accuracy (the runtime paths are type-generic / Float32-native).

using ROTIR
using FiniteDifferences
using LinearAlgebra
using Printf

const FDM = central_fdm(5, 1)
# `scale` floors the denominator so that genuinely-zero derivatives (e.g. the
# scale-invariant dx/drpole, or d/dld2 for laws that don't use ld2) report ~0
# instead of a spurious O(1) ratio of two numerical-noise vectors.
relerr(a, b; scale = 0.0) = norm(a .- b) / max(norm(a), norm(b), scale, eps())
status(r) = r < 1e-4 ? "✓" : (r < 1e-2 ? "~" : "✗")
npass = Ref(0); nfail = Ref(0); zygote_ran = Ref(false)
function check(label, r; tol = 1e-4)
    ok = r < tol
    ok ? (npass[] += 1) : (nfail[] += 1)
    @printf("  %-46s relerr=%.2e  %s\n", label, r, status(r))
end

# ─── Setup (Float64) ─────────────────────────────────────────────────────────
oifitsfiles = ["./demos/data/2011Sep02.lam_And_prepped.oifits"]
data_all = readoifits_multiepochs(oifitsfiles; T = Float64)
data = data_all[1, :]
tepochs = Float64.([d.mean_mjd for d in data]); tepochs .-= tepochs[1]
n = 3; tessels = tessellation_healpix(n, T = Float64)

base = (surface_type = 2, rpole = 1.37, tpole = 4800.0, ldtype = 3,
        ld1 = 0.23, ld2 = 0.1, inclination = 78.0, position_angle = 24.0,
        rotation_period = 54.8, beta = 0.15, frac_escapevel = 0.6, B_rot = 0.0)

sinθ = sin.(tessels.unit_spherical[:, 5, 2])
cosθ = cos.(tessels.unit_spherical[:, 5, 2])
rpole, fev, inc, PA, β, ld1, ld2, tpole = 1.37, 0.6, 78.0, 24.0, 0.15, 0.23, 0.1, 4800.0

# ─── (1a) von Zeipel leaf derivatives ────────────────────────────────────────
println("\n[1a] von Zeipel map derivatives (analytic vs FD)")
x, dx_drp, dx_dfev, dx_dβ, dx_dtp = vonzeipel_map_and_derivs(rpole, fev, β, tpole, sinθ, cosθ)
sc = 1e-3 * norm(x)   # derivative-agreement scale (dx/drpole ≡ 0: R is scale-invariant)
check("dx/drpole (≈0, scale-invariant)", relerr(dx_drp, FDM(r -> vonzeipel_map(r, fev, β, tpole, sinθ, cosθ), rpole); scale = sc))
check("dx/dfev",   relerr(dx_dfev, FDM(f -> vonzeipel_map(rpole, f,  β, tpole, sinθ, cosθ), fev);   scale = sc))
check("dx/dbeta",  relerr(dx_dβ,   FDM(b -> vonzeipel_map(rpole, fev, b, tpole, sinθ, cosθ), β);     scale = sc))
check("dx/dtpole", relerr(dx_dtp,  FDM(t -> vonzeipel_map(rpole, fev, β, t,     sinθ, cosθ), tpole); scale = sc))

# ─── (1b) LD leaf derivatives (all three laws) ───────────────────────────────
println("\n[1b] LD map derivatives (analytic vs FD)")
_, _, nz = project_geometry(rpole, fev, inc, PA, tessels, 0.0, base)
for lt in (1, 2, 3)
    ld, dld_dnz, dld_dld1, dld_dld2 = ld_and_derivs(nz, lt, ld1, ld2)
    ldsc = 1e-3 * norm(ld)   # d/dld2 ≡ 0 for ldtype 1 & 3 (they don't use ld2)
    check("ldtype=$lt  d/dld1", relerr(dld_dld1, FDM(a -> ld_weight(nz, lt, a, ld2), ld1); scale = ldsc))
    check("ldtype=$lt  d/dld2", relerr(dld_dld2, FDM(a -> ld_weight(nz, lt, ld1, a), ld2); scale = ldsc))
    check("ldtype=$lt  d/dnz",  relerr(dld_dnz,  FiniteDifferences.grad(FDM, z -> sum(ld_weight(z, lt, ld1, ld2)), nz)[1]); tol = 1e-3)
end

# ─── (2) forward consistency vs spheroid_parametric_f (now LD-aware) ─────────
println("\n[2] forward χ² consistency")
θ = [rpole, fev, inc, PA, β, ld1, ld2]
logπ = build_parametric_logπ(data, tessels, tepochs, base)   # :linear
params = merge(base, (rpole = rpole, frac_escapevel = fev, inclination = inc,
                      position_angle = PA, beta = β, ld1 = ld1, ld2 = ld2))
# NB: my differentiable path keeps ALL pixels (soft vw→0), whereas spheroid_parametric_f
# hard-thresholds at vw>0.01, so agreement is close but not exact — informational (tol 2e-2).
chi2_ref = spheroid_parametric_f(params, tessels, data, tepochs)
check("logπ ≈ -0.5·spheroid_parametric_f", abs(logπ(θ) - (-0.5*chi2_ref)) / abs(0.5*chi2_ref); tol = 2e-2)

# ─── (3) end-to-end Zygote gradient vs FD ────────────────────────────────────
println("\n[3] end-to-end ∇logπ (Zygote vs FD)  —  requires Zygote")
try
    @eval using Zygote
    for lt in (1, 2, 3)
        b = merge(base, (ldtype = lt,))
        lp = build_parametric_logπ(data, tessels, tepochs, b)
        g_ad = Zygote.gradient(lp, θ)[1]
        g_fd = FiniteDifferences.grad(FDM, lp, θ)[1]
        check("ldtype=$lt  full ∇ (7 params)", relerr(g_ad, g_fd); tol = 1e-3)
    end
    # ─── (4) tpole null-space gate + Planck ──────────────────────────────────
    println("\n[4] tpole degeneracy gate")
    θ8 = vcat(θ, tpole)
    lp_lin = build_parametric_logπ(data, tessels, tepochs, base;
                                   intensity_model = :linear, tpole_free = true)
    g_lin = Zygote.gradient(lp_lin, θ8)[1]
    check("linear ∇_tpole ≈ 0 (null-space)", abs(g_lin[8]) / (norm(g_lin) + eps()))
    λH = 1.6e-6  # H band (m)
    lp_pl = build_parametric_logπ(data, tessels, tepochs, base;
                                  intensity_model = :planck, band = λH, tpole_free = true)
    g_pl = Zygote.gradient(lp_pl, θ8)[1]
    g_pl_fd = FiniteDifferences.grad(FDM, lp_pl, θ8)[1]
    check("planck full ∇ (8 params)", relerr(g_pl, g_pl_fd); tol = 1e-3)
    @printf("  planck ∇_tpole = %.4e (should be ≠ 0)\n", g_pl[8])
    zygote_ran[] = true   # runtests.jl asserts this, so a skip can't masquerade as a pass
catch e
    @warn "Zygote section skipped" exception = (e, catch_backtrace())
end

@printf("\n=== %d passed, %d failed ===\n", npass[], nfail[])
