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


# ===========================================================================================
# [P5] The ParametricLayout and the model-atmosphere gradient path
# ===========================================================================================
# What this guards, in order of how badly it fails silently:
#
#  1. `logg_map` MUTATES its output, so Zygote refuses it outright ("Mutating arrays is not
#     supported"). Without its rrule the distance cannot reach a gradient-based fit at all —
#     `d` could sit in the parameter vector contributing nothing.
#  2. The layout must reproduce the historical 7-vector EXACTLY, or every existing caller of
#     `fit_parametric` silently maps names to the wrong slots.
#  3. `d` must be in the NULL SPACE under a provider that ignores logg, and OUT of it under
#     one that reads logg. If the first fails, something other than the atmosphere is leaking
#     the distance into the likelihood; if the second fails, the grid is not being indexed by
#     logg and the whole point is lost.
println("\n[P5] ParametricLayout + the model-atmosphere gradient path")

let L0 = parametric_layout(), L8 = parametric_layout(tpole_free = true)
    check("layout == legacy names",
          L0.names == ["rpole","omega","inc","PA","beta","ld1","ld2"] ? 0.0 : 1.0)
    check("layout(tpole_free) == legacy names",
          L8.names == ["rpole","omega","inc","PA","beta","ld1","ld2","tpole"] ? 0.0 : 1.0)
    check("legacy bounds unchanged",
          (default_parametric_bounds()[1] == [1e-3,0.0,0.0,-180.0,0.0,0.0,-1.0] &&
           default_parametric_bounds()[2] == [Inf,0.99,180.0,180.0,1.0,2.0,1.0]) ? 0.0 : 1.0)
    check("legacy free_indices unchanged",
          (parametric_free_indices(nothing) == collect(1:7) &&
           parametric_free_indices(["rpole","ld1"]) == [1,6]) ? 0.0 : 1.0)
    # ldtype = 0 drops both LD coefficients, because the provider owns the mu dependence and
    # a fitted coefficient would be perfectly unconstrained.
    Ld = parametric_layout(ldtype = 0, distance_free = true)
    check("ldtype=0 + distance_free == [rpole,omega,inc,PA,beta,d]",
          Ld.names == ["rpole","omega","inc","PA","beta","d"] ? 0.0 : 1.0)
    sp = default_star_params(2; ldtype = 0, rpole = 0.849, d = 16.8,
                             frac_escapevel = 0.92, rotation_period = 1/1.12)
    θd = layout_theta(Ld, sp)
    check("layout_theta/layout_merge round trip",
          layout_theta(Ld, layout_merge(Ld, sp, θd)) == θd ? 0.0 : 1.0)
end

let
    θl = collect(range(1e-3, π - 1e-3, length = 41))
    sl, cl = sin.(θl), cos.(θl)
    wl = [sin(3k) for k in 1:41]           # deterministic, non-constant weights
    fl(p) = dot(wl, logg_map(p[1], p[2], p[3], p[4], sl, cl))
    p0 = [0.849, 16.8, 0.92, 1/1.12]
    ana = FiniteDifferences.grad(FDM, fl, p0)[1]
    if zygote_ran[]
        gz = Zygote.gradient(fl, p0)[1]
        check("logg_map rrule vs FD (Zygote must not hit the mutation error)",
              relerr(gz, ana); tol = 1e-8)
    end
end

if zygote_ran[]
    let gpath = joinpath(@__DIR__, "..", "demos", "data", "betcas_H_korg.fits")
        if isfile(gpath)
            prov = TabulatedProvider(load_intensity_grid(gpath))
            bp0 = merge(base, (ldtype = 0, d = 16.8, rotation_period = 1/1.12,
                               tpole = 7208.0))
            Ld = parametric_layout(ldtype = 0, distance_free = true)
            θd = [0.849, 0.92, 19.9, -7.09, 0.25, 16.8]
            bnd = band_of(data[1])
            lpT = build_parametric_logπ(data, tessels, tepochs, bp0;
                                        provider = prov, layout = Ld, band = bnd)
            gT = Zygote.gradient(lpT, θd)[1]
            gF = FiniteDifferences.grad(FDM, lpT, θd)[1]
            check("provider path: full 6-param ∇logπ vs FD", relerr(gT, gF); tol = 1e-6)
            # The two null-space gates.
            lpP = build_parametric_logπ(data, tessels, tepochs, bp0;
                                        provider = PlanckProvider(), layout = Ld, band = bnd)
            gP = Zygote.gradient(lpP, θd)[1]
            check("∂logπ/∂d is EXACTLY zero under Planck (d enters only via logg)",
                  abs(gP[6]) / (norm(gP) + eps()))
            check("∂logπ/∂d is NON-zero under the atmosphere grid",
                  abs(gT[6]) / (norm(gT) + eps()) > 1e-8 ? 0.0 : 1.0)
            # Counting limb darkening twice must be refused, not quietly done.
            check("provider + ldtype≠0 is refused",
                  (try build_parametric_logπ(data, tessels, tepochs, base;
                                             provider = prov, layout = Ld, band = bnd)
                       1.0
                   catch; 0.0; end))
            check("provider without `band` is refused",
                  (try build_parametric_logπ(data, tessels, tepochs, bp0;
                                             provider = prov, layout = Ld)
                       1.0
                   catch; 0.0; end))
        else
            @warn "betcas_H_korg.fits absent — the provider-gradient section was skipped"
        end
    end
end

# ===========================================================================================
# [P5b] every sampler reaches the new parameters through ONE spec
# ===========================================================================================
# `_fit_hmc`, `_fit_pigeons` and `fit_parametric_ultranest` each carried their own copy of
# "build the rapid rotator's logπ, its bounds and its free set". `parametric_posterior_spec`
# is now that copy, once, in the core — which is also where `_box_transform` lives, and for
# the same reason: sibling package extensions cannot import from one another.
println("\n[P5b] the shared posterior spec, and the samplers' kwargs")

let gpath = joinpath(@__DIR__, "..", "demos", "data", "betcas_H_korg.fits")
    prov = isfile(gpath) ? TabulatedProvider(load_intensity_grid(gpath)) : PlanckProvider()
    bp0  = merge(base, (ldtype = 0, d = 16.8, rotation_period = 1/1.12, tpole = 7208.0))
    Ld   = parametric_layout(ldtype = 0, distance_free = true)
    θd   = [1.37, 0.6, 78.0, 24.0, 0.15, 16.8]
    bnd  = band_of(data[1])
    lp, lo, hi, idx, L2 = parametric_posterior_spec(data, tessels, tepochs, bp0;
                              free = ["rpole", "beta"], provider = prov,
                              layout = Ld, band = bnd)
    check("spec returns the layout it was given", L2 === Ld ? 0.0 : 1.0)
    check("spec's bounds are the layout's",
          (lo == Ld.lower && hi == Ld.upper) ? 0.0 : 1.0)
    check("spec resolves `free` by name", idx == [1, 5] ? 0.0 : 1.0)
    # The spec must be the SAME posterior `build_parametric_logπ` gives — if the two drift, a
    # sampler and a gradient fit silently optimise different models.
    direct = build_parametric_logπ(data, tessels, tepochs, bp0;
                                   provider = prov, layout = Ld, band = bnd)
    check("spec's logπ == build_parametric_logπ's",
          abs(lp(θd) - direct(θd)) / max(abs(lp(θd)), 1.0))
end

# Every entry point must ACCEPT provider/layout. `bootstrap_parametric` in particular took
# them and silently dropped them at one point, which would have bootstrapped the wrong model.
let want = (:provider, :layout)
    for (nm, fn) in (("fit_parametric", fit_parametric),
                     ("bootstrap_parametric", bootstrap_parametric))
        ms = methods(fn)
        ok = !isempty(ms) && all(w -> w in Base.kwarg_decl(first(ms)), want)
        # Zero methods means the Zygote extension is not loaded; that is a skip, not a failure.
        isempty(ms) ? (@warn "$nm has no methods (extension not loaded) — kwarg check skipped") :
                      check("$nm accepts provider and layout", ok ? 0.0 : 1.0)
    end
end

# The GUI derives its name->θ table from the layout instead of restating it. Assert the
# derivation reproduces the table it replaced, since the GUI itself cannot be loaded here.
let L = parametric_layout(; tpole_free = true)
    derived = Dict(L.fields[i] => L.names[i] for i in eachindex(L.names))
    check("GUI's PARAMETRIC_THETA derivation is unchanged",
          derived == Dict(:rpole=>"rpole", :frac_escapevel=>"omega", :inclination=>"inc",
                          :position_angle=>"PA", :beta=>"beta", :ld1=>"ld1", :ld2=>"ld2",
                          :tpole=>"tpole") ? 0.0 : 1.0)
    check("GUI's θfull order matches parametric_param_names",
          [derived[f] for f in L.fields] == parametric_param_names(; tpole_free=true) ?
          0.0 : 1.0)
end

# Printed LAST so it counts the [P5] section appended above; runtests.jl asserts nfail[].
@printf("\n=== %d passed, %d failed ===\n", npass[], nfail[])
