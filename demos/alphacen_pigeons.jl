#!/usr/bin/env julia
# Does Pigeons actually work? — parallel tempering on a problem with a KNOWN answer.
#
#     julia -t auto --project=bin demos/alphacen_pigeons.jl
#
# WHY THIS EXISTS. Every Pigeons run in this repository so far has FAILED. On β Cas it
# finished `round_trips = 0` with a global barrier Λ ≈ 7.95 against only 10 chains, the χ² at
# its posterior median (6779, χ²ᵣ = 5.38) was WORSE than VMLMB's in the same basin, and the β
# marginal sat at 0.397 where the local optimum there is 0.23 — so the sample was discarded.
# That is evidence that the SETUP was wrong, and no evidence at all about the sampler. Until
# Pigeons has been seen to get a right answer, "Pigeons ran and produced a posterior" is not
# something any test should assert.
#
# THE PROBLEM. α Cen A, one night of VLTI H-band data (324 V², 156 T3amp, 156 T3φ, baselines
# 10–128 m), as a limb-darkened sphere: `sphere_param_names(1) == ["radius", "ld1"]`. Exactly
# radius plus one limb-darkening coefficient.
#
# WHY TWO PARAMETERS IS THE POINT. The posterior is then a function on a plane, so it can be
# computed EXACTLY by quadrature on a grid — marginals, quantiles and the log-evidence. The
# comparison is therefore against truth, not against another sampler that might be wrong in
# the same way. Agreeing with Nautilus would not establish that either of them is right.
#
# And there is an external anchor as well: Kervella et al. (2017) give α Cen A
# θ_LD = 8.502 ± 0.038 mas, i.e. radius = 4.251 mas.
#
# WHAT WOULD COUNT AS WORKING. Three things, and the third is the one β Cas failed:
#   1. the marginals match the quadrature truth,
#   2. `stepping_stone` log Z matches the quadrature log-evidence,
#   3. `round_trips` is comfortably above zero — the ladder was actually traversed.
# A sampler can pass 1 and 2 on a unimodal posterior while failing 3, which is precisely how
# a multimodal problem then goes wrong silently.

using ROTIR, Statistics, Printf, LinearAlgebra
# ALL FIVE, because `ROTIRPigeonsExt` is triggered by the whole list — Project.toml says
# `ROTIRPigeonsExt = ["Pigeons", "Distributions", "LogDensityProblems", "ADTypes", "Zygote"]`.
# With `using Pigeons` alone `methods(ROTIR._fit_pigeons)` is EMPTY and
# `ROTIR.pigeons_available()` is false, which is how the first attempt got a MethodError
# listing argument types for a function that had no methods at all.
using Pigeons, Distributions, LogDensityProblems, ADTypes, Zygote

const OUT = String[]
# PRINTED AS IT GOES, and kept for the recap at the end. The buffer alone is not enough: the
# first run of this script died at the Pigeons call AFTER the quadrature had finished, and
# because nothing had been printed yet the 3 minutes of grid work went in the bin with the
# stack trace.
say(s) = (push!(OUT, s); println(s); flush(stdout); nothing)
# `@sprintf` is a macro and needs its format string at parse time, so it cannot be splatted
# into — `Printf.format` is the runtime equivalent and is what lets this be one helper rather
# than an `@sprintf` at every call site.
sayf(fmt, args...) = say(Printf.format(Printf.Format(fmt), args...))

# ── the model ───────────────────────────────────────────────────────────────────────────
const NSIDE_EXP = 4          # 3072 tessels, 6.1 ms/logπ — measured; nside 3 biases the
                             # diameter low by ~1% because the mesh under-resolves the limb
const FILE = joinpath(pkgdir(ROTIR), "demos", "data", "AlphaCenA.oifits")

# ── WHICH OBSERVABLES, AND WHY NOT T3φ ──────────────────────────────────────────────────
# `weights = (V², T3amp, T3φ)`. T3φ is OFF, matching `fit_sphere_ld`'s own default and for the
# reason its docstring gives with measured numbers: a limb-darkened sphere is CENTROSYMMETRIC,
# so its visibilities are real and its closure phases are identically 0° or 180°. The measured
# T3φ therefore cannot be improved by any (radius, ld1) — the model has no freedom to move them
# — while every calibration residual in them lands in the χ², making the likelihood
# overconfident about the parameters it CAN move. On RW Cep, adding T3φ moved the fitted
# diameter 2.93 → 2.56 mas and took χ²v2/n from 8.9 to 24.8.
#
# T3amp is KEPT: it is |V₁V₂V₃| and does carry diameter information, and the same RW Cep table
# shows V² + T3amp is the best of the three combinations.
#
# There is a numerical reason too, found while adding this: near zero closure phase `angle(t3)`
# is ill-conditioned, and the finite-difference check of ∇logπ degrades from 3e-9 (V² only) to
# 8e-6 once T3φ is in — on the unchanged code path. The analytic gradient is exactly linear in
# these weights to 3e-15, so it is the finite difference that suffers, but it is a sign that
# T3φ on a symmetric target is poorly behaved as well as uninformative.
const WEIGHTS = (1, 1, 0)

data = readoifits(FILE)[1]
tess = tessellation_healpix(NSIDE_EXP; T = Float64)      # Float64: the reference must not
                                                         # carry mesh precision into the answer
base = default_star_params(0; radius = 4.25, ldtype = 1, ld1 = 0.3)
logπ = build_sphere_logπ([data], tess, [0.0], base; weights = WEIGHTS)
const NAMES = sphere_param_names(1)
const LB, UB = default_sphere_bounds(1)

say("α Cen A, limb-darkened sphere: θ = $(NAMES)")
const NDATA = WEIGHTS[1]*length(data.v2) + (WEIGHTS[2]>0)*length(data.t3amp) +
              (WEIGHTS[3]>0)*length(data.t3phi)
sayf("data: %d V², %d T3amp, %d T3φ   weights = %s  → %d points fitted",
     length(data.v2), length(data.t3amp), length(data.t3phi), string(WEIGHTS), NDATA)
sayf("mesh: HEALPix nside_exp = %d (%d tessels), Float64", NSIDE_EXP, 12 * (2^NSIDE_EXP)^2)

# ── 1. the peak, and the curvature that sets the grid box ───────────────────────────────
# Coarse-to-fine rather than an optimiser: this has to be reproducible and cannot be allowed
# to land in a line-search failure, which `fit_parametric`'s own docstring warns about.
function refine_peak(logπ; r0 = 4.25, l0 = 0.2)
    r, l = r0, l0
    for (dr, dl) in ((0.25, 0.5), (0.05, 0.1), (0.01, 0.02), (0.002, 0.004))
        best = (-Inf, r, l)
        for rr in (r - 5dr):dr:(r + 5dr), ll in (l - 5dl):dl:(l + 5dl)
            (ll < LB[2] || ll > UB[2]) && continue
            v = logπ([rr, ll])
            isfinite(v) && v > best[1] && (best = (v, rr, ll))
        end
        _, r, l = best
    end
    return (r, l)
end

"""
    profile_halfwidth(logπ, θ, k; drop) -> h

The half-width in coordinate `k` at which logπ has fallen `drop` nats below the peak.

NOT A HESSIAN, and that is the point. The first version of this script used a finite-difference
Hessian and got σ(radius) = 9e-5 mas when the truth is ~50x larger — so the ±8σ grid box came
out flat (density at the edge equal to the peak) and three minutes of quadrature were void.

The reason is structural rather than a bad step choice: a TESSELLATED likelihood is only
piecewise smooth, because moving the radius slides the limb across a 3072-tessel mesh in
discrete jumps. A four-point second difference divides by 4h², amplifying that mesh noise by
~6e4 at h = 2e-3, which is the same order as the true curvature 1/σ² — so the Hessian is noise.
Doubling outward needs no such luck: it only requires logπ to decrease on average, and it
reports the width at a STATED contour instead of a curvature. `drop = 30` nats is ~7.7σ for a
Gaussian, so ±h brackets essentially all the mass without assuming the shape is Gaussian.
"""
function profile_halfwidth(logπ, θ, k; drop = 30.0, h0 = 1e-5, hmax = 1.0, lo = -Inf, hi = Inf)
    f0 = logπ(θ)
    h = h0
    while h < hmax
        a = copy(θ); a[k] = min(θ[k] + h, hi)
        b = copy(θ); b[k] = max(θ[k] - h, lo)
        va, vb = logπ(a), logπ(b)
        # The WORSE side: a coordinate against a bound stops falling on one side, and taking the
        # max there would run h to hmax and blow the box up.
        w = min(isfinite(va) ? va : -Inf, isfinite(vb) ? vb : -Inf)
        w < f0 - drop && return h
        h *= 1.6
    end
    return hmax
end

θhat = refine_peak(logπ)
σ = [profile_halfwidth(logπ, collect(θhat), 1; hmax = 0.5),
     profile_halfwidth(logπ, collect(θhat), 2; lo = LB[2], hi = UB[2], hmax = 1.0)] ./ 7.7
sayf("peak:  radius = %.5f mas  ld1 = %.5f   (diameter = %.4f mas)",
     θhat[1], θhat[2], 2 * θhat[1])
sayf("logπ at the peak = %.4f, i.e. χ² = %.2f, χ²ᵣ = %.3f",
     logπ(collect(θhat)), -2 * logπ(collect(θhat)),
     -2 * logπ(collect(θhat)) / (NDATA - 2))
sayf("profile σ (width at Δlogπ = 30, /7.7): radius %.6f, ld1 %.6f", σ[1], σ[2])

# ── 2. THE EXACT POSTERIOR, by quadrature ───────────────────────────────────────────────
# ±8σ in each direction on a 241 x 241 grid. Threaded over rows: the docstring for
# `_fit_pigeons` states the likelihood builds its own geometry per call and shares nothing
# mutable, which is what makes this safe — and it is the same property Pigeons relies on.
const NG = 201
const SPAN = 8.0

function grid_posterior(logπ, θhat, σ; ng = NG, span = SPAN)
    rs = range(θhat[1] - span * σ[1], θhat[1] + span * σ[1]; length = ng)
    ls = range(max(θhat[2] - span * σ[2], LB[2]),
               min(θhat[2] + span * σ[2], UB[2]); length = ng)
    L = fill(-Inf, ng, ng)
    Threads.@threads for i in 1:ng
        for j in 1:ng
            v = logπ([rs[i], ls[j]])
            L[i, j] = isfinite(v) ? v : -Inf
        end
    end
    return rs, ls, L
end

"""
Marginal quantiles and the log-evidence from a log-posterior on a regular grid.

Trapezoid in both directions, in a shifted exponential so the sum never overflows: the peak
logπ here is ~-2000 and `exp` of that is zero in Float64, so normalising BEFORE exponentiating
is not an optimisation but the only way the integral exists at all.
"""
function grid_summary(rs, ls, L)
    m = maximum(L)
    P = exp.(L .- m)                                   # unnormalised, peak 1
    dr = step(rs); dl = step(ls)
    w(n) = (v = fill(1.0, n); v[1] = v[end] = 0.5; v)   # trapezoid weights
    wr, wl = w(length(rs)), w(length(ls))
    Z = sum(P .* (wr * wl')) * dr * dl                  # ∫∫ exp(logπ - m)
    logZ = m + log(Z)
    # Marginals, each normalised to 1.
    pr = vec(sum(P .* (ones(length(rs)) * wl'), dims = 2)) .* dl
    pl = vec(sum(P .* (wr * ones(length(ls))'), dims = 1)) .* dr
    pr ./= sum(pr .* wr) * dr
    pl ./= sum(pl .* wl) * dl
    function quantiles(x, p)
        c = cumsum(p .* step(x)); c ./= c[end]
        [x[findfirst(>=(q), c)] for q in (0.16, 0.5, 0.84)]
    end
    return (logZ = logZ,
            r = quantiles(rs, pr), l = quantiles(ls, pl),
            mass_at_edge = max(maximum(P[1, :]), maximum(P[end, :]),
                               maximum(P[:, 1]), maximum(P[:, end])))
end

t_grid = @elapsed ((rs, ls, L) = grid_posterior(logπ, θhat, σ))
G = grid_summary(rs, ls, L)
say("")
sayf("── EXACT posterior by quadrature (%d x %d grid, ±%.0fσ, %d threads, %.1f s) ──",
     NG, NG, SPAN, Threads.nthreads(), t_grid)
sayf("radius = %.5f  [%.5f, %.5f]  (q16, median, q84)", G.r[2], G.r[1], G.r[3])
sayf("ld1    = %.5f  [%.5f, %.5f]", G.l[2], G.l[1], G.l[3])
sayf("log Z  = %.4f", G.logZ)
# A box that clips the posterior invalidates the quadrature, so say whether it did.
sayf("largest density on the box edge, relative to the peak: %.2e %s",
     G.mass_at_edge, G.mass_at_edge < 1e-6 ? "(negligible — the box contains the mass)" :
                                             "(!! WIDEN THE BOX)")

# ── 3. Pigeons: a PILOT first, to measure the barrier ───────────────────────────────────
# The β Cas lesson, applied: size `n_chains` against the global barrier Λ rather than
# guessing. The guidance is n_chains ≳ 2Λ, and β Cas was run at 10 chains against Λ ≈ 7.95.
"Λ, however this Pigeons version exposes it."
function barrier(pt)
    for f in (() -> Pigeons.global_barrier(pt),
              () -> Pigeons.global_barrier(pt.shared.tempering),
              () -> pt.shared.tempering.communication_barriers.globals[end])
        try
            v = f()
            v isa Real && isfinite(v) && return Float64(v)
        catch
        end
    end
    return NaN
end

# `ROTIR._fit_pigeons`, not `_fit_pigeons`: the samplers are internal hooks filled in by
# weak-dependency extensions (`ROTIR.pigeons_available()` is `!isempty(methods(_fit_pigeons))`),
# so they are reachable but deliberately unexported.
run_pigeons(; n_chains, n_rounds) = ROTIR._fit_pigeons([data], tess, [0.0], base;
    θ0 = collect(θhat), free = NAMES, model = :sphere,
    n_chains = n_chains, n_rounds = n_rounds, explorer = :slice,
    multithreaded = true, seed = 20261001, verb = true, weights = WEIGHTS)

t_pilot = @elapsed pilot = run_pigeons(n_chains = 10, n_rounds = 5)
Λ = barrier(pilot.result)
say("")
sayf("── Pigeons PILOT: 10 chains, 2^5 scans, %.1f s ──", t_pilot)
sayf("global barrier Λ = %.3f  → guidance wants n_chains ≳ %d", Λ, ceil(Int, 2Λ))
sayf("round trips = %d", pilot.round_trips)

nchains = isfinite(Λ) ? max(10, ceil(Int, 2Λ) + 2) : 20
# 2^8, not 2^9: each round doubles the work, so the last one costs as much as all the
# previous put together. 256 scans is ample for a two-parameter unimodal posterior, and the
# round-trip count below is what says whether it was.
t_pig = @elapsed pig = run_pigeons(n_chains = nchains, n_rounds = 8)
say("")
sayf("── Pigeons SIZED BY Λ: %d chains, 2^8 scans, %.1f s ──", nchains, t_pig)
sayf("round trips = %d   (0 means the ladder was never traversed)", pig.round_trips)
sayf("samples = %d", size(pig.samples, 1))
sayf("radius = %.5f  [%.5f, %.5f]", pig.median[1], pig.q16[1], pig.q84[1])
sayf("ld1    = %.5f  [%.5f, %.5f]", pig.median[2], pig.q16[2], pig.q84[2])
sayf("log Z  = %.4f   (stepping stone)", pig.logz)

# ── 4. the verdict, against truth and against the literature ────────────────────────────
# In units of the TRUE posterior width, which is the only scale on which "close" means
# anything: a 0.001 mas discrepancy is nothing if σ is 0.01 and fatal if σ is 0.0001.
σr_true = (G.r[3] - G.r[1]) / 2
σl_true = (G.l[3] - G.l[1]) / 2
dr = (pig.median[1] - G.r[2]) / σr_true
dl = (pig.median[2] - G.l[2]) / σl_true
wr = ((pig.q84[1] - pig.q16[1]) / 2) / σr_true
wl = ((pig.q84[2] - pig.q16[2]) / 2) / σl_true

say("")
say("──────────────────────── VERDICT ────────────────────────")
sayf("%-22s %12s %12s %10s", "", "quadrature", "Pigeons", "Δ / σ_true")
sayf("%-22s %12.5f %12.5f %10.2f", "radius [mas]", G.r[2], pig.median[1], dr)
sayf("%-22s %12.5f %12.5f %10.2f", "ld1", G.l[2], pig.median[2], dl)
sayf("%-22s %12.5f %12.5f %10.2f", "σ(radius)", σr_true, (pig.q84[1]-pig.q16[1])/2, wr)
sayf("%-22s %12.5f %12.5f %10.2f", "σ(ld1)", σl_true, (pig.q84[2]-pig.q16[2])/2, wl)
sayf("%-22s %12.4f %12.4f %10.4f", "log Z", G.logZ, pig.logz, pig.logz - G.logZ)
say("")
# Kervella et al. 2017, θ_LD = 8.502 ± 0.038 mas. NOT a pass/fail criterion for the sampler —
# the model here is one night, one band, a uniform-ld1 sphere and a 3072-tessel mesh — but a
# diameter far from it would mean the LIKELIHOOD is wrong, in which case agreeing with the
# quadrature would only prove the two agree about the wrong thing.
sayf("diameter = %.4f mas against Kervella+2017 θ_LD = 8.502 ± 0.038  (%.1fσ_lit)",
     2 * pig.median[1], (2 * pig.median[1] - 8.502) / 0.038)
say("")
ok_r = abs(dr) < 0.25
ok_l = abs(dl) < 0.25
ok_w = 0.75 < wr < 1.33 && 0.75 < wl < 1.33
ok_z = abs(pig.logz - G.logZ) < 1.0
ok_t = pig.round_trips >= 3
for (lbl, ok) in (("marginal medians match truth (<0.25σ)", ok_r && ok_l),
                  ("marginal WIDTHS match truth (within 33%)", ok_w),
                  ("log Z matches the quadrature (<1 nat)", ok_z),
                  ("the ladder was traversed (round trips ≥ 3)", ok_t))
    sayf("  %-45s %s", lbl, ok ? "PASS" : "FAIL")
end
say("")
say(all((ok_r, ok_l, ok_w, ok_z, ok_t)) ?
    "PIGEONS WORKS on this problem. The β Cas failure was the setup, not the sampler." :
    "PIGEONS DID NOT REPRODUCE THE KNOWN ANSWER — see which line failed above.")

println("\n" * "="^76)
for l in OUT; println(l); end
println("="^76)
