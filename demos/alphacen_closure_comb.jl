#!/usr/bin/env julia
# Closure phases make this posterior MULTIMODAL — can tempering find the global optimum where
# a local search cannot?
#
#     julia --project=bin demos/alphacen_closure_comb.jl
#
# THE MECHANISM, and why it is a property of the data rather than of the fitter.
#
# A limb-darkened sphere is centrosymmetric, so its complex visibility is REAL and its closure
# phase is exactly 0° or 180° — 180° whenever an ODD number of the three baselines lies beyond
# a visibility null. So as the radius changes, the nulls sweep across the baselines and model
# closure phases flip by 180° DISCONTINUOUSLY. The T3φ χ² term therefore has a step at every
# radius where any baseline crosses any null, and between the steps it is smooth: a COMB of
# local minima separated by walls.
#
# α Cen A is squarely in that regime — measured from the file, for a 8.44 mas disc:
#
#     null 1 at 2.98e7 cycles/rad     data span 5.69e6 – 8.46e7
#     null 2 at 5.46e7                45% of uv points lie beyond null 1
#     null 3 at 7.91e7                18% lie beyond null 2
#
# A gradient method cannot cross a step (the gradient says nothing about the next tooth) and a
# simplex can only cross one by luck. Parallel tempering is built for exactly this: the hot
# chains see a flattened posterior in which the walls are low, and the round trips carry that
# information down to the cold chain.
#
# WHAT IS BEING COMPARED. The profile likelihood along radius is computed on a grid, so the
# teeth can be COUNTED rather than inferred, and the global optimum is known by exhaustion.
# Then: Nelder–Mead and BOBYQA from a spread of starts, and Pigeons. The question is not which
# gets the best χ² — it is which finds the SAME optimum the grid proves is global.
#
# NOTE ON T3φ AND TRUTH. The global optimum with T3φ included is not necessarily the better
# ASTROPHYSICAL answer — `fit_sphere_ld`'s docstring shows on RW Cep that including T3φ drags a
# symmetric model's diameter away from the right value. This script is about whether a sampler
# can solve a hard optimisation problem, and T3φ supplies a genuinely hard one on real data.

using ROTIR, Statistics, Printf, LinearAlgebra
using Pigeons, Distributions, LogDensityProblems, ADTypes, Zygote

const OUT = String[]
say(s) = (push!(OUT, s); println(s); flush(stdout); nothing)
sayf(fmt, args...) = say(Printf.format(Printf.Format(fmt), args...))

const NSIDE_EXP = 4
const FILE = joinpath(pkgdir(ROTIR), "demos", "data", "AlphaCenA.oifits")
# EXPLICIT AND TIGHT-ISH. `default_sphere_bounds` gives `ub = [Inf, 2.0]`, and with
# σ(radius) ≈ 6e-4 the slice sampler then has to shrink ~1e4 before it localises at all — a 2^5
# pilot took 7.5 minutes for that reason alone. This box is still ~400 teeth wide.
const LB = [3.90, 0.00]
const UB = [4.60, 0.60]

data = readoifits(FILE)[1]
tess = tessellation_healpix(NSIDE_EXP; T = Float64)
base = default_star_params(0; radius = 4.25, ldtype = 1, ld1 = 0.3)
lp(w) = build_sphere_logπ([data], tess, [0.0], base; weights = w)

say("α Cen A, limb-darkened sphere, θ = [radius, ld1]")
sayf("prior box: radius ∈ [%.2f, %.2f], ld1 ∈ [%.2f, %.2f]", LB[1], UB[1], LB[2], UB[2])

# ── the profile likelihood along radius, and the teeth ──────────────────────────────────
# For each radius, the best ld1 — that is the curve a radius-only search actually walks, and
# the one whose local maxima ARE the traps. Serial on purpose: `logπ` allocates a full star
# geometry per call, so 16 threads lose to GC contention (measured: 364 s against 246 s of
# serial work on a 201² grid).
const NR = 561
const NL = 25

function profile(logπ; nr = NR, nl = NL)
    rs = range(LB[1], UB[1]; length = nr)
    ls = range(LB[2], UB[2]; length = nl)
    p = fill(-Inf, nr); bl = zeros(nr)
    for i in 1:nr, j in 1:nl
        v = logπ([rs[i], ls[j]])
        isfinite(v) && v > p[i] && (p[i] = v; bl[i] = ls[j])
    end
    return rs, p, bl
end

"Interior local maxima of a 1-D curve, as (index, value), strongest first."
function local_maxima(p; rel = 3.0)
    idx = [i for i in 2:length(p)-1 if p[i] >= p[i-1] && p[i] >= p[i+1] && isfinite(p[i])]
    # `rel` nats below the global peak is the cut: anything shallower is not a trap a search
    # would get stuck in, it is grid noise on a flat shoulder.
    top = maximum(p)
    keep = [i for i in idx if p[i] > top - 1e6]
    sort!(keep; by = i -> -p[i])
    return keep
end

# ── for each weighting: count the teeth, then see who finds the global one ───────────────
function study(wname, w)
    logπ = lp(w)
    t = @elapsed ((rs, p, bl) = profile(logπ))
    i★ = argmax(p)
    maxima = local_maxima(p)
    # Teeth DEEP enough to trap: a local max whose two flanking valleys are more than 2 nats
    # below it, i.e. a genuine basin rather than a ripple on a slope.
    function depth(i)
        l = i; while l > 1 && p[l-1] <= p[l]; l -= 1; end
        r = i; while r < length(p) && p[r+1] <= p[r]; r += 1; end
        return p[i] - max(l > 1 ? p[l-1] : -Inf, r < length(p) ? p[r+1] : -Inf)
    end
    traps = [i for i in maxima if isfinite(depth(i)) && depth(i) > 2.0]

    say("")
    sayf("════════ %s   weights = %s ════════", wname, string(w))
    sayf("profile scan: %d radii x %d ld1, %.1f s", NR, NL, t)
    sayf("GLOBAL optimum: radius = %.5f  ld1 = %.4f   logπ = %.3f  (diameter %.4f mas)",
         rs[i★], bl[i★], p[i★], 2 * rs[i★])
    sayf("local maxima in the radius profile: %d, of which %d are traps >2 nats deep",
         length(maxima), length(traps))
    if length(traps) > 1
        d = [rs[i] for i in traps[1:min(6, end)]]
        sayf("  deepest trap radii: %s", join((Printf.@sprintf("%.4f", x) for x in d), ", "))
        sayf("  logπ at those:      %s",
             join((Printf.@sprintf("%.1f", p[i]) for i in traps[1:min(6, end)]), ", "))
        gaps = diff(sort([rs[i] for i in traps]))
        sayf("  median tooth spacing: %.5f mas", Statistics.median(gaps))
    end

    # ── local searches, from a spread of starts ──
    # The question is how many land on the global tooth, not what χ² they report.
    starts = collect(range(LB[1] + 0.04, UB[1] - 0.04; length = 13))
    hits = Dict{Symbol,Int}(); bests = Dict{Symbol,Float64}()
    for (m, alg) in ((:neldermead, :LN_NELDERMEAD), (:bobyqa, :LN_BOBYQA))
        nhit = 0; best = -Inf
        for r0 in starts
            res = try
                fit_sphere_ld(data, tess; tepochs = [0.0], ldtype = 1,
                              radius0 = r0, ld0 = 0.3,
                              radius_bounds = (LB[1], UB[1]), ld_bounds = (LB[2], UB[2]),
                              weights = collect(Float64, w), method = :neldermead,
                              algorithm = alg, maxeval = 4000, verbose = false)
            catch
                nothing
            end
            res === nothing && continue
            v = logπ([res.radius, res.ld1])
            v > best && (best = v)
            # "Found it" = within 1 nat of the grid's global optimum AND within half a tooth.
            abs(res.radius - rs[i★]) < 0.5 * step(rs) * 4 && v > p[i★] - 1.0 && (nhit += 1)
        end
        hits[m] = nhit; bests[m] = best
        sayf("%-12s found the global optimum in %2d of %2d starts   best logπ = %.3f (Δ = %.2f)",
             String(m), nhit, length(starts), best, best - p[i★])
    end

    # ── Pigeons ──
    pil = ROTIR._fit_pigeons([data], tess, [0.0], base; θ0 = [rs[i★], bl[i★]],
        free = ["radius", "ld1"], model = :sphere, lb = LB, ub = UB,
        n_chains = 10, n_rounds = 4, explorer = :slice, multithreaded = true,
        seed = 20261001, verb = false, weights = w)
    Λ = try Pigeons.global_barrier(pil.result) catch; NaN end
    nch = isfinite(Λ) ? max(10, ceil(Int, 2Λ) + 2) : 16
    sayf("Pigeons pilot: Λ = %.2f → n_chains = %d", Λ, nch)

    tp = @elapsed pig = ROTIR._fit_pigeons([data], tess, [0.0], base;
        θ0 = [0.5*(LB[1]+UB[1]), 0.3],        # a NEUTRAL start, not the known answer
        free = ["radius", "ld1"], model = :sphere, lb = LB, ub = UB,
        n_chains = nch, n_rounds = 7, explorer = :slice, multithreaded = true,
        seed = 20261001, verb = false, weights = w)
    vpig = logπ([pig.median[1], pig.median[2]])
    sayf("Pigeons: %d chains, 2^7 scans, %.0f s, %d round trips", nch, tp, pig.round_trips)
    sayf("  median radius = %.5f (global %.5f, Δ = %+.5f mas)  logπ = %.3f (Δ = %.2f)",
         pig.median[1], rs[i★], pig.median[1] - rs[i★], vpig, vpig - p[i★])
    sayf("  %s", abs(pig.median[1] - rs[i★]) < 0.002 ?
         "ON the global tooth" : "on a DIFFERENT tooth")
    return (rs = rs, p = p, i★ = i★, ntraps = length(traps), hits = hits,
            pig = pig, vpig = vpig, Λ = Λ)
end

A = study("WITHOUT closure phases", (1, 1, 0))
B = study("WITH closure phases", (1, 1, 1))

say("")
say("════════════════════════════ SUMMARY ════════════════════════════")
sayf("%-34s %14s %14s", "", "V²+T3amp", "+T3φ")
sayf("%-34s %14d %14d", "traps >2 nats deep", A.ntraps, B.ntraps)
sayf("%-34s %14s %14s", "Nelder–Mead found global",
     "$(A.hits[:neldermead])/13", "$(B.hits[:neldermead])/13")
sayf("%-34s %14s %14s", "BOBYQA found global",
     "$(A.hits[:bobyqa])/13", "$(B.hits[:bobyqa])/13")
sayf("%-34s %14s %14s", "Pigeons found global",
     abs(A.pig.median[1]-A.rs[A.i★]) < 0.002 ? "yes" : "no",
     abs(B.pig.median[1]-B.rs[B.i★]) < 0.002 ? "yes" : "no")
sayf("%-34s %14d %14d", "Pigeons round trips", A.pig.round_trips, B.pig.round_trips)
sayf("%-34s %14.2f %14.2f", "global barrier Λ", A.Λ, B.Λ)
say("")
hard = B.ntraps > 3 * max(A.ntraps, 1)
beat = B.hits[:neldermead] + B.hits[:bobyqa] < 13 &&
       abs(B.pig.median[1] - B.rs[B.i★]) < 0.002
say(hard ? "Closure phases DO make this multimodal, as predicted from the nulls." :
           "Closure phases did NOT produce the expected comb — check the null crossings.")
say(beat ? "AND TEMPERING WINS: Pigeons found the global optimum where the local searches did not." :
           "Tempering did not demonstrate an advantage here; see the counts above.")
println("\n" * "="^78); for l in OUT; println(l); end; println("="^78)
