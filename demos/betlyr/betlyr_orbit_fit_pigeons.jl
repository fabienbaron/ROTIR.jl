#!/usr/bin/env julia
# betlyr_orbit_fit_pigeons.jl — β Lyrae orbital elements by parallel tempering.
#
#   julia --project=demos -t auto betlyr/betlyr_orbit_fit_pigeons.jl
#   N_ROUNDS=12 N_CHAINS=16 julia --project=demos -t auto betlyr/betlyr_orbit_fit_pigeons.jl
#
# Same model as betlyr_orbit_fit_ultranest.jl (shared via betlyr_model.jl);
# different sampler.
#
# ---------------------------------------------------------------------------
# Why Pigeons rather than UltraNest here
# ---------------------------------------------------------------------------
# Two reasons, one of them a hard constraint:
#
#  1. PARALLELISM. UltraNest is reached through PyCall, and PyCall is not thread-safe:
#     `pydecref_` calls Py_DecRef with no GIL check, so when Julia's GC runs a PyObject
#     finalizer on a worker thread the process segfaults. That happens regardless of
#     whether the likelihood itself touches Python. So with UltraNest the likelihood is
#     stuck on one core. Pigeons is native Julia and `multithreaded = true` parallelises
#     over the chain ladder — every scan steps all `n_chains` at once, so the speedup is
#     both larger and steadier than UltraNest's uneven batches (measured there: mean 14.5
#     rows per call, median 5).
#
#  2. MULTIMODALITY DIAGNOSTICS. The point of tempering is that hot chains roam freely and
#     swaps carry the cold chain between basins; `n_round_trips` then *tells you* whether
#     that happened. Nested sampling gives no equivalent statement. For an orbit fit — where
#     the (Ω, T0) and inclination branches are genuine, well-separated modes — that matters.
#
# Pigeons also returns log(Z) via stepping stone, so nothing is given up by not using
# nested sampling.
#
# ---------------------------------------------------------------------------
# The PyCall hazard does not disappear just by changing sampler
# ---------------------------------------------------------------------------
# `using ROTIR` loads PyCall and PyPlot unconditionally (src/oiplot_spheroid.jl), so
# PyObjects exist in this session too. With `multithreaded = true` the log-potential runs
# on worker threads, and a PyObject finalizer firing there would crash exactly as before.
# Mitigation, in order of importance:
#   * NOTHING is plotted before sampling — PyPlot is imported only at the end.
#   * A full GC is forced before `pigeons(...)` so transient PyObjects created while
#     reading the OIFITS are finalized on the main thread, where it is safe.
#   * MULTI=0 falls back to single-threaded if a crash still appears.
# Set PLOT=0 to skip plotting entirely and keep the session Python-free after loading.

using Pigeons, Distributions, LogDensityProblems, LinearAlgebra, Random, Printf, Statistics
using DelimitedFiles
include(joinpath(@__DIR__, "betlyr_model.jl"))

const N_ROUNDS = parse(Int, get(ENV, "N_ROUNDS", "10"))     # 2^N_ROUNDS scans
const N_CHAINS = parse(Int, get(ENV, "N_CHAINS", "12"))
const MULTI    = get(ENV, "MULTI", "1") == "1"
const DOPLOT   = get(ENV, "PLOT", "1") == "1"
const OUTDIR   = get(ENV, "OUTDIR", joinpath(@__DIR__, "results"))
const EXPLORER = Symbol(get(ENV, "EXPLORER", "slice"))      # :slice | :mala
const VERBOSE  = get(ENV, "VERBOSE", "1") == "1"            # print at every round boundary
mkpath(OUTDIR)

const IDX  = free_indices()
const NDIM = length(IDX)
describe_model(IDX)
@printf("\nPigeons: %d chains, 2^%d = %d scans, multithreaded=%s, %d threads\n",
        N_CHAINS, N_ROUNDS, 2^N_ROUNDS, MULTI, Threads.nthreads())

# ── Target ──────────────────────────────────────────────────────────────────
# Sampled in the unconstrained z space (see betlyr_model.jl): θ = lo + (hi−lo)·σ(z), with
# the log-Jacobian included so a uniform prior on the box stays uniform after the change of
# variables. This also lets the reference be a plain wide Gaussian, which the tempering
# ladder can anneal from cleanly.
const POS = free_positions(IDX)

# Differentiable form: `z_to_theta_ad` builds θ by comprehension instead of writing into a
# `Vector{Float64}`, so Duals promote through it. Used for both the value and the gradient
# so the two can never drift apart.
logpost_z(z) = begin
    c = chi2_total(z_to_theta_ad(z, POS))
    c >= 1e12 && return -Inf
    v = -0.5c + log_jacobian(z, IDX)
    isfinite(v) ? v : -Inf
end

struct BetLyrTarget end
(::BetLyrTarget)(z) = logpost_z(z)
LogDensityProblems.logdensity(::BetLyrTarget, z) = logpost_z(z)
LogDensityProblems.dimension(::BetLyrTarget) = NDIM

# ---------------------------------------------------------------------------
# Gradients (EXPLORER=mala)
# ---------------------------------------------------------------------------
# EXPLORER=mala uses a TRUE ANALYTIC gradient (betlyr_gradient.jl) — no ForwardDiff, no
# Zygote in the sampling loop. Cost of value+gradient, 10 free parameters:
#     analytic          2.7x the primal      <- used here
#     ForwardDiff       8.3x
#     Zygote           11.3x
#     central diffs    20x    (2n evaluations)
#     central_fdm(5,1) 193x   (adaptive step search; a correctness oracle, not a gradient)
#
# MEASURED, 8 chains x 6 rounds, all three on identical settings:
#
#     SliceSampler                 77 s   0 round trips   log Z = -173750.7
#     AutoMALA + ANALYTIC gradient 146 s  0 round trips   log Z = -171077.0   (174790 grad calls)
#     AutoMALA + ForwardDiff       528 s  0 round trips   log Z = -171077.0
#
# The analytic gradient makes AutoMALA 3.6x faster and closes the gap with slice from 6.9x to
# 1.9x — but slice still wins, and NEITHER completes a round trip. The posterior is very sharp
# (χ² ≈ 2e6) and spans ~7 orders of magnitude in |∂χ²/∂θ| (2e4 for inclination, 2e11 for Ṗ);
# slice's per-coordinate adaptive width copes, a diagonally-preconditioned MALA does not.
# Keep EXPLORER=slice unless the posterior becomes better conditioned.
#
# ---------------------------------------------------------------------------
# HOW TO GIVE PIGEONS A HAND-WRITTEN GRADIENT
# ---------------------------------------------------------------------------
# Declaring `capabilities = LogDensityOrder{1}()` and defining `logdensity_and_gradient` is
# NOT enough — Pigeons routes the target through `LogDensityProblemsAD.ADgradient` and
# ForwardDiffs `logdensity` regardless. Verified with a call counter: zero calls.
#
# The hook is the THREE-argument form, reached from
#     get_buffer(ad_buffers, :target, kind, path.target, replica)   (BufferedAD.jl)
# and BOTH arguments must be typed. Pigeons' own method is
#     ADgradient(kind::AbstractADType, log_potential, replica::Replica)
# so typing only the target leaves the call AMBIGUOUS and the dispatch is silently lost —
# that failure mode cost three attempts. The returned `BufferedAD` forwards
# `logdensity_and_gradient` straight to ours.
if EXPLORER === :mala
    @eval using ForwardDiff, ADTypes          # ForwardDiff only for test_gradient()
    include(joinpath(@__DIR__, "betlyr_gradient.jl"))
    const LDPAD = Pigeons.LogDensityProblemsAD
    LogDensityProblems.capabilities(::Type{BetLyrTarget}) = LogDensityProblems.LogDensityOrder{1}()
    LogDensityProblems.logdensity_and_gradient(::BetLyrTarget, z) =
        logpost_and_grad_z(z, IDX, POS)
    LDPAD.ADgradient(kind::ADTypes.AbstractADType, lp::BetLyrTarget, r::Pigeons.Replica; kwargs...) =
        Pigeons.BufferedAD(lp, r.recorders.buffers)
    test_gradient() || error("analytic gradient does not match ForwardDiff — refusing to sample")
else
    LogDensityProblems.capabilities(::Type{BetLyrTarget}) = LogDensityProblems.LogDensityOrder{0}()
end

Pigeons.default_reference(::BetLyrTarget) =
    DistributionLogPotential(MvNormal(zeros(NDIM), Diagonal(fill(9.0, NDIM))))
const Z0 = theta_to_z(THETA_LIT, IDX)
Pigeons.initialization(::BetLyrTarget, ::AbstractRNG, ::Int) = copy(Z0)

@printf("start: logpost(z0) = %.2f  (χ² = %.1f)\n", logpost_z(Z0), chi2_total(THETA_LIT))

# Finalize any PyObject left over from reading the OIFITS *here*, on the main thread,
# before worker threads start touching Julia's GC. See the header.
GC.gc(true); GC.gc(true)

# The backend is irrelevant for :mala — the ADgradient hook above intercepts before any AD
# runs — but AutoMALA requires the field, so it is set to something valid.
explorer = EXPLORER === :mala ?
    Pigeons.AutoMALA(default_autodiff_backend = ADTypes.AutoForwardDiff()) :
    Pigeons.SliceSampler()

kw = (target = BetLyrTarget(), n_chains = N_CHAINS, multithreaded = MULTI,
      checkpoint = true, explorer = explorer,
      record = [Pigeons.traces, Pigeons.round_trip])

# ---------------------------------------------------------------------------
# Live reporting (VERBOSE=1)
# ---------------------------------------------------------------------------
# Pigeons has no per-scan callback, and "every 50 scans" would not be the right granularity
# anyway: it works in DOUBLING ROUNDS (round r is 2^r scans) and adapts the ladder and the
# explorer BETWEEN rounds, so the round boundary is the only point where the sampler's state
# is meaningfully settled. Because the rounds double, reporting at each one is automatically
# dense early and sparse late — 12 reports for N_ROUNDS=12, the last covering half the run.
#
# Running one round at a time via `increment_n_rounds!` is not an approximation of a single
# `pigeons(n_rounds=N)` call: it IS that call, paused. Verified bit-for-bit — same log(Z),
# same round-trip count, same seed behaviour.
#
# The summary is a MEDIAN over the current round's cold-chain trace, not the single latest
# state: one MCMC draw wanders by a full posterior width and reads as noise.
function report_round!(pt, r)
    Z = collect(Pigeons.get_sample(pt))
    isempty(Z) && return
    S = reduce(vcat, [transpose(z_to_theta_ad(z[1:NDIM], POS)[IDX]) for z in Z])
    med = [median(S[:, k]) for k in 1:NDIM]
    θ = copy(THETA_LIT); for (k, p) in enumerate(IDX); θ[p] = med[k]; end
    @printf("r%-3d %7d scans | logZ %11.3f | trips %4s | χ²/n %8.2f | %s\n",
            r, 2^r, Pigeons.stepping_stone(pt), string(Pigeons.n_round_trips(pt)),
            sum(chi2_split(θ)) / NTOT,
            join([@sprintf("%s=%.4g", PNAMES[p], med[k]) for (k, p) in enumerate(IDX)], " "))
end

t0 = time()
if VERBOSE
    @printf("\nexplorer = %s\n", EXPLORER === :mala ? "AutoMALA (ForwardDiff)" : "SliceSampler")
    pt = pigeons(; kw..., n_rounds = 1, show_report = false)
    report_round!(pt, 1)
    for r in 2:N_ROUNDS
        global pt = pigeons(Pigeons.increment_n_rounds!(pt, 1))
        report_round!(pt, r)
    end
else
    pt = pigeons(; kw..., n_rounds = N_ROUNDS)
end
wall = time() - t0

# ── Results ─────────────────────────────────────────────────────────────────
# get_sample returns z vectors with the log potential appended as a trailing component.
Z = collect(Pigeons.get_sample(pt))
S = reduce(vcat, [transpose(z_to_theta(z[1:NDIM], IDX)[IDX]) for z in Z])   # nsamples × NDIM
logz = Pigeons.stepping_stone(pt)
@printf("\nwall %.1f s | log(Z) = %.3f | %d samples | round trips %s\n",
        wall, logz, size(S, 1), string(Pigeons.n_round_trips(pt)))

θ_fit = copy(THETA_LIT)
for (k, p) in enumerate(IDX); θ_fit[p] = median(S[:, k]); end

println("\n" * "="^78); println("RESULT"); println("="^78)
@printf("%-16s %14s %14s %22s\n", "parameter", "literature", "fitted", "68% interval")
for (k, p) in enumerate(IDX)
    @printf("%-16s %14.6f %14.6f  [%.6f, %.6f]\n", PNAMES[p], THETA_LIT[p], θ_fit[p],
            quantile(S[:, k], 0.16), quantile(S[:, k], 0.84))
end
for (lab, θ) in (("literature", THETA_LIT), ("fitted", θ_fit))
    c = chi2_split(θ)
    @printf("\n%-11s χ²ᵥ₂/n=%8.2f  χ²ₜ₃ₐ/n=%8.2f  χ²ₜ₃ₚ/n=%8.2f  total=%8.2f",
            lab, c[1]/NV2, c[2]/NT3A, c[3]/NT3P, sum(c)/NTOT)
end
println()
if 5 in IDX
    @printf("\ndP/dt = %+.4e d/d = %+.2f s/yr   (literature %+.2f s/yr)\n",
            θ_fit[5], θ_fit[5]*365.25*86400, DPDT_SY)
end
@printf("donor %.4f mas = %.1f R☉ diam | disc FWHM %.4f mas = %.1f R☉, ratio %.2f, PA %.0f°\n",
        θ_fit[6], phys_radius_rsun(θ_fit[6], D_PC),
        θ_fit[7], phys_radius_rsun(θ_fit[7], D_PC), θ_fit[8], θ_fit[9])

println("\nper-night reduced χ² :")
@printf("%3s %-46s %9s %9s %7s\n", "#", "night", "lit", "fitted", "n")
for i in 1:NEP
    d = DATA[i]
    rc(θ) = begin
        v2m, t3am, t3pm = model_observables(θ, i)
        (sum(abs2, (v2m .- d.v2)./d.v2_err) + sum(abs2, (t3am .- d.t3amp)./d.t3amp_err) +
         sum(abs2, mod360(t3pm .- d.t3phi)./d.t3phi_err)) / NDATA[i]
    end
    @printf("%3d %-46s %9.2f %9.2f %7d\n", i, first(NIGHTS[i], 46), rc(THETA_LIT), rc(θ_fit), NDATA[i])
end

writedlm(joinpath(OUTDIR, "pigeons_posterior.txt"), S)
writedlm(joinpath(OUTDIR, "pigeons_best.txt"), hcat(PNAMES[IDX], θ_fit[IDX]))

# Round trips are THE diagnostic for multimodality: a chain that never completes one has
# not travelled between reference and target, so it has not shown it can leave the basin it
# started in. Few or zero ⇒ raise n_rounds or n_chains before believing the intervals.
nrt = Pigeons.n_round_trips(pt)
println("\nround trips: $nrt")
nrt isa Number && nrt < 5 &&
    println("⚠ few round trips — the tails are not yet trustworthy; increase N_ROUNDS/N_CHAINS.")

# Plot LAST: PyPlot is only imported here, so no PyObject exists while worker threads run.
if DOPLOT
    @eval using PythonPlot
    rd(θ, t) = orbit_to_rotir_offset((i = θ[2], Ω = θ[3], ω = OMEGA_PERI, P = P_ORB,
                                      a = θ[1], e = E_ORB, T0 = θ[4], q = Q_BIN,
                                      dP = θ[5], dω = 0.0), t)
    fig, ax = subplots(figsize = (8, 8)); ax.set_aspect("equal")
    tt = θ_fit[4] .+ range(0, P_ORB, length = 400)
    for (θ, lab, sty) in ((THETA_LIT, "literature", "--"), (θ_fit, "fitted", "-"))
        xy = [rd(θ, t) for t in tt]
        ax.plot([-p[1] for p in xy], [p[2] for p in xy], sty, lw = 1.5, label = lab)
    end
    for i in 1:NEP
        p = rd(θ_fit, TMEAN[i]); ax.plot([-p[1]], [p[2]], "o", ms = 6, mfc = "none", color = "C1")
    end
    ax.plot([0], [0], "r*", ms = 18, label = "donor")
    ax.invert_xaxis(); ax.set_xlabel("ΔRA East (mas)"); ax.set_ylabel("ΔDec North (mas)")
    ax.legend(); ax.grid(alpha = 0.3); ax.set_title("β Lyrae relative orbit — Pigeons")
    fig.savefig(joinpath(OUTDIR, "orbit_fit_pigeons.png"), dpi = 130, bbox_inches = "tight")
    PyPlot.close(fig)
end
@info "results in $OUTDIR"
