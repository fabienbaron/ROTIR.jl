#!/usr/bin/env julia
# betlyr_orbit_fit_ultranest.jl — β Lyrae orbital elements by nested sampling (UltraNest, via PythonCall).
#
#   julia --project=demos betlyr/betlyr_orbit_fit_ultranest.jl
#   LIVE=400 FITDP=1 FREESIZE=1 OLDDATA=1 julia --project=demos betlyr/betlyr_orbit_fit_ultranest.jl
#
# The model lives in betlyr_model.jl, shared with the Pigeons driver
# (betlyr_orbit_fit_pigeons.jl). This file is only the sampler.
#
# Single-threaded. Under PyCall this was mandatory (a PyObject finalized on a worker thread
# segfaulted); PythonCall defers decrefs behind a GIL check, so that hazard is gone and a
# threaded session is safe. Python must still only be CALLED from a thread holding the GIL,
# so do not wrap the likelihood in `Threads.@threads` — use the Pigeons driver for cores.
#
# Knobs: LIVE=400  DATASET=old|new|all  FITDP=1  FREESIZE=1  DRYRUN=1  BENCH=1
#        STEP=1  NSTEPS=25  FRACREMAIN=0.01  OUTDIR=...

ENV["MPLBACKEND"] = get(ENV, "MPLBACKEND", "Agg")
using PythonCall, PythonPlot, Printf, Statistics, DelimitedFiles
include(joinpath(@__DIR__, "betlyr_model.jl"))

const nlive  = parse(Int, get(ENV, "LIVE", "400"))
const outdir = get(ENV, "OUTDIR", joinpath(@__DIR__, "results"))
mkpath(outdir)
const IDX = free_indices()
describe_model(IDX)

if get(ENV, "BENCH", "0") == "1"
    b(f, n=100) = (f(); minimum(begin s=time_ns(); f(); (time_ns()-s)/1e6 end for _=1:n))
    @printf("\nlikelihood %.3f ms\n", b(() -> chi2_total(THETA_LIT)))
end
if get(ENV, "DRYRUN", "0") == "1"
    @info "DRYRUN=1: stopping before sampling"; exit(0)
end

# ── UltraNest ───────────────────────────────────────────────────────────────
# Vectorized, following OITOOLS' fit_model_ultranest idiom. Three PythonCall details:
#   * DECLARED argument types (`::AbstractMatrix{<:Real}`) make PythonCall convert the numpy
#     batch into a Julia matrix; without them the argument stays a `Py`.
#   * Returns need `.to_numpy()`: a bare Julia array arrives as a `juliacall.VectorValue`,
#     on which UltraNest calls numpy methods.
#   * `pylist(names)`: UltraNest evaluates `names + [...]`, which a VectorValue rejects.
println("\n" * "="^78)
println("Direct fit of the orbital elements to the visibilities (UltraNest)")
println("="^78)
const Δ = HI[IDX] .- LO[IDX]
prior_1(u::AbstractVector{<:Real}) = u .* Δ .+ LO[IDX]
prior_v(U::AbstractMatrix{<:Real}) = Py(reduce(vcat, (u -> prior_1(u)').(eachrow(U)))).to_numpy()
function loglike_1(p::AbstractVector{<:Real})
    θ = copy(THETA_LIT); θ[IDX] .= p
    v = -0.5 * chi2_total(θ)
    return isfinite(v) ? v : -1e30
end
# THREADED over the batch. `vectorized = true` only means UltraNest hands over a whole batch
# per call — the evaluation itself was a serial broadcast, so with JULIA_NUM_THREADS=16 the
# other 15 threads sat idle. `chi2_total` is pure Julia and reads only const globals, so the
# rows are independent.
#
# This is safe under PythonCall specifically. X has ALREADY been converted to a Julia matrix
# (the ::AbstractMatrix annotation does that), and the result is converted back on this
# thread after the loop, so no worker thread ever calls Python — which is the one thing that
# still segfaults. Under PyCall the pattern was unsafe regardless, because a PyObject
# finalizer could fire on a worker thread; PythonCall defers decrefs behind a GIL check.
#
# THREADED=0 falls back to the serial broadcast.
const THREADED = get(ENV, "THREADED", "1") == "1"
function loglike_v(X::AbstractMatrix{<:Real})
    THREADED || return Py(loglike_1.(eachrow(X))).to_numpy()
    out = Vector{Float64}(undef, size(X, 1))
    Threads.@threads for k in axes(X, 1)
        @inbounds out[k] = loglike_1(view(X, k, :))
    end
    return Py(out).to_numpy()
end

ultranest = pyimport("ultranest")
# WRAPPED (circular) PARAMETERS. T0 spans exactly one period and disc PA exactly one
# symmetry period of the ellipse, so in both cases the two ends of the prior are the SAME
# point. Without telling UltraNest that, a mode sitting near a boundary is split in two and
# the region has to span the whole range to cover both halves — which is a large part of why
# efficiency collapses here (the T0 panel shows satellite clumps at both edges).
const WRAPPED = pylist([p in (4, 9) for p in IDX])   # 4 = T0_JD, 9 = disc_pa_deg
sampler = ultranest.ReactiveNestedSampler(pylist(PNAMES[IDX]), loglike_v;
                                          transform = prior_v, vectorized = true,
                                          wrapped_params = WRAPPED)
# STEP SAMPLER. With ~2000 points at chi2/n ~ 37 the likelihood is ~exp(-36000), so the
# posterior occupies a vanishing fraction of a 9-D box prior. Plain region rejection cannot
# compress that: UltraNest reports "0/2494 accepted in iteration 41" and effectively stalls.
# A slice sampler walks along the constrained region instead of rejecting into it, which is
# what UltraNest's own warning recommends. `RegionSliceSampler` with adaptive step counts is
# the configuration OITOOLS uses in fit_model_ultranest.
#
# STEP=0 restores plain rejection sampling (fine for well-conditioned, low-dimensional fits).
if get(ENV, "STEP", "1") == "1"
    ss = pyimport("ultranest.stepsampler")
    sampler.stepsampler = ss.RegionSliceSampler(nsteps = parse(Int, get(ENV, "NSTEPS", "25")),
                                                adaptive_nsteps = "move-distance")
end
res = sampler.run(min_num_live_points = nlive,
                  frac_remain = parse(Float64, get(ENV, "FRACREMAIN", "0.01")),
                  show_status = true)
sampler.print_results()

post = pyconvert(Matrix{Float64}, res["samples"])
θ_fit = copy(THETA_LIT); θ_fit[IDX] .= [median(post[:, k]) for k in 1:length(IDX)]

println("\n" * "="^78); println("RESULT"); println("="^78)
@printf("%-16s %14s %14s %22s\n", "parameter", "literature", "fitted", "68% interval")
for (k, p) in enumerate(IDX)
    @printf("%-16s %14.6f %14.6f  [%.6f, %.6f]\n", PNAMES[p], THETA_LIT[p], θ_fit[p],
            quantile(post[:, k], 0.16), quantile(post[:, k], 0.84))
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
writedlm(joinpath(outdir, "ultranest_posterior.txt"), post)
writedlm(joinpath(outdir, "ultranest_best.txt"), hcat(PNAMES[IDX], θ_fit[IDX]))

# ── Corner plot ─────────────────────────────────────────────────────────────
# `ultranest.plot.cornerplot`, the same entry point OITOOLS uses in `fit_model_ultranest`,
# so ROTIR and OITOOLS posteriors are drawn identically. It takes the UltraNest result
# object directly and forwards extra keywords to `corner.corner`, which is how the
# literature values get overplotted as `truths` — every panel then shows the posterior
# against the published element with no separate comparison needed.
#
# `ultranest.plot` imports scipy, which UltraNest's core sampler does not. ROTIR's
# CondaPkg.toml declares scipy (and corner) for exactly this reason; without it the fit
# succeeds and only the plot fails. Wrapped in try/catch so that cannot discard the run.
try
    pyimport("matplotlib").rcParams["font.size"] = 8
    pyimport("ultranest.plot").cornerplot(res;
        truths      = pylist(Float64.(THETA_LIT[IDX])),
        truth_color = "tab:red",
        color       = "black",
        contour_kwargs = pydict(Dict(
            "colors"     => pylist(["#0072B2", "#56B4E9", "#009E73", "#F0E442"]),
            "linestyles" => pylist(["-", "-", "-", "-"]))))
    cpath = joinpath(outdir, "betlyr_corner.png")
    pyplot.savefig(cpath, dpi = 150, bbox_inches = "tight")
    pyplot.close("all")
    @info "corner plot: $cpath  (red = literature values)"
catch err
    @warn "corner plot failed; fit results are still written" err
end

@info "results in $outdir"
