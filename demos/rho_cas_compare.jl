#!/usr/bin/env julia
# Compare every uncertainty method on ρ Cas: error bars, wall time, evidence.
#
#   julia --project=demos demos/rho_cas_compare.jl
#
# Run the producers first (any subset, in any order):
#   demos/rho_cas_pathfinder.jl   — minutes, approximate, finds modes
#   demos/rho_cas_pigeons.jl      — hours, exact, crosses modes, gives log(Z)
#   demos/rho_cas_ultranest.jl    — hours, exact, crosses modes, gives log(Z) ± error
#   demos/rho_cas_bootstrap.jl    — hours, frequentist, does *not* trust the error bars
#
# The four do not measure the same thing, and the point of the table is the disagreement:
#   • Pathfinder/Pigeons/UltraNest are Bayesian posteriors of the same model and prior, so
#     they should agree. Where they do not, the approximation (Pathfinder) or the sampling
#     (round trips, live points) is at fault.
#   • The bootstrap asks a different question — how much does the answer move if the
#     observations had been slightly different — and does not assume the quoted error bars
#     are right. With χ²ᵣ ≈ 6 the posterior widths are optimistic by roughly √6, so a
#     bootstrap noticeably wider than the posteriors is the expected, honest outcome.

using Printf, Statistics
include(joinpath(@__DIR__, "posterior_utils.jl"))

const METHODS = ["pathfinder", "pigeons", "ultranest", "bootstrap"]

results = Dict{String,Any}()
for m in METHODS
    path = joinpath(RESULTS_DIR, m * ".csv")
    isfile(path) || continue
    S, labels, meta = load_posterior(path)
    results[m] = (S = S, labels = labels, meta = meta)
end

isempty(results) && error("no results in $RESULTS_DIR — run the sampler scripts first")

labels = first(values(results)).labels
for (m, r) in results
    r.labels == labels || error("$m sampled $(r.labels), others sampled $labels")
end

# ── Table ───────────────────────────────────────────────────────────────────
println("\n", "="^92)
@printf("%-12s %8s %10s %12s", "method", "samples", "wall (s)", "log(Z)")
for l in labels; @printf("  %22s", l); end
println("\n", "="^92)

for m in METHODS
    haskey(results, m) || continue
    S, meta = results[m].S, results[m].meta
    wall = get(meta, "wall_seconds", "—")
    logz = get(meta, "logz", "—")
    @printf("%-12s %8d %10s %12s", m, size(S, 1),
            wall == "—" ? "—" : @sprintf("%.0f", parse(Float64, wall)),
            logz == "—" ? "—" : @sprintf("%.2f", parse(Float64, logz)))
    for j in eachindex(labels)
        v = view(S, :, j)
        med = median(v)
        σ = 0.5 * (quantile(v, 0.84) - quantile(v, 0.16))
        @printf("  %12.5f ±%8.5f", med, σ)
    end
    println()
end
println("="^92)

# ── Ratio of error bars, relative to the widest ─────────────────────────────
println("\nσ relative to the widest method (1.00 = widest):")
for j in eachindex(labels)
    σs = Dict(m => 0.5 * (quantile(view(results[m].S, :, j), 0.84) -
                          quantile(view(results[m].S, :, j), 0.16))
              for m in keys(results))
    σmax = maximum(values(σs))
    @printf("  %-8s", labels[j])
    for m in METHODS
        haskey(σs, m) && @printf("  %s %.2f", m, σs[m] / σmax)
    end
    println()
end

# ── Overlaid marginals ──────────────────────────────────────────────────────
try
    plt = pyimport("matplotlib.pyplot")
    fig, axes = plt.subplots(1, length(labels), figsize = (5 * length(labels), 4))
    axs = length(labels) == 1 ? [axes] : axes
    for j in eachindex(labels)
        for m in METHODS
            haskey(results, m) || continue
            axs[j].hist(Array{Float64}(view(results[m].S, :, j)); bins = 60, density = true,
                        histtype = "step", linewidth = 1.6, label = m)
        end
        axs[j].set_xlabel(labels[j]); axs[j].set_ylabel("density")
        j == 1 && axs[j].legend()
    end
    fig.suptitle("ρ Cas — marginal posteriors by method")
    path = joinpath(RESULTS_DIR, "compare_marginals.png")
    fig.savefig(path, dpi = 130, bbox_inches = "tight"); plt.close(fig)
    println("\nwrote ", path)
catch e
    @warn "comparison plot skipped" exception = e
end
