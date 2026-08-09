#!/usr/bin/env julia
# Map the χ² basins of the ρ Cas parametric fit.
#
#   julia --project=demos -t auto demos/rho_cas_basins.jl
#   NSTART=100 FREE=rpole,ld1 julia --project=demos -t auto demos/rho_cas_basins.jl
#
# The fit is multimodal: starts at different polar radii converge to different solutions,
#2.2 to 3.7 mas in diameter, with χ²ᵣ from 6 to 500. Every uncertainty in this package —
# bootstrap, posterior, analytic covariance — describes one basin, so which basin you are
# in matters more than the error bar inside it.
#
# This scans the *starting* polar radius (the other free parameters start at THETA0),
# records where each start converges, and groups the results. Embarrassingly parallel:
# NSTART independent fits.
#
# Pathfinder's per-path optima (demos/rho_cas_pathfinder.jl) answer a similar question by
# dispersing over every parameter rather than radius alone — use this when you want a
# controlled 1-D scan, that when you want coverage of the whole box.

using ROTIR, Zygote, Printf, Statistics
include(joinpath(@__DIR__, "rho_cas_model.jl"))
include(joinpath(@__DIR__, "posterior_utils.jl"))

const NSTART  = parse(Int, get(ENV, "NSTART", "60"))     # starting radii, log-spaced
const R_RANGE = (parse(Float64, get(ENV, "RMIN", "0.5")),
                 parse(Float64, get(ENV, "RMAX", "3.5")))

describe_model()
@printf("scanning %d starting radii over [%.2f, %.2f] mas (polar radius)\n",
        NSTART, R_RANGE[1], R_RANGE[2])

starts = exp.(range(log(R_RANGE[1]), log(R_RANGE[2]), length = NSTART))
res    = fill(NaN, NSTART, NDIM + 1)          # free parameters, then chi2r

tick = progress_counter(NSTART; label = "basins:")
Threads.@threads for i in eachindex(starts)
    θ0 = copy(THETA0)
    θ0[1] = starts[i]                          # rpole is always θ[1]
    try
        θ̂, chi2r, _ = fit_parametric(DATA, TESSELS, TEPOCHS, BASE;
                                     θ0 = θ0, free = FREE_NAMES, maxiter = 400)
        res[i, 1:NDIM] = θ̂[IFREE]
        res[i, end]    = chi2r
    catch e
        e isa InterruptException && rethrow(e)
    end
    tick()
end

ok = [all(isfinite, @view res[i, :]) for i in 1:NSTART]
@printf("\n%d/%d starts converged\n\n", count(ok), NSTART)

# Group solutions whose polar radius agrees to 1e-3 mas
order  = sortperm(res[ok, 1])
R      = res[ok, :][order, :]
S      = starts[ok][order]
groups = Vector{Vector{Int}}()
for i in axes(R, 1)
    if !isempty(groups) && abs(R[i, 1] - R[groups[end][end], 1]) < 1e-3
        push!(groups[end], i)
    else
        push!(groups, [i])
    end
end

print(" diameter(mas)")
for l in LABELS[2:end]; @printf(" %9s", l); end
println("      χ²ᵣ   starts   from radii")
for g in groups
    @printf("  %10.4f", 2R[g[1], 1])
    for j in 2:NDIM; @printf(" %9.4f", R[g[1], j]); end
    @printf(" %9.3f %6d   %.2f..%.2f\n", R[g[1], end], length(g), minimum(S[g]), maximum(S[g]))
end

best = argmin(R[:, end])
@printf("\nlowest χ²ᵣ = %.3f at diameter %.4f mas\n", R[best, end], 2R[best, 1])
println("Bootstrap and the samplers all need a starting point inside the basin you mean to")
println("quote — for the samplers that is THETA0 in demos/rho_cas_model.jl.")
