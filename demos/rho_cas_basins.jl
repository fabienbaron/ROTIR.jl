#!/usr/bin/env julia
# Map the χ² basins of the ρ Cas parametric fit.
#
#   julia --project=demos -t auto demos/rho_cas_basins.jl
#
# The fit is multimodal: on this data a start at 1.65 mas polar radius converges to a
# 2.51 mas diameter, and a start elsewhere lands on 5.2 mas. Every uncertainty in this
# package — bootstrap, posterior, analytic covariance — describes one basin, so which
# basin you are in matters more than the error bar inside it.
#
# This scans starting points, records where each converges, and groups the results. Cheap
# and embarrassingly parallel: NSTART independent fits.

using ROTIR, Zygote, Printf, Statistics
include(joinpath(@__DIR__, "posterior_utils.jl"))

const NSTART = parse(Int, get(ENV, "NSTART", "60"))   # starting radii, log-spaced
const R_RANGE = (0.3, 4.0)              # polar radius in mas

data    = readoifits_multiepochs(["./demos/data/rho_Cas_example.oifits"]; T = Float64)[1, :]
tepochs = [0.0]
tessels = tessellation_healpix(3, T = Float64)

base = (surface_type = 2, rpole = 1.25, tpole = 4000.0,
        ldtype = 3, ld1 = 1.75, ld2 = 0.0,
        inclination = 0.0, position_angle = 0.0, rotation_period = 1e6,
        beta = 0.08, frac_escapevel = 0.0, B_rot = 0.0)
free = ["rpole", "ld1"]

starts = exp.(range(log(R_RANGE[1]), log(R_RANGE[2]), length = NSTART))
res    = fill(NaN, NSTART, 3)           # rpole, ld1, chi2r

tick = progress_counter(NSTART; label = "basins:")
Threads.@threads for i in eachindex(starts)
    θ0 = [starts[i], 0.0, 0.0, 0.0, 0.08, 1.0, 0.0]
    try
        θ̂, chi2r, _ = fit_parametric(data, tessels, tepochs, base;
                                     θ0 = θ0, free = free, maxiter = 400)
        res[i, :] = [θ̂[1], θ̂[6], chi2r]
    catch e
        e isa InterruptException && rethrow(e)
    end
    tick()
end

ok = [all(isfinite, @view res[i, :]) for i in 1:NSTART]
@printf("%d/%d starts converged\n\n", count(ok), NSTART)

# Group solutions that agree to 1e-3 mas in radius
order   = sortperm(res[ok, 1])
R       = res[ok, :][order, :]
groups  = Vector{Vector{Int}}()
for i in axes(R, 1)
    if !isempty(groups) && abs(R[i, 1] - R[groups[end][end], 1]) < 1e-3
        push!(groups[end], i)
    else
        push!(groups, [i])
    end
end

println(" diameter (mas)     ld1       χ²ᵣ     basin reached from N starts")
for g in groups
    @printf("  %8.4f     %7.4f  %8.3f   %3d   (start radii %.2f..%.2f mas)\n",
            2R[g[1], 1], R[g[1], 2], R[g[1], 3], length(g),
            minimum(starts[ok][order][g]), maximum(starts[ok][order][g]))
end

best = argmin(R[:, 3])
@printf("\nlowest χ²ᵣ = %.3f at diameter %.4f mas, ld1 %.4f\n", R[best, 3], 2R[best, 1], R[best, 2])
println("Bootstrap and NUTS both need a starting point inside the basin you mean to quote.")
