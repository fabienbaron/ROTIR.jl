#!/usr/bin/env julia
# ρ Cas posterior approximation by multi-path Pathfinder.
#
#   julia --project=demos -t auto demos/rho_cas_pathfinder.jl
#   NRUNS=32 NDRAWS=4000 FREE=rpole,omega,inc,PA,beta,ld1 julia --project=demos -t auto demos/rho_cas_pathfinder.jl
#
# Pathfinder runs quasi-Newton optimisation from dispersed starting points, builds a
# normal approximation along each path, and importance-resamples the mixture. It is the
# cheapest of the four methods by a wide margin — minutes rather than hours — and with
# `nruns` dispersed starts it also *finds* the basins, which is what a single NUTS chain
# cannot do here.
#
# What it is not: an exact posterior. The approximation is a mixture of normals in the
# unconstrained space, so skewed or banana-shaped posteriors come out too symmetric, and
# the Pareto-k diagnostic below tells you when to distrust it. Use it to locate modes and
# to warm-start the expensive samplers; quote error bars from Pigeons/UltraNest.

using Pathfinder, Zygote, ADTypes, Transducers, Random, LinearAlgebra, Printf, Statistics
include(joinpath(@__DIR__, "rho_cas_model.jl"))
include(joinpath(@__DIR__, "posterior_utils.jl"))

const NRUNS  = parse(Int, get(ENV, "NRUNS", "16"))      # independent paths
const NDRAWS = parse(Int, get(ENV, "NDRAWS", "2000"))   # resampled draws
const SEED   = parse(Int, get(ENV, "SEED", "20260801"))
# Paths run sequentially by default, which wastes a many-core machine: 4 paths took 660 s
# here. Each path still threads internally through ROTIR's kernels, so expect contention —
# measure both on your machine. PARALLEL_PATHS=0 restores sequential.
const PARALLEL_PATHS = get(ENV, "PARALLEL_PATHS", "1") == "1"

describe_model()
@printf("Pathfinder: %d paths, %d draws, %d threads\n", NRUNS, NDRAWS, Threads.nthreads())

rng = Xoshiro(SEED)
inits = dispersed_starts(NRUNS; rng = rng)

result, wall = timed(() -> multipathfinder(logpost_z, NDRAWS;
                                           init = inits,
                                           nruns = NRUNS,
                                           adtype = ADTypes.AutoZygote(),
                                           executor = PARALLEL_PATHS ? Transducers.ThreadedEx() :
                                                                       Transducers.SequentialEx(),
                                           rng = rng); label = "pathfinder")

# draws are D × ndraws in the unconstrained space
Zdraws = Array(result.draws)
S = reduce(vcat, transpose.([z_to_theta(Zdraws[:, i]) for i in axes(Zdraws, 2)]))

@printf("\nwall time %.1f s   %d draws\n", wall, size(S, 1))
try
    @printf("Pareto k = %.3f  (>0.7 ⇒ the normal mixture is a poor fit; treat as mode-finding only)\n",
            result.psis_result.pareto_shape)
catch
end
summarise(S, LABELS; title = "Pathfinder approximation:")
save_posterior("pathfinder", S, LABELS; wall_seconds = wall, nruns = NRUNS)
corner_plot(S, LABELS, "pathfinder_corner.png"; title = "ρ Cas — Pathfinder (approximate)")

# Where did the individual paths end up? Distinct optima = distinct basins, which is the
# cheap version of demos/rho_cas_basins.jl and the reason to prefer multi-path here.
try
    modes = [z_to_theta(Array(r.fit_distribution.μ)) for r in result.pathfinder_results]
    M = reduce(vcat, transpose.(modes))
    println("\nper-path optima (physical units):")
    for j in 1:NDIM
        v = sort(M[:, j])
        @printf("  %-8s min %.5f  median %.5f  max %.5f   distinct(1e-3) %d\n",
                LABELS[j], v[1], median(v), v[end],
                length(unique(round.(v, digits = 3))))
    end
catch e
    @warn "per-path summary unavailable" exception = e
end
