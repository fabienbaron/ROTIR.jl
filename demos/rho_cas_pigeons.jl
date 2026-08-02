#!/usr/bin/env julia
# ρ Cas posterior by non-reversible parallel tempering (Pigeons.jl).
#
#   julia --project=demos -t auto demos/rho_cas_pigeons.jl
#   N_ROUNDS=12 N_CHAINS=16 EXPLORER=mala julia --project=demos -t auto demos/rho_cas_pigeons.jl
#   FREE=all7 N_ROUNDS=12 julia --project=demos -t auto demos/rho_cas_pigeons.jl
#
# Why tempering: the χ² of this fit is multimodal — demos/rho_cas_basins.jl finds several
# minima between 2.2 and 3.7 mas diameter. A single NUTS chain samples whichever basin it
# starts in and reports a confident, incomplete answer. Parallel tempering moves between
# basins through the annealed chain ladder, and returns log(Z) for model comparison.
#
# Gradients: Pigeons' AutoMALA/MALA explorers use them; SliceSampler does not. ROTIR's
# gradient costs ~3.6× a likelihood evaluation, so MALA is usually the better trade — but
# `EXPLORER=slice` is there because slice sampling is more robust on ragged posteriors.
# Set EXPLORER=auto (default), mala, or slice.

using Pigeons, Distributions, LogDensityProblems, Zygote, Random, Printf, Statistics
include(joinpath(@__DIR__, "rho_cas_model.jl"))
include(joinpath(@__DIR__, "posterior_utils.jl"))

const N_ROUNDS = parse(Int, get(ENV, "N_ROUNDS", "10"))     # 2^N_ROUNDS scans
const N_CHAINS = parse(Int, get(ENV, "N_CHAINS", "10"))
const EXPLORER = Symbol(get(ENV, "EXPLORER", "auto"))

describe_model()
@printf("Pigeons: %d chains, 2^%d = %d scans, explorer=%s, %d threads\n",
        N_CHAINS, N_ROUNDS, 2^N_ROUNDS, EXPLORER, Threads.nthreads())

# ── Target ──────────────────────────────────────────────────────────────────
struct RhoCasTarget end
(::RhoCasTarget)(z) = logpost_z(z)

LogDensityProblems.logdensity(::RhoCasTarget, z) = logpost_z(z)
LogDensityProblems.dimension(::RhoCasTarget) = NDIM

# Declaring order 1 (and supplying the gradient) is what lets Pigeons use AutoMALA/MALA;
# with order 0 it silently falls back to slice sampling and the gradients go unused.
LogDensityProblems.capabilities(::Type{RhoCasTarget}) = LogDensityProblems.LogDensityOrder{1}()
function LogDensityProblems.logdensity_and_gradient(::RhoCasTarget, z)
    v, back = Zygote.pullback(logpost_z, z)
    isfinite(v) || return (-Inf, zeros(eltype(z), length(z)))
    g = back(one(v))[1]
    return v, g === nothing ? zeros(eltype(z), length(z)) : g
end

# Reference: a wide Gaussian in z, which the tempering ladder anneals towards the target.
Pigeons.default_reference(::RhoCasTarget) =
    DistributionLogPotential(MvNormal(zeros(NDIM), 9.0 * I(NDIM)))
Pigeons.initialization(::RhoCasTarget, ::AbstractRNG, ::Int) = copy(Z0)

explorer = EXPLORER === :mala  ? Pigeons.AutoMALA() :
           EXPLORER === :slice ? Pigeons.SliceSampler() : nothing

# checkpoint=true writes state after every round, so a run that turns out to need more
# sampling can be extended instead of restarted:
#     pt = Pigeons.PT("results/all/latest")      # or the specific run directory
#     pt = Pigeons.increment_n_rounds!(pt, 2); pt = pigeons(pt)
# Work doubles each round, so total time ≈ 2× the last round — start lower than you think.
kw = (target = RhoCasTarget(), n_chains = N_CHAINS, n_rounds = N_ROUNDS,
      multithreaded = true, checkpoint = true,
      record = [Pigeons.traces, Pigeons.round_trip])

pt, wall = timed(() -> explorer === nothing ? pigeons(; kw...) :
                                              pigeons(; kw..., explorer = explorer))

# ── Results in physical units ───────────────────────────────────────────────
Z = Pigeons.get_sample(pt)
S = reduce(vcat, transpose.(z_to_theta.(Z)))          # nsamples × NDIM
logz = Pigeons.stepping_stone(pt)

@printf("\nwall time %.1f s   log(Z) = %.3f   %d samples\n", wall, logz, size(S, 1))
summarise(S, LABELS; title = "Pigeons posterior:")
save_posterior("pigeons", S, LABELS; wall_seconds = wall, logz = logz,
               n_chains = N_CHAINS, n_rounds = N_ROUNDS, explorer = EXPLORER)
corner_plot(S, LABELS, "pigeons_corner.png"; title = "ρ Cas — Pigeons")

# Round trips are the diagnostic that matters for multimodality: a chain that never
# completes one has not travelled between the reference and the target, so it has not
# demonstrated it can leave the basin it started in.
@printf("\nround trips: %s\n", string(Pigeons.n_round_trips(pt)))
println("Few or zero round trips ⇒ increase n_rounds or n_chains before trusting the tails.")
