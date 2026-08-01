#!/usr/bin/env julia
# NUTS sampling of the parametric rapid-rotator model on β Cas (CHARA/MIRC-X),
# using the Zygote-differentiable log-likelihood from parametric_gradient.jl:
#   build_parametric_logπ  →  LogDensityProblems wrapper (value + Zygote gradient)
#   →  Pathfinder warm-start (init draw + covariance metric)  →  AdvancedHMC NUTS.
#
# Sampled parameters: θ = [rpole, ω, inc, PA, β, ld1, ld2].
# tpole is omitted — exactly degenerate under the linear intensity proxy (see the
# ∇_tpole ≡ 0 gate in test/test_parametric_gradient.jl); use intensity_model=:planck
# (with a band wavelength) to make it identifiable.
#
# Sampling is done in an UNCONSTRAINED space z (rpole = exp(z₁), ω = ω_max·σ(z₂)) with
# the change-of-variables log-Jacobian, so NUTS can never propose rpole ≤ 0 or ω ≥ 1
# (the break-up singularity). The physical priors live on θ.

using ROTIR, Zygote, LogDensityProblems, AdvancedHMC, Pathfinder
using LinearAlgebra, Statistics, Printf, Base.Threads

# ── Data (β Cas) ────────────────────────────────────────────────────────────
oifile = joinpath(@__DIR__, "data",
                  "MEDIAN5.MIRCX_L2.2025Oct30.HD_432.MIRCX_IDL.bet_Cas.AVG10m.oifits")
# Float32 throughout: the parametric model is Float32-clean (no Float64 promotion — see
# test/test_parametric_gradient.jl), so NUTS samples ~2× faster. Run with threads for the
# ~10× kernel speedup:  julia -t auto --project=. demos/rapid_rotator_betCas_nuts.jl
@printf("Julia threads = %d  (run with `-t auto` for the threaded FT kernels)\n", Threads.nthreads())

data = readoifits_multiepochs([oifile]; T = Float32)[1, :]
tepochs = Float32.([d.mean_mjd for d in data]); tepochs .-= tepochs[1]
tessels = tessellation_healpix(3, T = Float32)

base = (surface_type = 2, rpole = 0.849f0, tpole = 7208.0f0, ldtype = 3,
        ld1 = 0.21f0, ld2 = 0.0f0, inclination = 19.9f0, position_angle = -7.09f0,
        rotation_period = 1f0/1.12f0, beta = 0.146f0, frac_escapevel = 0.92f0, B_rot = 0.0f0)

pnames = ["rpole", "ω", "inc", "PA", "β", "ld1", "ld2"]
θ0     = Float32[0.849, 0.92, 19.9, -7.09, 0.146, 0.21, 0.0]

# ── Weak physical priors + log-posterior on θ ───────────────────────────────
μ_pr = copy(θ0)
σ_pr = Float32[0.10, 0.05, 5.0, 5.0, 0.08, 0.10, 0.08]
logprior(θ) = -0.5 * sum(abs2, (θ .- μ_pr) ./ σ_pr)
logπθ = build_parametric_logπ(data, tessels, tepochs, base;
                              intensity_model = :linear, logprior = logprior)

# ── Unconstrained reparameterization  z ∈ ℝ⁷  →  θ ──────────────────────────
const ω_max = 0.999f0
logistic(x) = 1 / (1 + exp(-x))
z_to_θ(z) = [exp(z[1]), ω_max * logistic(z[2]), z[3], z[4], z[5], z[6], z[7]]
θ_to_z(θ) = [log(θ[1]), log(θ[2] / (ω_max - θ[2])), θ[3], θ[4], θ[5], θ[6], θ[7]]
function logabsjac(z)                     # log|d rpole/dz₁| + log|dω/dz₂|
    s = logistic(z[2])
    return z[1] + log(ω_max) + log(s) + log1p(-s)
end
logπz(z) = logπθ(z_to_θ(z)) + logabsjac(z)

# ── LogDensityProblems wrapper (value + Zygote gradient) ────────────────────
struct RotatorPosterior{F}
    logπ::F
    D::Int
end
LogDensityProblems.capabilities(::Type{<:RotatorPosterior}) = LogDensityProblems.LogDensityOrder{1}()
LogDensityProblems.dimension(m::RotatorPosterior) = m.D
# Force value AND gradient to eltype(z): Zygote.withgradient can return a Float64 value
# alongside a Float32 gradient, which trips AdvancedHMC's PhasePoint (DualValue type mismatch).
function LogDensityProblems.logdensity(m::RotatorPosterior, z)
    R = eltype(z)
    v = convert(R, m.logπ(z))
    return isfinite(v) ? v : R(-1f30)
end
function LogDensityProblems.logdensity_and_gradient(m::RotatorPosterior, z)
    R = eltype(z)
    r = Zygote.withgradient(m.logπ, z)
    v = convert(R, r.val)
    g = r.grad[1] === nothing ? zero(z) : convert(Vector{R}, r.grad[1])
    # Clean reject on non-finite state (Float32 leapfrog can blow up in the tails); returning
    # a finite -1e30 with zero gradient keeps NUTS from propagating NaN/Inf into its integrator.
    (isfinite(v) && all(isfinite, g)) || return R(-1f30), zero(z)
    return v, g
end

ℓ  = RotatorPosterior(logπz, length(θ0))
z0 = θ_to_z(θ0)
@printf("logπ(z0) = %.3f\n", LogDensityProblems.logdensity(ℓ, z0))

# ── Pathfinder warm-start (unconstrained space; robust to Float32 hiccups) ──
D = length(θ0)
z_init = collect(z0)
Σ_init = Matrix{Float32}(I, D, D)
println("\nRunning Pathfinder…")
try
    pf = pathfinder(ℓ; init = collect(z0), ndraws = 200)
    z_init = Float32.(pf.draws[:, 1])
    Σ_init = Float32.(Matrix(pf.fit_distribution.Σ))     # Laplace-ish covariance → mass matrix
    println("Pathfinder mode (θ):")
    θ_mode = z_to_θ(Float32.(pf.optim_solution))
    for j in eachindex(pnames)
        @printf("  %-6s %9.4f\n", pnames[j], θ_mode[j])
    end
catch e
    @warn "Pathfinder failed — using θ0 init + identity metric" exception = (e, catch_backtrace())
end

# ── AdvancedHMC NUTS (Pathfinder covariance as the initial metric) ──────────
n_adapt, n_samples = 300, 400   # demo scale; increase for production inference
metric     = DenseEuclideanMetric(Σ_init)
ham        = Hamiltonian(metric, ℓ)
ϵ0         = find_good_stepsize(ham, collect(z_init))
integrator = Leapfrog(ϵ0)
kernel     = HMCKernel(Trajectory{MultinomialTS}(integrator, GeneralisedNoUTurn()))
adaptor    = StanHMCAdaptor(MassMatrixAdaptor(metric), StepSizeAdaptor(0.8f0, integrator))  # 0.8f0: match Float32 integrator

println("\nRunning NUTS ($n_samples samples, $n_adapt adaptation)…")
zsamples, stats = sample(ham, kernel, collect(z_init), n_samples, adaptor, n_adapt;
                         drop_warmup = true, progress = false, verbose = false)

# ── Transform back to physical θ and summarize ──────────────────────────────
Θ = reduce(hcat, (z_to_θ(z) for z in zsamples))'      # n_samples × D
ndiv = count(s -> getproperty(s, :numerical_error), stats)
@printf("\nDivergences: %d / %d\n", ndiv, length(stats))
println("\nβ Cas posterior (mean ± std, 16/84 %):")
for j in 1:D
    q = quantile(Θ[:, j], (0.16, 0.84))
    @printf("  %-6s  %9.4f ± %.4f   [%.4f, %.4f]\n",
            pnames[j], mean(Θ[:, j]), std(Θ[:, j]), q[1], q[2])
end
