# NUTS sampling of the parametric rapid-rotator model.
#
# Included by ext/ROTIRHMCExt.jl. This is the demo in demos/rapid_rotator_betCas_nuts.jl made
# callable: `build_parametric_logπ` → a LogDensityProblems wrapper carrying value AND Zygote
# gradient → AdvancedHMC's NUTS, with the bounded parameters reparameterised so the sampler
# cannot propose an invalid one.
#
# WHY A HAMILTONIAN SAMPLER HERE AND NOT ANOTHER NESTED ONE.
#
# `:ultranest` and `:nautilus` explore a BOX with no gradient, and a ROTIR likelihood is
# expensive per call. NUTS uses the analytic gradient the parametric model already carries
# (src/parametric_gradient.jl), so it walks the posterior in a number of evaluations that grows
# far more slowly with dimension. What it does NOT give is an evidence — for model comparison
# the nested samplers are still the answer.

# `_box_transform` used to live here. It moved to src/parametric_fit.jl when the Pigeons
# backend needed it too: the two samplers must agree about what the unconstrained coordinates
# MEAN, or a posterior from one is not comparable with a posterior from the other — and
# sibling extensions cannot import from each other, so a shared piece has to be in the core.

"""
    ParametricPosterior

`LogDensityProblems` wrapper over an unconstrained log-posterior, with a Zygote gradient.

Value and gradient are both forced to `eltype(z)`: `Zygote.withgradient` will happily return a
Float64 value beside a Float32 gradient, and AdvancedHMC's phase point then fails on the type
mismatch. A non-finite state returns a finite `-1e30` with a ZERO gradient rather than
propagating NaN into the integrator, which is what turns a Float32 excursion into a rejected
step instead of a dead chain.
"""
struct ParametricPosterior{F}
    logπ::F
    D::Int
end
LogDensityProblems.capabilities(::Type{<:ParametricPosterior}) =
    LogDensityProblems.LogDensityOrder{1}()
LogDensityProblems.dimension(m::ParametricPosterior) = m.D
function LogDensityProblems.logdensity(m::ParametricPosterior, z)
    R = eltype(z)
    v = convert(R, m.logπ(z))
    return isfinite(v) ? v : R(-1e30)
end
function LogDensityProblems.logdensity_and_gradient(m::ParametricPosterior, z)
    R = eltype(z)
    r = Zygote.withgradient(m.logπ, z)
    v = convert(R, r.val)
    g = r.grad[1] === nothing ? zero(z) : convert(Vector{R}, r.grad[1])
    (isfinite(v) && all(isfinite, g)) || return R(-1e30), zero(z)
    return v, g
end

"""
    _fit_hmc(data_epochs, tessels, tepochs, base_params; θ0, free, lb, ub, ...) -> NamedTuple

NUTS over the free entries of the parametric θ vector.

Returns the same fields the nested samplers do — `median`, `q16`, `q84`, `mean`, `std`,
`samples` — plus `divergences` and `n_adapt`/`n_samples`, and `logz = NaN`: a Hamiltonian
sampler does not compute an evidence, and reporting a number there would be worse than
reporting that it has none.

`n_samples` and `n_adapt` are the two knobs. The defaults are a working scale, not a
production one: 300 adaptation steps is enough for the step size and a diagonal mass matrix to
settle on a smooth posterior, and 400 draws is enough to see the shape.
"""
function _fit_hmc(data_epochs, tessels, tepochs, base_params;
                  θ0, free = nothing, lb = nothing, ub = nothing,
                  tpole_free::Bool = false, intensity_model::Symbol = :linear, band = nothing,
                  κ = 50, GM = 1, n_samples::Int = 400, n_adapt::Int = 300,
                  target_accept::Real = 0.8, verb::Bool = false)
    T = eltype(tessels.unit_xyz)
    θfull = collect(T, θ0)
    dlb, dub = default_parametric_bounds(; tpole_free = tpole_free)
    lower = collect(T, lb === nothing ? dlb : lb)
    upper = collect(T, ub === nothing ? dub : ub)
    idx = parametric_free_indices(free; tpole_free = tpole_free)

    logπ = build_parametric_logπ(data_epochs, tessels, tepochs, base_params;
                                 intensity_model = intensity_model, band = band,
                                 κ = κ, GM = GM, tpole_free = tpole_free, logprior = nothing)

    # The reduced problem: θ = θ_frozen + S·θ_free, a constant scatter matrix, so Zygote
    # differentiates through it with no special handling and the frozen entries contribute
    # exactly zero gradient. Same construction as `fit_parametric`.
    θ_frozen = copy(θfull); θ_frozen[idx] .= zero(T)
    S = zeros(T, length(θfull), length(idx))
    for (j, i) in enumerate(idx); S[i, j] = one(T); end

    to_θ, to_z, ljac = _box_transform(lower[idx], upper[idx])
    logπz = z -> logπ(θ_frozen .+ S * to_θ(z)) + ljac(z)

    ℓ  = ParametricPosterior(logπz, length(idx))
    z0 = collect(T, to_z(θfull[idx]))
    isfinite(LogDensityProblems.logdensity(ℓ, z0)) ||
        error("_fit_hmc: the starting point has zero posterior density; check θ0 against the bounds")

    # A DIAGONAL metric, adapted. A dense one is what the demo takes from Pathfinder, and it is
    # better on a correlated posterior — but Pathfinder is a further dependency and a diagonal
    # metric adapted over `n_adapt` steps is what makes this runnable without one.
    # `DiagEuclideanMetric(D)` builds a Float64 metric, and the gradient here is whatever the
    # tessellation's element type is — typically Float32. AdvancedHMC then fails constructing a
    # PhasePoint from a Float32 DualValue and a Float64 one. Building the metric from a vector
    # of the right type is what keeps one float type through the whole integrator.
    metric     = AdvancedHMC.DiagEuclideanMetric(ones(T, length(idx)))
    ham        = AdvancedHMC.Hamiltonian(metric, ℓ)
    ϵ0         = AdvancedHMC.find_good_stepsize(ham, z0)
    integrator = AdvancedHMC.Leapfrog(ϵ0)
    kernel     = AdvancedHMC.HMCKernel(
        AdvancedHMC.Trajectory{AdvancedHMC.MultinomialTS}(integrator,
                                                          AdvancedHMC.GeneralisedNoUTurn()))
    adaptor = AdvancedHMC.StanHMCAdaptor(AdvancedHMC.MassMatrixAdaptor(metric),
                                         AdvancedHMC.StepSizeAdaptor(T(target_accept), integrator))

    zs, stats = AdvancedHMC.sample(ham, kernel, z0, n_samples, adaptor, n_adapt;
                                   drop_warmup = true, progress = false, verbose = verb)
    Θ = reduce(hcat, (to_θ(z) for z in zs))'          # n_samples x nfree
    ndiv = count(s -> hasproperty(s, :numerical_error) && s.numerical_error, stats)
    verb && Printf.@printf("NUTS: %d draws after %d adaptation, %d divergences\n",
                           size(Θ, 1), n_adapt, ndiv)
    ndiv > 0.05 * n_samples &&
        @warn "NUTS: $(ndiv) divergences in $(n_samples) draws — the posterior has geometry " *
              "the sampler is struggling with; raise target_accept or narrow the bounds."

    q(pr) = [Statistics.quantile(view(Θ, :, j), pr) for j in 1:length(idx)]
    return (median = q(0.5), q16 = q(0.16), q84 = q(0.84),
            mean = vec(Statistics.mean(Θ, dims = 1)),
            std = vec(Statistics.std(Θ, dims = 1)),
            samples = Matrix(Θ), logz = NaN, logzerr = NaN,
            divergences = ndiv, n_adapt = n_adapt, n_samples = n_samples, stats = stats)
end
