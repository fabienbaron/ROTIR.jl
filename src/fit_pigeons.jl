# The Pigeons.jl backend for ROTIR's parametric fits.
#
# Included by ext/ROTIRPigeonsExt.jl. Nothing here is public: `_fit_pigeons` is declared in the
# core package and reached through the same `method =` switch that selects `:nautilus`,
# `:ultranest` and `:hmc`.
#
# WHY A FOURTH SAMPLER, with three already here.
#
# Non-reversible parallel tempering (Syed et al. 2021) is the one that handles a MULTIMODAL
# posterior, which is the shape ROTIR's χ² actually has: `demos/rho_cas_basins.jl` finds
# several minima between 2.2 and 3.7 mas diameter on one star. Against that:
#
#   * **NUTS samples whichever basin it started in** and reports a tight, confident, wrong
#     interval. Nothing in its output says a second basin exists.
#   * **Nested sampling does explore globally**, but its cost climbs steeply with dimension and
#     it cannot use the analytic gradient at all.
#   * Tempering moves between basins along the annealed chain ladder, uses the gradient when
#     the explorer is MALA, and returns log(Z) by stepping stone — so nothing is given up.
#
# It is also pure Julia and genuinely parallel: `multithreaded = true` scales over the
# machine's threads, which `_fit_ultranest` cannot do at all (Python may only be called from a
# thread holding the GIL, so it is pinned to one core).
#
# The diagnostic that matters here is ROUND TRIPS, not the sample count. A run whose chains
# never traverse the ladder has not visited the other modes, and its posterior is the same
# single-basin answer NUTS would have given — which is why the count is returned.

"""
    PigeonsPosterior{F,R,ORD}

`LogDensityProblems` wrapper over an unconstrained log-posterior, for Pigeons.

Separate from [`ParametricPosterior`](@ref) despite the near-identical body, because Pigeons
needs three things AdvancedHMC does not: a reference distribution to anneal from, a starting
point per replica, and a declared differentiation order that changes with the explorer.

ALL THREE TRAVEL WITH THE VALUE — the reference and `z0` as fields, the order as a type
parameter — and that shape is the whole reason this type looks over-parameterised.

Pigeons dispatches `default_reference`, `initialization` and `capabilities` on the target, so
the obvious implementation is to `@eval` a method per run. That does not work: a method
defined inside `_fit_pigeons` is in a NEWER WORLD AGE than the `pigeons(...)` call two lines
below it, and Julia does not see it — the run dies with "Attempted to call an abstract
function" from `create_reference_log_potential`, which names neither `@eval` nor world age.
Fields and a type parameter are visible immediately and cost nothing.
"""
struct PigeonsPosterior{F,R,ORD}
    logπ::F
    D::Int
    reference::R
    z0::Vector{Float64}
end

PigeonsPosterior(logπ, D::Int, reference, z0::Vector{Float64}, order::Int) =
    PigeonsPosterior{typeof(logπ),typeof(reference),order}(logπ, D, reference, z0)

# Pigeons evaluates a log potential by CALLING it, not only through LogDensityProblems.
(m::PigeonsPosterior)(z) = LogDensityProblems.logdensity(m, z)

LogDensityProblems.dimension(m::PigeonsPosterior) = m.D
# Order 0 under a slice explorer, 1 under MALA. Declaring 1 with a slice explorer only wastes
# the gradients; declaring 0 with MALA silently drops it back to slice sampling.
LogDensityProblems.capabilities(::Type{<:PigeonsPosterior{F,R,ORD}}) where {F,R,ORD} =
    LogDensityProblems.LogDensityOrder{ORD}()
Pigeons.default_reference(m::PigeonsPosterior) = m.reference
Pigeons.initialization(m::PigeonsPosterior, ::Random.AbstractRNG, ::Int) = copy(m.z0)
function LogDensityProblems.logdensity(m::PigeonsPosterior, z)
    v = Float64(m.logπ(z))
    return isfinite(v) ? v : -1e30
end
function LogDensityProblems.logdensity_and_gradient(m::PigeonsPosterior, z)
    r = Zygote.withgradient(m.logπ, z)
    v = Float64(r.val)
    g = r.grad[1] === nothing ? zero(z) : convert(Vector{Float64}, r.grad[1])
    (isfinite(v) && all(isfinite, g)) || return -1e30, zero(Vector{Float64}(undef, length(z)))
    return v, g
end

# THE GRADIENT HOOK, and the reason a hand-written gradient reaches Pigeons at all.
#
# Declaring `capabilities = LogDensityOrder{1}()` and defining `logdensity_and_gradient` is NOT
# enough: Pigeons routes every target through `LogDensityProblemsAD.ADgradient` and
# differentiates `logdensity` itself regardless. The demos verified that with a call
# counter — zero calls to the hand-written method.
#
# The hook is the three-argument form reached from
#     get_buffer(ad_buffers, :target, kind, path.target, replica)      (Pigeons/BufferedAD.jl)
# and BOTH the kind and the target must be typed. Pigeons' own method is
#     ADgradient(kind::AbstractADType, log_potential, replica::Replica)
# so typing only the target leaves the call ambiguous and the dispatch is silently lost — a
# failure mode that cost three attempts in demos/betlyr/betlyr_orbit_fit_pigeons.jl and is
# recorded there too.
const _LDPAD = Pigeons.LogDensityProblemsAD
_LDPAD.ADgradient(kind::ADTypes.AbstractADType, lp::PigeonsPosterior,
                  r::Pigeons.Replica; kwargs...) =
    Pigeons.BufferedAD(lp, r.recorders.buffers)

"""
    _fit_pigeons(data_epochs, tessels, tepochs, base_params; θ0, free, lb, ub, ...) -> NamedTuple

Non-reversible parallel tempering over the free entries of the parametric θ vector.

Returns the fields every ROTIR sampler returns — `median`, `q16`, `q84`, `mean`, `std`,
`samples`, `logz`, `logzerr` — plus `round_trips`, `n_chains` and `n_rounds`.

`n_rounds` is a POWER: the run does `2^n_rounds` scans, doubling the work at each round and
using the previous round to tune the ladder. Ten rounds is a working scale; the cost of the
last round equals everything before it, so raising it by one doubles the total.

`explorer`:

  * `:slice` (default) — `SliceSampler`, no gradient. Per-coordinate adaptive widths, which is
    what copes when the posterior is badly conditioned. MEASURED on β Lyr: the orbit posterior
    spans seven orders of magnitude in |∂χ²/∂θ| (2e4 for inclination against 2e11 for Ṗ) and a
    diagonally-preconditioned MALA does not cope with that at all.
  * `:mala` — `AutoMALA` driven by the analytic Zygote gradient through the hook above. Worth
    it on a well-conditioned posterior, where a gradient costs ~3.6x a likelihood here.

`multithreaded = true` by default, which is most of the point: unlike `_fit_ultranest` this
scales over the machine. The likelihood builds its own geometry per call and shares nothing
mutable, so it is safe to evaluate in parallel.
"""
function _fit_pigeons(data_epochs, tessels, tepochs, base_params;
                      θ0, free = nothing, lb = nothing, ub = nothing,
                      tpole_free::Bool = false, intensity_model::Symbol = :linear,
                      band = nothing, κ = 50, GM = 1,
                      n_rounds::Int = 10, n_chains::Int = 10, explorer::Symbol = :slice,
                      multithreaded::Bool = true, reference_sigma::Real = 3.0,
                      seed::Union{Nothing,Integer} = nothing, verb::Bool = false)
    explorer in (:slice, :mala) ||
        throw(ArgumentError("_fit_pigeons: explorer must be :slice or :mala (got $explorer)"))
    T = eltype(tessels.unit_xyz)
    θfull = collect(T, θ0)
    dlb, dub = default_parametric_bounds(; tpole_free = tpole_free)
    lower = collect(T, lb === nothing ? dlb : lb)
    upper = collect(T, ub === nothing ? dub : ub)
    idx = parametric_free_indices(free; tpole_free = tpole_free)
    isempty(idx) && error("_fit_pigeons: nothing is free")

    logπ = build_parametric_logπ(data_epochs, tessels, tepochs, base_params;
                                 intensity_model = intensity_model, band = band,
                                 κ = κ, GM = GM, tpole_free = tpole_free, logprior = nothing)

    # The reduced problem, exactly as `_fit_hmc` builds it: θ = θ_frozen + S·θ_free with a
    # constant scatter matrix, so Zygote differentiates through it with no special handling.
    θ_frozen = copy(θfull); θ_frozen[idx] .= zero(T)
    S = zeros(T, length(θfull), length(idx))
    for (j, i) in enumerate(idx); S[i, j] = one(T); end

    to_θ, to_z, ljac = _box_transform(lower[idx], upper[idx])
    # Float64 THROUGHOUT, unlike the NUTS path. Pigeons' reference is a `MvNormal` over
    # Float64 and its replicas carry Float64 states; feeding a Float32 target into that ladder
    # mixes element types inside the explorer rather than at a boundary where it would fail
    # loudly. The likelihood still runs at the tessellation's own precision internally.
    logπz = z -> Float64(logπ(θ_frozen .+ S * to_θ(T.(z))) + ljac(T.(z)))

    D = length(idx)
    z0 = collect(Float64, to_z(θfull[idx]))
    # The reference the ladder anneals FROM: a wide normal on the unconstrained coordinates.
    # The box transform maps every bound to the whole line with the bulk within a few units,
    # so σ = 3 covers the prior without putting most of the reference's mass where the
    # likelihood is zero.
    σ2 = Float64(reference_sigma)^2
    reference = Pigeons.DistributionLogPotential(
        Distributions.MvNormal(zeros(D), LinearAlgebra.Diagonal(fill(σ2, D))))
    ℓ = PigeonsPosterior(logπz, D, reference, z0, explorer === :mala ? 1 : 0)
    isfinite(LogDensityProblems.logdensity(ℓ, z0)) ||
        error("_fit_pigeons: the starting point has zero posterior density; check θ0 " *
              "against the bounds")

    expl = explorer === :mala ?
           Pigeons.AutoMALA(default_autodiff_backend = ADTypes.AutoZygote()) :
           Pigeons.SliceSampler()
    kw = (target = ℓ, n_chains = n_chains, n_rounds = n_rounds, multithreaded = multithreaded,
          explorer = expl, record = [Pigeons.traces, Pigeons.round_trip],
          show_report = verb)
    pt = seed === nothing ? Pigeons.pigeons(; kw...) :
                            Pigeons.pigeons(; kw..., seed = Int(seed))

    # `get_sample` yields entries carrying the D parameters PLUS the log potential as a
    # trailing component, and the array does not broadcast elementwise the way a
    # Vector{Vector} does — hence the comprehension and the `1:D` slice.
    Z = collect(Pigeons.get_sample(pt))
    Θ = reduce(vcat, (transpose(to_θ(T.(z[1:D]))) for z in Z))
    logz = Pigeons.stepping_stone(pt)
    rt = try
        Pigeons.n_round_trips(pt)
    catch
        -1
    end
    if verb
        println(Printf.@sprintf("Pigeons: %d chains, 2^%d scans, explorer=%s, %d samples, ",
                                n_chains, n_rounds, explorer, size(Θ, 1)),
                Printf.@sprintf("log(Z) = %.4f, %d round trips", logz, rt))
    end
    # ROUND TRIPS are the multimodality diagnostic. A ladder nothing traverses has not
    # visited the other modes, and the posterior below is then the single-basin answer a
    # cheaper sampler would have given — worth saying rather than leaving in the return value.
    0 <= rt < 3 && @warn "Pigeons: only $(rt) round trip(s) — the chains barely traversed " *
                         "the temperature ladder, so a multimodal posterior may not have " *
                         "been explored. Raise n_rounds or n_chains."

    q(pr) = [Statistics.quantile(view(Θ, :, j), pr) for j in 1:D]
    return (median = q(0.5), q16 = q(0.16), q84 = q(0.84),
            mean = vec(Statistics.mean(Θ, dims = 1)),
            std = vec(Statistics.std(Θ, dims = 1)),
            samples = Matrix(Θ), logz = logz,
            # Pigeons quotes no evidence error; 1/sqrt(round trips) is the honest scale and
            # NaN is better than a number with no basis when the ladder was not traversed.
            logzerr = rt > 0 ? 1 / sqrt(rt) : NaN,
            round_trips = rt, n_chains = n_chains, n_rounds = n_rounds, result = pt)
end
