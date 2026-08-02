# bootstrap.jl
#
# Block bootstrap uncertainties for ROTIR's *parametric* fits.
#
# The resampling units, the replicate loop and the statistics all come from OITOOLS
# (`data_blocks`, `resample_blocks`, `bootstrap_driver`); what lives here is the part
# OITOOLS cannot know about: ROTIR data is a `Vector{OIdata}`, one entry per epoch, each
# with its own geometry, so a replicate has to resample every epoch.
#
# Scope is deliberate: this is for parameter fits (`fit_parametric`, or any closure that
# maps a resampled dataset to a parameter vector). It is NOT for image reconstruction —
# a bootstrap around a regularized inversion measures how much the data move the image
# while saying nothing about the regularization bias, and per-pixel percentiles from such
# a run read as error bars without being any.

using Random

# ─────────────────────────────────────────────────────────────────────────────
# Blocks over multiple epochs
# ─────────────────────────────────────────────────────────────────────────────

"""
    epoch_blocks(data_epochs; granularity=:config, mjd_digits=5) -> Vector{DataBlocks}

Partition each epoch of `data_epochs` into resampling blocks, one `DataBlocks` per epoch
(see OITOOLS' [`data_blocks`](@ref) for the granularities).

`:config` (default) gives one block per (MJD, telescope configuration) — the baseline for
V², the triangle for T3. `:epoch` gives one block per MJD, which inside a ROTIR epoch means
one block per exposure, since `readoifits_multiepochs` groups exposures into an epoch
without collapsing their individual MJDs. `:point` is the i.i.d. bootstrap and will
underestimate uncertainties on real (correlated) data.

MJDs are stored in `Float64` by `readoifits` whatever `T` is, so the block structure is the
same for `Float32` and `Float64` data.
"""
function epoch_blocks(data_epochs::AbstractVector; granularity::Symbol=:config,
                      mjd_digits::Int=5)
    return [data_blocks(d; granularity=granularity, mjd_digits=mjd_digits)
            for d in data_epochs]
end

_npoints(d) = d.nv2 + d.nt3amp + d.nt3phi + d.nvisamp + d.nvisphi + d.nflux
_nt3points(d) = d.nt3amp + d.nt3phi

# ─────────────────────────────────────────────────────────────────────────────
# Resampling a multi-epoch dataset
# ─────────────────────────────────────────────────────────────────────────────

function _apply_draw(d, b, mode::Symbol, rng::AbstractRNG)
    if mode === :weights
        return apply_block_weights(d, b, block_weights(length(b); rng=rng))
    else
        return apply_block_counts(d, b, block_counts(length(b), mode; rng=rng))
    end
end

"""
    resample_epochs(data_epochs, blocks; mode=:replacement, stratify=true, rng) -> Vector{OIdata}

One bootstrap replicate of a multi-epoch dataset: every epoch is resampled and the epochs
are returned in their original order, so they still line up with the geometry built for
them.

# Stratification
- `stratify=true` (default) draws blocks independently within each epoch, so every epoch
  keeps roughly its data volume. In ROTIR the epochs *are* the rotational phase coverage
  that constrains inclination, position angle and ω; a draw that empties one inflates the
  uncertainties for a reason that has nothing to do with data quality.
- `stratify=false` draws once over the pooled block set and slices the result per epoch, so
  an epoch can lose most or all of its data. Statistically the more literal bootstrap of
  "the dataset I happen to have"; warns when an epoch drops below a tenth of its points.

# Modes
See OITOOLS' [`block_counts`](@ref) and [`block_weights`](@ref). `:replacement` (default)
is the textbook nonparametric bootstrap; `:weights` is the multiplier bootstrap, which only
rescales error bars and therefore never drops a block; `:halfsample` keeps a random half.
"""
function resample_epochs(data_epochs::AbstractVector, blocks::AbstractVector{DataBlocks};
                         mode::Symbol=:replacement, stratify::Bool=true,
                         rng::AbstractRNG=Random.default_rng())
    length(blocks) == length(data_epochs) ||
        error("resample_epochs: got $(length(blocks)) block sets for " *
              "$(length(data_epochs)) epochs")

    if stratify
        return [_apply_draw(data_epochs[e], blocks[e], mode, rng)
                for e in eachindex(data_epochs)]
    end

    # Pooled draw: one multiplicity (or weight) per block over the whole dataset, then
    # sliced back per epoch. Blocks are laid out epoch after epoch.
    nb  = [length(b) for b in blocks]
    tot = sum(nb)
    draw = mode === :weights ? block_weights(tot; rng=rng) : block_counts(tot, mode; rng=rng)

    out = similar(data_epochs, Any)
    off = 0
    for e in eachindex(data_epochs)
        slice = draw[(off + 1):(off + nb[e])]
        out[e] = mode === :weights ?
                 apply_block_weights(data_epochs[e], blocks[e], slice) :
                 apply_block_counts(data_epochs[e], blocks[e], slice)
        off += nb[e]
        n0, n1 = _npoints(data_epochs[e]), _npoints(out[e])
        n1 < n0 ÷ 10 && @warn("resample_epochs: pooled draw left epoch $e with $n1 of " *
                              "$n0 points; rotational phase coverage is degraded. " *
                              "Use stratify=true unless you specifically want this.",
                              maxlog=3)
    end
    return [out[e] for e in eachindex(out)]
end

function resample_epochs(data_epochs::AbstractVector; granularity::Symbol=:config,
                         mjd_digits::Int=5, kwargs...)
    return resample_epochs(data_epochs,
                           epoch_blocks(data_epochs; granularity=granularity,
                                        mjd_digits=mjd_digits); kwargs...)
end

# ─────────────────────────────────────────────────────────────────────────────
# Result
# ─────────────────────────────────────────────────────────────────────────────

"""
    ParametricBootstrap

Outcome of [`bootstrap_parametric`](@ref): an OITOOLS `BootstrapResult` (median, 16/84
percentiles, covariance, per-replicate samples, …) plus the diagnostics that matter for a
*bounded* parametric fit. All fields of the wrapped result are reachable directly, e.g.
`b.sigma`, `b.median`, `b.samples`.

# Extra fields
- `natbound` — per parameter, how many kept replicates converged onto a box bound. A large
  count means the percentiles are shaped by the bounds, not by the data.
- `ndegenerate` — replicates whose resampled dataset lost all closure phases/amplitudes.
- `nmaxiter` — replicates whose fit stopped on the iteration limit. These are not
  minimisers, so they widen the distribution with optimiser noise; raise `maxiter` if the
  count is not near zero.
- `nblocks_per_epoch`, `stratify` — the resampling setup.
"""
struct ParametricBootstrap
    result            ::BootstrapResult
    natbound          ::Vector{Int}
    ndegenerate       ::Int
    nmaxiter          ::Int
    nblocks_per_epoch ::Vector{Int}
    stratify          ::Bool
end

function Base.getproperty(b::ParametricBootstrap, s::Symbol)
    s in fieldnames(ParametricBootstrap) && return getfield(b, s)
    return getproperty(getfield(b, :result), s)   # forward sigma/median/samples/...
end

Base.propertynames(b::ParametricBootstrap) =
    (fieldnames(ParametricBootstrap)..., propertynames(getfield(b, :result))...)

function Base.show(io::IO, b::ParametricBootstrap)
    show(io, b.result)
    @printf(io, "  blocks/epoch: %s  (stratified: %s)\n",
            string(b.nblocks_per_epoch), b.stratify ? "yes" : "no")
    b.ndegenerate > 0 &&
        @printf(io, "  %d replicate(s) lost all closure data\n", b.ndegenerate)
    b.nmaxiter > 0 &&
        @printf(io, "  WARNING %d replicate(s) stopped on the iteration limit\n", b.nmaxiter)
    if any(b.natbound .> 0)
        for (j, p) in enumerate(b.list_free_params)
            b.natbound[j] > 0 &&
                @printf(io, "  WARNING %-16s hit a bound in %d/%d kept replicates\n",
                        p, b.natbound[j], count(b.mask))
        end
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Driver
# ─────────────────────────────────────────────────────────────────────────────

_astuple(out::Tuple) = length(out) >= 3 ? (out[1], out[2], out[3]) :
                       (out[1], out[2], nothing)
_astuple(out)        = (out, NaN, nothing)

function _count_at_bound(samples, mask, lb, ub)
    npar = size(samples, 2)
    n = zeros(Int, npar)
    (lb === nothing && ub === nothing) && return n
    for j in 1:npar
        lo = lb === nothing ? -Inf : Float64(lb[j])
        hi = ub === nothing ?  Inf : Float64(ub[j])
        tol_lo = 1e-6 * max(1.0, abs(lo));  tol_hi = 1e-6 * max(1.0, abs(hi))
        for i in axes(samples, 1)
            mask[i] || continue
            x = samples[i, j]
            (isfinite(lo) && abs(x - lo) <= tol_lo) ||
                (isfinite(hi) && abs(hi - x) <= tol_hi) ? (n[j] += 1) : nothing
        end
    end
    return n
end

"""
    bootstrap_parametric(fitfun, data_epochs; kwargs...) -> ParametricBootstrap

Nonparametric block bootstrap of a parametric fit: refit `nboot` replicates in which the
blocks of data are resampled, and summarise the scatter of the best-fit parameters.

`fitfun(replicate_data) -> θ̂` is called with a resampled `Vector{OIdata}` and returns the
best-fit parameters; it may also return `(θ̂, chi2r)` or `(θ̂, chi2r, extra)`. Without a
`chi2r` the `chi2r_max` rejection is unavailable (the replicate χ²ᵣ is recorded as `NaN`).
The `using Zygote` method of this function supplies `fitfun` for you from
[`fit_parametric`](@ref).

Unlike the analytic or posterior uncertainties, this does not assume the quoted error bars
are correct or uncorrelated: correlated calibration errors and mis-stated errors appear as
extra scatter between replicates. It does assume the blocks are numerous and exchangeable.

# Keywords
- `nboot` — replicates (default 200; ≥1000 for stable percentiles)
- `mode` — `:replacement` (default), `:weights` or `:halfsample`. See
  `OITOOLS/demos/bootstrap_validation` for how each is calibrated against simulated truth
- `granularity`, `mjd_digits` — block definition, see [`epoch_blocks`](@ref)
- `stratify` — resample within each epoch (default `true`), see [`resample_epochs`](@ref)
- `lb`, `ub` — box bounds, used only to count replicates that converged onto a bound
- `list_free_params` — parameter names for the summary
- `x_opt` — full-data fit, if you already have it (otherwise `fitfun` is called once on the
  unresampled data)
- `parallel` — `:outer` (default) runs replicates on all threads; `:inner` runs them one at
  a time and leaves the threads to ROTIR's own kernels
- `sigma_clipping`, `chi2r_max`, `seed`, `verb` — passed to OITOOLS' `bootstrap_driver`
"""
function bootstrap_parametric(fitfun, data_epochs::AbstractVector;
    nboot            ::Int    = 200,
    mode             ::Symbol = :replacement,
    granularity      ::Symbol = :config,
    mjd_digits       ::Int    = 5,
    stratify         ::Bool   = true,
    lb                        = nothing,
    ub                        = nothing,
    list_free_params          = nothing,
    x_opt                     = nothing,
    parallel         ::Symbol = :outer,
    sigma_clipping            = nothing,
    chi2r_max                 = nothing,
    seed                      = nothing,
    verb             ::Bool   = true,
)
    parallel in (:outer, :inner) ||
        error("bootstrap_parametric: parallel must be :outer or :inner (got :$parallel)")

    # ── Full-data fit ────────────────────────────────────────────────────────
    if x_opt === nothing
        verb && println("Fitting the full dataset...")
        x_opt, = _astuple(fitfun(data_epochs))
    end
    npar  = length(x_opt)
    names = list_free_params === nothing ? ["p$j" for j in 1:npar] :
            collect(String, list_free_params)

    # ── Blocks ───────────────────────────────────────────────────────────────
    blocks = epoch_blocks(data_epochs; granularity=granularity, mjd_digits=mjd_digits)
    nb     = [length(b) for b in blocks]
    if verb
        @printf("Resampling %d epoch(s), %d blocks total (%s, %s)\n",
                length(blocks), sum(nb), string(granularity),
                stratify ? "stratified per epoch" : "pooled")
    end
    sum(nb) < 10 && granularity !== :point &&
        @warn "Only $(sum(nb)) resampling blocks: bootstrap uncertainties will be poorly " *
              "determined."

    # ── Replicates ───────────────────────────────────────────────────────────
    info = Vector{Any}(undef, nboot)
    fill!(info, nothing)

    driver_fitfun = function (_state, rng)
        rdata = resample_epochs(data_epochs, blocks; mode=mode, stratify=stratify, rng=rng)
        θ, chi2r, extra = _astuple(fitfun(rdata))
        nt3 = sum(_nt3points, rdata)
        return (θ, chi2r, (nt3 = nt3, npoints = sum(_npoints, rdata), fit = extra))
    end

    result = bootstrap_driver(driver_fitfun, collect(Float64, x_opt), names;
                              nboot=nboot, nworkers=(parallel === :inner ? 1 : 0),
                              info=info, sigma_clipping=sigma_clipping,
                              chi2r_max=chi2r_max, seed=seed, verb=verb,
                              mode=mode, granularity=granularity, nblocks=sum(nb))

    # ── ROTIR-specific diagnostics ───────────────────────────────────────────
    natbound = _count_at_bound(result.samples, result.mask, lb, ub)
    ndegen   = count(i -> info[i] !== nothing && info[i].nt3 == 0, 1:nboot)
    nmaxit   = count(1:nboot) do i
        f = info[i] === nothing ? nothing : info[i].fit
        f !== nothing && hasproperty(f, :hit_maxiter) && f.hit_maxiter
    end

    b = ParametricBootstrap(result, natbound, ndegen, nmaxit, nb, stratify)
    verb && (any(natbound .> 0) || ndegen > 0 || nmaxit > 0) && show(stdout, b)
    return b
end

# ─────────────────────────────────────────────────────────────────────────────
# Point estimator (methods provided by ext/ROTIRZygoteExt.jl)
# ─────────────────────────────────────────────────────────────────────────────

"""
    fit_parametric(data_epochs, tessels, tepochs, base_params; kwargs...) -> (θ̂, chi2r, info)

Best-fit parameters of the von Zeipel parametric model, by bounded minimization of
`-logπ` from [`build_parametric_logπ`](@ref).

θ = `[rpole, ω, inc, PA, β, ld1, ld2]`, with `tpole` appended when `tpole_free=true`
(identifiable only with `intensity_model=:planck`, see [`build_parametric_logπ`](@ref)).

!!! note "Requires an AD backend"
    ROTIR itself carries no AD dependency: the parametric model only defines
    ChainRulesCore rrules. `using Zygote` alongside ROTIR loads the extension that provides
    this method.
"""
function fit_parametric end

"""
    default_parametric_bounds(; tpole_free=false) -> (lb, ub)

Box bounds for [`fit_parametric`](@ref), in the θ order
`[rpole, ω, inc, PA, β, ld1, ld2]` (+ `tpole`).

`ω` stops at 0.99 rather than 1 — the rapid-rotator radius diverges at break-up — and the
limb-darkening and gravity-darkening coefficients are held in their physical ranges. These
are starting points, not physics: override them per target.

`ld1` runs to 2, not 1: for the Hestroffer power law `I ∝ μ^ld1`, cool extended stars sit
above 1 (ρ Cas fits at α ≈ 1.02), and a bound the solution wants to cross both biases the
answer and destabilises the bounded line search.
"""
function default_parametric_bounds(; tpole_free::Bool=false)
    # rpole is bounded strictly away from 0: the visibilities are normalised by the total
    # flux, so a zero radius gives 0/0 and the objective is NaN exactly at the bound. A
    # bounded line search projects trial steps onto the box, lands there, and dies.
    lb = [1e-3, 0.0,   0.0, -180.0, 0.0,  0.0, -1.0]
    ub = [Inf,  0.99, 180.0, 180.0, 1.0,  2.0,  1.0]
    if tpole_free
        push!(lb, 0.0); push!(ub, Inf)
    end
    return lb, ub
end

"""
    parametric_param_names(; tpole_free=false) -> Vector{String}

Names of the parametric θ entries, in order.
"""
parametric_param_names(; tpole_free::Bool=false) =
    tpole_free ? ["rpole", "omega", "inc", "PA", "beta", "ld1", "ld2", "tpole"] :
                 ["rpole", "omega", "inc", "PA", "beta", "ld1", "ld2"]

"""
    parametric_free_indices(free; tpole_free=false) -> Vector{Int}

Resolve which θ entries are fitted. `free` may be `nothing` (all of them), a vector of
indices, or a vector of names/symbols drawn from [`parametric_param_names`](@ref).

Fitting a subset is the normal case, not an exception. A resolved single star constrains
its angular size and its limb darkening and nothing else, so
`free = ["rpole", "ld1"]` with `omega` held at 0 *is* the sphere + limb-darkening model —
the rapid-rotator radius reduces to `rpole` everywhere at ω = 0, and with a uniform von
Zeipel map the inclination, position angle, β and tpole have no effect at all. Leaving
them free would let the optimizer wander in a flat subspace and pile replicates onto the
bounds, which is what `natbound` is there to reveal.
"""
function parametric_free_indices(free; tpole_free::Bool=false)
    names = parametric_param_names(; tpole_free=tpole_free)
    free === nothing && return collect(1:length(names))
    idx = Int[]
    for f in free
        if f isa Integer
            1 <= f <= length(names) ||
                error("parametric_free_indices: index $f outside 1:$(length(names))")
            push!(idx, Int(f))
        else
            j = findfirst(==(String(f)), names)
            j === nothing &&
                error("parametric_free_indices: unknown parameter \"$f\"; " *
                      "expected one of $(names)")
            push!(idx, j)
        end
    end
    allunique(idx) || error("parametric_free_indices: repeated parameter in $free")
    isempty(idx) && error("parametric_free_indices: no free parameters")
    return sort(idx)
end
