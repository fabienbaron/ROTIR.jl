# ultranest.jl
#
# Nested sampling of the parametric model through the Python `ultranest` package.
#
# Included by ext/ROTIRUltraNestExt.jl, which loads when PythonCall does. It was core until
# the GUI measured what `using PythonCall` costs — see the header of src/fit_ultranest.jl.
# The Python side is reached exactly as OITOOLS' `fit_model_ultranest` does it.
#
# When to reach for it instead of the optimiser + bootstrap:
#   - the χ² is multimodal (ρ Cas has several minima between 2.2 and 3.7 mas diameter) and
#     nested sampling explores all of them without a starting point;
#   - you want the evidence log(Z) to compare models — sphere vs oblate vs oblate with
#     gravity darkening — rather than just fit one;
#   - you would rather spend evaluations than tune warm starts.
# Against: it ignores the gradients the model can supply (a gradient costs ~3.6× a
# likelihood here), and the cost grows steeply with the number of free parameters.


"""
    fit_parametric_ultranest(data_epochs, tessels, tepochs, base_params; kwargs...)

Sample the posterior of the parametric von Zeipel model by nested sampling, using the
Python `ultranest` package (`pip install ultranest`).

Uniform priors over the box `[lb, ub]` for every free parameter — nested sampling needs
finite bounds, so `Inf` is rejected. Returns a NamedTuple with `θ_mean`, `θ_median`,
`θ_std`, `q16`/`q84`, `logz`, `logzerr`, `samples` (nsamples × nfree), `free`,
`list_free_params` and the raw `result` dictionary.

# Keywords
- `θ0` — values for the frozen entries (required; free entries are ignored)
- `free` — which entries to sample: `nothing` (all), indices, or names
- `lb`, `ub` — box bounds; default [`default_parametric_bounds`](@ref) with `rpole`'s
  infinite upper bound replaced by `rpole_max`
- `rpole_max` — finite upper bound for `rpole` in mas (default 10)
- `intensity_model`, `band`, `κ`, `GM`, `tpole_free` — passed to `build_parametric_logπ`
- `min_num_live_points` (400), `frac_remain` (1e-3), `use_stepsampler`, `nsteps` (400),
  `log_dir`, `resume`, `verb` — UltraNest settings. Above ~5 free parameters the region
  sampler struggles; set `use_stepsampler=true` for a slice sampler.

!!! note "Cost"
    One likelihood evaluation is ~0.13 s at `n=3`/Float64 on ρ Cas. Two free parameters
    need ~10–20k of them; seven with a step sampler can need 10⁵–10⁶. Budget accordingly.
"""
function fit_parametric_ultranest(data_epochs::AbstractVector, tessels, tepochs, base_params;
    θ0,
    free                     = nothing,
    lb                       = nothing,
    ub                       = nothing,
    rpole_max       ::Real   = 10.0,
    intensity_model ::Symbol = :linear,
    band                     = nothing,
    κ                        = 50,
    GM                       = 1,
    tpole_free      ::Bool   = false,
    logprior                 = nothing,
    min_num_live_points::Int = 400,
    frac_remain              = 1e-3,
    use_stepsampler ::Bool   = false,
    nsteps          ::Int    = 400,
    log_dir                  = nothing,
    resume          ::String = "overwrite",
    verb            ::Bool   = true,
)
    T   = eltype(tessels.unit_xyz)
    θ   = collect(T, θ0)
    idx = parametric_free_indices(free; tpole_free=tpole_free)
    names = parametric_param_names(; tpole_free=tpole_free)[idx]

    dlb, dub = default_parametric_bounds(; tpole_free=tpole_free)
    dub[1] = rpole_max                      # nested sampling cannot use an infinite prior
    lo = Float64.(collect(lb === nothing ? dlb : lb)[idx])
    hi = Float64.(collect(ub === nothing ? dub : ub)[idx])
    all(isfinite, lo) && all(isfinite, hi) ||
        error("fit_parametric_ultranest: every free parameter needs finite bounds " *
              "(got lb=$lo, ub=$hi)")
    all(hi .> lo) || error("fit_parametric_ultranest: need ub > lb, got $lo / $hi")

    logπ = build_parametric_logπ(data_epochs, tessels, tepochs, base_params;
                                 intensity_model=intensity_model, band=band,
                                 κ=κ, GM=GM, tpole_free=tpole_free, logprior=logprior)

    # VECTORISED likelihood and transform, following OITOOLS' fit_model_ultranest. UltraNest
    # hands over a whole BATCH of points as an n x d numpy array and expects arrays back, so
    # the Julia<->Python crossing is paid once per batch rather than once per live point.
    #
    # Three PythonCall details:
    #   * ARGUMENTS arrive as `Py`. Annotating the closure `::AbstractMatrix{<:Real}` makes
    #     PythonCall convert the numpy batch to a Julia matrix; without it, `collect` raises
    #     `MethodError: no method matching (Array{Float64})(::Py)`.
    #   * RETURN values need `.to_numpy()`. A bare Julia array reaches Python as a
    #     `juliacall.VectorValue`, on which UltraNest calls numpy methods (`.transpose`).
    #   * `names` must be a real Python list: UltraNest evaluates `names + [...]`, and a
    #     `VectorValue` does not support `+`.
    #
    # UltraNest applies `transform` itself and passes the likelihood the PHYSICAL
    # parameters, NOT the unit cube. Do NOT re-apply the transform inside the likelihood:
    # that double-maps every parameter (a radius of 1.465 with bounds [0.05, 20] becomes
    # 29.3) and the sampler then converges confidently on nonsense — tight error bars around
    # a chi2 many orders of magnitude worse than the true optimum.
    transform_v = let lo = lo, hi = hi
        (U::AbstractMatrix{<:Real}) ->
            Py(reduce(vcat, (u -> (lo .+ (hi .- lo) .* u)').(eachrow(U)))).to_numpy()
    end
    function loglike_1(params)
        θtry = copy(θ)
        θtry[idx] .= T.(collect(Float64, params))
        v = logπ(θtry)
        # NaN/Inf (e.g. a degenerate geometry) must not reach the sampler
        return isfinite(v) ? Float64(v) : -1e30
    end
    loglike_v = (X::AbstractMatrix{<:Real}) -> Py(loglike_1.(eachrow(X))).to_numpy()

    # Safe in a multithreaded session: PythonCall's finalizer enqueues the pointer
    # (`GC.enqueue`) and the queue is drained only by a thread holding the GIL
    # (`PyGILState_Check`, PythonCall/src/GC/GC.jl), so Py objects may be finalized on any
    # thread. Verified with a full UltraNest run under `julia -t 8`.
    #
    # Python must still only be CALLED from a thread holding the GIL. Driving the sampler
    # from one thread (as here) is fine; calling Python inside a `Threads.@threads` body
    # without `PythonCall.GIL.@lock` segfaults.
    ultranest = pyimport("ultranest")
    GC.gc(true); GC.gc(true)      # finalize stray PyObjects on the main thread first

    # `log_dir` makes UltraNest store points in HDF5, which needs h5py. Rather than dying
    # after the setup work — and, for a resumed run, after hours of sampling — fall back to
    # running without a log directory and say what was lost.
    if log_dir !== nothing
        try
            pyimport("h5py")
        catch
            @warn "log_dir requires the Python h5py package; continuing without it " *
                  "(no resume, no on-disk monitoring). Install it into the environment " *
                  "PythonCall is using, e.g. via CondaPkg."
            log_dir = nothing
        end
    end

    # `names` MUST be a Python list. PythonCall wraps a Julia Vector as a
    # `juliacall.VectorValue`, and UltraNest concatenates paramnames with a list
    # internally, which then raises `TypeError: unsupported operand type(s) for +`.
    sampler = ultranest.ReactiveNestedSampler(pylist(names), loglike_v;
                                              transform = transform_v, vectorized = true,
                                              log_dir = log_dir, resume = resume)
    if use_stepsampler
        stepsampler = pyimport("ultranest.stepsampler")
        sampler.stepsampler = stepsampler.SliceSampler(
            nsteps = nsteps,
            generate_direction = stepsampler.generate_mixture_random_direction)
    end

    res = sampler.run(min_num_live_points = min_num_live_points,
                      frac_remain = frac_remain,
                      show_status = verb)
    verb && sampler.print_results()

    # `res[...]` entries are `Py`; convert here rather than leaking them to callers.
    S = pyconvert(Matrix{Float64}, res["samples"])        # nsamples × nfree
    q(p) = [quantile(view(S, :, j), p) for j in eachindex(idx)]
    out = (θ_mean   = vec(mean(S, dims=1)),
           θ_median = q(0.5),
           θ_std    = vec(std(S, dims=1)),
           q16      = q(0.16),
           q84      = q(0.84),
           logz     = pyconvert(Float64, res["logz"]),
           logzerr  = pyconvert(Float64, res["logzerr"]),
           samples  = S,
           free     = idx,
           list_free_params = names,
           result   = res)

    if verb
        @printf("\nlog(Z) = %.3f ± %.3f\n", out.logz, out.logzerr)
        for (j, p) in enumerate(names)
            @printf("  %-8s %10.5f  +%.5f -%.5f\n", p, out.θ_median[j],
                    out.q84[j] - out.θ_median[j], out.θ_median[j] - out.q16[j])
        end
    end
    return out
end
