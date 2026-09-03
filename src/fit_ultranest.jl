# The UltraNest backend for ROTIR's parametric and shape fits.
#
# Included by ext/ROTIRUltraNestExt.jl, which loads when PythonCall does. Nothing here is
# public: `_fit_ultranest` is declared in the core package and reached through the same
# `method = :ultranest` switch that selects `:nautilus`.
#
# WHY THIS IS AN EXTENSION AND NOT CORE, which it was until the GUI measured the cost.
#
# `using PythonCall` invalidates the plot-construction code OITOOLS precompiles for its live
# canvas: MEASURED, `build_canvas` goes from 338 ms to 2477 ms with PythonCall loaded and
# nothing else changed. That is 1.2 s added to every ROTIR GUI start, paid by every session
# whether or not it ever samples anything. PythonCall also cannot go into a PackageCompiler
# build, and matplotlib's backend probe maps a second Qt into the process — the collision the
# GUI exists downstream of.
#
# OITOOLS made PythonCall weak for exactly these reasons. ROTIR keeping it hard undid that
# decision for anyone who loaded ROTIR, which is everyone using this GUI.
#
# `:nautilus` is the pure-Julia alternative and needs no Python at all; MEASURED against this
# backend on a λ And sphere + limb-darkening fit, the two agree to 1.4e-6 mas in radius.

function _fit_ultranest(chi2_of, names, lo, hi;
                        min_num_live_points::Int = 400, frac_remain = 1e-3,
                        use_stepsampler::Bool = false, nsteps::Int = 400, verb::Bool = false)
    # Safe in a multithreaded session: PythonCall's finalizer enqueues the pointer
    # (`GC.enqueue`) and the queue is drained only by a thread holding the GIL
    # (`PyGILState_Check`, PythonCall/src/GC/GC.jl), so Py objects may be finalized on any
    # thread. Verified with a full UltraNest run under `julia -t 8`.
    #
    # Python must still only be CALLED from a thread holding the GIL. Driving the sampler
    # from one thread (as here) is fine; calling Python inside a `Threads.@threads` body
    # without `PythonCall.GIL.@lock` segfaults.
    ultranest = pyimport("ultranest")
    # Finalize any stray PyObject on the main thread before sampling starts.
    GC.gc(true); GC.gc(true)
    # UltraNest applies `transform` itself and hands `loglike` the PHYSICAL parameters, not
    # the unit cube. Do NOT apply the transform a second time inside `loglike`: that maps a
    # radius of 1.465 to 0.05 + 19.95*1.465 = 29.3 mas, and the sampler converges confidently
    # on nonsense — D = 0.14 mas with a ±0.0002 error bar at χ²/n = 2e7, against a true
    # optimum of 2.93 mas at χ²/n = 4.8.
    # VECTORISED likelihood and transform, following OITOOLS' fit_model_ultranest. UltraNest
    # hands over a whole BATCH of points as an n x d numpy array and expects arrays back, so
    # the Julia<->Python crossing is paid once per batch instead of once per live point.
    #
    # Three PythonCall details:
    #   * ARGUMENTS arrive as `Py`. Annotating the closure `::AbstractMatrix{<:Real}` makes
    #     PythonCall convert the numpy batch to a Julia matrix; without it, `collect` raises
    #     `MethodError: no method matching (Array{Float64})(::Py)`.
    #   * RETURN values need `.to_numpy()`. A bare Julia array reaches Python as a
    #     `juliacall.VectorValue`, on which UltraNest calls numpy methods (`.transpose`).
    #   * `names` must be a real Python list: UltraNest evaluates `names + [...]`, and a
    #     `VectorValue` does not support `+`.
    transform_v = let lo = lo, hi = hi
        (U::AbstractMatrix{<:Real}) ->
            Py(reduce(vcat, (u -> (lo .+ (hi .- lo) .* u)').(eachrow(U)))).to_numpy()
    end
    loglike_1(p) = (c = chi2_of(collect(Float64, p)); isfinite(c) ? -0.5*Float64(c) : -1e30)
    loglike_v = (X::AbstractMatrix{<:Real}) -> Py(loglike_1.(eachrow(X))).to_numpy()

    sampler = ultranest.ReactiveNestedSampler(pylist(names), loglike_v;
                                              transform = transform_v, vectorized = true)
    if use_stepsampler
        ss = pyimport("ultranest.stepsampler")
        sampler.stepsampler = ss.SliceSampler(nsteps = nsteps,
                                  generate_direction = ss.generate_mixture_random_direction)
    end
    res = sampler.run(min_num_live_points = min_num_live_points,
                      frac_remain = frac_remain, show_status = verb)
    verb && sampler.print_results()
    # `res[...]` entries are `Py`; convert on the way out rather than leaking them to callers.
    S = pyconvert(Matrix{Float64}, res["samples"])
    q(pr) = [quantile(view(S, :, j), pr) for j in eachindex(names)]
    return (median = q(0.5), q16 = q(0.16), q84 = q(0.84),
            mean = vec(mean(S, dims = 1)), std = vec(std(S, dims = 1)),
            logz = pyconvert(Float64, res["logz"]),
            logzerr = pyconvert(Float64, res["logzerr"]), samples = S, result = res)
end
