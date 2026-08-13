# parametric_fit.jl — fitting purely parametric stars (no surface map, no regularizer).
#
# ROTIR's existing `fit_parametric` / `build_parametric_logπ` are hardcoded to the RAPID
# ROTATOR vector θ = [rpole, fev, inc, PA, β, ld1, ld2]. That is the right model for β Cas or
# Altair, and the wrong one for a slow rotator: a red supergiant or a Cepheid has no measurable
# oblateness and no gravity darkening, so `fev` and `β` are not merely unconstrained but
# degenerate, and the fit wanders.
#
# These wrappers fit the SHAPE parameters of a star whose surface is otherwise uniform, using
# NLopt's Nelder-Mead (`:LN_NELDERMEAD`) — derivative-free, bound-constrained, and the
# best-established simplex implementation available. No gradient is needed at 2–5 parameters,
# and the χ² surface here is smooth but not cheap, so robustness matters more than order.
#
# Typical use: get the angular size before a snapshot reconstruction, since the tessellation
# is built at a fixed `radius` and getting it wrong is not a small error — a 55 % radius error
# on RW Cep gave χ²/n ≈ 4000 versus ≈ 4 once fitted.


# ─── Shared UltraNest driver for the wrappers ────────────────────────────────
# Nested sampling over a uniform box prior on the same parameters the optimiser fits. Kept
# private: `fit_sphere_ld` / `fit_ellipsoid_ld` expose it through `method = :ultranest`, so
# the caller never has to know the Python side exists.
#
# What it buys over Nelder-Mead: UNCERTAINTIES and log(Z). The optimiser returns a point with
# no error bar, which is useless when the question is "how well is the limb-darkening
# coefficient constrained" — exactly the question RADFLAT exists to answer. log(Z) also lets
# you compare sphere vs ellipsoid honestly, which Δχ² cannot do when χ²/n ≫ 1.
#
# Against: no gradients are used, cost grows steeply with dimension, and PyCall is not
# thread-safe. Safe to call from a multithreaded session (see the note in the body).
function _fit_ultranest(chi2_of, names, lo, hi;
                        min_num_live_points::Int = 400, frac_remain = 1e-3,
                        use_stepsampler::Bool = false, nsteps::Int = 400, verb::Bool = false)
    # NOTE: no thread guard here any more. Under PyCall this had to fail fast, because
    # `pydecref_` called `Py_DecRef` straight from a Julia GC finalizer with no GIL check,
    # so a PyObject finalized on a worker thread killed the process with SIGSEGV part-way
    # through a run. PythonCall does not have that failure mode: its finalizer enqueues the
    # pointer (`GC.enqueue`), and the queue is only drained by a thread that actually holds
    # the GIL (`PyGILState_Check`, PythonCall/src/GC/GC.jl). Verified empirically -- a full
    # UltraNest run completes under `julia -t 8` with Py objects being finalized on worker
    # threads throughout.
    #
    # Still true: Python must only be CALLED from a thread holding the GIL. Driving the
    # sampler from one thread (as here) is fine; calling Python from inside a
    # `Threads.@threads` body without `PythonCall.GIL.@lock` still segfaults.
    ultranest = pyimport("ultranest")
    # Finalize any stray PyObject on the main thread before sampling starts.
    GC.gc(true); GC.gc(true)
    # UltraNest applies `transform` itself and hands `loglike` the PHYSICAL parameters, not
    # the unit cube. Applying the transform a second time inside `loglike` maps e.g. a radius
    # of 1.465 to 0.05 + 19.95*1.465 = 29.3 mas and the sampler converges confidently on
    # nonsense — it did exactly that here, returning D = 0.14 mas with a ±0.0002 error bar at
    # χ²/n = 2e7, against a true optimum of 2.93 mas at χ²/n = 4.8.
    # VECTORISED likelihood and transform, following OITOOLS' fit_model_ultranest. UltraNest
    # hands over a whole BATCH of points as an n x d numpy array and expects arrays back, so
    # the Julia<->Python crossing is paid once per batch instead of once per live point.
    #
    # Three PythonCall details, each of which PyCall used to hide:
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

"""
    parametric_chi2(params, tessels, data, tepochs; weights=[1,1,1]) -> χ²

χ² of a purely parametric star: build the geometry, generate its map from `params` alone, and
compare to `data`. No surface reconstruction and no regularizer are involved.

`weights` scales `[V², T3amp, T3φ]`. The default weights everything equally; see
[`fit_sphere_ld`](@ref) for why you usually do **not** want T3φ in a symmetric fit.

`uniform = true` forces a UNIFORM surface brightness (`tpole` everywhere) instead of calling
`parametric_temperature_map`, so the only structure left is limb darkening. For
`surface_type = 0` that is what the map already is, so it changes nothing; for the ellipsoid
it matters, because `temperature_map_vonZeipel_ellipsoid` is a placeholder gravity model
(`g ∝ 1/r²`, with an unverified `rpole = radius_x`) that you do not want silently folded into
a shape fit.
"""
function parametric_chi2(params, tessels, data, tepochs;
                         weights = [1.0, 1.0, 1.0], uniform::Bool = false)
    stars = create_star_multiepochs(tessels, params, tepochs)
    setup_oi!(data, stars)
    x = uniform ? fill(eltype(stars[1])(params.tpole), stars[1].npix) :
                  parametric_temperature_map(params, stars[1])
    c = 0.0
    for (i, s) in enumerate(stars)
        v2m, t3am, t3pm = cvis_to_obs(poly_to_cvis(x, s), data[i])
        d = data[i]
        c += weights[1]*sum(abs2, (v2m  .- d.v2)    ./ d.v2_err) +
             weights[2]*sum(abs2, (t3am .- d.t3amp) ./ d.t3amp_err) +
             weights[3]*sum(abs2, mod360(t3pm .- d.t3phi) ./ d.t3phi_err)
    end
    return c
end

"""
    fit_sphere_ld(data, tessels; kwargs...) -> NamedTuple

Fit a **round star with a uniform surface**, limb darkening only — the simplest parametric
model there is, and the right starting point for a non-rotating star (red supergiant,
Cepheid) before snapshot imaging.

Free parameters: `radius` (mas, POLAR radius = angular diameter / 2) and the limb-darkening
coefficients `ld1` (and `ld2` for `ldtype = 2`). `tpole` only sets the flux scale, which
cancels in the normalised visibilities, so it is held fixed.

# Why `weights = [1, 1, 0]` by default — T3φ is excluded

Any centrosymmetric model produces closure phases of exactly 0° or 180°. A spotted star does
not: RW Cep's T3φ has an rms of 118°. Including T3φ therefore does not improve the fit, it
drags the symmetric model toward a compromise that is *worse* in the amplitudes. Measured on
RW Cep, fitting a power-law LD disc:

| fitted on | D (mas) | α | χ²v2/n | χ²t3a/n |
|---|---|---|---|---|
| V² only | 3.00 | 1.80 | 8.6 | 6.1 |
| **V² + T3amp** | **2.93** | **1.58** | 8.9 | 5.4 |
| V² + T3amp + T3φ | 2.56 | 0.94 | 24.8 | 9.5 |

The large residual T3φ is not a failure — it *is* the asymmetry, and it is what the imaging
step then reconstructs. Pass `weights = [1,1,1]` only if you believe the star is symmetric.

# Accuracy

Validated against OITOOLS' analytic `visibility_ldpow` on RW Cep: at HEALPix nside 5 the
fitted diameter agrees to 0.003 % and α to 0.07 %; at nside 3, to 0.12 % and 1.2 %. The
discretisation-induced χ² is negligible at every resolution (rms ΔV² = 3.3e-4 at nside 3,
against typical error bars of 2e-3).

# Arguments

  - `ldtype` — 1 linear `1−u(1−μ)`, 2 quadratic, 3 Hestroffer `μ^α` (default 3)
  - `radius0`, `ld0` — starting guesses; `radius0` defaults to half the first-null estimate
  - `radius_bounds`, `ld_bounds` — box constraints
  - `fit_ld2` — also fit `ld2` (only meaningful for `ldtype = 2`)
  - `algorithm` — any derivative-free NLopt algorithm; `:LN_NELDERMEAD` (default),
    `:LN_BOBYQA` and `:LN_SBPLX` are all sensible here

Returns `(; radius, ld1, ld2, chi2, chi2_per_datum, params, status)`, where `params` is the
NamedTuple ready to hand to `create_star` / `create_star_multiepochs`.
"""
function fit_sphere_ld(data, tessels;
                       tepochs = zeros(length(data)),
                       ldtype::Int = 3, tpole::Real = 5000.0,
                       radius0 = nothing, ld0::Real = ldtype == 3 ? 1.0 : 0.3,
                       radius_bounds = (0.05, 20.0), ld_bounds = ldtype == 3 ? (0.0, 4.0) : (0.0, 1.0),
                       ld2_0::Real = 0.0, fit_ld2::Bool = false,
                       weights = [1.0, 1.0, 0.0],
                       method::Symbol = :neldermead,
                       algorithm::Symbol = :LN_NELDERMEAD,
                       maxeval::Int = 2000, xtol_rel::Real = 1e-6, ftol_rel::Real = 1e-8,
                       min_num_live_points::Int = 400, use_stepsampler::Bool = false,
                       verbose::Bool = false)
    data isa AbstractVector || (data = [data])
    # Starting radius from the shortest baseline if not supplied: V² ≈ 0.6 at B_min for a disc
    # comparable to λ/2B_min, which is good enough for bobyqa's trust region to find the rest.
    if radius0 === nothing
        B = sqrt.(data[1].uv[1,:].^2 .+ data[1].uv[2,:].^2)
        radius0 = 0.5 * 206264806.2 / (2 * maximum(B))       # half of λ/2B_max
        radius0 = clamp(radius0, radius_bounds[1]*2, radius_bounds[2]/2)
    end
    mk(r, l1, l2) = (surface_type = 0, radius = r, tpole = tpole, ldtype = ldtype,
                     ld1 = l1, ld2 = l2, inclination = 90.0, position_angle = 0.0,
                     rotation_period = 1.0)
    nfree = fit_ld2 ? 3 : 2
    x0 = fit_ld2 ? [radius0, ld0, ld2_0] : [radius0, ld0]
    lb = fit_ld2 ? [radius_bounds[1], ld_bounds[1], -1.0] : [radius_bounds[1], ld_bounds[1]]
    ub = fit_ld2 ? [radius_bounds[2], ld_bounds[2],  1.0] : [radius_bounds[2], ld_bounds[2]]
    nd = sum(d.nv2 + d.nt3amp + d.nt3phi for d in data)
    f(θ) = begin
        c = parametric_chi2(mk(θ[1], θ[2], fit_ld2 ? θ[3] : ld2_0), tessels, data, tepochs;
                            weights = weights, uniform = true)
        verbose && @info "fit_sphere_ld" radius=θ[1] ld1=θ[2] chi2=c
        isfinite(c) ? c : 1e30
    end
    if method === :ultranest
        nm = fit_ld2 ? ["radius", "ld1", "ld2"] : ["radius", "ld1"]
        un = _fit_ultranest(θ -> f(θ), nm, lb, ub;
                            min_num_live_points = min_num_live_points,
                            use_stepsampler = use_stepsampler, verb = verbose)
        m  = un.median
        l2 = fit_ld2 ? m[3] : ld2_0
        return (radius = m[1], ld1 = m[2], ld2 = l2,
                radius_err = (un.q84[1]-un.q16[1])/2, ld1_err = (un.q84[2]-un.q16[2])/2,
                chi2 = f(m), chi2_per_datum = f(m)/nd, params = mk(m[1], m[2], l2),
                logz = un.logz, logzerr = un.logzerr, samples = un.samples,
                q16 = un.q16, q84 = un.q84, status = :ultranest)
    end
    method === :neldermead || error("fit_sphere_ld: method must be :neldermead or :ultranest")
    opt = NLopt.Opt(algorithm, nfree)
    opt.lower_bounds  = lb
    opt.upper_bounds  = ub
    opt.min_objective = (θ, grad) -> f(θ)      # derivative-free: `grad` is never filled
    opt.maxeval       = maxeval
    opt.xtol_rel      = xtol_rel
    opt.ftol_rel      = ftol_rel
    # Nelder-Mead's initial simplex is built from the step sizes; leaving them at the default
    # makes it crawl when the parameters differ in scale (radius ~1 mas vs ld1 ~1).
    opt.initial_step  = fit_ld2 ? [0.1*radius0, 0.2, 0.2] : [0.1*radius0, 0.2]
    fopt, xopt, status = NLopt.optimize(opt, x0)
    l2 = fit_ld2 ? xopt[3] : ld2_0
    return (radius = xopt[1], ld1 = xopt[2], ld2 = l2, chi2 = fopt,
            chi2_per_datum = fopt/nd, params = mk(xopt[1], xopt[2], l2), status = status)
end

"""
    fit_ellipsoid_ld(data, tessels; kwargs...) -> NamedTuple

Fit an **oblate spheroid with a uniform surface**, limb darkening only — the next model up
from [`fit_sphere_ld`](@ref) when a star is measurably non-circular on the sky but you do not
want to commit to a gravity-darkening law.

Free parameters: equatorial radius `req` (mas), flattening `f = 1 − r_pol/r_eq`, inclination
and position angle (degrees), and `ld1`. The geometry uses `radius_x = radius_y = req` and
`radius_z = req·(1−f)`, i.e. rotational symmetry about the polar axis — the underlying
`surface_type = 1` allows three independent semi-axes, but fitting a triaxial ellipsoid to
interferometry is almost never justified.

!!! note "Uniform surface on purpose"
    The surface is held uniform (`uniform = true` in [`parametric_chi2`](@ref)), NOT
    gravity-darkened. `temperature_map_vonZeipel_ellipsoid` is a placeholder — `g ∝ 1/r²`
    with an unverified `rpole = radius_x` and a `# to check` in the source — so folding it
    into a shape fit would attribute its errors to the shape. If you want a physical
    gravity-darkened oblate star, use `surface_type = 2` (rapid rotator) and `fit_parametric`,
    which solves the Roche problem properly.

!!! warning "Degeneracies"
    Pole-on (`inc → 0`) the projection is circular whatever the flattening, so `f`, `inc` and
    `PA` all become unidentifiable together; and a nearly circular star leaves `PA` free. Check
    the fitted `f` against its uncertainty before believing an orientation.

Same `weights = [1, 1, 0]` default as [`fit_sphere_ld`](@ref), and for the same reason: a
uniform ellipsoid is still centrosymmetric in projection, so it predicts closure phases of
0° or 180° only.
"""
function fit_ellipsoid_ld(data, tessels;
                          tepochs = zeros(length(data)),
                          ldtype::Int = 3, tpole::Real = 5000.0,
                          req0 = nothing, f0::Real = 0.02, inc0 = nothing, pa0::Real = 45.0,
                          ld0 = nothing, seed_from_sphere::Bool = true,
                          req_bounds = (0.05, 20.0), f_bounds = (0.0, 0.6),
                          inc_bounds = (0.0, 90.0), pa_bounds = (0.0, 180.0),
                          ld_bounds = ldtype == 3 ? (0.0, 4.0) : (0.0, 1.0),
                          weights = [1.0, 1.0, 0.0],
                          method::Symbol = :neldermead,
                          algorithm::Symbol = :LN_NELDERMEAD,
                          maxeval::Int = 4000, xtol_rel::Real = 1e-6, ftol_rel::Real = 1e-8,
                          min_num_live_points::Int = 400, use_stepsampler::Bool = true,
                          verbose::Bool = false)
    data isa AbstractVector || (data = [data])
    # Seed from the SPHERE fit by default. Starting cold, Nelder-Mead walks into the pole-on
    # corner (inc → 0/180°, where the projection is circular and the flattening is
    # unidentifiable) and stalls with parameters pinned at bounds and a χ² WORSE than the
    # nested sphere — which is impossible for a converged fit, since the sphere is the f = 0
    # special case. Seeding puts the simplex in the right basin.
    chi2_sphere = Inf
    if seed_from_sphere && (req0 === nothing || ld0 === nothing)
        sph = fit_sphere_ld(data, tessels; tepochs = tepochs, ldtype = ldtype, tpole = tpole,
                            weights = weights, algorithm = algorithm)
        req0 === nothing && (req0 = sph.radius)
        ld0  === nothing && (ld0  = sph.ld1)
        chi2_sphere = sph.chi2
    end
    if req0 === nothing
        B = sqrt.(data[1].uv[1,:].^2 .+ data[1].uv[2,:].^2)
        req0 = clamp(0.5 * 206264806.2 / (2 * maximum(B)), req_bounds[1]*2, req_bounds[2]/2)
    end
    ld0 === nothing && (ld0 = ldtype == 3 ? 1.0 : 0.3)
    # Two exact two-fold degeneracies removed by the bounds above rather than left for the
    # optimiser to wander through: an oblate spheroid at inclination i looks identical at
    # 180° − i, and at position angle PA + 180°.
    inc0 === nothing && (inc0 = 60.0)
    mk(req, fl, inc, pa, l1) = (surface_type = 1, radius_x = req, radius_y = req,
                                radius_z = req*(1 - fl), tpole = tpole, ldtype = ldtype,
                                ld1 = l1, ld2 = 0.0, beta = 0.25,
                                inclination = inc, position_angle = pa, rotation_period = 1.0)
    nd = sum(d.nv2 + d.nt3amp + d.nt3phi for d in data)
    obj(θ) = begin
        c = parametric_chi2(mk(θ...), tessels, data, tepochs; weights = weights, uniform = true)
        verbose && @info "fit_ellipsoid_ld" req=θ[1] f=θ[2] inc=θ[3] PA=θ[4] ld1=θ[5] chi2=c
        isfinite(c) ? c : 1e30
    end
    lo5 = [req_bounds[1], f_bounds[1], inc_bounds[1], pa_bounds[1], ld_bounds[1]]
    hi5 = [req_bounds[2], f_bounds[2], inc_bounds[2], pa_bounds[2], ld_bounds[2]]
    if method === :ultranest
        # 5 parameters with two near-degenerate directions: the region sampler struggles, so
        # the slice step sampler is the default here (unlike the 2-parameter sphere).
        un = _fit_ultranest(θ -> obj(θ), ["req","flattening","inc","PA","ld1"], lo5, hi5;
                            min_num_live_points = min_num_live_points,
                            use_stepsampler = use_stepsampler, verb = verbose)
        m = un.median
        return (req = m[1], flattening = m[2], inclination = m[3], position_angle = m[4],
                ld1 = m[5], rpol = m[1]*(1-m[2]), chi2 = obj(m), chi2_per_datum = obj(m)/nd,
                chi2_sphere = chi2_sphere, params = mk(m...), logz = un.logz,
                logzerr = un.logzerr, samples = un.samples, q16 = un.q16, q84 = un.q84,
                status = :ultranest)
    end
    method === :neldermead || error("fit_ellipsoid_ld: method must be :neldermead or :ultranest")
    opt = NLopt.Opt(algorithm, 5)
    opt.lower_bounds  = lo5
    opt.upper_bounds  = hi5
    opt.min_objective = (θ, grad) -> obj(θ)
    opt.maxeval       = maxeval
    opt.xtol_rel      = xtol_rel
    opt.ftol_rel      = ftol_rel
    opt.initial_step  = [0.1*req0, 0.05, 10.0, 20.0, 0.2]
    fopt, xopt, status = NLopt.optimize(opt, [req0, f0, inc0, pa0, ld0])
    # The sphere is nested inside this model, so a converged ellipsoid can never fit worse.
    # If it does, the optimiser stalled — say so instead of returning a plausible-looking
    # parameter set.
    if isfinite(chi2_sphere) && fopt > chi2_sphere*(1 + 1e-8)
        @warn "fit_ellipsoid_ld did not converge: χ² ($(fopt)) exceeds the nested sphere fit " *
              "($(chi2_sphere)). Returned parameters are unreliable — check for values pinned " *
              "at bounds (inc, PA and flattening are jointly unidentifiable pole-on)." xopt
    end
    return (req = xopt[1], flattening = xopt[2], inclination = xopt[3], position_angle = xopt[4],
            ld1 = xopt[5], rpol = xopt[1]*(1-xopt[2]), chi2 = fopt, chi2_per_datum = fopt/nd,
            chi2_sphere = chi2_sphere, params = mk(xopt...), status = status)
end
