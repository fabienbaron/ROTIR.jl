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
# Against: no gradients are used, cost grows steeply with dimension, and PythonCall is not
# thread-safe — which is why `threaded` reaches `_fit_nautilus` only, never UltraNest. Safe to
# call from a multithreaded session either way.
"""
    parametric_posterior_spec(data_epochs, tessels, tepochs, base_params; kwargs...)
        -> (logπ, lower, upper, free_idx, layout)

The rapid rotator's log-posterior together with everything a sampler needs around it: box
bounds, the resolved free-parameter indices, and the [`ParametricLayout`](@ref) they all
refer to.

Factored out because `_fit_hmc` and `_fit_pigeons` carried a byte-identical copy of this and
`fit_parametric_ultranest` a third variant, so every new physical parameter had to be added
to all three — and a missed one does not fail loudly, it just indexes a vector of the wrong
length. Lives in the core rather than in an extension for the same reason
[`_box_transform`](@ref) does: sibling package extensions cannot import from one another.

`provider` and `layout` are what let a sampler reach the parameters the legacy 7-vector has
no room for — `d` above all, which is what makes a model-atmosphere `logg` fittable. See
[`build_parametric_logπ`](@ref).
"""
function parametric_posterior_spec(data_epochs, tessels, tepochs, base_params;
                                   free = nothing, intensity_model::Symbol = :linear,
                                   band = nothing, κ = 50, GM = 1,
                                   tpole_free::Bool = false, gravity_law = nothing,
                                   logprior = nothing, provider = nothing, layout = nothing,
                                   weights = (1, 1, 1))
    L = layout === nothing ? parametric_layout(; tpole_free = tpole_free) : layout
    lp = build_parametric_logπ(data_epochs, tessels, tepochs, base_params;
                               intensity_model = intensity_model, band = band,
                               κ = κ, GM = GM, tpole_free = tpole_free,
                               gravity_law = gravity_law, logprior = logprior,
                               provider = provider, layout = L, weights = weights)
    return lp, copy(L.lower), copy(L.upper), layout_free_indices(L, free), L
end

"""
    _box_transform(lo, hi) -> (to_θ, to_z, logabsjac)

An unconstrained reparameterisation of a box, and its log-Jacobian.

Every free parameter is mapped from ℝ so NUTS cannot step outside its bounds. That is not
convenience: a leapfrog trajectory routinely overshoots, and a sampler that has to reject those
proposals loses most of its steps at exactly the point where the posterior is informative.

Three cases, because the demo's hand-written pair is the one-star special case of them:

  * both bounds finite → a scaled logistic, so `θ` lives strictly inside `(lo, hi)`;
  * `lo` finite only    → `lo + exp(z)`, the positive-parameter case (`rpole`);
  * neither finite      → the identity, for an angle.

The log-Jacobian is the SUM over parameters, and it must be added to the log-posterior or the
sampler is walking a different distribution from the one it reports.
"""
function _box_transform(lo::AbstractVector, hi::AbstractVector)
    n = length(lo)
    kinds = [isfinite(lo[j]) && isfinite(hi[j]) ? :box :
             isfinite(lo[j])                    ? :pos : :free for j in 1:n]
    logistic(x) = 1 / (1 + exp(-x))
    function to_θ(z)
        return [kinds[j] === :box  ? lo[j] + (hi[j] - lo[j]) * logistic(z[j]) :
                kinds[j] === :pos  ? lo[j] + exp(z[j]) : z[j] for j in 1:n]
    end
    function to_z(θ)
        return [kinds[j] === :box  ? log((θ[j] - lo[j]) / max(hi[j] - θ[j], eps())) :
                kinds[j] === :pos  ? log(max(θ[j] - lo[j], eps())) : θ[j] for j in 1:n]
    end
    function logabsjac(z)
        s = zero(eltype(z))
        @inbounds for j in 1:n
            if kinds[j] === :box
                u = logistic(z[j])
                s += log(hi[j] - lo[j]) + log(u) + log1p(-u)
            elseif kinds[j] === :pos
                s += z[j]
            end
        end
        return s
    end
    return to_θ, to_z, logabsjac
end

"""
    _fit_nested(method, chi2_of, names, lo, hi; kwargs...) -> NamedTuple

Dispatch to whichever nested sampler `method` names, and give one useful error when it is not
available.

Both backends return the same fields (`median`, `q16`, `q84`, `logz`, `logzerr`, `samples`),
which is what lets every caller treat them interchangeably. UltraNest-only keywords
(`use_stepsampler`, `min_num_live_points`) and Nautilus-only ones (`n_live`, `n_eff`,
`threaded`) are filtered here rather than at each call site, so a caller can pass whatever it
has.

`threaded` IS FORWARDED, and it was not: the keyword existed on `_fit_nautilus` and this
function dropped it on the floor along with the rest of `kwargs`, so a caller that asked for a
serial likelihood got a threaded one and no warning. That is not a preference — a `chi2_of`
that is not re-entrant THROWS on every concurrent evaluation, and a caller whose objective
swallows exceptions then gets a converged-looking run built entirely out of sentinel values.
ROTIR's own kernels are all re-entrant now (`:t3`'s mutable scratch comes from
[`t3_checkout!`](@ref), not from the shared plan), so the default is `true`; an objective
supplied by a caller is the case `threaded = false` exists for.
"""
function _fit_nested(method::Symbol, chi2_of, names, lo, hi; verb::Bool = false,
                     min_num_live_points::Int = 400, use_stepsampler::Bool = false,
                     n_live::Int = 1000, n_eff::Real = 10_000, threaded::Bool = true,
                     kwargs...)
    if method === :ultranest
        ultranest_available() || error(
            "method = :ultranest needs PythonCall loaded: add `using PythonCall` to this " *
            "session. It is a weak dependency, because loading it costs every GUI session " *
            "1.2 s of recompilation whether or not anything samples. `method = :nautilus` " *
            "is the pure-Julia sampler and needs no Python at all.")
        return _fit_ultranest(chi2_of, names, lo, hi; verb = verb,
                              min_num_live_points = min_num_live_points,
                              use_stepsampler = use_stepsampler)
    elseif method === :nautilus
        nautilus_available() || error(
            "method = :nautilus needs Nautilus.jl loaded: add `using Nautilus` to this " *
            "session. It is a weak dependency, so ROTIR does not pull it in on its own.")
        return _fit_nautilus(chi2_of, names, lo, hi; verb = verb, n_live = n_live,
                             n_eff = n_eff, threaded = threaded)
    end
    error("unknown nested sampler :$method (expected :ultranest or :nautilus)")
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
                         weights = (1.0, 1.0, 1.0), uniform::Bool = false)
    stars = create_star_multiepochs(tessels, params, tepochs)
    # NO `setup_oi!`. This built the dense (nuv × npix) `polyft` for every evaluation, used it
    # for a single matvec and threw it away — 5 to 10 GiB of temporaries to produce a 189 MiB
    # matrix that is read once. The fused kernel computes the same visibilities without ever
    # forming it, and it was already in the package: `parametric_gradient.jl` and
    # `shape_gradient.jl` have used it all along, and only this function did not.
    #
    # MEASURED on polaris at level 3, one epoch: 784 ms → 114 ms with the scalar kernel and
    # 5.1 ms with `:turbo`. This is the objective every optimiser and every sampler calls, so
    # it is the cost of a fit.
    #
    # Imaging is the opposite case and rightly keeps `setup_oi!`: there the geometry is fixed
    # and the map varies over hundreds of iterations, so the matrix is worth building once.
    x = uniform ? fill(eltype(stars[1])(params.tpole), stars[1].npix) :
                  parametric_temperature_map(params, stars[1])
    # THE ACCUMULATOR FOLLOWS THE MODEL. `c = 0.0` and a `Vector{Float64}` of weights made
    # this Float64 whatever the mesh was, so every term of a Float32 model promoted and the
    # function returned Float64 from a Float32 star. The weights are a TUPLE now as well: a
    # `Vector{Float64}` default allocates on every call and forces the same promotion through
    # `weights[k]`, while a tuple of literals costs nothing and converts at the multiply.
    #
    # Small beside `fused_cvis`, which is where this function's time actually goes — but a χ²
    # whose type depends on a default keyword rather than on the model is a trap, and the
    # comparison against `-2 logπ` is cleaner when both sides are the mesh's own type.
    S = eltype(x)
    c = zero(S)
    for (i, s) in enumerate(stars)
        v2m, t3am, t3pm = cvis_to_obs(fused_cvis(x, s, data[i]), data[i])
        d = data[i]
        c += S(weights[1])*sum(abs2, (v2m  .- d.v2)    ./ d.v2_err) +
             S(weights[2])*sum(abs2, (t3am .- d.t3amp) ./ d.t3amp_err) +
             S(weights[3])*sum(abs2, mod360(t3pm .- d.t3phi) ./ d.t3phi_err)
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
                       # TRUE, and it was `false` for one release while the default `:t3`
                       # kernel was not re-entrant. Nautilus parallelises its likelihood calls
                       # and ROTIR's χ² threads internally, so this NESTS threading inside
                       # threading — which is fine now that `:t3`'s mutable scratch comes from
                       # a pool (`t3_checkout!`) rather than from the globally cached plan,
                       # and the schedules are `:dynamic`. Nesting is also worth having: 8
                       # concurrent transforms run 1.96x faster than the same 8 in a loop, on
                       # top of the 4.8-7.2x (HEALPix 3-5, 16 threads) that the inner
                       # threading already gives.
                       #
                       # Pass `threaded = false` for an objective of your own that is not
                       # re-entrant. The failure it protects against is NOT a crash: a
                       # `chi2_of` that throws under concurrency and swallows it returns the
                       # 1e30 sentinel, and the sampler reports the FLAT PRIOR as a converged
                       # posterior.
                       threaded::Bool = true,
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
    if method in (:ultranest, :nautilus)
        nm = fit_ld2 ? ["radius", "ld1", "ld2"] : ["radius", "ld1"]
        un = _fit_nested(method, θ -> f(θ), nm, lb, ub;
                         min_num_live_points = min_num_live_points, threaded = threaded,
                         use_stepsampler = use_stepsampler, verb = verbose)
        m  = un.median
        l2 = fit_ld2 ? m[3] : ld2_0
        return (radius = m[1], ld1 = m[2], ld2 = l2,
                radius_err = (un.q84[1]-un.q16[1])/2, ld1_err = (un.q84[2]-un.q16[2])/2,
                chi2 = f(m), chi2_per_datum = f(m)/nd, params = mk(m[1], m[2], l2),
                logz = un.logz, logzerr = un.logzerr, samples = un.samples,
                q16 = un.q16, q84 = un.q84, status = method)
    end
    method === :neldermead ||
        error("fit_sphere_ld: method must be :neldermead, :ultranest or :nautilus")
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

Fit a **triaxial ellipsoid with a uniform surface**, limb darkening only — the next model up
from [`fit_sphere_ld`](@ref) when a star is measurably non-circular on the sky but you do not
want to commit to a gravity-darkening law.

Free parameters: the three body-frame semi-axes `radius_x`, `radius_y`, `radius_z` (mas),
inclination and position angle (degrees), and `ld1`. That is the model's own parameterization —
`surface_type = 1` carries three independent radii and
[`ellipsoid_param_names`](@ref) lists them — so nothing here is a reparameterization of it.

!!! warning "This function used to fit an oblate spheroid, and that was wrong"
    Until this change it imposed `radius_x = radius_y = req` and `radius_z = req*(1 - f)`,
    fitting a 5-parameter OBLATE figure: a two-dimensional slice of the model's three-radius
    shape space, with the symmetry axis pinned to the spin axis. A tidally distorted star —
    the usual reason a giant is non-circular — is elongated TOWARDS ITS COMPANION, i.e.
    `radius_x > radius_y ≈ radius_z`, which that slice cannot represent at any parameter value.
    Asked for such a figure it drives `f` to its bound and reports a shape that is not the
    star's. The old `req0`, `f0`, `req_bounds` and `f_bounds` keywords are refused rather than
    silently reinterpreted.

!!! note "The spin phase is now load-bearing"
    An oblate figure is axisymmetric, so `rotation_period` and the epoch times did not affect
    it and the old code hardcoded `rotation_period = 1.0`. A triaxial figure presents a
    DIFFERENT projected ellipse at every rotation phase — which is exactly what makes the three
    radii separately identifiable — so the period is required whenever the epochs differ, and
    there is no default. `t0` sets the phase zero (`rotate_star` uses `ψ = 2π(t - t0)/P` and the
    schema has no phase-offset parameter). For a synchronised binary, set `t0` to the
    conjunction epoch so that phase 0 is the companion-facing orientation; otherwise the fitted
    axes are labelled against an arbitrary zero and only their RATIOS are meaningful.

    From a single rotation phase the projected outline is one ellipse — three numbers — against
    five shape-and-orientation parameters, so the radii are NOT separately identifiable. The
    fit still runs (the best-fitting projected ellipse is well defined) and warns.

!!! note "Uniform surface on purpose"
    The surface is held uniform (`uniform = true` in [`parametric_chi2`](@ref)), NOT
    gravity-darkened. `temperature_map_vonZeipel_ellipsoid` is a placeholder — `g ∝ 1/r²`
    with an unverified `rpole = radius_x` and a `# to check` in the source — so folding it
    into a shape fit would attribute its errors to the shape. If you want a physical
    gravity-darkened oblate star, use `surface_type = 2` (rapid rotator) and `fit_parametric`,
    which solves the Roche problem properly; for a tidally distorted star in a close binary the
    physically correct figure is `surface_type = 3` (Roche lobe), which has ONE shape parameter
    rather than three free radii.

!!! warning "Degeneracies"
    Pole-on (`inc → 0`) the projection is circular whatever the shape, so the radii, `inc` and
    `PA` become jointly unidentifiable; and a nearly spherical star leaves `PA` free. `inc` is
    bounded to [0°, 90°] and `PA` to [0°, 180°] because a figure at `180° − inc` and at
    `PA + 180°` are all but indistinguishable — all but, rather than exactly, because for a
    triaxial body flipping the inclination reverses the direction the long axis sweeps with
    phase, which a long enough phase baseline can in principle tell apart. Check every fitted
    radius against its uncertainty before believing an orientation.

`restarts` (default 6) sweeps the starting orientation and keeps the best χ², because this
objective has several basins and the bad ones look converged — see the comment at the optimiser
for the ζ And measurement that sets the default. Pass `restarts = 1` to reproduce a single run
from `inc0`/`pa0`.

A uniform ellipsoid is still centrosymmetric in projection, so it predicts closure phases of
0° or 180° only — hence the same `weights = [1, 1, 0]` default as [`fit_sphere_ld`](@ref). On a
spotted star the T3φ residual is the spot signal, and it is usually far larger than the shape
signal: ζ And's figure is 6 % distorted (0.14 mas, a third of a resolution element at CHARA/H)
against a ~20 % surface brightness contrast. A uniform-surface shape fit on such a star is
driven by the spots, which is an argument for a joint model + imaging fit rather than for
tightening the bounds.
"""
function fit_ellipsoid_ld(data, tessels;
                          tepochs = zeros(length(data)),
                          ldtype::Int = 3, tpole::Real = 5000.0,
                          r0 = nothing, inc0 = nothing, pa0::Real = 45.0,
                          ld0 = nothing, seed_from_sphere::Bool = true,
                          rotation_period = nothing, t0::Real = 0.0,
                          r_bounds = (0.05, 20.0),
                          inc_bounds = (0.0, 90.0), pa_bounds = (0.0, 180.0),
                          ld_bounds = ldtype == 3 ? (0.0, 4.0) : (0.0, 1.0),
                          weights = [1.0, 1.0, 0.0],
                          method::Symbol = :neldermead,
                          algorithm::Symbol = :LN_NELDERMEAD,
                          restarts::Int = 6,
                          maxeval::Int = 6000, xtol_rel::Real = 1e-6, ftol_rel::Real = 1e-8,
                          min_num_live_points::Int = 400, use_stepsampler::Bool = true,
                          # TRUE — see the note on `fit_sphere_ld`. Nested threading is
                          # safe since `:t3`'s scratch moved into a pool; before that it was
                          # `false` because the sampler's concurrent likelihood calls shared
                          # one plan's grids.
                          threaded::Bool = true,
                          verbose::Bool = false,
                          # REFUSED, not reinterpreted: these named the oblate parameterization
                          # this function used to impose. Accepting them silently would fit a
                          # different model than the caller's script says it fits.
                          req0 = nothing, f0 = nothing, req_bounds = nothing,
                          f_bounds = nothing)
    for (k, v) in (("req0", req0), ("f0", f0), ("req_bounds", req_bounds),
                   ("f_bounds", f_bounds))
        v === nothing || error(
            "fit_ellipsoid_ld: `$k` belonged to the oblate parameterization (radius_x = " *
            "radius_y = req, radius_z = req*(1-f)) that this function no longer imposes. It " *
            "now fits the model's own three semi-axes. Use `r0` (a scalar or a 3-tuple of " *
            "starting radii in mas) and `r_bounds`, which apply to all three.")
    end
    data isa AbstractVector || (data = [data])
    ts = collect(float.(tepochs)) .- float(t0)
    # THE PERIOD HAS NO DEFAULT when it matters. A triaxial figure's sky projection depends on
    # the rotation phase, so a wrong period is a wrong model rather than a wrong scale — and
    # the old hardcoded `rotation_period = 1.0` would silently turn a 20-day baseline into 20
    # full turns. When every epoch is the same instant the phase cannot matter, and only then
    # is a placeholder safe.
    P = rotation_period
    if P === nothing
        if length(unique(ts)) <= 1
            P = 1.0
        else
            error("fit_ellipsoid_ld: `rotation_period` is required for multi-epoch data. A " *
                  "triaxial figure presents a different projected ellipse at every rotation " *
                  "phase — that is what identifies the three radii — so the period sets the " *
                  "model, not just a time unit. Pass the photometric or orbital period in the " *
                  "same units as `tepochs` (days), and `t0` for the phase zero.")
        end
    end
    P > 0 || error("fit_ellipsoid_ld: `rotation_period` must be positive, got $P")
    phases = sort(unique(round.(mod.(ts ./ P, 1.0); digits = 6)))
    if length(phases) <= 1
        @warn "fit_ellipsoid_ld: all epochs fall at the same rotation phase, so the data " *
              "constrain ONE projected ellipse — three numbers against five shape and " *
              "orientation parameters. The returned radii are one point on a degenerate " *
              "surface, not a measurement of the figure." phase = only(phases)
    end
    # Seed from the SPHERE fit by default. Starting cold, Nelder-Mead walks into the pole-on
    # corner (inc → 0/180°, where the projection is circular and the shape is unidentifiable)
    # and stalls with parameters pinned at bounds and a χ² WORSE than the nested sphere —
    # which is impossible for a converged fit, since the sphere is the equal-radii special
    # case. Seeding puts the simplex in the right basin.
    chi2_sphere = Inf
    if seed_from_sphere && (r0 === nothing || ld0 === nothing)
        sph = fit_sphere_ld(data, tessels; tepochs = tepochs, ldtype = ldtype, tpole = tpole,
                            weights = weights, algorithm = algorithm)
        r0  === nothing && (r0  = sph.radius)
        ld0 === nothing && (ld0 = sph.ld1)
        chi2_sphere = sph.chi2
    end
    if r0 === nothing
        B = sqrt.(data[1].uv[1,:].^2 .+ data[1].uv[2,:].^2)
        r0 = clamp(0.5 * 206264806.2 / (2 * maximum(B)), r_bounds[1]*2, r_bounds[2]/2)
    end
    # A SCALAR START IS BROKEN DELIBERATELY. At exactly equal radii the body is a sphere, which
    # has no orientation at all: `inc` and `PA` sit on a flat direction and a derivative-free
    # simplex has nothing to descend. 2 % is far below any shape these data could confuse it
    # with and far above the mesh's own asphericity.
    rx0, ry0, rz0 = r0 isa Number ? (1.02*r0, float(r0), 0.98*r0) : float.(Tuple(r0))
    ld0 === nothing && (ld0 = ldtype == 3 ? 1.0 : 0.3)
    inc0 === nothing && (inc0 = 60.0)
    mk(rx, ry, rz, inc, pa, l1) = (surface_type = 1, radius_x = rx, radius_y = ry,
                                   radius_z = rz, tpole = tpole, ldtype = ldtype,
                                   ld1 = l1, ld2 = 0.0, beta = 0.25,
                                   inclination = inc, position_angle = pa,
                                   rotation_period = P)
    nd = sum(d.nv2 + d.nt3amp + d.nt3phi for d in data)
    obj(θ) = begin
        c = parametric_chi2(mk(θ...), tessels, data, ts; weights = weights, uniform = true)
        verbose && @info "fit_ellipsoid_ld" rx=θ[1] ry=θ[2] rz=θ[3] inc=θ[4] PA=θ[5] ld1=θ[6] chi2=c
        isfinite(c) ? c : 1e30
    end
    lo6 = [r_bounds[1], r_bounds[1], r_bounds[1], inc_bounds[1], pa_bounds[1], ld_bounds[1]]
    hi6 = [r_bounds[2], r_bounds[2], r_bounds[2], inc_bounds[2], pa_bounds[2], ld_bounds[2]]
    names6 = ["radius_x", "radius_y", "radius_z", "inc", "PA", "ld1"]
    # The shape numbers worth reading, as opposed to the three radii separately: how elongated
    # the equator is (tidal) and how flattened the figure is pole to equator (rotational).
    # Reported alongside the radii because those two have different physical causes and the
    # three radii on their own do not separate them.
    shape(x) = (equatorial_elongation = x[1]/x[2],
                polar_flattening = 1 - 2x[3]/(x[1] + x[2]),
                mean_radius = cbrt(x[1]*x[2]*x[3]))
    if method in (:ultranest, :nautilus)
        # 6 parameters with several near-degenerate directions: the region sampler struggles,
        # so the slice step sampler is the default here (unlike the 2-parameter sphere).
        # Nautilus has no equivalent knob — its neural bounds are what it uses instead — so the
        # keyword is simply not passed on to it.
        un = _fit_nested(method, θ -> obj(θ), names6, lo6, hi6;
                         min_num_live_points = min_num_live_points, threaded = threaded,
                         use_stepsampler = use_stepsampler, verb = verbose)
        m = un.median
        return (radius_x = m[1], radius_y = m[2], radius_z = m[3], inclination = m[4],
                position_angle = m[5], ld1 = m[6], shape(m)...,
                rotation_period = P, phases = phases,
                chi2 = obj(m), chi2_per_datum = obj(m)/nd, chi2_sphere = chi2_sphere,
                params = mk(m...), logz = un.logz, logzerr = un.logzerr,
                samples = un.samples, q16 = un.q16, q84 = un.q84, status = method)
    end
    method === :neldermead ||
        error("fit_ellipsoid_ld: method must be :neldermead, :ultranest or :nautilus")
    # RESTARTS ARE NOT A LUXURY HERE, and the measurement is on real data: ζ And at
    # HEALPix 3, twelve starts over `inc0` x `pa0`, gave THREE basins — the best at
    # chi2/n = 3.626 (equatorial/polar = 1.070, agreeing with Roettenbacher et al.'s 1.06),
    # a pole-on corner at 4.325 with `radius_z` pinned to its ceiling, and an edge-on one at
    # 5.561. Four of the twelve, including the default start, landed in the pole-on corner and
    # reported 0.41 against the 1.07 the other basin and the literature agree on. A single
    # Nelder-Mead run on this objective is a coin toss,
    # and the bad basin LOOKS converged — hence the deterministic sweep, and hence the
    # warnings above, which fire on exactly those corners.
    #
    # The sweep is over orientation only. The radii and `ld1` are seeded from the sphere, which
    # is a good start for all of them; it is `inc`/`PA` that select the basin.
    incs = restarts <= 1 ? [inc0] :
           collect(range(inc_bounds[1] + 0.2*(inc_bounds[2]-inc_bounds[1]),
                         inc_bounds[1] + 0.8*(inc_bounds[2]-inc_bounds[1]),
                         length = max(2, round(Int, sqrt(restarts)))))
    pas  = restarts <= 1 ? [pa0] :
           collect(range(pa_bounds[1], pa_bounds[2],
                         length = max(2, cld(restarts, length(incs))) + 1))[1:end-1]
    fopt = Inf; xopt = [rx0, ry0, rz0, inc0, pa0, ld0]; status = :none
    for i0 in incs, p0 in pas
        opt = NLopt.Opt(algorithm, 6)
        opt.lower_bounds  = lo6
        opt.upper_bounds  = hi6
        opt.min_objective = (θ, grad) -> obj(θ)
        opt.maxeval       = maxeval
        opt.xtol_rel      = xtol_rel
        opt.ftol_rel      = ftol_rel
        opt.initial_step  = [0.1*rx0, 0.1*ry0, 0.1*rz0, 10.0, 20.0, 0.2]
        f, x, st = NLopt.optimize(opt, [rx0, ry0, rz0, clamp(i0, inc_bounds...),
                                        clamp(p0, pa_bounds...), ld0])
        verbose && @info "fit_ellipsoid_ld restart" inc0=i0 pa0=p0 chi2=f
        f < fopt && ((fopt, xopt, status) = (f, x, st))
    end
    # The sphere is nested inside this model, so a converged ellipsoid can never fit worse.
    # If it does, the optimiser stalled — say so instead of returning a plausible-looking
    # parameter set.
    if isfinite(chi2_sphere) && fopt > chi2_sphere*(1 + 1e-8)
        @warn "fit_ellipsoid_ld did not converge: χ² ($(fopt)) exceeds the nested sphere fit " *
              "($(chi2_sphere)). Returned parameters are unreliable — check for values pinned " *
              "at bounds (the radii, inc and PA are jointly unidentifiable pole-on)." xopt
    end
    # AT-BOUND, MEASURED AGAINST THE RANGE, not multiplicatively. `x <= l*(1 + 1e-6)` can
    # never fire on a bound of 0 — which is `inc`'s and `PA`'s — and `x >= h*(1 - 1e-6)` missed
    # a radius sitting at 99.94 % of its ceiling on real data. A fraction of the interval
    # catches both.
    atbound = [n for (n, x, l, h) in zip(names6, xopt, lo6, hi6)
               if x <= l + 1e-3*(h - l) || x >= h - 1e-3*(h - l)]
    isempty(atbound) || @warn "fit_ellipsoid_ld: parameters pinned at a bound — the figure " *
        "the data want is outside the box, or outside this model." atbound
    # POLE-ON IS THE DEGENERATE CORNER and it is reached long before `inc` hits zero: the
    # projected outline of ANY ellipsoid is circular at inc = 0, so the radii stop being
    # identifiable smoothly as inc falls. A fit that walks here has usually found a way to
    # absorb something that is not shape — surface structure, most often — rather than a
    # nearly pole-on star.
    xopt[4] < 15.0 && @warn "fit_ellipsoid_ld: inclination fitted to $(round(xopt[4], digits=2))°, " *
        "where the projection is nearly circular and the three radii are weakly identifiable. " *
        "Treat the figure as unmeasured unless an independent inclination says otherwise."
    return (radius_x = xopt[1], radius_y = xopt[2], radius_z = xopt[3], inclination = xopt[4],
            position_angle = xopt[5], ld1 = xopt[6], shape(xopt)...,
            rotation_period = P, phases = phases,
            chi2 = fopt, chi2_per_datum = fopt/nd, chi2_sphere = chi2_sphere,
            params = mk(xopt...), status = status)
end
