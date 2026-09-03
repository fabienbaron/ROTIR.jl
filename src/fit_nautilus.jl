# The Nautilus.jl backend for ROTIR's parametric fits.
#
# Included by ext/ROTIRNautilusExt.jl. Nothing here is public: `_fit_nautilus` is declared in
# the core package and reached through the same `method = :nautilus` switch that already
# selects `:neldermead` and `:ultranest`.
#
# WHY A THIRD SAMPLER, when `_fit_ultranest` already works.
#
# Nautilus.jl is a Julia port of the `nautilus` package (Lange 2023, MNRAS 525, 3181):
# importance nested sampling with neural-network-boosted bounds. Three properties matter for
# ROTIR specifically:
#
#   * **Pure Julia.** `_fit_ultranest` goes through PythonCall into a conda environment, and
#     every likelihood batch crosses that boundary. It also cannot survive into a
#     PackageCompiler build, and it is the reason ROTIR carries a Python dependency at all.
#   * **Importance nested sampling reuses every likelihood evaluation** instead of discarding
#     dead points, so it buys far more evidence precision per χ² call. In ROTIR a χ² call is
#     `parametric_chi2` — building the geometry at every epoch and evaluating the polygon
#     FT — which is orders of magnitude more expensive than the analytic model χ² OITOOLS
#     measured this on. The fewer calls a sampler needs, the more it wins here.
#   * **It accepts a single free parameter**, which the Python package does not; a one-radius
#     fit is a shape ROTIR routinely wants.
#
# OITOOLS measured the two head to head on its own model fits (see its
# src/fit_model_nautilus.jl header): on a 6-parameter fit, 6.4x FEWER likelihood calls and 27x
# tighter evidence, with the best-fit parameters agreeing to three or four decimals. That
# measurement is on a different likelihood and does not transfer automatically — but the
# direction of the call-count difference is a property of importance sampling, not of that
# problem, and the call count is what ROTIR pays for.
#
# The return shape is IDENTICAL to `_fit_ultranest`'s, so `fit_sphere_ld`, `fit_orbit` and the
# GUI can switch samplers without knowing which one ran.

"""
    _fit_nautilus(chi2_of, names, lo, hi; ...) -> NamedTuple

Nested-sampling fit of `chi2_of` over the box `[lo, hi]`, using Nautilus.

Returns the same fields as [`_fit_ultranest`](@ref): `median`, `q16`, `q84`, `mean`, `std`,
`logz`, `logzerr`, `samples`, plus `result` (the sampler itself, for anyone who wants
`Nautilus.eta` or the raw weights).

`n_live` is the live-point count and `n_eff` the effective sample size the run targets — the
two knobs worth turning. `n_eff` sets how thick the posterior is: equal-weight resampling
cannot produce more distinct points than the effective sample size, so a thin corner plot means
raising this rather than asking for more samples.

**The defaults are deliberately modest.** MEASURED on a 2-parameter sphere + limb-darkening fit
to λ And (nside 8, one epoch): at `n_eff = 10_000` the run took 836 s against 4.4 s for
Nelder–Mead on the same objective, and the two agreed to 1.4e-6 mas in radius and 1e-4 in
`ld1`. A ROTIR likelihood rebuilds the geometry at every epoch and evaluates the polygon
FT — it is orders of magnitude dearer than the analytic model χ² Nautilus is usually driven
with — so nested sampling here is bought for the POSTERIOR and the EVIDENCE, not for the point
estimate. If all you want is the best fit, `:neldermead` gets it in seconds.

Threaded by default. `parametric_chi2` builds its own geometry per call and shares nothing
mutable between calls, so the likelihood is safe to evaluate in parallel — which is exactly
the property `_fit_ultranest` does NOT have, because Python may only be called from a thread
holding the GIL.
"""
function _fit_nautilus(chi2_of, names, lo, hi;
                       n_live::Int = 500, n_eff::Real = 2000, f_live::Real = 0.01,
                       n_networks::Int = 4, discard_exploration::Bool = false,
                       n_like_max::Int = typemax(Int), timeout::Real = Inf,
                       equal_weight_boost::Real = 1, threaded::Bool = true,
                       seed::Union{Nothing,Integer} = nothing, verb::Bool = false)
    ndims = length(names)
    ndims == length(lo) == length(hi) ||
        throw(ArgumentError("names, lo and hi must have the same length"))
    all(isfinite, lo) && all(isfinite, hi) ||
        throw(ArgumentError("nautilus needs a FINITE box; an infinite bound has no prior " *
                            "volume to sample. Narrow it, or use :neldermead."))

    L, U = collect(Float64, lo), collect(Float64, hi)
    # The unit cube to physical parameters. Nautilus hands the transform a unit vector and the
    # likelihood the PHYSICAL parameters that come out of it — the same convention UltraNest
    # uses, and applying the transform a second time inside the likelihood is the mistake that
    # makes a sampler converge confidently on nonsense.
    prior_transform(u) = L .+ (U .- L) .* u

    # Track the best point seen. Nautilus reports the posterior, not the maximum, and a fit
    # wants both: the median for the parameter estimate, the argmax for the χ² that goes in a
    # table beside it.
    best = Ref((-Inf, Float64[]))
    lk = ReentrantLock()
    function loglikelihood(θ)
        c = chi2_of(collect(Float64, θ))
        ll = isfinite(c) ? -0.5 * Float64(c) : -1e30
        if ll > best[][1]
            # Under `threaded = true` several evaluations land at once, so the comparison and
            # the store have to be one operation. The unlocked pre-check keeps the lock off the
            # common path, where the point is not an improvement.
            lock(lk) do
                ll > best[][1] && (best[] = (ll, collect(Float64, θ)))
            end
        end
        return ll
    end

    smplr = Nautilus.Sampler(prior_transform, loglikelihood;
                             n_dim = ndims, n_live = n_live, n_networks = n_networks,
                             threaded = threaded,
                             seed = seed === nothing ? rand(UInt32) : seed)
    finished = Nautilus.run!(smplr; f_live = f_live, n_eff = n_eff, n_like_max = n_like_max,
                             timeout = timeout, discard_exploration = discard_exploration,
                             verbose = verb)
    finished || @warn "nautilus stopped on n_like_max or timeout before converging; the " *
                      "evidence and posterior are usable but less accurate than requested."

    logz = Nautilus.log_z(smplr)
    neff = Nautilus.n_eff(smplr)
    # Nautilus quotes no evidence error of its own; 1/sqrt(N_eff) is the standard estimate and
    # is what OITOOLS' backend reports, so the two are comparable.
    logzerr = neff > 0 ? 1 / sqrt(neff) : NaN

    # Equal-weight samples, so `samples` means the same thing as UltraNest's and a corner plot
    # can take it directly without per-sample weights.
    S, _, _ = Nautilus.posterior(smplr; equal_weight = true,
                                 equal_weight_boost = equal_weight_boost)
    isempty(S) && error("nautilus returned no posterior samples; the run failed")

    q(pr) = [Statistics.quantile(view(S, :, j), pr) for j in 1:ndims]
    if verb
        Printf.@printf("logZ = %.4f ± %.4f   (%d samples, %d calls, N_eff = %.0f)\n",
                       logz, logzerr, size(S, 1), smplr.n_like, neff)
    end
    return (median = q(0.5), q16 = q(0.16), q84 = q(0.84),
            mean = vec(Statistics.mean(S, dims = 1)),
            std = vec(Statistics.std(S, dims = 1)),
            logz = logz, logzerr = logzerr, samples = S,
            best = best[][2], best_chi2 = -2 * best[][1], result = smplr)
end
