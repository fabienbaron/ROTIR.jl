# ultranest.jl
#
# Nested sampling of the parametric model through the Python `ultranest` package.
#
# Why this lives in core rather than in an extension: nested sampling never asks for a
# gradient, so it needs no AD backend, and PyCall is already a ROTIR dependency. The
# Python side is reached exactly as OITOOLS' `fit_model_ultranest` does it.
#
# When to reach for it instead of the optimiser + bootstrap:
#   - the χ² is multimodal (ρ Cas has several minima between 2.2 and 3.7 mas diameter) and
#     nested sampling explores all of them without a starting point;
#   - you want the evidence log(Z) to compare models — sphere vs oblate vs oblate with
#     gravity darkening — rather than just fit one;
#   - you would rather spend evaluations than tune warm starts.
# Against: it ignores the gradients the model can supply (a gradient costs ~3.6× a
# likelihood here), and the cost grows steeply with the number of free parameters.

using PyCall

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

    # UltraNest works on the unit hypercube; the transform makes it a uniform box prior.
    transform(cube) = lo .+ (hi .- lo) .* collect(Float64, cube)

    function loglike(cube)
        θtry = copy(θ)
        θtry[idx] .= T.(transform(cube))
        v = logπ(θtry)
        # NaN/Inf (e.g. a degenerate geometry) must not reach the sampler
        return isfinite(v) ? Float64(v) : -1e30
    end

    ultranest = pyimport("ultranest")
    sampler = ultranest.ReactiveNestedSampler(names, loglike, transform;
                                              log_dir=log_dir, resume=resume)
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

    S = Array{Float64}(res["samples"])                    # nsamples × nfree
    q(p) = [quantile(view(S, :, j), p) for j in eachindex(idx)]
    out = (θ_mean   = vec(mean(S, dims=1)),
           θ_median = q(0.5),
           θ_std    = vec(std(S, dims=1)),
           q16      = q(0.16),
           q84      = q(0.84),
           logz     = res["logz"],
           logzerr  = res["logzerr"],
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
