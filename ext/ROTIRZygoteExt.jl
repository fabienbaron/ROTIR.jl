module ROTIRZygoteExt

# Point estimator for the parametric von Zeipel model.
#
# ROTIR core defines the model and its ChainRulesCore rrules but takes no AD dependency,
# so this is where an actual gradient — and therefore an actual fit — becomes available.
# `using Zygote` alongside ROTIR is all it takes.

using ROTIR
using Zygote
using OptimPackNextGen
using LinearAlgebra
using Printf

import ROTIR: fit_parametric, bootstrap_parametric, default_parametric_bounds,
              parametric_param_names, parametric_free_indices, build_parametric_logπ

"""
    fit_parametric(data_epochs, tessels, tepochs, base_params; kwargs...) -> (θ̂, chi2r, info)

See the docstring in ROTIR core. Minimizes `-logπ` with `OptimPackNextGen.vmlmb` under box
bounds, using a Zygote gradient.

Bounds replace the exp/sigmoid reparametrization of the NUTS demo: a sampler needs an
unconstrained space, a bounded optimizer does not.

# Keywords
- `θ0` — starting point; required (full length, including the frozen entries)
- `free` — which entries to fit: `nothing` (all), indices, or names such as
  `["rpole", "ld1"]`. Frozen entries keep their `θ0` value. See
  [`parametric_free_indices`](@ref)
- `lb`, `ub` — box bounds (default [`default_parametric_bounds`](@ref))
- `intensity_model`, `band`, `κ`, `GM`, `tpole_free` — passed to `build_parametric_logπ`
- `logprior` — optional `θ -> log p(θ)` added to the objective. The reported χ²ᵣ is always
  the likelihood part alone
- `maxiter`, `gtol`, `mem`, `verb` — optimizer settings

Returns `(θ̂, chi2r, info)` — `θ̂` at full length including frozen entries — where `info`
holds `evaluations`, `gnorm`, `gratio` (‖g‖ relative to the start), `ndof`, `free` and
`hit_maxiter`. There is deliberately no `converged` flag: vmlmb prints its stop reason
without returning it, and ‖g‖ is not a usable substitute (a highly resolved star has a
χ² steep enough that ‖g‖ stays O(10) at a minimum whose σ is ~1e-3 mas).
"""
function fit_parametric(data_epochs::AbstractVector, tessels, tepochs, base_params;
    θ0,
    free                     = nothing,
    lb                       = nothing,
    ub                       = nothing,
    intensity_model ::Symbol = :linear,
    band                     = nothing,
    κ                        = 50,
    GM                       = 1,
    tpole_free      ::Bool   = false,
    logprior                 = nothing,
    maxiter         ::Int    = 200,
    gtol                     = (0.0, 1e-6),
    xtol                     = (0.0, 1e-9),
    ftol                     = (0.0, 1e-12),
    mem             ::Int    = 7,
    verb            ::Bool   = false,
)
    T  = eltype(tessels.unit_xyz)     # vonzeipel_map_and_derivs needs one shared type
    θ  = collect(T, θ0)
    dlb, dub = default_parametric_bounds(; tpole_free=tpole_free)
    length(θ) == length(dlb) ||
        error("fit_parametric: θ0 has $(length(θ)) entries, expected $(length(dlb)) " *
              "for tpole_free=$tpole_free")
    lower = collect(T, lb === nothing ? dlb : lb)
    upper = collect(T, ub === nothing ? dub : ub)

    idx = parametric_free_indices(free; tpole_free=tpole_free)

    # Likelihood only: the prior (if any) is added to the objective, never to χ².
    logπ = build_parametric_logπ(data_epochs, tessels, tepochs, base_params;
                                 intensity_model=intensity_model, band=band,
                                 κ=κ, GM=GM, tpole_free=tpole_free, logprior=nothing)
    # Scale the objective by 1/npoints. χ² here is O(10⁶) and its gradient O(10⁵); at that
    # magnitude vmlmb's bound-constrained line search trips its own `stp == ls.stp`
    # assertion. A uniform rescale moves neither the minimiser nor the relative gradient
    # tolerance, and the reported χ²ᵣ below is computed from `logπ`, not from this.
    npts  = sum(d -> d.nv2 + d.nt3amp + d.nt3phi, data_epochs)
    fscal = one(T) / max(npts, 1)
    full = logprior === nothing ? (θ -> -fscal * logπ(θ)) :
                                  (θ -> -fscal * (logπ(θ) + logprior(θ)))

    # Reduced objective: θ = θ_frozen + S·z, a constant scatter matrix, so Zygote
    # differentiates through it with no special handling and the frozen entries
    # contribute exactly zero gradient.
    θ_frozen = copy(θ);  θ_frozen[idx] .= zero(T)
    S = zeros(T, length(θ), length(idx))
    for (j, i) in enumerate(idx); S[i, j] = one(T); end
    objective = z -> full(θ_frozen .+ S * z)

    neval  = Ref(0)
    gnorm0 = Ref(zero(T))
    function fg!(x, g)
        neval[] += 1
        val, back = Zygote.pullback(objective, x)
        g .= back(one(eltype(x)))[1]
        neval[] == 1 && (gnorm0[] = norm(g))
        return val
    end

    # xtol/ftol are tightened well below vmlmb's defaults (1e-7 / 1e-8): every replicate
    # starts at the full-data θ̂, and with the default relative step test the optimiser
    # declares victory before it has moved, which collapses the bootstrap scatter to zero.
    ẑ = OptimPackNextGen.vmlmb(fg!, θ[idx]; lower=lower[idx], upper=upper[idx],
                                mem=min(mem, length(idx)), maxiter=maxiter,
                                blmvm=false, gtol=gtol, xtol=xtol, ftol=ftol, verb=verb)
    θ̂ = copy(θ);  θ̂[idx] .= ẑ

    # vmlmb prints its stop reason but does not return it, and ‖g‖ is not a usable proxy
    # here: a highly resolved star gives a χ² so steep in diameter that ‖g‖ stays O(10) at
    # a minimum whose σ is ~1e-3 mas. vmlmb usually stops on its x/f tests instead. So do
    # not claim convergence — report the gradient and flag the one failure mode that is
    # unambiguous: running out of iterations.
    gfin  = Zygote.gradient(objective, ẑ)[1]
    gnorm = norm(gfin)
    ndof  = max(npts - length(idx), 1)
    chi2  = -2 * logπ(θ̂)
    info  = (evaluations = neval[], gnorm = gnorm, ndof = ndof, free = idx,
             gratio = gnorm0[] > 0 ? gnorm / gnorm0[] : NaN,
             hit_maxiter = neval[] >= maxiter)
    return θ̂, chi2 / ndof, info
end

"""
    bootstrap_parametric(data_epochs, tessels, tepochs, base_params; θ0, kwargs...)
        -> ParametricBootstrap

Block bootstrap of the parametric fit: fit the full dataset, then refit `nboot` resampled
replicates and summarise the parameter scatter.

Every replicate is warm-started at the full-data θ̂ and run to convergence, so the spread
measures the data rather than the optimizer's stopping rule.

Fit keywords (`lb`, `ub`, `intensity_model`, `band`, `κ`, `GM`, `tpole_free`, `logprior`,
`maxiter`, `gtol`, `mem`) go to [`fit_parametric`](@ref); the rest go to the core
[`bootstrap_parametric`](@ref) method.
"""
function bootstrap_parametric(data_epochs::AbstractVector, tessels, tepochs, base_params;
    θ0,
    free                     = nothing,
    lb                       = nothing,
    ub                       = nothing,
    intensity_model ::Symbol = :linear,
    band                     = nothing,
    κ                        = 50,
    GM                       = 1,
    tpole_free      ::Bool   = false,
    logprior                 = nothing,
    maxiter         ::Int    = 200,
    gtol                     = (0.0, 1e-6),
    mem             ::Int    = 7,
    verb            ::Bool   = true,
    kwargs...
)
    dlb, dub = default_parametric_bounds(; tpole_free=tpole_free)
    lower = lb === nothing ? dlb : lb
    upper = ub === nothing ? dub : ub
    idx   = parametric_free_indices(free; tpole_free=tpole_free)
    names = parametric_param_names(; tpole_free=tpole_free)[idx]
    fitkw = (free=idx, lb=lower, ub=upper, intensity_model=intensity_model, band=band,
             κ=κ, GM=GM, tpole_free=tpole_free, logprior=logprior, maxiter=maxiter,
             gtol=gtol, mem=mem)

    # Full-data fit first: it seeds every replicate and warms up Zygote's pullback
    # compilation on this thread, before any replicate task is spawned.
    verb && println("Fitting the full dataset...")
    θ̂, chi2r, info = fit_parametric(data_epochs, tessels, tepochs, base_params;
                                    θ0=θ0, verb=verb, fitkw...)
    if verb
        @printf("Full-data fit: χ²ᵣ = %.4f  (%d evaluations, ‖g‖/‖g₀‖ = %.2e)\n",
                chi2r, info.evaluations, info.gratio)
        @printf("  free: %s = %s\n", join(names, ", "),
                string(round.(Float64.(θ̂[idx]), sigdigits=5)))
        info.hit_maxiter && @warn "Full-data fit stopped at maxiter=$maxiter"
    end

    # Replicates warm-start at the full-data θ̂ and report only the free entries, which is
    # what the driver summarises.
    fitfun = function (rdata)
        θr, cr, ir = fit_parametric(rdata, tessels, tepochs, base_params;
                                    θ0=θ̂, verb=false, fitkw...)
        return θr[idx], cr, ir
    end

    return bootstrap_parametric(fitfun, data_epochs;
                                x_opt=θ̂[idx], lb=lower[idx], ub=upper[idx],
                                list_free_params=names, verb=verb, kwargs...)
end

end # module
