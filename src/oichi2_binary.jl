# Binary star forward model
# Computes combined complex visibilities for two stars using phase-shifted polygon FTs.
# Star 1 is at the origin; star 2 is offset by its orbital position.

"""
    orbit_to_rotir_offset(bparams, tepoch_jd) -> (offset_x, offset_y)

Convert orbital position to ROTIR's projected coordinate frame.
Returns the secondary's offset relative to the primary in mas, in ROTIR's (West, North) frame.

The orbital code (`binary_orbit_abs`) returns positions where x=North, y=East.
ROTIR's projected coordinates have proj_west=West, proj_north=North.
"""
function orbit_to_rotir_offset(bparams, tepoch_jd)
    x1, y1, z1, x2, y2, z2 = binary_orbit_abs(bparams, tepoch_jd)
    dx_north = x2 - x1
    dy_east  = y2 - y1
    offset_x = -dy_east   # West = -East
    offset_y = dx_north   # North
    return offset_x, offset_y
end

"""
    orbit_to_rotir_offset(bparams, tepochs::AbstractVector) -> (offsets_x, offsets_y)

Vectorised over epochs: same `(West, North)` offsets in mas as the scalar method, returned
as two arrays.

This is not just convenience — it is the form that can be differentiated affordably.
`binary_orbit_abs` is scalar-only and routes through the true anomaly and the individual
component radii; calling it in a loop makes reverse-mode AD tape one closure per epoch, and
the per-call overhead then dominates completely (measured: **512× the primal** for 55
epochs, versus ~4–13× for every other stage of the visibility model). Broadcasting instead
tapes a handful of array operations regardless of how many epochs there are.

The relative orbit is used directly rather than differencing the two absolute positions:
the components sit at true anomalies `υ` and `υ+π` with radii `r₂ = D/(1+q)` and
`r₁ = qD/(1+q)`, which sum to `D`, so `r₂ − r₁` is exactly the relative vector and `q`
drops out. That also removes the `υ ∈ [0, π]` branch in `binary_orbit_abs`, which is a
discrete test that AD would have to step around.
"""
function orbit_to_rotir_offset(bparams, tepochs::AbstractVector)
    # Integer literals and `pi` (an Irrational) throughout, never `180.0`: both specialise to
    # whatever float type the elements already are, so a Float32 parameter set stays Float32
    # instead of being silently widened here.
    Ω   = bparams.Ω * pi / 180
    inc = bparams.i * pi / 180
    E   = compute_eccentric_anomaly(bparams, tepochs)
    ω0  = bparams.ω * pi / 180
    dω  = hasproperty(bparams, :dω) ? bparams.dω : zero(typeof(bparams.ω))
    # Function barrier on the ω branch. `omega_at` returns a SCALAR when dω = 0 and a
    # VECTOR otherwise, so calling it inline gives ω the type `Union{Float64,Vector}` and
    # every broadcast below then infers as `Any` — JET flags each one as a runtime dispatch,
    # and reverse-mode AD has to box the lot. Dispatching to `_relative_offsets` first means
    # each branch is separately specialised and fully typed; both return the same
    # `Tuple{Vector,Vector}`, so this method is itself type-stable.
    return dω == 0 ? _relative_offsets(bparams, Ω, inc, E, ω0) :
                     _relative_offsets(bparams, Ω, inc, E,
                                       ω0 .+ (dω * pi / 180) .* (tepochs .- bparams.T0))
end

function _relative_offsets(bparams, Ω, inc, E, ω)
    a  = bparams.a
    e  = bparams.e
    β  = sqrt(1 - e^2)
    cE = cos.(E);  sE = sin.(E)
    cΩ = cos(Ω);   sΩ = sin(Ω);  ci = cos(inc)
    cω = cos.(ω);  sω = sin.(ω)
    L1 =  cΩ .* cω .- sΩ .* sω .* ci
    M1 =  sΩ .* cω .+ cΩ .* sω .* ci
    L2 = -cΩ .* sω .- sΩ .* cω .* ci
    M2 = -sΩ .* sω .+ cΩ .* cω .* ci
    north = a .* (L1 .* cE .+ β .* L2 .* sE .- e .* L1)
    east  = a .* (M1 .* cE .+ β .* M2 .* sE .- e .* M1)
    return -east, north                        # West = −East
end

"""
    binary_phase_shift(uv, offset_x, offset_y) -> Vector{Complex}

Compute per-baseline phase shift for a star displaced by (offset_x, offset_y) mas
in ROTIR's projected frame (West, North).

Uses the same kx/ky sign convention as the polygon FT in `setup_polyft_single`.
"""
# The phase argument kx·Δx + ky·Δy runs to many radians at long baselines, so its precision
# is set here, not downstream. Follow the inputs.
function binary_phase_shift(uv, offset_x, offset_y;
                            T = float(real(promote_type(eltype(uv), typeof(offset_x),
                                                        typeof(offset_y)))))
    C = T(180 * 3600000)
    kx = uv[1,:] .* T(-pi / C)
    ky = uv[2,:] .* T( pi / C)
    return cis.(-T(pi) .* (kx .* offset_x .+ ky .* offset_y))
end

"""
    binary_cvis(x1, star1, x2, star2, phase_shift; intensity_model=:linear, band=nothing)
        -> Vector{Complex}

Compute combined complex visibilities for a binary system.
Star 1 is at the origin; star 2's visibilities are multiplied by `phase_shift`.
Each map is weighted by its own soft visibility × limb-darkening map (`ldmap`, built by
`create_star`), so the two components' fluxes — and hence their flux ratio in the joint
normalization — are limb-darkened consistently with the single-star path.

`intensity_model = :linear` (default) uses the maps directly as surface brightness, the
Rayleigh–Jeans proxy. `:planck` treats them as *temperature* maps and converts them with
[`intensity`](@ref) at wavelength `band` (metres).

The flux ratio is where this matters most: for a 25300 K / 20585 K pair the linear proxy
misstates it by −4.1 % in H and −13.1 % in V. The non-dimensional Planck form
(λ⁵ and all constants dropped) stays exact here because both components are evaluated at
the same `band`, so the discarded prefactor cancels from the ratio as well as from the
flux normalization.

`occultation` folds in mutual occultation: pass `true` to compute it on the fly from the
two components' `center_offsets`, or a precomputed `(w1, w2)` pair from
[`occultation_weights`](@ref) to reuse it across a scan, or `:exact`/`:soft` to pick the
method. The default `false` sums the two
components unconditionally, which is wrong once their disks overlap on the sky — see
`check_binary_overlap` for the epochs where that happens.
"""
function binary_cvis(x1, star1, x2, star2, phase_shift;
                     intensity_model::Symbol = :linear, band = nothing,
                     occultation = false, data = nothing)
    I1 = intensity_model === :linear ? x1 : intensity(x1, intensity_model, band)
    I2 = intensity_model === :linear ? x2 : intensity(x2, intensity_model, band)

    ow1, ow2 = if occultation === false || occultation === nothing
        (nothing, nothing)
    elseif occultation === true
        occultation_weights(star1, star2)[1:2]
    elseif occultation isa Symbol
        occultation_weights(star1, star2; method=occultation)[1:2]
    else
        occultation[1], occultation[2]
    end

    # MATRIX-FREE when `setup_oi!` has not been run, exactly as `observables` chooses for a
    # single star — and it needs `data` to do it, since the uv points are what the transform is
    # evaluated at. Whichever `POLYFT_BACKEND[]` is selected then applies to a binary too:
    # :nufft (the default), :turbo or :scalar.
    #
    # The dense route stays for the case it was built for — imaging, where the geometry is
    # fixed and the same matrix is reused over hundreds of iterations.
    F1, flux1, F2, flux2 = if data === nothing || !isempty(star1.polyft)
        indx1 = star1.index_quads_visible
        xw1 = I1[indx1] .* star1.vis_weights[indx1] .* star1.ldmap[indx1]  # soft vis × LD
        ow1 === nothing || (xw1 = xw1 .* ow1[indx1])                       # mutual occultation
        indx2 = star2.index_quads_visible
        xw2 = I2[indx2] .* star2.vis_weights[indx2] .* star2.ldmap[indx2]
        ow2 === nothing || (xw2 = xw2 .* ow2[indx2])
        (star1.polyft * xw1, dot(star1.polyflux, xw1),
         star2.polyft * xw2, dot(star2.polyflux, xw2))
    else
        f1, fl1 = fused_cvis_parts(x1, star1, data; intensity_model = intensity_model,
                                   band = band, extra_weights = ow1)
        f2, fl2 = fused_cvis_parts(x2, star2, data; intensity_model = intensity_model,
                                   band = band, extra_weights = ow2)
        (f1, fl1, f2, fl2)
    end

    return (F1 .+ F2 .* phase_shift) ./ (flux1 + flux2)
end

"""
    binary_observables(x1, star1, x2, star2, data, phase_shift;
                       intensity_model=:linear, band=nothing) -> (v2, t3amp, t3phi)

Compute model observables (V2, T3amp, T3phi) for a binary system.
Uses `cvis_to_obs` (shared with single-star path) for the cvis→observables step.
See [`binary_cvis`](@ref) for the intensity-model keywords.
"""
function binary_observables(x1, star1, x2, star2, data, phase_shift;
                            intensity_model::Symbol = :linear, band = nothing,
                            occultation = false)
    cvis = binary_cvis(x1, star1, x2, star2, phase_shift;
                       intensity_model=intensity_model, band=band, occultation=occultation,
                       data=data)
    return cvis_to_obs(cvis, data)
end

"""
    binary_chi2_f(x1, star1, x2, star2, data, phase_shift; verbose=false,
                  intensity_model=:linear, band=nothing) -> Float

Compute chi-squared for a binary model against interferometric data.
See [`binary_cvis`](@ref) for the intensity-model keywords.
"""
function binary_chi2_f(x1, star1, x2, star2, data, phase_shift; verbose::Bool=false,
                       intensity_model::Symbol = :linear, band = nothing,
                       occultation = false)
    v2_model, t3amp_model, t3phi_model = binary_observables(x1, star1, x2, star2, data, phase_shift;
                                                            intensity_model=intensity_model, band=band,
                                                            occultation=occultation)
    chi2_v2 = sum(abs2, (v2_model .- data.v2) ./ data.v2_err)
    chi2_t3amp = sum(abs2, (t3amp_model .- data.t3amp) ./ data.t3amp_err)
    chi2_t3phi = sum(abs2, mod360(t3phi_model .- data.t3phi) ./ data.t3phi_err)
    if verbose
        printstyled(@sprintf("V2: %.4f ", chi2_v2/data.nv2), color=:red)
        printstyled(@sprintf("T3A: %.4f ", chi2_t3amp/data.nt3amp), color=:blue)
        printstyled(@sprintf("T3P: %.4f\n", chi2_t3phi/data.nt3phi), color=:green)
    end
    return chi2_v2 + chi2_t3amp + chi2_t3phi
end

"""
    binary_chi2_fg(x1, g1, star1, x2, g2, star2, data, phase_shift;
                   verbose=false, occultation=false) -> chi2

χ² of a binary against one epoch, **with the gradient with respect to both surface maps** —
the imaging counterpart of [`binary_chi2_f`](@ref).

`g1` and `g2` are overwritten with `∂χ²/∂x1` and `∂χ²/∂x2`. Both must be as long as their
component's map; tessels outside `index_quads_visible` get a zero, exactly as
`spheroid_chi2_fg` does, because a tessel the epoch cannot see does not enter the transform.

# The derivative

Write each component's weighted map as `xwₖ = xₖ[visible] · vis_weights · ldmap` (times the
occultation weights when they are in play), its unnormalized transform as `Fₖ = Pₖ·xwₖ` and
its flux as `flₖ = polyfluxₖ·xwₖ`. The model is

    V = (F₁ + s∘F₂) / (fl₁ + fl₂)        with `s` the secondary's phase shift

so, with `N = fl₁ + fl₂` and `a` the adjoint source from `cvis_to_chi2_fg`,

    ∂χ²/∂xw₁ = [ ℜ(P₁ᵀa)        − S·polyflux₁/N ] / N
    ∂χ²/∂xw₂ = [ ℜ(P₂ᵀ(s∘a))    − S·polyflux₂/N ] / N,   S = xw₁·ℜ(P₁ᵀa) + xw₂·ℜ(P₂ᵀ(s∘a))

which is the single-star expression with one extra term: the two components share `S`,
because they share the flux they are normalized by. Everything that couples them is in that
one scalar — brighten the primary and the secondary's *relative* contribution falls, and it
is `S` that carries it.

# Why the adjoint comes from OITOOLS rather than from `spheroid_chi2_fg`'s expressions

That function builds `ℜ(Pᵀa)` for V², T3amp and T3φ as nine separate indexed matrix products
with the adjoint written out inline. Duplicating those for two components — and a phase shift
on one of them — is nine more chances to get a conjugate backwards, and a conjugation error
is not uniformly wrong: it leaves flux-like directions nearly right while inverting
phase-like ones. `cvis_to_chi2_fg` returns exactly that adjoint source in one sweep, under
the convention `g_params = ℜ(Jᵀ·g_cvis)` with `J = ∂V/∂p` — no conjugation — which is the
convention `cvis_chi2`'s rrule already documents and `test_reflection.jl` already pins.
Checked against `spheroid_chi2_fg` on a single component: the two agree to 2e-9, the
accumulation-order difference between the two summations of the same χ².

`intensity_model` is deliberately absent: `:planck` makes the map a *temperature* and the
brightness a nonlinear function of it, and this derivative is the linear one. The single-star
`spheroid_chi2_fg` has the same restriction for the same reason.

`occultation` takes the same values as [`binary_cvis`](@ref). The weights depend on the
geometry, not on the maps, so they enter as constants — the derivative is exact with them on.
"""
@views function binary_chi2_fg(x1, g1, star1, x2, g2, star2, data, phase_shift;
                               verbose::Bool = false, occultation = false)
    T = real(float(promote_type(eltype(x1), eltype(x2))))
    i1 = star1.index_quads_visible
    i2 = star2.index_quads_visible
    length(g1) == length(x1) || throw(DimensionMismatch("g1 has $(length(g1)) entries for a " *
                                                        "map of $(length(x1))"))
    length(g2) == length(x2) || throw(DimensionMismatch("g2 has $(length(g2)) entries for a " *
                                                        "map of $(length(x2))"))

    ow1, ow2 = if occultation === false || occultation === nothing
        (nothing, nothing)
    elseif occultation === true
        occultation_weights(star1, star2)[1:2]
    elseif occultation isa Symbol
        occultation_weights(star1, star2; method = occultation)[1:2]
    else
        occultation[1], occultation[2]
    end

    w1 = star1.vis_weights[i1] .* star1.ldmap[i1]
    w2 = star2.vis_weights[i2] .* star2.ldmap[i2]
    ow1 === nothing || (w1 = w1 .* ow1[i1])
    ow2 === nothing || (w2 = w2 .* ow2[i2])
    xw1 = x1[i1] .* w1
    xw2 = x2[i2] .* w2

    F1 = star1.polyft * xw1;  fl1 = dot(star1.polyflux, xw1)
    F2 = star2.polyft * xw2;  fl2 = dot(star2.polyflux, xw2)
    N = fl1 + fl2
    cvis = (F1 .+ F2 .* phase_shift) ./ N

    chi2, a = cvis_to_chi2_fg(_ascomplexfloat(cvis), data; weights = OI_DEFAULT_WEIGHTS)

    gs1 = real(transpose(star1.polyft) * a)
    gs2 = real(transpose(star2.polyft) * (a .* phase_shift))
    S = dot(xw1, gs1) + dot(xw2, gs2)

    g1 .= zero(T);  g2 .= zero(T)
    g1[i1] .= w1 .* (gs1 .- S .* star1.polyflux ./ N) ./ N
    g2[i2] .= w2 .* (gs2 .- S .* star2.polyflux ./ N) ./ N

    if verbose
        # The split costs one more pass over the data and is computed only to be printed:
        # `cvis_to_chi2_fg` returns the total, and a per-observable trace is the only way to
        # see a reconstruction being dragged by one of the three.
        v2m, t3am, t3pm = cvis_to_obs(cvis, data)
        printstyled(@sprintf("V2: %.4f ",
                             sum(abs2, (v2m .- data.v2) ./ data.v2_err) / max(data.nv2, 1)),
                    color = :red)
        printstyled(@sprintf("T3A: %.4f ",
                             sum(abs2, (t3am .- data.t3amp) ./ data.t3amp_err) /
                             max(data.nt3amp, 1)), color = :blue)
        printstyled(@sprintf("T3P: %.4f ",
                             sum(abs2, mod360(t3pm .- data.t3phi) ./ data.t3phi_err) /
                             max(data.nt3phi, 1)), color = :green)
        printstyled(@sprintf("Flux: %.4f + %.4f\n", fl1, fl2), color = :normal)
    end
    return chi2
end

"""
    binary_crit_allepochs_fg(x, g, stars1, stars2, data, phase_shifts;
                             regularizers1=[], regularizers2=[], epochs_weights=[],
                             verbose=false) -> crit

The criterion a binary reconstruction minimises: [`binary_chi2_fg`](@ref) summed over epochs,
plus a regularizer on each component.

`x` is the two maps CONCATENATED, `[x1; x2]`, split at `stars1[1].npix` — one vector because
that is what VMLMB optimises, and one split point because the two components need not share a
tessellation.

Two regularizer lists, not one, and neither defaults to the other. A regularizer entry carries
a precomputed structure built from *its* star — `radflat_bins` bins by that component's radii,
`orthold_direction` is the degenerate direction of that component's geometry — so handing the
primary's list to the secondary regularizes the secondary against the primary's shape. Sharing
is only safe for the purely combinatorial ones (`sobel`, `tv`) and only when the two
tessellations match, which is a condition the caller knows and this function does not.

A regularizer's pixel subset (element 4 of an entry) indexes ITS OWN component's map, not the
concatenated vector: each list is evaluated on its half. So `1:npix` means all of that
component, on either side.

`epochs_weights` behaves as in `spheroid_crit_allepochs_fg`: applied to the per-epoch terms
before the sum, and to the gradient with them.
"""
function binary_crit_allepochs_fg(x, g, stars1, stars2, data, phase_shifts;
                                  regularizers1 = [], regularizers2 = [], epochs_weights = [],
                                  verbose = false, T = eltype(x))
    nepochs = length(data)
    length(stars1) == nepochs && length(stars2) == nepochs ||
        throw(DimensionMismatch("$(length(stars1)) primary and $(length(stars2)) secondary " *
                                "geometries for $(nepochs) epochs"))
    length(phase_shifts) == nepochs ||
        throw(DimensionMismatch("$(length(phase_shifts)) phase shifts for $(nepochs) epochs"))
    n1 = stars1[1].npix
    n2 = stars2[1].npix
    length(x) == n1 + n2 ||
        throw(DimensionMismatch("x has $(length(x)) entries for $(n1) + $(n2) tessels"))

    chi2_t = zeros(T, nepochs)
    # One gradient buffer per epoch, summed after the loop: the epochs are threaded and they
    # all write the same tessels, so a shared accumulator would be a race. Same shape as
    # `spheroid_crit_allepochs_fg`, for the same reason.
    ge = [zeros(T, n1 + n2) for _ in 1:nepochs]
    Threads.@threads for i in 1:nepochs
        chi2_t[i] = binary_chi2_fg(view(x, 1:n1), view(ge[i], 1:n1), stars1[i],
                                   view(x, n1+1:n1+n2), view(ge[i], n1+1:n1+n2), stars2[i],
                                   data[i], phase_shifts[i]; verbose = verbose)
    end
    w = _epoch_weights(epochs_weights, nepochs, T)
    f = w === nothing ? sum(chi2_t) : sum(chi2_t .* w)
    g[:] .= w === nothing ? sum(ge) : sum(w[i] .* ge[i] for i in 1:nepochs)
    if verbose
        printstyled(@sprintf("Total χ²: %.4f\n", f), color = :white)
    end

    if !isempty(regularizers1)
        rg = zeros(T, n1)
        f += spheroid_regularization(view(x, 1:n1), rg; regularizers = regularizers1,
                                     verbose = verbose)
        g[1:n1] .+= rg
    end
    if !isempty(regularizers2)
        rg = zeros(T, n2)
        f += spheroid_regularization(view(x, n1+1:n1+n2), rg; regularizers = regularizers2,
                                     verbose = verbose)
        g[n1+1:n1+n2] .+= rg
    end
    return f
end

"""
    binary_reconstruct_oi(x_start, data, stars1, stars2, phase_shifts;
                          maxiter, regularizers1, regularizers2, callback, ...) -> x

Reconstruct BOTH surface maps of a binary by VMLMB — [`image_reconstruct_oi`](@ref) for two
components.

`x_start` and the result are the two maps concatenated; `split_binary_map` takes them apart.

`phase_shifts[i]` is the secondary's phase ramp at epoch `i`, from
[`binary_phase_shift`](@ref) on that epoch's offset. Precomputed by the caller and constant
through the run, which is the point: imaging fixes the geometry, so the separation is not
being fitted here — it comes from the orbit, or from a fitted offset, and this reconstructs
the surfaces at that separation.

The dense `polyft` route is what runs, not the matrix-free one: the geometry is fixed and the
same two matrices are reused over hundreds of iterations, which is the case they were built
for. Run `setup_oi!(data, stars1)` and `setup_oi!(data, stars2)` before calling.

`callback` follows `image_reconstruct_oi` exactly, including that `x` belongs to the optimiser
and must be copied to be kept, and that the count is evaluations rather than iterations.
"""
function binary_reconstruct_oi(x_start, data, stars1, stars2, phase_shifts;
                               epochs_weights = [], verbose = true, lower = 0, upper = Inf,
                               maxiter = 100, regularizers1 = [], regularizers2 = [],
                               callback = nothing, callback_every::Int = 25)
    nevals = Ref(0)
    crit = function (x, g)
        f = binary_crit_allepochs_fg(x, g, stars1, stars2, data, phase_shifts;
                                     regularizers1 = regularizers1,
                                     regularizers2 = regularizers2,
                                     epochs_weights = epochs_weights, verbose = verbose)
        nevals[] += 1
        if callback !== nothing && (nevals[] == 1 || nevals[] % max(1, callback_every) == 0)
            try
                callback(x, nevals[], f)
            catch err
                @warn "binary_reconstruct_oi: callback failed, continuing without it" exception = err
                callback = nothing
            end
        end
        return f
    end
    return OptimPackNextGen.vmlmb(crit, x_start; verb = verbose, lower = lower, upper = upper,
                                  maxiter = maxiter, blmvm = false, gtol = (0, 1e-8))
end

"""
    split_binary_map(x, stars1) -> (x1, x2)

Take a concatenated binary map apart, splitting at the primary's tessel count.

Views, not copies: the caller usually wants to plot or measure the halves, and the whole
vector is the one VMLMB owns.
"""
split_binary_map(x, stars1) = split_binary_map(x, _npix_of(stars1))
split_binary_map(x, n1::Integer) = (view(x, 1:n1), view(x, n1+1:length(x)))
_npix_of(s) = s isa AbstractVector ? first(s).npix : s.npix
