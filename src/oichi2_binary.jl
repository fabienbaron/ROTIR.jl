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
                     occultation = false)
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

    indx1 = star1.index_quads_visible
    xw1 = I1[indx1] .* star1.vis_weights[indx1] .* star1.ldmap[indx1]  # soft visibility × limb darkening
    ow1 === nothing || (xw1 = xw1 .* ow1[indx1])                       # mutual occultation
    flux1 = dot(star1.polyflux, xw1)
    F1 = star1.polyft * xw1

    indx2 = star2.index_quads_visible
    xw2 = I2[indx2] .* star2.vis_weights[indx2] .* star2.ldmap[indx2]  # soft visibility × limb darkening
    ow2 === nothing || (xw2 = xw2 .* ow2[indx2])
    flux2 = dot(star2.polyflux, xw2)
    F2 = star2.polyft * xw2

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
                       intensity_model=intensity_model, band=band, occultation=occultation)
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
