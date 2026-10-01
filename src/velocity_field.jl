# velocity_field.jl
# ---------------------------------------------------------------------------
# The per-tessel line-of-sight velocity, in km/s.
#
# WHY. Nothing in the interferometric path has ever carried a velocity: the geometry
# struct has no such field, `rotate_star` computes a rotation PHASE from
# `rotation_period` but never differentiates it, and the only velocity in the package is
# the loose scalar `vsini` threaded through the (dead) Doppler-imaging code. Without a
# velocity field there are no spectrally-resolved observables at all — no differential
# visibility across a line, no line profile from the same surface the visibilities come
# from.
#
# THE ONE IDENTITY THIS RESTS ON. A rigid rotation about the spin axis has
#
#     v = Ω n̂ × r
#
# in ANY frame, because a rotation matrix maps cross products to cross products. ROTIR
# already has both pieces in the sky frame: `rot_vertex`'s third ROW is the spin axis in
# sky coordinates (geometry.jl:182-184, and the comment at :234 says so), and
# `proj_west`/`proj_north` are the sky-frame x and y. So the line-of-sight component needs
# no body-frame reconstruction and, more importantly, no new geometry derivatives — it is
# a function of the very quantities `project_geometry` already returns with a hand-written
# rrule, so Zygote composes it for free.
#
#     v_rad = Ω sin i (cos(PA)·x + sin(PA)·y)        [positive = receding]
#
# from v_rad = −(Ω n̂ × r)_z with n̂ = (−sin i sin PA, sin i cos PA, cos i). Positive is
# RECEDING, the spectroscopic convention, so λ_obs = λ₀(1 + v_rad/c) and `vgamma` adds
# with its usual sign. Note ROTIR's +z points at the observer (`center_offsets` is
# documented "(West, North, toward-observer)", and a visible tessel has nz > 0), hence the
# minus sign in that step.
#
# WHY POSITIONS AND NOT NORMALS. `src/di.jl:454` computes its velocity from
# `normals[:,1:2]` rather than from positions. For a SPHERE the two are parallel and it
# does not matter. For an oblate rotator — the entire point of surface_type 2 — it does:
# the equator bulges, so an equatorial tessel sits at a larger radius and moves faster,
# and the surface normal of an oblate figure is not radial. Using positions is what makes
# V(θ) follow the Roche shape, which is also what PMOIRED does (rotastar.py:160-163).
#
# THE CENTROID IS EXACT, NOT AN APPROXIMATION. For rigid rotation v_rad is a LINEAR
# function of position, so its mean over a flat quadrilateral equals its value at the
# quad's centroid, which is the mean of the four vertices. Averaging the projected
# vertices is therefore the exact area-weighted mean velocity of the tessel, not a
# midpoint approximation — and it is the right quantity, because the polygon Fourier
# transform integrates over that same quad. (This exactness is lost when `B_rot ≠ 0`, see
# below.)
# ---------------------------------------------------------------------------

"Speed of light in km/s. Not 3e5: `src/di.jl:371` used that and it is wrong in the 4th digit."
const _C_KMS = 299792.458

"""
    los_velocity_from_proj(pw, pn, incl_deg, pa_deg, period_days, d_pc;
                           vgamma = 0, B_rot = 0, cosθ = nothing) -> Vector

Line-of-sight velocity per tessel in km/s, positive RECEDING, from the projected
sky-frame vertex coordinates.

`pw`, `pn` are `(npix, 4)` projected West/North in mas, exactly as
`stellar_geometry.proj_west`/`.proj_north` and as `project_geometry` returns them — which
is what lets this differentiate with respect to the geometry parameters through their
existing rrule rather than needing one of its own.

`d_pc` converts mas to a length; without it an angular rate is not a velocity.

## Differential rotation

`B_rot ≠ 0` applies a solar-type law in the colatitude θ,

    Ω(θ) = Ω₀ (1 − B_rot cos²θ)

so `rotation_period` is the EQUATORIAL period and `B_rot > 0` makes the poles lag. It
needs `cosθ` (body-frame colatitude cosines, `star.vertices_spherical[:,5,3]`-adjacent —
see [`los_velocity`](@ref)).

Two honest caveats. First, the centroid identity above stops being exact, because Ω now
varies across a tessel; evaluating at the centre is then an approximation, good while the
mesh resolves the latitudinal shear. Second, and more fundamental: a differentially
rotating surface has **no rotational potential**, so its shape does not follow ROTIR's
Roche factor at all — the figure would have to come from `g_eff·ds = 0` instead, and the
gravity-darkening law would need Zorec et al. (2017) rather than ELR (the reasoning is
spelled out at the end of `src/gravity_darkening.jl`). So `B_rot ≠ 0` is a kinematic
perturbation on a shape that is not self-consistent with it. Useful for asking whether the
data want shear at all; not a substitute for the real model.
"""
function los_velocity_from_proj(pw::AbstractMatrix, pn::AbstractMatrix,
                                incl_deg, pa_deg, period_days, d_pc;
                                vgamma = zero(eltype(pw)),
                                B_rot = zero(eltype(pw)),
                                cosθ = nothing)
    # The MESH sets the precision, not the parameters — `compute_radii` and `rotate_star`
    # both do this via `convert_params`, and promoting here instead would make a Float32 mesh
    # silently run in Float64 whenever a caller passed an unconverted Float64 angle.
    T = float(eltype(pw))
    deg = T(π) / T(180)
    si  = sin(T(incl_deg) * deg)
    sp, cp = sincos(T(pa_deg) * deg)
    Ω = T(2π) / (T(period_days) * T(_DAY_S))            # rad/s, equatorial

    # mas -> m, then m/s -> km/s. One factor, applied once.
    scale = T(d_pc) * T(_MASPC_TO_M) / T(1000)

    # Centroid of each quad = exact mean position, since v_rad is linear in position.
    xc = dropdims(sum(pw, dims = 2), dims = 2) ./ T(4)
    yc = dropdims(sum(pn, dims = 2), dims = 2) ./ T(4)

    v = (Ω * si * scale) .* (cp .* xc .+ sp .* yc)
    if !iszero(B_rot)
        cosθ === nothing && throw(ArgumentError(
            "los_velocity_from_proj: B_rot = $B_rot needs `cosθ`, the body-frame " *
            "colatitude cosines, to evaluate Ω(θ) = Ω₀(1 − B_rot cos²θ)"))
        v = v .* (one(T) .- T(B_rot) .* T.(cosθ) .^ 2)
    end
    return v .+ T(vgamma)
end

"""
    los_velocity(star, star_params) -> Vector

[`los_velocity_from_proj`](@ref) for an existing `stellar_geometry`, reading `inclination`,
`position_angle`, `rotation_period`, `d` and the optional `vgamma`/`B_rot` from
`star_params`.

Errors without a distance: an angular rotation rate is not a velocity, and silently
substituting one would put a plausible-looking but meaningless number into a line profile.
"""
function los_velocity(star, star_params)
    has_physical_scale(star_params) || throw(ArgumentError(
        "los_velocity: `star_params` carries no distance `d`, so mas/s cannot be turned " *
        "into km/s. Add `d` (pc) — see src/stellar_physics.jl."))
    vγ = hasproperty(star_params, :vgamma) ? star_params.vgamma : 0.0
    B  = hasproperty(star_params, :B_rot)  ? star_params.B_rot  : 0.0
    # Body-frame colatitude of each tessel centre; slot 5 is the centre, index 2 is θ.
    cθ = iszero(B) ? nothing : cos.(@view star.vertices_spherical[:, 5, 2])
    return los_velocity_from_proj(star.proj_west, star.proj_north,
                                  star_params.inclination, star_params.position_angle,
                                  star_params.rotation_period, star_params.d;
                                  vgamma = vγ, B_rot = B, cosθ = cθ)
end

"""
    doppler_lambda(λ₀, v_rad) -> λ

Observed wavelength(s) for a rest wavelength `λ₀` and a receding-positive velocity in km/s:

    λ = λ₀ (1 + v_rad/c)

Broadcasts, so `λ₀` may be a scalar and `v_rad` a per-tessel vector — which is the shape
the spectral cube needs, one shifted wavelength per tessel per channel.

The first-order form is deliberate and sufficient here: the relativistic correction is
O((v/c)²), i.e. 3e-7 at 300 km/s, far below the 1e-4 a Doppler-imaging grid resolves.
What is NOT optional is doing the resampling in log λ rather than linear λ once a spectrum
is interpolated — `src/di.jl:371` shifted in linear λ and its own comment at :323 admits
the error varies across the window.
"""
@inline function doppler_lambda(λ0, v_rad)
    # `_C_KMS` is a Float64 const and `1` an Int, so writing this as `λ0 .* (1 .+ v_rad/_C_KMS)`
    # promotes a Float32 model to Float64 for the whole wavelength vector. It is invisible from
    # outside — `provider_intensity` narrows its OUTPUT back to the input eltype — but the
    # intermediate λ allocates at double width and the grid query runs in Float64.
    T = float(promote_type(eltype(λ0), eltype(v_rad)))
    return λ0 .* (one(T) .+ v_rad ./ T(_C_KMS))
end

"""
    velocity_field_summary(star, star_params) -> NamedTuple

`(vmin, vmax, vspan, vsini_proj)` of the visible hemisphere, in km/s.

`vsini_proj` is half the peak-to-peak line-of-sight velocity over visible tessels — the
quantity a rotational broadening measurement actually constrains. Comparing it against
[`projected_veq`](@ref), which comes from the geometry alone, is the cross-check that the
velocity field and the surface agree; they are computed by completely different routes.
"""
function velocity_field_summary(star, star_params)
    v = los_velocity(star, star_params)
    vis = star.index_quads_visible
    vv = @view v[vis]
    lo, hi = extrema(vv)
    return (vmin = lo, vmax = hi, vspan = hi - lo, vsini_proj = (hi - lo) / 2)
end
