# binary_geometry.jl
# ---------------------------------------------------------------------------
# Two stars in ONE frame.
#
# `create_star` builds a star at its own origin and orients it with
# `rotate_star(2πt/P_rot, inclination, position_angle)` — a spin about a fixed sky axis.
# That is right for a single rotating star, but it leaves nothing tying a Roche lobe's
# tidal bulge to the line of centres, so:
#
#   * the bulge does not follow the companion around the orbit, and
#   * there is no way to ask "where is the companion, relative to this surface element?",
#     which is exactly what mutual irradiation needs.
#
# `demos/spica_binary_roche.jl` used to work around this by passing t = 0 for every
# epoch, freezing both the tidal shape and its orientation. This file replaces that
# workaround: it derives one body→sky rotation from the orbit itself and applies it to
# both components.
#
# ---------------------------------------------------------------------------
# Frames
# ---------------------------------------------------------------------------
# ROTIR sky frame  : (1) West, (2) North, (3) toward the observer.  RIGHT-handed
#                    (East-left/North-up on the sky makes (East, North, toward-obs)
#                    LEFT-handed, hence the West convention).
# orbits.jl frame  : (x = North, y = East, z = +receding).
#   ⇒ sky = (-y, x, -z)                                    [`sky_of_orbit`]
# which is `orbit_to_rotir_offset` (src/oichi2_binary.jl) plus the z component.
#
# Binary co-rotating frame:
#   x̂ = (r₂ − r₁)/D   (primary → secondary)
#   ẑ = orbital angular momentum
#   ŷ = ẑ × x̂
#
# BOTH Roche lobes use this SAME x̂. Reading the potentials:
#   compute_potential_primary   puts the companion at +D x̂ (geometry_rochelobe.jl:135)
#   compute_potential_secondary puts the companion at −D x̂ (geometry_rochelobe.jl:148,
#       confirmed by solve_R_L1 using θ = −π/2 for the secondary at :176)
# i.e. the primary's body +x̂ points AT the secondary and the secondary's body +x̂ points
# AWAY from the primary — which is the same direction in space. One matrix orients both.
# ---------------------------------------------------------------------------

"""
    sky_of_orbit(x, y, z) -> (west, north, toward_observer)

Map a vector from the orbital-code frame of `src/orbits.jl` (x = North, y = East,
z = positive receding) into ROTIR's sky frame (West, North, toward observer).
"""
@inline sky_of_orbit(x, y, z) = (-y, x, -z)

"""
    binary_frame(bparams, tepoch_jd) -> (R, D_mas, r_rel_sky)

Orientation of the binary at time `tepoch_jd` (JD).

* `R` — 3×3 matrix whose **columns** are `x̂, ŷ, ẑ` of the co-rotating binary frame
  expressed in ROTIR's sky frame, so `v_sky = R * v_body`. ROTIR stores vertices as row
  vectors, so the corresponding right-multiply is `xyz * transpose(R)`.
* `D_mas` — instantaneous centre-to-centre separation (mas).
* `r_rel_sky` — the secondary's offset from the primary in the sky frame (mas),
  `= R * [D_mas, 0, 0]`.

`R` is built from `binary_orbit_abs` and `compute_coeff` directly, so it cannot drift out
of sync with the orbit used everywhere else (`orbit_to_rotir_offset`, `plot2d_binary`).
`bparams` may be either a `binaryparameters` NamedTuple or a Roche `star_params`
NamedTuple — both carry `i, Ω, ω, P, a, e, T0, q, dP`.
"""
function binary_frame(bparams, tepoch_jd)
    x1, y1, z1, x2, y2, z2 = binary_orbit_abs(bparams, Float64(tepoch_jd))
    r_rel_sky = collect(sky_of_orbit(x2 - x1, y2 - y1, z2 - z1))   # primary -> secondary, mas
    D_mas = norm(r_rel_sky)
    D_mas > 0 || error("binary_frame: degenerate orbit, D = 0 at t = $tepoch_jd")
    x̂ = r_rel_sky ./ D_mas

    # Orbital angular momentum from the (time-independent) Thiele-Innes basis:
    # the relative position is r ∝ cos(ν) e₁ + sin(ν) e₂, so h ∝ e₁ × e₂.
    Ω = bparams.Ω * pi / 180.0
    i = bparams.i * pi / 180.0
    ω = bparams.ω * pi / 180.0
    L1, M1, N1, L2, M2, N2 = compute_coeff(Ω, i, ω)
    e1 = collect(sky_of_orbit(L1, M1, -N1))   # compute_xyz_abs returns z = -astro_z
    e2 = collect(sky_of_orbit(L2, M2, -N2))
    ẑ = cross(e1, e2); ẑ ./= norm(ẑ)
    ŷ = cross(ẑ, x̂)

    R = hcat(x̂, ŷ, ẑ)          # columns = basis vectors: v_sky = R * v_body
    return R, D_mas, r_rel_sky
end

"""
    create_binary_star(tessels, star_params, bparams, tepoch_jd;
                       secondary=false, T=Float32, κ=T(50),
                       volume_conserving=false, omega=nothing)

One component of a binary, deformed by the *shared* instantaneous separation `D(t)` and
oriented by the *shared* [`binary_frame`](@ref) rather than by a free spin phase.

Identical to [`create_star`](@ref) except for those two substitutions, so everything
downstream (normals, soft visibility, projection, `ldmap`, `setup_oi!`, `polyft`) is
unchanged. Works for any `surface_type`; only type 3 (Roche) actually uses `D`.

`volume_conserving=true` holds the component's *volume* fixed as `D(t)` varies (PHOEBE's
choice, and the physically right invariant on an orbital timescale) instead of holding
`rpole` fixed. It costs a root solve per call — see [`roche_omega_for_volume`](@ref) and
[`roche_omega_table`](@ref) for the cached version. `omega` short-circuits it with a
precomputed equipotential.

The vertices are centred on this component. Its offset from the primary lives in
`center_offsets`, set by [`create_binary_geometry`](@ref).
"""
function create_binary_star(tessels::tessellation, star_params, bparams, tepoch_jd;
                            secondary=false, T=Float32, κ=T(50),
                            volume_conserving=false, omega=nothing)
    R, D_mas, _ = binary_frame(bparams, tepoch_jd)
    p = convert_params(T, star_params)
    D = p.surface_type == 3 ? T(D_mas / p.a) : nothing

    if p.surface_type == 3 && omega === nothing && volume_conserving
        # Root-solve in Float64: the nested Halley stack needs the headroom.
        omega = roche_omega_for_volume(tessels, star_params, Float64(D); secondary=secondary)
    end

    r, xyz = compute_radii(tessels, star_params, tepoch_jd;
                           secondary=secondary, T=T, D=D, omega=omega)

    # v_sky = R * v_body, and ROTIR stores vertices as row vectors  ⇒  right-multiply by Rᵀ.
    npix = tessels.npix
    Rt = T.(transpose(R))
    xyz = reshape(reshape(xyz, (npix * 5, 3)) * Rt, (npix, 5, 3))

    return finish_star(xyz, r, tessels, star_params, tepoch_jd; T=T, κ=κ)
end

"""
    create_binary_geometry(tessels1, params1, tessels2, params2, bparams, tepoch_jd;
                           T=Float32, κ=T(50), volume_conserving=false,
                           omega1=nothing, omega2=nothing) -> (star1, star2)

Both components at one epoch, in one frame. The secondary is built with `secondary=true`
(its potential and gravity use the inverted mass ratio — `params2.q` must already be
`M₁/M₂`, as in `demos/spica_binary_roche.jl`).

The secondary's sky-frame offset from the primary is stored in `star2.center_offsets`
(West, North, toward-observer; mas), with `star1.center_offsets = [0,0,0]`. That single
field is the whole interface used by [`handle_reflection`](@ref) and by the binary
plotting/animation code to place the two meshes relative to each other — the vertex
arrays themselves stay centred on their own component, so `polyft` and
`binary_phase_shift` are unaffected.
"""
function create_binary_geometry(tessels1::tessellation, params1,
                                tessels2::tessellation, params2,
                                bparams, tepoch_jd;
                                T=Float32, κ=T(50), volume_conserving=false,
                                omega1=nothing, omega2=nothing)
    star1 = create_binary_star(tessels1, params1, bparams, tepoch_jd;
                               secondary=false, T=T, κ=κ,
                               volume_conserving=volume_conserving, omega=omega1)
    star2 = create_binary_star(tessels2, params2, bparams, tepoch_jd;
                               secondary=true, T=T, κ=κ,
                               volume_conserving=volume_conserving, omega=omega2)
    _, _, r_rel_sky = binary_frame(bparams, tepoch_jd)
    star1.center_offsets = zeros(T, 3)
    star2.center_offsets = T.(r_rel_sky)
    return star1, star2
end

"""
    create_binary_geometry_multiepochs(tessels1, params1, tessels2, params2, bparams,
                                       tepochs_jd; kwargs...) -> (stars1, stars2)

[`create_binary_geometry`](@ref) over a vector of epochs (JD). Returns two vectors of
`stellar_geometry`, ready for `setup_oi!(data, stars1)` / `setup_oi!(data, stars2)`.
"""
function create_binary_geometry_multiepochs(tessels1::tessellation, params1,
                                            tessels2::tessellation, params2,
                                            bparams, tepochs_jd; kwargs...)
    nepochs = length(tepochs_jd)
    stars1 = Array{stellar_geometry}(undef, nepochs)
    stars2 = Array{stellar_geometry}(undef, nepochs)
    for i in 1:nepochs
        stars1[i], stars2[i] = create_binary_geometry(tessels1, params1, tessels2, params2,
                                                      bparams, tepochs_jd[i]; kwargs...)
    end
    return stars1, stars2
end

# ---------------------------------------------------------------------------
# Volume-conserving Roche shape
# ---------------------------------------------------------------------------

"""
    roche_mesh_volume(tessels, star_params, D, Ω; secondary=false) -> V

Volume (in units of a³) enclosed by the `Ω` equipotential, evaluated on `tessels`:
`V = (1/3) Σᵢ rᵢ³ ΔΩᵢ` with `rᵢ` the radius at each tessel centre and `ΔΩᵢ` its solid
angle on the unit sphere.

This is deliberately the volume of the *mesh ROTIR actually integrates over*, not the
exact analytic volume. It is what makes [`roche_omega_for_volume`](@ref) well-posed —
the same discretisation appears on both sides of `V(Ω) = V_target`, so the conserved
quantity is exactly the thing the forward model sees. It is also ~5× cheaper than
[`roche_volume`](@ref) and uses the same safe `rpole` seed as `update_roche_radii`.
"""
function roche_mesh_volume(tessels::tessellation, star_params, D, Ω; secondary=false)
    p = star_params
    pf = secondary ? compute_potential_secondary : compute_potential_primary
    async_ratio = p.rotation_period / p.P
    q = p.q
    r0 = p.rpole / p.a          # the pole is the minimum radius: always inside the surface
    dΩ = tessel_solid_angles(tessels)
    V = 0.0
    @inbounds for i in 1:tessels.npix
        θ = Float64(tessels.unit_spherical[i, 5, 2])
        φ = Float64(tessels.unit_spherical[i, 5, 3])
        r = solve_radius(r0, Ω, D, θ, φ, q, async_ratio, pf; verbose=false, secondary=secondary)
        V += r^3 * dΩ[i]
    end
    return V / 3
end

"""
    tessel_solid_angles(tessels) -> Vector

Solid angle subtended by each tessel on the unit sphere, from the two-triangle area of
its unit-sphere corners. For HEALPix these are all equal (`≈ 4π/npix`); for the
latitude/longitude scheme they are not, hence the explicit computation.
"""
function tessel_solid_angles(tessels::tessellation)
    n = tessels.npix
    u = tessels.unit_xyz
    A = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        e1 = (u[i,2,1]-u[i,1,1], u[i,2,2]-u[i,1,2], u[i,2,3]-u[i,1,3])
        e2 = (u[i,3,1]-u[i,1,1], u[i,3,2]-u[i,1,2], u[i,3,3]-u[i,1,3])
        e3 = (u[i,4,1]-u[i,1,1], u[i,4,2]-u[i,1,2], u[i,4,3]-u[i,1,3])
        A[i] = (_cross_norm(e1, e2) + _cross_norm(e2, e3)) / 2
    end
    # Flat quads slightly under-cover the sphere; renormalise so Σ ΔΩ = 4π exactly.
    return A .* (4pi / sum(A))
end

"""
    roche_omega_for_volume(tessels, star_params, D; secondary=false, V_target=nothing,
                           D_ref=1.0, rtol=1e-6) -> Ω

Surface equipotential giving this component the same volume at separation `D` as it has
at `D_ref` with its nominal `rpole`.

`update_roche_radii` normally holds `rpole` fixed and lets Ω follow `D`, which makes the
star's *volume* breathe over an eccentric orbit. Physically the volume is the invariant
on an orbital timescale — PHOEBE re-solves Ω from a fixed equivalent-volume radius
(`universe.py:1421-1439`) — and the difference between the two choices *is* the
ellipsoidal variation this whole exercise is about.

Volume is monotonically decreasing in Ω, so the solve is a bracket plus Brent on
[`roche_mesh_volume`](@ref). Each evaluation costs `npix` Halley solves, so prefer
[`roche_omega_table`](@ref) when sweeping an orbit.

**Reference phase.** `D_ref = 1.0` means the conserved volume is the one the star would
have at separation `D = a`. Codes differ here and the choice shifts the absolute size
slightly (not the modulation): ELISa also references `components_distance = 1.0`, whereas
Wilson–Devinney and PHOEBE legacy reference **periastron** (`D = 1 − e`). PHOEBE 2's
`requiv` is phase-independent. Pass `D_ref = 1 - bparams.e` to match the WD convention.
"""
function roche_omega_for_volume(tessels::tessellation, star_params, D; secondary=false,
                                V_target=nothing, D_ref=1.0, rtol=1e-6)
    p = star_params
    async_ratio = p.rotation_period / p.P
    q = p.q
    pot_fun = secondary ? compute_potential_secondary : compute_potential_primary

    Ω_ref, _ = pot_fun(p.rpole / p.a, D_ref, 0.0, 0.0, q, async_ratio)
    Vt = V_target === nothing ?
         roche_mesh_volume(tessels, p, D_ref, Ω_ref; secondary=secondary) : V_target

    f(Ω) = roche_mesh_volume(tessels, p, D, Ω; secondary=secondary) - Vt

    Ω0, _ = pot_fun(p.rpole / p.a, D, 0.0, 0.0, q, async_ratio)
    lo, hi = Ω0, Ω0
    step = 0.02 * abs(Ω0) + 1e-3
    for _ in 1:60                     # V(lo) > Vt : star too big, so lower Ω
        f(lo) > 0 && break
        lo -= step; step *= 1.3
    end
    step = 0.02 * abs(Ω0) + 1e-3
    for _ in 1:60                     # V(hi) < Vt
        f(hi) < 0 && break
        hi += step; step *= 1.3
    end
    if !(f(lo) > 0 && f(hi) < 0)
        @warn "roche_omega_for_volume: could not bracket Ω (D=$D, secondary=$secondary); " *
              "falling back to the fixed-rpole equipotential"
        return Ω0
    end
    return brent_root(f, lo, hi; tol=rtol * abs(Ω0))
end

"""
    roche_omega_table(tessels, star_params, bparams; secondary=false, npoints=15) -> Ω(D)

Precompute the volume-conserving equipotential on a grid of separations spanning
`[1−e, 1+e]` (units of a) and return a linear-interpolating closure `Ω(D)`. One table
amortises the solve across a whole movie or epoch set; the interpolation error is far
below the mesh resolution.
"""
function roche_omega_table(tessels::tessellation, star_params, bparams;
                           secondary=false, npoints=15)
    e = bparams.e
    Ds = collect(range(1 - e, 1 + e, length=npoints))
    Ωs = [roche_omega_for_volume(tessels, star_params, D; secondary=secondary) for D in Ds]
    function Ωof(D)
        D <= Ds[1] && return Ωs[1]
        D >= Ds[end] && return Ωs[end]
        k = searchsortedlast(Ds, D)
        w = (D - Ds[k]) / (Ds[k+1] - Ds[k])
        return (1 - w) * Ωs[k] + w * Ωs[k+1]
    end
    return Ωof
end

# ---------------------------------------------------------------------------
# Mutual-occultation diagnostic
# ---------------------------------------------------------------------------

"""
    projected_separation(bparams, tepoch_jd) -> ρ

Sky-projected centre-to-centre separation (mas) — the quantity that decides whether the
two components overlap on the image plane, as opposed to the 3-D separation returned by
[`binary_frame`](@ref).
"""
function projected_separation(bparams, tepoch_jd)
    ox, oy = orbit_to_rotir_offset(bparams, Float64(tepoch_jd))
    return hypot(ox, oy)
end

"""
    check_binary_overlap(star1, star2, bparams, tepochs_jd; verbose=true) -> Vector{Bool}

Warn where the two components overlap on the sky.

`binary_cvis` adds the two components' visibilities with **no mutual occultation**: the
hidden part of the far star still contributes flux. That is harmless while the disks are
well separated and wrong once they overlap. Spica sits right at the limit
(R₁+R₂ ≈ 0.43 a versus a projected separation of ≈ 0.44 a at conjunction), so this check
is not academic — any conclusion drawn from a flagged epoch carries the caveat.

Returns one `Bool` per epoch, `true` where `ρ < R₁ + R₂`.
"""
function check_binary_overlap(star1, star2, bparams, tepochs_jd; verbose=true)
    s1 = star1 isa AbstractVector ? star1 : [star1]
    s2 = star2 isa AbstractVector ? star2 : [star2]
    ts = tepochs_jd isa AbstractVector ? tepochs_jd : [tepochs_jd]
    n = length(ts)
    flags = falses(n)
    for k in 1:n
        g1 = s1[min(k, length(s1))]
        g2 = s2[min(k, length(s2))]
        # Largest on-sky extent of each mesh (projected radius, mas).
        r1 = sqrt(maximum(g1.proj_west .^ 2 .+ g1.proj_north .^ 2))
        r2 = sqrt(maximum(g2.proj_west .^ 2 .+ g2.proj_north .^ 2))
        ρ = projected_separation(bparams, ts[k])
        flags[k] = ρ < r1 + r2
        if verbose && flags[k]
            @warn "Epoch $k (JD $(round(ts[k], digits=4))): components overlap on the sky " *
                  "(ρ = $(round(ρ, digits=4)) mas < R₁+R₂ = $(round(r1 + r2, digits=4)) mas). " *
                  "binary_cvis has no mutual occultation, so the far component's hidden " *
                  "face is still contributing flux at this epoch."
        end
    end
    if verbose && !any(flags)
        @info "check_binary_overlap: no epoch overlaps on the sky " *
              "(closest approach $(round(minimum(projected_separation.(Ref(bparams), ts)), digits=4)) mas)."
    end
    return flags
end
