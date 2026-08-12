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
    # ẑ below is the orbit-plane normal, which is set by Ω and i alone — apsidal motion
    # rotates e1/e2 *within* the plane and leaves e1×e2 invariant. `omega_at` is used
    # anyway so this stays consistent with binary_orbit_abs if the derivation changes.
    ω = omega_at(bparams, tepoch_jd)
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
                       secondary=false, T=eltype(tessels), κ=T(50),
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
                            secondary=false, T=eltype(tessels), κ=T(50),
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
                           T=promote_type(eltype(tessels1), eltype(tessels2)), κ=T(50),
                           volume_conserving=false,
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
                                T=promote_type(eltype(tessels1), eltype(tessels2)), κ=T(50),
                                volume_conserving=false,
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
    async_ratio = synchronicity(p)
    q = p.q
    r0 = p.rpole / p.a          # the pole is the minimum radius: always inside the surface
    dΩ = tessel_solid_angles(tessels)
    # Angles follow the mesh; the ACCUMULATOR stays Float64 regardless. Same rule as the χ²:
    # per-element quantities take the data's type, a summed scalar takes the wider one.
    V = 0.0
    @inbounds for i in 1:tessels.npix
        θ = tessels.unit_spherical[i, 5, 2]
        φ = tessels.unit_spherical[i, 5, 3]
        r = solve_radius(r0, Ω, D, θ, φ, q, async_ratio, pf; verbose=false, secondary=secondary)
        V += r^3 * dΩ[i]
    end
    return V / 3
end

"""
    tessel_solid_angles(tessels) -> Vector{eltype(tessels)}

Solid angle subtended by each tessel on the unit sphere, from the two-triangle area of
its unit-sphere corners. For HEALPix these are all equal (`≈ 4π/npix`); for the
latitude/longitude scheme they are not, hence the explicit computation.

Returns the **mesh's own float type**. The renormalising sum is accumulated in `Float64`
and the scale factor narrowed back, so a Float32 mesh gives `Σ ΔΩ = 4π` to ~3e-8 relative
rather than being silently widened to `Float64`.
"""
function tessel_solid_angles(tessels::tessellation)
    n = tessels.npix
    u = tessels.unit_xyz
    Tel = float(eltype(tessels))
    A = Vector{Tel}(undef, n)
    @inbounds for i in 1:n
        e1 = (u[i,2,1]-u[i,1,1], u[i,2,2]-u[i,1,2], u[i,2,3]-u[i,1,3])
        e2 = (u[i,3,1]-u[i,1,1], u[i,3,2]-u[i,1,2], u[i,3,3]-u[i,1,3])
        e3 = (u[i,4,1]-u[i,1,1], u[i,4,2]-u[i,1,2], u[i,4,3]-u[i,1,3])
        A[i] = (_cross_norm(e1, e2) + _cross_norm(e2, e3)) / 2
    end
    # Flat quads slightly under-cover the sphere; renormalise so Σ ΔΩ = 4π.
    # The normalising sum is accumulated WIDE and the scale factor is then narrowed back to
    # the mesh type — `A .* (4pi / sum(A))` would have promoted the whole array to Float64,
    # because `4pi` is Float64, silently undoing the element type chosen just above.
    return A .* Tel(4π / sum(Float64, A))
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
    async_ratio = synchronicity(p)
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
    projected_limb_profile(star; nbins=72) -> (r_max, nbins)

Azimuthal profile of a component's projected outline: the maximum projected radius of any
vertex, binned by position angle about the component's own projected centre.

A tidally distorted lobe is not a circle on the sky, so a single radius would either
over- or under-occult depending on orientation. `nbins = 72` (5° bins) resolves the ~10%
tidal elongation of a Spica-like lobe with margin.
"""
function projected_limb_profile(star; nbins::Int = 72)
    T = float(real(eltype(star.proj_west)))
    rmax = zeros(T, nbins)
    pw = star.proj_west; pn = star.proj_north
    @inbounds for i in 1:star.npix, v in 1:4
        x = -pw[i, v]; y = pn[i, v]                       # East, North
        r = hypot(x, y)
        b = mod(floor(Int, (atan(y, x) + pi) / (2pi) * nbins), nbins) + 1
        rmax[b] = max(rmax[b], r)
    end
    # Empty bins happen whenever nbins is comparable to the vertex count (a coarse mesh).
    # Fill them by interpolating between the nearest occupied neighbours — filling with
    # the global maximum instead would put radial spikes in the profile.
    if any(iszero, rmax)
        occupied = findall(!iszero, rmax)
        isempty(occupied) && error("projected_limb_profile: no vertices?")
        @inbounds for b in 1:nbins
            iszero(rmax[b]) || continue
            # nearest occupied bin on each side, wrapping
            lo = hi = 0
            for d in 1:nbins
                l = mod(b - d - 1, nbins) + 1; h = mod(b + d - 1, nbins) + 1
                lo == 0 && !iszero(rmax[l]) && (lo = l)
                hi == 0 && !iszero(rmax[h]) && (hi = h)
                lo != 0 && hi != 0 && break
            end
            rmax[b] = (rmax[lo] + rmax[hi]) / 2
        end
    end
    return rmax
end

"""
    convex_hull_2d(xs, ys) -> (hx, hy)

Convex hull of a 2-D point set, counter-clockwise, by Andrew's monotone chain
(O(n log n)).

Used for the exact occultation path: `polygon_convex_clip_area` intersects half-planes,
which only reproduces the clip polygon if that polygon is genuinely convex. An
azimuthally-sampled outline is *not* — with a coarse mesh it is spiky, and the half-plane
intersection then collapses to nothing. The hull of the projected vertices is convex by
construction and is also the geometrically correct silhouette, since the projection of a
convex body is convex.

(For a lobe close to filling its Roche surface the true silhouette has an L1 cusp and is
slightly non-convex; the hull then very slightly over-occults near the tip. Detached
systems, where the silhouette really is convex, are exact.)
"""
function convex_hull_2d(xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real})
    T = float(promote_type(eltype(xs), eltype(ys)))
    n = length(xs)
    n < 3 && return (collect(T.(xs)), collect(T.(ys)))
    idx = sortperm(collect(zip(T.(xs), T.(ys))))
    cross3(o, a, b) = (xs[a]-xs[o])*(ys[b]-ys[o]) - (ys[a]-ys[o])*(xs[b]-xs[o])
    hull = Int[]
    for i in idx                                   # lower hull
        while length(hull) >= 2 && cross3(hull[end-1], hull[end], i) <= 0; pop!(hull); end
        push!(hull, i)
    end
    lower = length(hull) + 1
    for i in Iterators.reverse(idx)                # upper hull
        while length(hull) >= lower && cross3(hull[end-1], hull[end], i) <= 0; pop!(hull); end
        push!(hull, i)
    end
    pop!(hull)                                     # last point repeats the first
    return (T.(xs[hull]), T.(ys[hull]))
end

"""
    limb_radius(profile, φ) -> r

Projected limb radius at position angle `φ` (radians, measured as `atan(North, East)`),
by linear interpolation between the bin centres of [`projected_limb_profile`](@ref).
"""
@inline function limb_radius(profile::AbstractVector{<:Real}, φ::Real)
    n = length(profile)
    u = (mod(φ + pi, 2pi) / (2pi)) * n - 0.5      # position in bin-centre coordinates
    b0 = floor(Int, u)
    w = u - b0
    i0 = mod(b0, n) + 1
    i1 = mod(b0 + 1, n) + 1
    return (1 - w) * profile[i0] + w * profile[i1]
end

"""
    silhouette_polygon(star; nbins=72) -> (sx, sy)

The component's projected outline as an explicit closed polygon in `(East, North)`,
**relative to its own projected centre** and wound counter-clockwise (as
`polygon_convex_clip_area` requires).

Built from [`projected_limb_profile`](@ref): vertex `k` sits at the bin-centre angle with
that bin's maximum projected vertex radius. For a detached lobe the true silhouette is
convex — the projection of a convex body always is — so this is its convex hull sampled
azimuthally; for a near-fillout lobe with an L1 cusp it degrades gracefully to the
star-shaped outline rather than over-occulting the way a strict hull would.
"""
function silhouette_polygon(star; nbins::Int = 72)
    prof = projected_limb_profile(star; nbins=nbins)
    sx = Vector{eltype(prof)}(undef, nbins); sy = similar(sx)
    @inbounds for k in 1:nbins
        φ = -pi + (k - 0.5) * 2pi / nbins
        sx[k] = prof[k] * cos(φ); sy[k] = prof[k] * sin(φ)
    end
    return sx, sy, prof
end

"""
    occultation_weights(star1, star2; method=:exact, softness=nothing, nbins=72)
        -> (w1, w2, occulting)

Soft mutual-occultation factors in `[0,1]`, one per tessel, to multiply into each
component's `vis_weights`.

`binary_cvis` otherwise sums the two components unconditionally, so the hidden face of the
far star keeps contributing flux once the disks overlap on the sky. This supplies the
missing factor, in the same differentiable sigmoid style as `soft_visibility`.

The two components are treated as disjoint convex bodies (the same assumption the
radiosity solver makes), so whichever centre has the larger toward-observer coordinate is
in front *everywhere*: only the far star is occulted, and its returned weights are

    wᵢ = σ( (dᵢ − R_limb(φᵢ)) / width )

with `dᵢ`, `φᵢ` the tessel centre's projected distance and position angle from the near
star's projected centre. `R_limb` follows the near star's actual (possibly tidally
distorted) outline via [`projected_limb_profile`](@ref).

`method`:

* `:exact` (default) — the **fraction of each far tessel's projected area** that survives,
  by clipping its quad against the silhouette polygon
  ([`polygon_convex_clip_area`](@ref)). No free parameters, and the answer no longer
  depends on the mesh: a tessel straddling the limb is partially occulted, as it should
  be. Only tessels that actually straddle pay for the clip — the rest are settled by a
  cheap per-tessel circle test — so this costs about the same as `:soft`.
* `:soft` — the tessel *centre* is either in or out, smoothed by a sigmoid of width
  `softness` (a fraction of the limb radius, defaulting to the mesh scale
  `√(4π/npix)/2`). Cheaper to differentiate through (C∞ rather than C0), but its accuracy
  is tied to the tessellation: at Spica's closest approach it gives 5.8 % of the
  secondary's flux occulted at nside 3/2 against a converged 9.8 %, a 40 % error, because
  the straddling band is a large fraction of a coarse mesh.

**Fast path:** an O(1) bounding-circle test on the projected separation against the sum of
the two maximum projected radii. When the disks are disjoint — the overwhelmingly common
case — both weight vectors are `nothing`, which callers treat as "no occultation" with no
allocation. `occulting` is `0` (none), `1` (star 1 occults star 2) or `2`.
"""
function occultation_weights(star1, star2; method::Symbol = :exact,
                             softness = nothing, nbins::Int = 72)
    method in (:exact, :soft) ||
        error("occultation_weights: method must be :exact or :soft (got $(method))")
    # Geometry follows the meshes; only the epoch (an absolute JD) is forced wide elsewhere.
    G  = float(promote_type(eltype(star1), eltype(star2)))
    o1 = G.(star1.center_offsets)
    o2 = G.(star2.center_offsets)
    # projected (East, North) centres and the maximum projected radius of each component
    c1 = (-o1[1], o1[2]); c2 = (-o2[1], o2[2])
    rad(s) = sqrt(maximum(s.proj_west .^ 2 .+ s.proj_north .^ 2))
    r1 = rad(star1); r2 = rad(star2)
    ρ = hypot(c2[1] - c1[1], c2[2] - c1[2])
    ρ >= r1 + r2 && return (nothing, nothing, 0)          # O(1) reject: disks disjoint

    # Larger toward-observer coordinate is nearer.
    near_is_1 = o1[3] > o2[3]
    near, far, cn, cf = near_is_1 ? (star1, star2, c1, c2) : (star2, star1, c2, c1)

    prof = projected_limb_profile(near; nbins=nbins)
    # The exact path clips against a true convex hull (half-plane intersection is only
    # valid for a convex clip region); the soft path uses the azimuthal profile, which
    # needs no convexity.
    sx, sy = if method === :exact
        convex_hull_2d(vec(-near.proj_west), vec(near.proj_north))
    else
        (G[], G[])
    end
    nc = length(sx)
    Rmax = maximum(prof); Rmin = minimum(prof)
    w = softness === nothing ? sqrt(4pi / near.npix) / 2 : G(softness)
    dx = cf[1] - cn[1]; dy = cf[2] - cn[2]     # far centre relative to near centre

    T = eltype(far.proj_west)
    n = far.npix
    wf = Vector{T}(undef, n)
    pw = far.proj_west; pn = far.proj_north
    # Scratch for the clipper: a 4-gon clipped by an nc-gon has ≤ 4+nc vertices.
    #
    # DELIBERATELY Float64 regardless of the mesh type, unlike everything else here. This is
    # the standard "exact predicates on inexact input" pattern: Sutherland–Hodgman decides
    # inside/outside from the sign of a cross product, and at a grazing intersection that
    # sign is the difference of two nearly equal products. In Float32 the sign can flip
    # spuriously, which does not perturb the clipped area slightly — it changes the vertex
    # COMBINATORICS and can drop or duplicate a whole corner. Widening the coordinates is
    # lossless and costs a few hundred bytes of scratch per call.
    Ax = Vector{Float64}(undef, nc + 8); Ay = similar(Ax)
    Bx = similar(Ax); By = similar(Ax)
    qx = Vector{Float64}(undef, 4); qy = Vector{Float64}(undef, 4)

    @inbounds for i in 1:n
        # the tessel's four projected corners, relative to the near star's centre
        cx = 0.0; cy = 0.0
        for v in 1:4
            qx[v] = -Float64(pw[i, v]) + dx
            qy[v] =  Float64(pn[i, v]) + dy
            cx += qx[v]; cy += qy[v]
        end
        cx /= 4; cy /= 4
        d = hypot(cx, cy)

        if method === :soft
            wf[i] = T(sigmoid((d - limb_radius(prof, atan(cy, cx))) / (w * limb_radius(prof, atan(cy, cx)))))
            continue
        end

        # cheap per-tessel accept/reject: circumradius of the quad about its centroid
        rt = 0.0
        for v in 1:4; rt = max(rt, hypot(qx[v] - cx, qy[v] - cy)); end
        if d - rt >= Rmax
            wf[i] = one(T); continue           # wholly outside the silhouette
        elseif d + rt <= Rmin
            wf[i] = zero(T); continue          # wholly inside it (Rmin disc ⊂ silhouette)
        end
        # straddling: exact clipped-area fraction
        atot = _shoelace(qx, qy, 4)
        if atot <= 0
            wf[i] = one(T)
        else
            ahid = polygon_convex_clip_area(qx, qy, 4, sx, sy, nc, Ax, Ay, Bx, By)
            wf[i] = T(clamp(1 - ahid / atot, 0, 1))
        end
    end
    ones_near = ones(T, near.npix)
    return near_is_1 ? (ones_near, wf, 1) : (wf, ones_near, 2)
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
