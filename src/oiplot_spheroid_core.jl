# =======================================================================================
# Backend-agnostic geometry for the spheroid plots
# =======================================================================================
# Everything here answers a question about the MESH — where the pole is, how the body frame
# maps to the sky, what a parallel looks like once occluded — and draws nothing. It is shared
# by the matplotlib front-end (src/oiplot_spheroid.jl) and the Makie one
# (src/oiplot_spheroid_makie.jl) so the two cannot disagree about the geometry they render.
#
# That split follows OITOOLS, which pushed the toolkit-free half of its plotting into the core
# package for the same reason (see ~/SOFTWARE/OITOOLS.jl/ext/OITOOLSMakieExt.jl:28-31). None of
# these functions touches a plotting library: the only external calls are `svd`, `norm`
# (LinearAlgebra), `mean` and `sortperm`, all of which the parent module already has.
#
# Names are underscore-prefixed and NOT exported — they are implementation detail shared
# between two front-ends, not public API.

function _polar_radius(star, star_params)
    if star_params !== nothing && hasproperty(star_params, :surface_type)
        stype = star_params.surface_type
        if stype == 0
            return star_params.radius
        elseif stype == 1
            return star_params.radius_z
        elseif stype == 2
            return star_params.rpole
        end
    end
    return maximum(sqrt.(star.vertices_xyz[:,:,1].^2 .+ star.vertices_xyz[:,:,2].^2 .+ star.vertices_xyz[:,:,3].^2))
end

"""
    _mesh_rotation(star) -> R   (3x3, body frame -> sky frame)

Recover the body->sky rotation directly from the mesh, with no angles involved.

`stellar_geometry` stores each tessel centre twice: in `vertices_spherical` as body-frame
(r, θ, φ), and in `vertices_xyz` after rotation to the sky. Those two point sets differ by
exactly the rotation we want, so an orthogonal Procrustes fit recovers it.

This matters because it captures the rotation PHASE ψ as well as inclination and position
angle. Reconstructing angles from the spin axis alone would place the pole correctly but
leave the meridians at an arbitrary phase.
"""
function _mesh_rotation(star)
    sp = star.vertices_spherical[:, 5, :]
    r, θ, φ = Float64.(sp[:, 1]), Float64.(sp[:, 2]), Float64.(sp[:, 3])
    B = hcat(r .* sin.(θ) .* cos.(φ), r .* sin.(θ) .* sin.(φ), r .* cos.(θ))
    S = Float64.(star.vertices_xyz[:, 5, :])
    F = svd(B' * S)                     # minimise ||B*R - S|| over orthogonal R
    R = F.U * F.Vt
    if det(R) < 0                       # reject a reflection, keep it a proper rotation
        R = F.U * Diagonal([1.0, 1.0, -1.0]) * F.Vt
    end
    return R
end

"""
    _spin_axis(star, star_params, inclination, position_angle) -> (north, south)

The stellar rotation axis in sky coordinates, as the two pole positions.

With `inclination`/`position_angle` given (degrees) the axis is analytic. Otherwise it is
recovered from the mesh by averaging the tessels of extreme COLATITUDE in the body frame.

Selecting on colatitude is load-bearing: do NOT take the first and last few pixels by
index. ROTIR tessellations are HEALPix NESTED, where the leading pixels lie in base pixel 0,
a diamond straddling mid-latitudes — only under RING ordering are they the pole. Indexing
puts the axis 76.5° off for an nside=3 sphere at inclination 60°/PA 30°, with the north
arrow pointing *away* from the observer. Colatitude selection is ordering-independent and
matches the analytic axis to 0.02° (the residual is Float32 mesh discretisation).
"""
function _spin_axis(star, star_params, inclination, position_angle)
    R = _polar_radius(star, star_params)
    if !isnan(inclination) && !isnan(position_angle)
        inc_rad = inclination * π / 180
        PA_rad  = position_angle * π / 180
        spin = [-sin(PA_rad) * sin(inc_rad), cos(PA_rad) * sin(inc_rad), cos(inc_rad)]
    else
        θ = star.vertices_spherical[:, 5, 2]
        k = min(4, star.npix)
        inorth = partialsortperm(θ, 1:k)              # smallest colatitude = north pole
        isouth = partialsortperm(θ, 1:k, rev = true)  # largest  colatitude = south pole
        spin = Float64.(vec(mean(star.vertices_xyz[inorth, 5, :], dims = 1)) .-
                        vec(mean(star.vertices_xyz[isouth, 5, :], dims = 1)))
        spin ./= norm(spin)
    end
    return R .* spin, -R .* spin
end

"""
    _mesh_surface_field(star) -> (dx, dy, dz, r, nz)

Sample the star's own surface as a scattered field on the body-frame unit sphere.

`stellar_geometry` stores, per tessel centre, the body-frame direction (θ, φ) and the
deformed radius `r(θ, φ)` in `vertices_spherical`, plus the SKY-frame normal in `normals`.
Together those are everything a graticule needs — the shape and where it stops being
visible — for *any* surface, including ones with no closed form.

`(dx, dy, dz)` are the body-frame unit vectors, kept as three separate vectors rather than
an `npix × 3` matrix so the neighbour search below reads contiguous memory.
"""
function _mesh_surface_field(star)
    θ = Float64.(star.vertices_spherical[:, 5, 2])
    φ = Float64.(star.vertices_spherical[:, 5, 3])
    s = sin.(θ)
    return (s .* cos.(φ), s .* sin.(φ), cos.(θ),
            Float64.(star.vertices_spherical[:, 5, 1]), Float64.(star.normals[:, 3]))
end

"""
    _interp_surface(field, θ, φ; k=8) -> (r, nz)

Interpolate the mesh surface field at an arbitrary body-frame direction.

Modified Shepard (Franke–Little) interpolation over the `k` nearest tessel centres, with
the cutoff radius `R` set by the (k+1)-th. The weight `((R-d)/(R·d))²` vanishes *smoothly*
at `d = R`, so a neighbour entering or leaving the set contributes nothing at the moment it
does. Plain inverse-distance weighting over a fixed `k` is discontinuous in its derivative
wherever the neighbour set changes, and on a graticule that shows up as faint kinks along
each curve.

Distances are squared chords `2(1 - u·v)` on the unit sphere — monotonic in angle, free of
trigonometry, and immune to the φ wrap that plagues differencing in (θ, φ).

Weights are normalised, so a constant field is reproduced exactly: a sphere interpolates to
a sphere, whatever the mesh.
"""
function _interp_surface(field, θq::Real, φq::Real; k::Int=8)
    dx, dy, dz, rs, nzs = field
    n = length(rs)
    kk = min(k, n - 1)
    # Promote FIRST. The mesh is Float32, so a query taken straight off it would be
    # evaluated in Float32 trig while the field was built in Float64 — a ~3e-8 disagreement,
    # enough to miss the exact-hit shortcut below on the very samples that define the field.
    θ = Float64(θq); φ = Float64(φq)
    st = sin(θ)
    v1 = st * cos(φ); v2 = st * sin(φ); v3 = cos(θ)
    idx = zeros(Int, kk + 1)
    dd  = fill(Inf, kk + 1)
    @inbounds for i in 1:n
        d = 2 * (1 - (dx[i] * v1 + dy[i] * v2 + dz[i] * v3))
        d < dd[kk+1] || continue
        j = kk + 1
        while j > 1 && dd[j-1] > d
            dd[j] = dd[j-1]; idx[j] = idx[j-1]; j -= 1
        end
        dd[j] = d; idx[j] = i
    end
    dd[1] <= 1e-12 && return (rs[idx[1]], nzs[idx[1]])   # query sits on a sample
    R = sqrt(max(dd[kk+1], eps()))
    wsum = 0.0; rsum = 0.0; nsum = 0.0
    @inbounds for j in 1:kk
        di = sqrt(dd[j])
        w = ((R - di) / (R * di))^2
        wsum += w; rsum += w * rs[idx[j]]; nsum += w * nzs[idx[j]]
    end
    wsum > 0 || return (rs[idx[1]], nzs[idx[1]])          # degenerate: all neighbours tied
    return (rsum / wsum, nsum / wsum)
end

"""
    _mesh_body_curve(field, θs, φs; k=8) -> (body_pts, vis)

A graticule curve read off the mesh: `npoints × 3` body-frame points, plus the sky-frame
visibility of each.

`vis` comes from the interpolated normal, not from `z > 0`. Those agree only for a sphere:
on any distorted body the silhouette is the set of points whose normal is perpendicular to
the line of sight, which is a plane through the origin tilted away from `z = 0` (2.6° for
the near-lobe-filling Roche star in the docs). Using the mesh normal makes the curves
terminate at exactly the boundary `draw_limb!` draws, since that hull is built from
`index_quads_visible`, which is thresholded on the same quantity.
"""
function _mesh_body_curve(field, θs, φs; k::Int=8)
    npts = length(θs)
    body_pts = Matrix{Float64}(undef, npts, 3)
    vis = Vector{Bool}(undef, npts)
    @inbounds for m in 1:npts
        θ = Float64(θs[m]); φ = Float64(φs[m])
        r, nz = _interp_surface(field, θ, φ; k=k)
        st, ct = sin(θ), cos(θ)
        body_pts[m, 1] = r * st * cos(φ)
        body_pts[m, 2] = r * st * sin(φ)
        body_pts[m, 3] = r * ct
        vis[m] = nz > 0
    end
    return body_pts, vis
end

"""
Extract contiguous visible segments as Nx2 matrices for PolyCollection.

Visibility defaults to the front hemisphere `z > 0`, which is exact for a sphere and a
good approximation elsewhere. Pass `vis` to override it — the mesh path supplies the sign
of the interpolated surface normal, which is the true silhouette test on a distorted body.
"""
function _visible_segments(sky_pts, offset_west, offset_north; vis=nothing)
    segments = Vector{Matrix{Float64}}()
    seg_start = 0
    npts = size(sky_pts, 1)
    for k in 1:npts
        if vis === nothing ? sky_pts[k, 3] > 0 : vis[k]
            if seg_start == 0; seg_start = k; end
        else
            if seg_start > 0 && k - seg_start >= 2
                rng = seg_start:k-1
                push!(segments, hcat(-(sky_pts[rng, 1] .+ offset_west), sky_pts[rng, 2] .+ offset_north))
            end
            seg_start = 0
        end
    end
    if seg_start > 0 && npts - seg_start + 1 >= 2
        rng = seg_start:npts
        push!(segments, hcat(-(sky_pts[rng, 1] .+ offset_west), sky_pts[rng, 2] .+ offset_north))
    end
    return segments
end

# Recover (ntheta, nphi) from a lat/long mesh. `tessellation_latlong` lays pixels out
# theta-major, phi-minor (`ilong_range = (i-1)*nphi+1 : i*nphi`), so the colatitude of the
# tessel centres is constant across the first ring and changes at pixel nphi+1. That first
# change point IS nphi, and ntheta follows. This is why the dimensions need not be stored
# on `tessellation`, which carries only `npix`.
function _latlong_dims(star)
  θ = star.vertices_spherical[:, 5, 2]
  k = findfirst(i -> !isapprox(θ[i], θ[1]; atol = 1e-6), 2:star.npix)
  nphi = k === nothing ? star.npix : k    # index into 2:npix ⇒ ring length
  ntheta, r = divrem(star.npix, nphi)
  r == 0 || error("plot_mollweide: npix=$(star.npix) is not divisible by the recovered " *
                  "ring length nphi=$nphi — mesh is not a regular lat/long grid.")
  return ntheta, nphi
end

"""
    graticule_segments(star; nlat=5, nlon=8, offset_west=0.0, offset_north=0.0,
                       inclination=NaN, position_angle=NaN, npoints=200,
                       star_params=nothing) -> Vector{Matrix{Float64}}

Parallels and meridians as visibility-clipped polylines in plot coordinates
(`x = East = -proj_west`), one `N×2` matrix per contiguous visible run.

This is the whole of the graticule calculation and it draws nothing, so both plotting
back-ends share it rather than each carrying a copy of the projection, the occlusion test
and the segment splitting.

Two paths. With `star_params` giving an analytic `surface_type` (0 sphere, 1 ellipsoid,
2 rapid rotator) the curve is evaluated in closed form. Otherwise — a Roche lobe, or any
star whose parameters were not supplied — it is interpolated from the MESH itself through
`_mesh_body_curve`, so the graticule follows the real deformed surface rather than a sphere
that only approximates it.
"""
function graticule_segments(star; nlat=5, nlon=8, offset_west=0.0, offset_north=0.0,
                            inclination=NaN, position_angle=NaN, npoints=200,
                            star_params=nothing)

    # Determine surface model
    use_exact = star_params !== nothing && hasproperty(star_params, :surface_type)
    stype = use_exact ? star_params.surface_type : -1

    # Read the surface off the mesh whenever there is no closed form to use.
    mesh_field = (!use_exact || stype ∉ (0, 1, 2)) ? _mesh_surface_field(star) : nothing

    # Build rotation matrix (same convention as rotate_star).
    #
    # Orientation precedence: explicit angles, else the MESH. `star_params` supplies the
    # SHAPE only and carries no orientation authority.
    #
    # Do not be tempted to read `inclination`/`position_angle` off `star_params`: that is
    # right for a single star but WRONG for a binary component, because
    # `create_binary_geometry` orients
    # both components by the shared `binary_frame` built from the ORBIT (i, Ω, ω), and
    # ignores each star's own inclination/position_angle entirely. For the Spica-like test
    # binary those two answers differ by 102.6°.
    #
    # The mesh is the one source that is always right, because it IS the star being drawn,
    # whatever built it — and `_mesh_rotation` recovers the rotation phase too, which no
    # pair of angles can express without also knowing the epoch and period.
    if isnan(inclination) || isnan(position_angle)
        R = _mesh_rotation(star)
    else
        # Rotation angle from the epoch (needs the period, which only star_params carries)
        rotation_angle = (star_params !== nothing && hasproperty(star_params, :rotation_period)) ?
            360.0 * star.t / star_params.rotation_period : 0.0
        R = rot_vertex(rotation_angle * π / 180,
                       inclination * π / 180, position_angle * π / 180)
    end

    # The mesh normals are SKY-frame vectors, so they only describe visibility for the
    # orientation the mesh actually has. Usually that IS `R` — either it came from the mesh,
    # or the caller passed the star's own angles, which reconstruct the same rotation. When
    # it does not (deliberately overridden angles), fall back to the z > 0 hemisphere clip
    # rather than clipping against a different star's silhouette.
    mesh_orient = mesh_field === nothing || norm(R .- _mesh_rotation(star)) < 1e-3

    graticule_lines = Vector{Matrix{Float64}}()

    # Latitude circles (constant colatitude)
    θ_targets = collect(range(π/(nlat+1), stop=π*nlat/(nlat+1), length=nlat))
    ϕ_range = collect(range(-π, stop=π, length=npoints))
    for θ0 in θ_targets
        vis = nothing
        if stype == 0
            r = star_params.radius
            body_pts = hcat(r .* sin(θ0) .* cos.(ϕ_range),
                            r .* sin(θ0) .* sin.(ϕ_range),
                            fill(r * cos(θ0), npoints))
        elseif stype == 1
            body_pts = hcat(star_params.radius_x .* sin(θ0) .* cos.(ϕ_range),
                            star_params.radius_y .* sin(θ0) .* sin.(ϕ_range),
                            fill(star_params.radius_z * cos(θ0), npoints))
        elseif stype == 2
            arg = star_params.frac_escapevel * sin(θ0)
            r0 = abs(arg) < 1e-5 ? star_params.rpole : star_params.rpole * f_rapid_rot(arg)
            body_pts = hcat(r0 .* sin(θ0) .* cos.(ϕ_range),
                            r0 .* sin(θ0) .* sin.(ϕ_range),
                            fill(r0 * cos(θ0), npoints))
        else
            body_pts, vis = _mesh_body_curve(mesh_field, fill(θ0, npoints), ϕ_range)
            mesh_orient || (vis = nothing)
        end
        sky_pts = body_pts * R
        append!(graticule_lines, _visible_segments(sky_pts, offset_west, offset_north; vis=vis))
    end

    # Longitude lines (constant azimuth)
    ϕ_targets = collect(range(0, stop=2π*(1 - 1/nlon), length=nlon))
    θ_range = collect(range(0, stop=π, length=npoints))
    for ϕ0 in ϕ_targets
        vis = nothing
        if stype == 0
            r = star_params.radius
            body_pts = hcat(r .* sin.(θ_range) .* cos(ϕ0),
                            r .* sin.(θ_range) .* sin(ϕ0),
                            r .* cos.(θ_range))
        elseif stype == 1
            body_pts = hcat(star_params.radius_x .* sin.(θ_range) .* cos(ϕ0),
                            star_params.radius_y .* sin.(θ_range) .* sin(ϕ0),
                            star_params.radius_z .* cos.(θ_range))
        elseif stype == 2
            ω = star_params.frac_escapevel
            args = ω .* sin.(θ_range)
            r_vals = star_params.rpole .* f_rapid_rot.(args)
            r_vals[abs.(args) .< 1e-5] .= star_params.rpole
            body_pts = hcat(r_vals .* sin.(θ_range) .* cos(ϕ0),
                            r_vals .* sin.(θ_range) .* sin(ϕ0),
                            r_vals .* cos.(θ_range))
        else
            body_pts, vis = _mesh_body_curve(mesh_field, θ_range, fill(ϕ0, npoints))
            mesh_orient || (vis = nothing)
        end
        sky_pts = body_pts * R
        append!(graticule_lines, _visible_segments(sky_pts, offset_west, offset_north; vis=vis))
    end

    return graticule_lines
end

# ── shared with BOTH plotting back-ends ──────────────────────────────────────────────────
#
# These live here rather than beside the code that uses them because a sibling extension
# cannot import from another one: a helper defined in ROTIRMakieExt is invisible to
# ROTIRGUIExt however it is imported, and duplicating it is how the two come to disagree.

"""
    _map_range(values; vmin = nothing, vmax = nothing) -> (lo, hi)

The colour range for a map, with the guard both back-ends need: a uniform map gives
`lo == hi`, which Makie rejects as a colorrange and matplotlib renders as a single flat
colour with no indication that the scale collapsed.
"""
function _map_range(projmap; vmin = nothing, vmax = nothing)
    pmin = vmin === nothing ? minimum(projmap) : vmin
    pmax = vmax === nothing ? maximum(projmap) : vmax
    if pmax - pmin < 1.0
        # SYMMETRIC about the value, not open-ended upwards. Widening only the top puts every
        # value on the floor of the colormap, so a uniform-temperature star — a legitimate
        # model, and what a non-rotating rapid rotator IS — renders solid black and reads as a
        # failed draw rather than as a flat map.
        mid = (pmin + pmax) / 2
        half = max(abs(mid) * 0.005, 0.5)
        pmin, pmax = mid - half, mid + half
    end
    return pmin, pmax
end

"The largest distance from the body centre to any vertex, used to frame a 3-D view."
_body_radius(star) = maximum(sqrt.(sum(abs2, star.vertices_xyz[:, 1:4, :], dims = 3)))
