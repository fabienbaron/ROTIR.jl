using PythonPlot, PythonCall, LaTeXStrings, Statistics

# Suppress matplotlib warning about tight_layout + equal-aspect axes
# (deferred to first plot call via set_oiplot_defaults)

const _oiplot_initialized = Ref(false)

"""
    set_oiplot_defaults(; compact=<current>)

Apply ROTIR's matplotlib style. This **delegates to `OITOOLS.set_oiplot_defaults`** so a
figure produced by ROTIR and one produced by OITOOLS look the same in the same document —
they routinely appear side by side in a paper, and a serif/sans or 12pt/14pt mismatch
between them is immediately visible.

Delegating to OITOOLS keeps a single source of truth for the rcParams — a private copy
drifts (font size, minor tick count, `markeredgewidth`) and the two packages' figures stop
matching. It also picks up OITOOLS' `compact` mode for free, worth passing
`compact=true` when stacking panels.

Every plotting entry point in this file calls this first. That matters for more than
style: the rcParams are global, so a function that does *not* set them inherits whatever
the last caller left behind, and its font then depends on what you happened to plot
earlier in the session.
"""
function set_oiplot_defaults(; kwargs...)
    if !_oiplot_initialized[]
        pyimport("warnings").filterwarnings("ignore", message=".*tight_layout.*")
        _oiplot_initialized[] = true
    end
    OITOOLS.set_oiplot_defaults(; kwargs...)
end

const global oiplot_colors=["black", "gold","chartreuse","blue","red", "pink","lightgray","darkorange","darkgreen","aqua",
"fuchsia","saddlebrown","dimgray","darkslateblue","violet","indigo","blue","dodgerblue",
"sienna","olive","purple","darkorchid","tomato","darkturquoise","steelblue","seagreen","darkgoldenrod","darkseagreen","salmon","slategray","lime","coral","maroon","mistyrose","sandybrown","tan","olivedrab"]

const global oiplot_markers=["o","s","v","P","*","x","^","D","p",1,"<","H","X","4",4,"_","1",6,"8","d",9]

############################################################
#
# Imaging on spheroids
#
############################################################

"""
Physical width of `ax` in inches, from the axes rectangle and the figure size. Uses
`get_position` rather than `get_window_extent` so no renderer is needed — this is called
while the figure is still being built.
"""
function _axes_width_inches(ax)
    pos = ax.get_position()
    return pyconvert(Float64, pos.width) *
           pyconvert(Float64, ax.figure.get_size_inches()[0])
end

function _axes_height_inches(ax)
    pos = ax.get_position()
    return pyconvert(Float64, pos.height) *
           pyconvert(Float64, ax.figure.get_size_inches()[1])
end

"""
    draw_compass(ax, axis_max; size_frac=0.12, fontsize=nothing, color="black")

Draw E/N compass arrows in the lower-right corner of a 2D sky-plane plot.
East points left (astronomical convention when x-axis is inverted).

`fontsize=nothing` scales the labels with the axes' physical width, so the compass stays
legible in a small multi-panel subplot instead of colliding with its own arrows. Pass a
number to pin it.
"""
function draw_compass(ax, axis_max; size_frac=0.12, fontsize=nothing, color="black")
    # The compass is a figure annotation, like a scale bar: it should stay a roughly FIXED
    # PHYSICAL size rather than shrink with the panel. Sizing the arrows in data units
    # (`size_frac * axis_max`) while the labels are in points is what makes the compass
    # illegible in plot2d_allepochs, where three panels share one figure — the arrow shrinks
    # 3x but the text cannot, so N and E collide with their own arrows.
    #
    # So: derive the arrow length in POINTS from size_frac, clamp it to a legible physical
    # range, then convert back to data units. A large single-panel plot lands mid-range and
    # is unchanged; a small subplot gets a proportionally larger — and readable — compass.
    h_pts     = max(_axes_height_inches(ax) * 72, eps())
    arrow_pts = clamp(0.5 * size_frac * h_pts, 22.0, 40.0)
    s         = (arrow_pts / h_pts) * (2 * axis_max)
    margin    = 0.20 * axis_max
    cx = -(axis_max - margin)
    cy = -(axis_max - margin)
    if fontsize === nothing
        fontsize = clamp(0.34 * arrow_pts, 7.0, 13.0)
    end
    # Offset the labels by their own size in DATA units, so the gap tracks the text.
    gap = 0.75 * fontsize * (2 * axis_max) / h_pts
    ax.annotate("N", xy=(cx - gap, cy + s), xytext=(cx - gap, cy),
        arrowprops=Dict("arrowstyle" => "-|>", "lw" => 1.5, "color" => color),
        fontsize=fontsize, color=color, ha="center", va="bottom", fontweight="bold", zorder=8)
    ax.annotate("E", xy=(cx + s, cy - gap), xytext=(cx, cy - gap),
        arrowprops=Dict("arrowstyle" => "-|>", "lw" => 1.5, "color" => color),
        fontsize=fontsize, color=color, ha="right", va="center", fontweight="bold", zorder=8)
end

"""
    _padded_norm(lo, hi, range; cfloor=0.15, cceil=0.95) -> matplotlib Normalize

Map `[lo, hi]` onto the colormap fraction `[cfloor, cceil]` rather than `[0, 1]`, by
widening the normalisation limits at both ends.

Both ends matter. The default background is WHITE and the top of `gist_heat` is pure
white, so without a ceiling the hottest part of a surface renders in exactly the background
colour. On a binary that is fatal: the two components share one scale, so the hotter one
sits at the top of the range and its disk disappears completely. Padding only the floor
("avoiding black pixels on a dark background") guards the end that is not the problem.
"""
function _padded_norm(lo, hi, range; cfloor = 0.15, cceil = 0.95)
    span = cceil - cfloor
    return matplotlib.colors.Normalize(vmin = lo - cfloor       / span * range,
                                       vmax = hi + (1 - cceil)  / span * range)
end

"""
    draw_limb!(ax, star; offset_west=0.0, offset_north=0.0, color="black",
               linewidth=0.8, alpha=1.0, zorder=4)

Draw the projected stellar limb — the silhouette of the tessellated surface on the sky.

This is what keeps a star visible regardless of its colours. A surface that maps to the
pale end of a colormap (a hot component on `gist_heat`) is indistinguishable from a white
background, so the disk vanishes and any decoration drawn on it — rotation axis, spin
arrow — appears to float in empty space with nothing to reference it against.

The outline is the CONVEX HULL of the projected vertices, not `silhouette_polygon`. Both
describe the same boundary, but the latter bins vertices by azimuth and keeps the largest
radius per bin: wherever no vertex happens to fall near the true edge of its bin the radius
under-estimates, and the outline notches inward. On a coarse mesh those notches are plainly
visible. The hull is exact for a convex body and has no resolution parameter at all.

(`silhouette_polygon` remains the right tool for occultation, where a radius-vs-azimuth
profile is what the calculation needs.)

Vertices are taken in `x = -proj_west`, the same convention as the plots, so only the
offsets need applying.
"""
function draw_limb!(ax, star; offset_west=0.0, offset_north=0.0, color="black",
                    linewidth=0.8, alpha=1.0, zorder=4)
    vis = star.index_quads_visible
    hx, hy = convex_hull_2d(vec(-star.proj_west[vis, :]), vec(star.proj_north[vis, :]))
    ax.plot(Float64.(vcat(hx, hx[1])) .- offset_west,
            Float64.(vcat(hy, hy[1])) .+ offset_north,
            "-", color=color, linewidth=linewidth, alpha=alpha, zorder=zorder)
end

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
    draw_rotation_axis(ax, star; inclination=NaN, position_angle=NaN,
        arrow_frac=0.3, color="black", linewidth=1.5, offset_west=0.0, offset_north=0.0)

Draw the projected stellar rotation axis on a 2D sky-plane plot.
Shows the pole-to-pole line (dashed) extended beyond the limb, with an arrow at the north pole.

When `inclination` (degrees from LOS) and `position_angle` (degrees, N through E) are given,
the axis is computed analytically. Otherwise it is estimated from the tessellation vertices.
"""
function draw_rotation_axis(ax, star; arrow_frac=0.3, color="black", linewidth=1.5,
    offset_west=0.0, offset_north=0.0, inclination=NaN, position_angle=NaN,
    alpha=1.0, inside_alpha=0.5, star_params=nothing)
    north, south = _spin_axis(star, star_params, inclination, position_angle)
    delta = north .- south
    north_tip = north .+ arrow_frac .* delta
    south_tip = south .- arrow_frac .* delta
    points = hcat(south_tip, south, north, north_tip)  # 3×4
    seg_alpha = [alpha, inside_alpha, alpha]  # outside, inside, outside
    for j in 1:3
        p1 = points[:, j]
        p2 = points[:, j+1]
        ax.plot([-(p1[1]+offset_west), -(p2[1]+offset_west)],
                [p1[2]+offset_north, p2[2]+offset_north],
                "--", color=color, linewidth=linewidth, alpha=seg_alpha[j], zorder=6)
    end
    ax.annotate("", xy=(-(north_tip[1]+offset_west), north_tip[2]+offset_north),
        xytext=(-(north[1]+offset_west), north[2]+offset_north),
        arrowprops=Dict("arrowstyle" => "-|>", "shrinkA" => 0, "shrinkB" => 0,
            "color" => color, "lw" => 0, "mutation_scale" => 15),
        zorder=6)
end

"""
    draw_rotation_arrow(ax, star; pole="N", radius_frac=0.15, offset_frac=0.05,
        color="black", linewidth=1.5, inclination=NaN, position_angle=NaN)

Draw a curved arrow around the rotation axis showing the sense of rotation.
A 290° ellipse centered at the pole, solid in front (z>0), dashed behind.
Prograde rotation = counter-clockwise around north pole.
When `inclination`/`position_angle` are given, the axis is computed analytically.
"""
function draw_rotation_arrow(ax, star; pole="N", radius_frac=0.07, offset_frac=0.15,
    color="black", linewidth=1.5, npoints=200, offset_west=0.0, offset_north=0.0,
    inclination=NaN, position_angle=NaN, star_params=nothing)
    north, south = _spin_axis(star, star_params, inclination, position_angle)
    axis = north .- south
    axis_len = norm(axis)
    axis_hat = axis ./ axis_len
    center = pole == "N" ? north .+ offset_frac .* axis : south .- offset_frac .* axis
    r = radius_frac * axis_len
    ref = abs(axis_hat[3]) < 0.9 ? [0.0, 0.0, 1.0] : [1.0, 0.0, 0.0]
    e1 = cross(axis_hat, ref); e1 ./= norm(e1)
    e2 = cross(axis_hat, e1)
    # Reverse sweep for retrograde rotation (negative rotation_period)
    if star_params !== nothing && hasproperty(star_params, :rotation_period) && star_params.rotation_period < 0
        e1, e2 = e2, e1
    end
    θ = collect(range(0, 300π/180, length=npoints))
    pts = hcat([center .+ r .* (cos(t) .* e1 .+ sin(t) .* e2) for t in θ]...)
    x2d = -(pts[1, :] .+ offset_west); y2d = pts[2, :] .+ offset_north; z = pts[3, :]
    front = z .> 0
    xf = copy(x2d); xf[.!front] .= NaN; yf = copy(y2d); yf[.!front] .= NaN
    ax.plot(xf, yf, "-", color=color, linewidth=linewidth, zorder=7)
    xb = copy(x2d); xb[front] .= NaN; yb = copy(y2d); yb[front] .= NaN
    ax.plot(xb, yb, "--", color=color, linewidth=linewidth, zorder=1)
    zord_tip = z[end] > 0 ? 7 : 1
    # The head is symmetric about the line xytext -> xy, so that line must be the curve's
    # TANGENT at the endpoint, not a chord back to an earlier sample. The chord to
    # `end-3` spans 4.5 deg of arc and sits 1.37 deg off the tangent, which rotates the
    # head slightly against the curve it terminates. The curve is parametric, so use the
    # analytic derivative d/dt [cos t*e1 + sin t*e2] and project it the same way as pts.
    dpt = -sin(θ[end]) .* e1 .+ cos(θ[end]) .* e2
    tx, ty = -dpt[1], dpt[2]
    tnorm = hypot(tx, ty)
    tx /= tnorm; ty /= tnorm
    # Size the head against the LOOP, not at a fixed 15 points. The loop's radius `r` is in
    # data units and comes out around 15 pt at ordinary figure scales — so a fixed 15 pt
    # head is as long as the circle's radius and chords straight across the arc. The head
    # is a straight triangle, so spanning ~60 deg of a tight curve makes it look attached
    # at the wrong angle even when its axis is exactly tangent (which it is).
    ylo, yhi = pyconvert(NTuple{2,Float64}, ax.get_ylim())
    h_pts       = max(_axes_height_inches(ax) * 72, eps())
    data_per_pt = abs(yhi - ylo) / h_pts
    # matplotlib's "-|>" draws a head whose LENGTH is 0.4 * mutation_scale points, so pick
    # the length first (a fraction of the loop radius) and back out the scale.
    head_len    = 0.45 * r                                    # data units
    head_pts    = clamp(head_len / data_per_pt / 0.4, 6.0, 30.0)
    # Two settings make the head sit ON the arc rather than beside it:
    #
    #   lw = 0      — the annotate must contribute the HEAD ONLY. With a non-zero width it
    #                 also strokes its own shaft, and that shaft is the STRAIGHT segment
    #                 xytext -> xy. A straight chord against a curving arc leaves a kink at
    #                 the join, which is what made the arrow look unstitched. The plotted
    #                 curve is already the shaft, so the annotate need not draw one.
    #   shrink = 0  — the 2-point default pulls the head back off the curve's endpoint,
    #                 leaving a gap.
    #
    # draw_rotation_axis uses the same pair for the same reason.
    # "-|>" places the arrow's TIP at `xy`. Pointing `xy` at the curve's last point
    # therefore lays the head BACKWARD along the arc, burying it in the line it is meant
    # to terminate. The base belongs on the endpoint and the head should continue past it,
    # so put `xytext` on the curve and `xy` one head-length further along the tangent.
    ax.annotate("", xy=(x2d[end] + head_len * tx, y2d[end] + head_len * ty),
        xytext=(x2d[end], y2d[end]),
        arrowprops=Dict("arrowstyle" => "-|>", "shrinkA" => 0, "shrinkB" => 0,
            "color" => color, "lw" => 0, "mutation_scale" => head_pts),
        zorder=zord_tip)
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
    draw_graticules(ax, star; star_params=nothing, nlat=5, nlon=8, inclination=NaN, position_angle=NaN, ...)

Draw latitude/longitude graticule lines on the star surface using parametric curves.
Generates smooth curves in the body frame, rotates them with the same Euler rotation
as the star (rot_vertex), and clips to the visible side.
Renders via matplotlib PolyCollection for efficiency.

When `star_params` is provided and the surface has a closed form, that form is used:
- Type 0 (sphere): radius
- Type 1 (triaxial ellipsoid): radius_x, radius_y, radius_z
- Type 2 (rapid rotator): rpole, frac_escapevel via f_rapid_rot

Everything else — Roche lobes (type 3), and any star drawn without `star_params` — is
interpolated from the MESH's own `r(θ, φ)`; see `_interp_surface`.

That path exists because a Roche lobe has no axisymmetric description. The previous
fallback fitted a biaxial ellipsoid to the mesh-averaged polar and equatorial radii, which
is symmetric about the spin axis and therefore structurally unable to show the tidal
teardrop pointing at the companion — the one feature such a figure is drawn to show. The
mesh needs no `D`, no `q`, and no potential solve, and stays correct for surface types that
do not exist yet.
"""
function draw_graticules(ax, star; nlat=5, nlon=8, color="black", linewidth=0.8, alpha=0.5,
    offset_west=0.0, offset_north=0.0, inclination=NaN, position_angle=NaN,
    npoints=200, star_params=nothing, limb=true, zorder=5)
    collections = pyimport("matplotlib.collections")

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

    # The projected limb. Without it the meridians and parallels float with no boundary,
    # and on a filled plot the silhouette is invisible whenever the surface maps to the
    # pale end of the colormap (a hot component on gist_heat against a white background
    # disappears entirely). `silhouette_polygon` already works in x = -proj_west, the same
    # convention used here, so its output needs only the offsets applied.
    limb && draw_limb!(ax, star; offset_west=offset_west, offset_north=offset_north,
                       color=color, linewidth=linewidth, alpha=alpha, zorder=zorder-1)

    if !isempty(graticule_lines)
        ax.add_collection(collections.PolyCollection(graticule_lines, closed=false,
            ec=color, fc="none", linewidths=[linewidth], alpha=alpha, zorder=zorder))
    end
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

function set_tick_spacing(ax, axis_max)
  if (ceil(axis_max) <= 3.0)
    long_tick = 1.0; short_tick = 0.1;
  elseif (ceil(axis_max) <= 6.0)
    long_tick = 2.0; short_tick = 0.2;
  else
    long_tick = 3.0; short_tick = 0.5;
  end
  ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(long_tick));
  ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(long_tick));
  ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(short_tick));
  ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(short_tick));
end

"""
    rgba(cmap, x)                 -> NTuple{4,Float64}
    rgba(cmap, xs::AbstractArray) -> Vector{Vector{Float64}}

Evaluate a matplotlib colormap and bring the result back to Julia as something matplotlib
will accept.

Two PythonCall behaviours bite here:

1. **No automatic conversion.** `cmap(x)` returns a `Py`, so `cmap.(xs)` gives a
   `Vector{Py}`, which matplotlib turns into a numpy *object* array rather than a float
   array. Colour lookup then fails, or silently produces a dtype nothing downstream can use.
2. **A Julia `Matrix` is not an `n × 4` RGBA array.** Passing one gives
   `ValueError: Invalid RGBA argument: 0.1` — PythonCall hands numpy a Fortran-ordered view
   and `to_rgba_array` reads it as a sequence of scalars. Transposing does not help.

A `Vector` of length-4 `Vector`s is unambiguous under both conventions, so that is what this
returns. (`Vector{NTuple{4,Float64}}` also works for colours, but a vector of tuples can
convert to a numpy *structured* array in other contexts — see the `_xyz_rows` note below —
so the nested-vector form is the safer habit.)
"""
rgba(cmap, x::Real) = pyconvert(NTuple{4,Float64}, cmap(x))
rgba(cmap, xs::AbstractArray) =
    [collect(pyconvert(NTuple{4,Float64}, cmap(x))) for x in xs]

"""
    add_tessel_collection!(ax, star, colours; plotmesh=false, zorder=2,
                           offset_west=0.0, offset_north=0.0)

Draw a star's visible tessels as a single `matplotlib.collections.PolyCollection`.

One collection instead of one `add_patch` per tessel: at 3072 tessels that is the
difference between a plot that takes seconds and one that takes milliseconds, which is
what makes the per-frame rendering in `animation.jl` practical.

`colours` must be indexed by *position within* `star.index_quads_visible` (i.e. the same
ordering as `tmap[star.index_quads_visible]`).
"""
function add_tessel_collection!(ax, star, colours; plotmesh=false, zorder=2,
                                offset_west=0.0, offset_north=0.0,
                                edgecolors=nothing, linewidths=nothing,
                                indices=nothing)
  collections = pyimport("matplotlib.collections")
  # `indices` defaults to the visible (front-facing) tessels. Passing the complement draws
  # the far hemisphere instead, which is what a see-through wireframe needs.
  vis = indices === nothing ? star.index_quads_visible : indices
  verts = Vector{Matrix{Float64}}(undef, length(vis))
  @inbounds for (i, idx) in enumerate(vis)
    # plot coords: x = East = -West
    verts[i] = hcat(-(Float64.(star.proj_west[idx, :]) .+ offset_west),
                     Float64.(star.proj_north[idx, :]) .+ offset_north)
  end
  # Non-mesh mode strokes each polygon in its own face colour: a zero-width edge leaves
  # antialiasing seams that read as a spurious grid.
  # Mid grey for the mesh, not "lightgrey": it has to read against BOTH ends of the
  # colormap. A light edge vanishes against a correctly-lit pale surface and
  # `plotmesh=true` then silently does nothing. 0.45 grey holds on both.
  ec = edgecolors === nothing ? (plotmesh ? "0.45" : colours) : edgecolors
  lw = linewidths === nothing ? (plotmesh ? 0.25 : 0.35)      : linewidths
  pc = collections.PolyCollection(verts, facecolors=colours, edgecolors=ec,
                                  linewidths=lw, rasterized=false, zorder=zorder)
  ax.add_collection(pc)
  return pc
end

function plot2Dquad(star,i) # plots the ith quad projected onto the imaging plane
  proj_west = star.proj_west;
  proj_north = star.proj_north;
  #plots the nth quad in the 2D plane, using ABCD
  # this can be used to debug lots of stuff...
  fig = figure("Test counter",figsize=(10,10),facecolor="White");
  scatter(-proj_west[i,:], proj_north[i,:]);  # -proj_west = East (astronomical convention)
  annotate("A", xy=[-proj_west[i,1];proj_north[i,1]], xycoords="data");
  annotate("B", xy=[-proj_west[i,2];proj_north[i,2]], xycoords="data");
  annotate("C", xy=[-proj_west[i,3];proj_north[i,3]], xycoords="data");
  annotate("D", xy=[-proj_west[i,4];proj_north[i,4]], xycoords="data");
  pyplot.draw()
  return 1
end

function plot3d(star_temperature_map,star) # this plots the temperature map
  set_oiplot_defaults()
  # 3D view in ROTIR's internal frame: x₁=West on sky, x₂=North on sky, x₃=toward observer
  corners_xyz = star.vertices_xyz[:,1:4,:];
  Art3D = pyimport("mpl_toolkits.mplot3d.art3d")
  Poly3DCollection = Art3D.Poly3DCollection
  fig2 = figure("Spheroid plot",figsize=(10,10),facecolor="White");
  ax = subplot(projection="3d")
  xlabel("West (mas)"); ylabel("North (mas)"); zlabel("toward obs.");
  axis_max = maximum(sqrt.(star.vertices_xyz[:,:,1].^2 .+ star.vertices_xyz[:,:,2].^2 .+ star.vertices_xyz[:,:,3].^2))*1.5;
  xlim([axis_max,-axis_max]);
  ylim([-axis_max,axis_max]);
  zlim(bottom=-axis_max,top=axis_max);

  # One Poly3DCollection for the whole surface rather than one per tessel, and the map
  # normalisation hoisted out of the loop (it was recomputing minimum/maximum of the whole
  # temperature map for every pixel — O(npix²)).
  #
  # `_xyz_rows` note: the vertices must be a nested Vector of length-3 Vectors, NOT
  # `collect(zip(x,y,z))`. A `Vector{NTuple{3,Float32}}` converts to a numpy STRUCTURED array
  # (`dtype([('f0','<f4'),('f1','<f4'),('f2','<f4')])`) under PythonCall, and matplotlib then
  # fails with "Cannot cast array data ... to dtype('float64')".
  tmin, tmax = extrema(star_temperature_map)
  trange = tmax > tmin ? tmax - tmin : one(tmax)
  cmap = get_cmap("gist_heat")
  # Poly3DCollection needs an EXPLICIT pylist nesting — unlike the 2-D PolyCollection, which
  # accepts Julia arrays. It calls `np.asarray(verts, copy=None)` on anything array-like, and
  # `copy=None` is a numpy-2 API, so a Julia nested Vector (or a numpy array) fails with
  # "NoneType copy mode not allowed" while a genuine Python list takes a different code path.
  faces = pylist([pylist([pylist([Float64(corners_xyz[i,k,1]), Float64(corners_xyz[i,k,2]),
                                  Float64(corners_xyz[i,k,3])])
                          for k in axes(corners_xyz, 2)]) for i in 1:star.npix])
  colours = pylist([pylist(c) for c in
                    rgba(cmap, [(star_temperature_map[i] - tmin)/trange for i in 1:star.npix])])
  ax.add_collection3d(Poly3DCollection(faces, edgecolor="none", facecolor=colours))
  ax.set_aspect("equal")
  pyplot.draw()
end


"""
    plot2d(tmap, star; intensity=false, intensity_model=:linear, band=nothing,
           vmin=nothing, vmax=nothing, kwargs...) -> (fig, ax)

Plot a temperature map on the projected sky plane (East left, North up).

`intensity = false` colours by temperature. `intensity = true` colours by apparent
surface brightness, i.e. `I(T) × ldmap`, where `I` is selected by `intensity_model`:
`:linear` (`I = T`, the Rayleigh–Jeans proxy) or `:planck` at wavelength `band` (metres).
Note the two keywords are different things — `intensity` decides *whether* limb darkening
and the brightness conversion are applied at all, `intensity_model` decides *which*
conversion.

`vmin` / `vmax` pin the colour scale. Leave them `nothing` for the usual per-call
autoscale; set them when rendering a sequence of frames, or a change of scale between
frames will masquerade as a change in the star.
"""
function plot2d(tmap, star; intensity = false, figtitle ="", plotmesh=false, pad = 0.5,
    colormap="gist_heat", xlim=Float64[], ylim=Float64[], background="white", flipx=false,
    compass=true, rotation_axis=false, rotation_arrow=false, graticules=false,
    contours=Float64[], contour_color="gray", contour_labels=true, contour_fontsize=10,
    limb=true, limb_color="black", limb_linewidth=0.8,
    inclination=NaN, position_angle=NaN, star_params=nothing,
    intensity_model::Symbol = :linear, band = nothing,
    vmin=nothing, vmax=nothing,
    graticule_kwargs=(;))
  # Plot temperature map onto the projected 2D image plane (= observer view)
  # Convention: East left, North up (astronomical standard)
  set_oiplot_defaults()
  axdiv= pyimport("mpl_toolkits.axes_grid1.axes_divider")
  facecolor="White"
  if background=="black"
    facecolor="Black"
  end
  fig = figure("Epoch image",figsize=(10,10),facecolor="White")
  clf();
  ax = gca();
  title(figtitle)
  ax.set_facecolor(facecolor)
  ax.set_aspect("equal", adjustable="box")
  axis_max = maximum(sqrt.(star.vertices_xyz[:,:,1].^2 .+ star.vertices_xyz[:,:,2].^2 .+ star.vertices_xyz[:,:,3].^2))+pad;
  if flipx==false
    ax.set_xlim([axis_max,-axis_max]);
  else
    ax.set_xlim([-axis_max,axis_max]);
  end
  ax.set_ylim([-axis_max,axis_max]);
  if intensity == true
    Imap = intensity_model === :linear ? tmap : ROTIR.intensity(tmap, intensity_model, band)
    projmap = Imap[star.index_quads_visible] .* star.ldmap[star.index_quads_visible]
  else
    projmap = tmap[star.index_quads_visible];
  end
  pmin = vmin === nothing ? minimum(projmap) : vmin
  pmax = vmax === nothing ? maximum(projmap) : vmax
  prange = pmax - pmin
  if prange < 1.0; prange = max(abs(pmax) * 0.01, 1.0); end
  norm_plot = _padded_norm(pmin, pmax, prange; cfloor = 0.08)
  colours = rgba(get_cmap(colormap), norm_plot.(projmap))

  add_tessel_collection!(ax, star, colours; plotmesh=plotmesh, zorder=2)
  visible = star.index_quads_visible
  xlabel(L"x $\leftarrow$ E (mas)", fontsize=20)
  ylabel(L"y $\rightarrow$ N (mas)", fontsize=20)
  # Contours from triangulated visible polygon centers
  if !isempty(contours)
    cx = vec(mean(-star.proj_west[visible, :], dims=2))
    cy = vec(mean(star.proj_north[visible, :], dims=2))
    cs = ax.tricontour(cx, cy, Float64.(projmap), levels=Float64.(sort(contours)), colors=contour_color, zorder=4)
    if contour_labels
      ax.clabel(cs, inline=true, fontsize=contour_fontsize, fmt="%.0f K", colors=contour_color)
    end
  end
  # Decorations: graticules (z=5) < pole line (z=6) < spin arrow (z=7) < compass (z=8)
  # Limb first: it is what keeps the disk visible when the surface maps to the pale end
  # of the colormap. draw_graticules would also draw it, so suppress the duplicate.
  if limb && !graticules; draw_limb!(ax, star; color=limb_color, linewidth=limb_linewidth); end
  if graticules; draw_graticules(ax, star; inclination=inclination, position_angle=position_angle, star_params=star_params, limb=limb, color=limb_color, graticule_kwargs...); end
  if rotation_axis; draw_rotation_axis(ax, star, inclination=inclination, position_angle=position_angle, star_params=star_params); end
  if rotation_arrow; draw_rotation_arrow(ax, star, inclination=inclination, position_angle=position_angle, star_params=star_params); end
  if compass; draw_compass(ax, axis_max); end
  cmap=ColorMap(colormap)
  divider = axdiv.make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.07)
  cb=colorbar(matplotlib.cm.ScalarMappable(norm=norm_plot,cmap=cmap), cax=cax)
  if pmin < pmax; cb.ax.set_ylim(pmin, pmax); end
  return fig, ax
  end

"""
    plot2d_binary(tmap1, tmap2, star1, star2, bparams, tepoch;
                  intensity=false, intensity_model=:linear, band=nothing,
                  vmin=nothing, vmax=nothing, ax=nothing, colorbar_on=true, ...) -> (fig, ax)

Plot a binary system on the 2D sky plane with correct occlusion ordering and orbital
offset. `tepoch` is the observation time in JD; the secondary's offset relative to the
primary comes from the orbital elements in `bparams`. The farther star (larger z =
receding) is drawn behind the nearer one.

`intensity`, `intensity_model` and `band` work exactly as in [`plot2d`](@ref). The
intensity model matters more here than for a single star: it sets the *relative*
brightness of two components at different temperatures, and the linear proxy misstates
that ratio by ~4 % in H and ~13 % in V.

`vmin` / `vmax` pin the shared colour scale across both components — essential when
rendering a frame sequence. Pass an existing `ax` (and `colorbar_on=false`) to draw into
a multi-panel figure instead of creating a new one.
"""
function plot2d_binary(tmap1, tmap2, star1, star2, bparams, tepoch;
    intensity=false, plotmesh=false, colormap="gist_heat", pad=1.0, background="white",
    compass=true, rotation_axis=false, rotation_arrow=false, graticules=false, figtitle="",
    inclination1=NaN, position_angle1=NaN, inclination2=NaN, position_angle2=NaN,
    star_params1=nothing, star_params2=nothing,
    limb=true, limb_color="black", limb_linewidth=0.8,
    intensity_model::Symbol = :linear, band = nothing,
    vmin=nothing, vmax=nothing, ax=nothing, colorbar_on=true, axis_max=nothing,
    graticule_kwargs=(;))
  set_oiplot_defaults()
  axdiv = pyimport("mpl_toolkits.axes_grid1.axes_divider")
  facecolor = background == "black" ? "Black" : "White"
  if ax === nothing
    fig = figure("Binary epoch", figsize=(10,10), facecolor="White")
    clf()
    ax = gca()
  else
    fig = ax.get_figure()
  end
  ax.set_title(figtitle)
  ax.set_facecolor(facecolor)
  ax.set_aspect("equal", adjustable="box")
  # Compute orbital offset of secondary relative to primary (West, North) in mas
  offset_west, offset_north = orbit_to_rotir_offset(bparams, tepoch)
  # Axis limits encompassing both stars at their orbital positions
  # Star 1 at origin, star 2 offset by (offset_west, offset_north)
  # Plot coords: x = -proj_west (East left), y = proj_north (North up)
  r1 = maximum(sqrt.(star1.vertices_xyz[:,:,1].^2 .+ star1.vertices_xyz[:,:,2].^2))
  r2 = maximum(sqrt.(star2.vertices_xyz[:,:,1].^2 .+ star2.vertices_xyz[:,:,2].^2))
  east_offset  = -offset_west  # East = -West
  north_offset = offset_north
  if axis_max === nothing
    axis_max = max(r1, r2 + abs(east_offset), r2 + abs(north_offset),
                   abs(east_offset) + r2, abs(north_offset) + r2) + pad
  end
  ax.set_xlim([axis_max, -axis_max])
  ax.set_ylim([-axis_max, axis_max])
  # Shared color normalization across both stars
  # Pad the bottom of the range so the coolest temperature maps to ~0.15
  # instead of 0.0 (pure black in gist_heat), ensuring both stars are visible.
  if intensity
    I1 = intensity_model === :linear ? tmap1 : ROTIR.intensity(tmap1, intensity_model, band)
    I2 = intensity_model === :linear ? tmap2 : ROTIR.intensity(tmap2, intensity_model, band)
    projmap1 = I1[star1.index_quads_visible] .* star1.ldmap[star1.index_quads_visible]
    projmap2 = I2[star2.index_quads_visible] .* star2.ldmap[star2.index_quads_visible]
  else
    projmap1 = tmap1[star1.index_quads_visible]
    projmap2 = tmap2[star2.index_quads_visible]
  end
  tmin = vmin === nothing ? min(minimum(projmap1), minimum(projmap2)) : vmin
  tmax = vmax === nothing ? max(maximum(projmap1), maximum(projmap2)) : vmax
  trange = tmax - tmin
  if trange < 1.0; trange = max(tmax * 0.01, 1.0); end
  norm_plot = _padded_norm(tmin, tmax, trange; cfloor = 0.15)
  colours1 = rgba(get_cmap(colormap), norm_plot.(projmap1))
  colours2 = rgba(get_cmap(colormap), norm_plot.(projmap2))
  # Determine z-ordering: farther star (larger z = receding) drawn first (behind)
  _, _, z1, _, _, z2 = binary_orbit_abs(bparams, tepoch)
  zord1 = z1 > z2 ? 2 : 3
  zord2 = z1 > z2 ? 3 : 2
  # Star 1 at origin; star 2 shifted by the orbital offset
  add_tessel_collection!(ax, star1, colours1; plotmesh=plotmesh, zorder=zord1)
  add_tessel_collection!(ax, star2, colours2; plotmesh=plotmesh, zorder=zord2,
                         offset_west=offset_west, offset_north=offset_north)
  ax.set_xlabel(L"x $\leftarrow$ E (mas)", fontsize=20)
  ax.set_ylabel(L"y $\rightarrow$ N (mas)", fontsize=20)
  # Decorations: limb (z=4) < graticules (z=5) < pole line (z=6) < spin arrow (z=7) < compass (z=8)
  #
  # The limb matters more here than for a single star. The two components sit on ONE shared
  # colour scale, so the hotter one lands at the pale end of the colormap and — on a white
  # background — becomes invisible, taking its rotation axis and spin arrow with it into
  # apparently empty space. Drawn here rather than inside the graticule block so it appears
  # whether or not graticules were asked for.
  # Each component's decorations sit just above ITS OWN surface, NOT above both. Drawing
  # every limb at one fixed zorder paints the FAR star's outline over the NEAR star, which
  # reads as the near star being transparent. Fractional zorders interleave them:
  #   far surface 2 < far limb 2.4 < near surface 3 < near limb 3.4
  #
  # The limb is drawn HERE in both cases, and `draw_graticules` is told not to draw its
  # own (`limb = false`). Its internal one sits at `zorder - 1`, which is right for the
  # single-star default of 5 but lands at `zord - 0.5` here — i.e. UNDER that component's
  # own surface, where it is invisible. That is what silently removed the outline from
  # every binary plotted with `graticules = true`. Drawing it here also gets it
  # `limb_color`/`limb_linewidth` rather than the graticule's colour and alpha.
  if limb
    draw_limb!(ax, star1; color=limb_color, linewidth=limb_linewidth, zorder=zord1+0.4)
    draw_limb!(ax, star2; offset_west=offset_west, offset_north=offset_north,
               color=limb_color, linewidth=limb_linewidth, zorder=zord2+0.4)
  end
  if graticules
    draw_graticules(ax, star1; inclination=inclination1, position_angle=position_angle1,
        star_params=star_params1, limb=false, zorder=zord1+0.5, graticule_kwargs...)
    draw_graticules(ax, star2; offset_west=offset_west, offset_north=offset_north,
        inclination=inclination2, position_angle=position_angle2, star_params=star_params2,
        limb=false, zorder=zord2+0.5, graticule_kwargs...)
  end
  if rotation_axis
    draw_rotation_axis(ax, star1, inclination=inclination1, position_angle=position_angle1, star_params=star_params1)
    draw_rotation_axis(ax, star2, offset_west=offset_west, offset_north=offset_north,
        inclination=inclination2, position_angle=position_angle2, star_params=star_params2)
  end
  if rotation_arrow
    draw_rotation_arrow(ax, star1, inclination=inclination1, position_angle=position_angle1, star_params=star_params1)
    draw_rotation_arrow(ax, star2, offset_west=offset_west, offset_north=offset_north,
        inclination=inclination2, position_angle=position_angle2, star_params=star_params2)
  end
  if compass; draw_compass(ax, axis_max); end
  # Colorbar — use the padded norm so colors match the patches
  if colorbar_on
    cmap = ColorMap(colormap)
    divider = axdiv.make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.07)
    cb = colorbar(matplotlib.cm.ScalarMappable(norm=norm_plot, cmap=cmap), cax=cax)
    if tmin < tmax; cb.ax.set_ylim(tmin, tmax); end  # clip colorbar ticks to actual temperature range
  end
  return fig, ax
end

"""
    plot2d_wireframe(star; compass=true, rotation_axis=false, hidden=true,
                     front_color="black", hidden_color="lightgrey", linewidth=0.5)

Wireframe view of the tessellation projected onto the sky plane.

With `hidden=true` the far-side tessels are drawn too, in `hidden_color`, so the mesh reads
as a transparent globe and the near/far hemispheres are told apart by edge weight rather
than by one hiding the other. This needs the near side UNFILLED: an opaque white fill looks
identical on a white background but occludes everything behind it, leaving no back edges to
see.

`hidden=false` restores the opaque near-side-only view.
"""
function plot2d_wireframe(star; compass=true, rotation_axis=false, hidden=true,
                          front_color="black", hidden_color="lightgrey", linewidth=0.5)
  # Global rcParams: without this the font depends on what was plotted before.
  set_oiplot_defaults()
  # Wireframe view of the tessellation projected onto the sky plane
  fig = figure("Epoch tmap",figsize=(10,10),facecolor="White")
  ax = fig.add_axes([0.1,0.1,0.85,0.85]);
  xlabel(L"x $\leftarrow$ E (mas)", fontsize=14);
  ylabel(L"y $\rightarrow$ N (mas)", fontsize=14);
  axis_max = maximum(sqrt.(star.vertices_xyz[:,:,1].^2 .+ star.vertices_xyz[:,:,2].^2 .+ star.vertices_xyz[:,:,3].^2))*1.5;
  ax.set_xlim([axis_max,-axis_max]);
  ax.set_ylim([-axis_max,axis_max]);
  # One PolyCollection instead of nquads_visible separate Polygon patches: matplotlib
  # re-validates and re-transforms every artist it holds, so the per-patch loop cost grew
  # with the mesh (noticeable from nside 4 upward, and it is the whole plot at nside 6).
  if hidden
    # Far side first and lighter, near side unfilled on top of it.
    back = setdiff(1:star.npix, star.index_quads_visible)
    isempty(back) || add_tessel_collection!(ax, star, "none"; indices=back,
        edgecolors=hidden_color, linewidths=linewidth, zorder=1)
    add_tessel_collection!(ax, star, "none"; edgecolors=front_color,
        linewidths=linewidth, zorder=2)
  else
    add_tessel_collection!(ax, star, "white"; edgecolors=front_color,
        linewidths=linewidth, zorder=2)
  end
  ax.tick_params(axis="both", which="both", labelsize=15, width=1, length=5);
  set_tick_spacing(ax, axis_max)
  ax.tick_params(axis="both", which="minor", width=1, length=5);
  for spine in ["top", "bottom", "left", "right"]
    ax.spines[spine].set_linewidth(1);
  end
  ax.set_aspect("equal")
  if compass; draw_compass(ax, axis_max, color="black"); end
  if rotation_axis; draw_rotation_axis(ax, star); end
  ax.plot();
  pyplot.draw();
  return;
end

"""
    plot2d_allepochs(tmap, star; plotmesh=false, tepochs=[], colormap="gist_heat",
                     arr_box=nothing, ncols=nothing, compass=true)

Grid of sky-plane maps, one panel per epoch, on a shared colour scale.

The grid and figure size follow the number of epochs, so three epochs give a 1x3 row
rather than a half-empty 2x3 block with the axis labels floating below it.

`ncols` pins the column count; `arr_box` is also accepted, in the legacy two-digit
`<rows><cols>` form, for existing scripts.
"""
function plot2d_allepochs(tmap, star; plotmesh=false, tepochs = [], colormap="gist_heat",
                          arr_box=nothing, ncols=nothing, compass=true)
    set_oiplot_defaults()
    nepochs = length(star)
    # Layout: explicit ncols, else the legacy arr_box, else auto (max 4 across).
    if ncols !== nothing
        nc = ncols
        nr = cld(nepochs, nc)
    elseif arr_box !== nothing
        nr, nc = divrem(arr_box, 10)
    else
        nc = min(nepochs, 4)
        nr = cld(nepochs, nc)
    end
    fig = figure("Temperature map -- All epochs",
                 figsize=(4.2*nc + 1.2, 4.6*nr + 0.8), facecolor="White")
    if plotmesh == true
      meshcolor = "grey"
    else
      meshcolor = "none"
    end
    minT = minimum(tmap);
    maxT = maximum(tmap);

    # Handle uniform maps
    if minT == maxT
      minT = 0.95*minT
      maxT = 1.05*maxT
    end

    for t=1:nepochs
      ax = fig.add_subplot(nr, nc, t)
      if tepochs !=[]
        title("Epoch $t — $(tepochs[t])")
      end
      axis_max = maximum(sqrt.(star[t].vertices_xyz[:,:,1].^2 .+ star[t].vertices_xyz[:,:,2].^2 .+ star[t].vertices_xyz[:,:,3].^2))*1.5;
      ax.set_xlim([axis_max,-axis_max])
      ax.set_ylim([-axis_max,axis_max])
      projmap = (tmap[star[t].index_quads_visible].-minT)./(maxT-minT);
      # One PolyCollection per epoch rather than nepochs x nquads_visible patches — this
      # was the worst offender of the three, being a product of two loops.
      cols = rgba(get_cmap(colormap), view(projmap, 1:star[t].nquads_visible))
      add_tessel_collection!(ax, star[t], cols; edgecolors=meshcolor,
                             linewidths = plotmesh ? 0.2 : 0.0, zorder=2)
      ax.plot();
      # Tick styling comes from set_oiplot_defaults (i.e. OITOOLS' rcParams). Hardcoding
      # labelsize/width/length here would desynchronise these panels from plot2d's ticks.
      set_tick_spacing(ax, axis_max)
      for spine in ["top", "bottom", "left", "right"]
        ax.spines[spine].set_linewidth(1);
      end
      ax.set_aspect("equal");
      if compass && t == 1; draw_compass(ax, axis_max, color="black"); end
    end
    # supxlabel/supylabel place the shared labels relative to the ACTUAL axes, so they
    # follow the grid instead of floating at a fixed figure fraction below an empty row.
    fig.supxlabel(L"x $\leftarrow$ E (mas)")
    fig.supylabel(L"y $\rightarrow$ N (mas)")
    fig.tight_layout()
    fig.canvas.draw();
    return;
end

"""
    plot_rv(bparams; K1, K2, γ=0.0, rv_data1=nothing, rv_data2=nothing, figtitle="Radial Velocities")

Plot radial velocity model curves vs orbital phase, optionally overlaying data.
`rv_data1`/`rv_data2` should be Nx3 matrices with columns [JD, RV_km/s, error_km/s].
"""
function plot_rv(bparams; K1::Float64, K2::Float64, γ::Float64=0.0,
    rv_data1=nothing, rv_data2=nothing, figtitle="Radial Velocities")
  set_oiplot_defaults()
  phases_plot = collect(range(0, 1, length=500))
  tepochs_rv = bparams.T0 .+ phases_plot .* bparams.P
  rv1_model, rv2_model = binary_RV(bparams, tepochs_rv, K1=K1, K2=K2, γ=γ)

  fig, ax = subplots(1, 1, figsize=(10, 6))
  ax.plot(phases_plot, rv1_model, "b-", linewidth=2, label="Primary model")
  ax.plot(phases_plot, rv2_model, "r-", linewidth=2, label="Secondary model")

  if rv_data1 !== nothing
    phi1 = mod.(rv_data1[:,1] .- bparams.T0, bparams.P) ./ bparams.P
    if size(rv_data1, 2) >= 3
      ax.errorbar(phi1, rv_data1[:,2], yerr=rv_data1[:,3], fmt="o", color="blue", ms=4, label="Primary data", zorder=5)
    else
      ax.scatter(phi1, rv_data1[:,2], color="blue", marker="o", s=30, label="Primary data", zorder=5)
    end
  end
  if rv_data2 !== nothing
    phi2 = mod.(rv_data2[:,1] .- bparams.T0, bparams.P) ./ bparams.P
    if size(rv_data2, 2) >= 3
      ax.errorbar(phi2, rv_data2[:,2], yerr=rv_data2[:,3], fmt="s", color="red", ms=4, label="Secondary data", zorder=5)
    else
      ax.scatter(phi2, rv_data2[:,2], color="red", marker="s", s=30, label="Secondary data", zorder=5)
    end
  end

  ax.axhline(γ, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)
  ax.set_xlabel("Orbital Phase")
  ax.set_ylabel("Radial Velocity (km/s)")
  ax.set_xlim(0, 1)
  ax.legend()
  ax.set_title(figtitle)
  ax.grid(true, alpha=0.3)
  tight_layout()
  return fig, ax
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

function plot_mollweide(tmap, star; ntheta = nothing, nphi = nothing, kwargs...)
  if star.tessellation_type == 0
    mollplot_temperature_healpix(tmap; kwargs...)
  else
    # `mollplot_temperature_longlat` needs the grid dimensions, which `stellar_geometry`
    # does not store. Recover them from the mesh unless the caller pins them explicitly.
    nt, np = (ntheta === nothing || nphi === nothing) ? _latlong_dims(star) : (ntheta, nphi)
    mollplot_temperature_longlat(tmap, nt, np; kwargs...)
  end
  return
end

function mollplot_temperature_healpix(tmap; visible_pixels = [], vmin = -Inf, vmax = Inf, incl=90.0, colormap="gist_heat", figtitle="Mollweide", mask_unobserved=true, bad_color="lightgray", lon_color="white", lat_color="black")
  set_oiplot_defaults()
  np = pyimport("numpy")
  xsize = 2000
  ysize = div(xsize,2)
  theta = collect(range(pi, stop=0.0, length=ysize))
  phi   = collect(range(-pi, stop=pi, length=xsize))
  longitude = collect(range(-179.999, stop=179.999, length=xsize))/180*pi
  latitude = collect(range(-89.999, stop=89.999, length=ysize))/180*pi
  # project the map to a rectangular matrix xsize x ysize
  nside = npix2nside(length(tmap))
  PHI = [i for j in theta, i in phi]
  THETA = [j for j in theta, i in phi]
  grid_pix = reshape(ang2pix_nest(nside, vec(THETA), vec(PHI)), size(PHI))
  grid_map = Float64.(tmap[grid_pix])
  # Mask unobserved pixels
  if mask_unobserved && visible_pixels != []
    vis_set = Set(visible_pixels)
    for idx in eachindex(grid_pix)
      if !(grid_pix[idx] in vis_set)
        grid_map[idx] = NaN
      end
    end
  end
  fig = figure(figtitle, figsize=(10, 7))
  fig.clear();
  ax = subplot(111,projection="mollweide")
  # rasterized makes the map bitmap while the labels remain vectorial
  # NOTE: Matplotlib Mollweide has longitude increasing to the right (math convention).
  # Astronomical convention has RA increasing to the left. For stellar surface maps
  # the current orientation shows the star as if unwrapped; negate `longitude` if the
  # observer-facing (sub-observer) convention is preferred. Verify with an asymmetric map.
  if visible_pixels == []
    if (vmin == -Inf)
      vmin = minimum(tmap);
    end
    if (vmax == Inf)
      vmax = maximum(tmap);
    end
  else
    if (vmin == -Inf)
      vmin = minimum(tmap[visible_pixels]);
    end
    if (vmax == Inf)
      vmax = maximum(tmap[visible_pixels]);
    end
  end
  cmap_obj = matplotlib.colormaps[colormap]   # matplotlib.cm.get_cmap removed in 3.9
  cmap_obj.set_bad(color=bad_color)
  moll = pcolormesh(longitude, latitude, grid_map, vmin=vmin, vmax=vmax, rasterized=true, cmap=cmap_obj)
  # graticule
  ax.set_longitude_grid(20);
  ax.set_latitude_grid(20);
  ax.set_longitude_grid_ends(90);
  spacing = 0.04
  subplots_adjust(bottom=spacing, top=1-spacing, left=spacing, right=1-spacing);
  grid(true)
  if incl != 90.0
    ax.axhline(-incl * pi/180, c="black", ls="-.")
  end
  ticks = collect(range(vmin, stop=vmax, length=7));
  cb = colorbar(moll, orientation="horizontal", shrink=.6, pad=0.05, ticks=ticks)
  cb.ax.xaxis.labelpad=5
  cb.ax.xaxis.set_label_text("Temperature (K)")
  # workaround for issue with viewers, see colorbar docstring
  cb.solids.set_edgecolor("face")
  # `lon_color`/`lat_color` colour BOTH the tick labels and the corresponding graticule
  # lines: meridians are the x-axis grid, parallels the y-axis grid. Setting only
  # tick_params (as this did) left the graticule at matplotlib's default grey, so the
  # keywords appeared to do nothing to the grid they are named after.
  ax.tick_params(axis="x", labelsize=12, colors=lon_color)
  ax.tick_params(axis="y", labelsize=12, colors=lat_color)
  ax.xaxis.grid(true, color=lon_color, alpha=0.6)
  ax.yaxis.grid(true, color=lat_color, alpha=0.6)
  return
end

function mollplot_temperature_longlat(tmap, ntheta, nphi; visible_pixels = [], vmin = -Inf, vmax = Inf, colormap="gist_heat", figtitle="Mollweide", incl=90.0)
  set_oiplot_defaults()
  xsize = 2000
  ysize = div(xsize,2)
  theta = collect(range(pi, stop=0, length=ysize))
  phi   = collect(range(-pi, stop=pi, length=xsize))
  longitude = collect(range(-180.0, stop=180.0, length=xsize))/180.0*pi
  latitude = collect(range(90.0, stop=-90.0, length=ysize))/180.0*pi
  # project the map to a rectangular matrix xsize x ysize
  PHI = [i for j in theta, i in phi]
  THETA = [j for j in theta, i in phi]
  grid_pix = longlat_ang2pix(ntheta, nphi, THETA, PHI);
  grid_map = tmap[circshift(grid_pix,(0,Int(xsize/2)))]
  fig = figure(figtitle, figsize=(10, 7))
  clf();
  ax = subplot(111,projection="mollweide",title=figtitle)
  if visible_pixels == []
    if (vmin == -Inf)
      vmin = minimum(tmap);
    end
    if (vmax == Inf)
      vmax = maximum(tmap);
    end
  else
    if (vmin == -Inf)
      vmin = minimum(tmap[visible_pixels]);
    end
    if (vmax == Inf)
      vmax = maximum(tmap[visible_pixels]);
    end
  end
  moll = pcolormesh(longitude, latitude, grid_map, vmin=vmin, vmax=vmax, rasterized=true, cmap=colormap)
  # graticule
  ax.set_longitude_grid(30)
  ax.set_latitude_grid(30)
  ax.set_longitude_grid_ends(90)
  spacing = 0.04
  subplots_adjust(bottom=spacing, top=1-spacing, left=spacing, right=1-spacing)
  grid(true)
  if incl != 90.0
    ax.axhline(-incl * pi/180, c="black", ls="-.")
  end

  ticks = collect(range(vmin, stop=vmax, length=7));
  cb = colorbar(moll, orientation="horizontal", shrink=.6, pad=0.05, ticks=ticks)
  cb.ax.xaxis.labelpad=5
  cb.ax.xaxis.set_label_text("Temperature (K)")
  # workaround for issue with viewers, see colorbar docstring
  cb.solids.set_edgecolor("face")
  ax.tick_params(axis="x", labelsize=15)
  ax.tick_params(axis="y", labelsize=15)
  return
end
