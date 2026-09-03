# The live plot surfaces: draw once, then only ever push data.
#
# WHY THIS IS NOT JUST `plot2d_makie`.
#
# GLMakie allocates a plot's GPU buffers the moment the plot is inserted into a scene that is
# already on screen. Under QMLMakie the GL context belongs to Qt's render thread and is current
# only inside Qt's render callback, so an insertion from a QML callback — which is where every
# user action arrives — allocates with no context bound:
#
#     gl_renderobject: glGenBuffers returned invalid id. OpenGL Context active?
#
# It does not always bite. With nothing else competing for the context the redraw often gets a
# usable one; opening a ComboBox popup or a file dialog gives that surface its own context, and
# from then on every insertion fails. That is why "the first plot works and nothing after it
# does", and why a scripted test that never opens a dialog can pass on a machine where the GUI
# is plainly broken. (This is OITOOLS' measured finding, in src/gui/livecanvas.jl; ROTIR
# inherits the constraint unchanged because it is a property of QMLMakie, not of either
# package.)
#
# So: create every plot ONCE, before the window exists, and afterwards only assign to
# Observables. Updating an Observable uploads at the next render, on the render thread, with a
# valid context. Nothing here inserts a plot, a Legend or a Colorbar after `loadqml`.
#
# Two consequences worth stating, because both look like bugs otherwise:
#
#   * The legend, if there ever is one, must be HAND-DRAWN into a decoration-free Axis. A
#     `Makie.Legend` fixes its entry count at construction.
#   * A colorbar is hidden by collapsing its column (`colsize!(…, Fixed(0))`), not by a
#     `visible` attribute, which it does not have.
#
# The offline builders in src/oiplot_spheroid_makie.jl create plots freely — they build a fresh
# Figure that is not on screen, so the constraint does not apply — and they remain what the
# comparison against matplotlib checks. This file must produce the same geometry and the same
# colours as they do; two implementations that drift is a failure worth guarding against.

# Typography for the live canvases, and minor ticks.
#
# Larger than the offline defaults on purpose: those are set for a figure that ends up in a
# paper at a fixed size, and these are read on screen inside a panel that is a fraction of the
# window. The scale follows the screen (`live_plot_scale`), so a HiDPI panel gets the same
# apparent size rather than the same pixel count.
#
# FIVE minor intervals per major tick, which puts a subtick every 0.2 of a labelled step — the
# reading a sky map actually gets used for is "how far is that spot from the limb", and that is
# an interpolation between ticks.
const LIVE_MINOR_INTERVALS = 5

"Apply the live typography to an axis. Called on every redraw, since the scale can change."
function style_live_axis!(ax; scale = live_plot_scale())
    ax.xlabelsize[] = 17 * scale
    ax.ylabelsize[] = 17 * scale
    ax.xticklabelsize[] = 14 * scale
    ax.yticklabelsize[] = 14 * scale
    ax.titlesize[] = 18 * scale
    ax.xminorticksvisible[] = true
    ax.yminorticksvisible[] = true
    ax.xminorticks[] = Makie.IntervalsBetween(LIVE_MINOR_INTERVALS)
    ax.yminorticks[] = Makie.IntervalsBetween(LIVE_MINOR_INTERVALS)
    ax.xminorgridvisible[] = false
    ax.yminorgridvisible[] = false
    return ax
end

"""
    SkyCanvas

The sky-plane tessel view: one `poly!` whose vertices and colours are Observables.

`limb` and `graticule` are separate line plots rather than part of the polygon set, because
they are toggled independently and a `lines!` fed an empty vector draws nothing — which is how
a decoration is switched off without inserting or removing a plot.
"""
struct SkyCanvas
    figure::Any
    axis::Any
    polys::Makie.Observable{Vector{Vector{Makie.Point2f}}}
    colors::Makie.Observable{Vector{Makie.RGBAf}}
    strokecolors::Makie.Observable{Vector{Makie.RGBAf}}
    polyplot::Any
    limb::Makie.Observable{Vector{Makie.Point2f}}
    limbplot::Any
    grat::Makie.Observable{Vector{Makie.Point2f}}
    gratplot::Any
    compass::Makie.Observable{Vector{Makie.Point2f}}
    compassplot::Any
    compasstext::Makie.Observable{Vector{Makie.Point2f}}
    compasslabels::Makie.Observable{Vector{String}}
    compasstextplot::Any
    cbarlimits::Makie.Observable{Tuple{Float32,Float32}}
    cbarlabel::Makie.Observable{String}
    colormap::Makie.Observable{Any}
    colorbar::Any
    # The idle state. An empty canvas with a 0-1 colorbar beside it reads as a plot that failed
    # rather than as one with nothing to show, so both of these exist from the start: a text
    # placeholder (empty string draws nothing) and the colorbar's column, which is COLLAPSED
    # rather than hidden — a Colorbar has no `visible` attribute.
    message::Makie.Observable{String}
    messageplot::Any
    # What was last drawn. Kept so a colormap change can redo the value -> RGBA mapping (the
    # colours are computed on the Julia side, so the colormap Observable feeds only the
    # colorbar), and so the zoom clamp has the data extent to measure against.
    lastvalues::Base.RefValue{Vector{Float64}}
    laststroke::Base.RefValue{Bool}
    homespan::Base.RefValue{Float64}
    # The spin axis and the rotation arrow, as two more pre-built polylines. They are toggled
    # from the view bar, and a `lines!` fed an empty vector draws nothing — which is how a
    # decoration is switched off without inserting or removing a plot.
    axis3d::Makie.Observable{Vector{Makie.Point2f}}
    axisplot::Any
    spin::Makie.Observable{Vector{Makie.Point2f}}
    spinplot::Any
end

"""
    build_sky_canvas(fig) -> SkyCanvas

Create every plot the sky view will ever need. **Call before the window is shown.**
"""
function build_sky_canvas(fig)
    ax = Makie.Axis(fig[1, 1]; xlabel = "x ← E (mas)", ylabel = "y → N (mas)",
                    aspect = Makie.DataAspect())
    # East to the LEFT, the astronomical convention, and done with REVERSED LIMITS rather than
    # `xreversed = true` — `xlims!(ax, lo, hi)` with lo < hi sets `xreversed = false`, so the
    # flag was being cleared by the first draw and every sky view came out mirrored. The
    # offline layer (`style_sky_axis!`) has always reversed the limits; this now matches it.
    ax.xgridvisible = false
    ax.ygridvisible = false

    polys   = Makie.Observable(Vector{Makie.Point2f}[])
    colors  = Makie.Observable(Makie.RGBAf[])
    strokes = Makie.Observable(Makie.RGBAf[])
    # `color` and `strokecolor` are BOTH Vector{RGBAf}, never numbers. Handing Makie raw values
    # plus a colormap works offline but changes the colour Observable's element type the first
    # time a map with a different range arrives, which forces a plot rebuild — the one thing
    # this file exists to avoid. Stroking each polygon in its own face colour is also what
    # closes the antialiasing seams between tessels; see `add_tessel_collection_makie!`.
    polyplot = Makie.poly!(ax, polys; color = colors, strokecolor = strokes,
                           strokewidth = 0.35)

    limb  = Makie.Observable(Makie.Point2f[])
    limbp = Makie.lines!(ax, limb; color = :black, linewidth = 0.8)
    grat  = Makie.Observable(Makie.Point2f[])
    gratp = Makie.lines!(ax, grat; color = (:black, 0.55), linewidth = 0.6)
    comp  = Makie.Observable(Makie.Point2f[])
    compp = Makie.lines!(ax, comp; color = :black, linewidth = 1.4)
    ctpos = Makie.Observable(Makie.Point2f[])
    ctlab = Makie.Observable(String[])
    ctp   = Makie.text!(ax, ctpos; text = ctlab, fontsize = 12, font = :bold,
                        align = (:center, :center))

    lims  = Makie.Observable((0.0f0, 1.0f0))
    label = Makie.Observable("T (K)")
    cmap  = Makie.Observable{Any}(_padded_cmap("gist_heat"))
    cbar  = Makie.Colorbar(fig[1, 2]; colormap = cmap, colorrange = lims, label = label)
    axs   = Makie.Observable(Makie.Point2f[])
    axsp  = Makie.lines!(ax, axs; color = (:black, 0.8), linewidth = 1.4, linestyle = :dash)
    spn   = Makie.Observable(Makie.Point2f[])
    spnp  = Makie.lines!(ax, spn; color = :black, linewidth = 1.4)

    msg   = Makie.Observable("")
    msgp  = Makie.text!(ax, [Makie.Point2f(0, 0)]; text = msg, fontsize = 15,
                        color = Makie.RGBAf(0.45, 0.5, 0.55, 1), align = (:center, :center))

    c = SkyCanvas(fig, ax, polys, colors, strokes, polyplot, limb, limbp, grat, gratp,
                  comp, compp, ctpos, ctlab, ctp, lims, label, cmap, cbar, msg, msgp,
                  Ref(Float64[]), Ref(false), Ref(1.0), axs, axsp, spn, spnp)
    idle!(c, "no model — add one on the Model tab")
    return c
end

"""
    show_map!(canvas, star, values; ...)

Push one epoch's surface map into the sky view.

`values` is per VISIBLE tessel, in the order `tessel_polygons` produces — the same contract
`add_tessel_collection_makie!` has, and getting it wrong scrambles the map with no error.
"""
function show_map!(c::SkyCanvas, star, values;
                   colorrange = nothing, offset_west = 0.0, offset_north = 0.0,
                   limb = true, graticules = false, compass = true, plotmesh = false,
                   rotation_axis = false, rotation_arrow = false,
                   graticule_nlat::Int = 5, graticule_nlon::Int = 8,
                   graticule_color::Symbol = :black,
                   star_params = nothing, pad = 0.5, title = "")
    busy!(c)
    cr = colorrange === nothing ? _map_range(values) : colorrange
    c.lastvalues[] = collect(Float64, values)
    c.laststroke[] = plotmesh
    cols = map_colors(values, c.colormap[], cr)
    c.polys[]  = tessel_polygons(star; offset_west, offset_north)
    c.colors[] = cols
    c.strokecolors[] = plotmesh ?
        fill(Makie.RGBAf(0.45, 0.45, 0.45, 1), length(cols)) : cols
    c.cbarlimits[] = (Float32(cr[1]), Float32(cr[2]))

    amax = sky_axis_max(star; pad)
    if limb
        vis = star.index_quads_visible
        hx, hy = convex_hull_2d(vec(-star.proj_west[vis, :] .- offset_west),
                                vec( star.proj_north[vis, :] .+ offset_north))
        c.limb[] = _closed_ring(hx, hy)
    else
        c.limb[] = Makie.Point2f[]
    end
    c.grat[] = graticules ?
        _flatten_segments(graticule_segments(star; star_params = star_params,
                                             nlat = graticule_nlat, nlon = graticule_nlon,
                                             offset_west = offset_west,
                                             offset_north = offset_north)) : Makie.Point2f[]
    c.gratplot.color[] = (graticule_color, 0.55)
    # The spin axis and the rotation arrow, from the same geometry the offline layer uses —
    # `_spin_axis` is in the shared core, so the two cannot draw different axes.
    c.axis3d[] = rotation_axis ?
        _axis_polyline(star, star_params, offset_west, offset_north) : Makie.Point2f[]
    c.spin[] = rotation_arrow ?
        _spin_polyline(star, star_params, offset_west, offset_north) : Makie.Point2f[]

    if compass
        cpts, cpos, clab = _compass_geometry(amax)
        c.compass[] = cpts; c.compasstext[] = cpos; c.compasslabels[] = clab
    else
        c.compass[] = Makie.Point2f[]; c.compasstext[] = Makie.Point2f[]
        c.compasslabels[] = String[]
    end
    style_live_axis!(c.axis)
    c.cbarlabel[] = c.cbarlabel[]              # keep the label; the size follows the theme
    c.colorbar.labelsize[] = 16 * live_plot_scale()
    c.colorbar.ticklabelsize[] = 13 * live_plot_scale()
    c.homespan[] = 2amax
    Makie.xlims!(c.axis, amax, -amax)          # reversed: East to the left
    Makie.ylims!(c.axis, -amax, amax)
    isempty(title) || (c.axis.title = title)
    return c
end

# ── zoom ─────────────────────────────────────────────────────────────────────────────────
#
# The sky view is a picture of a star of known angular size, so unlimited zoom is not a
# feature: zoomed far enough out the star is one pixel, far enough in it is one tessel, and
# both look like the plot broke. The bounds are relative to the FRAMED extent
# (`homespan`), so they mean the same thing for a 0.5 mas star and a 20 mas one.
const ZOOM_MIN_SPAN = 0.05      # 20x in
const ZOOM_MAX_SPAN = 4.0       # 4x out

"""
    reset_zoom!(canvas) -> canvas

Back to the framing `show_map!` chose. Bound to a right-click, as in the OITOOLS GUI.
"""
function reset_zoom!(c::SkyCanvas)
    a = c.homespan[] / 2
    Makie.xlims!(c.axis, a, -a); Makie.ylims!(c.axis, -a, a)   # reversed, as `show_map!` sets it
    return c
end

"""
    clamp_zoom!(canvas) -> Bool

Pull the view back inside [`ZOOM_MIN_SPAN`, `ZOOM_MAX_SPAN`] if scrolling has taken it out,
and report whether it had to. Called from a `finallimits` listener, so it also catches a
drag-zoom.

Rescales about the CURRENT CENTRE rather than the origin: a clamp that recentred would yank
the view sideways at the moment the user hit the bound, which reads as the plot fighting back.
"""
function clamp_zoom!(c::SkyCanvas)
    h = c.homespan[]
    h > 0 || return false
    fl = c.axis.finallimits[]
    wx, wy = Float64(fl.widths[1]), Float64(fl.widths[2])
    (isfinite(wx) && isfinite(wy) && wx > 0 && wy > 0) || return false
    lo, hi = ZOOM_MIN_SPAN * h, ZOOM_MAX_SPAN * h
    span = max(wx, wy)
    (lo <= span <= hi) && return false
    f = clamp(span, lo, hi) / span
    cx = Float64(fl.origin[1]) + wx / 2
    cy = Float64(fl.origin[2]) + wy / 2
    # HIGH then LOW on x, or the clamp silently un-reverses the axis and the sky flips
    # east-west the first time somebody scrolls past a bound.
    Makie.xlims!(c.axis, cx + wx * f / 2, cx - wx * f / 2)
    Makie.ylims!(c.axis, cy - wy * f / 2, cy + wy * f / 2)
    return true
end

# The spin axis, as ONE polyline with the tips included. The plot convention is
# `x = -West`, the same negation `tessel_polygons` applies.
function _axis_polyline(star, star_params, ow, on; arrow_frac = 0.3)
    north, south = _spin_axis(star, star_params, NaN, NaN)
    d = north .- south
    ntip = north .+ arrow_frac .* d
    stip = south .- arrow_frac .* d
    P(p) = Makie.Point2f(-(p[1] + ow), p[2] + on)
    return [P(stip), P(south), P(north), P(ntip)]
end

# The 300-degree spin circle about the north pole, split at z = 0 with a NaN break so the half
# behind the star lifts the pen — the same treatment `draw_rotation_arrow_makie!` gives it.
function _spin_polyline(star, star_params, ow, on; radius_frac = 0.07, offset_frac = 0.15,
                        npoints = 80)
    north, south = _spin_axis(star, star_params, NaN, NaN)
    ax = north .- south
    alen = LinearAlgebra.norm(ax)
    alen > 0 || return Makie.Point2f[]
    ahat = ax ./ alen
    centre = north .+ offset_frac .* ax
    r = radius_frac * alen
    ref = abs(ahat[3]) < 0.9 ? [0.0, 0.0, 1.0] : [1.0, 0.0, 0.0]
    e1 = LinearAlgebra.cross(ahat, ref); e1 ./= LinearAlgebra.norm(e1)
    e2 = LinearAlgebra.cross(ahat, e1)
    out = Makie.Point2f[]
    θs = collect(range(0, 300π/180, length = npoints))
    for θ in θs
        p = centre .+ r .* (cos(θ) .* e1 .+ sin(θ) .* e2)
        push!(out, p[3] > 0 ? Makie.Point2f(-(p[1] + ow), p[2] + on) :
                              Makie.Point2f(NaN, NaN))
    end
    # The arrowhead, as a "V" appended to the same polyline rather than a `scatter!` marker:
    # this canvas may not insert a plot, and the sense of rotation is the whole point of the
    # arrow — an arc without a head says nothing.
    #
    # Built on the TANGENT at the end of the arc, not on a chord back to an earlier sample:
    # over 300 degrees a chord points visibly the wrong way.
    θe  = θs[end]
    pe  = centre .+ r .* (cos(θe) .* e1 .+ sin(θe) .* e2)
    tan = -sin(θe) .* e1 .+ cos(θe) .* e2
    tan ./= LinearAlgebra.norm(tan)
    nrm = LinearAlgebra.cross(ahat, tan)
    h   = 0.45r
    P(q) = Makie.Point2f(-(q[1] + ow), q[2] + on)
    a1 = pe .- h .* tan .+ 0.5h .* nrm
    a2 = pe .- h .* tan .- 0.5h .* nrm
    append!(out, (Makie.Point2f(NaN, NaN), P(a1), P(pe), P(a2)))
    return out
end

# A `lines!` draws ONE polyline, so a set of disjoint curves is joined with NaN breaks rather
# than by inserting one plot per curve — inserting is exactly what is forbidden here. Makie
# treats a non-finite point as a pen lift.
function _flatten_segments(segs)
    out = Makie.Point2f[]
    for s in segs
        isempty(s) && continue
        # Each segment is an N x 2 MATRIX — column 1 x, column 2 y — which is what
        # matplotlib's LineCollection takes and therefore what `graticule_segments` returns.
        # Iterating it as a list of points walks it column-major instead.
        append!(out, (Makie.Point2f(s[i, 1], s[i, 2]) for i in axes(s, 1)))
        push!(out, Makie.Point2f(NaN, NaN))          # break, not a join to the next run
    end
    return out
end

# The hull comes back as two coordinate VECTORS (that is `convex_hull_2d`'s contract), and the
# ring is closed by repeating the first vertex — `lines!` draws an open polyline.
function _closed_ring(hx, hy)
    isempty(hx) && return Makie.Point2f[]
    pts = [Makie.Point2f(hx[i], hy[i]) for i in eachindex(hx)]
    push!(pts, Makie.Point2f(hx[1], hy[1]))
    return pts
end

# The compass as ONE polyline with a NaN break between its two arms, plus two label positions.
# Drawn from geometry rather than through `draw_compass_makie!` for the same reason as
# everything else in this file: that function inserts plots.
function _compass_geometry(amax; size_frac = 0.09)
    L  = size_frac * 2amax
    # NEGATIVE x, because the axis is reversed: +x is East and East is drawn to the LEFT, so
    # the bottom-RIGHT corner of the screen is at negative x. Placing it at +0.72 put the
    # compass in the bottom-left, where the offline layer puts nothing.
    ox, oy = -0.72amax, -0.80amax
    pts = [Makie.Point2f(ox, oy), Makie.Point2f(ox, oy + L), Makie.Point2f(NaN, NaN),
           Makie.Point2f(ox, oy), Makie.Point2f(ox + L, oy)]
    pos = [Makie.Point2f(ox, oy + 1.25L), Makie.Point2f(ox + 1.3L, oy)]
    return pts, pos, ["N", "E"]
end

"""
    StarCanvas

The rotatable 3-D view: one `mesh!` in an `Axis3`, vertices and colours as Observables.

An `Axis3`, not an `LScene`: this view is read against a SCALE — how big the star is, how far
apart the two components of a binary are — and an LScene has no ticks to read it off. Axis3 is
just as rotatable (left-drag orbits) and `perspectiveness = 0` makes it orthographic, which is
what the observer at infinity requires and what keeps it agreeing with the sky view beside it.
"""
struct StarCanvas
    figure::Any
    scene::Any
    mesh::Makie.Observable{Any}
    colors::Makie.Observable{Vector{Makie.RGBAf}}
    meshplot::Any
    mesh2::Makie.Observable{Any}
    colors2::Makie.Observable{Vector{Makie.RGBAf}}
    meshplot2::Any
    orbit::Makie.Observable{Vector{Makie.Point3f}}
    orbitplot::Any
    colormap::Makie.Observable{Any}
    cbarlimits::Makie.Observable{Tuple{Float32,Float32}}
    colorbar::Any
    message::Makie.Observable{String}
    messageplot::Any
end

"""
    build_star_canvas(fig) -> StarCanvas

Two meshes, always: the second is the binary's secondary and stays empty for a single star.
Creating it up front rather than when a binary is first shown is the whole point — inserting
it later is the failure this file is organised around. An empty mesh draws nothing.
"""
function build_star_canvas(fig)
    # Labels in mas, on all three axes. The mesh frame is (West, North, toward-observer), which
    # is right-handed — see `star_mesh` — so viewed down +z it puts East on the left, matching
    # the sky view.
    #
    # The opening camera is the observer's, tilted 10 degrees off it.
    #
    # Straight down the line of sight is the view the orthographic panel already gives, so a
    # 3-D scene that opens there looks like a flat picture with axes round it and gives no
    # reason to drag. FORTY-FIVE degrees: the z axis gets the same visible extent as x and y,
    # the body reads as a solid rather than a disc, and the near hemisphere — the one the data
    # constrain — still fills most of the frame.
    ls = Makie.Axis3(fig[1, 1]; aspect = :data, perspectiveness = 0f0,
                     xlabel = "x → W (mas)", ylabel = "y → N (mas)", zlabel = "z → obs (mas)",
                     azimuth = pi/2, elevation = pi/2 - deg2rad(45))

    # The colour vector must match the mesh's vertex count, EMPTY INCLUDED. A 3-vertex
    # placeholder mesh with a zero-length colour vector fails to build a render object at all
    # ("Failed to resolve gl_renderobject", with `scaled_color = RGBAf[]`), and under QMLMakie
    # that surfaces only as "exception in render" — the plot never appears and nothing says
    # why. Fully transparent, so the placeholder is invisible even where it is drawn.
    m1 = Makie.Observable{Any}(_empty_mesh())
    c1 = Makie.Observable(_empty_colors())
    p1 = Makie.mesh!(ls, m1; color = c1, shading = Makie.NoShading)
    m2 = Makie.Observable{Any}(_empty_mesh())
    c2 = Makie.Observable(_empty_colors())
    p2 = Makie.mesh!(ls, m2; color = c2, shading = Makie.NoShading)
    orb = Makie.Observable(Makie.Point3f[])
    op  = Makie.lines!(ls, orb; color = (:grey35, 0.8), linewidth = 1.2)

    cmap = Makie.Observable{Any}(_padded_cmap("gist_heat"))
    lims = Makie.Observable((0.0f0, 1.0f0))
    cbar = Makie.Colorbar(fig[1, 2]; colormap = cmap, colorrange = lims, label = "T (K)")
    # In a Label rather than in the scene: an LScene has no data-space anchor that stays put as
    # the camera turns, so a `text!` inside it would drift off-screen on the first drag.
    msg  = Makie.Observable("")
    # A Label under the axis rather than text inside it: an Axis3's data space rotates with the
    # camera, so a `text!` anchored in it would drift off-screen on the first drag.
    msgp = Makie.Label(fig[2, 1], msg; fontsize = 15, color = Makie.RGBAf(0.45, 0.5, 0.55, 1),
                       tellwidth = false)
    c = StarCanvas(fig, ls, m1, c1, p1, m2, c2, p2, orb, op, cmap, lims, cbar, msg, msgp)
    idle!(c, "no model — add one on the Model tab")
    return c
end

# A degenerate single triangle rather than a zero-vertex mesh: GeometryBasics rejects the
# latter, and Makie needs SOMETHING to allocate buffers for at build time. All three vertices
# coincide at the origin, so it covers no pixels.
#
# TWO of them, 3-D and 2-D, and the distinction is load-bearing. Makie fixes a mesh plot's
# vertex type at construction and CONVERTS every later assignment to it — a 2-D mesh handed to
# a plot seeded with a 3-D one raises `Cannot convert Mesh{2,...} to Mesh{3,...}` inside the
# compute graph, where QMLMakie's render callback swallows it and prints only
# "exception in render". The Mollweide is a flat projection and must be seeded 2-D.
_empty_mesh() = Makie.GeometryBasics.Mesh(
    [Makie.Point3f(0, 0, 0), Makie.Point3f(0, 0, 0), Makie.Point3f(0, 0, 0)],
    [Makie.TriangleFace(1, 2, 3)])
_empty_mesh2d() = Makie.GeometryBasics.Mesh(
    [Makie.Point2f(0, 0), Makie.Point2f(0, 0), Makie.Point2f(0, 0)],
    [Makie.TriangleFace(1, 2, 3)])
_empty_colors() = fill(Makie.RGBAf(0, 0, 0, 0), 3)

"""
    show_star3d!(canvas, star, values; colorrange = nothing)

Push a single star's full surface into the 3-D view. The secondary mesh is emptied, so
switching from a binary back to one component leaves nothing behind.
"""
function show_star3d!(c::StarCanvas, star, values; colorrange = nothing)
    busy!(c)
    cr = colorrange === nothing ? _map_range(values) : colorrange
    msh, vals = star_mesh(star; values = values)
    c.mesh[]   = msh
    c.colors[] = map_colors(vals, c.colormap[], cr)
    c.mesh2[]  = _empty_mesh()
    c.colors2[] = _empty_colors()
    c.orbit[]  = Makie.Point3f[]
    c.cbarlimits[] = (Float32(cr[1]), Float32(cr[2]))
    _recenter_cam!(c, _body_radius(star))
    return c
end

"""
    show_binary3d!(canvas, star1, v1, star2, v2, bparams, tepoch; colorrange = nothing)

Both components at their full 3-D separation, plus the relative orbit track.
"""
function show_binary3d!(c::StarCanvas, star1, v1, star2, v2, bparams, tepoch;
                        colorrange = nothing)
    busy!(c)
    lo = colorrange === nothing ? min(minimum(v1), minimum(v2)) : colorrange[1]
    hi = colorrange === nothing ? max(maximum(v1), maximum(v2)) : colorrange[2]
    hi - lo < 1 && (hi = lo + max(abs(hi) * 0.01, 1.0))
    x1, y1, z1, x2, y2, z2 = binary_orbit_abs(bparams, tepoch)
    off = (-(y2 - y1), x2 - x1, -(z2 - z1))       # (West, North, toward observer)
    m1, c1 = star_mesh(star1; values = v1)
    m2, c2 = star_mesh(star2; offset = off, values = v2)
    c.mesh[]    = m1; c.colors[]  = map_colors(c1, c.colormap[], (lo, hi))
    c.mesh2[]   = m2; c.colors2[] = map_colors(c2, c.colormap[], (lo, hi))
    c.orbit[]   = relative_orbit_track(bparams)
    c.cbarlimits[] = (Float32(lo), Float32(hi))
    sep = sqrt(off[1]^2 + off[2]^2 + off[3]^2)
    _recenter_cam!(c, max(sep, _body_radius(star1) + _body_radius(star2)))
    return c
end

# Frame the body. Called on every new map: a star three times the size of the last one would
# otherwise open off-screen, and the user cannot tell that from a failed draw.
#
# `limits!` rather than a camera move, because an Axis3 derives its camera FROM its limits —
# and setting them is also what puts sensible numbers on the mas ticks. The camera ANGLE is
# left alone on purpose: it is the one piece of view state the user has set by dragging, and
# resetting it on every parameter edit would fight them.
function _recenter_cam!(c::StarCanvas, radius)
    r = isfinite(radius) && radius > 0 ? Float64(radius) * 1.05 : 1.0
    Makie.limits!(c.scene, (-r, r), (-r, r), (-r, r))
    return c
end

"""
    Chi2Canvas

Per-epoch χ², as a grouped bar chart: V², T3 amplitude and T3 phase side by side.

Reduced rather than absolute, because the point of the panel is to see WHICH epoch the model
fits badly, and epochs rarely have the same number of points.
"""
struct Chi2Canvas
    figure::Any
    axis::Any
    message::Makie.Observable{String}
    messageplot::Any
    barplot::Any
    unity::Any
end

function build_chi2_canvas(fig)
    ax = Makie.Axis(fig[1, 1]; xlabel = "epoch", ylabel = "χ²ᵣ")
    # Seeded with one zero-height bar rather than empty. `barplot!` with `dodge` computes
    # `maximum(dodge)` to work out how many groups there are, and `maximum` of an empty vector
    # throws — so an empty chart cannot be BUILT at all, which is fatal here: the plot has to
    # exist before the window does.
    xs, ys, gp = _chi2_placeholder()
    bp = Makie.barplot!(ax, xs, ys; dodge = gp, color = gp,
                        colormap = [Makie.RGBAf(0.80, 0.16, 0.16, 1),
                                    Makie.RGBAf(0.16, 0.35, 0.80, 1),
                                    Makie.RGBAf(0.13, 0.55, 0.24, 1)],
                        colorrange = (1, 3))
    # χ²ᵣ = 1 is the line the eye actually looks for, so it is drawn once and never touched.
    un = Makie.hlines!(ax, [1.0]; color = (:black, 0.45), linestyle = :dash, linewidth = 1)
    msg  = Makie.Observable("")
    msgp = Makie.text!(ax, [Makie.Point2f(1, 0.5)]; text = msg, fontsize = 14,
                       color = Makie.RGBAf(0.45, 0.5, 0.55, 1), align = (:center, :center))
    c = Chi2Canvas(fig, ax, msg, msgp, bp, un)
    idle!(c, "no dataset and no model")
    return c
end

"""
    show_chi2!(canvas, breakdowns)

Push one [`chi2_breakdown`](@ref) per epoch. Non-finite reduced values — an epoch with no
closure phases, say — are dropped rather than plotted as zero, which would read as a perfect fit
to data that is not there.

`Makie.update!` in ONE call, not three Observable assignments. `barplot!` resolves its compute
graph after each input changes, so setting `x` and then `y` leaves an instant where they have
different lengths, and the recipe raises `DimensionMismatch` from inside the graph — where
QMLMakie swallows it and prints only "exception in render". The batched update makes the three
arrive as a single transaction.
"""
function show_chi2!(c::Chi2Canvas, breakdowns)
    style_live_axis!(c.axis)
    xs = Float32[]; ys = Float32[]; gp = Int[]
    for (i, b) in enumerate(breakdowns), (g, v) in enumerate((b.v2r, b.t3ampr, b.t3phir))
        isfinite(v) || continue
        push!(xs, Float32(i)); push!(ys, Float32(v)); push!(gp, g)
    end
    if isempty(ys)
        # Same reason as at build time: an empty `dodge` makes the recipe throw on the next
        # recompute. A single zero-height bar reads as "nothing to show" and costs one pixel.
        px, py, pg = _chi2_placeholder()
        Makie.update!(c.barplot, px[], py[]; dodge = pg[], color = pg[])
        idle!(c, "no χ² yet — load a dataset and set a model")
        return c
    end
    busy!(c)
    Makie.update!(c.barplot, xs, ys; dodge = gp, color = gp)
    Makie.ylims!(c.axis, 0, max(1.2, 1.1maximum(ys)))
    Makie.xlims!(c.axis, 0.4, length(breakdowns) + 0.6)
    return c
end

_chi2_placeholder() = (Makie.Observable(Float32[1]), Makie.Observable(Float32[0]),
                       Makie.Observable(Int[1]))

"""
    MollCanvas

The whole-surface Mollweide view: one `mesh!` over the projected grid.

The graticule and the bounding ellipse are static — they depend only on the projection, never
on the data — so they are drawn once at build time and never updated. That is not merely an
optimisation here: it is one fewer thing that could need a plot inserted later.
"""
struct MollCanvas
    figure::Any
    axis::Any
    mesh::Makie.Observable{Any}
    colors::Makie.Observable{Vector{Makie.RGBAf}}
    meshplot::Any
    colormap::Makie.Observable{Any}
    cbarlimits::Makie.Observable{Tuple{Float32,Float32}}
    colorbar::Any
    message::Makie.Observable{String}
    messageplot::Any
end

function build_moll_canvas(fig)
    ax = Makie.Axis(fig[1, 1]; aspect = Makie.DataAspect())
    Makie.hidedecorations!(ax); Makie.hidespines!(ax)
    m = Makie.Observable{Any}(_empty_mesh2d())
    c = Makie.Observable(_empty_colors())
    mp = Makie.mesh!(ax, m; color = c, shading = Makie.NoShading)
    # Static graticule: meridians and parallels every 20 degrees, one polyline each way with
    # NaN breaks, traced through the same projection as the data so they cannot drift from it.
    lon = Makie.Point2f[]
    for λ0 in range(-π, π, length = 19)
        append!(lon, (Makie.Point2f(mollweide_xy(λ0, φ)...)
                      for φ in range(-π/2, π/2, length = 120)))
        push!(lon, Makie.Point2f(NaN, NaN))
    end
    Makie.lines!(ax, lon; color = Makie.RGBAf(1, 1, 1, 0.75), linewidth = 0.6)
    lat = Makie.Point2f[]
    for φ0 in range(-π/2, π/2, length = 10)
        append!(lat, (Makie.Point2f(mollweide_xy(λ, φ0)...)
                      for λ in range(-π, π, length = 240)))
        push!(lat, Makie.Point2f(NaN, NaN))
    end
    Makie.lines!(ax, lat; color = Makie.RGBAf(0, 0, 0, 0.5), linewidth = 0.6)
    # The bounding ellipse, and LABELS on both graticules — the matplotlib Mollweide has them
    # and a projection without them is a shape rather than a map. All three are traced through
    # the same projection as the data, so they cannot drift from it.
    Makie.lines!(ax, vcat([Makie.Point2f(mollweide_xy(-π + 1e-9, φ)...)
                           for φ in range(-π/2, π/2, length = 180)],
                          [Makie.Point2f(mollweide_xy(π - 1e-9, φ)...)
                           for φ in range(π/2, -π/2, length = 180)]);
                 color = (:black, 0.6), linewidth = 1.0)
    lonlab = [Makie.Point2f(mollweide_xy(deg2rad(λ), 0.0)...) for λ in -160:20:160]
    Makie.text!(ax, lonlab; text = ["$(λ)°" for λ in -160:20:160],
                color = Makie.RGBAf(1, 1, 1, 0.9), fontsize = 11,
                align = (:center, :center))
    # Latitudes just outside the ellipse at their own latitude, not at a fixed longitude: the
    # ellipse narrows towards the poles, so a fixed longitude walks the high labels into the map.
    latlab = [Makie.Point2f(mollweide_xy(-π, deg2rad(φ))[1] - 0.06,
                            mollweide_xy(-π, deg2rad(φ))[2]) for φ in (-80,-60,-40,-20,20,40,60,80)]
    Makie.text!(ax, latlab; text = ["$(φ)°" for φ in (-80,-60,-40,-20,20,40,60,80)],
                color = :black, fontsize = 11, align = (:right, :center))

    cmap = Makie.Observable{Any}(_padded_cmap("gist_heat"))
    lims = Makie.Observable((0.0f0, 1.0f0))
    cbar = Makie.Colorbar(fig[2, 1]; colormap = cmap, colorrange = lims, vertical = false,
                          label = "T (K)", width = Makie.Relative(0.6))
    msg  = Makie.Observable("")
    msgp = Makie.text!(ax, [Makie.Point2f(0, 0)]; text = msg, fontsize = 15,
                       color = Makie.RGBAf(0.45, 0.5, 0.55, 1), align = (:center, :center))
    cc = MollCanvas(fig, ax, m, c, mp, cmap, lims, cbar, msg, msgp)
    idle!(cc, "no model — add one on the Model tab")
    return cc
end

"""
    show_mollweide!(canvas, tmap, star; ...)

Push a whole-surface map into the Mollweide view.

`visible_pixels` is empty by default and nothing is masked — the same default as
`plot_mollweide_makie`, and for the same reason: one epoch's `index_quads_visible` is the
front hemisphere at that phase, and masking with it would grey out half a map the other epochs
constrain perfectly well.
"""
function show_mollweide!(c::MollCanvas, tmap, star = nothing; visible_pixels = Int[],
                         colorrange = nothing, nlon = 288, nlat = 144, mask_unobserved = true,
                         bad_color = Makie.RGBAf(0.83, 0.83, 0.83, 1))
    busy!(c)
    X, Y, vals = mollweide_grid(tmap, star; nlon, nlat, visible_pixels, mask_unobserved)
    finite = filter(isfinite, vals)
    cr = colorrange === nothing ?
        _map_range(isempty(finite) ? Float64.(tmap) : finite) : colorrange
    ny, nx = size(vals)
    pts  = [Makie.Point2f(X[i, j], Y[i, j]) for j in 1:nx for i in 1:ny]
    cmap = c.colormap[]
    lo, hi = cr; s = hi > lo ? 1 / (hi - lo) : 0.0
    cols = Vector{Makie.RGBAf}(undef, length(pts))
    @inbounds for j in 1:nx, i in 1:ny
        v = vals[i, j]
        cols[(j - 1) * ny + i] = isfinite(v) ?
            Makie.RGBAf(Makie.to_color(cmap[clamp((v - lo) * s, 0, 1)])) : bad_color
    end
    faces = Makie.TriangleFace{Int}[]
    sizehint!(faces, 2 * (nx - 1) * (ny - 1))
    @inbounds for j in 1:nx-1, i in 1:ny-1
        a = (j - 1) * ny + i; b = a + 1; cc = j * ny + i; d = cc + 1
        push!(faces, Makie.TriangleFace(a, b, d)); push!(faces, Makie.TriangleFace(a, d, cc))
    end
    c.mesh[]   = Makie.GeometryBasics.Mesh(pts, faces)
    c.colors[] = cols
    c.cbarlimits[] = (Float32(cr[1]), Float32(cr[2]))
    return c
end

"""
    idle!(canvas, message)  /  busy!(canvas)

Put a canvas into, or take it out of, its EMPTY state.

An axis with nothing in it and a 0–1 colorbar beside it reads as a plot that failed rather than
as one with nothing to show — which is the single most common way an idle GUI looks broken. So
an empty canvas says what it is waiting for, and its colorbar column is COLLAPSED
(`colsize!(…, Fixed(0))`) rather than hidden: a Colorbar has no `visible` attribute, and
collapsing its column is the supported way to make it go away and come back.

Neither call inserts a plot. The message is a `text!` created at build time whose string is
empty most of the time, and a layout size is not a plot at all — so both are legal after the
window exists, which is the constraint everything in this file is arranged around.
"""
function idle!(c, message::AbstractString)
    c.message[] = String(message)
    _colorbar_visible!(c, false)
    return c
end

function busy!(c)
    isempty(c.message[]) || (c.message[] = "")
    _colorbar_visible!(c, true)
    return c
end

# BOTH mechanisms, because neither is enough alone. Collapsing the row or column
# (`Fixed(0)`) is what makes the layout give the space back to the plot, but the Colorbar
# still DRAWS: it renders into its own block scene and its bar is a protrusion, so a
# zero-width column leaves a full-height bar sitting beside the axis (measured — the
# reference in the OITOOLS notes is about reclaiming the space, not about hiding the block).
# `blockscene.visible` is what actually hides it. `Auto()` restores whatever the layout would
# have chosen on its own.
_colorbar_visible!(c::Chi2Canvas, ::Bool) = c            # no colorbar to collapse
function _colorbar_visible!(c, on::Bool)
    if c isa MollCanvas
        Makie.rowsize!(c.figure.layout, 2, on ? Makie.Auto() : Makie.Fixed(0))
    else
        Makie.colsize!(c.figure.layout, 2, on ? Makie.Auto() : Makie.Fixed(0))
    end
    try
        c.colorbar.blockscene.visible[] = on
    catch
        # A Makie version without `blockscene` still gets the layout collapse, which is the
        # half that matters for the plot's size.
    end
    return c
end



# ── colormaps ────────────────────────────────────────────────────────────────────────────
#
# The same list the OITOOLS GUI offers, minus the ones that make no sense for a temperature:
# `gist_heat` first because it is what every ROTIR figure in the papers uses, `viridis` for
# anyone who needs a perceptually uniform ramp, and the reversed pairs because a dark-on-light
# map prints better than a light-on-dark one.
const SURFACE_COLORMAPS = ["gist_heat"    => "gist_heat",
                           "gist_heat_r"  => "gist_heat_r",
                           "viridis"      => "viridis",
                           "inferno"      => "inferno",
                           "magma"        => "magma",
                           "gist_gray"    => "gist_gray",
                           "gist_gray_r"  => "gist_gray_r"]

colormap_names() = join(first.(SURFACE_COLORMAPS), "\n")

"""
    set_colormap!(canvas, name) -> Bool

Switch a canvas's colormap by its button label, and RE-COLOUR what is already drawn.

Reassigning `c.colormap` alone is not enough for the sky and Mollweide canvases: their colours
are `Vector{RGBAf}` computed on the Julia side (which is what closes the seams between tessels
and what keeps the colour Observable's element type stable — see the note at the top of this
file), so the colormap Observable feeds only the colorbar. The mapping has to be redone, which
is why every canvas keeps the values it last drew.

Unknown names are ignored rather than thrown: this arrives from QML, and a typo there should
not take the window down.
"""
function set_colormap!(c, name::AbstractString)
    i = findfirst(p -> first(p) == String(name), SURFACE_COLORMAPS)
    i === nothing && return false
    spec = last(SURFACE_COLORMAPS[i])
    cmap = endswith(spec, "_r") ?
        Makie.cgrad(_padded_cmap(spec[1:end-2]), rev = true) : _padded_cmap(spec)
    c.colormap[] = cmap
    _recolor!(c)
    return true
end

# Re-run the value -> RGBA mapping with the new colormap. Each canvas keeps its last values for
# exactly this: nothing else in the session can supply them, and re-deriving them would mean
# rebuilding the geometry.
function _recolor!(c::SkyCanvas)
    isempty(c.lastvalues[]) && return c
    lo, hi = c.cbarlimits[]
    cols = map_colors(c.lastvalues[], c.colormap[], (Float64(lo), Float64(hi)))
    c.colors[] = cols
    c.strokecolors[] = c.laststroke[] ?
        fill(Makie.RGBAf(0.45, 0.45, 0.45, 1), length(cols)) : cols
    return c
end
_recolor!(c::StarCanvas)  = c    # the mesh keeps numeric colours; the Observable is enough
_recolor!(c::MollCanvas)  = c    # redrawn from `lastmap` by the caller when it matters
_recolor!(c::Chi2Canvas)  = c    # categorical colours, not a ramp

# ── the posterior ────────────────────────────────────────────────────────────────────────
#
# TWO AXES, FIXED, and not a corner plot — which is what a reader expects and what the
# plot-once constraint forbids. A corner grid has D×D panels and D is the number of free
# parameters, decided when the user ticks a box; creating axes after the window exists is the
# GL failure this whole file is arranged around, and pre-building an 8×8 grid against the
# chance that eight parameters are freed would cost sixty-four axes at every start to show ten.
#
# So: one MARGINAL, for the parameter a selector picks, and one PAIR scatter for two. That is
# what the corner plot is read for anyway — is this marginal a single narrow peak or something
# else, and are these two degenerate with each other — and in a panel this size it is legible,
# which a squashed 8×8 is not.

struct PostCanvas
    figure::Any
    axis1::Any                                    # the marginal
    axis2::Any                                    # the pair scatter
    dens::Makie.Observable{Vector{Makie.Point2f}} # density polyline
    densplot::Any
    band::Makie.Observable{Vector{Makie.Point2f}} # the 16-84 interval, as a filled band
    bandplot::Any
    medline::Makie.Observable{Vector{Makie.Point2f}}
    medplot::Any
    pair::Makie.Observable{Vector{Makie.Point2f}}
    pairplot::Any
    pairmed::Makie.Observable{Vector{Makie.Point2f}}
    pairmedplot::Any
    message::Makie.Observable{String}
    messageplot::Any
end

function build_post_canvas(fig)
    ax1 = Makie.Axis(fig[1, 1]; xlabel = "parameter", ylabel = "posterior density")
    ax2 = Makie.Axis(fig[1, 2]; xlabel = "x", ylabel = "y")
    Makie.colsize!(fig.layout, 1, Makie.Relative(0.5))

    # The band FIRST, so the density line and the median draw over it.
    band = Makie.Observable([Makie.Point2f(0, 0), Makie.Point2f(0, 0), Makie.Point2f(0, 0)])
    bandplot = Makie.poly!(ax1, band; color = Makie.RGBAf(0.20, 0.45, 0.80, 0.20),
                           strokewidth = 0)
    dens = Makie.Observable([Makie.Point2f(0, 0)])
    densplot = Makie.lines!(ax1, dens; color = Makie.RGBAf(0.10, 0.25, 0.55, 1), linewidth = 2)
    medline = Makie.Observable([Makie.Point2f(0, 0), Makie.Point2f(0, 0)])
    medplot = Makie.lines!(ax1, medline; color = Makie.RGBAf(0.75, 0.15, 0.15, 1),
                           linewidth = 1.5, linestyle = :dash)

    pair = Makie.Observable([Makie.Point2f(0, 0)])
    pairplot = Makie.scatter!(ax2, pair; markersize = 4,
                              color = Makie.RGBAf(0.10, 0.25, 0.55, 0.35))
    pairmed = Makie.Observable([Makie.Point2f(0, 0)])
    pairmedplot = Makie.scatter!(ax2, pairmed; marker = :cross, markersize = 14,
                                 color = Makie.RGBAf(0.75, 0.15, 0.15, 1))

    msg  = Makie.Observable("")
    msgp = Makie.text!(ax1, [Makie.Point2f(0, 0)]; text = msg, fontsize = 15,
                       color = Makie.RGBAf(0.45, 0.5, 0.55, 1), align = (:center, :center))
    c = PostCanvas(fig, ax1, ax2, dens, densplot, band, bandplot, medline, medplot,
                   pair, pairplot, pairmed, pairmedplot, msg, msgp)
    idle!(c, "no posterior yet — run a sampler")
    return c
end

function idle!(c::PostCanvas, msg::AbstractString)
    # EMPTY, not a point at the origin. `lines!` and `scatter!` both accept an empty vector,
    # and a single placeholder point is not invisible: the median marker is a cross at
    # markersize 14, so an idle panel drew a red cross at (0, 0) that reads as data. The band
    # keeps its degenerate triangle — `poly!` is the one primitive here that an empty vector
    # is not obviously safe for, and a zero-area triangle draws nothing.
    c.dens[]     = Makie.Point2f[]
    c.band[]     = [Makie.Point2f(0, 0), Makie.Point2f(0, 0), Makie.Point2f(0, 0)]
    c.medline[]  = Makie.Point2f[]
    c.pair[]     = Makie.Point2f[]
    c.pairmed[]  = Makie.Point2f[]
    c.message[]  = String(msg)
    c.messageplot.visible[] = true
    for ax in (c.axis1, c.axis2)
        Makie.limits!(ax, -1, 1, -1, 1)
        ax.xlabel = ""; ax.ylabel = ""
    end
    c.axis1.xlabel = "parameter"; c.axis1.ylabel = "posterior density"
    return c
end

_recolor!(::PostCanvas) = nothing
_colorbar_visible!(::PostCanvas, ::Bool) = nothing   # two axes, no colorbar

"""
    show_posterior!(canvas, samples, names, i, j, med, q16, q84)

Draw the marginal of parameter `i` and the (`i`, `j`) pair.

The density is a HISTOGRAM drawn as a polyline rather than a kernel estimate: a KDE has a
bandwidth, the bandwidth changes the shape, and a posterior panel that quietly smooths a
bimodal marginal into one bump would defeat the reason this exists. Bin count follows
Freedman–Diaconis, floored so that a short chain still shows something.
"""
function show_posterior!(c::PostCanvas, samples::AbstractMatrix, names, i::Int, j::Int,
                         med, q16, q84)
    if isempty(samples) || size(samples, 1) < 2
        idle!(c, "this fit has no posterior — it was an optimiser, not a sampler")
        return c
    end
    i = clamp(i, 1, size(samples, 2)); j = clamp(j, 1, size(samples, 2))
    x = view(samples, :, i)
    lo, hi = extrema(x)
    if !(hi > lo)
        idle!(c, "the marginal is a single value")
        return c
    end
    n = length(x)
    iqr = Statistics.quantile(x, 0.75) - Statistics.quantile(x, 0.25)
    h = iqr > 0 ? 2 * iqr / cbrt(n) : (hi - lo) / 20
    nb = clamp(h > 0 ? ceil(Int, (hi - lo) / h) : 20, 8, 60)
    edges = range(lo, hi; length = nb + 1)
    counts = zeros(Float64, nb)
    for v in x
        k = clamp(floor(Int, (v - lo) / (hi - lo) * nb) + 1, 1, nb)
        counts[k] += 1
    end
    counts ./= (n * (edges[2] - edges[1]))          # a density, so the y scale means something
    # A step outline: two points per bin, so the bars are visible as bars.
    pts = Makie.Point2f[]
    push!(pts, Makie.Point2f(edges[1], 0))
    for k in 1:nb
        push!(pts, Makie.Point2f(edges[k], counts[k]))
        push!(pts, Makie.Point2f(edges[k + 1], counts[k]))
    end
    push!(pts, Makie.Point2f(edges[end], 0))
    c.dens[] = pts
    ymax = max(maximum(counts), eps(Float32))

    # The 16-84 band as an explicit quad, and the median as a line: both are what a reader
    # takes off a marginal, and neither survives being left to the eye on a bare histogram.
    a = isfinite(q16) ? Float32(q16) : Float32(lo)
    b = isfinite(q84) ? Float32(q84) : Float32(hi)
    c.band[] = [Makie.Point2f(a, 0), Makie.Point2f(b, 0),
                Makie.Point2f(b, ymax), Makie.Point2f(a, ymax)]
    mv = isfinite(med) ? Float32(med) : Float32(Statistics.median(x))
    c.medline[] = [Makie.Point2f(mv, 0), Makie.Point2f(mv, ymax)]

    c.pair[] = [Makie.Point2f(samples[k, i], samples[k, j]) for k in 1:size(samples, 1)]
    c.pairmed[] = [Makie.Point2f(mv, Statistics.median(view(samples, :, j)))]

    c.message[] = ""
    c.messageplot.visible[] = false
    style_live_axis!(c.axis1); style_live_axis!(c.axis2)
    c.axis1.xlabel = String(names[i]); c.axis1.ylabel = "posterior density"
    c.axis2.xlabel = String(names[i]); c.axis2.ylabel = String(names[j])
    pad = 0.05f0 * (hi - lo)
    Makie.xlims!(c.axis1, lo - pad, hi + pad)
    Makie.ylims!(c.axis1, 0, 1.10f0 * ymax)
    Makie.autolimits!(c.axis2)
    return c
end
