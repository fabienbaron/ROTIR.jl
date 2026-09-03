# =======================================================================================
# Makie drawing layer for the spheroid plots
# =======================================================================================
# The Makie counterpart of src/oiplot_spheroid.jl. Both draw the same mesh, described by
# src/oiplot_spheroid_core.jl, so the two front-ends cannot disagree about geometry.
#
# NAMES CARRY A `_makie` SUFFIX. Both back-ends can be loaded at once — PythonPlot for the
# existing demos, Makie for the GUI — and a bare `plot2d` cannot have two implementations.
# OITOOLS solved it the same way (`uvplot_makie`, `imdisp_makie`, `plot_v2_makie`). When the
# matplotlib layer is retired the suffix goes and these take the plain names.
#
# Everything here BUILDS A FIGURE AND RETURNS IT. Nothing is interactive, nothing mentions a
# backend: `using CairoMakie` gets vector output with no GPU, `using GLMakie` gets a window,
# and the GUI reaches it through QMLMakie. That separation is what lets a headless script or a
# PackageCompiler build plot at all — which the matplotlib layer cannot do.
#
# The GUI does NOT call these. It pre-creates its plots once and then assigns Observables,
# because under QMLMakie inserting a plot into a live scene allocates GPU buffers with no GL
# context bound. See the plan's "create plots once" section.

# ---------------------------------------------------------------------------------------
# Colour
# ---------------------------------------------------------------------------------------
"""
    _padded_cmap(name; cfloor=0.08, cceil=0.95) -> ColorGradient

The matplotlib layer widens its `Normalize` so data lands in colormap fraction
`[cfloor, cceil]` rather than `[0,1]` (`_padded_norm`). The reason is not cosmetic: the
background is white and the top of `gist_heat` is also white, so an unpadded scale renders the
hottest surface in exactly the background colour — fatal on a binary, where both components
share one scale and the hotter one disappears.

In Makie the same thing is a colormap SLICE, which is simpler and leaves `colorrange` free to
carry the true data range — so the colorbar shows real temperatures with no separate clipping
step (matplotlib needed `cb.ax.set_ylim(pmin, pmax)` for that).
"""
_padded_cmap(name; cfloor = 0.08, cceil = 0.95) =
    Makie.cgrad(Makie.cgrad(Symbol(name))[range(cfloor, cceil, length = 256)])

# Font sizes in points, matched to the matplotlib layer's `set_oiplot_defaults()` at equal
# rendered size. Makie's own Axis defaults are smaller relative to the axis box, which makes an
# A/B against the matplotlib twin read as a different plot rather than as the same one. They
# live here, as plain constants, rather than in a Makie `Theme`: a theme is global state, and
# the GUI has to be able to build a figure without having reached into the caller's.
const PLOT_TITLESIZE     = 19
const PLOT_LABELSIZE     = 18
const PLOT_TICKLABELSIZE = 15
const PLOT_LEGENDSIZE    = 14

_axis(cell; kwargs...) = Makie.Axis(cell; titlesize = PLOT_TITLESIZE,
                                   xlabelsize = PLOT_LABELSIZE, ylabelsize = PLOT_LABELSIZE,
                                   xticklabelsize = PLOT_TICKLABELSIZE,
                                   yticklabelsize = PLOT_TICKLABELSIZE, kwargs...)

# Colorbar ticks as matplotlib writes them: N evenly spaced values with SHORT labels. Makie's
# default formatter prints the raw Float64 (6646.0439...), which is unreadable and, on a
# temperature map, spuriously precise.
function _bar_ticks(lo, hi; n = 7)
    v = collect(range(lo, hi, length = n))
    return (v, [Printf.@sprintf("%.0f", x) for x in v])
end

"""
    map_colors(values, cmap, colorrange) -> Vector{RGBAf}

Resolve values to explicit colours. Two reasons not to hand Makie the raw numbers:

 1. **Seams.** `poly!` antialiases each polygon against its neighbour, leaving pale hairlines
    that read as a spurious grid over the surface. The matplotlib layer avoids it by stroking
    every polygon in its OWN face colour, and that needs the colour as a value.
 2. **The GUI.** Swapping the element type of a `color` Observable forces Makie to rebuild the
    plot, which under QMLMakie means allocating GPU buffers with no context bound. A canvas
    that is always fed `Vector{RGBAf}` never changes type.
"""
function map_colors(values, cmap, colorrange)
    lo, hi = colorrange
    s = hi > lo ? 1 / (hi - lo) : 0.0
    return [Makie.RGBAf(Makie.to_color(cmap[clamp((v - lo) * s, 0, 1)])) for v in values]
end

"Value range for a map, with the degenerate/uniform guard the matplotlib layer uses."

"""
    _surface_values(tmap, star; intensity, intensity_model, band) -> Vector

The value actually drawn per visible tessel. `intensity` decides WHETHER limb darkening is
applied; `intensity_model` decides which temperature→brightness conversion is used. They are
different questions — with `intensity=false` a temperature map is band-independent, so
`:linear` and `:planck` must give identical output.
"""
function _surface_values(tmap, star; intensity::Bool = false,
                         intensity_model::Symbol = :linear, band = nothing)
    vis = star.index_quads_visible
    intensity || return Float64.(tmap[vis])
    Imap = intensity_model === :linear ? tmap : ROTIR.intensity(tmap, intensity_model, band)
    return Float64.(Imap[vis]) .* Float64.(star.ldmap[vis])
end

# ---------------------------------------------------------------------------------------
# The core primitive
# ---------------------------------------------------------------------------------------
"""
    tessel_polygons(star; offset_west=0.0, offset_north=0.0, indices=nothing)

The visible tessels as `Vector{Vector{Point2f}}`, ready for `poly!`.

`vertices_xyz` is a POLYGON SOUP — 4 corners per tessel with nothing shared between
neighbours (src/geometry.jl:24-42) — which is exactly what `poly!` wants, so this is a direct
translation of the matplotlib `PolyCollection` and not an approximation of it.

Plot convention: `x = East = -proj_west`. `indices` defaults to the front-facing tessels;
passing the complement draws the far hemisphere, which is what a see-through wireframe needs.
"""
function tessel_polygons(star; offset_west = 0.0, offset_north = 0.0, indices = nothing)
    vis = indices === nothing ? star.index_quads_visible : indices
    polys = Vector{Vector{Makie.Point2f}}(undef, length(vis))
    @inbounds for (i, idx) in enumerate(vis)
        polys[i] = [Makie.Point2f(-(Float64(star.proj_west[idx, c]) + offset_west),
                                    Float64(star.proj_north[idx, c]) + offset_north)
                    for c in 1:4]
    end
    return polys
end

"""
    add_tessel_collection_makie!(ax, star, values; ...) -> Makie.Poly

Draw the tessels into an existing axis. `values` is one number per drawn tessel, in the order
`tessel_polygons` produces — i.e. indexed by POSITION WITHIN `index_quads_visible`, not by
pixel. Getting that wrong scrambles the map with no error.

`plotmesh` strokes the mesh in mid grey rather than each polygon's own face colour. 0.45 grey
is deliberate: the edge has to read against BOTH ends of the colormap, and a light edge
vanishes against a pale surface, making `plotmesh=true` silently do nothing.
"""
function add_tessel_collection_makie!(ax, star, values;
                                      colormap = _padded_cmap("gist_heat"),
                                      colorrange = nothing, plotmesh = false,
                                      offset_west = 0.0, offset_north = 0.0,
                                      strokecolor = nothing, strokewidth = nothing,
                                      indices = nothing)
    polys = tessel_polygons(star; offset_west, offset_north, indices)
    cr = colorrange === nothing ? _map_range(values) : colorrange
    cols = map_colors(values, colormap, cr)
    # Non-mesh mode strokes each polygon in its OWN face colour — that is what closes the
    # antialiasing seams. Mid grey for the mesh, not "lightgrey": the edge has to read against
    # BOTH ends of the colormap, and a light edge vanishes against a pale surface, which would
    # make plotmesh=true silently do nothing.
    sc = strokecolor === nothing ? (plotmesh ? Makie.RGBAf(0.45, 0.45, 0.45, 1) : cols) :
                                   strokecolor
    sw = strokewidth === nothing ? (plotmesh ? 0.25 : 0.35) : strokewidth
    return Makie.poly!(ax, polys; color = cols, strokecolor = sc, strokewidth = sw)
end

# ---------------------------------------------------------------------------------------
# Decorations
# ---------------------------------------------------------------------------------------
"""
    draw_limb_makie!(ax, star; ...)

Convex hull of the visible corners. NOT `silhouette_polygon`: azimuthal binning notches
inward on a coarse mesh, which reads as a dented star.
"""
function draw_limb_makie!(ax, star; offset_west = 0.0, offset_north = 0.0,
                          color = :black, linewidth = 0.8, alpha = 1.0)
    vis = star.index_quads_visible
    hx, hy = convex_hull_2d(vec(-star.proj_west[vis, :]), vec(star.proj_north[vis, :]))
    pts = [Makie.Point2f(x - offset_west, y + offset_north)
           for (x, y) in zip(vcat(hx, hx[1]), vcat(hy, hy[1]))]
    return Makie.lines!(ax, pts; color = (color, alpha), linewidth = linewidth)
end

"""
    sky_axis_max(star; pad=0.5)

Half-width of the square view, from the mesh's own extent. Matches the matplotlib layer so
the two produce identically framed figures.
"""
sky_axis_max(star; pad = 0.5) =
    maximum(sqrt.(star.vertices_xyz[:, :, 1].^2 .+ star.vertices_xyz[:, :, 2].^2 .+
                  star.vertices_xyz[:, :, 3].^2)) + pad

"East-left, North-up, equal aspect — the astronomical convention, as reversed limits."
function style_sky_axis!(ax, amax; flipx = false)
    ax.aspect = Makie.DataAspect()
    ax.xlabel = "x ← E (mas)"
    ax.ylabel = "y → N (mas)"
    # No grid: the matplotlib twin has none, and over a filled disk it reads as an artefact of
    # the map rather than as an axis decoration.
    ax.xgridvisible = false
    ax.ygridvisible = false
    Makie.xlims!(ax, flipx ? (-amax, amax) : (amax, -amax))
    Makie.ylims!(ax, -amax, amax)
    return ax
end

"""
    draw_rotation_axis_makie!(ax, star; ...)

The spin axis as three dashed segments — south tip, the interior run at reduced alpha, north
tip — with an arrowhead at the north end.

The interior segment is faded on purpose: at full alpha a dashed line straight across the disk
reads as a feature ON the surface rather than an axis passing behind it.

The arrowhead's rotation is computed in SCREEN space. Plot coordinates negate west
(`x = -proj_west`) and the axis is then reversed for East-left, and those two cancel — so the
screen direction is `(delta_west, delta_north)`, not its negation. Getting this wrong points
the arrow at the south pole, which looks plausible enough to miss.
"""
function draw_rotation_axis_makie!(ax, star; arrow_frac = 0.3, color = :black,
                                   linewidth = 1.5, offset_west = 0.0, offset_north = 0.0,
                                   inclination = NaN, position_angle = NaN,
                                   alpha = 1.0, inside_alpha = 0.5, star_params = nothing)
    north, south = _spin_axis(star, star_params, inclination, position_angle)
    delta = north .- south
    ntip  = north .+ arrow_frac .* delta
    stip  = south .- arrow_frac .* delta
    P(p)  = Makie.Point2f(-(p[1] + offset_west), p[2] + offset_north)
    for (p1, p2, a) in ((stip, south, alpha), (south, north, inside_alpha),
                        (north, ntip, alpha))
        Makie.lines!(ax, [P(p1), P(p2)]; color = (color, a), linewidth = linewidth,
                     linestyle = :dash)
    end
    Makie.scatter!(ax, [P(ntip)]; marker = :utriangle, markersize = 12, color = color,
                   rotation = atan(delta[2], delta[1]) - π/2)
    return ax
end

"""
    draw_rotation_arrow_makie!(ax, star; pole="N", ...)

A 300° circular arrow around the pole, showing the sense of rotation.

The arc is split at `z = 0`: the half in front of the star is solid, the half behind it dashed,
so the arrow reads as encircling the pole rather than lying flat on it. The matplotlib layer
does this by NaN-masking and two fixed zorders; here it is two `lines!` calls in draw order.

The head is placed on the TANGENT at the endpoint, never on a chord back to an earlier sample:
over a 300° arc a chord points visibly the wrong way. Makie's screen-space markers also remove
the `mutation_scale` back-computation the matplotlib version needs to keep the head from
chording straight across the arc.
"""
function draw_rotation_arrow_makie!(ax, star; pole = "N", radius_frac = 0.07,
                                    offset_frac = 0.15, color = :black, linewidth = 1.5,
                                    npoints = 100, offset_west = 0.0, offset_north = 0.0,
                                    inclination = NaN, position_angle = NaN,
                                    star_params = nothing)
    north, south = _spin_axis(star, star_params, inclination, position_angle)
    axis = north .- south
    alen = LinearAlgebra.norm(axis)
    ahat = axis ./ alen
    centre = pole == "N" ? north .+ offset_frac .* axis : south .- offset_frac .* axis
    r = radius_frac * alen
    ref = abs(ahat[3]) < 0.9 ? [0.0, 0.0, 1.0] : [1.0, 0.0, 0.0]
    e1 = LinearAlgebra.cross(ahat, ref); e1 ./= LinearAlgebra.norm(e1)
    e2 = LinearAlgebra.cross(ahat, e1)
    θ  = collect(range(0, 300π/180, length = npoints))
    P3 = [centre .+ r .* (cos(t) .* e1 .+ sin(t) .* e2) for t in θ]
    xs = [-(p[1] + offset_west) for p in P3]
    ys = [ p[2] + offset_north  for p in P3]
    zs = [ p[3]                 for p in P3]
    front = zs .> 0
    # Two passes, NaN-broken so each stays one plot object.
    for (keep, style) in ((front, :solid), (.!front, :dash))
        pts = [k ? Makie.Point2f(x, y) : Makie.Point2f(NaN, NaN)
               for (k, x, y) in zip(keep, xs, ys)]
        Makie.lines!(ax, pts; color = color, linewidth = linewidth, linestyle = style)
    end
    dpt = -sin(θ[end]) .* e1 .+ cos(θ[end]) .* e2      # tangent, not a chord
    Makie.scatter!(ax, [Makie.Point2f(xs[end], ys[end])]; marker = :utriangle,
                   markersize = 12, color = color,
                   rotation = atan(dpt[2], dpt[1]) - π/2)
    return ax
end

"""
    draw_graticules_makie!(ax, star; ...)

Parallels and meridians, plus (by default) the limb. The geometry comes from
`graticule_segments` in the core, shared with the matplotlib layer — this function only draws.

Segments are joined into ONE line plot with NaN separators rather than one `lines!` per run:
Makie breaks a polyline at NaN, so a graticule with a hundred occlusion-clipped runs stays a
single plot object. That matters for the GUI, where the plot pool is fixed at construction.
"""
function draw_graticules_makie!(ax, star; nlat = 5, nlon = 8, color = :black,
                                linewidth = 0.8, alpha = 0.5,
                                offset_west = 0.0, offset_north = 0.0,
                                inclination = NaN, position_angle = NaN,
                                npoints = 200, star_params = nothing, limb = true)
    segs = graticule_segments(star; nlat, nlon, offset_west, offset_north,
                              inclination, position_angle, npoints, star_params)
    limb && draw_limb_makie!(ax, star; offset_west, offset_north,
                             color = color, linewidth = linewidth, alpha = alpha)
    isempty(segs) && return ax
    pts = Makie.Point2f[]
    for s in segs
        append!(pts, (Makie.Point2f(s[i, 1], s[i, 2]) for i in axes(s, 1)))
        push!(pts, Makie.Point2f(NaN, NaN))          # break, not a join to the next run
    end
    Makie.lines!(ax, pts; color = (color, alpha), linewidth = linewidth)
    return ax
end

"""
    draw_compass_makie!(ax, amax; size_frac=0.12, fontsize=13, color=:black)

N/E arrows. Placed at the same data coordinates as the matplotlib version so the two figures
agree; note that with the reversed x-axis (East-left) increasing x draws to the LEFT, which is
why the E arrowhead is `:ltriangle` and not `:rtriangle` — Makie markers are screen-space and
do not flip with the axis.

The matplotlib version derives the arrow length in POINTS and converts back to data units, so
a compass in a small subplot does not shrink below legibility (`_axes_height_inches`). Makie
has screen-space markers and per-axis font sizes, so that back-computation is not reproduced
here; if `plot2d_allepochs_makie` ever needs it, scale `size_frac` per panel instead.
"""
function draw_compass_makie!(ax, amax; size_frac = 0.09, fontsize = 12, color = :black)
    s      = size_frac * amax
    margin = 0.20 * amax
    cx, cy = -(amax - margin), -(amax - margin)
    gap    = 0.06 * amax
    # N: +y. E: +x, which the reversed axis renders leftward.
    Makie.lines!(ax, [Makie.Point2f(cx, cy), Makie.Point2f(cx, cy + s)];
                 color = color, linewidth = 1.5)
    Makie.scatter!(ax, [Makie.Point2f(cx, cy + s)]; marker = :utriangle,
                   markersize = 11, color = color)
    Makie.text!(ax, cx - gap, cy + s; text = "N", color = color, fontsize = fontsize,
                font = :bold, align = (:center, :bottom))
    Makie.lines!(ax, [Makie.Point2f(cx, cy), Makie.Point2f(cx + s, cy)];
                 color = color, linewidth = 1.5)
    Makie.scatter!(ax, [Makie.Point2f(cx + s, cy)]; marker = :ltriangle,
                   markersize = 11, color = color)
    Makie.text!(ax, cx + s, cy - gap; text = "E", color = color, fontsize = fontsize,
                font = :bold, align = (:right, :center))
    return ax
end

# ---------------------------------------------------------------------------------------
# Entry points
# ---------------------------------------------------------------------------------------
"""
    plot2d_makie(tmap, star; ...) -> (fig, ax)

Sky-plane view of the surface map. The Makie counterpart of [`plot2d`](@ref); same framing,
same colour treatment, same East-left convention.

`intensity=true` multiplies the limb-darkening map into the colour, giving the physically
apparent disk; `intensity=false` colours by raw temperature with a hard limb edge.
"""
function plot2d_makie(tmap, star; intensity = false, figtitle = "", plotmesh = false,
                      pad = 0.5, colormap = "gist_heat", background = :white, flipx = false,
                      compass = true, graticules = false, limb = true,
                      rotation_axis = false, rotation_arrow = false,
                      limb_color = :black, limb_linewidth = 0.8,
                      inclination = NaN, position_angle = NaN, star_params = nothing,
                      intensity_model::Symbol = :linear, band = nothing,
                      vmin = nothing, vmax = nothing, colorbar = true,
                      figsize = (760, 700), graticule_kwargs = (;))
    vals = _surface_values(tmap, star; intensity, intensity_model, band)
    pmin, pmax = _map_range(vals; vmin, vmax)
    amax = sky_axis_max(star; pad)

    fig = Makie.Figure(size = figsize, backgroundcolor = :white)
    ax  = _axis(fig[1, 1]; title = figtitle, backgroundcolor = background)
    style_sky_axis!(ax, amax; flipx)

    cmap = _padded_cmap(colormap)
    add_tessel_collection_makie!(ax, star, vals; colormap = cmap,
                                 colorrange = (pmin, pmax), plotmesh)
    limb && draw_limb_makie!(ax, star; color = limb_color, linewidth = limb_linewidth)
    graticules && draw_graticules_makie!(ax, star; star_params = star_params,
                                         inclination, position_angle, limb = false,
                                         graticule_kwargs...)
    rotation_axis && draw_rotation_axis_makie!(ax, star; star_params, inclination,
                                               position_angle)
    rotation_arrow && draw_rotation_arrow_makie!(ax, star; star_params, inclination,
                                                 position_angle)
    compass && draw_compass_makie!(ax, amax)
    # colorrange is the TRUE data range — the padding lives in the colormap slice — so the
    # colorbar needs no clipping step.
    colorbar && Makie.Colorbar(fig[1, 2]; colormap = cmap, colorrange = (pmin, pmax))
    return fig, ax
end

"""
    star_mesh(star; offset=(0,0,0), values=nothing) -> (mesh, vertexcolors)

The full closed surface as a `GeometryBasics.Mesh`, ready for `mesh!`.

Two things here are deliberate and both are wrong if changed:

**Vertices are NOT welded.** `vertices_xyz` is a polygon soup — 4 corners per tessel, nothing
shared with the neighbours (src/geometry.jl:24-42) — and it stays that way, so each tessel
carries one flat colour on all four of its corners. Welding would make Makie interpolate
across tessel boundaries and quietly smooth the map that the whole reconstruction is about.

**The sky frame is used as-is, West included.** `(West, North, toward-observer)` is already
RIGHT-handed — `Ŵ = N̂ × ẑ_obs`, since `Ê × N̂` points AWAY from the observer — so no axis is
flipped here. Negating x to "make it East" mirrors the star and reverses its apparent spin.
Handedness is also exactly why nothing is flipped: with the camera on +z and North up, screen
right is `N̂ × ẑ_obs` = West, which puts East on the left. That is the astronomical convention
`plot2d_makie` draws by reversing its x-axis, reached here through the geometry instead.
"""
function star_mesh(star; offset = (0.0, 0.0, 0.0), values = nothing)
    npix  = star.npix
    verts = Vector{Makie.Point3f}(undef, 4npix)
    faces = Vector{Makie.TriangleFace{Int}}(undef, 2npix)
    ox, oy, oz = offset
    @inbounds for i in 1:npix
        b = 4(i - 1)
        for c in 1:4
            verts[b + c] = Makie.Point3f( star.vertices_xyz[i, c, 1] + ox,
                                          star.vertices_xyz[i, c, 2] + oy,
                                          star.vertices_xyz[i, c, 3] + oz)
        end
        # Quad (1,2,3,4) -> two triangles sharing the 1-3 diagonal, the same split
        # `finish_star` assumes when it takes the AC x BD cross product for the normal.
        faces[2i - 1] = Makie.TriangleFace(b + 1, b + 2, b + 3)
        faces[2i]     = Makie.TriangleFace(b + 1, b + 3, b + 4)
    end
    msh = Makie.GeometryBasics.Mesh(verts, faces)
    values === nothing && return msh, nothing
    vcol = Vector{eltype(values)}(undef, 4npix)
    @inbounds for i in 1:npix, c in 1:4
        vcol[4(i - 1) + c] = values[i]
    end
    return msh, vcol
end

"""
    scene3d(fig, cell; axistype=:lscene, title="", radius=1.0) -> (parent, axis-like)

The container a 3-D surface is drawn into.

`:lscene` is the default because it is the one the mouse behaves properly in: left-drag
orbits freely about the target, scroll zooms, right-drag pans, and there are no ticks or box
in the way of the star. `:axis3` keeps Makie's boxed axis with labelled mas ticks — better
when the point is to read a size off the plot rather than to look at the shape.

Either way the camera starts where the observer actually is: on +z, looking back at the
origin, North up — which, the sky frame being right-handed, puts East on the left. That is
exactly the view `plot2d_makie` draws, so the 3-D view opens as the 2-D one and the mouse
takes it from there.

The projection is ORTHOGRAPHIC, not Makie's default perspective. The observer is effectively
at infinity, so perspective would both disagree with the 2-D view at the opening camera and
inflate whichever component happens to be nearer — reading as a size difference between the
two stars that is not there.
"""
function scene3d(fig, cell = fig[1, 1]; axistype::Symbol = :lscene, title = "",
                 radius::Real = 1.0)
    if axistype === :axis3
        ax = Makie.Axis3(cell; title = title, aspect = :data,
                         xlabel = "x → W (mas)", ylabel = "y → N (mas)",
                         zlabel = "z → obs (mas)",
                         # Forty-five degrees off the line of sight. Straight down it is what
                         # the 2-D view already shows, and a 3-D scene opening there reads as
                         # a flat picture with axes round it; at exactly elevation = pi/2 the
                         # Axis3 up-vector is degenerate as well.
                         azimuth = pi/2, elevation = pi/2 - deg2rad(45))
        return ax, ax
    end
    axistype === :lscene || throw(ArgumentError("axistype must be :lscene or :axis3, got :$axistype"))
    ls = Makie.LScene(cell; show_axis = false)
    Makie.cam3d!(ls.scene; projectiontype = Makie.Orthographic)
    d = 3.2 * radius
    Makie.update_cam!(ls.scene, Makie.Vec3f(0, 0, d), Makie.Vec3f(0, 0, 0), Makie.Vec3f(0, 1, 0))
    isempty(title) || Makie.Label(fig[0, 1], title; fontsize = 16, font = :bold,
                                  tellwidth = false)
    return ls, ls
end

"""
    plot3d_makie(tmap, star; ...) -> (fig, ax)

The surface as a real 3-D object, rotatable with the mouse under GLMakie.

Unlike [`plot2d_makie`](@ref) this draws EVERY tessel, front and back: occlusion comes from
the depth buffer instead of from `index_quads_visible`, which is what makes the far side
appear correctly as the star is turned.

`shading = NoShading` throughout, and that is physics, not economy. The colours already are
the emergent surface brightness — gravity darkening, spots, limb darkening if `intensity` is
on — so letting Makie add a light would multiply a second, invented illumination on top of
the map being displayed.
"""
function plot3d_makie(tmap, star; colormap = "gist_heat", figtitle = "",
                      intensity = false, intensity_model::Symbol = :linear, band = nothing,
                      vmin = nothing, vmax = nothing, colorbar = true, plotmesh = false,
                      axistype::Symbol = :lscene, figsize = (760, 700))
    vals = _surface_values3d(tmap, star; intensity, intensity_model, band)
    pmin, pmax = _map_range(vals; vmin, vmax)
    cmap = _padded_cmap(colormap)
    msh, vcol = star_mesh(star; values = vals)

    fig = Makie.Figure(size = figsize, backgroundcolor = :white)
    parent, ax = scene3d(fig, fig[1, 1]; axistype, title = figtitle,
                         radius = _body_radius(star))
    Makie.mesh!(ax, msh; color = vcol, colormap = cmap, colorrange = (pmin, pmax),
                shading = Makie.NoShading)
    plotmesh && Makie.wireframe!(ax, msh; color = Makie.RGBAf(0.45, 0.45, 0.45, 1),
                                 linewidth = 0.4, depth_shift = -1f-3)
    colorbar && Makie.Colorbar(fig[1, 2]; colormap = cmap, colorrange = (pmin, pmax))
    return fig, ax
end

"""
    plot3d_binary_makie(tmap1, tmap2, star1, star2, bparams, tepoch; ...) -> (fig, ax)

Both components in one rotatable 3-D scene, placed at their orbital separation.

The separation is the full 3-D one, not the sky-projected pair `plot2d_binary_makie` uses:
the line-of-sight term is what makes an eclipse readable as one star passing in FRONT of the
other once the view is tilted away from the observer's. `binary_orbit_abs` returns
`(North, East, receding)`, so it is remapped to the mesh frame `(West, North, toward)` —
mixing the two up puts the secondary on the wrong side and, at conjunction, behind instead
of in front.

The primary sits at the origin; only the relative vector matters for the picture, and
holding it fixed keeps the camera from drifting as the epoch is stepped.
"""
function plot3d_binary_makie(tmap1, tmap2, star1, star2, bparams, tepoch;
                             colormap = "gist_heat", figtitle = "", intensity = false,
                             intensity_model::Symbol = :linear, band = nothing,
                             vmin = nothing, vmax = nothing, colorbar = true,
                             plotmesh = false, orbit = true,
                             orbit_color = Makie.RGBAf(0.35, 0.35, 0.35, 0.8),
                             axistype::Symbol = :lscene, figsize = (860, 700))
    v1 = _surface_values3d(tmap1, star1; intensity, intensity_model, band)
    v2 = _surface_values3d(tmap2, star2; intensity, intensity_model, band)
    tmin = vmin === nothing ? min(minimum(v1), minimum(v2)) : vmin
    tmax = vmax === nothing ? max(maximum(v1), maximum(v2)) : vmax
    tmax - tmin < 1.0 && (tmax = tmin + max(abs(tmax) * 0.01, 1.0))
    cmap = _padded_cmap(colormap; cfloor = 0.15)

    x1, y1, z1, x2, y2, z2 = binary_orbit_abs(bparams, tepoch)
    off = (-(y2 - y1), x2 - x1, -(z2 - z1))   # (West, North, toward observer)

    fig = Makie.Figure(size = figsize, backgroundcolor = :white)
    sep = sqrt(off[1]^2 + off[2]^2 + off[3]^2)
    parent, ax = scene3d(fig, fig[1, 1]; axistype, title = figtitle,
                         radius = max(sep, _body_radius(star1) + _body_radius(star2)))
    for (st, vals, o) in ((star1, v1, (0.0, 0.0, 0.0)), (star2, v2, off))
        msh, vcol = star_mesh(st; offset = o, values = vals)
        Makie.mesh!(ax, msh; color = vcol, colormap = cmap, colorrange = (tmin, tmax),
                    shading = Makie.NoShading)
        plotmesh && Makie.wireframe!(ax, msh; color = Makie.RGBAf(0.45, 0.45, 0.45, 1),
                                     linewidth = 0.4, depth_shift = -1f-3)
    end
    orbit && Makie.lines!(ax, relative_orbit_track(bparams); color = orbit_color,
                          linewidth = 1.2)
    colorbar && Makie.Colorbar(fig[1, 2]; colormap = cmap, colorrange = (tmin, tmax))
    return fig, ax
end

"""
    relative_orbit_track(bparams; npoints=360) -> Vector{Point3f}

The secondary's path around the primary over one period, in the mesh frame.

Drawn as a closed curve rather than sampled per epoch: with apsidal motion (`dω`) the track
does not quite close, and showing that is more honest than hiding it.
"""
function relative_orbit_track(bparams; npoints::Int = 360)
    P  = bparams.P
    t0 = hasproperty(bparams, :T0) ? bparams.T0 : 0.0
    pts = Vector{Makie.Point3f}(undef, npoints)
    @inbounds for (k, f) in enumerate(range(0, 1, length = npoints))
        x1, y1, z1, x2, y2, z2 = binary_orbit_abs(bparams, t0 + f * P)
        pts[k] = Makie.Point3f(-(y2 - y1), x2 - x1, -(z2 - z1))
    end
    return pts
end

# Same question as `_surface_values`, different answer: the 3-D view draws ALL tessels, so
# the value vector must be indexed by PIXEL, not by position within `index_quads_visible`.
function _surface_values3d(tmap, star; intensity::Bool = false,
                           intensity_model::Symbol = :linear, band = nothing)
    intensity || return Float64.(tmap)
    Imap = intensity_model === :linear ? tmap : ROTIR.intensity(tmap, intensity_model, band)
    return Float64.(Imap) .* Float64.(star.ldmap)
end


"""
    plot2d_wireframe_makie(star; ...) -> (fig, ax)

Quad edges only. With `hidden=false` the far hemisphere is drawn too, faintly — the
"transparent globe" view that shows how the mesh wraps.
"""
function plot2d_wireframe_makie(star; hidden = true, pad = 0.5, figtitle = "",
                                color = :black, linewidth = 0.4, flipx = false,
                                compass = true, figsize = (700, 700))
    amax = sky_axis_max(star; pad)
    fig = Makie.Figure(size = figsize, backgroundcolor = :white)
    ax  = _axis(fig[1, 1]; title = figtitle)
    style_sky_axis!(ax, amax; flipx)
    if !hidden
        far = setdiff(1:star.npix, star.index_quads_visible)
        Makie.poly!(ax, tessel_polygons(star; indices = far); color = :transparent,
                    strokecolor = (color, 0.18), strokewidth = linewidth)
    end
    Makie.poly!(ax, tessel_polygons(star); color = :transparent,
                strokecolor = color, strokewidth = linewidth)
    compass && draw_compass_makie!(ax, amax)
    return fig, ax
end

"""
    plot2d_binary_makie(tmap1, tmap2, star1, star2, bparams, tepoch; ...) -> (fig, ax)

Both components at their orbital positions on ONE shared colour scale.

The shared scale is why the padding floor is 0.15 here against 0.08 for a single star: the
hotter component lands at the pale end of the colormap, and against a white background it
would otherwise vanish — taking its limb and graticules into apparently empty space with it.

Depth comes from the orbit (`binary_orbit_abs` z, `+z` receding), and layering is simply DRAW
ORDER: far surface, far limb, far graticules, then the near component's three, then the
compass. The matplotlib layer needs fractional zorders (2, 2.4, 2.5 / 3, 3.4, 3.5) to express
the same thing, because painting every limb at one fixed z draws the far star's outline over
the near one.
"""
function plot2d_binary_makie(tmap1, tmap2, star1, star2, bparams, tepoch;
                             intensity = false, plotmesh = false, colormap = "gist_heat",
                             pad = 1.0, background = :white, compass = true,
                             graticules = false, figtitle = "",
                             inclination1 = NaN, position_angle1 = NaN,
                             inclination2 = NaN, position_angle2 = NaN,
                             star_params1 = nothing, star_params2 = nothing,
                             limb = true, limb_color = :black, limb_linewidth = 0.8,
                             intensity_model::Symbol = :linear, band = nothing,
                             vmin = nothing, vmax = nothing, amax = nothing,
                             colorbar = true, figsize = (760, 700), graticule_kwargs = (;))
    v1 = _surface_values(tmap1, star1; intensity, intensity_model, band)
    v2 = _surface_values(tmap2, star2; intensity, intensity_model, band)
    tmin = vmin === nothing ? min(minimum(v1), minimum(v2)) : vmin
    tmax = vmax === nothing ? max(maximum(v1), maximum(v2)) : vmax
    tmax - tmin < 1.0 && (tmax = tmin + max(abs(tmax) * 0.01, 1.0))

    offset_west, offset_north = orbit_to_rotir_offset(bparams, tepoch)
    if amax === nothing
        r1 = maximum(sqrt.(star1.vertices_xyz[:,:,1].^2 .+ star1.vertices_xyz[:,:,2].^2))
        r2 = maximum(sqrt.(star2.vertices_xyz[:,:,1].^2 .+ star2.vertices_xyz[:,:,2].^2))
        amax = max(r1, r2 + abs(offset_west), r2 + abs(offset_north)) + pad
    end

    fig = Makie.Figure(size = figsize, backgroundcolor = :white)
    ax  = _axis(fig[1, 1]; title = figtitle, backgroundcolor = background)
    style_sky_axis!(ax, amax)
    cmap = _padded_cmap(colormap; cfloor = 0.15)

    # Farther component first. z is +receding, so the LARGER z is behind.
    _, _, z1, _, _, z2 = binary_orbit_abs(bparams, tepoch)
    far_first = z1 > z2
    comps = ((star1, v1, 0.0, 0.0, star_params1, inclination1, position_angle1),
             (star2, v2, offset_west, offset_north, star_params2, inclination2, position_angle2))
    for (st, vals, ow, on, sp, inc, pa) in (far_first ? comps : reverse(comps))
        add_tessel_collection_makie!(ax, st, vals; colormap = cmap,
                                     colorrange = (tmin, tmax), plotmesh,
                                     offset_west = ow, offset_north = on)
        limb && draw_limb_makie!(ax, st; offset_west = ow, offset_north = on,
                                 color = limb_color, linewidth = limb_linewidth)
        graticules && draw_graticules_makie!(ax, st; star_params = sp, inclination = inc,
                                             position_angle = pa, offset_west = ow,
                                             offset_north = on, limb = false,
                                             graticule_kwargs...)
    end
    compass && draw_compass_makie!(ax, amax)
    colorbar && Makie.Colorbar(fig[1, 2]; colormap = cmap, colorrange = (tmin, tmax))
    return fig, ax
end

"""
    plot2d_allepochs_makie(tmap, stars; ...) -> fig

One panel per epoch on a shared colour scale.

Note the scale is plain min–max over the WHOLE map with no padding and no colorbar, matching
the matplotlib version — this view is for comparing epochs against each other, not for reading
absolute temperatures off a bar.
"""
function plot2d_allepochs_makie(tmap, stars; plotmesh = false, tepochs = [],
                                colormap = "gist_heat", ncols = nothing, compass = true,
                                pad = 0.5, panelsize = 300)
    nep = length(stars)
    nc  = ncols === nothing ? min(nep, 4) : ncols
    nr  = cld(nep, nc)
    minT, maxT = extrema(tmap)
    minT == maxT && ((minT, maxT) = (0.95minT, 1.05maxT))   # uniform-map guard

    fig  = Makie.Figure(size = (panelsize * nc + 60, panelsize * nr + 60),
                        backgroundcolor = :white)
    cmap = _padded_cmap(colormap)
    amax = maximum(sky_axis_max(s; pad) for s in stars)     # one frame for every panel
    for t in 1:nep
        r, c = fldmod1(t, nc)
        ax = _axis(fig[r, c];
                        title = isempty(tepochs) ? "epoch $t" :
                                Printf.@sprintf("t = %.3f", tepochs[t]))
        style_sky_axis!(ax, amax)
        add_tessel_collection_makie!(ax, stars[t], Float64.(tmap[stars[t].index_quads_visible]);
                                     colormap = cmap, colorrange = (minT, maxT), plotmesh)
        # Shrink the compass with the panel. Makie font sizes are in points and do not scale
        # with the axis, so a compass sized only in data units keeps full-size text in a small
        # panel and the labels collide with their own arrows — the same failure the matplotlib
        # layer avoids by deriving the arrow length in points (`_axes_height_inches`).
        # First panel only, as the matplotlib version does: the orientation is the same in
        # every panel, and eight copies of it is eight pieces of noise.
        compass && t == 1 && draw_compass_makie!(ax, amax; size_frac = 0.10,
                                       fontsize = clamp(13 * 3 / max(nc, nr), 7, 13))
        # Only edge panels carry labels; a grid of repeated axis text is unreadable.
        r < nr && (ax.xlabel = ""; Makie.hidexdecorations!(ax; grid = false))
        c > 1  && (ax.ylabel = ""; Makie.hideydecorations!(ax; grid = false))
    end
    return fig
end

"""
    plot_rv_makie(bparams; K1, K2, γ=0.0, rv_data1=nothing, rv_data2=nothing, ...) -> (fig, ax)

Radial-velocity curves against orbital phase, with optional data overlaid. `rv_data*` is
`[time value]` or `[time value error]`; a third column switches scatter to error bars.
"""
function plot_rv_makie(bparams; K1::Real, K2::Real, γ::Real = 0.0,
                       rv_data1 = nothing, rv_data2 = nothing,
                       figtitle = "Radial Velocities", figsize = (900, 540))
    ph = collect(range(0, 1, length = 500))
    rv1, rv2 = binary_RV(bparams, bparams.T0 .+ ph .* bparams.P; K1 = K1, K2 = K2, γ = γ)

    fig = Makie.Figure(size = figsize, backgroundcolor = :white)
    ax  = _axis(fig[1, 1]; title = figtitle, xlabel = "Orbital Phase",
                ylabel = "Radial Velocity (km/s)")
    # The systemic velocity, as a reference: the two curves cross ON it, and without it the
    # crossings look like they should be at zero even when γ is not.
    Makie.hlines!(ax, [γ]; color = (:grey60, 0.9), linestyle = :dash, linewidth = 1)
    Makie.lines!(ax, ph, rv1; color = :blue, linewidth = 2, label = "Primary model")
    Makie.lines!(ax, ph, rv2; color = :red,  linewidth = 2, label = "Secondary model")
    for (d, col, mk, lab) in ((rv_data1, :blue, :circle, "Primary data"),
                              (rv_data2, :red,  :rect,   "Secondary data"))
        d === nothing && continue
        φ = mod.(d[:, 1] .- bparams.T0, bparams.P) ./ bparams.P
        if size(d, 2) >= 3
            Makie.errorbars!(ax, φ, d[:, 2], d[:, 3]; color = col)
        end
        Makie.scatter!(ax, φ, d[:, 2]; color = col, marker = mk, markersize = 9, label = lab)
    end
    Makie.axislegend(ax; position = :rt, labelsize = PLOT_LEGENDSIZE)
    return fig, ax
end

# ---------------------------------------------------------------------------------------
# Mollweide
# ---------------------------------------------------------------------------------------
"""
    mollweide_xy(λ, φ) -> (x, y)

Mollweide (equal-area) projection of longitude `λ` and latitude `φ`, both in radians.

Makie has no projected axis and GeoMakie is not a dependency, so this is done by hand. The
projection needs an auxiliary angle θ solving `2θ + sin 2θ = π sin φ`, which has no closed
form; Newton converges in a handful of steps everywhere except the poles, where `θ = ±π/2`
exactly and the derivative vanishes — hence the explicit branch.
"""
function mollweide_xy(λ::Real, φ::Real)
    abs(φ) ≥ π/2 - 1e-12 && return (0.0, sign(φ) * sqrt(2))
    θ = float(φ)
    πsφ = π * sin(φ)
    for _ in 1:20
        den = 2 + 2cos(2θ)
        abs(den) < 1e-12 && break
        δ = (2θ + sin(2θ) - πsφ) / den
        θ -= δ
        abs(δ) < 1e-12 && break
    end
    return (2sqrt(2)/π * λ * cos(θ), sqrt(2) * sin(θ))
end

"""
    mollweide_grid(tmap, star=nothing; nlon=720, nlat=360, visible_pixels=[], mask_unobserved=true)

Resample a surface map onto a regular (lon, lat) grid and project it, returning `(X, Y, values)`
— matrices of projected coordinates plus the sampled map, with unobserved pixels set to `NaN`.

Handles both tessellations. HEALPix goes through `ang2pix_nest`, exactly as the matplotlib
layer does. The long-lat grid is binned directly here rather than through `longlat_ang2pix`,
which ignores its `THETA`/`PHI` arguments entirely and tiles purely by index — so it silently
requires `nlat` and `nlon` to divide `ntheta` and `nphi`, and leaves the caller to fix the
orientation afterwards (the matplotlib layer's `circshift` by half the width is that fix).
Binning from the angles themselves has neither constraint and puts tessel column 1 at
longitude 0, which is where that `circshift` was putting it.

720x360 is 259k quads by default — lower than matplotlib's 2000x1000 because this feeds a mesh
rather than a `pcolormesh` raster, and is already past visibly smooth.
"""
function mollweide_grid(tmap, star = nothing; nlon::Int = 720, nlat::Int = 360,
                        visible_pixels = Int[], mask_unobserved::Bool = true)
    lon = collect(range(-π + 1e-6, π - 1e-6, length = nlon))
    lat = collect(range(-π/2 + 1e-6, π/2 - 1e-6, length = nlat))
    if star === nothing || star.tessellation_type == 0
        nside = npix2nside(length(tmap))
        THETA = [π/2 - φ for φ in lat, _ in lon]        # colatitude
        PHI   = [λ for _ in lat, λ in lon]
        pix   = reshape(ang2pix_nest(nside, vec(THETA), vec(PHI)), size(THETA))
    else
        nt, np = _latlong_dims(star)
        # Tessel numbering runs φ fastest inside θ (src/tessellation_latlong.jl:37-45), with
        # θ-row 1 at the North pole.
        pix = [ (clamp(1 + floor(Int, (π/2 - φ) / π * nt), 1, nt) - 1) * np +
                 clamp(1 + floor(Int, mod(λ, 2π) / (2π) * np), 1, np)
                for φ in lat, λ in lon ]
    end
    vals  = Float64.(tmap[pix])
    if mask_unobserved && !isempty(visible_pixels)
        vset = Set(visible_pixels)
        @inbounds for i in eachindex(pix)
            pix[i] in vset || (vals[i] = NaN)
        end
    end
    X = similar(vals); Y = similar(vals)
    @inbounds for j in axes(vals, 2), i in axes(vals, 1)
        X[i, j], Y[i, j] = mollweide_xy(lon[j], lat[i])
    end
    return X, Y, vals
end

"""
    plot_mollweide_makie(tmap, star=nothing; ...) -> (fig, ax)

Whole-surface view in a Mollweide equal-area projection.

Drawn as a `mesh!` over the projected grid rather than a heatmap: the projected cells are not
axis-aligned, so no rectangular primitive fits them.

`visible_pixels` is EMPTY by default and nothing is masked, matching `plot_mollweide`. Passing
`star` selects the tessellation, it does not select a coverage mask — `index_quads_visible` is
one epoch's front hemisphere, so masking with it would grey out half a map that the other
epochs constrain perfectly well. To show coverage, pass the union of the epochs' visible sets
explicitly; those pixels are then painted `bad_color` rather than dropped, so a gap reads as a
gap instead of as absence.
"""
function plot_mollweide_makie(tmap, star = nothing; visible_pixels = Int[],
                              vmin = nothing, vmax = nothing, colormap = "gist_heat",
                              figtitle = "Mollweide", mask_unobserved = true, incl = 90.0,
                              bad_color = Makie.RGBAf(0.83, 0.83, 0.83, 1),
                              lon_color = Makie.RGBAf(1, 1, 1, 0.85),
                              lat_color = Makie.RGBAf(0, 0, 0, 0.55),
                              labels = true, nlon = 720, nlat = 360, colorbar = true,
                              figsize = (960, 620))
    X, Y, vals = mollweide_grid(tmap, star; nlon, nlat, visible_pixels, mask_unobserved)
    finite = filter(isfinite, vals)
    pmin, pmax = _map_range(isempty(finite) ? Float64.(tmap) : finite; vmin, vmax)
    cmap = _padded_cmap(colormap)

    ny, nx = size(vals)
    pts  = [Makie.Point2f(X[i, j], Y[i, j]) for j in 1:nx for i in 1:ny]
    cols = Vector{Makie.RGBAf}(undef, length(pts))
    lo, hi = pmin, pmax; s = hi > lo ? 1/(hi - lo) : 0.0
    @inbounds for j in 1:nx, i in 1:ny
        v = vals[i, j]
        cols[(j-1)*ny + i] = isfinite(v) ?
            Makie.RGBAf(Makie.to_color(cmap[clamp((v - lo)*s, 0, 1)])) : bad_color
    end
    faces = Makie.TriangleFace{Int}[]
    sizehint!(faces, 2*(nx-1)*(ny-1))
    @inbounds for j in 1:nx-1, i in 1:ny-1
        a = (j-1)*ny + i; b = a + 1; c = j*ny + i; d = c + 1
        push!(faces, Makie.TriangleFace(a, b, d)); push!(faces, Makie.TriangleFace(a, d, c))
    end

    fig = Makie.Figure(size = figsize, backgroundcolor = :white)
    ax  = _axis(fig[1, 1]; title = figtitle, aspect = Makie.DataAspect())
    Makie.hidedecorations!(ax); Makie.hidespines!(ax)
    Makie.mesh!(ax, Makie.GeometryBasics.Mesh(pts, faces); color = cols,
                shading = Makie.NoShading)

    # Graticule, bounding ellipse and labels are all traced through the SAME projection as the
    # data, so they cannot drift from it — which is the one thing a hand-rolled projection can
    # get wrong invisibly. Meridians every 20 deg and parallels every 20 deg, matching
    # `set_longitude_grid(20)` / `set_latitude_grid(20)`.
    for λ0 in range(-π, π, length = 19)
        Makie.lines!(ax, [Makie.Point2f(mollweide_xy(λ0, φ)...)
                          for φ in range(-π/2, π/2, length = 180)];
                     color = lon_color, linewidth = 0.6)
    end
    for φ0 in range(-π/2, π/2, length = 10)
        Makie.lines!(ax, [Makie.Point2f(mollweide_xy(λ, φ0)...)
                          for λ in range(-π, π, length = 360)];
                     color = lat_color, linewidth = 0.6)
    end
    Makie.lines!(ax, [Makie.Point2f(mollweide_xy(λ, φ)...)
                      for (λ, φ) in vcat([(-π + 1e-9, φ) for φ in range(-π/2, π/2, length=180)],
                                          [( π - 1e-9, φ) for φ in range(π/2, -π/2, length=180)])];
                 color = (:black, 0.6), linewidth = 1.0)
    if labels
        for λ0 in range(-160, 160, step = 20)
            x, y = mollweide_xy(deg2rad(λ0), 0.0)
            Makie.text!(ax, x, y; text = "$(Int(λ0))°", color = lon_color,
                        fontsize = PLOT_TICKLABELSIZE - 2, align = (:center, :center))
        end
        for φ0 in (-80, -60, -40, -20, 20, 40, 60, 80)
            # Just outside the bounding ellipse at that latitude, not at a fixed longitude:
            # the ellipse narrows towards the poles, so a fixed longitude walks the high
            # labels INTO the map.
            xe, y = mollweide_xy(-π, deg2rad(φ0))
            Makie.text!(ax, xe - 0.06, y; text = "$(φ0)°", color = :black,
                        fontsize = PLOT_TICKLABELSIZE, align = (:right, :center))
        end
    end
    # The sub-observer inclination, drawn as the parallel at -incl the way `axhline` does.
    if incl != 90.0
        Makie.lines!(ax, [Makie.Point2f(mollweide_xy(λ, -deg2rad(incl))...)
                          for λ in range(-π, π, length = 360)];
                     color = :black, linestyle = :dashdot, linewidth = 1.2)
    end
    colorbar && Makie.Colorbar(fig[2, 1]; colormap = cmap, colorrange = (pmin, pmax),
                               vertical = false, label = "Temperature (K)",
                               labelsize = PLOT_LABELSIZE, ticklabelsize = PLOT_TICKLABELSIZE,
                               width = Makie.Relative(0.6), ticks = _bar_ticks(pmin, pmax))
    return fig, ax
end
