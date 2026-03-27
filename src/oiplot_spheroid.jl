using PyPlot,PyCall, LaTeXStrings, Statistics

# Suppress matplotlib warning about tight_layout + equal-aspect axes
# (deferred to first plot call via set_oiplot_defaults)

const _oiplot_initialized = Ref(false)

function set_oiplot_defaults()
    if !_oiplot_initialized[]
        pyimport("warnings").filterwarnings("ignore", message=".*tight_layout.*")
        _oiplot_initialized[] = true
    end
    rc = PyDict(pyimport("matplotlib")."rcParams")
    rc["font.family"] = "serif"
    rc["font.size"] = 14
    rc["xtick.major.size"] = 6
    rc["ytick.major.size"] = 6
    rc["xtick.minor.size"] = 6
    rc["ytick.minor.size"] = 6
    rc["xtick.major.width"] = 1
    rc["ytick.major.width"] = 1
    rc["xtick.minor.width"] = 1
    rc["ytick.minor.width"] = 1
    rc["lines.markeredgewidth"] = 1
    rc["legend.numpoints"] = 1
    rc["legend.handletextpad"] = 0.3
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
    draw_compass(ax, axis_max; size_frac=0.12, fontsize=12, color="black")

Draw E/N compass arrows in the lower-right corner of a 2D sky-plane plot.
East points left (astronomical convention when x-axis is inverted).
"""
function draw_compass(ax, axis_max; size_frac=0.12, fontsize=12, color="black")
    s = size_frac * axis_max
    margin = 0.20 * axis_max
    cx = -(axis_max - margin)
    cy = -(axis_max - margin)
    gap = 0.15 * s  # small offset so E and N labels don't crowd the corner
    ax.annotate("N", xy=(cx - gap, cy + s), xytext=(cx - gap, cy),
        arrowprops=Dict("arrowstyle" => "-|>", "lw" => 1.5, "color" => color),
        fontsize=fontsize, color=color, ha="center", va="bottom", fontweight="bold", zorder=8)
    ax.annotate("E", xy=(cx + s, cy - gap), xytext=(cx, cy - gap),
        arrowprops=Dict("arrowstyle" => "-|>", "lw" => 1.5, "color" => color),
        fontsize=fontsize, color=color, ha="right", va="center", fontweight="bold", zorder=8)
end

"""
    draw_rotation_axis(ax, star; inclination=NaN, position_angle=NaN,
        arrow_frac=0.3, color="black", linewidth=1.5, offset_west=0.0, offset_north=0.0)

Draw the projected stellar rotation axis on a 2D sky-plane plot.
Shows the pole-to-pole line extended beyond the limb, with an arrow at the north pole.
Parts behind the star (z < 0) are drawn dashed.

When `inclination` (degrees from LOS) and `position_angle` (degrees, N through E) are given,
the axis is computed analytically. Otherwise it is estimated from the tessellation vertices.
"""
function draw_rotation_axis(ax, star; arrow_frac=0.3, color="black", linewidth=1.5,
    offset_west=0.0, offset_north=0.0, inclination=NaN, position_angle=NaN)
    if !isnan(inclination) && !isnan(position_angle)
        inc_rad = inclination * π / 180
        PA_rad  = position_angle * π / 180
        spin = [-sin(PA_rad) * sin(inc_rad), cos(PA_rad) * sin(inc_rad), cos(inc_rad)]
        R = maximum(sqrt.(star.vertices_xyz[:,:,1].^2 .+ star.vertices_xyz[:,:,2].^2 .+ star.vertices_xyz[:,:,3].^2))
        north = R .* spin
        south = -R .* spin
    else
        np = star.npix
        north = vec(mean(star.vertices_xyz[1:min(4,np), 5, :], dims=1))
        south = vec(mean(star.vertices_xyz[max(1,np-3):np, 5, :], dims=1))
    end
    delta = north .- south
    north_tip = north .+ arrow_frac .* delta
    south_tip = south .- arrow_frac .* delta
    points = hcat(south_tip, south, north, north_tip)  # 3×4
    for j in 1:3
        p1 = points[:, j]
        p2 = points[:, j+1]
        z_avg = (p1[3] + p2[3]) / 2
        ls = z_avg > 0 ? "-" : "--"
        zord = z_avg > 0 ? 6 : 1
        ax.plot([-(p1[1]+offset_west), -(p2[1]+offset_west)],
                [p1[2]+offset_north, p2[2]+offset_north],
                ls, color=color, linewidth=linewidth, zorder=zord)
    end
    zord_tip = north_tip[3] > 0 ? 6 : 1
    ax.annotate("", xy=(-(north_tip[1]+offset_west), north_tip[2]+offset_north),
        xytext=(-(north[1]+offset_west), north[2]+offset_north),
        arrowprops=Dict("arrowstyle" => "-|>", "color" => color, "lw" => linewidth),
        zorder=zord_tip)
end

"""
    draw_rotation_arrow(ax, star; pole="N", radius_frac=0.15, offset_frac=0.05,
        color="black", linewidth=1.5, inclination=NaN, position_angle=NaN)

Draw a curved arrow around the rotation axis showing the sense of rotation.
A 270° ellipse centered at the pole, solid in front (z>0), dashed behind.
Prograde rotation = counter-clockwise around north pole.
When `inclination`/`position_angle` are given, the axis is computed analytically.
"""
function draw_rotation_arrow(ax, star; pole="N", radius_frac=0.15, offset_frac=0.05,
    color="black", linewidth=1.5, npoints=200, offset_west=0.0, offset_north=0.0,
    inclination=NaN, position_angle=NaN)
    if !isnan(inclination) && !isnan(position_angle)
        inc_rad = inclination * π / 180
        PA_rad  = position_angle * π / 180
        spin = [-sin(PA_rad) * sin(inc_rad), cos(PA_rad) * sin(inc_rad), cos(inc_rad)]
        R = maximum(sqrt.(star.vertices_xyz[:,:,1].^2 .+ star.vertices_xyz[:,:,2].^2 .+ star.vertices_xyz[:,:,3].^2))
        north = R .* spin
        south = -R .* spin
    else
        north = star.vertices_xyz[1, 1, 1:3]
        south = star.vertices_xyz[end, 3, 1:3]
    end
    axis = north .- south
    axis_len = norm(axis)
    axis_hat = axis ./ axis_len
    center = pole == "N" ? north .+ offset_frac .* axis : south .- offset_frac .* axis
    r = radius_frac * axis_len
    ref = abs(axis_hat[3]) < 0.9 ? [0.0, 0.0, 1.0] : [1.0, 0.0, 0.0]
    e1 = cross(axis_hat, ref); e1 ./= norm(e1)
    e2 = cross(axis_hat, e1)
    θ = collect(range(0, 3π/2, length=npoints))
    pts = hcat([center .+ r .* (cos(t) .* e1 .+ sin(t) .* e2) for t in θ]...)
    x2d = -(pts[1, :] .+ offset_west); y2d = pts[2, :] .+ offset_north; z = pts[3, :]
    front = z .> 0
    xf = copy(x2d); xf[.!front] .= NaN; yf = copy(y2d); yf[.!front] .= NaN
    ax.plot(xf, yf, "-", color=color, linewidth=linewidth, zorder=7)
    xb = copy(x2d); xb[front] .= NaN; yb = copy(y2d); yb[front] .= NaN
    ax.plot(xb, yb, "--", color=color, linewidth=linewidth, zorder=1)
    zord_tip = z[end] > 0 ? 7 : 1
    ax.annotate("", xy=(x2d[end], y2d[end]), xytext=(x2d[end-3], y2d[end-3]),
        arrowprops=Dict("arrowstyle" => "-|>", "color" => color, "lw" => linewidth),
        zorder=zord_tip)
end

"""
    draw_graticules(ax, star; nlat=5, nlon=8, inclination=NaN, position_angle=NaN, ...)

Draw latitude/longitude graticule lines on the star surface using parametric curves.
Generates smooth curves in the body frame, rotates them with the same Euler rotation
as the star (rot_vertex), and z-clips to the front hemisphere.
Renders via matplotlib PolyCollection for efficiency.

Semi-axes are extracted automatically from the tessellation's body-frame radii.
"""
function draw_graticules(ax, star; nlat=5, nlon=8, color="black", linewidth=0.8, alpha=0.5,
    offset_west=0.0, offset_north=0.0, inclination=NaN, position_angle=NaN,
    rotation_angle=0.0, npoints=200)
    collections = pyimport("matplotlib.collections")

    # Extract effective semi-axes from body-frame tessellation
    all_colat = star.vertices_spherical[:, 1:4, 2]   # colatitude (body frame)
    all_r     = star.vertices_spherical[:, 1:4, 1]    # radius (body frame)

    pole_mask = (all_colat .< 0.3) .| (all_colat .> π - 0.3)
    r_pole = any(pole_mask) ? mean(all_r[pole_mask]) : mean(all_r)
    eq_mask = abs.(all_colat .- π/2) .< 0.3
    r_eq = any(eq_mask) ? mean(all_r[eq_mask]) : mean(all_r)

    # Build rotation matrix (same convention as rotate_star)
    inc_rad = isnan(inclination) ? 0.0 : inclination * π / 180
    PA_rad  = isnan(position_angle) ? 0.0 : position_angle * π / 180
    ψ_rad   = rotation_angle * π / 180
    R = rot_vertex(ψ_rad, inc_rad, PA_rad)

    graticule_lines = Vector{Matrix{Float64}}()

    # Latitude circles (constant colatitude)
    θ_targets = collect(range(π/(nlat+1), stop=π*nlat/(nlat+1), length=nlat))
    ϕ_range = collect(range(-π, stop=π, length=npoints))
    for θ0 in θ_targets
        body_pts = hcat(r_eq .* sin(θ0) .* cos.(ϕ_range),
                        r_eq .* sin(θ0) .* sin.(ϕ_range),
                        fill(r_pole * cos(θ0), npoints))
        sky_pts = body_pts * R
        append!(graticule_lines, _visible_segments(sky_pts, offset_west, offset_north))
    end

    # Longitude lines (constant azimuth)
    ϕ_targets = collect(range(0, stop=2π*(1 - 1/nlon), length=nlon))
    θ_range = collect(range(0, stop=π, length=npoints))
    for ϕ0 in ϕ_targets
        body_pts = hcat(r_eq .* sin.(θ_range) .* cos(ϕ0),
                        r_eq .* sin.(θ_range) .* sin(ϕ0),
                        r_pole .* cos.(θ_range))
        sky_pts = body_pts * R
        append!(graticule_lines, _visible_segments(sky_pts, offset_west, offset_north))
    end

    if !isempty(graticule_lines)
        ax.add_collection(collections.PolyCollection(graticule_lines, closed=false,
            ec=color, fc="none", linewidths=[linewidth], alpha=alpha, zorder=5))
    end
end

"""Extract contiguous front-hemisphere (z>0) segments as Nx2 matrices for PolyCollection."""
function _visible_segments(sky_pts, offset_west, offset_north)
    segments = Vector{Matrix{Float64}}()
    seg_start = 0
    npts = size(sky_pts, 1)
    for k in 1:npts
        if sky_pts[k, 3] > 0
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
  PyPlot.draw()
  return 1
end

function plot3d(star_temperature_map,star) # this plots the temperature map
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

  for i=1:star.npix
      verts = (collect(zip(corners_xyz[i, :, 1], corners_xyz[i, :, 2], corners_xyz[i, :, 3])),);
      color = get_cmap("gist_heat")((star_temperature_map[i] - minimum(star_temperature_map)) / (maximum(star_temperature_map) - minimum(star_temperature_map)));
      ax.add_collection3d(Poly3DCollection(verts, edgecolor="none", facecolor=color));
  end
  ax.set_aspect("equal")
  PyPlot.draw()
end


function plot2d(tmap, star; intensity = false, figtitle ="", plotmesh=false, pad = 0.5,
    colormap="gist_heat", xlim=Float64[], ylim=Float64[], background="white", flipx=false,
    compass=true, rotation_axis=false, rotation_arrow=false, graticules=false,
    inclination=NaN, position_angle=NaN)
  # Plot temperature map onto the projected 2D image plane (= observer view)
  # Convention: East left, North up (astronomical standard)
  set_oiplot_defaults()
  patches = pyimport("matplotlib.patches")
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
  projmap = tmap[star.index_quads_visible];
  if intensity == true
    projmap .*=star.ldmap[star.index_quads_visible]
  end
  # Normalize with padded floor so minimum maps to cfloor (~dark red, not black)
  pmin = minimum(projmap); pmax = maximum(projmap)
  prange = pmax - pmin
  if prange < 1.0; prange = max(abs(pmax) * 0.01, 1.0); end
  cfloor = 0.08
  vmin_padded = pmin - cfloor / (1.0 - cfloor) * prange
  norm_plot = matplotlib.colors.Normalize(vmin=vmin_padded, vmax=pmax)
  colours = get_cmap(colormap).(norm_plot.(projmap))

  visible = star.index_quads_visible
  for i=1:star.nquads_visible
  idx = visible[i]
  p = patches.Polygon(hcat(-star.proj_west[idx,:],star.proj_north[idx,:]),closed=true,edgecolor= (plotmesh == true) ? "lightgrey" : colours[i],facecolor=colours[i],fill=true,rasterized=false, zorder=2)
  ax.add_patch(p);
  end
  xlabel(L"x $\leftarrow$ E (mas)", fontsize=20)
  ylabel(L"y $\rightarrow$ N (mas)", fontsize=20)
  # Decorations: graticules (z=5) < pole line (z=6) < spin arrow (z=7) < compass (z=8)
  if graticules; draw_graticules(ax, star, inclination=inclination, position_angle=position_angle); end
  if rotation_axis; draw_rotation_axis(ax, star, inclination=inclination, position_angle=position_angle); end
  if rotation_arrow; draw_rotation_arrow(ax, star, inclination=inclination, position_angle=position_angle); end
  if compass; draw_compass(ax, axis_max); end
  cmap=ColorMap(colormap)
  divider = axdiv.make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.07)
  cb=colorbar(matplotlib.cm.ScalarMappable(norm=norm_plot,cmap=cmap), cax=cax)
  if pmin < pmax; cb.ax.set_ylim(pmin, pmax); end
  return fig, ax
  end

"""
    plot2d_binary(tmap1, tmap2, star1, star2, bparams, tepoch; ...)

Plot a binary system on the 2D sky plane with correct occlusion and orbital offset.
`tepoch` is the observation time in JD; the secondary's offset relative to the primary
is computed from the orbital elements in `bparams`.
The farther star (larger z = receding) is drawn behind the nearer one.
"""
function plot2d_binary(tmap1, tmap2, star1, star2, bparams, tepoch;
    intensity=false, plotmesh=false, colormap="gist_heat", pad=1.0, background="white",
    compass=true, rotation_axis=false, rotation_arrow=false, graticules=false, figtitle="",
    inclination1=NaN, position_angle1=NaN, inclination2=NaN, position_angle2=NaN)
  set_oiplot_defaults()
  patches = pyimport("matplotlib.patches")
  axdiv = pyimport("mpl_toolkits.axes_grid1.axes_divider")
  facecolor = background == "black" ? "Black" : "White"
  fig = figure("Binary epoch", figsize=(10,10), facecolor="White")
  clf()
  ax = gca()
  title(figtitle)
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
  axis_max = max(r1, r2 + abs(east_offset), r2 + abs(north_offset),
                 abs(east_offset) + r2, abs(north_offset) + r2) + pad
  ax.set_xlim([axis_max, -axis_max])
  ax.set_ylim([-axis_max, axis_max])
  # Shared color normalization across both stars
  # Pad the bottom of the range so the coolest temperature maps to ~0.15
  # instead of 0.0 (pure black in gist_heat), ensuring both stars are visible.
  projmap1 = tmap1[star1.index_quads_visible]
  projmap2 = tmap2[star2.index_quads_visible]
  if intensity
    projmap1 = projmap1 .* star1.ldmap[star1.index_quads_visible]
    projmap2 = projmap2 .* star2.ldmap[star2.index_quads_visible]
  end
  tmin = min(minimum(projmap1), minimum(projmap2))
  tmax = max(maximum(projmap1), maximum(projmap2))
  trange = tmax - tmin
  if trange < 1.0; trange = max(tmax * 0.01, 1.0); end
  cfloor = 0.15  # minimum colormap fraction (avoids black pixels on dark background)
  vmin_padded = tmin - cfloor / (1.0 - cfloor) * trange  # tmin maps to cfloor in colormap
  norm_plot = matplotlib.colors.Normalize(vmin=vmin_padded, vmax=tmax)
  colours1 = get_cmap(colormap).(norm_plot.(projmap1))
  colours2 = get_cmap(colormap).(norm_plot.(projmap2))
  # Determine z-ordering: farther star (larger z = receding) drawn first (behind)
  _, _, z1, _, _, z2 = binary_orbit_abs(bparams, tepoch)
  zord1 = z1 > z2 ? 2 : 3
  zord2 = z1 > z2 ? 3 : 2
  # Star 1 at origin
  vis1 = star1.index_quads_visible
  for i=1:star1.nquads_visible
    idx = vis1[i]
    ec = plotmesh ? "lightgrey" : colours1[i]
    p = patches.Polygon(hcat(-star1.proj_west[idx,:], star1.proj_north[idx,:]),
      closed=true, edgecolor=ec, facecolor=colours1[i], fill=true, rasterized=false, zorder=zord1)
    ax.add_patch(p)
  end
  # Star 2 shifted by orbital offset
  vis2 = star2.index_quads_visible
  for i=1:star2.nquads_visible
    idx = vis2[i]
    ec = plotmesh ? "lightgrey" : colours2[i]
    p = patches.Polygon(hcat(-(star2.proj_west[idx,:] .+ offset_west), star2.proj_north[idx,:] .+ offset_north),
      closed=true, edgecolor=ec, facecolor=colours2[i], fill=true, rasterized=false, zorder=zord2)
    ax.add_patch(p)
  end
  xlabel(L"x $\leftarrow$ E (mas)", fontsize=20)
  ylabel(L"y $\rightarrow$ N (mas)", fontsize=20)
  # Decorations: graticules (z=5) < pole line (z=6) < spin arrow (z=7) < compass (z=8)
  if graticules
    draw_graticules(ax, star1, inclination=inclination1, position_angle=position_angle1)
    draw_graticules(ax, star2, offset_west=offset_west, offset_north=offset_north,
        inclination=inclination2, position_angle=position_angle2)
  end
  if rotation_axis
    draw_rotation_axis(ax, star1, inclination=inclination1, position_angle=position_angle1)
    draw_rotation_axis(ax, star2, offset_west=offset_west, offset_north=offset_north,
        inclination=inclination2, position_angle=position_angle2)
  end
  if rotation_arrow
    draw_rotation_arrow(ax, star1, inclination=inclination1, position_angle=position_angle1)
    draw_rotation_arrow(ax, star2, offset_west=offset_west, offset_north=offset_north,
        inclination=inclination2, position_angle=position_angle2)
  end
  if compass; draw_compass(ax, axis_max); end
  # Colorbar — use the padded norm so colors match the patches
  cmap = ColorMap(colormap)
  divider = axdiv.make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.07)
  cb = colorbar(matplotlib.cm.ScalarMappable(norm=norm_plot, cmap=cmap), cax=cax)
  if tmin < tmax; cb.ax.set_ylim(tmin, tmax); end  # clip colorbar ticks to actual temperature range
  return fig, ax
end

function plot2d_wireframe(star; compass=true, rotation_axis=false)
  # Wireframe view of the tessellation projected onto the sky plane
  patches = pyimport("matplotlib.patches")
  fig = figure("Epoch tmap",figsize=(10,10),facecolor="White")
  ax = fig.add_axes([0.1,0.1,0.85,0.85]);
  xlabel(L"x $\leftarrow$ E (mas)", fontsize=14);
  ylabel(L"y $\rightarrow$ N (mas)", fontsize=14);
  axis_max = maximum(sqrt.(star.vertices_xyz[:,:,1].^2 .+ star.vertices_xyz[:,:,2].^2 .+ star.vertices_xyz[:,:,3].^2))*1.5;
  ax.set_xlim([axis_max,-axis_max]);
  ax.set_ylim([-axis_max,axis_max]);
  visible = star.index_quads_visible
  for i=1:star.nquads_visible
  idx = visible[i]
  p = patches.Polygon(hcat(-star.proj_west[idx,:],star.proj_north[idx,:]),closed=true,edgecolor="black", facecolor="white",rasterized=false)
  ax.add_patch(p);
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
  PyPlot.draw();
  return;
end

function plot2d_allepochs(tmap, star; plotmesh=false, tepochs = [], colormap="gist_heat", arr_box=23, compass=true)
    nepochs = length(star)
    patches = pyimport("matplotlib.patches")
    fig = figure("Temperature map -- All epochs",figsize=(15,10),facecolor="White")
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
      ax = subplot(arr_box*10+t)
      if tepochs !=[]
        title("Epoch $t $(tepochs[t])",fontweight="bold")
      end
      axis_max = maximum(sqrt.(star[t].vertices_xyz[:,:,1].^2 .+ star[t].vertices_xyz[:,:,2].^2 .+ star[t].vertices_xyz[:,:,3].^2))*1.5;
      ax.set_xlim([axis_max,-axis_max])
      ax.set_ylim([-axis_max,axis_max])
      projmap = (tmap[star[t].index_quads_visible].-minT)./(maxT-minT);
      visible = star[t].index_quads_visible
      for i=1:star[t].nquads_visible
        idx = visible[i]
        p = patches.Polygon(hcat(-star[t].proj_west[idx, :],star[t].proj_north[idx, :]),
        closed=true,edgecolor=meshcolor,facecolor=get_cmap(colormap)(projmap[i]),fill="true",rasterized=false)
        ax.add_patch(p);
      end
      ax.plot();
      ax.tick_params(axis="both", which="both", labelsize=15, width=2, length=10);
      set_tick_spacing(ax, axis_max)
      ax.tick_params(axis="both", which="minor", width=1, length=5);
      for spine in ["top", "bottom", "left", "right"]
        ax.spines[spine].set_linewidth(1);
      end
      ax.set_aspect("equal");
      if compass && t == 1; draw_compass(ax, axis_max, color="black"); end
    end
    fig.canvas.draw();
    fig.text(0.5, 0.04, L"x $\leftarrow$ E (mas)", ha="center", va="center", fontweight="bold", fontsize=15);
    fig.text(0.06, 0.5, L"y $\rightarrow$ N (mas)", ha="center", va="center", rotation="vertical", fontweight="bold", fontsize=15);
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

function plot_mollweide(tmap, star; kwargs...)
  if star.tessellation_type == 0
    mollplot_temperature_healpix(tmap,kwargs...)
  else
    # Longitude and latitude need to be provided in kwargs
    #... or could we do otherwise and recompute from theta/phi?
    mollplot_temperature_longlat(tmap,kwargs...)
  end
  return
end

function mollplot_temperature_healpix(tmap; visible_pixels = [], vmin = -Inf, vmax = Inf, incl=90.0, colormap="gist_heat", figtitle="Mollweide")
  xsize = 2000
  ysize = div(xsize,2)
  theta = collect(range(pi, stop=0.0, length=ysize))
  phi   = collect(range(-pi, stop=pi, length=xsize))
  longitude = collect(range(-180, stop=180, length=xsize))/180*pi
  latitude = collect(range(-90, stop=90, length=ysize))/180*pi
  # project the map to a rectangular matrix xsize x ysize
  nside = npix2nside(length(tmap))
  PHI = [i for j in theta, i in phi]
  THETA = [j for j in theta, i in phi]
  grid_pix = reshape(ang2pix_nest(nside, vec(THETA), vec(PHI)), size(PHI))
  grid_map = tmap[grid_pix]
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
  moll = pcolormesh(longitude, latitude, grid_map, vmin=vmin, vmax=vmax, rasterized=true, cmap=colormap)
  # graticule
  ax.set_longitude_grid(20);
  ax.set_latitude_grid(20);
  ax.set_longitude_grid_ends(90);
  spacing = 0.04
  subplots_adjust(bottom=spacing, top=1-spacing, left=spacing, right=1-spacing);
  grid(true)
  if incl != 90.0
    ax.axhline(-incl * pi/180, c=:black, ls="-.")
  end
  ticks = collect(range(vmin, stop=vmax, length=7));
  cb = colorbar(moll, orientation="horizontal", shrink=.6, pad=0.05, ticks=ticks)
  cb.ax.xaxis.labelpad=5
  cb.ax.xaxis.set_label_text("Temperature (K)")
  # workaround for issue with viewers, see colorbar docstring
  cb.solids.set_edgecolor("face")
  ax.tick_params(axis="x", labelsize=12)
  ax.tick_params(axis="y", labelsize=12)
  return
end

function mollplot_temperature_longlat(tmap, ntheta, nphi; visible_pixels = [], vmin = -Inf, vmax = Inf, colormap="gist_heat", figtitle="Mollweide", incl=90.0)
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
    ax.axhline(-incl * pi/180, c=:black, ls="-.")
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
