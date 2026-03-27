using PyPlot,PyCall, LaTeXStrings, Statistics

function set_oiplot_defaults()
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
    draw_compass(ax, axis_max; size_frac=0.12, fontsize=12, color="white")

Draw E/N compass arrows in the lower-right corner of a 2D sky-plane plot.
East points left (astronomical convention when x-axis is inverted).
"""
function draw_compass(ax, axis_max; size_frac=0.12, fontsize=12, color="white")
    s = size_frac * axis_max
    margin = 0.20 * axis_max
    # Lower-right corner (x-axis runs [+max, -max] = East left)
    cx = -(axis_max - margin)
    cy = -(axis_max - margin)
    # North arrow (up)
    ax.annotate("N", xy=(cx, cy + s), xytext=(cx, cy),
        arrowprops=Dict("arrowstyle" => "-|>", "lw" => 1.5, "color" => color),
        fontsize=fontsize, color=color, ha="center", va="bottom", fontweight="bold")
    # East arrow (toward positive x = left on sky, since xlim=[max, -max])
    ax.annotate("E", xy=(cx + s, cy), xytext=(cx, cy),
        arrowprops=Dict("arrowstyle" => "-|>", "lw" => 1.5, "color" => color),
        fontsize=fontsize, color=color, ha="right", va="center", fontweight="bold")
end

"""
    draw_rotation_axis(ax, star; arrow_frac=0.3, color="cyan", linewidth=1.5)

Draw the projected stellar rotation axis on a 2D sky-plane plot.
Shows the pole-to-pole line extended beyond the limb, with an arrow at the north pole.
Parts behind the star (z < 0) are drawn dashed.
"""
function draw_rotation_axis(ax, star; arrow_frac=0.3, color="cyan", linewidth=1.5)
    north = star.vertices_xyz[1, 1, 1:3]
    south = star.vertices_xyz[end, 3, 1:3]
    delta = north .- south
    north_tip = north .+ arrow_frac .* delta
    south_tip = south .- arrow_frac .* delta
    # Project to 2D sky: East = -xyz[1] (proj_west), North = xyz[2] (proj_north)
    points = hcat(south_tip, south, north, north_tip)  # 3×4
    for j in 1:3
        p1 = points[:, j]
        p2 = points[:, j+1]
        z_avg = (p1[3] + p2[3]) / 2
        ls = z_avg > 0 ? "-" : "--"
        zord = z_avg > 0 ? 11 : 1
        ax.plot([-p1[1], -p2[1]], [p1[2], p2[2]], ls, color=color, linewidth=linewidth, zorder=zord)
    end
    # Arrow head at north tip
    zord_tip = north_tip[3] > 0 ? 11 : 1
    ax.annotate("", xy=(-north_tip[1], north_tip[2]), xytext=(-north[1], north[2]),
        arrowprops=Dict("arrowstyle" => "-|>", "color" => color, "lw" => linewidth),
        zorder=zord_tip)
end

"""
    draw_rotation_arrow(ax, star; pole="N", radius_frac=0.15, offset_frac=0.4,
        color="cyan", linewidth=1.5)

Draw a curved arrow around the rotation axis showing the sense of rotation.
A 270° ellipse centered along the pole, solid in front (z>0), dashed behind.
Prograde rotation = counter-clockwise around north pole.
"""
function draw_rotation_arrow(ax, star; pole="N", radius_frac=0.15, offset_frac=0.4,
    color="cyan", linewidth=1.5, npoints=200)
    north = star.vertices_xyz[1, 1, 1:3]
    south = star.vertices_xyz[end, 3, 1:3]
    axis = north .- south
    axis_len = norm(axis)
    axis_hat = axis ./ axis_len
    center = pole == "N" ? north .+ offset_frac .* axis : south .- offset_frac .* axis
    r = radius_frac * axis_len
    # Two orthogonal vectors perpendicular to the rotation axis
    ref = abs(axis_hat[3]) < 0.9 ? [0.0, 0.0, 1.0] : [1.0, 0.0, 0.0]
    e1 = cross(axis_hat, ref); e1 ./= norm(e1)
    e2 = cross(axis_hat, e1)
    # 270° arc (prograde = CCW around north pole in 3D)
    θ = collect(range(0, 3π/2, length=npoints))
    pts = hcat([center .+ r .* (cos(t) .* e1 .+ sin(t) .* e2) for t in θ]...)
    x2d = -pts[1, :]; y2d = pts[2, :]; z = pts[3, :]
    # Solid in front, dashed behind, using NaN to break the line
    front = z .> 0
    xf = copy(x2d); xf[.!front] .= NaN; yf = copy(y2d); yf[.!front] .= NaN
    ax.plot(xf, yf, "-", color=color, linewidth=linewidth, zorder=11)
    xb = copy(x2d); xb[front] .= NaN; yb = copy(y2d); yb[front] .= NaN
    ax.plot(xb, yb, "--", color=color, linewidth=linewidth, zorder=1)
    # Arrow tip at end of arc
    zord_tip = z[end] > 0 ? 11 : 1
    ax.annotate("", xy=(x2d[end], y2d[end]), xytext=(x2d[end-3], y2d[end-3]),
        arrowprops=Dict("arrowstyle" => "-|>", "color" => color, "lw" => linewidth),
        zorder=zord_tip)
end

"""
    draw_graticules(ax, star; nlat=5, nlon=8, color="gray", linewidth=0.5, alpha=0.5, front_only=true)

Draw latitude/longitude lines on the star surface by finding contour crossings
in the tessellation. Works for any surface type (sphere, ellipsoid, Roche lobe).
Uses the body-frame spherical coordinates from `vertices_spherical`.
"""
function draw_graticules(ax, star; nlat=5, nlon=8, color="gray", linewidth=0.5, alpha=0.5, front_only=true)
    visible = star.index_quads_visible
    quad_edges = ((1,2), (2,3), (3,4), (4,1))
    function _contour_segments(field_idx, targets; periodic=false)
        for target in targets
            for i in visible
                crossings = Vector{Vector{Float64}}()
                for (j1, j2) in quad_edges
                    v1 = star.vertices_spherical[i, j1, field_idx]
                    v2 = star.vertices_spherical[i, j2, field_idx]
                    if periodic
                        d1 = mod(v1 - target + π, 2π) - π
                        d2 = mod(v2 - target + π, 2π) - π
                    else
                        d1 = v1 - target; d2 = v2 - target
                    end
                    if d1 * d2 < 0 && abs(d1 - d2) < π
                        t = d1 / (d1 - d2)
                        pt = star.vertices_xyz[i, j1, :] .+ t .* (star.vertices_xyz[i, j2, :] .- star.vertices_xyz[i, j1, :])
                        push!(crossings, pt)
                    end
                end
                if length(crossings) == 2
                    p1, p2 = crossings
                    z_avg = (p1[3] + p2[3]) / 2
                    if !front_only || z_avg > 0
                        zord = z_avg > 0 ? 10 : 1
                        ax.plot([-p1[1], -p2[1]], [p1[2], p2[2]], "-",
                            color=color, linewidth=linewidth, alpha=alpha, zorder=zord)
                    end
                end
            end
        end
    end
    # Latitude lines (constant colatitude θ)
    θ_targets = collect(range(π/(nlat+1), stop=π*nlat/(nlat+1), length=nlat))
    _contour_segments(2, θ_targets, periodic=false)
    # Longitude lines (constant azimuth φ)
    φ_targets = collect(range(0, stop=2π - 2π/nlon, length=nlon))
    _contour_segments(3, φ_targets, periodic=true)
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
    colormap="gist_heat", xlim=Float64[], ylim=Float64[], background="black", flipx=false,
    compass=true, rotation_axis=false, rotation_arrow=false, graticules=false)
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
  axis("equal")
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
  colours = get_cmap(colormap).(projmap/maximum(projmap))

  for i=1:star.nquads_visible
  p = patches.Polygon(hcat(-star.proj_west[i,:],star.proj_north[i,:]),closed=true,edgecolor= (plotmesh == true) ? "lightgrey" : colours[i],facecolor=colours[i],fill=true,rasterized=false)
  ax.add_patch(p);
  end
  xlabel(L"x $\leftarrow$ E (mas)", fontsize=20)
  ylabel(L"y $\rightarrow$ N (mas)", fontsize=20)
  compass_color = background == "black" ? "white" : "black"
  if compass; draw_compass(ax, axis_max, color=compass_color); end
  if rotation_axis; draw_rotation_axis(ax, star); end
  if rotation_arrow; draw_rotation_arrow(ax, star); end
  if graticules; draw_graticules(ax, star); end
  cmap=ColorMap(colormap)
  projmap ./= maximum(projmap)  # TODO: this is for intensity -- we will have to rewrite it properly for temperature
  norm = matplotlib.colors.Normalize(vmin=minimum(projmap), vmax=maximum(projmap))
  divider = axdiv.make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.07)
  cb=colorbar(matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap), cax=cax)
  end

"""
    plot2d_binary(tmap1, tmap2, star1, star2, bparams, tepoch; ...)

Plot a binary system on the 2D sky plane with correct occlusion.
The farther star (larger z = receding) is drawn behind the nearer one.
Stars must have been created with `create_binary` (projections include orbital offsets).
"""
function plot2d_binary(tmap1, tmap2, star1, star2, bparams, tepoch;
    plotmesh=false, colormap="gist_heat", pad=1.0, background="black",
    compass=true, rotation_axis=false, rotation_arrow=false, graticules=false, figtitle="")
  set_oiplot_defaults()
  patches = pyimport("matplotlib.patches")
  axdiv = pyimport("mpl_toolkits.axes_grid1.axes_divider")
  facecolor = background == "black" ? "Black" : "White"
  fig = figure("Binary epoch", figsize=(10,10), facecolor="White")
  clf()
  ax = gca()
  title(figtitle)
  ax.set_facecolor(facecolor)
  axis("equal")
  # Axis limits encompassing both stars
  all_x = vcat(vec(star1.vertices_xyz[:,:,1]), vec(star2.vertices_xyz[:,:,1]))
  all_y = vcat(vec(star1.vertices_xyz[:,:,2]), vec(star2.vertices_xyz[:,:,2]))
  axis_max = max(maximum(abs.(all_x)), maximum(abs.(all_y))) + pad
  ax.set_xlim([axis_max, -axis_max])
  ax.set_ylim([-axis_max, axis_max])
  # Shared color normalization
  projmap1 = tmap1[star1.index_quads_visible]
  projmap2 = tmap2[star2.index_quads_visible]
  tmin = min(minimum(projmap1), minimum(projmap2))
  tmax = max(maximum(projmap1), maximum(projmap2))
  if tmin == tmax; tmax = tmin * 1.05 + 1.0; end
  colours1 = get_cmap(colormap).((projmap1 .- tmin) ./ (tmax - tmin))
  colours2 = get_cmap(colormap).((projmap2 .- tmin) ./ (tmax - tmin))
  # Determine z-ordering: farther star (larger z = receding) drawn first (behind)
  _, _, z1, _, _, z2 = binary_orbit_abs(bparams, tepoch)
  zord1 = z1 > z2 ? 2 : 3
  zord2 = z1 > z2 ? 3 : 2
  for i=1:star1.nquads_visible
    ec = plotmesh ? "lightgrey" : colours1[i]
    p = patches.Polygon(hcat(-star1.proj_west[i,:], star1.proj_north[i,:]),
      closed=true, edgecolor=ec, facecolor=colours1[i], fill=true, rasterized=false, zorder=zord1)
    ax.add_patch(p)
  end
  for i=1:star2.nquads_visible
    ec = plotmesh ? "lightgrey" : colours2[i]
    p = patches.Polygon(hcat(-star2.proj_west[i,:], star2.proj_north[i,:]),
      closed=true, edgecolor=ec, facecolor=colours2[i], fill=true, rasterized=false, zorder=zord2)
    ax.add_patch(p)
  end
  xlabel(L"x $\leftarrow$ E (mas)", fontsize=20)
  ylabel(L"y $\rightarrow$ N (mas)", fontsize=20)
  compass_color = background == "black" ? "white" : "black"
  if compass; draw_compass(ax, axis_max, color=compass_color); end
  if rotation_axis; draw_rotation_axis(ax, star1); draw_rotation_axis(ax, star2); end
  if rotation_arrow; draw_rotation_arrow(ax, star1); draw_rotation_arrow(ax, star2); end
  if graticules; draw_graticules(ax, star1); draw_graticules(ax, star2); end
  # Colorbar
  cmap = ColorMap(colormap)
  norm_cb = matplotlib.colors.Normalize(vmin=tmin, vmax=tmax)
  divider = axdiv.make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.07)
  cb = colorbar(matplotlib.cm.ScalarMappable(norm=norm_cb, cmap=cmap), cax=cax)
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
  for i=1:star.nquads_visible
  p = patches.Polygon(hcat(-star.proj_west[i,:],star.proj_north[i,:]),closed=true,edgecolor="black", facecolor="white",rasterized=false)
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
      for i=1:star[t].nquads_visible
        p = patches.Polygon(hcat(-star[t].proj_west[i, :],star[t].proj_north[i, :]),
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
