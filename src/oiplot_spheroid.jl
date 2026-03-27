using PyPlot,PyCall, LaTeXStrings, Statistics

function set_oiplot_defaults()
    PyDict(pyimport("matplotlib")."rcParams")["font.family"]=["serif"]
    PyDict(pyimport("matplotlib")."rcParams")["font.size"]=[14]
    PyDict(pyimport("matplotlib")."rcParams")["xtick.major.size"]=[6]
    PyDict(pyimport("matplotlib")."rcParams")["ytick.major.size"]=[6]
    PyDict(pyimport("matplotlib")."rcParams")["xtick.minor.size"]=[6]
    PyDict(pyimport("matplotlib")."rcParams")["ytick.minor.size"]=[6]
    PyDict(pyimport("matplotlib")."rcParams")["xtick.major.width"]=[1]
    PyDict(pyimport("matplotlib")."rcParams")["ytick.major.width"]=[1]
    PyDict(pyimport("matplotlib")."rcParams")["xtick.minor.width"]=[1]
    PyDict(pyimport("matplotlib")."rcParams")["ytick.minor.width"]=[1]
    PyDict(pyimport("matplotlib")."rcParams")["lines.markeredgewidth"]=[1]
    PyDict(pyimport("matplotlib")."rcParams")["legend.numpoints"]=[1]
    PyDict(pyimport("matplotlib")."rcParams")["legend.handletextpad"]=[0.3]
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

function plot3d_vertices(star)
# 3D view in ROTIR's internal frame: x₁=West on sky, x₂=North on sky, x₃=toward observer
center_xyz = star.vertices_xyz[:,5,:]
corners_xyz = star.vertices_xyz[:,1:4,:]
fig = figure("Center of tessels",figsize=(10,10),facecolor="White");
ax = subplot(projection="3d")
xlabel("West (mas)"); ylabel("North (mas)"); zlabel("toward obs.");
axis_max = maximum(sqrt.(star.vertices_xyz[:,:,1].^2 .+ star.vertices_xyz[:,:,2].^2 .+ star.vertices_xyz[:,:,3].^2))*1.5;
xlim([axis_max,-axis_max]);
ylim([-axis_max,axis_max]);
zlim([-axis_max,axis_max]);
plot3D(center_xyz[:,1],center_xyz[:,2],center_xyz[:,3], ".", color="red");
for i=1:4
  plot3D(corners_xyz[:, i, 1],corners_xyz[:,i, 2],corners_xyz[:,i, 3], ".", color="blue");
end
ax.set_aspect("equal")
PyPlot.draw()
end

function plot2d(tmap, star; intensity = false, figtitle ="", plotmesh=false, pad = 0.5,
    colormap="gist_heat", xlim=Float64[], ylim=Float64[], background="black", flipx=false,
    compass=true, rotation_axis=false)
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
  cmap=ColorMap(colormap)
  projmap ./= maximum(projmap)  # TODO: this is for intensity -- we will have to rewrite it properly for temperature
  norm = matplotlib.colors.Normalize(vmin=minimum(projmap), vmax=maximum(projmap))
  divider = axdiv.make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.07)
  cb=colorbar(matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap), cax=cax)
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
