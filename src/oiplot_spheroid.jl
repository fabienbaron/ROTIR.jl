# gather common display tasks
using PyCall
using PyPlot
using LaTeXStrings

############################################################
#
# Imaging on spheroids
#
############################################################


function plot2Dquad(star_geometry,i) # plots the ith quad projected onto the imaging plane
  projx = star_geometry.projx;
  projy = star_geometry.projy;
  #plots the nth quad in the 2D plane, using ABCD
  # this can be used to debug lots of stuff...
  fig = figure("Test counter",figsize=(10,10),facecolor="White");
  scatter(projx[i,:], projy[i,:]);
  annotate("A", xy=[projx[i,1];projy[i,1]], xycoords="data");
  annotate("B", xy=[projx[i,2];projy[i,2]], xycoords="data");
  annotate("C", xy=[projx[i,3];projy[i,3]], xycoords="data");
  annotate("D", xy=[projx[i,4];projy[i,4]], xycoords="data");
  PyPlot.draw()
  return 1
end

function plot3d_temperature(star_temperature_map,star_geometry) # this plots the temperature map
#  quads_visible = star_geometry.quads_visible;
  #patches = pyimport("matplotlib.mplot3d")
  corners_xyz = star_geometry.vertices_xyz[:,1:4,:];
  Art3D = pyimport("mpl_toolkits.mplot3d.art3d")
  Poly3DCollection = Art3D.Poly3DCollection
  fig2 = figure("Spheroid plot",figsize=(10,10),facecolor="White");
  ax = fig2.add_subplot(111, projection="3d")
  ax = gca(); xlabel("x"); ylabel("y"); zlabel("z"); grid("on");
  axis_max = maximum(sqrt.(star_geometry.vertices_xyz[:,1:4,1].^2 +
    star_geometry.vertices_xyz[:,1:4,2].^2 +
    star_geometry.vertices_xyz[:,1:4,3].^2));
  ax.set_xlim([-axis_max,axis_max]);
  ax.set_ylim([-axis_max,axis_max]);
  ax.set_zlim([-axis_max,axis_max]);
  for i=1:star_geometry.npix
      verts = (collect(zip(corners_xyz[i, :, 1], corners_xyz[i, :, 2], corners_xyz[i, :, 3])),);
      color = get_cmap("gist_heat")(star_temperature_map[i]/maximum(star_temperature_map));
      ax.add_collection3d(Poly3DCollection(verts, edgecolor="none", facecolor=color));
  end
  # maybe repeat the same here with another color for hidden polygons
  PyPlot.draw()
end

function plot3d_temperature_binary(star_temperature_map1,star_temperature_map2,star_geometry1,star_geometry2) # this plots the temperature map
  #  quads_visible = star_geometry.quads_visible;
  #patches = pyimport("matplotlib.mplot3d")
  corners_xyz1 = star_geometry1.vertices_xyz[:,1:4,:];
  corners_xyz2 = star_geometry2.vertices_xyz[:,1:4,:];
  Art3D = pyimport("mpl_toolkits.mplot3d.art3d")
  Poly3DCollection = Art3D[:Poly3DCollection]
  fig2 = figure("Spheroid plot",figsize=(10,10),facecolor="White", projection="3d");
  ax = gca();
  xlabel("x"); ylabel("y"); zlabel("z");
  grid("on");
  axis_max = maximum(sqrt.(star_geometry2.vertices_xyz[:,1:4,1].^2 +
    star_geometry2.vertices_xyz[:,1:4,2].^2 +
    star_geometry2.vertices_xyz[:,1:4,3].^2));
  ax.set_xlim([-axis_max,axis_max]);
  ax.set_ylim([-axis_max,axis_max]);
  ax.set_zlim([-axis_max,axis_max]);
  for i=1:star_geometry1.npix
    # if (quads_visible[i] > 0)
    # star 1
    verts1 = (collect(zip(corners_xyz1[i, :, 1], corners_xyz1[i, :, 2], corners_xyz1[i, :, 3])),);
    color1 = get_cmap("gist_heat")(star_temperature_map1[i]/maximum(star_temperature_map2));
    ax.add_collection3d(Poly3DCollection(verts1, edgecolor="none", facecolor=color1));
    # star 2
    verts2 = (collect(zip(corners_xyz2[i, :, 1], corners_xyz2[i, :, 2], corners_xyz2[i, :, 3])),);
    color2 = get_cmap("gist_heat")(star_temperature_map2[i]/maximum(star_temperature_map2));
    ax.add_collection3d(Poly3DCollection(verts2, edgecolor="none", facecolor=color2));
    # end
  end
    # maybe repeat the same here with another color for hidden polygons
  PyPlot.draw()
end

function plot3d_vertices(star_geometry)
center_xyz = star_geometry.vertices_xyz[:,5,:]
corners_xyz = star_geometry.vertices_xyz[:,1:4,:]
fig = figure("Center of healpixels",figsize=(10,10),facecolor="White");
ax = gca(projection="3d");
axis("equal")
xlabel("x"); ylabel("y"); zlabel("z");
axis_max = maximum(sqrt(star_geometry.vertices_xyz[:,:,1].^2 + star_geometry.vertices_xyz[:,:,2].^2 + star_geometry.vertices_xyz[:,:,3].^2));
xlim([-axis_max,axis_max]);
ylim([-axis_max,axis_max]);
zlim([-axis_max,axis_max]);
plot3D(center_xyz[:,1],center_xyz[:,2],center_xyz[:,3], ".", color="red");
for i=1:4
  plot3D(corners_xyz[:, i, 1],corners_xyz[:,i, 2],corners_xyz[:,i, 3], ".", color="blue");
end
PyPlot.draw()
end

function plot2d_temperature(star_map, star_geometry; plotmesh=true, colormap="gist_heat", xlim=Float64[], ylim=Float64[]) # this plots the temperature map onto the projected 2D image plane (= observer view)
set_oiplot_defaults()
patches = pyimport("matplotlib.patches")
fig = figure("Epoch image",figsize=(10,10),facecolor="Black")
ax = fig.add_axes([0.05,0.05,0.85,0.85])
if plotmesh == true
  meshcolor = "grey"
else
  meshcolor = "none"
end
axis("equal")
if xlim == []
  xlim = [minimum(star_geometry.projx),maximum(star_geometry.projx)]
end
if ylim == []
  ylim = [minimum(star_geometry.projy),maximum(star_geometry.projy)]
end
ax.set_xlim(xlim)
ax.set_ylim(ylim)
projmap = star_map[star_geometry.index_quads_visible];
for i=1:star_geometry.nquads_visible
  p = patches.Polygon(hcat(star_geometry.projx[i,:],star_geometry.projy[i,:]),
  closed=true,edgecolor=meshcolor,facecolor=get_cmap(colormap)(projmap[i]/maximum(projmap)),fill=true,rasterized=false)
  ax.add_patch(p);
end
xlabel("x ← E (mas)")
ylabel("y → N (mas)")
ax.plot();
PyPlot.draw()
end


function plot2d_intensity(star_map, star_geometry; plotmesh=true, colormap="gist_heat") # this plots the temperature map onto the projected 2D image plane (= observer view)
# still missing the actual intensity (includes LD)
patches = pyimport("matplotlib.patches")
fig = figure("Epoch image",figsize=(10,10),facecolor="White")
ax = fig.add_axes([0.05,0.05,0.85,0.85])
if plotmesh == true
  meshcolor = "grey"
else
  meshcolor = "none"
end
xlabel("x");
ylabel("y");
axis("equal");
ax.set_xlim([maximum(star_geometry.projx),minimum(star_geometry.projx)])
ax.set_ylim([minimum(star_geometry.projy),maximum(star_geometry.projy)])
projmap = (star_map.*star_geometry.ldmap)[star_geometry.index_quads_visible]
for i=1:star_geometry.nquads_visible
  p = patches.Polygon(hcat(-star_geometry.projx[i,:],star_geometry.projy[i,:]),
  closed=true,edgecolor=meshcolor,facecolor=get_cmap("gist_heat")(projmap[i]),fill="true",rasterized=false)
  ax.add_patch(p);
end
ax.plot();
PyPlot.draw()
end


function plot2d_temperature_savefig(star_map,star_geometry;file_loc="./image",
    labels=false,iteration="0",omega="?",rotational_vel="?",LD="?",rotation_period="?",
    plotmesh=true,plot_lim=[-3,3],colormap="Blues_r",longlat=true,color_perc=0.75,starname="")
patches = pyimport("matplotlib.patches")
ioff()
fig = figure(figsize=(10,10),facecolor="White")
ax = fig.add_axes([0.1,0.1,0.85,0.85]);
gca().set_axis_on();
xlabel("East (mas)", fontweight="bold", fontsize=15);
ylabel("North (mas)", fontweight="bold", fontsize=15);
plot_lim_r = -1*plot_lim;
ax.set_xlim(plot_lim_r);
ax.set_ylim(plot_lim);
projmap = (star_map.*star_geometry.ldmap)[star_geometry.index_quads_visible];
if (plotmesh == false)
for i=1:star_geometry.nquads_visible
    p = patches.Polygon(hcat(star_geometry.projx[i,:],star_geometry.projy[i,:]),closed=true,
        edgecolor="none",facecolor=get_cmap(colormap)(projmap[i]/maximum(projmap)*color_perc),fill=true,rasterized=false);
    ax.add_patch(p);
end
end
# calculate number of latitudes/longitudes & make grid lines
if ((longlat==true) & (plotmesh == true))
    meshcolor_arr = repeat(["none"],star_geometry.npix,1);
    #=nlongitude = length(star_geometry.vertices_spherical[:,2,5][find(star_geometry.vertices_spherical[:,2,5] .== star_geometry.vertices_spherical[1,2,5])]);
    nlatitude = Int(star_geometry.npix/nlongitude);
    indx = 0;
    for i = 1:nlatitude
        for j = 1:nlongitude
            indx += 1;
            if (j%10 == 0)
                meshcolor_arr[indx] = "grey";
            end
        end
    end=#
    meshcolor_arr[Int.(collect(linspace(1,star_geometry.npix,star_geometry.npix))) %10 .== 0] = "grey";
    meshcolor_proj = meshcolor_arr[star_geometry.index_quads_visible];
    for i=1:star_geometry.nquads_visible
        p = patches.Polygon(hcat(star_geometry.projx[i,:],star_geometry.projy[i,:]),closed=true,
            edgecolor="none",facecolor=get_cmap(colormap)(projmap[i]/maximum(projmap)*color_perc),fill=true,rasterized=false);
        ax.add_patch(p);
        p = patches.Polygon(hcat(star_geometry.projx[i,1:2],star_geometry.projy[i,1:2]),closed=true,
            edgecolor=meshcolor_proj[i],facecolor=get_cmap(colormap)(projmap[i]/maximum(projmap)*color_perc),fill=true,rasterized=false);
        ax.add_patch(p);
    end
end

# adds labels to graph
if (labels == true)
    if (length(string(omega)) > 5)
        ax.text(plot_lim[2]-plot_lim[2]/6.,plot_lim[1]+plot_lim[2]/6.,L"$\omega$ = "*string(omega)[1:6],fontsize=20);
    else
        ax.text(plot_lim[2]-plot_lim[2]/6.,plot_lim[1]+plot_lim[2]/6.,L"$\omega$ = "*string(omega),fontsize=20);
    end
    ax.text(plot_lim[2]-plot_lim[2]/6.,plot_lim[2]-plot_lim[2]/6.,string(rotational_vel)[1:6]*" revolutions/day",fontsize=20);
    ax.text(plot_lim[2]-plot_lim[2]/6.,plot_lim[2]-plot_lim[2]/4.,"Rotation period = "*string(rotation_period)[1:6]*" days",fontsize=20);
    ax.text(plot_lim[1]+plot_lim[2]/2.,plot_lim[1]+plot_lim[2]/6.,string(starname),fontsize=20);
    if (length(string(LD)) > 3)
        ax.text(plot_lim[2]-plot_lim[2]/6.,plot_lim[1]+plot_lim[2]/3.,"LD = "*string(LD)[1:4],fontsize=20);
    else
        ax.text(plot_lim[2]-plot_lim[2]/6.,plot_lim[1]+plot_lim[2]/3.,"LD = "*string(LD),fontsize=20);
    end
end

ax.tick_params(axis="both", which="both", labelsize=15, width=2, length=10);
ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(1));
ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(1));
ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(0.2));
ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(0.2));
ax.tick_params(axis="both", which="minor", width=2, length=5);
# make axes thicker
ax.spines["top"].set_linewidth(2);
ax.spines["bottom"].set_linewidth(2);
ax.spines["left"].set_linewidth(2);
ax.spines["right"].set_linewidth(2);

ax[:plot]();
PyPlot.draw()
if (length(string(iteration)) == 1)
    savefig(file_loc*"000"*string(iteration)*".png");
elseif (length(string(iteration)) == 2)
    savefig(file_loc*"00"*string(iteration)*".png");
elseif (length(string(iteration)) == 3)
    savefig(file_loc*"0"*string(iteration)*".png");
else
    savefig(file_loc*string(iteration)*".png");
end
close()
end

function plot2d_wire(star_geometry) # this plots the temperature map onto the projected 2D image plane (= observer view)
patches = pyimport("matplotlib.patches")
fig = figure("Epoch image",figsize=(10,10),facecolor="White")
ax = fig.add_axes([0.05,0.05,0.85,0.85])
xlabel("x");
ylabel("y");
axis_max = maximum(sqrt(star_geometry.vertices_xyz[:,1:4,1].^2 +star_geometry.vertices_xyz[:,1:4,2].^2 + star_geometry.vertices_xyz[:,1:4,3].^2));
axis("equal")
ax.set_xlim([-axis_max,axis_max]);
ax.set_ylim([-axis_max,axis_max]);
for i=1:star_geometry.nquads_visible
p = patches.Polygon(hcat(star_geometry.projx[i,:],star_geometry.projy[i,:]),closed=true,edgecolor="black", facecolor="white",rasterized=false)
ax.add_patch(p);
end
ax.plot();
PyPlot.draw()
end

function plot2d_temperature_allepochs(star_map, star_geometry; plotmesh=false, tepochs = [], colormap="gist_heat")
patches = pyimport("matplotlib.patches")
fig = figure("Temperature map -- All epochs",figsize=(15,10),facecolor="Black")
if plotmesh == true
  meshcolor = "grey"
else
  meshcolor = "none"
end
subplots_adjust(hspace=0.0)
rows = ceil(Int64,sqrt(length(star_geometry)))
cols = rows
for t=1:length(star_geometry)
  fig.add_subplot(rows,cols,t)
  if tepochs !=[]
    title("Epoch $t $(tepochs[t])") # Give the most recent axis a title
  end
  ax = gca();
  axis("equal")
  ax.set_xlim([2,-2])
  ax.set_ylim([-2,2])
  visible_pixels = sometimes_visible(star_geometry);
  minT = minimum(star_map[visible_pixels]);
  maxT = maximum(star_map[visible_pixels]);
  projmap = (star_map[star_geometry[t].index_quads_visible].-minT)./(maxT-minT);

  for i=1:star_geometry[t].nquads_visible
    p = patches.Polygon(hcat(-star_geometry[t].projx[i,:],star_geometry[t].projy[i,:]),
    closed=true,edgecolor=meshcolor,facecolor=get_cmap(colormap)(projmap[i]),fill="true",rasterized=false)
    ax.add_patch(p);
  end
  ax.plot();
end
fig.canvas.draw() # Update the figure
suptitle("All epochs")
end

function plot2d_intensity_allepochs(star_map, star_geometry; plotmesh=false, tepochs = [])
patches = pyimport("matplotlib.patches")
fig = figure("Intensity map -- All epochs",figsize=(15,10),facecolor="White")
if plotmesh == true
  meshcolor = "grey"
else
  meshcolor = "none"
end

for t=1:length(star_geometry)
  ax= subplot(230+t) # Create the 1st axis of a 2x2 arrax of axes
  if tepochs !=[]
    title("Epoch $t $(tepochs[t])") # Give the most recent axis a title
  end
  axis("equal")
  ax.set_xlim([2,-2])
  ax.set_ylim([-2,2])
  visible_pixels = sometimes_visible(star_geometry);
  minT = minimum(star_map[visible_pixels]);
  maxT = maximum(star_map[visible_pixels]);
  projmap = ( (star_map.-minT)./(maxT-minT).*star_geometry[t].ldmap)[star_geometry[t].index_quads_visible];
  for i=1:star_geometry[t].nquads_visible
    p = patches.Polygon(hcat(-star_geometry[t].projx[i,:],star_geometry[t].projy[i,:]),
    closed=true,edgecolor=meshcolor,facecolor=get_cmap("gist_heat")(projmap[i]),fill="true",rasterized=false)
    ax.add_patch(p);
  end
  ax.plot();
end
fig.canvas.draw() # Update the figure
suptitle("All epochs")
end

function mollplot_temperature_healpix(image; visible_pixels = [], vmin = -Inf, vmax = Inf, colormap="gist_heat", figtitle="Mollweide")
  xsize = 2000
  ysize = div(xsize,2)
  theta = collect(range(pi, stop=0.0, length=ysize))
  phi   = collect(range(-pi, stop=pi, length=xsize))
  longitude = collect(range(-180, stop=180, length=xsize))/180*pi
  latitude = collect(range(-90, stop=90, length=ysize))/180*pi
  # project the map to a rectangular matrix xsize x ysize
  nside = npix2nside(length(image))
  PHI = [i for j in theta, i in phi]
  THETA = [j for j in theta, i in phi]
  grid_pix = reshape(ang2pix_nest(nside, vec(THETA), vec(PHI)), size(PHI))
  grid_map = image[grid_pix]
  fig = figure(figtitle, figsize=(10, 7))
  clf();
  ax = subplot(111,projection="mollweide")
  # rasterized makes the map bitmap while the labels remain vectorial
  # flip longitude to the astro convention
  if visible_pixels == []
      if (vmin == -Inf)
          vmin = minimum(image);
      end
      if (vmax == Inf)
        vmax = maximum(image);
      end
  else
    if (vmin == -Inf)
         vmin = minimum(image[visible_pixels]);
    end
    if (vmax == Inf)
     vmax = maximum(image[visible_pixels]);
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
# class ThetaFormatterShiftPi(GeoAxes.ThetaFormatter):
#     """Shifts labelling by pi
#     Shifts labelling from -180,180 to 0-360"""
#     def __call__(self, x, pos=None):
#         if x != 0:
#             x *= -1
#         if x < 0:
#             x += 2*np.pi
#         return GeoAxes.ThetaFormatter.__call__(self, x, pos)
#
#
# ax[:xaxis][:set_major_formatter](ThetaFormatterShiftPi(60))
# colorbar
ticks = collect(range(vmin, stop=vmax, length=7));
cb = colorbar(moll, orientation="horizontal", shrink=.6, pad=0.05, ticks=ticks)
cb.ax.xaxis.labelpad=5
cb.ax.xaxis.set_label_text("Temperature")

# workaround for issue with viewers, see colorbar docstring
cb.solids.set_edgecolor("face")
ax.tick_params(axis="x", labelsize=15)
ax.tick_params(axis="y", labelsize=15)
end

function mollplot_temperature_longlat(image, ntheta, nphi; visible_pixels = [], colormap="gist_heat", figtitle="Mollweide")
  size = 2000
  ysize = div(xsize,2)
  theta = collect(range(pi, stop=0, length=ysize))
  phi   = collect(range(-pi, stop=pi, length=xsize))
  longitude = collect(range(-180.0, stop=180.0, length=xsize))/180.0*pi
  latitude = collect(range(90.0, stop=-90.0, length=ysize))/180.0*pi
  # project the map to a rectangular matrix xsize x ysize
  PHI = [i for j in theta, i in phi]
  THETA = [j for j in theta, i in phi]
  grid_pix = longlat_ang2pix(ntheta, nphi, THETA, PHI); # for long lat scheme
  #grid_map = image[grid_pix]
  grid_map = image[circshift(grid_pix,(0,Int(xsize/2)))]
  fig = figure(figtitle, figsize=(10, 7))
  clf();
  ax = subplot(111,projection="mollweide",title=title)
  # rasterized makes the map bitmap while the labels remain vectorial
  # flip longitude to the astro convention
  if visible_pixels == []
      vmin = minimum(image);
      vmax = maximum(image);
  else
     vmin = minimum(image[visible_pixels]);
     vmax = maximum(image[visible_pixels]);
  end

  moll = pcolormesh(longitude, latitude, grid_map, vmin=vmin, vmax=vmax, rasterized=true, cmap=colormap)
  # graticule
  ax.set_longitude_grid(30)
  ax.set_latitude_grid(30)
  ax.set_longitude_grid_ends(90)
  spacing = 0.04
  subplots_adjust(bottom=spacing, top=1-spacing, left=spacing, right=1-spacing)
  grid(true)

  ticks = collect(range(vmin, stop=vmax, length=7));
  cb = colorbar(moll, orientation="horizontal", shrink=.6, pad=0.05, ticks=ticks)
  cb.ax.xaxis.labelpad=5
  cb.ax.xaxis.set_label_text("Temperature")
  # workaround for issue with viewers, see colorbar docstring
  cb.solids.set_edgecolor("face")
  ax.tick_params(axis="x", labelsize=15)
  ax.tick_params(axis="y", labelsize=15)
end

#note: this requires healpy
function healpyplot(x; vmin = -1e9, vmax = 1e9)
  if vmin==-1e9
    vmin = minimum(x)
  end
  if vmax==1e9
    vmax = maximum(x)
  end
hp=pyimport("healpy.visufunc")
hp.mollview(x, nest=true, fig=1, min=vmin, max=vmax, rot=0, flip="geo", cbar=false, cmap="gist_heat", title="")
hp.graticule(dpar=10, dmer=10,  force=true, verbose=false, lw=0.25)
hp.cartview(x, nest=true, fig=2, min=vmin, max=vmax, rot=0, flip="geo", cbar=false, cmap="gist_heat", title="")
hp.graticule(dpar=10, dmer=10,  force=true, verbose=false, lw=0.25)
hp.orthview(x, nest=true, fig=3, min=vmin, max=vmax, rot=0, flip="geo", cbar=false, cmap="gist_heat", title="")
hp.graticule(dpar=10, dmer=10,  force=true, verbose=false, lw=0.25)
end
