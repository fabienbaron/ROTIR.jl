# gather common display tasks
using PyCall
using PyPlot
using LazySets
using Statistics
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
# star_temperature_map =secondary_temperature_map; star_geometry= secondary_geom
# star_temperature_map = primary_temperature_map; star_geometry= primary_geom
  corners_xyz = star_geometry.vertices_xyz[:,1:4,:];
  Art3D = pyimport("mpl_toolkits.mplot3d.art3d")
  Poly3DCollection = Art3D.Poly3DCollection
  fig2 = figure("Spheroid plot",figsize=(10,10),facecolor="White");
  ax = subplot(projection="3d")
  #ax = Axes3D(fig2)
  xlabel("x"); ylabel("y"); zlabel("z");
  #grid("off");
  axis_max = maximum(sqrt.(star_geometry.vertices_xyz[:,:,1].^2 .+ star_geometry.vertices_xyz[:,:,2].^2 .+ star_geometry.vertices_xyz[:,:,3].^2))*1.5;
  xlim([axis_max,-axis_max]);
  ylim([-axis_max,axis_max]);
  zlim(bottom=-axis_max,top=axis_max);
 
  for i=1:star_geometry.npix
#    if (quads_visible[i] > 0)
      verts = (collect(zip(corners_xyz[i, :, 1], corners_xyz[i, :, 2], corners_xyz[i, :, 3])),);
      color = get_cmap("gist_heat")((star_temperature_map[i] - minimum(star_temperature_map)) / (maximum(star_temperature_map) - minimum(star_temperature_map)));
      ax.add_collection3d(Poly3DCollection(verts, edgecolor="none", facecolor=color));
#    end
  end
  ax.set_aspect("equal")
  PyPlot.draw()
end

function plot3d_temperature_binary(star_temperature_map1,star_temperature_map2,star_geometry1,star_geometry2) # this plots the temperature map
  #  quads_visible = star_geometry.quads_visible;
  #patches = pyimport("matplotlib.mplot3d")
  corners_xyz1 = star_geometry1.vertices_xyz[:,1:4,:];
  corners_xyz2 = star_geometry2.vertices_xyz[:,1:4,:];
  Art3D = pyimport("mpl_toolkits.mplot3d.art3d")
  Poly3DCollection = Art3D.Poly3DCollection
  # set up figure
  fig2 = figure("Spheroid plot",figsize=(10,10),facecolor="White")
  ax = Axes3D(fig2)
  #ax = gca(projection="3d");
  xlabel("x"); ylabel("y"); zlabel("z");
  grid("on");
  star1_radiusmax = maximum(sqrt.(star_geometry1.vertices_xyz[:,:,1].^2 .+ star_geometry1.vertices_xyz[:,:,2].^2 .+ star_geometry1.vertices_xyz[:,:,3].^2));
  star2_radiusmax = maximum(sqrt.(star_geometry2.vertices_xyz[:,:,1].^2 .+ star_geometry2.vertices_xyz[:,:,2].^2 .+ star_geometry2.vertices_xyz[:,:,3].^2));
  axis_max = (star1_radiusmax + star2_radiusmax)*1.5;
  #axis_max = maximum(sqrt.((star_geometry2.vertices_xyz[:,:,1] .- star_geometry1.vertices_xyz[:,:,1]).^2 .+ (star_geometry2.vertices_xyz[:,:,2] .- star_geometry1.vertices_xyz[:,:,2]).^2.^2 .+ (star_geometry2.vertices_xyz[:,:,3] .- star_geometry1.vertices_xyz[:,:,3]).^2.^2))/2.0;
  xlim([axis_max,-axis_max]);
  ylim([-axis_max,axis_max]);
  zlim(bottom=-axis_max,top=axis_max);

  # set up maximum temperature
  max_temp = maximum([maximum(star_temperature_map1),maximum(star_temperature_map2)])*1.1;
  for i=1:star_geometry1.npix
    # if (quads_visible[i] > 0)
    # star 1
    verts1 = (collect(zip(corners_xyz1[i, :, 1], corners_xyz1[i, :, 2], corners_xyz1[i, :, 3])),);
    color1 = get_cmap("gist_heat")(star_temperature_map1[i]/max_temp);
    ax.add_collection3d(Poly3DCollection(verts1, edgecolor="none", facecolor=color1));
    # star 2
    verts2 = (collect(zip(corners_xyz2[i, :, 1], corners_xyz2[i, :, 2], corners_xyz2[i, :, 3])),);
    color2 = get_cmap("gist_heat")(star_temperature_map2[i]/max_temp);
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
axis_max = maximum(sqrt.(star_geometry.vertices_xyz[:,:,1].^2 .+ star_geometry.vertices_xyz[:,:,2].^2 .+ star_geometry.vertices_xyz[:,:,3].^2))*1.5;
xlim([axis_max,-axis_max]);
ylim([-axis_max,axis_max]);
zlim([-axis_max,axis_max]);
plot3D(center_xyz[:,1],center_xyz[:,2],center_xyz[:,3], ".", color="red");
for i=1:4
  plot3D(corners_xyz[:, i, 1],corners_xyz[:,i, 2],corners_xyz[:,i, 3], ".", color="blue");
end
PyPlot.draw()
end

function plot2d_temperature(star_map, star_geometry, stellar_parameters, oblate_consts; plotmesh=false, colormap="gist_heat", offset=[0.0,0.0,0.0], fig=[], ax=[], draw_polelines=true, poleline_frac=0.35, draw_rotation_arrow=true, rotation_arrow_axis="S", rotation_arrow_axis_fraction=1.0, draw_graticules=true, draw_convex_hull=true) # this plots the temperature map onto the projected 2D image plane (= observer view)
    # still missing the actual intensity (includes LD)
    patches = pyimport("matplotlib.patches")
    collections = pyimport("matplotlib.collections")
    if fig == []
        fig = figure("Epoch image",figsize=(10,10),facecolor="White")
    end
    if ax == []
        ax = fig.add_axes([0.1,0.1,0.85,0.85]);
    end
    ax.set_aspect("equal")
    if plotmesh == true
    meshcolor = "grey"
    else
    meshcolor = "none"
    end
    xlabel("East (mas)", fontweight="bold", fontsize=15);
    ylabel("North (mas)", fontweight="bold", fontsize=15);
    axis_max = maximum(sqrt.(star_geometry.vertices_xyz[:,:,1].^2 .+ star_geometry.vertices_xyz[:,:,2].^2 .+ star_geometry.vertices_xyz[:,:,3].^2))*1.5;
    ax.set_xlim([axis_max,-axis_max])
    ax.set_ylim([-axis_max,axis_max])
    #projmap = star_map[star_geometry.index_quads_visible];
    projmap = (star_map[star_geometry.index_quads_visible] .- minimum(star_map))/(maximum(star_map) - minimum(star_map));
    #projmap = ((star_map[star_geometry.index_quads_visible].-minimum(star_map[star_geometry.index_quads_visible]))./(maximum(star_map[star_geometry.index_quads_visible]) - minimum(star_map[star_geometry.index_quads_visible]))).*star_geometry.ldmap[star_geometry.index_quads_visible];
    for i=1:star_geometry.nquads_visible
    p = patches.Polygon(hcat(-star_geometry.projx[i,:],star_geometry.projy[i,:]),
    closed=true,edgecolor=meshcolor,facecolor=get_cmap(colormap)(projmap[i]),fill=true,rasterized=false, zorder=9)
    ax.add_patch(p);
    #  println(i);
    #  readline()
    #  PyPlot.draw()
    end

    PA = stellar_parameters.position_angle
    i = stellar_parameters.inclination
    rot_angle = stellar_parameters.selfrotangle
    Rx = α -> [1.0 0.0 0.0; 0.0 cos(α) -sin(α); 0.0 sin(α) cos(α)]
    Ry = α -> [cos(α) 0.0 sin(α); 0.0 1.0 0.0; -sin(α) 0.0 cos(α)]
    Rz = α -> [cos(α) -sin(α) 0.0; sin(α) cos(α) 0.0; 0.0 0.0 1.0]

    function compute_polelines(star_geometry)
        north = star_geometry.vertices_xyz[1,1,1:3]
        south = star_geometry.vertices_xyz[end,3,1:3]
        Δcoords = north .- south
        north_end = north .+ Δcoords*poleline_frac
        south_end = south .- Δcoords*poleline_frac
        north_line = [north north_end]
        south_line = [south south_end]
        return north_line, south_line
    end

    function plot_polelines(ax, north_poleline, south_poleline)
        ## Plot polelines
        znorth = north_poleline[3, 2]
        zsouth = south_poleline[3, 2]
        if znorth > 0; nzorder = 10; else; nzorder = 1; end;
        if zsouth > 0; szorder = 10; else; szorder = 1; end;
        Δcoords = north_poleline[:, 2] - north_poleline[:, 1]
        ax.arrow(-north_poleline[1, 1], north_poleline[2, 1], -Δcoords[1], Δcoords[2], color="k", zorder=nzorder, linewidth=1, length_includes_head=true, head_width=0.075, head_length=0.1125)
        ax.plot(-south_poleline[1, :], south_poleline[2, :], "k-", zorder=szorder)
    end

    north_poleline, south_poleline = compute_polelines(star_geometry)
    Δcoords = north_poleline[:, 2] - north_poleline[:, 1]
    if (i in [0.0 180.0]) == false && draw_polelines == true
        plot_polelines(ax, north_poleline, south_poleline)
    end

    function compute_rotation_arrow(c, r, i, PA)
        if 80.0 < i <= 90.0 
            i = 80.0;
        elseif 90.0 < i < 100.0
            i = 100.0
        end

        θ = collect(range(-pi/4, stop=-7pi/4, length=201))
        x = r*cos.(θ)
        y = r*sin.(θ)
        z = zeros(length(θ))
        coords = [x y z]'

        # Rz = [cos(-PA*pi/180) -sin(-PA*pi/180) 0.0; sin(-PA*pi/180) cos(-PA*pi/180) 0.0; 0.0 0.0 1.0]
        # Rx = [1.0 0.0 0.0; 0.0 cos(-i*pi/180) -sin(-i*pi/180); 0.0 sin(-i*pi/180) cos(-i*pi/180)]

        for n=1:length(θ)
            coords[:, n] = Rx(-i*pi/180) * coords[:, n]
            coords[:, n] = Rz(-PA*pi/180) * coords[:, n]
        end

        coords[1, :] .+= -c[1]
        coords[2, :] .+= c[2]
        # coords[3, :] .+= zc

        tip_angle = atan((coords[2, end] - coords[2, end-1]), (coords[1, end] - coords[1, end-1]))
        tip_radius = poleline_frac/10
        tip_coords = coords[1:2, end]
        tip_coords[1] += tip_radius/2*cos(tip_angle)
        tip_coords[2] += tip_radius/2*sin(tip_angle)
        return coords, tip_coords, tip_radius, tip_angle
    end

    function plot_rotation_arrow(ax, rotation_line_coords, tip_coords, tip_radius, tip_angle)
        if rotation_line_coords[3, end] > 0.0; zorder = 10; else; zorder = 1; end;
        ax.plot(rotation_line_coords[1, rotation_line_coords[3, :] .> 0.0], rotation_line_coords[2, rotation_line_coords[3, :] .> 0.0], "-", color="k", linewidth=1, zorder=10)
        ax.plot(rotation_line_coords[1, rotation_line_coords[3, :] .< 0.0], rotation_line_coords[2, rotation_line_coords[3, :] .< 0.0], "--", color="k", linewidth=1, zorder=1)
        ax.add_patch(patches.RegularPolygon((tip_coords[1], tip_coords[2]), 3, tip_radius, tip_angle-pi/2, color="k", zorder=zorder))
    end

    if rotation_arrow_axis == "N"; 
        rotation_ellipse_center = star_geometry.vertices_xyz[1,1,1:3] .+ rotation_arrow_axis_fraction*Δcoords; 
    elseif rotation_arrow_axis == "S"; 
        rotation_ellipse_center = star_geometry.vertices_xyz[end,3,1:3] .- rotation_arrow_axis_fraction*Δcoords;
    end

    rotation_line_coords, tip_coords, tip_radius, tip_angle = compute_rotation_arrow(rotation_ellipse_center, poleline_frac, i, PA)
    if draw_rotation_arrow == true
        plot_rotation_arrow(ax, rotation_line_coords, tip_coords, tip_radius, tip_angle)
    end

    function compute_graticules(oblate_consts, incl, PA, rot_angle)
        # Rz = [cos(-PA*pi/180) -sin(-PA*pi/180) 0.0; sin(-PA*pi/180) cos(-PA*pi/180) 0.0; 0.0 0.0 1.0]
        # Rx = [1.0 0.0 0.0; 0.0 cos((90-incl)*pi/180) -sin((90-incl)*pi/180); 0.0 sin((90-incl)*pi/180) cos((90-incl)*pi/180)]
        # Ry = [cos(rot_angle*pi/180) 0.0 sin(rot_angle*pi/180); 0.0 1.0 0.0; -sin(rot_angle*pi/180) 0.0 cos(rot_angle*pi/180)]

        # Rx = α -> [1.0 0.0 0.0; 0.0 cos(α) -sin(α); 0.0 sin(α) cos(α)]
        # Ry = α -> [cos(α) 0.0 sin(α); 0.0 1.0 0.0; -sin(α) 0.0 cos(α)]
        # Rz = α -> [cos(α) -sin(α) 0.0; sin(α) cos(α) 0.0; 0.0 0.0 1.0]

        a, b, c = oblate_consts

        r = ones(200) * b
        θs = vec([-pi/3 -pi/6 0.0 pi/6 pi/3])
        ϕ = collect(range(-pi/2, stop=3pi/2, length=200))
        latlong_lines = []
        for i=eachindex(θs)
            θ = ones(200) * θs[i]
            coords_rθϕ = [r θ ϕ]'
            sortixs = sortperm(ϕ)
            x = a * cos.(θ) .* cos.(ϕ)
            y = c * sin.(θ)
            z = b * cos.(θ) .* sin.(ϕ)
            coords_xyz = [x y z]'
            # coords_xyz = coords_xyz[:, sortixs]
            for j=axes(coords_xyz, 2)
                coords_xyz[:, j] = Rx((90.0-i)*pi/180) * coords_xyz[:, j]
                coords_xyz[:, j] = Rz(-PA*pi/180) * coords_xyz[:, j]
            end
            if length(coords_xyz[1, coords_xyz[3, :] .> 0.0]) > 0
                coords_xyz = coords_xyz[:, coords_xyz[3, :] .> 0.0]
                push!(latlong_lines, hcat(coords_xyz[1, :], coords_xyz[2, :]))
            end
        end

        θ = collect(range(-pi/2, stop=pi/2, length=200))
        ϕs = collect(range(0.0, stop=2pi, step=pi/4))[1:end-1]
        for i=eachindex(ϕs)
            ϕ  = ones(200) * ϕs[i]
            x = a * cos.(θ) .* cos.(ϕ)
            y = c * sin.(θ)
            z = b * cos.(θ) .* sin.(ϕ)
            coords_xyz = [x y z]'
            for j=axes(coords_xyz, 2)
                coords_xyz[:, j] = Ry(rot_angle*pi/180) * coords_xyz[:, j]
                coords_xyz[:, j] = Rx((90.0-i)*pi/180) * coords_xyz[:, j]
                coords_xyz[:, j] = Rz(-PA*pi/180) * coords_xyz[:, j]
            end
            if length(coords_xyz[1, coords_xyz[3, :] .> 0.0]) > 0
                coords_xyz = coords_xyz[:, coords_xyz[3, :] .> 0.0]
                push!(latlong_lines, hcat(coords_xyz[1, :], coords_xyz[2, :]))
            end

            ## Longitude labels
            # x = a * cos(ϕs[i])
            # y = 0.0
            # z = b * sin(ϕs[i])
            # coords_xyz = [x y z]'
            # coords_xyz = Ry * coords_xyz
            # coords_xyz = Rx * coords_xyz
            # coords_xyz = Rz * coords_xyz
            # ϕ_label = mod(ϕs[i] - pi, 2pi)
            # ϕ_label = round(Int, ϕ_label*180/pi)
            # if z >= 0
            #     ax.text(coords_xyz[1], coords_xyz[2], "$(ϕ_label)°", zorder=13)
            # end
        end
        return latlong_lines
    end

    function plot_graticules(ax, graticules)
        ax.add_collection(collections.PolyCollection(graticules, closed=false, ec="k", fc="none", zorder=13))
    end

    graticules = compute_graticules(oblate_consts, i, PA, rot_angle)
    if draw_graticules == true
        plot_graticules(ax, graticules)
    end

    function compute_convex_hull(star_geometry)
        xpts = vec(star_geometry.vertices_xyz[:, :, 1])
        ypts = vec(star_geometry.vertices_xyz[:, :, 2])
        pts = [vec([xpts[i] ypts[i]]) for i=eachindex(xpts)]
        chull = convex_hull(pts)
        chull = mapreduce(permutedims, vcat, chull)
        return chull
    end

    function plot_convex_hull(ax, chull)
        chull_polygon = patches.Polygon(hcat(-chull[:, 1], chull[:, 2]), closed=true, fill=false, ec="k", rasterized=false, zorder=13)
        ax.add_patch(chull_polygon)
    end
    
    chull = compute_convex_hull(star_geometry)
    ## Plot Convex Hull
    if draw_convex_hull == true
        plot_convex_hull(ax, chull)
    end

    ax.tick_params(axis="both", which="both", labelsize=15, width=2, length=10);
    if (ceil(axis_max) <= 3.0)
        long_tick = 1.0;
        short_tick = 0.1;
    elseif ((ceil(axis_max) > 3.0) & (ceil(axis_max) <= 6.0))
        long_tick = 2.0;
        short_tick = 0.2;
    elseif (ceil(axis_max) > 6.0)
        long_tick = 3.0;
        short_tick = 0.5;
    end
    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(long_tick));
    ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(long_tick));
    ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(short_tick));
    ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(short_tick));
    ax.tick_params(axis="both", which="minor", width=2, length=5);
    # make axes thicker
    ax.spines["top"].set_linewidth(2);
    ax.spines["bottom"].set_linewidth(2);
    ax.spines["left"].set_linewidth(2);
    ax.spines["right"].set_linewidth(2);

    ax.set_xlim(ax.get_xlim().*[1.1, 1.1])
    ax.set_ylim(ax.get_ylim().*[1.1, 1.1])

    ax.plot();
    PyPlot.draw()
end

function plot2d_temperature_cmap(star_map, star_geometry; plotmesh=true, minmaxT=[], colormap="gist_heat",offset=[0.0,0.0,0.0]) # this plots the temperature map onto the projected 2D image plane (= observer view)
  # still missing the actual intensity (includes LD)
  patches = pyimport("matplotlib.patches")
  fig = figure("Epoch image",figsize=(12,10),facecolor="White")
  ax = fig.add_axes([0.1,0.1,0.85,0.85]);
  if plotmesh == true
    meshcolor = "grey"
  else
    meshcolor = "none"
  end
  xlabel("East (mas)", fontweight="bold", fontsize=15);
  ylabel("North (mas)", fontweight="bold", fontsize=15);
  axis_max = maximum(sqrt.(star_geometry.vertices_xyz[:,:,1].^2 .+ star_geometry.vertices_xyz[:,:,2].^2 .+ star_geometry.vertices_xyz[:,:,3].^2))*1.5;
  ax.set_xlim([axis_max,-axis_max])
  ax.set_ylim([-axis_max,axis_max])
  #projmap = star_map[star_geometry.index_quads_visible];
  projmap = (star_map[star_geometry.index_quads_visible] .- minimum(star_map))/(maximum(star_map) - minimum(star_map));
  #projmap = ((star_map[star_geometry.index_quads_visible].-minimum(star_map[star_geometry.index_quads_visible]))./(maximum(star_map[star_geometry.index_quads_visible]) - minimum(star_map[star_geometry.index_quads_visible]))).*star_geometry.ldmap[star_geometry.index_quads_visible];

  #p_list = [];

  for i=1:star_geometry.nquads_visible
    p = patches.Polygon(hcat(-star_geometry.projx[i,:],star_geometry.projy[i,:]),
    closed=true,edgecolor=meshcolor,facecolor=get_cmap(colormap)(projmap[i]),fill=true,rasterized=false)
    ax.add_patch(p);
    #col = PatchCollection(p)
    #append!(p_list, star_map[star_geometry.index_quads_visible][i]);
  #  println(i);
  #  readline()
  #  PyPlot.draw()
  end

  ax.tick_params(axis="both", which="both", labelsize=15, width=2, length=10);
  if (ceil(axis_max) <= 3.0)
    long_tick = 1.0;
    short_tick = 0.1;
  elseif ((ceil(axis_max) > 3.0) & (ceil(axis_max) <= 6.0))
    long_tick = 2.0;
    short_tick = 0.2;
  elseif (ceil(axis_max) > 6.0)
    long_tick = 3.0;
    short_tick = 0.5;
  end
  ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(long_tick));
  ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(long_tick));
  ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(short_tick));
  ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(short_tick));
  ax.tick_params(axis="both", which="minor", width=2, length=5);
  # make axes thicker
  ax.spines["top"].set_linewidth(2);
  ax.spines["bottom"].set_linewidth(2);
  ax.spines["left"].set_linewidth(2);
  ax.spines["right"].set_linewidth(2);

  if (minmaxT == [])
    vmin = minimum(star_map[star_geometry.index_quads_visible]);
    vmax = maximum(star_map[star_geometry.index_quads_visible]);
  else
    vmin = minmaxT[1];
    vmax = minmaxT[2];
  end

  #moll = pcolormesh(star_map[star_geometry.index_quads_visible], vmin=vmin, vmax=vmax, rasterized=true, cmap=colormap)
  sm = plt.cm.ScalarMappable(cmap=colormap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
  sm.set_array([])
  plt.colorbar(sm)

  ax.plot();
  PyPlot.draw()
end




function plot2d_temperature_savefig(star_map,star_geometry;file_loc="./",labels="none",
  iteration="0",omega="?",rotational_vel="?",LD="?",rotation_period="?",
  plotmesh=true,colormap="Blues_r",longlat=true,starname="")

  patches = pyimport("matplotlib.patches")
  ioff()
  fig = figure(figsize=(10,10),facecolor="White")
  ax = fig.add_axes([0.1,0.1,0.85,0.85]);
  gca().set_axis_on();
  xlabel("East (mas)", fontweight="bold", fontsize=15);
  ylabel("North (mas)", fontweight="bold", fontsize=15);
  axis_max = maximum(sqrt.(star_geometry.vertices_xyz[:,:,1].^2 .+ star_geometry.vertices_xyz[:,:,2].^2 .+ star_geometry.vertices_xyz[:,:,3].^2))*1.5;
  ax.set_xlim([axis_max,-axis_max])
  ax.set_ylim([-axis_max,axis_max])
  projmap = ((star_map.*star_geometry.ldmap)[star_geometry.index_quads_visible] .- minimum((star_map.*star_geometry.ldmap)[star_geometry.index_quads_visible]))/(maximum((star_map.*star_geometry.ldmap)[star_geometry.index_quads_visible]) - minimum((star_map.*star_geometry.ldmap)[star_geometry.index_quads_visible]));
  if (plotmesh == false)
    for i=1:star_geometry.nquads_visible
      p = patches.Polygon(hcat(-star_geometry.projx[i,:],star_geometry.projy[i,:]),closed=true,
        edgecolor="none",facecolor=get_cmap(colormap)(projmap[i]),fill=true,rasterized=false);
      ax.add_patch(p);
    end
  end
  # calculate number of latitudes/longitudes & make grid lines
  if ((longlat==true) & (plotmesh == true))
    meshcolor_arr = repeat(["none"],star_geometry.npix,1);
    meshcolor_arr[Int.(collect(range(1,stop=star_geometry.npix,length=star_geometry.npix))) .%10 .== 0] .= "grey";
    meshcolor_proj = meshcolor_arr[star_geometry.index_quads_visible];
    for i=1:star_geometry.nquads_visible
      p = patches.Polygon(hcat(-star_geometry.projx[i,:],star_geometry.projy[i,:]),closed=true,
        edgecolor="none",facecolor=get_cmap(colormap)(projmap[i]),fill=true,rasterized=false);
      ax.add_patch(p);
      p = patches.Polygon(hcat(-star_geometry.projx[i,1:2],star_geometry.projy[i,1:2]),closed=true,
        edgecolor=meshcolor_proj[i],facecolor=get_cmap(colormap)(projmap[i]),fill=true,rasterized=false);
      ax.add_patch(p);
    end
  end

  # adds labels to graph
  if (labels == "none")
    println("Plotting with no labels on graph");
  elseif (labels == "rapid_rotator")
    println("Plotting with labels used for rapid rotators");
    if (length(string(omega)) > 5)
      ax.text(axis_max-axis_max/6.0,-axis_max+axis_max/6.0,L"$\omega$ = "*string(omega)[1:6],fontsize=20);
    else
      ax.text(axis_max-axis_max/6.0,-axis_max+axis_max/6.0,L"$\omega$ = "*string(omega),fontsize=20);
    end
    ax.text(axis_max-axis_max/6.0,axis_max-axis_max/6.0,string(rotational_vel)[1:6]*" revolutions/day",fontsize=20);
    ax.text(axis_max-axis_max/6.0,axis_max-axis_max/4.0,"Rotation period = "*string(rotation_period)[1:6]*" days",fontsize=20);
    ax.text(-axis_max+axis_max/2.0,-axis_max+axis_max/6.0,string(starname),fontsize=20);
    if (length(string(LD)) > 3)
      ax.text(axis_max-axis_max/6.0,-axis_max+axis_max/3.0,"LD = "*string(LD)[1:4],fontsize=20);
    else
      ax.text(axis_max-axis_max/6.0,-axis_max+axis_max/3.0,"LD = "*string(LD),fontsize=20);
    end
  elseif (labels == "spherical")
    println("Plotting with labels used for 'spherical' star");
    ax.text(axis_max-axis_max/6.0,axis_max-axis_max/6.0,string(rotational_vel)[1:6]*" revolutions/day",fontsize=20);
    ax.text(axis_max-axis_max/6.0,axis_max-axis_max/4.0,"Rotation period = "*string(rotation_period)[1:6]*" days",fontsize=20);
    ax.text(-axis_max+axis_max/2.0,-axis_max+axis_max/6.0,string(starname),fontsize=20);
    if (length(string(LD)) > 3)
      ax.text(axis_max-axis_max/6.0,-axis_max+axis_max/3.0,"LD = "*string(LD)[1:4],fontsize=20);
    else
      ax.text(axis_max-axis_max/6.0,-axis_max+axis_max/3.0,"LD = "*string(LD),fontsize=20);
    end
  end

  ax.tick_params(axis="both", which="both", labelsize=15, width=2, length=10);
  if (ceil(axis_max) <= 2.5)
    long_tick = 1.0;
    short_tick = 0.1;
  elseif ((ceil(axis_max) > 2.5) & (ceil(axis_max) <= 6.0))
    long_tick = 2.0;
    short_tick = 0.2;
  elseif (ceil(axis_max) > 6.0)
    long_tick = 3.0;
    short_tick = 0.5;
  end
  ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(long_tick));
  ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(long_tick));
  ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(short_tick));
  ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(short_tick));
  ax.tick_params(axis="both", which="minor", width=2, length=5);
  # make axes thicker
  ax.spines["top"].set_linewidth(2);
  ax.spines["bottom"].set_linewidth(2);
  ax.spines["left"].set_linewidth(2);
  ax.spines["right"].set_linewidth(2);

  ax.plot();
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
  ax = fig.add_axes([0.1,0.1,0.85,0.85]);
  xlabel("X (mas)", fontweight="bold", fontsize=15);
  ylabel("Y (mas)", fontweight="bold", fontsize=15);
  axis_max = maximum(sqrt.(star_geometry.vertices_xyz[:,:,1].^2 .+ star_geometry.vertices_xyz[:,:,2].^2 .+ star_geometry.vertices_xyz[:,:,3].^2))*1.5;
  axis("equal")
  ax.set_xlim([axis_max,-axis_max]);
  ax.set_ylim([-axis_max,axis_max]);
  for i=1:star_geometry.nquads_visible
  p = patches.Polygon(hcat(-star_geometry.projx[i,:],star_geometry.projy[i,:]),closed=true,edgecolor="black", facecolor="white",rasterized=false)
  ax.add_patch(p);
  end
  ax.tick_params(axis="both", which="both", labelsize=15, width=2, length=10);
  if (ceil(axis_max) <= 3.0)
    long_tick = 1.0;
    short_tick = 0.1;
  elseif ((ceil(axis_max) > 3.0) & (ceil(axis_max) <= 6.0))
    long_tick = 2.0;
    short_tick = 0.2;
  elseif (ceil(axis_max) > 6.0)
    long_tick = 3.0;
    short_tick = 0.5;
  end
  ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(long_tick));
  ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(long_tick));
  ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(short_tick));
  ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(short_tick));
  ax.tick_params(axis="both", which="minor", width=2, length=5);
  # make axes thicker
  ax.spines["top"].set_linewidth(2);
  ax.spines["bottom"].set_linewidth(2);
  ax.spines["left"].set_linewidth(2);
  ax.spines["right"].set_linewidth(2);
  ax.plot();
  PyPlot.draw()
end


function plot2d_temperature_allepochs(star_map, star_geometry; plotmesh=false, tepochs = [], colormap="gist_heat",arr_box=23)
    patches = pyimport("matplotlib.patches")
    fig = figure("Temperature map -- All epochs",figsize=(15,10),facecolor="White")
    if plotmesh == true
      meshcolor = "grey"
    else
      meshcolor = "none"
    end
    #subplots_adjust(hspace=0.0)
    #rows = ceil(Int64,sqrt(length(star_geometry)))
    #cols = rows

    # minT = minimum(minimum.(star_map));
    # maxT = maximum(maximum.(star_map));
    minT = minimum(star_map);
    maxT = maximum(star_map);

    for t=1:length(star_geometry)
      #fig.add_subplot(rows,cols,t)
      ax = subplot(arr_box*10+t) # Create the 1st axis of a 2x2 arrax of axes
      if tepochs !=[]
        title("Epoch $t $(tepochs[t])",fontweight="bold") # Give the most recent axis a title
      end
      ax = gca();
      axis("equal")
      axis_max = maximum(sqrt.(star_geometry[t].vertices_xyz[:,:,1].^2 .+ star_geometry[t].vertices_xyz[:,:,2].^2 .+ star_geometry[t].vertices_xyz[:,:,3].^2))*1.5;
      ax.set_xlim([axis_max,-axis_max])
      ax.set_ylim([-axis_max,axis_max])
      visible_pixels = sometimes_visible(star_geometry);
      projmap = (star_map[star_geometry[t].index_quads_visible].-minT)./(maxT-minT);
      for i=1:star_geometry[t].nquads_visible
        p = patches.Polygon(hcat(-star_geometry[t].projx[i, :],star_geometry[t].projy[i, :]),
        closed=true,edgecolor=meshcolor,facecolor=get_cmap(colormap)(projmap[i]),fill="true",rasterized=false)
        ax.add_patch(p);
      end
      ax.plot();
      ax.tick_params(axis="both", which="both", labelsize=15, width=2, length=10);
      if (ceil(axis_max) <= 3.0)
        long_tick = 1.0;
        short_tick = 0.1;
      elseif ((ceil(axis_max) > 3.0) & (ceil(axis_max) <= 6.0))
        long_tick = 2.0;
        short_tick = 0.2;
      elseif (ceil(axis_max) > 6.0)
        long_tick = 3.0;
        short_tick = 0.5;
      end
      ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(long_tick));
      ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(long_tick));
      ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(short_tick));
      ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(short_tick));

      ax.tick_params(axis="both", which="minor", width=2, length=5);
      # make axes thicker
      ax.spines["top"].set_linewidth(2);
      ax.spines["bottom"].set_linewidth(2);
      ax.spines["left"].set_linewidth(2);
      ax.spines["right"].set_linewidth(2);
    end
    fig.canvas.draw() # Update the figure
    #suptitle("All epochs")

    fig.text(0.5, 0.04, "East (mas)", ha="center", va="center",fontweight="bold", fontsize=15)
    fig.text(0.06, 0.5, "North (mas)", ha="center", va="center", rotation="vertical", fontweight="bold", fontsize=15)
end



# use of star_map has temperature maps in a (n,m) format where n is the pixels and m is the # of epochs
function plot2d_temperature_allcompound_epochs(star_maps, star_geometry; plotmesh=false, tepochs = [], minmaxT = [], colormap="gist_heat",arr_box=23)
  patches = pyimport("matplotlib.patches")
  fig = figure("Temperature map -- All epochs",figsize=(15,10),facecolor="White")
  if plotmesh == true
    meshcolor = "grey"
  else
    meshcolor = "none"
  end

  if (minmaxT == [])
    minT = minimum(minimum.(star_maps));
    maxT = maximum(maximum.(star_maps));
  else
    minT = minmaxT[1];
    maxT = minmaxT[2];
  end

  # minT = minimum(star_maps);
  # maxT = maximum(star_maps);

  for t=1:length(star_geometry)
    ax= subplot(arr_box*10+t) # Create the 1st axis of a 2x2 arrax of axes
    if tepochs !=[]
      epoch_title = string(tepochs[t])[1:10];
      title("Epoch $t "*epoch_title,fontweight="bold") # Give the most recent axis a title
    end
    axis_max = maximum(sqrt.(star_geometry[t].vertices_xyz[:,:,1].^2 .+ star_geometry[t].vertices_xyz[:,:,2].^2 .+ star_geometry[t].vertices_xyz[:,:,3].^2))*1.5;
    ax.set_xlim([axis_max,-axis_max])
    ax.set_ylim([-axis_max,axis_max])
    visible_pixels = sometimes_visible(star_geometry);
    projmap = (star_maps[:,t][star_geometry[t].index_quads_visible].-minT)./(maxT-minT);
    for i=1:star_geometry[t].nquads_visible
      p = patches.Polygon(hcat(-star_geometry[t].projx[i,:],star_geometry[t].projy[i,:]),
      closed=true,edgecolor=meshcolor,facecolor=get_cmap(colormap)(projmap[i]),fill="true",rasterized=false)
      ax.add_patch(p);
    end
    ax.plot();
    ax.tick_params(axis="both", which="both", labelsize=15, width=2, length=10);
    if (ceil(axis_max) <= 3.0)
      long_tick = 1.0;
      short_tick = 0.1;
    elseif ((ceil(axis_max) > 3.0) & (ceil(axis_max) <= 6.0))
      long_tick = 2.0;
      short_tick = 0.2;
    elseif (ceil(axis_max) > 6.0)
      long_tick = 3.0;
      short_tick = 0.5;
    end
    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(long_tick));
    ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(long_tick));
    ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(short_tick));
    ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(short_tick));
    ax.tick_params(axis="both", which="minor", width=2, length=5);
    # make axes thicker
    ax.spines["top"].set_linewidth(2);
    ax.spines["bottom"].set_linewidth(2);
    ax.spines["left"].set_linewidth(2);
    ax.spines["right"].set_linewidth(2);
  end
  fig.canvas.draw() # Update the figure
  #suptitle("All epochs")

  fig.text(0.5, 0.04, "East (mas)", ha="center", va="center",fontweight="bold", fontsize=15)
  fig.text(0.06, 0.5, "North (mas)", ha="center", va="center", rotation="vertical", fontweight="bold", fontsize=15)

  cbar_ax = fig.add_axes([0.92, 0.15, 0.05, 0.7])
  sm = plt.cm.ScalarMappable(cmap=colormap, norm=plt.Normalize(vmin=minT, vmax=maxT))
  sm.set_array([])
  plt.colorbar(sm,cax=cbar_ax,orientation="vertical",shrink=0.6)
end


function plot2d_temperature_allepochs_poleline(star_map, star_geometry; plotmesh=false, tepochs = [], colormap="gist_heat",arr_box=23)
  patches = pyimport("matplotlib.patches")
  fig = figure("Temperature map -- All epochs",figsize=(15,10),facecolor="White")
  if plotmesh == true
    meshcolor = "grey"
  else
    meshcolor = "none"
  end

  # minT = minimum(minimum.(star_map));
  # maxT = maximum(maximum.(star_map));
  minT = minimum(star_map);
  maxT = maximum(star_map);

  for t=1:length(star_geometry)
    ax= subplot(arr_box*10+t) # Create the 1st axis of a 2x2 arrax of axes
    if tepochs !=[]
      title("Epoch $t $(tepochs[t])",fontweight="bold") # Give the most recent axis a title
    end
    axis_max = maximum(sqrt.(star_geometry[t].vertices_xyz[:,:,1].^2 .+ star_geometry[t].vertices_xyz[:,:,2].^2 .+ star_geometry[t].vertices_xyz[:,:,3].^2))*1.5;
    ax.set_xlim([axis_max,-axis_max])
    ax.set_ylim([-axis_max,axis_max])
    visible_pixels = sometimes_visible(star_geometry);
    projmap = (star_map[star_geometry[t].index_quads_visible].-minT)./(maxT-minT);
    # TBD: Still needs to have pole lines partially hidden if behind star (add conditional statement?)
    for i=1:star_geometry[t].nquads_visible
        north = star_geometry[t].vertices_xyz[1,1,1:2]
        south = star_geometry[t].vertices_xyz[end,3,1:2]
        #lines = Any[collect(zip(-[north[1], south[1]],[north[2], south[2]]))];
        lines = Any[collect(zip(-[1.5*north[1], north[1]],[1.5*north[2], north[2]]))];
        push!(lines,collect(zip(-[south[1], 1.5*south[1]],[south[2], 1.5*south[2]])));
        #poleline = matplotlib.collections.LineCollection(Any[collect(zip(-[north[1], south[1]],[north[2], south[2]]))])
        poleline = matplotlib.collections.LineCollection(lines);
        ax.add_collection(poleline);
      p = patches.Polygon(hcat(-star_geometry[t].projx[i,:],star_geometry[t].projy[i,:]),
      closed=true,edgecolor=meshcolor,facecolor=get_cmap(colormap)(projmap[i]),fill="true",rasterized=false)
      ax.add_patch(p);

    end
    ax.plot();
    ax.tick_params(axis="both", which="both", labelsize=15, width=2, length=10);
    if (ceil(axis_max) <= 3.0)
      long_tick = 1.0;
      short_tick = 0.1;
    elseif ((ceil(axis_max) > 3.0) & (ceil(axis_max) <= 6.0))
      long_tick = 2.0;
      short_tick = 0.2;
    elseif (ceil(axis_max) > 6.0)
      long_tick = 3.0;
      short_tick = 0.5;
    end
    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(long_tick));
    ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(long_tick));
    ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(short_tick));
    ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(short_tick));
    ax.tick_params(axis="both", which="minor", width=2, length=5);
    # make axes thicker
    ax.spines["top"].set_linewidth(2);
    ax.spines["bottom"].set_linewidth(2);
    ax.spines["left"].set_linewidth(2);
    ax.spines["right"].set_linewidth(2);
  end
  fig.canvas.draw() # Update the figure
  #suptitle("All epochs")

  fig.text(0.5, 0.04, "East (mas)", ha="center", va="center",fontweight="bold", fontsize=15)
  fig.text(0.06, 0.5, "North (mas)", ha="center", va="center", rotation="vertical", fontweight="bold", fontsize=15)
end


function plot2d_temperature_allepochs_cmap(star_map, star_geometry; plotmesh=false, tepochs = [], minmaxT = [], colormap="gist_heat",arr_box=23)
  patches = pyimport("matplotlib.patches")
  fig = figure("Temperature map -- All epochs",figsize=(15,10),facecolor="White")
  if plotmesh == true
    meshcolor = "grey"
  else
    meshcolor = "none"
  end

  # minT = minimum(minimum.(star_map));
  # maxT = maximum(maximum.(star_map));
  if (minmaxT == [])
    minT = minimum(star_map);
    maxT = maximum(star_map);
  else
    minT = minmaxT[1];
    maxT = minmaxT[2];
  end

  for t=1:length(star_geometry)
    ax= subplot(arr_box*10+t) # Create the 1st axis of a 2x2 arrax of axes
    if tepochs !=[]
      epoch_title = string(tepochs[t])[1:10];
      title("Epoch $t "*epoch_title,fontweight="bold") # Give the most recent axis a title
    end
    axis_max = maximum(sqrt.(star_geometry[t].vertices_xyz[:,:,1].^2 .+ star_geometry[t].vertices_xyz[:,:,2].^2 .+ star_geometry[t].vertices_xyz[:,:,3].^2))*1.5;
    ax.set_xlim([axis_max,-axis_max])
    ax.set_ylim([-axis_max,axis_max])
    visible_pixels = sometimes_visible(star_geometry);
    projmap = (star_map[star_geometry[t].index_quads_visible].-minT)./(maxT-minT);
    for i=1:star_geometry[t].nquads_visible
      p = patches.Polygon(hcat(-star_geometry[t].projx[i,:],star_geometry[t].projy[i,:]),
      closed=true,edgecolor=meshcolor,facecolor=get_cmap(colormap)(projmap[i]),fill="true",rasterized=false)
      ax.add_patch(p);
    end
    ax.plot();
    ax.tick_params(axis="both", which="both", labelsize=15, width=2, length=10);
    if (ceil(axis_max) <= 3.0)
      long_tick = 1.0;
      short_tick = 0.1;
    elseif ((ceil(axis_max) > 3.0) & (ceil(axis_max) <= 6.0))
      long_tick = 2.0;
      short_tick = 0.2;
    elseif (ceil(axis_max) > 6.0)
      long_tick = 3.0;
      short_tick = 0.5;
    end
    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(long_tick));
    ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(long_tick));
    ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(short_tick));
    ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(short_tick));
    ax.tick_params(axis="both", which="minor", width=2, length=5);
    # make axes thicker
    ax.spines["top"].set_linewidth(2);
    ax.spines["bottom"].set_linewidth(2);
    ax.spines["left"].set_linewidth(2);
    ax.spines["right"].set_linewidth(2);
  end
  fig.canvas.draw() # Update the figure
  #suptitle("All epochs")

  fig.text(0.5, 0.04, "East (mas)", ha="center", va="center",fontweight="bold", fontsize=15)
  fig.text(0.06, 0.5, "North (mas)", ha="center", va="center", rotation="vertical", fontweight="bold", fontsize=15)

  #fig.subplots_adjust(right=0.8)
  cbar_ax = fig.add_axes([0.92, 0.15, 0.05, 0.7])
  sm = plt.cm.ScalarMappable(cmap=colormap, norm=plt.Normalize(vmin=minT, vmax=maxT))
  sm.set_array([])
  plt.colorbar(sm,cax=cbar_ax,orientation="vertical",shrink=0.6)
  #plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9);
  #cax = plt.axes([0.85, 0.1, 0.075, 0.8])
  #plt.colorbar(cax=cax)
end


# use if star_map array has multiple individual maps within it
function plot2d_intensity_allepochs(star_map, star_geometry; plotmesh=false, tepochs = [], colormap="gist_heat",arr_box=23)
  patches = pyimport("matplotlib.patches")
  fig = figure("Intensity map -- All epochs",figsize=(15,10),facecolor="White")
  if plotmesh == true
    meshcolor = "grey"
  else
    meshcolor = "none"
  end

  minT = minimum(star_map);
  maxT = maximum(star_map);

  for t=1:length(star_geometry)
    ax= subplot(arr_box*10+t) # Create the 1st axis of a 2x2 arrax of axes
    if tepochs !=[]
      title("Epoch $t $(tepochs[t])",fontweight="bold") # Give the most recent axis a title
    end
    axis("equal")
    axis_max = maximum(sqrt.(star_geometry[t].vertices_xyz[:,:,1].^2 .+ star_geometry[t].vertices_xyz[:,:,2].^2 .+ star_geometry[t].vertices_xyz[:,:,3].^2))*1.5;
    ax.set_xlim([axis_max,-axis_max])
    ax.set_ylim([-axis_max,axis_max])
    visible_pixels = sometimes_visible(star_geometry);
    projmap = ((star_map[star_geometry[t].index_quads_visible].-minT)./(maxT-minT).*star_geometry[t].ldmap[star_geometry[t].index_quads_visible]);
    for i=1:star_geometry[t].nquads_visible
      p = patches.Polygon(hcat(-star_geometry[t].projx[i,:],star_geometry[t].projy[i,:]),
      closed=true,edgecolor=meshcolor,facecolor=get_cmap(colormap)(projmap[i]),fill="true",rasterized=false)
      ax.add_patch(p);
    end
    ax.plot();
    ax.tick_params(axis="both", which="both", labelsize=15, width=2, length=10);
    if (ceil(axis_max) <= 3.0)
      long_tick = 1.0;
      short_tick = 0.1;
    elseif ((ceil(axis_max) > 3.0) & (ceil(axis_max) <= 6.0))
      long_tick = 2.0;
      short_tick = 0.2;
    elseif (ceil(axis_max) > 6.0)
      long_tick = 3.0;
      short_tick = 0.5;
    end
    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(long_tick));
    ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(long_tick));
    ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(short_tick));
    ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(short_tick));
    ax.tick_params(axis="both", which="minor", width=2, length=5);
    # make axes thicker
    ax.spines["top"].set_linewidth(2);
    ax.spines["bottom"].set_linewidth(2);
    ax.spines["left"].set_linewidth(2);
    ax.spines["right"].set_linewidth(2);
  end
  fig.canvas.draw() # Update the figure
  #suptitle("All epochs")

  fig.text(0.5, 0.04, "East (mas)", ha="center", va="center",fontweight="bold", fontsize=15)
  fig.text(0.06, 0.5, "North (mas)", ha="center", va="center", rotation="vertical", fontweight="bold", fontsize=15)
end


function plot2d_intensity_allepochs_cmap(star_map, star_geometry; plotmesh=false, tepochs = [], minmaxT = [], colormap="gist_heat",arr_box=23)
  patches = pyimport("matplotlib.patches")
  fig = figure("Intensity map -- All epochs",figsize=(15,10),facecolor="White")
  if plotmesh == true
    meshcolor = "grey"
  else
    meshcolor = "none"
  end

  if (minmaxT == [])
    minT = minimum(star_map);
    maxT = maximum(star_map);
  else
    minT = minmaxT[1];
    maxT = minmaxT[2];
  end

  vmax = deepcopy(minT);
  vmin = deepcopy(maxT);

  for t=1:length(star_geometry)
    ax= subplot(arr_box*10+t) # Create the 1st axis of a 2x2 arrax of axes
    if tepochs !=[]
      epoch_title = string(tepochs[t])[1:10];
      title("Epoch $t "*epoch_title,fontweight="bold") # Give the most recent axis a title
    end
    axis_max = maximum(sqrt.(star_geometry[t].vertices_xyz[:,:,1].^2 .+ star_geometry[t].vertices_xyz[:,:,2].^2 .+ star_geometry[t].vertices_xyz[:,:,3].^2))*1.5;
    ax.set_xlim([axis_max,-axis_max])
    ax.set_ylim([-axis_max,axis_max])
    visible_pixels = sometimes_visible(star_geometry);
    projmap = ((star_map[star_geometry[t].index_quads_visible].-minT)./(maxT-minT).*star_geometry[t].ldmap[star_geometry[t].index_quads_visible]);

    if (maximum(star_map[star_geometry[t].index_quads_visible].*star_geometry[t].ldmap[star_geometry[t].index_quads_visible]) > vmax)
      vmax = maximum(star_map[star_geometry[t].index_quads_visible].*star_geometry[t].ldmap[star_geometry[t].index_quads_visible]);
    end

    if (minimum(star_map[star_geometry[t].index_quads_visible].*star_geometry[t].ldmap[star_geometry[t].index_quads_visible]) < vmin)
      vmin = minimum(star_map[star_geometry[t].index_quads_visible].*star_geometry[t].ldmap[star_geometry[t].index_quads_visible]);
    end

    for i=1:star_geometry[t].nquads_visible
      p = patches.Polygon(hcat(-star_geometry[t].projx[i,:],star_geometry[t].projy[i,:]),
      closed=true,edgecolor=meshcolor,facecolor=get_cmap(colormap)(projmap[i]),fill="true",rasterized=false)
      ax.add_patch(p);
      #println(star_map[star_geometry[t].index_quads_visible][i]*star_geometry[t].ldmap[star_geometry[t].index_quads_visible][i])
    end
    ax.plot();
    ax.tick_params(axis="both", which="both", labelsize=15, width=2, length=10);
    if (ceil(axis_max) <= 3.0)
      long_tick = 1.0;
      short_tick = 0.1;
    elseif ((ceil(axis_max) > 3.0) & (ceil(axis_max) <= 6.0))
      long_tick = 2.0;
      short_tick = 0.2;
    elseif (ceil(axis_max) > 6.0)
      long_tick = 3.0;
      short_tick = 0.5;
    end
    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(long_tick));
    ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(long_tick));
    ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(short_tick));
    ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(short_tick));
    ax.tick_params(axis="both", which="minor", width=2, length=5);
    # make axes thicker
    ax.spines["top"].set_linewidth(2);
    ax.spines["bottom"].set_linewidth(2);
    ax.spines["left"].set_linewidth(2);
    ax.spines["right"].set_linewidth(2);
  end
  fig.canvas.draw() # Update the figure
  #suptitle("All epochs")

  cbar_ax = fig.add_axes([0.92, 0.15, 0.05, 0.7])
  sm = plt.cm.ScalarMappable(cmap=colormap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
  sm.set_array([])
  plt.colorbar(sm,cax=cbar_ax,orientation="vertical",shrink=0.6)

  fig.text(0.5, 0.04, "East (mas)", ha="center", va="center",fontweight="bold", fontsize=15)
  fig.text(0.06, 0.5, "North (mas)", ha="center", va="center", rotation="vertical", fontweight="bold", fontsize=15)
end



function mollplot_temperature_healpix(image; visible_pixels = [], vmin = -Inf, vmax = Inf, incl=90.0, colormap="gist_heat", figtitle="Mollweide")
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
  if incl != 90.0
    ax.axhline(-incl * pi/180, c=:black, ls="-.")
  end
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
cb.ax.xaxis.set_label_text("Temperature (K)")

# workaround for issue with viewers, see colorbar docstring
cb.solids.set_edgecolor("face")
ax.tick_params(axis="x", labelsize=15)
ax.tick_params(axis="y", labelsize=15)
end

function mollplot_temperature_longlat(image, ntheta, nphi; visible_pixels = [], vmin = -Inf, vmax = Inf, colormap="gist_heat", figtitle="Mollweide", incl=90.0)
  xsize = 2000
  ysize = div(xsize,2)

  ## DJ: had to set xsize=ntheta and ysize=nphi to get plotting to work for large
  ##  number of latitudes and longitudes
  # xsize = ntheta
  # ysize = nphi

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
  # if visible_pixels == []
  #     vmin = minimum(image);
  #     vmax = maximum(image)*1.01;
  # else
  #    vmin = minimum(image[visible_pixels]);
  #    vmax = maximum(image[visible_pixels])*1.01;
  # end
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
end

function mollplot_temperature_longlat_visblock(image, ntheta, nphi, star_epoch_geom; vmin = -1000.0, vmax = 1000.0, colormap="gist_heat", figtitle="Mollweide")
  xsize = 2000
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

  pixels_vis = Array{Bool}(undef, star_epoch_geom.npix);
  pixels_hidden = Array{Bool}(undef, star_epoch_geom.npix);
  pixels_vis[:] .= false;
  pixels_hidden[:] .= true;
  pixels_vis[star_epoch_geom.index_quads_visible] .= true;
  pixels_hidden[star_epoch_geom.index_quads_visible] .= false;

  grid_map_vis = (image.*pixels_vis)[circshift(grid_pix,(0,Int(xsize/2)))];
  grid_map_hidden = (image.*pixels_hidden)[circshift(grid_pix,(0,Int(xsize/2)))];
  fig = figure(figtitle, figsize=(10, 7))
  clf();
  ax = subplot(111,projection="mollweide",title=title)
  # rasterized makes the map bitmap while the labels remain vectorial
  # flip longitude to the astro convention
  # if visible_pixels == []
  #     vmin = minimum(image);
  #     vmax = maximum(image)*1.01;
  # else
  #    vmin = minimum(image[visible_pixels]);
  #    vmax = maximum(image[visible_pixels])*1.01;
  # end
  vmin = minimum(image[star_epoch_geom.index_quads_visible]);
  vmax = maximum(image[star_epoch_geom.index_quads_visible]);

  moll = pcolormesh(longitude, latitude, grid_map_vis, vmin=vmin, vmax=vmax, rasterized=true, cmap=colormap)
  moll2 = pcolormesh(longitude, latitude, grid_map_hidden, rasterized=true, cmap="spring_r")
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
  cb.ax.xaxis.set_label_text("Temperature (K)")
  # workaround for issue with viewers, see colorbar docstring
  cb.solids.set_edgecolor("face")
  ax.tick_params(axis="x", labelsize=15)
  ax.tick_params(axis="y", labelsize=15)
end


function plot2d_temperature_binary(star_map1, star_map2, star_epoch_geom1, star_epoch_geom2; zorder = 1, plotmesh=false,labels=false, poleline_frac = 0.35, draw_orbit_arrow=true, draw_convex_hull=true, draw_graticules=true)
    patches = pyimport("matplotlib.patches")
    collections = pyimport("matplotlib.collections")
    fig = figure("Epoch image",figsize=(10,10),facecolor="White")
    
    if plotmesh == true
      meshcolor = "grey"
    else
      meshcolor = "none"
    end
    xlabel("East (mas)", fontweight="bold", fontsize=15);
    ylabel("West (mas)", fontweight="bold", fontsize=15);

    # find what the axes are
    star1_radiusmax = maximum(sqrt.(star_epoch_geom1.vertices_xyz[:,:,1].^2 .+ star_epoch_geom1.vertices_xyz[:,:,2].^2 .+ star_epoch_geom1.vertices_xyz[:,:,3].^2));
    star2_radiusmax = maximum(sqrt.(star_epoch_geom2.vertices_xyz[:,:,1].^2 .+ star_epoch_geom2.vertices_xyz[:,:,2].^2 .+ star_epoch_geom2.vertices_xyz[:,:,3].^2));
    axis_max = (star1_radiusmax + star2_radiusmax)*1.5;
    
    axis_max = 3.0
    ax = gca();
    ax.set_aspect("equal")
    ax.set_xlim([axis_max,-axis_max])
    ax.set_ylim([-axis_max,axis_max])


    projmap1 = star_map1[star_epoch_geom1.index_quads_visible];
    projmap2 = star_map2[star_epoch_geom2.index_quads_visible];
    max_temperature = maximum([projmap1;projmap2])*1.1; # multiply maximum temperature for best color scheme

    if (zorder == 1) # Primary is in front
      for i=1:star_epoch_geom1.nquads_visible
          # star 1
          p = patches.Polygon(hcat(star_epoch_geom1.projx[i,:], star_epoch_geom1.projy[i,:]),
          closed=true,edgecolor=meshcolor,facecolor=get_cmap("gist_heat")(projmap1[i]/max_temperature),fill=true,rasterized=false, zorder=2);
          ax.add_patch(p);
      end
      for i=1:star_epoch_geom2.nquads_visible
          # star 2
          p = patches.Polygon(hcat(star_epoch_geom2.projx[i,:], star_epoch_geom2.projy[i,:]),
          closed=true,edgecolor=meshcolor,facecolor=get_cmap("gist_heat")(projmap2[i]/max_temperature),fill=true,rasterized=false, zorder=3)
          ax.add_patch(p);
      end
      else
      for i=1:star_epoch_geom2.nquads_visible
          # star 2
          p = patches.Polygon(hcat(star_epoch_geom2.projx[i,:], star_epoch_geom2.projy[i,:]),
          closed=true,edgecolor=meshcolor,facecolor=get_cmap("gist_heat")(projmap2[i]/max_temperature),fill=true,rasterized=false, zorder=2)
          ax.add_patch(p);
      end
      for i=1:star_epoch_geom1.nquads_visible
          # star 1
          p = patches.Polygon(hcat(star_epoch_geom1.projx[i,:], star_epoch_geom1.projy[i,:]),
          closed=true,edgecolor=meshcolor,facecolor=get_cmap("gist_heat")(projmap1[i]/max_temperature),fill=true,rasterized=false, zorder=3);
          ax.add_patch(p);
      end
    end
    ax.tick_params(axis="both", which="both", labelsize=15, width=2, length=10);
    if (ceil(axis_max) <= 3.0)
    long_tick = 1.0;
    short_tick = 0.1;
    elseif ((ceil(axis_max) > 3.0) & (ceil(axis_max) <= 6.0))
    long_tick = 2.0;
    short_tick = 0.2;
    elseif (ceil(axis_max) > 6.0)
    long_tick = 3.0;
    short_tick = 0.5;
    end
    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(long_tick));
    ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(long_tick));
    ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(short_tick));
    ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(short_tick));
    ax.tick_params(axis="both", which="minor", width=2, length=5);
    # make axes thicker
    ax.spines["top"].set_linewidth(2);
    ax.spines["bottom"].set_linewidth(2);
    ax.spines["left"].set_linewidth(2);
    ax.spines["right"].set_linewidth(2);
end

# this plots the temperature map onto the projected 2D image plane (= observer view)
# Adds polelines. TBD: Needs more centering of polelines
function plot2d_temperature_binary_poleline(star_map1, star_map2, star_geometry1, star_geometry2, bparameters, tepoch; plotmesh=false,labels=false)
  rotation_period = bparameters.binary_period;
  fillout_factor = bparameters.fillout_factor;
  
  # still missing the actual intensity (includes LD)
  patches = pyimport("matplotlib.patches")
  fig = figure("Epoch image",figsize=(10,10),facecolor="White")
  ax = fig.add_axes([0.1,0.1,0.85,0.85])
  if plotmesh == true
    meshcolor = "grey"
  else
    meshcolor = "none"
  end
  xlabel("X (mas)", fontweight="bold", fontsize=15);
  ylabel("Y (mas)", fontweight="bold", fontsize=15);

  # find what the axes are
  star1_radiusmax = maximum(sqrt.(star_geometry1.vertices_xyz[:,:,1].^2 .+ star_geometry1.vertices_xyz[:,:,2].^2 .+ star_geometry1.vertices_xyz[:,:,3].^2));
  star2_radiusmax = maximum(sqrt.(star_geometry2.vertices_xyz[:,:,1].^2 .+ star_geometry2.vertices_xyz[:,:,2].^2 .+ star_geometry2.vertices_xyz[:,:,3].^2));
  axis_max = (star1_radiusmax + star2_radiusmax)*1.5;
  ax.set_xlim([axis_max,-axis_max])
  ax.set_ylim([-axis_max,axis_max])

  projmap1 = star_map1[star_geometry1.index_quads_visible];
  projmap2 = star_map2[star_geometry2.index_quads_visible];
  max_temperature = maximum([projmap1;projmap2])*1.1; # multiply maximum temperature for best color scheme

  # Adds arrows to plotting
  a1 = star_geometry2.projx[1,:][2]
  b1 = star_geometry2.projy[1,:][2]
  a2 = star_geometry2.projx[end,:][2]
  b2 = star_geometry2.projy[end,:][2]
  sx = a2-a1
  sy = b2-b1
  c = -sx*0.25
  d = -sy*0.25
  arrow(a1, b1, c, d,
  head_width=0.008,
  width=0.005,
  head_length=0.007,
  overhang=0.2,
  facecolor="black")
  arrow(a2, b2, -c, -d,
  head_width=0.008,
  width=0.005,
  head_length=0.000,
  overhang=0.0,
  facecolor="black")

  aa1 = star_geometry1.projx[1,:][2]
  bb1 = star_geometry1.projy[1,:][2]
  aa2 = star_geometry1.projx[end,:][2]
  bb2 = star_geometry1.projy[end,:][2]
  sxx = aa2-aa1
  syy = bb2-bb1
  cc = -sxx*.25
  dd = -syy*.25
  arrow(aa1, bb1, cc, dd,
  head_width=0.008,
  width=0.005,
  head_length=0.007,
  overhang=0.2,
  head_starts_at_zero="true",
  facecolor="black")

  arrow(aa2, bb2, -cc, -dd,
  head_width=0.008,
  width=0.005,
  head_length=0.000,
  overhang=0.0,
  head_starts_at_zero="true",
  facecolor="black")

  x1, y1, z1, x2, y2, z2 = binary_orbit_abs(bparameters,tepoch);
  if (z1 > z2) # plots furthest star first
    for i=1:star_geometry1.nquads_visible
      # star 1
      p = patches.Polygon(hcat(star_geometry1.projx[i,:], star_geometry1.projy[i,:]),
      closed=true,edgecolor=meshcolor,facecolor=get_cmap("gist_heat")(projmap1[i]/max_temperature),fill=true,rasterized=false);
      ax.add_patch(p);
    end
    for i=1:star_geometry2.nquads_visible
      # star 2
      p = patches.Polygon(hcat(star_geometry2.projx[i,:], star_geometry2.projy[i,:]),
      closed=true,edgecolor=meshcolor,facecolor=get_cmap("gist_heat")(projmap2[i]/max_temperature),fill=true,rasterized=false)
      ax.add_patch(p);
    end
  else
    for i=1:star_geometry2.nquads_visible
      # star 2
      p = patches.Polygon(hcat(star_geometry2.projx[i,:], star_geometry2.projy[i,:]),
      closed=true,edgecolor=meshcolor,facecolor=get_cmap("gist_heat")(projmap2[i]/max_temperature),fill=true,rasterized=false)
      ax.add_patch(p);
    end
    for i=1:star_geometry1.nquads_visible
      # star 1
      p = patches.Polygon(hcat(star_geometry1.projx[i,:], star_geometry1.projy[i,:]),
      closed=true,edgecolor=meshcolor,facecolor=get_cmap("gist_heat")(projmap1[i]/max_temperature),fill=true,rasterized=false);
      ax.add_patch(p);
    end
  end

  ax.tick_params(axis="both", which="both", labelsize=15, width=2, length=10);
  if (ceil(axis_max) <= 3.0)
    long_tick = 1.0;
    short_tick = 0.1;
  elseif ((ceil(axis_max) > 3.0) & (ceil(axis_max) <= 6.0))
    long_tick = 2.0;
    short_tick = 0.2;
  elseif (ceil(axis_max) > 6.0)
    long_tick = 3.0;
    short_tick = 0.5;
  end
  ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(long_tick));
  ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(long_tick));
  ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(short_tick));
  ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(short_tick));
  ax.tick_params(axis="both", which="minor", width=2, length=5);
  # make axes thicker
  ax.spines["top"].set_linewidth(2);
  ax.spines["bottom"].set_linewidth(2);
  ax.spines["left"].set_linewidth(2);
  ax.spines["right"].set_linewidth(2);

  ax.plot();
  PyPlot.draw()

end


# #note: this requires healpy
# function healpyplot(x; vmin = -1e9, vmax = 1e9)
#   if vmin==-1e9
#     vmin = minimum(x)
#   end
#   if vmax==1e9
#     vmax = maximum(x)
#   end
# hp=pyimport("healpy.visufunc")
# hp.mollview(x, nest=true, fig=1, min=vmin, max=vmax, rot=0, flip="geo", cbar=false, cmap="gist_heat", title="")
# hp.graticule(dpar=10, dmer=10,  force=true, verbose=false, lw=0.25)
# hp.cartview(x, nest=true, fig=2, min=vmin, max=vmax, rot=0, flip="geo", cbar=false, cmap="gist_heat", title="")
# hp.graticule(dpar=10, dmer=10,  force=true, verbose=false, lw=0.25)
# hp.orthview(x, nest=true, fig=3, min=vmin, max=vmax, rot=0, flip="geo", cbar=false, cmap="gist_heat", title="")
# hp.graticule(dpar=10, dmer=10,  force=true, verbose=false, lw=0.25)
# end
