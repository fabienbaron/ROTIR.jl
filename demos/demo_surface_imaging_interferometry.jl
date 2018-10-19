using ROTIR
oifitsfiles = ["./data/2011Sep02.lam_And_prepped.oifits", "./data/2011Sep06.lam_And_prepped.oifits",
"./data/2011Sep10.lam_And_prepped.oifits","./data/2011Sep14.lam_And_prepped.oifits",
"./data/2011Sep19.lam_And_prepped.oifits","./data/2011Sep24.lam_And_prepped.oifits"];

# Import OIFITS
nepochs, tepochs, data = readoifits_multiepochs(oifitsfiles)#, redundance_chk=true);
#
# tot_nv2 = sum([data[i].nv2 for i=1:nepochs])
# tot_nt3amp = sum([data[i].nt3amp for i=1:nepochs]);
# tot_nt3phi = sum([data[i].nt3phi for i=1:nepochs]);
# tot_pts = tot_nv2 + tot_nt3amp + tot_nt3phi;

# Define stellar parameters for each epoch
stellar_parameters = Array{starparameters}(undef, nepochs);
# starparams = [1.37131,  # milliarcseconds (radius at pole)
#     4800.,              # Kelvin (at pole)
#     0.,                 # unitless; fractional rotational velocity
#     [3,0.22886],        # limb darkening,first coefficient is for LD law type, then LD coefficients
#     0.08,               # exponent for von Zeipel law
#     -1.43075,           # 2nd constant for rotational velocity
#     71.0962,            # degrees; inclination
#     18.6033,            # degrees; position_angle
#     51.8513];           # days; rotation_period

starparams = [1.3795,4800.,0.,[3,0.22886],0.08,-1.43075,78.0,-21.5,54.0]; # Rob Park's parameters (used for VR)

for i=1:nepochs
    stellar_parameters[i]=starparameters(starparams[1],starparams[2],starparams[3],starparams[4],starparams[5],
        starparams[6],starparams[7],starparams[8],360.0/starparams[9]*(tepochs[i]-tepochs[1]),starparams[9]);
end

# Create 3D geometry from parameters

# Healpix
# n=3;star_epoch_geom = create_geometry( healpix_round_star(n,stellar_parameters[1].radius), stellar_parameters);

# Longitude/latitude
const ntheta = 20;#40;
const nphi = 40;#80;
a,b,c = oblate_const(stellar_parameters[1]);
#star_epoch_geom = create_geometry(latlong_ellipsoid_star(ntheta,nphi, a, b, c), stellar_parameters);
star_epoch_geom = create_geometry(latlong_ellipsoid_star(ntheta,nphi, a, b, c),stellar_parameters,diff_rot=false, ntheta=ntheta, nphi=nphi);


# Might want to have a small iteration temperature map and then
# optimize parameters. Try to see if you can do small iterations on NLopt

# Initial temperature map
temperature_map_start = calc_tempmap_vZ(stellar_parameters[1],star_epoch_geom[1]);#4800*ones(star_epoch_geom[1].npix);
# Setup Fourier transform (polygons -> complex visibilities)
polyflux, polyft = setup_polygon_ft(data, star_epoch_geom);
# Define epoch weights and function to minimize
epochs_weights = ones(Float64, nepochs)/nepochs;

# Define neighbors for healpix or longlat
#neighbors,south_neighbors,west_neighbors,south_neighbors_reverse,west_neighbors_reverse = tv_neighbours_healpix(n);
neighbors,south_neighbors,west_neighbors,south_neighbors_reverse,west_neighbors_reverse = tv_neighbors_longlat(ntheta,nphi);

# Define regularization weight
mu = 0.01;

# Define imaging criterion (will be minimized)
#crit_imaging = (x,g)->chi2_allepochs_fg(x, g, epochs_weights, polyflux, polyft, data); # chi red = 6300 for 20 iterations INVESTIGATE WHERE SLOWS DOWN
crit_imaging = (x,g)->spheroid_chi2_allepochs_tv_fg(x, g, mu, south_neighbors,west_neighbors, south_neighbors_reverse, west_neighbors_reverse, epochs_weights, polyflux, polyft, data);
#crit_imaging = (x,g)->chi2_allepochs_l2_fg(x, g, mu, epochs_weights, polyflux, polyft, data); # chi red = for 20 iterations
# Optimize criterion
using OptimPackNextGen
temperature_map = vmlmb(crit_imaging, temperature_map_start, verb=true, lower=0, maxiter=150, blmvm=false);

# Plot results
plot2d_temperature_allepochs(temperature_map, star_epoch_geom, tepochs = tepochs);
plot2d_intensity_allepochs(temperature_map, star_epoch_geom, tepochs = tepochs);

temperature_map_withhiddenblack=copy(temperature_map);
temperature_map_withhiddenblack[never_visible(star_epoch_geom)].=0; # or -1.6375e+30; # this constant is the Healpix constant for hidden vals
plot3d_temperature(temperature_map_withhiddenblack,star_epoch_geom[1]);

mollplot_temp_longlat(temperature_map_withhiddenblack, ntheta, nphi); # long/lat pix scheme
#mollplot_temperature(temperature_map_withhiddenblack, false); # for Healpix


# save figures and make movie
file_loc = "./lam_And_demo"
for i = 1:nepochs
    plot2d_temperature_savefig(temperature_map,star_epoch_geom[i],file_loc=file_loc,iteration=i,plot_lim=[-2,2],plotmesh=false,color_map="gist_heat");
end

run(`ffmpeg -f image2 -r 6 -pattern_type glob -i 'lam_And_demo*.png' lam_And_demo_all.mp4`) # make a movie
run(`ffmpeg -i lam_And_demo_all.mp4 -filter:v "setpts=5*PTS" lam_And_demo_all2.mp4`) # make video longer

using NLopt
function crit_modeling(params::Vector, dummy::Vector)
#  radius = x[1]
#  inclination = x[2]
#  position_angle = x[3]
#  rotation_period = x[4]
#  limbdarkening = x[5]
  println(params);
  for i=1:nepochs
      stellar_parameters[i]=starparameters(params[1],params[2],params[3],params[4],params[5],360./params[6]*(tepochs[i]-tepochs[1]),params[6],params[7],params[8],params[9],params[10]);
  end
  ntheta = 20; nphi = 40; # For lambda And
  a,b,c = oblate_const(stellar_parameters[1]);
  star_epoch_geom = create_geometry(latlong_ellipsoid_star(ntheta,nphi, a, b, c),stellar_parameters,diff_rot=true, ntheta=ntheta, nphi=nphi);
  polyflux_new, polyft_new = setup_polygon_ft(data, star_epoch_geom);
  f = chi2_allepochs_f(temperature_map, epochs_weights, polyflux_new, polyft_new, data);
  return f
end

opt = Opt(:LN_NELDERMEAD, 7)
#lower_bounds!(opt, [1.26798, 65., 10., 45., 0.1]) #[1.36798,70.2921,19.3993,52.7254,0.21786]; # chi of 6208.91577961831 for 2011 dataset
#upper_bounds!(opt, [1.46798, 85., 30., 75., 0.4])
lower_bounds!(opt, [1.26798, 67., 16., 49., 0.15, -2.5, -3]) # needs updating
upper_bounds!(opt, [1.46798, 73., 22., 55., 0.31, 3, 3]) # needs updating
#xtol_rel!(opt,1e-3)
tol = [1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3];
xtol_abs!(opt,1e-3)
min_objective!(opt, crit_modeling)
(minchi2,params,ret) = optimize(opt, params)
#println("got $minchi2 at $params_opt (returned $ret)")
println("got $minchi2 at $params (returned $ret)")


# To save reconstruction
using HDF5
h5open("reconstruction.h5", "w") do file
  write(file, "temperature", temperature_map)
  write(file, "starparams", params)
end

# To resume
temperature_map = h5read("reconstruction.h5", "temperature")
params = h5read("reconstruction.h5", "starparams")


# Save parameters for ParaView
using WriteVTK
npix = ntheta*nphi;
nxyz = length(star_epoch_geom[1].vertices_xyz[:,1,1:4]);
x_val = reshape(star_epoch_geom[1].vertices_xyz[:,1,1:4]',1,nxyz);
y_val = reshape(star_epoch_geom[1].vertices_xyz[:,2,1:4]',1,nxyz);
z_val = reshape(star_epoch_geom[1].vertices_xyz[:,3,1:4]',1,nxyz);
lam_Andpos = Array{Float64}(3,nxyz);
lam_Andpos[1,:] = x_val; lam_Andpos[2,:] = y_val; lam_Andpos[3,:] = z_val;

celltype = VTKCellTypes.VTK_QUAD;
cells = MeshCell[];

for i=1:npix
    # Define connectivity of cell.
    indx = collect(((i-1)*4+1):(i*4));

    # Define cell.
    cell = MeshCell(celltype, indx);

    push!(cells, cell);
end

vtkfile = vtk_grid("lam_And",lam_Andpos,cells);

vtk_cell_data(vtkfile,temperature_map,"Temperature");
outfiles = vtk_save(vtkfile);


# Time series
pvd = paraview_collection("multi_epoch_lamAnd");
# make vtk setup
for i=1:nepochs
    x_val = reshape(star_epoch_geom[i].vertices_xyz[:,1,1:4]',1,nxyz);
    y_val = reshape(star_epoch_geom[i].vertices_xyz[:,2,1:4]',1,nxyz);
    z_val = reshape(star_epoch_geom[i].vertices_xyz[:,3,1:4]',1,nxyz);
    lam_Andpos = Array{Float64}(3,nxyz);
    lam_Andpos[1,:] = x_val; lam_Andpos[2,:] = y_val; lam_Andpos[3,:] = z_val;
    #println(lam_Andpos);
    vtkfile = vtk_grid("lam_And_epoch"*string(i),lam_Andpos,cells);
    vtk_cell_data(vtkfile,temperature_map,"Temperature");
    epoch_time = tepochs[i]-tepochs[1];
    collection_add_timestep(pvd, vtkfile, epoch_time);
end
vtk_save(pvd)


# More epochs
new_tepochs = collect(linspace(0.,tepochs[6]-tepochs[1],22)); new_nepochs = length(new_tepochs);
new_stellar_parameters = Array{starparameters}(new_nepochs);
for i=1:new_nepochs
    new_stellar_parameters[i]=starparameters(starparams[1],starparams[2],starparams[3],starparams[4],starparams[5],
        starparams[6],starparams[7],starparams[8],360./starparams[9]*(new_tepochs[i]-new_tepochs[1]),starparams[9]);
end
new_star_epoch_geom = create_geometry(latlong_ellipsoid_star(ntheta,nphi, a, b, c),new_stellar_parameters,diff_rot=false, ntheta=ntheta, nphi=nphi);

# Time series
pvd = paraview_collection("multi_epoch_lamAnd_finertime");
# make vtk setup
for i=1:new_nepochs
    x_val = reshape(new_star_epoch_geom[i].vertices_xyz[:,1,1:4]',1,nxyz);
    y_val = reshape(new_star_epoch_geom[i].vertices_xyz[:,2,1:4]',1,nxyz);
    z_val = reshape(new_star_epoch_geom[i].vertices_xyz[:,3,1:4]',1,nxyz);
    lam_Andpos = Array{Float64}(3,nxyz);
    lam_Andpos[1,:] = x_val; lam_Andpos[2,:] = y_val; lam_Andpos[3,:] = z_val;
    #println(lam_Andpos);
    vtkfile = vtk_grid("lam_And_epoch_finer"*string(i),lam_Andpos,cells);
    vtk_cell_data(vtkfile,temperature_map,"Temperature");
    epoch_time = new_tepochs[i]-new_tepochs[1];
    collection_add_timestep(pvd, vtkfile, epoch_time);
end
vtk_save(pvd)
