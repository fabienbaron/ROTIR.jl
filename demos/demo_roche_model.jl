#include("./rotir_mod.jl")
#using rotir_module
include("./lib/geometry.jl");
include("./lib/readoifits.jl");
include("./lib/oichi2_spheroid.jl");
include("./lib/oiplot.jl");
include("./lib/oistars.jl");
include("./lib/geometry_rochelobe.jl");

"""
Code still in progress
TO DO:
1. Add restrictions for fill factor and async ratio -- DONE
2. Further testing on elliptical orbits (with final testing with different fill factors)
3. Make update_binary function work for this occasion 
"""
# using OptimPack
# oifitsfiles = ["./DATA/2011Sep02.lam_And_prepped.oifits", "./DATA/2011Sep06.lam_And_prepped.oifits",
# "./DATA/2011Sep10.lam_And_prepped.oifits","./DATA/2011Sep14.lam_And_prepped.oifits",
# "./DATA/2011Sep19.lam_And_prepped.oifits","./DATA/2011Sep24.lam_And_prepped.oifits"];
# nepochs, tepochs, data = readoifits_multiepochs(oifitsfiles)#, redundance_chk=true);

nstars = 2; stellar_parameters = Array{starparameters}(nstars);
nepochs = 12; binary_parameters = Array{binaryparameters}(nepochs);

# parameters for beta Lyrae
# parameters from Zhao et al. 2008
# http://iopscience.iop.org/article/10.1086/592146/pdf
# and Mennickent and Djurašević 2013
# https://academic.oup.com/mnras/article/432/1/799/1142291

starparams = Array{Any}(nstars,9);
starparams[1,:] = [0.239, 13800., 0., [3,0.22886], 0.25, -1.43075, 0., 0., 12.9414];
starparams[2,:] = [0.094, 30200., 0., [3,0.22886], 0.25, 0., 0., 0., 12.9414];
# starparams[1,:] = [0.239, 5000., 0., [3,0.22886], 0.25, -1.43075, 0., 0., 12.9414];
# starparams[2,:] = [0.094, 10000., 0., [3,0.22886], 0.25, 0., 0., 0., 12.9414];
# MUST ADD RESTRICTIONS IN FUNCTIONS FOR FILL FACTOR AND ASYNC RATIO
binaryparams = [91.96, 251.87, 45., 12.9414, 0., 0.865, 0., 0.226, 1., 0.6005975]; 
# 0.993 # 0.6005975

tepochs = collect(linspace(0.,binaryparams[4],Int(ceil(binaryparams[4]))));

# only works for 2 stars, for now
for i=1:nepochs
    for j=1:nstars
        stellar_parameters[j]=starparameters(starparams[j,1],starparams[j,2],starparams[j,3],starparams[j,4],starparams[j,5],
            starparams[j,6],starparams[j,7],starparams[j,8],360./starparams[j,9]*(tepochs[i]-tepochs[1]),starparams[j,9]);
    end
    binary_parameters[i] = binaryparameters(stellar_parameters[1],stellar_parameters[2],binaryparams[1],binaryparams[2],binaryparams[3],
        binaryparams[4],binaryparams[5],binaryparams[6],binaryparams[7],binaryparams[8],binaryparams[9],binaryparams[10]);
end

const ntheta = 40;
const nphi = 80;
star_epoch_geom1 = Array{Any}(nepochs);
star_epoch_geom2 = Array{Any}(nepochs);


# star 1
star_base_geom1 = update_star(latlong_rochelobe(ntheta,nphi,binary_parameters[1]),binary_parameters[1].star1, ntheta=ntheta, nphi=nphi);
# star 2 (binary flag for secondary)
star_base_geom2 = update_star(latlong_rochelobe(ntheta,nphi,binary_parameters[1];secondary=true),binary_parameters[1].star2, ntheta=ntheta, 
    nphi=nphi,secondary=true,bparameters=binary_parameters[1]);

# for multiple epochs
for i=1:nepochs
    star_epoch_geom1[i] = update_star(latlong_rochelobe(ntheta,nphi,binary_parameters[i]),binary_parameters[i].star1, ntheta=ntheta, nphi=nphi)
    star_epoch_geom2[i] = update_star(latlong_rochelobe(ntheta,nphi,binary_parameters[i];secondary=true),binary_parameters[i].star2, ntheta=ntheta, 
        nphi=nphi,secondary=true,bparameters=binary_parameters[i])
end


# star 1 temperature map
temperature_map1 = compute_teff_vonzeipel(binary_parameters[1].star1.radius,binary_parameters[1].star1.temperature,star_base_geom1.vertices_spherical[:,1,5],
    binary_parameters[1].separation,star_base_geom1.vertices_spherical[:,2,5],star_base_geom1.vertices_spherical[:,3,5],binary_parameters[1].async_ratio,
    binary_parameters[1].mass_ratio,binary_parameters[1].star1.beta);
# star 2 temperature map
temperature_map2 = compute_teff_vonzeipel(binary_parameters[1].star2.radius,binary_parameters[1].star2.temperature,star_base_geom2.vertices_spherical[:,1,5],
    binary_parameters[1].separation,star_base_geom2.vertices_spherical[:,2,5],star_base_geom2.vertices_spherical[:,3,5],binary_parameters[1].async_ratio,
    binary_parameters[1].mass_ratio,binary_parameters[1].star2.beta);
#plot3d_temperature(temperature_map1,latlong_rochelobe(ntheta,nphi,binary_parameters[1].star1.radius,binary_parameters[1]));
#plot3d_temperature(temperature_map2,star_base_geom2);
plot3d_temperature_binary(temperature_map1,temperature_map2,star_base_geom1,star_base_geom2);


# Save parameters for ParaView
using WriteVTK
npix = ntheta*nphi;
nxyz = length(star_base_geom1.vertices_xyz[:,1,1:4]);
# Primary
x_val1 = reshape(star_base_geom1.vertices_xyz[:,1,1:4]',1,nxyz);
y_val1 = reshape(star_base_geom1.vertices_xyz[:,2,1:4]',1,nxyz);
z_val1 = reshape(star_base_geom1.vertices_xyz[:,3,1:4]',1,nxyz);
primary_pos = Array{Float64}(3,nxyz);
primary_pos[1,:] = x_val1; primary_pos[2,:] = y_val1; primary_pos[3,:] = z_val1;

# Secondary
x_val2 = reshape(star_base_geom2.vertices_xyz[:,1,1:4]',1,nxyz);
y_val2 = reshape(star_base_geom2.vertices_xyz[:,2,1:4]',1,nxyz);
z_val2 = reshape(star_base_geom2.vertices_xyz[:,3,1:4]',1,nxyz);
secondary_pos = Array{Float64}(3,nxyz);
secondary_pos[1,:] = x_val2; secondary_pos[2,:] = y_val2; secondary_pos[3,:] = z_val2;

celltype = VTKCellTypes.VTK_QUAD;
cells = MeshCell[];

for i=1:npix
    # Define connectivity of cell.
    indx = collect(((i-1)*4+1):(i*4));

    # Define cell.
    cell = MeshCell(celltype, indx);

    push!(cells, cell);
end

vtkfile1 = vtk_grid("binary_primaryff1_tempalt",primary_pos,cells);
vtk_cell_data(vtkfile1,temperature_map1,"Temperature");
outfiles1 = vtk_save(vtkfile1);

vtkfile2 = vtk_grid("binary_secondaryff1_tempalt",secondary_pos,cells);
vtk_cell_data(vtkfile2,temperature_map2,"Temperature");
outfiles2 = vtk_save(vtkfile2);


# make binary as one system
binary_pos = hcat(primary_pos,secondary_pos);
temperature_map_comb = vcat(temperature_map1,temperature_map2);

# Redefine cells
cells = MeshCell[];

for i=1:npix*2
    # Define connectivity of cell.
    indx = collect(((i-1)*4+1):(i*4));

    # Define cell.
    cell = MeshCell(celltype, indx);

    push!(cells, cell);
end

vtkfile = vtk_grid("binary_ff1_tempalt",binary_pos,cells); # save cells to vtk file
vtk_cell_data(vtkfile,temperature_map_comb,"Temperature"); # save temperature data to cells
outfiles = vtk_save(vtkfile); # save data file


pvd = paraview_collection("binary_orbit_1_tempalt");
for i = 1:nepochs
    x_val_temp1 = reshape(star_epoch_geom1[i].vertices_xyz[:,1,1:4]',1,nxyz);
    y_val_temp1 = reshape(star_epoch_geom1[i].vertices_xyz[:,2,1:4]',1,nxyz);
    z_val_temp1 = reshape(star_epoch_geom1[i].vertices_xyz[:,3,1:4]',1,nxyz);

    star1_pos = Array{Float64}(3,nxyz);
    star1_pos[1,:] = x_val_temp1; star1_pos[2,:] = y_val_temp1; star1_pos[3,:] = z_val_temp1;

    x_val_temp2 = reshape(star_epoch_geom2[i].vertices_xyz[:,1,1:4]',1,nxyz);
    y_val_temp2 = reshape(star_epoch_geom2[i].vertices_xyz[:,2,1:4]',1,nxyz);
    z_val_temp2 = reshape(star_epoch_geom2[i].vertices_xyz[:,3,1:4]',1,nxyz);

    star2_pos = Array{Float64}(3,nxyz);
    star2_pos[1,:] = x_val_temp2; star2_pos[2,:] = y_val_temp2; star2_pos[3,:] = z_val_temp2;

    binary_pos_multepoch = hcat(star1_pos,star2_pos);

    vtkfile = vtk_grid("binary_orbit_1_tempalt"*string(i),binary_pos_multepoch,cells);

    vtk_cell_data(vtkfile,temperature_map_comb,"Temperature");
    collection_add_timestep(pvd, vtkfile, tepochs[i]);
end

vtk_save(pvd)