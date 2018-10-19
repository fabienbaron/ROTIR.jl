include("./lib/geometry.jl");
include("./lib/oiplot.jl");
include("./lib/oistars.jl")

# make vectors
nepochs = 600; nepochs_2 = 2*nepochs; tot_nepochs = 2*nepochs+nepochs_2;
ld_grid = collect(linspace(0.0,0.0,nepochs)); # no limb darkening, no rapid rotation
ld_grid = vcat(ld_grid,collect(linspace(0.2,0.2,nepochs_2))); # introduce rapid rotation
ld_grid = vcat(ld_grid,collect(linspace(0.2,0.2,nepochs))); # have rapid rotation
oblate = collect(linspace(1.e-5,1.e-5,nepochs)); # no limb darkening, no rapid rotation
oblate = vcat(oblate,collect(logspace(-5,-0.6057970289817278,nepochs_2))); # introduce rapid rotation
oblate = vcat(oblate,collect(linspace(0.24785801713586308,0.24785801713586308,nepochs))); # have rapid rotation
omega, R_pole, R_equ = calc_omega(1.634,oblate); # true omega for Atlair is 0.923
rotational_vel, rotation_period = calc_rotspin(1.634,R_equ,omega,1.791); # make omega variable (i.e., third term)
time_grid = collect(linspace(0.,20.,tot_nepochs)); # days

# put stellar parameters of desired star
stellar_parameters = Array{starparameters}(tot_nepochs);
params = repmat([1.479,     # milliarcseconds (at pole)
    8450.,                  # Kelvin (at pole)
    omega[1],               # unitless; fractional rotational velocity
    [3,ld_grid[1]],         # limb darkening,first coefficient is for LD law type, then LD coefficients
    0.190,                  # exponent for von Zeipel law
    0.,                     # 2nd constant for rotational velocity
    57.2,                   # degrees
    -61.8,                  # degrees
    rotation_period[1]],    # days
    1,tot_nepochs);
for i = 1:tot_nepochs
    params[3,i] = omega[i];
    params[4,i][2] = ld_grid[i];
    params[9,i] = rotation_period[i];
end
rotation_angle=zeros(tot_nepochs);
for i=1:tot_nepochs
    if i>1
        rotation_angle[i] = rotation_angle[i-1] + 360./params[9,i]*(time_grid[i]-time_grid[i-1]);
    end
    stellar_parameters[i]=starparameters(params[1,i],params[2,i],params[3,i],params[4,i],params[5,i],params[6,i],params[7,i],params[8,i],rotation_angle[i],params[9,i]);
end
for i = 1:tot_nepochs
    if ((i%60 == 0) & (i > 61))
        println("i = $i")
        println("Amount of rotation in 60 frames: \t",(stellar_parameters[i].selfrotangle-stellar_parameters[i-60].selfrotangle)/360.)
        println(360./params[9,i]*(time_grid[i]-time_grid[i-60]))
        println(time_grid[i]-time_grid[i-60])
        println(" ")
    end
end

using WriteVTK
ntheta = 40; nphi = 80;
npix = ntheta*nphi;
nxyz = length(latlong_rapidrot_star(ntheta,nphi, stellar_parameters[1]).vertices_xyz[:,1,1:4]);

celltype = VTKCellTypes.VTK_QUAD;
cells = MeshCell[];
for i=1:npix
    # Define connectivity of cell.
    indx = collect(((i-1)*4+1):(i*4));

    # Define cell.
    cell = MeshCell(celltype, indx);
    push!(cells, cell);
end

pvd = paraview_collection("multi_epoch_star_spinup");
tic();
for i = 1:tot_nepochs
    star_epoch_geom = update_star(latlong_rapidrot_star(ntheta,nphi, stellar_parameters[i]), stellar_parameters[i]);
    x_val = reshape(star_epoch_geom.vertices_xyz[:,1,1:4]',1,nxyz);
    y_val = reshape(star_epoch_geom.vertices_xyz[:,2,1:4]',1,nxyz);
    z_val = reshape(star_epoch_geom.vertices_xyz[:,3,1:4]',1,nxyz);
    star_pos = Array{Float64}(3,nxyz);
    star_pos[1,:] = x_val; star_pos[2,:] = y_val; star_pos[3,:] = z_val;
    vtkfile = vtk_grid("multi_epoch_star_spinup"*string(i),star_pos,cells);


    temperature_map = calc_tempmap_vZ(stellar_parameters[i],star_epoch_geom);
    vtk_cell_data(vtkfile,temperature_map,"Temperature");
    collection_add_timestep(pvd, vtkfile, time_grid[i]);
    println("frames left: ",tot_nepochs-i); # countdown timer till finish
end
toc();

vtk_save(pvd)
