include("./lib/geometry.jl");
include("./lib/readoifits.jl");
include("./lib/oichi2_spheroid.jl");
include("./lib/oiplot.jl");
include("./lib/oistars.jl");
using OptimPack

oifitsfiles = ["./DATA/2011Sep02.lam_And_prepped.oifits", "./DATA/2011Sep06.lam_And_prepped.oifits",
"./DATA/2011Sep10.lam_And_prepped.oifits","./DATA/2011Sep14.lam_And_prepped.oifits",
"./DATA/2011Sep19.lam_And_prepped.oifits","./DATA/2011Sep24.lam_And_prepped.oifits"];

nepochs, tepochs, data = readoifits_multiepochs(oifitsfiles)#, redundance_chk=true);

# get the total number of v2, t3amp, t3phi points from all epochs
nv2 = zeros(nepochs); nt3amp = zeros(nepochs); nt3phi = zeros(nepochs);
for i=1:nepochs
    nv2[i] = data[i].nv2;
    nt3amp[i] = data[i].nt3amp;
    nt3phi[i] = data[i].nt3phi;
end
tot_nv2 = sum(nv2); tot_nt3amp = sum(nt3amp); tot_nt3phi = sum(nt3phi);
tot_pts = tot_nv2 + tot_nt3amp + nt3phi;
npts = sum(tot_pts);

stellar_parameters = Array{starparameters}(nepochs);

# Initial parameters
# starparams = [1.375,#1.37131,  # milliarcseconds (at pole)
#     4800.,              # Kelvin (at pole)
#     0.,                 # unitless; fractional rotational velocity
#     [3,0.22886],        # limb darkening,first coefficient is for LD law type, then LD coefficients
#     0.08,               # exponent for von Zeipel law
#     -1.43075,           # 2nd constant for rotational velocity
#     71.0962,            # degrees; inclination
#     18.6033,            # degrees; position_angle
#     51.8513];           # days; rotation_period

starparams = [1.3795,4800.,0.,[3,0.22886],0.08,-1.43075,78.0,21.5,54.0]; # Rob Park's parameters

for i=1:nepochs
    stellar_parameters[i]=starparameters(starparams[1],starparams[2],starparams[3],starparams[4],starparams[5],
        starparams[6],starparams[7],starparams[8],360./starparams[9]*(tepochs[i]-tepochs[1]),starparams[9]);
end

# Create 3D geometry from parameters

# Healpix
# n=3;star_epoch_geom = create_geometry( healpix_round_star(n,stellar_parameters[1].radius), stellar_parameters);

# Longitude/latitude
const ntheta = 20; const nphi = 40; # for lambda And (Nyquist sampled for its angular diameter)
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


# If doing total variation, then the following will do an L-curve to find the best weight value
# Can skip if you know the best weight value
tv_weights = logspace(-15,15,61);
#tv_weights = logspace(-5,5,11);
lcurve_chi2 = zeros(length(tv_weights));
lcurve_reg = zeros(length(tv_weights));
@time begin
for i=1:length(tv_weights)
    temperature_map = reconstruct(temperature_map_start, polyflux, polyft,tv_weights[i],south_neighbors,west_neighbors,south_neighbors_reverse,west_neighbors_reverse);
    gradient = Array{Float64}(size(temperature_map));
    for t=1:3 # make sure we converged
        temperature_map = reconstruct(temperature_map, polyflux, polyft,tv_weights[i],south_neighbors,west_neighbors,south_neighbors_reverse,west_neighbors_reverse);
    end
    lcurve_chi2[i] = chi2_allepochs_tv_fg(temperature_map, gradient, tv_weights[i], south_neighbors,west_neighbors, south_neighbors_reverse, west_neighbors_reverse, epochs_weights, polyflux, polyft, data);
    lcurve_reg[i] = regularization(temperature_map,gradient,tv_weights[i],south_neighbors, west_neighbors, south_neighbors_reverse, west_neighbors_reverse);
    print_with_color(200,"i = $i \n",bold=true)
end
end

# Plot of regularization vs. chi2
loglog(lcurve_reg, lcurve_chi2); plot(lcurve_reg, lcurve_chi2,"o"); xlabel(L"\mu"*"*Regularization"); ylabel(L"\chi^{2}"); # could use xticks([array of numbers],[array of string])
# Plot of reordered regularization vs. chi2
loglog(lcurve_reg[sortperm(lcurve_reg)],lcurve_chi2[sortperm(lcurve_reg)]); plot(lcurve_reg[sortperm(lcurve_reg)],lcurve_chi2[sortperm(lcurve_reg)],"o"); xlabel(L"\mu"*"*Regularization"); ylabel(L"\chi^{2}");
# Plot of mu vs. chi2 (reordered)
loglog(tv_weights[sortperm(lcurve_reg)],lcurve_chi2[sortperm(lcurve_reg)]); plot(tv_weights[sortperm(lcurve_reg)],lcurve_chi2[sortperm(lcurve_reg)],"o"); xlabel(L"\mu"); ylabel(L"\chi^{2}");
# Plot of mu vs. chi2 (without reordering)
loglog(tv_weights,lcurve_chi2); plot(tv_weights,lcurve_chi2,"o"); xlabel(L"\mu"); ylabel(L"\chi^{2}");

# Define regularization weight
mu = tv_weights[27]; # longlat
#mu = tv_weights[26]; # healpix

function crit_modeling(params::Vector, dummy::Vector)
    println(params);
    for i=1:nepochs
        stellar_parameters[i]=starparameters(starparams[1],starparams[2],starparams[3],starparams[4],starparams[5],
            starparams[6],starparams[7],starparams[8],360./starparams[9]*(tepochs[i]-tepochs[1]),starparams[9]);
    end    
    a,b,c = oblate_const(stellar_parameters[1]);
    star_epoch_geom = create_geometry(latlong_ellipsoid_star(ntheta,nphi, a, b, c),stellar_parameters,diff_rot=false, ntheta=ntheta, nphi=nphi);
    temperature_map_start = calc_tempmap_vZ(stellar_parameters[1],star_epoch_geom[1]);
    polyflux_new, polyft_new = setup_polygon_ft(data, star_epoch_geom);
    #crit_imaging = (x,g)->chi2_allepochs_fg(x, g, epochs_weights, polyflux_new, polyft_new, data);
    #crit_imaging = (x,g)->chi2_allepochs_l2_fg(x, g, mu, epochs_weights, polyflux_new, polyft_new, data);
    crit_imaging = (x,g)->chi2_allepochs_tv_fg(x, g, mu, south_neighbors,west_neighbors, south_neighbors_reverse, west_neighbors_reverse, epochs_weights, polyflux_new, polyft_new, data);
    temperature_map = OptimPack.vmlmb(crit_imaging, temperature_map_start, verb=true, lower=0, maxiter=150, blmvm=false);

    f = chi2_allepochs_f(temperature_map, epochs_weights, polyflux_new, polyft_new, data);
    return f
end

# Grid search with just angular size and rotation period
#ang_radius_grid = collect(linspace(1.3,1.45,11));
#ang_radius_grid = collect(linspace(1.33,1.42,11));
#ang_radius_grid = collect(linspace(1.37131-1.37131/50.,1.37131+1.37131/50.,25));
ang_radius_grid = collect(linspace(1.325,1.385,121));
#period_grid = collect(linspace(45.,70.,11)); 
#period_grid = collect(linspace(32.5,57.5,11)); # then try 26
period_grid = collect(linspace(20.,62.,25));
#LD_grid = collect(linspace(0.1,0.4,21));
LD_grid = collect(linspace(0.175,0.225,51));


chi2_map_angrad_period = Array{Float64}(length(ang_radius_grid),length(period_grid));
chi2_map_angrad_LD = Array{Float64}(length(ang_radius_grid),length(LD_grid)); 

starparams_arr = Array{Any}(length(starparams),length(ang_radius_grid),length(LD_grid));

for i=1:length(LD_grid)
    starparams_arr[:,:,i] = repmat(starparams,1,length(ang_radius_grid));
end
for i=1:length(ang_radius_grid)
    starparams_arr[1,i,:] = ang_radius_grid[i];
    for j=1:length(LD_grid)
        starparams_arr[4,:,j] = repmat([[3,LD_grid[j]]],length(ang_radius_grid),1);
    end
end

# Angular radius vs period
# try map to make this faster
@time begin
for i = 1:length(ang_radius_grid)
    starparams[1] = ang_radius_grid[i];
    for j = 1:length(period_grid)
        print_with_color(200,"i = $i, j = $j \n",bold=true);
        starparams[9] = period_grid[j];
        chi2_map_angrad_period[i,j] = crit_modeling(starparams,Array{Float64}(1));
        print_with_color(200,"chi2 = ",chi2_map_angrad_period[i,j]," \n",bold=true);
    end
    println("chi2 for all period = ", chi2_map_angrad_period)
end
end

# Angular radius vs LD
@time begin
for i = 71:73
    #starparams[1] = ang_radius_grid[i];
    for j = 1:length(LD_grid)
        print_with_color(200,"i = $i, j = $j \n",bold=true);
        #starparams[4] = [3,LD_grid[j]];
        #print_with_color(200,"angular radius = ",starparams[1],", LD = ",starparams[4],"\n",bold=true);
        print_with_color(200,"angular radius = ",starparams_arr[1,i,j],", LD = ",starparams_arr[4,i,j],"\n",bold=true);
        chi2_map_angrad_LD[i,j] = crit_modeling(starparams_arr[:,i,j],Array{Float64}(1));
        print_with_color(200,"chi2 = ",chi2_map_angrad_LD[i,j]," \n",bold=true);
    end
    #println("chi2 for all period = ", chi2_map_angrad_LD)
end
end

"""
v0 -> initial test 11x11 grid with only 20 iterations
v1 -> initial test with 11x11 grid having 150 iterations
v2 -> test with 11x51 grid having 150 iterations and reverting back to original rotational matrix
v3 -> test using Rob Park's original parameters
v4 -> try a different range for the period (see if it yields better results) [45.,70.] -> [32.5,57.5]
    change the angular diameter grid to something smaller [1.3,1.45] -> [1.3,1.45]
v5 -> try original range for period but smaller angular diameter grid and now using l2 norm with Rob's 
    parameters (there's overfitting with just positivity)
v6 -> constrain angular diameter grid even more but do wider range for the period
"""

"""
v0 -> initial test for angular radius and LD with large grid size (using Rob's parameters)
v1 -> refine parameters for angular diameter and LD
"""
nversion = "1";

# To save chi2 values
using HDF5
h5open("chi2_map_angrad_LDv"*nversion*"pt23.h5", "w") do file
  write(file, "chi2_map_angrad_LD", chi2_map_angrad_LD)
  write(file, "ang_radius_grid", ang_radius_grid)
  write(file, "LD_grid", LD_grid)
end

# To resume
chi2_map = h5read("chi2_map_angrad_LDv"*nversion*"pt23.h5", "chi2_map_angrad_LD")
