using OITOOLS
include("../src/ROTIR.jl"); using Main.ROTIR;
# LOAD DATA
oifitsfile="./data/rhoCas_merged_mircx_except_2023Oct28.oifits"
nepochs = 1
tepochs = [0.0]
data = [readoifits(oifitsfile, filter_bad_data=true, use_vis=false )[1,1]];

# SETUP STAR MODEL PARAMETERS
stellar_parameters = Array{starparameters}(undef, nepochs);
starparams = [2.1/2,  # milliarcseconds (radius at pole)
     7200.,              # Kelvin (at pole)
     0.,                 # unitless; fractional rotational velocity
     [1, 0.0, 0.0],        # limb darkening,first coefficient is for LD law type, then LD coefficients
     0.0,               # exponent for von Zeipel law
     0.,           # 2nd constant for rotational velocity
     90.,            # degrees; inclination
     0.,            # degrees; position_angle
     0.];           # days; rotation_period

for i=1:nepochs
    stellar_parameters[i]=starparameters(starparams[1],starparams[2],starparams[3],starparams[4],starparams[5],
        starparams[6],starparams[7],starparams[8],0.0,starparams[9]);
end

# SETUP 3D GEOMETRY (HEALPIX)
n=4; star_epoch_geom = create_geometry( healpix_round_star(n,radius=stellar_parameters[1].radius), stellar_parameters);
polyflux, polyft = setup_polygon_ft(data, star_epoch_geom);
temperature_start = 7200*ones(star_epoch_geom[1].npix); # temperature map -- initial guess
f0_start = 0#1e-3 # background flux fraction -- initial guess
tvinfo = tv_neighbours_healpix(n);
indx_visible = sometimes_visible(star_epoch_geom)
regularizers = [["tv",5e-4, tvinfo,1:length(temperature_start)]]#,["tv2",5e-6, tvinfo,1:length(temperature_start)] ];
#regularizers = [["tv2",5e-6, tvinfo,1:length(temperature_start)] ];
xx_start = vcat(f0_start,temperature_start) # All parameters to be optimized
maxiter = 50
#xx_sol = optimize_bg_flux(xx_start, polyflux, polyft, data);
xx_sol =  spheroid_oi_reconstruct(xx_start, data, polyflux, polyft, regularizers = regularizers, verb = true, maxiter=maxiter);

# Get observables and chi2s from temperature map + background flux
v2_model, t3amp_model, t3phi_model = observables(xx_sol, polyflux[1], polyft[1], data[1]);
chi2_v2, chi2_t3amp, chi2_t3phi = chi2s(xx_sol, polyflux[1], polyft[1], data[1], verbose = true)

f0 = xx_sol[1]; # background flux
temperature = xx_sol[2:end]; # temperature map
plot2d_temperature(temperature, star_epoch_geom[1], figtitle="H-band", xlim=[1.5, -1.5], ylim=[-1.5, 1.5]); # Note: Temp = intensity since no LDD
