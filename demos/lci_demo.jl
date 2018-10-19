#
# Simple light curve inversion from phase data
#
using ROTIR
# Import light curve
#lcifile = "5110407_intphaselc_pcs_0001.txt";
lcifile = "./data/KIC_flux.txt"
lcidata = read_lci_absolute(lcifile);

# Define stellar parameters for each epoch
stellar_parameters = Array{starparameters}(undef, lcidata.nepochs);
starparams = [1.,   # radius (fixed for LCI)
          3000., # temperature Kelvin (at pole)
          0.,  # frac_escapevel::Float64 # unitless; fractional rotational velocity
          [3,0.22886], # Hestroffer limb darkening coefficient
          0.08,               # exponent for von Zeipel law
          0,           # 2nd constant for rotational velocity
          60.0,            # degrees; inclination
          0,            # degrees; position_angle
          0,    # degrees: selfrotation angle, defined by phase
          0           # days; rotation_period (unused here if we want to use the phase)
          ];
for i=1:lcidata.nepochs
    stellar_parameters[i]=starparameters(starparams[1],starparams[2],starparams[3],starparams[4],starparams[5],starparams[6],starparams[7],starparams[8],360.0.*lcidata.phase[i],starparams[10]);
end

# Create 3D geometry from parameters
n = 3; #Healpix tesselation level
star_epoch_geom = create_geometry( healpix_round_star(n,stellar_parameters[1].radius), stellar_parameters);
# Initial temperature map
temperature_map_start = 3000.0*ones(star_epoch_geom[1].npix); #start from lowest possible temperature
# Setup Fourier transform (polygons -> complex visibilities)
polyflux = setup_lci(star_epoch_geom);
somvis = sometimes_visible(star_epoch_geom); # lists all pixels at least visible once
tvinfo = tv_neighbours_healpix(n)
regularizers = [["tv", 0.01, tvinfo]];
#regularization λ=0.0002, B=6000;  # bias factor from Harmonλ=0.0002

temperature_map = lci_reconstruct(temperature_map_start, lcidata, polyflux, somvis);

lciplot_vs_model(lcidata, modelflux_lci(temperature_map, polyflux));

mollplot_temperature_healpix(temperature_map);

temperature_map[never_visible(star_epoch_geom)] .= mean(temperature_map[sometimes_visible(star_epoch_geom)]);
mollplot_temp_longlat(temperature_map,ntheta,nphi);
#plot2d_temperature(temperature_map, star_epoch_geom[1]);
