include("../src/ROTIR.jl"); using Main.ROTIR;
lcifile = "./data/kplr005110407_LC_CBV_Q02.txt"
lcidata = read_lci_absolute(lcifile);
# Define stellar parameters common to all epochs
base_parameters = [1.,   # radius (fixed for LCI)
          5200., # temperature Kelvin (at pole)
          0.,  # frac_escapevel::Float64 # unitless; fractional rotational velocity
          [2,0.7248,0.1941], # Logarithmic law for limb darkening, coefficients from RMR
          0.,               # exponent for von Zeipel law
          0,           # 2nd constant for rotational velocity
          60.0,            # degrees; inclination
          0,            # degrees; position_angle
          NaN,    # degrees: selfrotation angle, defined by phase
          3.4693           # days; rotation_period (unused here if we want to use the phase)
          ];
#lcidata.phase = mod.(lcidata.mjd,starparams[10]);
epochs_parameters = setup_stellar_parameters_lci(base_parameters, lcidata);

# Create 3D geometry from parameters
n = 4; #Healpix tesselation level
epochs_geometry = create_geometry( healpix_round_star(n), epochs_parameters);
polyflux, visible_pixels, hidden_pixels = setup_lci(epochs_geometry);

# Setup regularization
λ=0.0004; B=6000.0;  # bias regularization (see e.g. Harmon et al., 2000)
x_start= ones(epochs_geometry[1].npix);
regularizers = [["bias", λ, B,visible_pixels]];

x = lci_reconstruct(x_start, lcidata, polyflux, visible_pixels, regularizers = regularizers, lowtemp= 0, maxiter = 500, relative = false);
lciplot_vs_model_phase(lcidata, modelflux_lci(temperature_map, polyflux));

temperature_map[hidden_pixels] .= mean(temperature_map[visible_pixels]);
mollplot_temperature_healpix(temperature_map);
