#
# Simple light curve inversion demo code for relative flux, with fixed regularization
#
using ROTIR
using Statistics
#include("../src/ROTIR.jl"); using Main.ROTIR;

# Import light curve
lcifile = "./data/5110407_intphaselc_pcs_0001.txt";
ΔT_target = 1000;
Tphot_target = 5200.;
T_low = 3000.0 ; #lower accceptable temperature, for optimization purpose

lcidata = read_lci_relative(lcifile);
lcidata.flux ./=mean(lcidata.flux) #renorm
lcidata.fluxerr ./= 20. # this is an example of data with unknown error bars...


stellar_parameters = Array{starparameters}(undef, lcidata.nepochs);
starparams = [1.,   # radius (fixed for LCI)
          T_low, # temperature Kelvin (at pole)
          0.,  # frac_escapevel::Float64 # unitless; fractional rotational velocity
          [3,0.2], # Hestroffer limb darkening coefficient
          0.08,               # exponent for von Zeipel law
          0,           # 2nd constant for rotational velocity
          30.0,            # degrees; inclination
          0,            # degrees; position_angle
          0,    # degrees: selfrotation angle, defined by phase
          0           # days; rotation_period (unused here if we want to use the phase)
          ];
for i=1:lcidata.nepochs
    stellar_parameters[i]=starparameters(starparams[1],starparams[2],starparams[3],starparams[4],starparams[5],starparams[6],starparams[7],starparams[8],360.0.*lcidata.phase[i],starparams[10]);
end

# Create 3D geometry from parameters
n = 4; #Healpix tesselation level
star_epoch_geom = create_geometry( healpix_round_star(n), stellar_parameters);
polyflux = setup_lci(star_epoch_geom);
visible_pixels = sometimes_visible(star_epoch_geom); # lists all pixels at least visible once
hidden_pixels =  never_visible(star_epoch_geom);

# Initial temperature map
temperature_map_start = T_low*ones(star_epoch_geom[1].npix); #start from lowest possible temperature
#temperature_map_start = polyflux'*lcidata.flux # alternative way to start

λ=0.1; B=100.0;  # bias regularization (see e.g. Harmon et al., 2000)
regularizers = [["bias", λ, B,visible_pixels]];

temperature_map = lci_reconstruct(temperature_map_start, lcidata, polyflux, visible_pixels, maxiter = 500, regularizers = regularizers, relative = true, verb = true);

lciplot_vs_model_phase(lcidata, modelflux_lci_rel(temperature_map, polyflux));
temperature_map[hidden_pixels] .= mean(temperature_map[visible_pixels]);
mollplot_temperature_healpix(temperature_map);

end
