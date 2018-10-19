#
# Light curve inversion demo code for relative flux, with regularization optimization
#
using ROTIR
#include("../src/ROTIR.jl"); using Main.ROTIR;

using Statistics

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

for B in 10.0.^range(0,stop=10,length=10)
    for λ in 10.0.^range(-8,stop=4,length=15)
        regularizers = [["bias", λ, B,visible_pixels]];
        temperature_map = lci_reconstruct(temperature_map_start, lcidata, polyflux, visible_pixels, maxiter = 500, regularizers = regularizers, relative = true, verb = false);
        T_ratio=minimum(temperature_map[visible_pixels])/mean(temperature_map[visible_pixels]);
        chi2_r = chi2_lci_f(temperature_map, polyflux, lcidata)/length(lcidata.flux)
        #if abs(chi2_r - 1)<3./sqrt(length(lcidata.flux)) # this test may not be reliable when flux uncertainty is missing
            if abs(T_ratio-(Tphot_target-ΔT_target)/Tphot_target)<.1
                    temperature_map[hidden_pixels] .= mean(temperature_map[visible_pixels]);
                    mollplot_temperature_healpix(temperature_map, title = "B=$(B) - λ=$(λ)");
                    printstyled("B=$(B) \t λ=$(λ) \t chi2_r =$(chi2_r)  \t Tratio=$(T_ratio) Tratiotarget=$((Tphot_target-ΔT_target)/Tphot_target)\n", color=:red)
            else
                print("B=$(B) \t λ=$(λ) \t chi2_r =$(chi2_r)  \t Tratio=$(T_ratio) Tratiotarget=$((Tphot_target-ΔT_target)/Tphot_target)\n");
            end

    end
end





end
