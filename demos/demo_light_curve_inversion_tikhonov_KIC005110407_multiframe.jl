#
# Simple light curve inversion with fixed regularization
# Method: matrix inversion
#
using Statistics
using LinearAlgebra
using SparseArrays
include("../src/ROTIR.jl"); using Main.ROTIR;
lcifile = "./data/kplr005110407_LC_CBV_Q02.txt"
lcidata = read_lci_absolute(lcifile);

# Period of rotation
period = 3.4693; # based on previous analysis

# For another star, to determine an approximate period with Lomb Scargle, uncomment (requires LombScargle package)
# using LombScargle
#pgram = lombscargle(lcidata.mjd, lcidata.flux, lcidata.fluxerr, maximum_frequency=1.0)
#period = findmaxperiod(pgram); plot(freqpower(pgram)...);
#println("Period determined by Lomb Scargle: $(period) days.")

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
          period           # days; rotation_period (unused here if we want to use the phase)
          ];
#lcidata.phase = mod.(lcidata.mjd,starparams[10]);
epochs_parameters = setup_stellar_parameters_lci(base_parameters, lcidata);
# Create 3D geometry from parameters
n = 4; #Healpix tesselation level
epochs_geometry = create_geometry( healpix_round_star(n), epochs_parameters);
polyflux, visible_pixels, hidden_pixels = setup_lci(epochs_geometry);

# Split data into epochs
nframes, Y, E, W, H = split_lcidata_by_period(lcidata,polyflux, period);
# Setup Harmon & Total Variation regularization
C, ∇s, ∇w = setup_regularization_matrices(epochs_geometry);

Tphot_wanted = 5200;
ΔT_wanted = 1000;
T_low = Tphot_wanted-2*ΔT_wanted # Threshold lowq flux
T_hi = Tphot_wanted+ΔT_wanted # Threshold hi flux

#
# Linear inversion method (faster)
#
B=5000.; λ=0.001; μ=0.00001;
xstart = [];
using PyPlot
npix = epochs_geometry[1].npix
x= Array{Float64}(undef, npix, nframes);
for i=1:nframes
 if i==1
     xstart= ones(npix)
 else
     xstart = x[:,i-1]
 end
 println("Reconstructing frame $i")
 x[:,i]=lci_linear_inversion_frame(xstart, H[i], W[i], Y[i], C, ∇s, ∇w, B, λ , μ);
 x[hidden_pixels,i] .= median(x[visible_pixels,i])
 println("Chi2: ", (H[i]*x[:,i]-Y[i])'*W[i]*(H[i]*x[:,i]-Y[i])/length(Y[i]))
 mollplot_temperature_healpix(x[:,i]/mean(x[visible_pixels,i])*Tphot_wanted, visible_pixels = visible_pixels)
 savefig("lci$i")
end
