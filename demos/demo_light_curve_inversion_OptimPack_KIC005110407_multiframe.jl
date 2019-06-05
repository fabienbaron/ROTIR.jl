#
# Light curve inversion from phase data with absolute flux - Movie version
#
using ROTIR
using Statistics
using PyPlot
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

# Split data into imaging epochs
nframes, Y, E, W, H = split_lcidata_by_period(lcidata,polyflux, 0.5*period); # we use half period, but any reasonable number should work

# Setup Harmon & Total Variation regularization
#C, ∇s, ∇w = setup_regularization_matrices(epochs_geometry);

# Reconstruction
Tphot_wanted = 5200;
ΔT_wanted = 1000;
T_low = Tphot_wanted-2*ΔT_wanted; # Threshold lowq flux
T_hi = Tphot_wanted+ΔT_wanted; # Threshold hi flux

B=5000.; λ=0.001; μ=0.00001;

xstart = [];

npix = epochs_geometry[1].npix
x_start= ones(npix*nframes);
B=5000.; λ=0.001; μ=0.00001;
regularizers = [repeat([[["bias", λ, B,visible_pixels], ["tv",μ, tv_neighbours_healpix(n), 1:npix]]], nframes);[[["temporal_tvsq",1e-3]]]];

# Reconstruction: x is the reconstructed temperature movie
xopt = lci_reconstruct_mutitemporal(x_start, Y, E, H, visible_pixels, verb = true, maxiter = 200, regularizers = regularizers);
x = rescale_temperature(xopt, Tphot_wanted, visible_pixels);

# Export movie as mollweide plot then mp4
vmin = minimum(x[visible_pixels,:]);
vmax = maximum(x[visible_pixels,:]);
for i=1:nframes
    x[hidden_pixels,i] .= median(x[visible_pixels,i]);
    mollplot_temperature_healpix(x[:,i], visible_pixels = visible_pixels, vmin = vmin, vmax = vmax);
    #readline(stdin)
    savefig("lci$i");
end

run(`ffmpeg -f image2 -pattern_type glob -framerate 5 -i 'lci%d.png' lci_all.mp4`)
run(`ls -1v *.png | xargs -I {} echo "file '{}'" > list.txt')
run(`ffmpeg -r 1/5 -f concat -i list.txt  lcimovie.mp4`)
