#Generate fake light-curve
include("./lib/lci.jl");
include("./lib/lciplot.jl");
include("./lib/geometry.jl");
include("./lib/oiplot.jl");
include("./lib/oistars.jl");

nepochs = 400;
phase = linspace(0,1,nepochs);

stellar_parameters = Array{starparameters}(nepochs);
params = [1.,   # radius (fixed for LCI)
          3000., # temperature
          0.,   # fractional critical velocity (assume 0 for non-rapid rotators)
          60.0, # inclination angle
          0.0, # position angle
          0,    # rotation angle (for LCI we use phase instead)
          0.0,  # rotation period (for Kepler LCI, does not matter)
          0.2,  # Hestroffer limb darkening coefficient
          0.08, # von Zeipel exponent
          0,    # differential rotation coefficient 1
          0     # differential rotation coefficient 2
          ];
for i=1:nepochs
    stellar_parameters[i]=starparameters(params[1],params[2],params[3],params[4],params[5],360.0.*phase[i],params[7],params[8],params[9],params[10],params[11]);
end

#n = 4; #Healpix tesselation level
#star_epoch_geom = create_geometry( healpix_round_star(n,stellar_parameters[1].radius), stellar_parameters);
ntheta = 40; nphi = 80;
a,b,c = oblate_const(stellar_parameters[1]);
star_epoch_geom = create_geometry(latlong_ellipsoid_star(ntheta,nphi, a, b, c),stellar_parameters,ntheta=ntheta,nphi=nphi);
polyflux = setup_lci(star_epoch_geom);

include("lci_fakemap.jl");

flux = polyflux*x_true;
fluxerr = 5e-4*maximum(flux)*ones(nepochs);
flux += fluxerr.*randn(nepochs);
lcidata  = LCI(phase, flux, fluxerr, length(phase));
write_lci("KIC_flux.txt", lcidata)
