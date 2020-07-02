include("../src/ROTIR.jl"); using Main.ROTIR
# make vectors
using PyPlot
file_loc = "./rapid";
nepochs = 100;
inclination = collect(range(-90.0,90.0,length=nepochs)); # no limb darkening, no rapid rotation
posangle = collect(range(0,0,length=nepochs)); # no limb darkening, no rapid rotation
tepochs = collect(range(0.,0.,length=nepochs)); # days

# put stellar parameters of desired star
starparams = [1.,     # milliarcseconds (at pole)
    8000.,                  # Kelvin (at pole)
    0,               # unitless; fractional rotational velocity
    [3,.2],         # limb darkening,first coefficient is for LD law type, then LD coefficients
    0,                  # exponent for von Zeipel law
    0.,                     # 2nd constant for rotational velocity
    0,                   # inclination degrees
    0,                  # position angle degrees
    1.0];

stellar_parameters = Array{starparameters}(undef, nepochs);
for i=1:nepochs
    stellar_parameters[i]=starparameters(starparams[1],starparams[2],starparams[3],starparams[4],starparams[5],
        starparams[6],inclination[i],posangle[i],360.0/starparams[9]*(tepochs[i]-tepochs[1]),starparams[9]);
end

ntheta = 40; nphi = 80;
star_epoch_geom = create_geometry(latlong_ellipsoid_star(ntheta,nphi, 1.0, 1.0, 1.0),stellar_parameters,diff_rot=false, ntheta=ntheta, nphi=nphi);
temperature_map = 8000*ones(star_epoch_geom[1].npix);
plot2d_temperature(temperature_map, star_epoch_geom[1])
for i=1:nepochs
    plot2d_temperature(temperature_map, star_epoch_geom[i])
    savefig(string("wire_longlat_",string(i,pad=3),".png"))
end
run(`ffmpeg -f image2 -pattern_type glob -framerate 2 -i 'wire_longlat*.png' longlat.avi`)

star_epoch_geom = create_geometry( healpix_round_star(4,radius=1.0), stellar_parameters);
temperature_map = 8000*ones(star_epoch_geom[1].npix);
for i=1:nepochs
    plot2d_temperature(temperature_map, star_epoch_geom[i])
    savefig(string("wire_healpix_",string(i,pad=3),".png"))
end
run(`ffmpeg -f image2 -pattern_type glob -framerate 2 -i 'wire_healpix*.png' healpix.avi`)
