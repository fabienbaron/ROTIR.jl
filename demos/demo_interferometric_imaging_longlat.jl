using ROTIR

# LOAD DATA
oifitsfiles = ["./data/2011Sep02.lam_And_prepped.oifits", "./data/2011Sep06.lam_And_prepped.oifits",
"./data/2011Sep10.lam_And_prepped.oifits","./data/2011Sep14.lam_And_prepped.oifits",
"./data/2011Sep19.lam_And_prepped.oifits","./data/2011Sep24.lam_And_prepped.oifits"];
nepochs, tepochs, data = readoifits_multiepochs(oifitsfiles);

# SETUP STAR MODEL PARAMETERS
stellar_parameters = Array{starparameters}(undef, nepochs);
starparams = [1.37131,  # milliarcseconds (radius at pole)
     4800.,              # Kelvin (at pole)
     0.,                 # unitless; fractional rotational velocity
     [3,0.22886],        # limb darkening,first coefficient is for LD law type, then LD coefficients
     0.08,               # exponent for von Zeipel law
     -1.43075,           # 2nd constant for rotational velocity
     71.0962,            # degrees; inclination
     18.6033,            # degrees; position_angle
     51.8513];           # days; rotation_period

for i=1:nepochs
    stellar_parameters[i]=starparameters(starparams[1],starparams[2],starparams[3],starparams[4],starparams[5],
        starparams[6],starparams[7],starparams[8],360.0/starparams[9]*(tepochs[i]-tepochs[1]),starparams[9]);
end

# SETUP 3D GEOMETRY (Longitude/latitude)
const ntheta = 20; const nphi = 40;
a,b,c = oblate_const(stellar_parameters[1]);
star_epoch_geom = create_geometry(latlong_ellipsoid_star(ntheta,nphi, a, b, c),stellar_parameters,diff_rot=false, ntheta=ntheta, nphi=nphi);
polyflux, polyft = setup_polygon_ft(data, star_epoch_geom);

# SETUP INITIAL MAP
temperature_map_start = calc_tempmap_vZ(stellar_parameters[1],star_epoch_geom[1]);#4800*ones(star_epoch_geom[1].npix);

# SETUP REGULARIZATION
tvinfo = tv_neighbors_longlat(ntheta,nphi);
regularizers = [["tv", 0.01, tvinfo,1:length(temperature_map_start)]];

# RECONSTRUCTION
temperature_map =  spheroid_oi_reconstruct(temperature_map_start, data, polyflux, polyft, regularizers = regularizers, verb = true);

#
# PLOTS
#

# Ortho plots
plot2d_temperature_allepochs(temperature_map, star_epoch_geom, tepochs = tepochs);
plot2d_intensity_allepochs(temperature_map, star_epoch_geom, tepochs = tepochs);

# Mollweide plot
temperature_map_withhiddenblack=copy(temperature_map);
temperature_map_withhiddenblack[never_visible(star_epoch_geom)].=0;
plot3d_temperature(temperature_map_withhiddenblack,star_epoch_geom[1]);
mollplot_temperature_longlat(temperature_map_withhiddenblack, ntheta, nphi); # long/lat pix scheme
