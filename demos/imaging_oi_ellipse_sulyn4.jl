using ROTIR
# LOAD DATA
oifitsfiles=["./data/SU_Lyn.oifits"]
data_all = readoifits_multiepochs(oifitsfiles; filter_bad_data=true, T=Float32);
data = data_all[1, :]; # select first wavelength bin, all epochs
nepochs = length(data)
tepochs = Float32.([d.mean_mjd for d in data])
tepochs = tepochs .- tepochs[1]; # First epoch set as t=0

# To use a Healpix scheme
n=3; 
tessels = tessellation_healpix(n)

## Rapid rotator
star_params = (
surface_type   = 1f0      ,   # Round:0, Ellipsoid: 1, Rapid Rotator:2, Roche: 3
radius_x              = 1.75f0,   # milliarcseconds (radius at pole)
radius_y              = 1.75f0,
radius_z              =  1.6f0, 
tpole          =   3500.0f0, #  # Kelvin (at pole)
ldtype         =      3f0,   # LD type  1: Linear 2: quadratic 3: power (Hestroffer)
ld1            =  0.26f0, # limb darkening,first coefficient is for LD law type, then LD coefficients
ld2            =    0.0f0,   # second ld coeff, used if needed
inclination    = 76f0,  # degrees; inclination
position_angle =      49f0,  # degrees; position_angle
rotation_period=    254.8f0,  # rotation period in days
beta  =0.08f0
)

stars = create_star_multiepochs(tessels, star_params, tepochs);
tmap_start = 1000*ones(Float32, stars[1].npix)
setup_oi!(data, stars)

# SETUP REGULARIZATION
regularizers_1 = [["tv", 0.1, tv_neighbours_healpix(n),1:length(tmap_start)]];
tmap_1 =  image_reconstruct_oi(tmap_start, data, stars, maxiter = 500, lower=1000, regularizers = regularizers_1, verbose = true);

regularizers_2 = [["tv", 0.1, tv_neighbours_healpix_visible(n, stars),1:length(tmap_start)]];
tmap_2 =  image_reconstruct_oi(tmap_start, data, stars, maxiter = 500, lower=1000, regularizers = regularizers_2, verbose = true);


# Should be a slight difference in limb-darkening

plot2d(tmap_1, stars[1])


plot2d(tmap_2, stars[1])


# Example of upsampling
tmap, stars = upsample_map_stars(tmap_2, stars, star_params, tepochs)
setup_oi!(data, stars)
n = n+1
regularizers = [["tv", 0.1, tv_neighbours_healpix_visible(n, stars),1:length(tmap)]];
tmap =  image_reconstruct_oi(tmap, data, stars, maxiter = 500, lower=1000, regularizers = regularizers, verbose = true);
