using ROTIR

# LOAD DATA
oifitsfiles = ["./data/2011Sep02.lam_And_prepped.oifits", "./data/2011Sep06.lam_And_prepped.oifits",
"./data/2011Sep10.lam_And_prepped.oifits","./data/2011Sep14.lam_And_prepped.oifits",
"./data/2011Sep19.lam_And_prepped.oifits","./data/2011Sep24.lam_And_prepped.oifits"];
data_all = readoifits_multiepochs(oifitsfiles, warn=false, verbose=false, T=Float32);
data = data_all[1, :]; # select first wavelength bin, all epochs
nepochs = length(data)
tepochs = [d.mean_mjd for d in data]
tepochs = tepochs .- tepochs[1]; # First epoch set as t=0

# To use a Healpix scheme
n=3;
tessels = tessellation_healpix(n)

## Rapid rotator
star_params = (
              surface_type   = 2,        # Round:0, Ellipsoid: 1, Rapid Rotator:2, Roche: 3
              rpole          = 1.37131,   # milliarcseconds (radius at pole)
              tpole          = 4800.0,    # Kelvin (at pole)
              ldtype         = 3,         # LD type  1: Linear 2: quadratic 3: power (Hestroffer)
              ld1            = 0.22886,   # limb darkening, first coefficient is for LD law type, then LD coefficients
              ld2            = 0.0,       # second ld coeff, used if needed
              inclination    = 78.0962,   # degrees; inclination
              position_angle = 24.0,      # degrees; position_angle
              rotation_period= 54.8,      # rotation period in days
              beta           = 0.08,      # exponent for von Zeipel law
              frac_escapevel = 0.0,       # unitless; fractional rotational velocity
              B_rot          = 0.0        # 2nd constant for rotational velocity
              )

stars = create_star_multiepochs(tessels, star_params, tepochs);
tmap_start = parametric_temperature_map(star_params,stars[1]);
setup_oi!(data, stars)

# SETUP REGULARIZATION
regularizers = [["tv2", 1e-5, tv_neighbours_healpix(n), 1:length(tmap_start)]];

# RECONSTRUCTION
tmap =  image_reconstruct_oi(tmap_start, data, stars, maxiter = 1000, regularizers = regularizers, verbose = true);
crit = image_reconstruct_oi_crit(tmap, data, stars, regularizers = [], verbose = true)
chi2 = image_reconstruct_oi_chi2(tmap, data, stars, verbose = true)
plot2d_allepochs(tmap, stars)
plot_mollweide(tmap, stars[1])
