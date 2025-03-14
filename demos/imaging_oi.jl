include("../src/ROTIR.jl"); using Main.ROTIR
include("../src/oitools-ext.jl")

using BenchmarkTools; 
# LOAD DATA
oifitsfiles = ["./data/2011Sep02.lam_And_prepped.oifits", "./data/2011Sep06.lam_And_prepped.oifits",
"./data/2011Sep10.lam_And_prepped.oifits","./data/2011Sep14.lam_And_prepped.oifits",
"./data/2011Sep19.lam_And_prepped.oifits","./data/2011Sep24.lam_And_prepped.oifits"];
nepochs, tepochs, data = readoifits_multiepochs(oifitsfiles);
tepochs = tepochs .- tepochs[1]; # First epoch set as t=0
data = dataF32.(data)
tepochs = Float32.(tepochs)

# To use a Healpix scheme
n=3; 
tessels = tessellation_healpix(n)

## Rapid rotator
star_params = (
              surface_type   = 2      ,   # Round:0, Ellipsoid: 1, Rapid Rotator:2, Roche: 3
              rpole          = 1.37131f0,   # milliarcseconds (radius at pole)
              tpole          =   4800f0, #  # Kelvin (at pole)
              ldtype         =      3,   # LD type  1: Linear 2: quadratic 3: power (Hestroffer)
              ld1            =  0.22886f0, # limb darkening,first coefficient is for LD law type, then LD coefficients
              ld2            =    0.0f0,   # second ld coeff, used if needed
              inclination    =   4.65f0,  # degrees; inclination
              position_angle =      24f0,  # degrees; position_angle
              rotation_period= 55.3044f0,  # rotation period in days
              beta           =    0.08f0,  # exponent for von Zeipel law
              frac_escapevel =      0.001f0, # unitless; fractional rotational velocity
              B_rot =               0.  # 2nd constant for rotational velocity              
              )

stars = create_star_multiepochs(tessels, star_params, tepochs);
tmap_start = parametric_temperature_map(star_params,stars[1]);
setup_oi!(data, stars)

# SETUP REGULARIZATION
regularizers = [["tv2", 0.01f0, tv_neighbours_healpix(n),1:length(tmap_start)]];

# RECONSTRUCTION
tmap =  spheroid_oi_reconstruct(tmap_start, data, stars, regularizers = regularizers, verbose = true);
