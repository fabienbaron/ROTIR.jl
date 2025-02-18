include("../src/ROTIR.jl"); using Main.ROTIR
using BenchmarkTools; 
# LOAD DATA
oifitsfiles = ["./data/2011Sep02.lam_And_prepped.oifits", "./data/2011Sep06.lam_And_prepped.oifits",
"./data/2011Sep10.lam_And_prepped.oifits","./data/2011Sep14.lam_And_prepped.oifits",
"./data/2011Sep19.lam_And_prepped.oifits","./data/2011Sep24.lam_And_prepped.oifits"];
nepochs, tepochs, data = readoifits_multiepochs(oifitsfiles);
tepochs = tepochs .- tepochs[1]; # First epoch set as t=0


# To use the latitude/longitude scheme
# ntheta=50
# nphi=50
# @btime tessellation_latlong(ntheta,nphi)

# To use a Healpix scheme
n=3; 
tessels = tessellation_healpix(n)

## Rapid rotator
star_params = (
              surface_type   = 2      ,   # Round:0, Ellipsoid: 1, Rapid Rotator:2, Roche: 3
              rpole          = 1.37131,   # milliarcseconds (radius at pole)
              tpole          =   4800.0, #  # Kelvin (at pole)
              ldtype         =      3,   # LD type  1: Linear 2: quadratic 3: power (Hestroffer)
              ld1            =  0.22886, # limb darkening,first coefficient is for LD law type, then LD coefficients
              ld2            =    0.0,   # second ld coeff, used if needed
              inclination    = 78.0962,  # degrees; inclination
              position_angle =      24,  # degrees; position_angle
              rotation_period=    54.8,  # rotation period in days
              beta           =    0.08,  # exponent for von Zeipel law
              frac_escapevel =      0.9, # unitless; fractional rotational velocity
              B_rot =               0.  # 2nd constant for rotational velocity              
              )

stars = create_star_multiepochs(tessels, star_params, tepochs);
star_map = temperature_map_vonZeipel_rapid_rotator(star_params,stars[1]);
plot2d_temperature_allepochs(star_map, stars)
