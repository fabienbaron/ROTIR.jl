include("../src/ROTIR.jl"); using Main.ROTIR
using BenchmarkTools; 
# LOAD DATA
oifitsfiles = ["./data/2011Sep02.lam_And_prepped.oifits", "./data/2011Sep06.lam_And_prepped.oifits",
"./data/2011Sep10.lam_And_prepped.oifits","./data/2011Sep14.lam_And_prepped.oifits",
"./data/2011Sep19.lam_And_prepped.oifits","./data/2011Sep24.lam_And_prepped.oifits"];
nepochs, tepochs, data = readoifits_multiepochs(oifitsfiles);
tepochs = tepochs .- tepochs[1]; # t=0 corresponds to face on?

# # SETUP STAR MODEL PARAMETERS         
star_params = (
              surface_type  = 0      ,   # Round:0, Ellipsoid: 1, Rapid Rotator:2, Roche: 3
              radius        = 1.37131,   # milliarcseconds (radius at pole)
              temperature    =   4800.0, #  # Kelvin (at pole)
              ldtype         =      3,   # LD type
              ld1            =  0.22886, # limb darkening,first coefficient is for LD law type, then LD coefficients
              ld2            =    0.0,   # second ld coeff, used if needed
              inclination    = 78.0962,  # degrees; inclination
              position_angle =      24,  # degrees; position_angle
              rotation_period=    54.8,  # rotation period in days
              beta           =    0.08,  # exponent for von Zeipel law
              frac_escapevel =      0., # unitless; fractional rotational velocity
              B_rot =               0.  # 2nd constant for rotational velocity              
              )

n=3; tesselation = healpix_sphere(n);
chi2=p->spheroid_surface_f(p, tesselation, data, tepochs); 
chi2(star_params)
  
