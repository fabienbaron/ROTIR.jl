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

star1_parameters = (
              surface_type  = 2      ,   # Round:0, Ellipsoid: 1, Rapid Rotator:2, Roche: 3
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

# Binary parameters for a single visible star (no orbit)
binary_parameters = ( star1 = star1_parameters, 
                    d = 77.0, # distance (parsecs)
                #    i, # degrees
                #    Ω, # degrees
                #    ω, # degrees
                     P = 10.0, # days
                     a = 10.0, # in milliarcseconds (aka semi-major axis)
                #    e, # unitless
                #    T0, # JD; time of periastron
                    q=0.61, # unitless, q = Mass secondary/Mass primary
                    fillout_factor_primary = 1.0, # unitless; value of the potential at Roche lobe divided by value of potential at the surface
                #    dP, # days/day; linear change of the binary period
                #    dω # periapsis change (degrees/day)
                );

#if binary but single, extract star params, then compute D
#stars = create_star_multiepochs(tessels, star_params, tepochs);

# WORK IN PROGRESS
