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

# To use the latitude/longitude scheme
# ntheta=50
# nphi=50
# @btime tessellation_latlong(ntheta,nphi)

# To use a Healpix scheme
n=3; 
tessels = tessellation_healpix(n)

# Binary parameters for a single visible star
roche_parameters = (  surface_type  = 3,  # Round:0, Ellipsoid: 1, Rapid Rotator:2, Roche: 3
                      rpole          = 1.37131f0,   # milliarcseconds (radius at pole)
                      tpole          =   4800.0f0, #  # Kelvin (at pole)
                      ldtype         =      3,   # LD type  1: Linear 2: quadratic 3: power (Hestroffer)
                      ld1            =  0.22886f0, # limb darkening,first coefficient is for LD law type, then LD coefficients
                      ld2            =    0.0f0,   # second ld coeff, used if needed
                      inclination    = 78.0962f0,  # degrees; inclination
                      position_angle =      24.0f0,  # degrees; position_angle
                      rotation_period=    54.8f0,  # rotation period in days
                      beta           =    0.08f0,  # exponent for von Zeipel law
# Now Roche parameters
                      d              =    77.0f0, # distance (parsecs)
                      q              =    0.61f0, # unitless, q = Mass secondary/Mass primary
                      fillout_factor_primary = 1.0f0, # unitless; value of the potential at Roche lobe divided by value of potential at the surface
# And orbital parameters                      
                      i = 0.0f0, # degrees
                      Ω = 0.0f0, # degrees
                      ω = 0.0f0, # degrees
                      P = 10.0f0, # days
                      a = 10.0f0, # semi-major axis in milliarcseconds
                      e = 0.0f0, # unitless
                      T0 = 0.0f0, # JD; time of periastron
                      dP = 0.0f0, # days/day; linear change of the binary period
                      dω = 0.0f0 # periapsis change (degrees/day)
                );

#if binary but single, extract star params, then compute D
stars = create_star_multiepochs(tessels, roche_parameters, tepochs);

# Create a single map based on the first epoch
star_maps = temperature_map_vonZeipel_roche_single(roche_parameters,stars[1], tepochs[1]);

# In the future, maybe create as many maps as epochs (useful if interactions)
#star_maps = temperature_map_vonZeipel_roche_single(roche_parameters,stars, tepochs);

# Setup the temperature-to-flux vector and the temperature-to-visility matrix
setup_oi!(data, stars)

v2_model, t3amp_model, t3phi_model = observables(star_maps, stars[1], data[1]);
chi2v2, chi2t3amp, chi2t3phi = chi2s(star_maps, stars[1], data[1])


# plot2d_temperature_allepochs(star_maps, stars) #to update for variable maps
