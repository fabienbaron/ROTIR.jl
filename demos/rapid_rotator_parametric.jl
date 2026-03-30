using ROTIR, BenchmarkTools
# LOAD DATA
oifitsfiles = ["./data/2011Sep02.lam_And_prepped.oifits", "./data/2011Sep06.lam_And_prepped.oifits",
"./data/2011Sep10.lam_And_prepped.oifits","./data/2011Sep14.lam_And_prepped.oifits",
"./data/2011Sep19.lam_And_prepped.oifits","./data/2011Sep24.lam_And_prepped.oifits"];
data_all = readoifits_multiepochs(oifitsfiles; T=Float32);
data = data_all[1, :]; # select first wavelength bin, all epochs
nepochs = length(data)
tepochs = Float32.([d.mean_mjd for d in data])
tepochs = tepochs .- tepochs[1]; # First epoch set as t=0

# To use a Healpix scheme
n=3; 
tessels = tessellation_healpix(n)

## Rapid rotator
star_params = (
              surface_type   = 2      ,   # Round:0, Ellipsoid: 1, Rapid Rotator:2, Roche: 3
              rpole          = 1.37131f0,   # milliarcseconds (radius at pole)
              tpole          =   4800.0f0, #  # Kelvin (at pole)
              ldtype         =      3,   # LD type  1: Linear 2: quadratic 3: power (Hestroffer)
              ld1            =  0.22886f0, # limb darkening,first coefficient is for LD law type, then LD coefficients
              ld2            =    0.0f0,   # second ld coeff, used if needed
              inclination    = 78.0962f0,  # degrees; inclination
              position_angle =      24f0,  # degrees; position_angle
              rotation_period=    54.8f0,  # rotation period in days
              beta           =    0.08f0,  # exponent for von Zeipel law
              frac_escapevel =      0.9f0, # unitless; fractional rotational velocity
              B_rot =               0f0  # 2nd constant for rotational velocity              
              )

chi2_parametric_surface=p->spheroid_parametric_f(p, tessels, data, tepochs) 
@btime chi2_parametric_surface(star_params)
