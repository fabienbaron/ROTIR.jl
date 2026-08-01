using ROTIR
# LOAD DATA
oifitsfiles = ["./data/MEDIAN5.MIRCX_L2.2025Oct30.HD_432.MIRCX_IDL.bet_Cas.AVG10m.oifits"]
data_all = readoifits_multiepochs(oifitsfiles; T=Float32);
data = data_all[1, :]; # select first wavelength bin, all epochs
nepochs = length(data)
tepochs = Float32.([d.mean_mjd for d in data])
tepochs = tepochs .- tepochs[1]; # First epoch set as t=0


# To use a Healpix scheme
n=3; 
tessels = tessellation_healpix(n, T=Float32)

## Rapid rotator
star_params = (
              surface_type   = 2      ,   # Round:0, Ellipsoid: 1, Rapid Rotator:2, Roche: 3
              rpole          = 0.849f0,   # milliarcseconds (radius at pole)
              tpole          =   7208f0, #  # Kelvin (at pole)
              ldtype         =      3,   # LD type  1: Linear 2: quadratic 3: power (Hestroffer)
              ld1            =   0.21f0, # limb darkening,first coefficient is for LD law type, then LD coefficients
              ld2            =    0f0,   # second ld coeff, used if needed
              inclination    = 19.9f0,  # degrees; inclination
              position_angle =      -7.09f0,  # degrees; position_angle
              rotation_period= 1/1.12f0,  # rotation period in days
              beta           =    0.146f0,  # exponent for von Zeipel law
              frac_escapevel =      0.92f0, # unitless; fractional rotational velocity
              B_rot =               0f0  # 2nd constant for rotational velocity              
              )

stars = create_star_multiepochs(tessels, star_params, tepochs);
tmap_start = parametric_temperature_map(star_params,stars[1]);
setup_oi!(data, stars)

# RECONSTRUCTION AND MODEL FITTING
tmap_start = convert(Vector{Float32}, tmap_start)
θ_start = Float32.([0.849, 0.92, 19.9, -7.09])   # [rpole, omega, inc, PA] for rapid rotator
#                 [rpole, omega, inc,   PA]
lower = Float32.([0.5,   0.0,   0.0, -180.0])
upper = Float32.([2.0,   0.99, 90.0,  180.0])
r_weight = Float32(1e-4)
κ_input = Float32(50.0)

tmap_final, θ_final = joint_reconstruct_oi(
    tmap_start, θ_start, data, tessels, star_params, tepochs;
    maxiter_xmap=200, maxiter_θ=100, nouter=3,
    reg_weight=r_weight,  κ=κ_input,
    θ_lower=lower,
    θ_upper=upper,
)

crit = image_reconstruct_oi_crit(tmap_final, data, stars, regularizers = [], verbose = true)
chi2 = image_reconstruct_oi_chi2(tmap_final, data, stars, verbose = true)

# tpole renormalization
tmap_final = rescale_temperature_tpole(tmap_final, stars, star_params)

plot2d(tmap_final, stars[1], contours=[6000, 6500, 7000, 7500])
# How to display non-measured pixels in mollweide
plot_mollweide(tmap_final, stars[1], visible_pixels=sometimes_visible(stars), bad_color="lightgray")