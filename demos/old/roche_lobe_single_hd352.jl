using ROTIR
#include("../src/oiplot_spheroid_geomakie.jl")
using ParameterHandling;
# LOAD DATA
oifitsfiles=["./data/MIRCX_L2.2023Sep08.5_Cet.MIRCX_IDL.TEST.SPLIT.oifits"]
data_all = readoifits_multiepochs(oifitsfiles; T=Float32);
data = data_all[1, :]; # select first wavelength bin, all epochs
nepochs = length(data)
tepochs = Float32.([d.mean_mjd for d in data])
tepochs = tepochs .- tepochs[1]; # First epoch set as t=0

# To use the latitude/longitude scheme
# ntheta=50
# nphi=50
# @btime tessellation_latlong(ntheta,nphi)

# To use a Healpix scheme
n=4; 
tessels = tessellation_healpix(n)


# Binary parameters for a single visible star
roche_parameters = (  surface_type  = 3,  # Round:0, Ellipsoid: 1, Rapid Rotator:2, Roche: 3
                     rpole          = 0.6502f0,   # milliarcseconds (radius at pole)
                     tpole           = 4000.0,      # Fixed pole temperature (K)
                      ldtype         =   3,   # LD type  1: Linear 2: quadratic 3: power (Hestroffer)
                      ld1            =   0.22886f0, # limb darkening,first coefficient is for LD law type, then LD coefficients
                      ld2            =   0.0f0,   # second ld coeff, used if needed
                      inclination    =   72.0f0,  # degrees; inclination
                      position_angle =   0.0f0,  # degrees; position_angle
                      rotation_period=   95.5f0,  # rotation period in days
                      beta           =   0.0969f0,  # exponent for von Zeipel law
# Now Roche parameters
                      d              =   337.8799f0, # distance (parsecs)
                      q              =    1.40f0, # unitless, q = Mass secondary/Mass primary
                      fillout_factor_primary = 1.0, # if negative, won't be used and rpole will be used to define potential
                                                   # unitless, value of the potential at Roche lobe divided by value of potential at the surface
# And orbital parameters                      
                      i = 72.0f0, # degrees
                      Ω = 0.0f0, # degrees
                      ω = 21.2f0, # degrees
                      P = 96.4371f0, # days
                      a = 1.667f0, # semi-major axis in milliarcseconds
                      e = 0.02f0, # unitless
                      T0 = 2453602.1470, # JD; time of periastron
                      dP = 0.0f0, # days/day; linear change of the binary period
                      dω = 0.0f0 # periapsis change (degrees/day)
                );
           



# Define a corresponding free-mask: true means the parameter is free to vary.
free_mask = (
    surface_type           = false,
    rpole                  = false,
    tpole                  = false,
    ldtype                 = false,
    ld1                    = true,
    ld2                    = false,
    inclination            = true,
    position_angle         = true,
    rotation_period        = false,
    beta                   = true,
    d                      = true,
    q                      = true,
    fillout_factor_primary = false,
    i                      = true,
    Ω                      = true,
    ω                      = true,
    P                      = false,
    a                      = true,
    e                      = false,
    T0                     = false,
    dP                     = false,
    dω                     = false
)

#if binary but single, extract star params, then compute D
stars = create_star_multiepochs(tessels, roche_parameters, tepochs);

# Create a single map based on the first epoch
tmap = parametric_temperature_map(roche_parameters,stars[1]);

# Examples of plots
#plot2d_wireframe(stars[1])
#plot2d_allepochs(tmap, stars)
#plot3d_vertices(stars[1])
#plot3d(tmap, stars[1])
#plot_mollweide(tmap, stars[1])
##plot3d_makie(tmap, stars[1]) # Work in progress
##plot2d_makie(tmap, stars[1]) # Work in progress
# In the future, maybe create as many maps as epochs (useful if interactions)
#tmaps = temperature_map_vonZeipel_roche_single(roche_parameters,stars, tepochs);

# Setup the temperature-to-flux vector and the temperature-to-visility matrix
setup_oi!(data, stars)

# To get the observables are a given epoch
v2_model, t3amp_model, t3phi_model = observables(tmap, stars[1], data[1]);

# individual chi2
chi2v2, chi2t3amp, chi2t3phi = chi2s(tmap, stars[1], data[1], verbose = true);

# Total chi2 summed over all epochs
chi2_epochs = spheroid_chi2_allepochs_f(tmap, stars, data)

# Now let's time the parameters_to_chi2 chain
chi2_parametric_surface=p->spheroid_parametric_f(p, tessels, data, tepochs) 

# Test
chi2_parametric_surface(roche_parameters)

# Work in progress
# Use ParameterHandling to transform NamedTuple into parameter vector
# Helper function: extract only free parameters into a vector.
function flatten_free(params::NamedTuple, mask::NamedTuple)
    free_vals = Float32[]
    for key in keys(params)
        if mask[key]
            push!(free_vals, params[key])
        end
    end
    return free_vals
end
# Helper function: reconstruct full parameters from the free vector.
function reconstruct_full(free_vec::Vector{Float64}, params::NamedTuple, mask::NamedTuple)
    new_params = Dict{Symbol,Any}()
    i = 1
    for key in keys(params)
        if mask[key]
            new_params[key] = free_vec[i]
            i += 1
        else
            new_params[key] = params[key]  # keep fixed value
        end
    end
    return (; new_params...)  # convert to NamedTuple
end


p_free = flatten_free(roche_parameters, free_mask)

lbounds = [0.20, # rpole
0.0, #ld1
40.0, #inclination
0.0, #positionangle
0, #beta
100, #d
0.0,#q
40.0, #i
 0., #Ω,
 0, #ω
 0.25]  #a

ubounds = [2, # rpole
1.0, #ld1
90.0, #inclination
90.0, #positionangle
1, #beta
500, #d
3,#q
90.0, #i
 180, #Ω,
 180, #ω
 2]  #a

# Define your objective (chi²) function.
# In your context, spheroid_parametric_f accepts a full parameter set.
function chi2_obj(free_vec)
    # Reconstruct the full parameter set, inserting the fixed values.
    current_params = reconstruct_full(free_vec, roche_parameters, free_mask)
    # Call your model function (e.g., spheroid_parametric_f) with the full parameters.
    return spheroid_parametric_f(current_params, tessels, data, tepochs)
end
chi2_parametric_surface = (p,g)->spheroid_surface_f(p, tessels, data, tepochs);
using NLopt

# Set up NLopt with only the free parameters.
fitter = :LN_NELDERMEAD
opt = Opt(fitter, length(p_free))
min_objective!(opt, (p, grad) -> chi2_obj(p))
xtol_rel!(opt, 1e-3)
lower_bounds!(opt, lbounds)
upper_bounds!(opt, ubounds)

# Run the optimizer starting from the free parameters.
min_f, min_free, ret = optimize(opt, p_free)

# Reconstruct the full parameter set after optimization.
best_params = reconstruct_full(min_free, roche_parameters, free_mask)

println("""
    Best chi² value       : $min_f
    Best parameters       : $best_params
    Solution status       : $ret
    Number of evaluations : $(NLopt.numevals(opt))
""")

