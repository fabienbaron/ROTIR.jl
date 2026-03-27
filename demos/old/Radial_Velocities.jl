include("../src/ROTIR.jl"); using Main.ROTIR; using PyPlot

# star 1 (primary) model
star1parameters = starparameters(.93/2,   # milliarcseconds (radius at pole)
    25300.0,                 # Kelvin (at pole)
    0.0,                    # unitless; fractional rotational velocity
    [3.0,0.15],            # limb darkening,first coefficient is for LD law type, then LD coefficients
    0.205,                   # exponent for von Zeipel law
    0.0,                    # 2nd constant for rotational velocity
    0.0,                    # degrees; inclination
    0.0,                    # degrees; position_angle
    0.0,
    7.226);               # days; rotation_period

# star 2 (secondary)
star2parameters = starparameters(0.57/2, 20585.0, 0.0, [3,0.15], 0.205, 0., 0., 0., 0., 7.226);

# Binary parameters
binary_parameters = binaryparameters(star1parameters, star2parameters,
    77.0,                  # parsecs; distance
    116.0,      # degrees; orbital inclination
    309.938,                 # degrees; longitude of ascending node
    233.0,                   # degrees; Argument of periapsis (periastron)
    4.0145,                # days; binary period
    1.53,                  # milliarcseconds; separation (aka semi-major axis)
    0.123,                    # unitless; eccentricity
    2454189.40,                   # JD; time of periastron
    0.6188,                  # unitless; mass ratio
    [1.0, 1.0],         # unitless; fillout factor 
    0.0,                    # days/day; linear change of the binary period
    0.0); # degrees/day apsidal motion period


nepochs = 1000
tepochs = binary_parameters.T0 .+ collect(range(0.0, stop=binary_parameters.P, length=1000))
phases = mod.(tepochs .- binary_parameters.T0, binary_parameters.P)/binary_parameters.P

rad1, rad2 = binary_RV(binary_parameters, tepochs, K1 = 123.9, K2 = 198.8, γ=0.0)
fig = figure("Radial Velocities", figsize=(10,10),facecolor = "White")
xlim(0.0, 1.0)
ylim(-250.0, 250.0)
scatter(phases, rad1, label = "Primary Star")
scatter(phases, rad2, label = "Secondary Star")

using DelimitedFiles
data1 = readdlm("./all_rv_1_ORIG.txt")
t1 = data1[:,1] 
ϕ1 = mod.(t1.- binary_parameters.T0, binary_parameters.P)/binary_parameters.P
rv1 = data1[:,2] 
data2 = readdlm("./all_rv_2_ORIG.txt")  
t2 = data2[:,1] 
ϕ2 = mod.(t2.- binary_parameters.T0, binary_parameters.P)/binary_parameters.P
rv2 = data2[:,2]

scatter(ϕ1, rv1, label ="RV1 data")
scatter(ϕ2, rv2, label ="RV2 data")
legend()

