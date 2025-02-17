using OITOOLS
using OIFITS
using PyPlot
using Statistics
include("../rotir/src/ROTIR.jl"); using Main.ROTIR;
include("../rotir/src/keplerorbit.jl")

oifitsfile = "2007_2012_2015_Spica.oifits";
data = (readoifits(oifitsfile))[1,1]

data.uv_mjd #displays the information stored in the data


plot(sort(data.uv_mjd))
#plots epochs in a line

sorted_mjd = sort(data.uv_mjd)

#get indices that sort original MJD values
sorted_indices = sortperm(data.uv_mjd)

#calculating derivatives for all indices
differences = diff(sorted_indices)

derivatives = []
for i in 1:length(differences) -1
    derivative = differences[i + 1] - differences[i]
    push!(derivatives, derivative)
end

#find all non-zero values 
#(i.e., where the obs day changes & obs_eppchs line jumps)
jump_positions = findall(derivatives .!=0)

jump_mjds = sorted_mjd[jump_positions] #calculating the corresponding MJD values at those jump points

# sort the values in jump_mjds and keep only the items with a 
# difference of exactly one integer (no decimals), you can use the round
#  function to round the values to the nearest integer and then use the 
# unique function to remove any duplicate values.

sorted_jump_mjds = unique(round.(Int, jump_mjds))

#re-read oifit file
# data = (readoifits("2007_2012_2015_Spica.oifits",mjd_range = ))[1]


# observedmjds_2007 = [54221.0, 54228, 54232]
# observedmjds_2012 = [56087.0, 56088, 56089, 56090]
# observedmjds_2015 = [57128.0, 57129, 57130, 57131]

observed_mjds = [54221.0, 54228.0, 54232.0, 56087.0, 56088.0, 56089.0, 56090.0, 57128.0, 57129.0, 57130.0, 57131.0]

nepochs = length(observed_mjds)
nstars = 2; stellar_parameters = Array{starparameters}(undef,nstars);
binary_parameters = Array{binaryparameters}(undef,12);


# parameters for Spica
# Stellar Parameters
starparams = Array{Any}(undef,nstars,9);

# star 1 (primary)
starparams[1,:] = [37.16,   # milliarcseconds (radius at pole)
    25900.0,                 # Kelvin (at pole)
    0.448,                    # unitless; fractional rotational velocity
    [3.0,0.15],            # limb darkening,first coefficient is for LD law type, then LD coefficients
    0.25,                   # exponent for von Zeipel law
    0.0,                    # 2nd constant for rotational velocity
    0.0,                    # degrees; inclination
    0.0,                    # degrees; position_angle
    46.07];               # days; rotation_period
# star 2 (secondary)
starparams[2,:] = [18.062, 20850.0, 0.229, [3,0.15], 0.25, 0., 0., 0., 52.53];


# Binary parameters
binaryparams = [116.0,      # degrees; orbital inclination
    309.938,                 # degrees; longitude of ascending node
    255.0,                   # degrees; Argument of periapsis (periastron)
    4.0145,                # days; binary period
    77.0,                  # parsecs; distance
    1.54,                  # milliarcseconds; separation (aka semi-major axis)
    0.123,                    # unitless; eccentricity
    0.6188,                  # unitless; mass ratio
    1.8,                    # unitless; async ratio
    [0.2,1.0],#0.6005975,         # unitless; fillout factor [DO NOT KNOW YET]
    0.0,                    # days/day; linear change of the binary period
    2440678.008];                   # JD; time of periastron


tepochs = collect(range(0.0,stop=4.0145,length=1000)); 




star1parameters = starparameters(starparams[1,1], starparams[1,2], starparams[1,3], starparams[1,4], starparams[1,5], starparams[1,6], starparams[1,7], starparams[1,8],360.0/starparams[1,9]*(1000.0), starparams[1,9])
star2parameters = starparameters(starparams[2,1], starparams[2,2], starparams[2,3], starparams[2,4], starparams[2,5], starparams[2,6], starparams[2,7], starparams[2,8], 360.0/starparams[2,9]*(1000.0), starparams[2,9])

binary_parameters = binaryparameters(star1parameters, star2parameters ,binaryparams[1],binaryparams[2],binaryparams[3],binaryparams[4],
binaryparams[5],binaryparams[6],binaryparams[7],binaryparams[8],binaryparams[9],binaryparams[10],binaryparams[11],binaryparams[12]);

#plotting full orbit
xyz = binary_orbit_py_alt(binary_parameters,tepochs)

x = xyz[1]
y = xyz[2]
scatter(x,y, s = 4)


#setting up xyz positions for the desired epochs
tepochs2 = collect(range(0.0,stop=4.0145,length=length(observed_mjds)));


xyz2 = binary_orbit_py_alt(binary_parameters, tepochs2)

#will be points for markers @ the epochs
x2=xyz2[1]
y2=xyz2[2]

PyPlot.hold(true)

scatter(x2, y2, c = "red", marker = "x", s = 60)
PyPlot.hold(true)
scatter(0,0, c = "red", marker = "o", s = 60)
