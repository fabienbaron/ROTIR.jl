using OITOOLS
using OIFITS
using PyPlot
using Statistics
using Interpolations
include("../rotir-lexi/src/ROTIR.jl"); using Main.ROTIR;
# include("../rotir/src/oistars.jl")
# include("../rotir/src/keplerorbit.jl")
# include("../rotir/src/geometry_rochelobe.jl")


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

# observedmjds_2007 = [54221.0, 54228, 54232]
# observedmjds_2012 = [56087.0, 56088, 56089, 56090]
# observedmjds_2015 = [57128.0, 57129, 57130, 57131]

observed_mjds = copy(sorted_jump_mjds)
#[54221.0, 54228.0, 54232.0, 56087.0, 56088.0, 56089.0, 56090.0, 57128.0, 57129.0, 57130.0, 57131.0]

nepochs = length(observed_mjds)
nstars = 2; stellar_parameters = Array{starparameters}(undef,nstars);

binary_parameters = Array{binaryparameters}(undef,12);


# parameters for Spica
# Stellar Parameters
starparams = Array{Any}(undef,nstars,9);

# star 1 (primary) model
starparams[1,:] = [37.16,   # milliarcseconds (radius at pole)
    25300.0,                 # Kelvin (at pole)
    0.448,                    # unitless; fractional rotational velocity
    [3.0,0.15],            # limb darkening,first coefficient is for LD law type, then LD coefficients
    0.25,                   # exponent for von Zeipel law
    0.0,                    # 2nd constant for rotational velocity
    0.0,                    # degrees; inclination
    0.0,                    # degrees; position_angle
    46.07];               # days; rotation_period
# star 2 (secondary)
starparams[2,:] = [18.062, 20585.0, 0.229, [3,0.15], 0.25, 0., 0., 0., 52.53];


# # Binary parameters checks model
# binaryparams = [63.1,      # degrees; orbital inclination
#     309.938,                 # degrees; longitude of ascending node
#     255.0,                   # degrees; Argument of periapsis (periastron)
#     4.0145,                # days; binary period
#     77.0,                  # parsecs; distance
#     2.44810371190724e-9,                  # milliarcseconds; separation (aka semi-major axis)****
#     0.133,                    # unitless; eccentricity
#     0.6307,                  # unitless; mass ratio
#     1.8,                    # unitless; async ratio
#     [0.2,1.0],#0.6005975,         # unitless; fillout factor [DO NOT KNOW YET]
#     0.0,                    # days/day; linear change of the binary period
#     2454000.0];                   # JD; time of periastron

#orbit model:
binaryparams = [114.751,      # degrees; orbital inclination
    122.294,                 # degrees; longitude of ascending node
    224.303,                   # degrees; Argument of periapsis (periastron)
    4.0145,                # days; binary period
    77.0,                  # parsecs; distance
    1.58891,                  # milliarcseconds; separation (aka semi-major axis)****
    0.0980014,                    # unitless; eccentricity
    0.6307,                  # unitless; mass ratio
    1.8,                    # unitless; async ratio
    [0.2,1.0],#0.6005975,         # unitless; fillout factor [DO NOT KNOW YET]
    0.0,                    # days/day; linear change of the binary period
    2454000.0];   

# # star 1 (primary)
# starparams[1,:] = [37.16,   # milliarcseconds (radius at pole)
#     25900.0,                 # Kelvin (at pole)
#     0.448,                    # unitless; fractional rotational velocity
#     [3.0,0.15],            # limb darkening,first coefficient is for LD law type, then LD coefficients
#     0.25,                   # exponent for von Zeipel law
#     0.0,                    # 2nd constant for rotational velocity
#     0.0,                    # degrees; inclination
#     0.0,                    # degrees; position_angle
#     46.07];               # days; rotation_period
# # star 2 (secondary)
# starparams[2,:] = [18.062, 20850.0, 0.229, [3,0.15], 0.25, 0., 0., 0., 52.53];


# # Binary parameters
# binaryparams = [116.0,      # degrees; orbital inclination
#     309.938,                 # degrees; longitude of ascending node
#     255.0,                   # degrees; Argument of periapsis (periastron)
#     4.0145,                # days; binary period
#     77.0,                  # parsecs; distance
#     1.54,                  # milliarcseconds; separation (aka semi-major axis)
#     0.123,                    # unitless; eccentricity
#     0.6188,                  # unitless; mass ratio
#     1.8,                    # unitless; async ratio
#     [0.2,1.0],#0.6005975,         # unitless; fillout factor [DO NOT KNOW YET]
#     0.0,                    # days/day; linear change of the binary period
#     2440678.008];                   # JD; time of periastron


tepochs = collect(range(0.0,stop=4.0145,length=1000)); 




star1parameters = starparameters(starparams[1,1], starparams[1,2], starparams[1,3], starparams[1,4], starparams[1,5], starparams[1,6], starparams[1,7], starparams[1,8],360.0/starparams[1,9]*(1000.0), starparams[1,9])
star2parameters = starparameters(starparams[2,1], starparams[2,2], starparams[2,3], starparams[2,4], starparams[2,5], starparams[2,6], starparams[2,7], starparams[2,8], 360.0/starparams[2,9]*(1000.0), starparams[2,9])

binary_parameters = binaryparameters(star1parameters, star2parameters ,binaryparams[1],binaryparams[2],binaryparams[3],binaryparams[4],
binaryparams[5],binaryparams[6],binaryparams[7],binaryparams[8],binaryparams[9],binaryparams[10],binaryparams[11],binaryparams[12]);


# Call the binary_orbit_abs function and store the result in a variable
# orb_coords= binary_orbit_abs_alt_vec(binary_parameters, tepochs)
orb_coords= binary_orbit_abs_alt_vec(binary_parameters, tepochs)

# Extract the coordinates from the tuple
x1, y1, z1, x2, y2, z2  = orb_coords
#create 2d plot of orbits
fig = figure("Orbit", figsize=(10,10),facecolor = "White")
ax = fig.add_axes([0.1,0.1,0.85,0.85]);
ax.set_aspect("equal")
xlabel("Δ RA (mas)", fontweight="bold", fontsize=15);
ylabel("Δ Dec (mas)", fontweight="bold", fontsize=15);
ax.set_xlim(2.0, -2.0)
ax.set_ylim(-2.0, 2.0)
ax.plot(x1, z1, label ="Secondary Star")
ax.plot(x2, z2, label ="Primary Star")

#setting up xyz positions for the desired epochs
tepochs2 = collect(range(0.0,stop=4.0145,length=length(observed_mjds)));


epoch = binary_orbit_abs_alt_vec(binary_parameters, tepochs2)

scatter(epoch[1], epoch[3])
scatter(epoch[4], epoch[6])

#calculate Center of Mass of orbit and add to orbit plot



#estimate initial positions by taking the average of all points of the observed epochs
x1, y1, z1, x2, y2, z2 = epoch
avgx1 = mean(x1)
avgy1 = mean(y1)
avgx2 = mean(x2)
avgy2 = mean(y2)
#declare intial positions
r1 = [avgx1, avgy1]
r2 = [avgx2, avgy2]


#computing masses 
m1, m2 = compute_masses(binary_parameters)

#CM calculation
# CMx, CMy = (m1 * r1 + m2 *r2) / (m1 + m2)

# Define the positions of the objects at a specific epoch
epoch_positions = (r1, r2)  # r1 and r2 are vectors representing the positions

# Calculate the center of mass
CMx, CMy = (m1 * epoch_positions[1] + m2 * epoch_positions[2]) / (m1 + m2)

#add point to plot
scatter(CMx, CMy, marker = "x", c= "black", s= 50, label= "Center of Mass")
legend()

rel_coords = compute_xy_simple(binary_parameters, tepochs)
x, y, theta, rho = rel_coords
#plot orbit
fig = figure("Orbit", figsize=(10,10),facecolor = "White")
ax = fig.add_axes([0.1,0.1,0.85,0.85]);
ax.set_aspect("equal")
xlabel("Δ RA (mas)", fontweight="bold", fontsize=15);
ylabel("Δ Dec (mas)", fontweight="bold", fontsize=15);
ax.set_xlim(2.0, -2.0)
ax.set_ylim(-2.0, 2.0)
ax.plot(x, y, label = "Secondary Star")
scatter(0,0, marker = "o", c= "Blue", s=7000, label = "Primary Star")

 epochs = compute_xy_simple(binary_parameters, tepochs2)
 x, y, theta, rho = epochs

 ax.scatter(x, y, s = 300)


#take difference of orbit positions to plot orbit relative to primary star
#full orbits
relx = y1 - y2
rely = x1 - x2
#plot orbit
fig = figure("Orbit", figsize=(10,10),facecolor = "White")
ax = fig.add_axes([0.1,0.1,0.85,0.85]);
ax.set_aspect("equal")
xlabel("Δ RA (mas)", fontweight="bold", fontsize=15);
ylabel("Δ Dec (mas)", fontweight="bold", fontsize=15);
ax.set_xlim(2.0, -2.0)
ax.set_ylim(-2.0, 2.0)
ax.plot(relx, rely, label = "Secondary Star")
#epochs
ex1, ey1, ez1, ex2, ey2, ez2 = epoch
relex = ey1 - ey2
reley = ex1 - ex2
#plot epochs
scatter(relex, reley, s=100)

scatter(0,0, marker = "o", c= "black", s=7000,)
legend()


