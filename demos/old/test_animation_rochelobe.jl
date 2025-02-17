include("../src/ROTIR.jl"); using Main.ROTIR; using PyPlot

# star 1 (primary) model
star1parameters = starparameters(.8/2,   # milliarcseconds (radius at pole)
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

# Binary parameters Monnier 
binary_parameters = binaryparameters(star1parameters, star2parameters,
    77.0,                  # parsecs; distance
    114.751,      # degrees; orbital inclination
    122.294+0*180,                 # degrees; longitude of ascending node
    224.303-180,                   # degrees; Argument of periapsis (periastron)
    4.0145898,                # days; binary period
    1.58891,                  # milliarcseconds; separation (aka semi-major axis)
    0.0980014,                    # unitless; eccentricity
    1.31443,                   # JD; time of periastron
    0.6188,                  # unitless; mass ratio
    [1.0, 1.0],         # unitless; fillout factor 
    0.0,                    # days/day; linear change of the binary period
    0.0); # degrees/day apsidal motion period


nepochs = 20
tepochs = binary_parameters.T0 .+ collect(range(0.0,stop=4.0145,length=nepochs)); 
epochphases = mod.(tepochs .- binary_parameters.T0, binary_parameters.P)

#
ntheta = 20
nphi = 20
prim_base = latlong_ellipsoid_star(ntheta,nphi,1.0,1.0,1.0);
sec_base = latlong_ellipsoid_star(ntheta,nphi,1.0,1.0,1.0);


n=3; #768 
prim_base = healpix_round_star(n,radius=1.0)
sec_base = healpix_round_star(n,radius=1.0)

for i=1:nepochs
    primary_geom, secondary_geom = update_roche_radii(prim_base, sec_base, binary_parameters, tepochs[i]);
    primary_temperature_map = compute_teff_vonzeipel(binary_parameters, primary_geom, tepochs[i]);
    secondary_temperature_map = compute_teff_vonzeipel(binary_parameters, secondary_geom, tepochs[i], secondary = true);  
#    plot2d_wire(primary_geom)
#    plot2d_wire(secondary_geom)
#    plot3d_temperature(primary_temperature_map, primary_geom)
#    plot3d_temperature(secondary_temperature_map, secondary_geom)
    plot2d_temperature_binary(primary_temperature_map, secondary_temperature_map, primary_geom, secondary_geom, plotmesh = true)
    savefig("frame_$i.png")
    close()
end 

run(`ffmpeg -y -framerate 10 -i frame_%d.png output.mp4`)


# Plot orbit
#
x1 = zeros(nepochs)
x2 = zeros(nepochs)
y1 = zeros(nepochs)
y2 = zeros(nepochs)
z1 = zeros(nepochs)
z2 = zeros(nepochs)

for ii=1:nepochs
    x1[ii], y1[ii], z1[ii], x2[ii], y2[ii], z2[ii] = binary_orbit_rel(binary_parameters,tepochs[ii])

   # x1[ii], y1[ii], z1[ii], x2[ii], y2[ii], z2[ii] = binary_orbit_rel_alt(binary_parameters,tepochs[ii])


    #scatter([x1[ii],x2[ii]], [z1[ii],z2[ii]] )
end

plot(x1, z1, x2, z2)

xlim(2,-2)
ylim(-2,2)
#set_aspect("equal")

plot(x2-x1, z2-z1)
plot(y2-y1, z2-z1)
plot(y1, z1, y2, z2)


