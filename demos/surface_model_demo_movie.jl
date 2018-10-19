include("./lib/geometry.jl");
include("./lib/oiplot.jl");
include("./lib/oistars.jl")

"""
Warning! Running this code will save 3600 images to the folder in which you are
running this code!
"""


# make vectors
file_loc = "./rapid";
nepochs = 600; nepochs_2 = 2*nepochs; tot_nepochs = 4*nepochs+nepochs_2;
ld_grid = collect(linspace(0.0,0.0,nepochs)); # no limb darkening, no rapid rotation
ld_grid = vcat(ld_grid,collect(linspace(0.,1.0,nepochs))); # raise up LD
ld_grid = vcat(ld_grid,collect(linspace(1.0,0.2,nepochs))); # put LD down to 0.2
ld_grid = vcat(ld_grid,collect(linspace(0.2,0.2,nepochs_2))); # introduce rapid rotation
ld_grid = vcat(ld_grid,collect(linspace(0.2,0.2,nepochs))); # have rapid rotation
oblate = collect(linspace(1.e-5,1.e-5,nepochs)); # no limb darkening, no rapid rotation
oblate = vcat(oblate,collect(linspace(1.e-5,1.e-5,nepochs))); # raise up LD
oblate = vcat(oblate,collect(linspace(1.e-5,1.e-5,nepochs))); # put LD down to 0.2
oblate = vcat(oblate,collect(logspace(-5,-0.6057970289817278,nepochs_2))); # introduce rapid rotation
oblate = vcat(oblate,collect(linspace(0.24785801713586308,0.24785801713586308,nepochs))); # have rapid rotation
omega, R_pole, R_equ = calc_omega(1.634,oblate); # true omega for Atlair is 0.923
rotational_vel, rotation_period = calc_rotspin(1.634,R_equ,omega,1.791); # make omega variable (i.e., third term)
time_grid = collect(linspace(0.,20.,tot_nepochs)); # days

# put stellar parameters of desired star
stellar_parameters = Array{starparameters}(tot_nepochs);
params = repmat([1.479,     # milliarcseconds (at pole)
    8450.,                  # Kelvin (at pole)
    omega[1],               # unitless; fractional rotational velocity
    [3,ld_grid[1]],         # limb darkening,first coefficient is for LD law type, then LD coefficients
    0.190,                  # exponent for von Zeipel law
    0.,                     # 2nd constant for rotational velocity
    57.2,                   # degrees
    -61.8,                  # degrees
    rotation_period[1]],    # days
    1,tot_nepochs);
for i = 1:tot_nepochs
    params[3,i] = omega[i];
    params[4,i][2] = ld_grid[i];
    params[9,i] = rotation_period[i];
end
rotation_angle=zeros(tot_nepochs);
for i=1:tot_nepochs
    if i>1
        rotation_angle[i] = rotation_angle[i-1] + 360./params[9,i]*(time_grid[i]-time_grid[i-1]);
    end
    stellar_parameters[i]=starparameters(params[1,i],params[2,i],params[3,i],params[4,i],params[5,i],params[6,i],params[7,i],params[8,i],rotation_angle[i],params[9,i]);
end
for i = 1:tot_nepochs
    if ((i%60 == 0) & (i > 61))
        println("i = $i")
        println("Amount of rotation in 60 frames: \t",(stellar_parameters[i].selfrotangle-stellar_parameters[i-60].selfrotangle)/360.)
        println(360./params[9,i]*(time_grid[i]-time_grid[i-60]))
        println(time_grid[i]-time_grid[i-60])
        println(" ")
    end
end

ntheta = 40; nphi = 80;
# Save figures! Warning, will take a few hours
tic();
for i = 1:tot_nepochs
    #a,b,c = oblate_const(stellar_parameters[i]);
    #star_epoch_geom = update_star(latlong_ellipsoid_star(ntheta,nphi, a, b, c), stellar_parameters[i]);
    star_epoch_geom = update_star(latlong_rapidrot_star(ntheta,nphi, stellar_parameters[i]), stellar_parameters[i]);
    temperature_map = calc_tempmap_vZ(stellar_parameters[i],star_epoch_geom);
    plot2d_temperature_savefig(temperature_map, star_epoch_geom,file_loc=file_loc,labels=true,iteration=i,
        omega=omega[i],rotational_vel=rotational_vel[i],LD=ld_grid[i],rotation_period=rotation_period[i]);
    println(tot_nepochs-i) # countdown timer till finish
end
toc();

run(`ffmpeg -f image2 -pattern_type glob -framerate 60 -i 'rapid*.png' rapid_all.mp4`)
