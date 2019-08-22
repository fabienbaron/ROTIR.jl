include("geometry_healpix.jl")
include("geometry_longlat.jl")

mutable struct base_geometry  # this is only affected by stellar parameters (e.g. Roche parameters)
  npix::Int64
  # Healpix vertices
  vertices_xyz::Array{Float64,3}
  vertices_spherical::Array{Float64,3}
end

mutable struct epoch_geometry # typically one per epoch, rotation and projection of the base geometry
  npix::Int64
  # Healpix vertices
  vertices_xyz::Array{Float64,3}
  vertices_spherical::Array{Float64,3}
  normals::Array{Float64,2}
  # Healpix projection onto the 2D imaging plane
 # quads_visible::Array{Bool,1} not necessary anymore
  index_quads_visible::Array{Int64,1}
  nquads_visible::Int64
  projx::Array{Float64,2}
  projy::Array{Float64,2}
  # Limb-darkening map
  ldmap::Array{Float64,1}
end

function oblate_const(stellar_parameters) # Approximate a rapid rotator by an oblate spheroid
    # Get oblate part using Gerard's approximation
    if (stellar_parameters.frac_escapevel >= 1.e-10)
        a = b = 3.0*stellar_parameters.radius.*cos((pi + acos(stellar_parameters.frac_escapevel*sin(pi/2.)))/3.)./
            (stellar_parameters.frac_escapevel*sin(pi/2.));
        c = stellar_parameters.radius;
    elseif (stellar_parameters.frac_escapevel <= 1.e-10)
        a = b = c = stellar_parameters.radius;
    end
    return a,b,c
end

# von Zeipel law
function calc_tempmap_vZ(stellar_parameters,star_epoch_geom; GM = 1.0)
    r_pole = stellar_parameters.radius;
    r_theta = sqrt.(star_epoch_geom.vertices_xyz[:,5,1].^2 + star_epoch_geom.vertices_xyz[:,5,2].^2 + star_epoch_geom.vertices_xyz[:,5,3].^2);
    theta = star_epoch_geom.vertices_spherical[:,5,2];
    teff_pole = stellar_parameters.temperature;

    omega_crit = sqrt.(8.0*GM/(27.0*r_pole^3));
    omega = stellar_parameters.frac_escapevel*omega_crit;
    g_r_theta = -GM./(r_theta.^2) + r_theta.*(omega*sin.(theta)).^2;
    g_theta_theta = r_theta.*(omega^2).*sin.(theta).*cos.(theta);
    g_theta = sqrt.(g_r_theta.^2 + g_theta_theta.^2);

    g_r_pole = -GM/(r_pole.^2); # second term is zero
    g_theta_pole = 0.;
    g_pole = sqrt.(g_r_pole.^2 + g_theta_pole.^2);

    # Teff(theta) = T_pole*(g(theta)/g_pole)^(beta)
    star_map = teff_pole*((g_theta/g_pole).^stellar_parameters.beta)
    return star_map
end

function calc_rotspin(R_pole,R_equ,omega_c,Mass)
    omega_k = sqrt.(8.0*((R_equ./R_pole).^3)/27.0).*omega_c;
    G = 6.67e-8; M_sun = 2.e33; R_sun = 7.e10;
    v_crit = sqrt.((2.0/3.0)*G*Mass*M_sun/(R_pole*R_sun))*(1.e-5); # km/s
    velocity = omega_c.*2.0*R_equ./(3.0*R_pole)*v_crit; # km/s
    rotation_period = 2.0*pi*R_equ.*R_sun*(1.e-5)./velocity; # s
    rotation_period /= (60.0*60.0*24.0); # day
    ang_vel = velocity./(R_equ*1.e-5*60.0*60.0*24.0); # degrees/day
    rotational_vel = ang_vel*(pi/180.0); # rotations/day
    return rotational_vel, rotation_period
end

function calc_omega(R_pole,oblate)
    R_equ = (1.0+oblate).*R_pole;
    omega_0 = 1.0 - R_pole./R_equ;
    omega = sqrt.(27.0*omega_0.*((1.0-omega_0).^2)/4.0);
    return omega, R_pole, R_equ
end

function rot_vertex(angle_r1, angle_r2, angle_r3) # new rotational matrix
dcm = Array{Float64}(undef, 3,3)
c1 = cos(angle_r1)
s1 = sin(angle_r1)
c2 = cos(angle_r2)
s2 = sin(angle_r2)
c3 = cos(angle_r3)
s3 = sin(angle_r3)
dcm[1,1] = -s1*c2*s3 + c1*c3;
dcm[1,2] = s1*c3*c2 + c1*s3;
dcm[1,3] = -s1*s2;

dcm[2,1] = -c1*c2*s3 - s1*c3;
dcm[2,2] = c1*c3*c2 - s1*s3;
dcm[2,3] = -c1*s2;

dcm[3,1] = -s2*s3;
dcm[3,2] = s2*c3;
dcm[3,3] = c2;

return dcm
end

function rot_vertex_old(angle_r1, angle_r2, angle_r3) # FB rotational matrix
dcm = Array{Float64}(undef, 3,3)
c1 = cos(angle_r1)
s1 = sin(angle_r1)
c2 = cos(angle_r2)
s2 = sin(angle_r2)
c3 = cos(angle_r3)
s3 = sin(angle_r3)
dcm[1,1] = -s1*c2*s3 + c1*c3;
dcm[1,2] = s1*c3*c2 + c1*s3;
dcm[1,3] = s1*s2;

dcm[2,1] = -c1*c2*s3 - s1*c3;
dcm[2,2] = c1*c3*c2 - s1*s3;
dcm[2,3] = c1*s2;

dcm[3,1] = s2*s3;
dcm[3,2] = -s2*c3;
dcm[3,3] = c2;

return dcm
end



function omega_rotation(A_rot, B_rot, latitude)
    #omega = A_rot + B_rot*((sin(pi/2. - latitude)).^2) + C_rot*((sin(pi/2. - latitude)).^4);
    omega = A_rot + B_rot*((sin(pi/2. - latitude)).^2);
    return omega
end

function update_star(star_base_geom, sparameters; diff_rot=false, ntheta=36, nphi=72, secondary=false,bparameters=Array{Any})
  npix = star_base_geom.npix;
  vertices_xyz = deepcopy(star_base_geom.vertices_xyz);
  vertices_spherical = deepcopy(star_base_geom.vertices_spherical);

  #
  # Rotate the star according to stellar parameters
  #
  # Differential rotation check
  if (diff_rot == false) # block rotation
      compound_rotation = rot_vertex(sparameters.selfrotangle*pi/180., sparameters.inclination*pi/180., sparameters.position_angle*pi/180.);
      vertices_xyz = reshape(reshape(vertices_xyz,(npix*5,3))*compound_rotation, (npix,5,3));
  else
    day_sep = sparameters.selfrotangle/360.0*sparameters.rotation_period;
    A_rot = 360.0/sparameters.rotation_period;
    B_rot = sparameters.B_rot;
    indx_sphere = 1;
    for i = 1:ntheta
      # CALCULATE SELFROTANGLE HERE
      latitude = vertices_spherical[indx_sphere,2,5]; # in astronomy coordinates
      #omega_rot = omega_rotation(A_rot, B_rot, C_rot, latitude);
      omega_rot = omega_rotation(A_rot, B_rot, latitude);
      new_selfrotation = omega_rot*day_sep;
      compound_rotation = rot_vertex(new_selfrotation*pi/180.0, sparameters.inclination*pi/180.0, sparameters.position_angle*pi/180.0);
      ilong_range = ((i-1)*nphi+1):(i*nphi);
      for j = ilong_range
          vertices_xyz[j,:,:] = compound_rotation*vertices_xyz[j,:,:];
      end
      indx_sphere = i*nphi + 1;
    end
  end
  #
  # Determine visible pixels based on whether their normals point toward us (z>0)
  #
  vecAC = vertices_xyz[:, 3, :]-vertices_xyz[:, 1, :];
  vecBD = vertices_xyz[:, 4, :]-vertices_xyz[:, 2, :];

  normals = Array{Float64}(undef, npix,3);  # cross.() needs arrays of arrays
  normals[:,1] = vecAC[:,2].*vecBD[:,3] - vecAC[:,3].*vecBD[:,2];
  normals[:,2] = vecAC[:,3].*vecBD[:,1] - vecAC[:,1].*vecBD[:,3];
  normals[:,3] = vecAC[:,1].*vecBD[:,2] - vecAC[:,2].*vecBD[:,1];
  normals ./= sqrt.(sum(abs2, normals, dims=2))
  #quads_visible = zeros(Bool,npix);
  #quads_visible[normals[:,3].> 0] .= true;
  #index_quads_visible = findall(quads_visible);
  index_quads_visible = findall(normals[:,3].>0);
  nquads_visible = length(index_quads_visible);

  # 2D the projection onto the (x,y) observing plane
  projx = vertices_xyz[index_quads_visible, 1:4, 1];
  projy = vertices_xyz[index_quads_visible, 1:4, 2];

  # Limb-darkening map
  ldmap = zeros(Float64,npix);
  μ = abs.(normals[index_quads_visible,3])
  if (sparameters.ld[1] == 1) # 1: quadratic
    ldmap[index_quads_visible] = 1.0 .- sparameters.ld[2]*(1.0 .-μ) - sparameters.ld[3]*(1.0.-μ.^2)
  elseif (sparameters.ld[1] == 2) # 2: logarithmic
    ldmap[index_quads_visible] = 1.0 .- sparameters.ld[2]*(1.0 .-μ) - sparameters.ld[3]*μ.*log.(μ)
  elseif (sparameters.ld[1] == 3)  # 3; Hestroffer
    ldmap[index_quads_visible] = μ.^sparameters.ld[2]
  end
  return epoch_geometry(npix, vertices_xyz, vertices_spherical, normals, index_quads_visible, nquads_visible, projx,  projy, ldmap);
end


function never_visible(star_epoch_geom)
  # return the list of pixels which are never visible at any epochs
  hidden = Array{Bool}(undef, star_epoch_geom[1].npix)
  hidden[:] .= true;
  for t=1:length(star_epoch_geom)
  hidden[star_epoch_geom[t].index_quads_visible] .= false;
  end
  return findall(hidden.==true)
end

function sometimes_visible(star_epoch_geom)
  # return the list of pixels which are at least visible during one epoch
  sometimes = Array{Bool}(undef, star_epoch_geom[1].npix)
  sometimes[:] .= false;
  for t=1:length(star_epoch_geom)
  sometimes[star_epoch_geom[t].index_quads_visible] .= true;
  end
  return findall(sometimes.==true)
end


function create_geometry(star_base_geom, stellar_parameters; kwargs...)
nepochs = size(stellar_parameters,1);
npix = star_base_geom.npix
star_epoch_geom = Array{epoch_geometry}(undef, nepochs);
println("Creating geometry for $(nepochs) epochs x $(npix) tessels");
for i=1:nepochs
  star_epoch_geom[i] = update_star(star_base_geom, stellar_parameters[i]; kwargs...);
end
return star_epoch_geom
end


#= experimental stuff
using QHull
function hull(star_epoch_geom)
println("Computing convex hull");
ch = chull(hcat(vec(star_epoch_geom.projx),vec(star_epoch_geom.projy)))
hull=[ch.vertices;ch.vertices[1]]; #convex hull, closed
scatter(ch.points[:,1], ch.points[:,2])
plot(ch.points[hull,1], ch.points[hull,2])

bounding_box = [minimum(star_epoch_geom.projx), minimum(star_epoch_geom.projy), maximum(star_epoch_geom.projx), maximum(star_epoch_geom.projy)]
end
=#



# function triangle_orientation(a,b,c)
# #CrossProductZ(a,b) = a[1] * b[2] - a[2] * b[1]
# #CrossProductZ(b,c) = b[1] * c[2] - b[2] * c[1]
# #CrossProductZ(c,a) = c[1] * a[2] - c[2] * a[1]
# #  return CrossProductZ(a, b) + CrossProductZ(b, c) + CrossProductZ(c, a)
# return a[1] * b[2] - a[2] * b[1] + b[1] * c[2] - b[2] * c[1] + c[1] * a[2] - c[2] * a[1]
# end
#
# function sort2D_quad_counter(quad)
# #given four points, return them sorted counterclockwise
# #initial order
# a = quad[1,:];
# b = quad[2,:];
# c = quad[3,:];
# d = quad[4,:];
#
# if triangle_orientation(a, b, c) > 0.0
#         # Triangle abc is already clockwise.  Where does d fit?
#         if triangle_orientation(a, c, d) > 0.0
#             return quad;
#         elseif triangle_orientation(a, b, d) > 0.0
#             return quad[[1,2,4,3],:];
#         else
#           return quad[[4,2,3,1],:];
#         end
# elseif triangle_orientation(a, c, d) > 0.0
#         # Triangle abc is counterclockwise, i.e. acb is clockwise.
#         # Also, acd is clockwise.
#         if triangle_orientation(a, b, d) > 0.0
#           return quad[[1,3,2,4],:];
#         else
#           return quad[[2,1,3,4],:];
#         end
# else
#   # Triangle abc is counterclockwise, and acd is counterclockwise.
#   # Therefore, abcd is counterclockwise.
#     return quad[[3,2,1,4],:];
# end
#
# end
