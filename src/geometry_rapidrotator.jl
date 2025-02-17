function f_rapid_rot(x)
   return 3*cos.((pi .+ acos.(x))/3)./x;
end


@views function update_radii_rapidrot(tessels::tessellation, star_parameters)
  # Return radius
  rpole = star_parameters.rpole;
  ω = star_parameters.frac_escapevel;
  r = rpole * f_rapid_rot(ω*sin.(tessels.unit_spherical[:,:,2])); 
  # Fix for ω sin \theta  -->0
  r[r .== Inf] .= rpole;
  # Rewrite pole radius values
  if tessels.tessellation_type==1 # Longitude/Latitude
    # Overwrite pole radius values with exact values
    # top of star
    r[1:tessels.nphi,1,1] .= rpole;
    r[1:tessels.nphi,4,1] .= rpole;
    # bottom of star
    r[(end-tessels.nphi+1):end,2,1] .= rpole;
    r[(end-tessels.nphi+1):end,3,2] .= rpole;
  end
  return r
end

function oblate_const(stellar_parameters) # Approximate a rapid rotator by an oblate spheroid
    # Get oblate part using Gerard's approximation
    if (stellar_parameters.frac_escapevel >= 1.e-10)
      a = b = 3.0*stellar_parameters.radius.*cos((pi + acos(stellar_parameters.frac_escapevel*sin(pi/2.0)))/3.0)./
        (stellar_parameters.frac_escapevel*sin(pi/2.));
      c = stellar_parameters.radius;
   elseif (stellar_parameters.frac_escapevel <= 1.e-10)
     a = b = c = stellar_parameters.radius;
   end
   return a,b,c
end
 

function calc_rotspin(rpole,R_equ,omega_c,Mass)
    omega_k = sqrt.(8.0*((R_equ./rpole).^3)/27.0).*omega_c;
    G = 6.67e-8; M_sun = 2.e33; R_sun = 7.e10;
    v_crit = sqrt.((2.0/3.0)*G*Mass*M_sun/(rpole*R_sun))*(1.e-5); # km/s
    velocity = omega_c.*2.0*R_equ./(3.0*rpole)*v_crit; # km/s
    rotation_period = 2.0*pi*R_equ.*R_sun*(1.e-5)./velocity; # s
    rotation_period /= (60.0*60.0*24.0); # day
    ang_vel = velocity./(R_equ*1.e-5*60.0*60.0*24.0); # degrees/day
    rotational_vel = ang_vel*(pi/180.0); # rotations/day
    return rotational_vel, rotation_period
end

# this is the fractional angular velocity (Keplerian angular velocity)
function calc_omega(rpole,oblate)
    R_equ = (1.0+oblate).*rpole;
    omega_0 = 1.0 - rpole./R_equ;
    omega = sqrt.(27.0*omega_0.*((1.0-omega_0).^2)/4.0);
    return omega, rpole, R_equ
end

# function calc_grelmap_vZ(stellar_parameters,star_epoch_geom; offsets = [0.0,0.0,0.0], GM = 1.0)
#     delx = offsets[1]; dely = offsets[2]; delz = offsets[3];
#     rpole = stellar_parameters.radius;
#     r_theta = sqrt.((star_epoch_geom.vertices_xyz[:,5,1] .- delx).^2 + (star_epoch_geom.vertices_xyz[:,5,2] .- dely).^2 + (star_epoch_geom.vertices_xyz[:,5,3] .- delz).^2);
#     theta = star_epoch_geom.vertices_spherical[:,5,2];
#     teff_pole = stellar_parameters.temperature;

#     omega_crit = sqrt.(8.0*GM/(27.0*rpole^3));
#     omega = stellar_parameters.frac_escapevel*omega_crit;
#     g_r_theta = -GM./(r_theta.^2) + r_theta.*(omega*sin.(theta)).^2;
#     g_theta_theta = r_theta.*(omega^2).*sin.(theta).*cos.(theta);
#     g_theta = sqrt.(g_r_theta.^2 + g_theta_theta.^2);

#     g_rpole = -GM/(rpole.^2); # second term is zero
#     g_theta_pole = 0.0;
#     g_pole = sqrt.(g_rpole.^2 + g_theta_pole.^2);

#     return g_theta / g_pole
# end

# von Zeipel law
@views function calc_tempmap_vZ(stellar_parameters,star_epoch_geom; offsets = [0.0,0.0,0.0], GM = 1.0)
  delx = offsets[1]; 
  dely = offsets[2]; 
  delz = offsets[3];
  rpole = stellar_parameters.radius;
  r_theta = sqrt.((star_epoch_geom.vertices_xyz[:,5,1] .- delx).^2 + (star_epoch_geom.vertices_xyz[:,5,2] .- dely).^2 + (star_epoch_geom.vertices_xyz[:,5,3] .- delz).^2);
  theta = star_epoch_geom.vertices_spherical[:,5,2];
  teff_pole = stellar_parameters.tpole;

  omega_crit = sqrt.(8.0*GM/(27.0*rpole^3));
  omega = stellar_parameters.frac_escapevel*omega_crit;
  g_r_theta = -GM./(r_theta.^2) + r_theta.*(omega*sin.(theta)).^2;
  g_theta_theta = r_theta.*(omega^2).*sin.(theta).*cos.(theta);
  g_theta = sqrt.(g_r_theta.^2 + g_theta_theta.^2);

  g_rpole = -GM/(rpole.^2); # second term is zero
  g_theta_pole = 0.0;
  g_pole = sqrt.(g_rpole.^2 + g_theta_pole.^2);
  star_map = teff_pole*((g_theta/g_pole).^stellar_parameters.beta)
  return star_map
end

# Modified von Zeipel law aka "Espinosa Lara-Rieutord law"
# https://www.aanda.org/articles/aa/pdf/2011/09/aa17252-11.pdf
# TBD!!
#=function calc_tempmap_ELR(stellar_parameters,star_epoch_geom)
    # rpole/R(theta)
    star_radius_ratio = stellar_parameters.radius./(sqrt.(star_epoch_geom.vertices_xyz[:,1,5].^2 +
        star_epoch_geom.vertices_xyz[:,2,5].^2 + star_epoch_geom.vertices_xyz[:,3,5].^2));

    # teff(theta) = t_pole*(g(theta)/g_pole)^(beta_vZ)
    star_map = stellar_parameters.temperature.*((star_radius_ratio.^4) -
        2.*8./27.*star_radius_ratio.*stellar_parameters.frac_escapevel.*(sin.(star_epoch_geom.vertices_spherical[:,2,5])).^2 +
        ((8.*star_radius_ratio/27.).^2).*((stellar_parameters.frac_escapevel.*sin.(star_epoch_geom.vertices_spherical[:,2,5])).^4) +
        ((8.*star_radius_ratio.*sin.(star_epoch_geom.vertices_spherical[:,2,5]).*cos.(star_epoch_geom.vertices_spherical[:,2,5])/27.).^2).*
        (stellar_parameters.frac_escapevel.^4)).^(stellar_parameters.beta_vZ/2.0);
    return star_map
end=#


#=function calc_tempmap_ELR(stellar_parameters,star_epoch_geom)
    # r_tilde = R/R_equator; unitless
    r_tilde = stellar_parameters.radius./maximum(stellar_parameters.radius);



    star_radius_ratio = stellar_parameters.radius./(sqrt.(star_epoch_geom.vertices_xyz[:,1,5].^2 +
        star_epoch_geom.vertices_xyz[:,2,5].^2 + star_epoch_geom.vertices_xyz[:,3,5].^2));

    # teff(theta) = t_pole*(g(theta)/g_pole)^(beta_vZ)
    star_map = stellar_parameters.temperature.*((star_radius_ratio.^4) -
        2.0*8.0/27.*star_radius_ratio.*stellar_parameters.frac_escapevel.*(sin.(star_epoch_geom.vertices_spherical[:,2,5])).^2 +
        ((8.*star_radius_ratio/27.).^2).*((stellar_parameters.frac_escapevel.*sin.(star_epoch_geom.vertices_spherical[:,2,5])).^4) +
        ((8.*star_radius_ratio.*sin.(star_epoch_geom.vertices_spherical[:,2,5]).*cos.(star_epoch_geom.vertices_spherical[:,2,5])/27.).^2).*
        (stellar_parameters.frac_escapevel.^4)).^(stellar_parameters.beta_vZ/2.0);
    return star_map
end=#
