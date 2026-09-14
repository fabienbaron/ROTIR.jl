function f_rapid_rot(x)
   return 3*cos.((pi .+ acos.(x))/3)./x;
end


@views function update_radii_rapidrot(tessels::tessellation, star_parameters)
  # Return radius
  rpole = star_parameters.rpole;
  ω = star_parameters.frac_escapevel;
  arg = ω*sin.(tessels.unit_spherical[:,:,2]);
  r = rpole * f_rapid_rot(arg);
  # Fix for ω sin(θ) → 0 (poles): f_rapid_rot is 0/0, limit is 1 → r = rpole
  # Tolerance must exceed Float32 epsilon (~1.2e-7) to catch sin(Float32(π)) ≈ -8.7e-8
  r[abs.(arg) .< 1e-5] .= rpole;
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
      a = b = 3.0*stellar_parameters.rpole.*cos((pi + acos(stellar_parameters.frac_escapevel*sin(pi/2.0)))/3.0)./
        (stellar_parameters.frac_escapevel*sin(pi/2.));
      c = stellar_parameters.rpole;
   elseif (stellar_parameters.frac_escapevel <= 1.e-10)
     a = b = c = stellar_parameters.rpole;
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

# function calc_grelmap_vZ(stellar_parameters,star; offsets = [0.0,0.0,0.0], GM = 1.0)
#     delx = offsets[1]; dely = offsets[2]; delz = offsets[3];
#     rpole = stellar_parameters.rpole;
#     r_theta = sqrt.((star.vertices_xyz[:,5,1] .- delx).^2 + (star.vertices_xyz[:,5,2] .- dely).^2 + (star.vertices_xyz[:,5,3] .- delz).^2);
#     theta = star.vertices_spherical[:,5,2];
#     teff_pole = stellar_parameters.tpole;

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
@views function temperature_map_vonZeipel_rapid_rotator(stellar_parameters, star; offsets = [0.0,0.0,0.0], GM = 1.0, T=eltype(star))
  p = convert_params(T, stellar_parameters)
  toff= T.(offsets)'
  GM = T(GM)
  r_theta = sqrt.(dropdims(sum(abs2, (star.vertices_xyz[:,5,:] .- toff), dims=2), dims=2));
  theta = star.vertices_spherical[:,5,2];
  omega_crit = T(sqrt(8*GM/(27*p.rpole^3)));
  omega = p.frac_escapevel*omega_crit;
  g_r_theta = -GM./(r_theta.^2) + r_theta.*(omega*sin.(theta)).^2;
  g_theta_theta = omega^2*r_theta.*sin.(theta).*cos.(theta);
  g_theta = sqrt.(g_r_theta.^2 + g_theta_theta.^2);
  g_pole = GM/p.rpole^2
  star_map = p.tpole*(g_theta/g_pole).^p.beta
  return star_map
end



# function omega_rotation(A_rot, B_rot, latitude)
#   #omega = A_rot - B_rot*((sin(pi/2. - latitude)).^2) - C_rot*((sin(pi/2. - latitude)).^4);
#   omega = A_rot - B_rot*((sin.(pi/2.0 - latitude)).^2); #A_rot - B_rot*((cos.(latitude)).^2)
#   return omega
# end


# =========================================================================================
# The Espinosa Lara & Rieutord law on the mesh
# =========================================================================================
# Two dead `calc_tempmap_ELR` sketches used to sit here. Both expanded `g_eff²` by hand and
# neither carried the latitudinal flux factor that IS the ELR result, so both were von Zeipel
# with algebra errors. The real law lives in `src/gravity_darkening.jl`; this is its mesh-side
# entry point, written to mirror `temperature_map_vonZeipel_rapid_rotator` line for line so
# the only difference between the two maps is the factor `F_ω(θ)/F_ω(0)`.
#
# `q = ω̃² r̃³` is taken from the MESH here — `r̃ = r(θ)/R_e` with `r(θ)` the tessel's own
# radius — rather than from the Roche shape factor the gradient path uses. On a mesh built by
# `compute_radii` the two agree identically, since that is where the shape factor comes from;
# taking it from the mesh keeps the flux factor consistent with the `g_eff` computed from the
# same radii, which is what matters if a caller ever hands in a perturbed surface.

"""
    temperature_map_ELR_rapid_rotator(star_params, star; offsets, GM=1, T) -> Vector

The Espinosa Lara & Rieutord (2011) temperature map on a rapid rotator's tessellation.

    T_eff(θ) = tpole · [ (F_ω(θ)/F_ω(0)) · (g_eff(θ)/g_pole) ]^β

Same arguments and same return as [`temperature_map_vonZeipel_rapid_rotator`](@ref), which
it reduces to as `frac_escapevel → 0`. At `beta = 0.25` this is their eq. (31); `beta` is
free so the two laws can be compared at equal parameter count. See
[`elr_flux_factor`](@ref) for the flux factor and `src/gravity_darkening.jl` for the
derivation.
"""
@views function temperature_map_ELR_rapid_rotator(stellar_parameters, star;
                                                  offsets = [0.0,0.0,0.0], GM = 1.0,
                                                  T=eltype(star))
  p = convert_params(T, stellar_parameters)
  toff = T.(offsets)'
  GM = T(GM)
  r_theta = sqrt.(dropdims(sum(abs2, (star.vertices_xyz[:,5,:] .- toff), dims=2), dims=2));
  theta = star.vertices_spherical[:,5,2];
  omega_crit = T(sqrt(8*GM/(27*p.rpole^3)));
  omega = p.frac_escapevel*omega_crit;
  g_r_theta = -GM./(r_theta.^2) + r_theta.*(omega*sin.(theta)).^2;
  g_theta_theta = omega^2*r_theta.*sin.(theta).*cos.(theta);
  g_theta = sqrt.(g_r_theta.^2 + g_theta_theta.^2);
  g_pole = GM/p.rpole^2
  # The ELR correction. `ω̃ = Ω/Ω_k` is the Roche rotation rate, NOT `omega` above (which
  # carries GM and the radius); `R_e` is the equatorial radius the same shape factor gives.
  fev = T(p.frac_escapevel)
  ωt  = elr_omega(fev)
  Re  = p.rpole * f_rapid_rot_and_deriv(fev)[1]
  # The pole's flux factor is the normalisation that makes `tpole` the polar temperature.
  # Its `q` is `(8/27) fev²`, since the shape factor is 1 on the axis.
  Fp  = elr_flux_factor(first(elr_q_and_deriv(fev, zero(T))), zero(T))
  Fth = similar(r_theta)
  @inbounds for i in eachindex(r_theta)
      Fth[i] = elr_flux_factor(ωt^2 * (r_theta[i] / Re)^3, T(theta[i]))
  end
  return p.tpole .* ((Fth ./ Fp) .* (g_theta ./ g_pole)) .^ p.beta
end

"""
    temperature_map_rapid_rotator(star_params, star; offsets, GM=1, T) -> Vector

The rapid rotator's temperature map under whichever gravity-darkening law `star_params`
names, from its `gravity_law` field: `:vonzeipel` (the default, and what a model without the
field means) or `:elr`. See [`gravity_law_spec`](@ref) for the laws and
`src/gravity_darkening.jl` for why the choice matters.
"""
function temperature_map_rapid_rotator(stellar_parameters, star; kwargs...)
    law = gravity_law_name(stellar_parameters)
    return law === :elr ?
        temperature_map_ELR_rapid_rotator(stellar_parameters, star; kwargs...) :
        temperature_map_vonZeipel_rapid_rotator(stellar_parameters, star; kwargs...)
end

