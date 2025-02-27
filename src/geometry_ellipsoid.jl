@views function temperature_map_vonZeipel_ellipsoid(stellar_parameters,star; offsets = [0.0,0.0,0.0], T=Float32)
    toff= T.(offsets)'
    GM = T(GM)
    r_theta = sqrt.(dropdims(sum(abs2, (star.vertices_xyz[:,5,:] .- toff), dims=2), dims=2));
    r_pole = 1*radius_x # to check, is it really x?
    theta = star.vertices_spherical[:,5,2];
    teff_pole = stellar_parameters.tpole;
    g_theta = 1 ./(r_theta.^2) 
    g_pole = 1/rpole^2
    star_map = teff_pole*(g_theta/g_pole).^stellar_parameters.beta
    return star_map
  end
  