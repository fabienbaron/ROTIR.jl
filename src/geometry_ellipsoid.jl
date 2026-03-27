@views function temperature_map_vonZeipel_ellipsoid(stellar_parameters,star; offsets = [0.0,0.0,0.0], T=Float32)
    p = convert_params(T, stellar_parameters)
    toff= T.(offsets)'
    r_theta = sqrt.(dropdims(sum(abs2, (star.vertices_xyz[:,5,:] .- toff), dims=2), dims=2));
    rpole = p.radius_x # to check, is it really x?
    theta = star.vertices_spherical[:,5,2];
    g_theta = 1 ./(r_theta.^2)
    g_pole = 1/rpole^2
    star_map = p.tpole*(g_theta/g_pole).^p.beta
    return star_map
  end
  