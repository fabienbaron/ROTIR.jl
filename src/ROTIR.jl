module ROTIR
using OITOOLS
using Statistics
using LinearAlgebra

include("oistars.jl");
include("geometry.jl");
include("oichi2_spheroid.jl");
include("oiplot_spheroid.jl");

# Re-export OITOOLS functions so users only need `using ROTIR`
export OIdata, readoifits, readoifits_multiepochs, readfits, writefits

# Tessellation
export tessellation
export tessellation_healpix, tessellation_latlong
export nside2npix, npix2n
export healpix_round_star, healpix_ellipsoid_star, tv_neighbours_healpix, tv_neighbours_healpix_visible
export latlong_round_star, latlong_ellipsoid_star, tv_neighbors_longlat, latlong_rochelobe
export upsample_map_stars, downsample_map_stars

# Stellar/binary parameters
export starparameters, binaryparameters
export update_star, update_binary

# Geometry: stars and binaries
export create_star, create_star_multiepochs, create_binary
export rotate_single_star, compute_separation
export oblate_const

# Geometry: Roche lobe
export update_roche_radii, get_surface_potential, update_roche_geom
export compute_potential_primary, compute_potential_secondary, solve_radius
export compute_dpotential_primary_L1, solve_R_L1, newton_raphson
export radius_eggleton, radius_leahy

# Geometry: rapid rotators
export temperature_map_vonZeipel_rapid_rotator, calc_grelmap_vZ
export compute_teff_vonzeipel
export calc_rotspin, calc_omega

# Geometry: temperature maps
export temperature_map_vonZeipel_roche_single, compute_gravity_primary, compute_gravity_secondary
export temperature_map_vonZeipel_ellipsoid
export parametric_temperature_map, spheroid_parametric_f

# Orbits
export compute_eccentric_anomaly, compute_true_anomaly
export compute_E_NR, compute_coeff, compute_xyz_rel
export binary_orbit_rel, binary_orbit_abs, binary_RV, binary_proj_plane

# OI chi2 and reconstruction
export setup_oi!, setup_polygon_ft, setup_polyflux_single, setup_polyft_single, setup_polyft_single_alt
export observables, chi2s, mod360
export spheroid_chi2_f, spheroid_chi2_fg, spheroid_chi2_fg_alt
export spheroid_chi2_allepochs_fg, spheroid_chi2_allepochs_f
export spheroid_total_variation, spheroid_crit_multiepochs_fg
export spheroid_l2_fg, spheroid_harmon_bias_fg, spheroid_regularization
export image_reconstruct_oi, image_reconstruct_oi_crit, image_reconstruct_oi_chi2, image_reconstruct_oi_chi2_fg
export multires_reconstruct_oi
export cvis_to_v2, poly_to_cvis, poly_to_flux, cvis_to_t3

# Plotting
export plot2d, plot2d_wireframe, plot2d_allepochs
export plot3d, plot3d_vertices
export plot_mollweide
export sometimes_visible, never_visible, invisible_neighbours, with_invisible_neighbours, without_invisible_neighbours

# Utilities
export make_circ_spot, make_spot_move
export rl1, max_rpole
export lci_linear_inversion_frame, lci_reconstruct_mutitemporal, rescale_temperature

end
