module ROTIR
import OITOOLS: OIdata,
    readoifits, readoifits_multiepochs, readfits, writefits,
    set_data_filter, filter_data,
    plot_v2_residuals, plot_t3amp_residuals, plot_t3phi_residuals, plot_residuals
using Statistics
using LinearAlgebra
using Printf
using PrecompileTools

include("oistars.jl");
include("soft_visibility.jl");

# Convert all AbstractFloat fields of a NamedTuple to type T in one shot.
# Integer fields (surface_type, ldtype) pass through unchanged.
convert_params(::Type{T}, p::NamedTuple) where T = NamedTuple{keys(p)}(map(v -> v isa AbstractFloat ? T(v) : v, values(p)))

include("geometry.jl");
include("oichi2_spheroid.jl");
include("oichi2_binary.jl");
include("fused_polyft.jl");
include("shape_gradient.jl");
include("oiplot_spheroid.jl");

# Re-export OITOOLS functions so users only need `using ROTIR`
export OIdata, readoifits, readoifits_multiepochs, readfits, writefits
export set_data_filter, filter_data
export plot_v2_residuals, plot_t3amp_residuals, plot_t3phi_residuals, plot_residuals

# Tessellation
export tessellation
export tessellation_healpix, tessellation_latlong
export nside2npix, npix2n
export healpix_round_star, healpix_ellipsoid_star, tv_neighbors_healpix, tv_neighbors_healpix_visible
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
export solve_R_L1, solve_R_L2, solve_R_L3, solve_lagrange_points
export newton_raphson, brent_root
export radius_eggleton, radius_leahy, rpole_to_fillout
export roche_volume, roche_area, roche_equivalent_radius, romberg_integrate

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
export observables, cvis_to_obs, chi2s, mod360
export spheroid_chi2_f, spheroid_chi2_fg
export spheroid_chi2_allepochs_fg, spheroid_chi2_allepochs_f
export spheroid_total_variation, spheroid_crit_multiepochs_fg
export spheroid_l2_fg, spheroid_harmon_bias_fg, spheroid_regularization
export image_reconstruct_oi, image_reconstruct_oi_crit, image_reconstruct_oi_chi2, image_reconstruct_oi_chi2_fg
export multires_reconstruct_oi
export cvis_to_v2, poly_to_cvis, poly_to_flux, cvis_to_t3

# Binary forward model
export binary_phase_shift, binary_cvis, binary_observables, binary_chi2_f, orbit_to_rotir_offset

# Soft visibility
export sigmoid, dsigmoid, soft_visibility

# Fused two-pass polyft (matrix-free forward/adjoint)
export compute_polyflux_and_cvis!, compute_adjoint_cvis!, compute_adjoint_vertices!
export precompute_k2_inv_im, fused_spheroid_chi2_fg

# Shape gradients (joint shape + map optimization)
export rotation_matrix, dR_dinc, dR_dPA
export projected_vertices_and_derivs, shape_chi2_fg!, joint_reconstruct_oi

# Plotting
export plot2d, plot2d_wireframe, plot2d_allepochs
export plot3d
export plot_mollweide
export draw_compass, draw_rotation_axis, draw_rotation_arrow, draw_graticules
export plot_rv, plot2d_binary
export sometimes_visible, never_visible, invisible_neighbors, with_invisible_neighbors, without_invisible_neighbors

# Utilities
export make_circ_spot, make_spot_move
export rl1, max_rpole
export lci_linear_inversion_frame, lci_reconstruct_mutitemporal, rescale_temperature

include("precompile.jl")

end
