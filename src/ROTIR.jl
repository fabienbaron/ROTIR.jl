module ROTIR
using OITOOLS
using Statistics
using LinearAlgebra

include("oistars.jl");
include("geometry.jl");
include("oichi2_spheroid.jl");
include("oiplot_spheroid.jl");

#include("di.jl");
#include("di_dynamical.jl");
#include("reconstruct.jl");

export tessellation
export tessellation_healpix, tessellation_latlong
export update_roche_radii, get_surface_potential, compute_potential_primary, compute_potential_secondary,solve_radius, update_roche_geom
export radius_eggleton, radius_leahy
export compute_true_anomaly, compute_E_NR,compute_coeff, compute_xyz_rel # debug 
export compute_dpotential_primary_L1,  solve_R_L1, newton_raphson
export rotate_single_star, compute_separation
export temperature_map_vonZeipel_roche_single, compute_gravity_primary, compute_gravity_secondary
export temperature_map_vonZeipel_ellipsoid
# Export stellar/binary parameters for initializing and updating
export starparameters, binaryparameters
export compute_eccentric_anomaly, compute_true_anomaly
export update_star, update_binary
export binary_orbit_rel, binary_orbit_abs, binary_RV, binary_proj_plane

# Export reading in oifits files, geometry initialization, and regularization
export readoifits_multiepochs, readoifits #, OIdata, LCI
export oblate_const
export healpix_round_star, healpix_ellipsoid_star, tv_neighbours_healpix
export latlong_round_star, latlong_ellipsoid_star, tv_neighbors_longlat, latlong_rochelobe
export create_star, create_star_multiepochs, create_binary

# Export for rapid rotators
export temperature_map_vonZeipel_rapid_rotator, calc_grelmap_vZ
export compute_teff_vonzeipel
export calc_rotspin, calc_omega

# Export criterion calculations for spheroids
export setup_oi!, setup_polygon_ft, observables, chi2s, mod360
export spheroid_chi2_f,spheroid_chi2_fg,spheroid_chi2_fg_alt,spheroid_chi2_allepochs_fg,spheroid_chi2_allepochs_f
export spheroid_total_variation,spheroid_crit_multiepochs_fg,spheroid_l2_fg,spheroid_harmon_bias_fg,spheroid_regularization
export image_reconstruct_oi, image_reconstruct_oi_crit, image_reconstruct_oi_chi2
export parametric_temperature_map, spheroid_parametric_f
# Export plotting routines
export plot2d_wireframe, plot2d_allepochs
export sometimes_visible, never_visible, invisible_neighbours
export plot_mollweide
export plot3d, plot3d_vertices, plot2d

export nside2npix, npix2n # Healpix
export lci_linear_inversion_frame,lci_reconstruct_mutitemporal,rescale_temperature
export readfits, writefits
# # Export Doppler imaging routines
# export modelGrid, globalLineProfiles, observedLineProfiles, jd2phase
# export setup_di, read_di_profiles, bin_spectrum, setup_stellar_parameters_di, di_reconstruct
# export make_circ_spot_DI, make_meridional_band_DI, make_equatorial_band_DI
# export calculate_global_profiles, calculate_global_profiles_unnormalized
# export compute_grid_atlas9_SPECTRUM, read_spectrum_grid, setup_model_grid
# export gaussian_instrumental_kernel, instrumental_broadening, doppler_shift_spectrum
# export calculate_chi2_f
# export plot_profiles, mollplot_healpix, mollplot_longlat
# export make_circ_spot_DI_spotfill, calculate_global_profiles_spotfill, di_reconstruct_spotfill
# export calculate_μv

# #Export Dynamical Doppler routines
# export split_didata_by_period, di_reconstruct_multitemporal
# export evolve_frame

# # Wrapped for LCI, DI, OI reconstructions
# export reconstruct

# # Misc. Spectrum utilities
# include("specutils.jl")
# export normalize_UVES, fit_rv, fit_continuum

# testing
export cvis_to_v2, poly_to_cvis, poly_to_flux, cvis_to_t3
export make_circ_spot, make_spot_move


export rl1, max_rpole

#include("lci.jl") # Light curve inversion main functions
#export LCI
#export lci_reconstruct, read_lci_relative,read_lci_absolute,read_lci_absolute_mjd,write_lci,modelflux_lci,modelflux_lci_rel
#export setup_lci, chi2_lci_f, chi2_lci_bias_fg, setup_stellar_parameters_lci, split_lcidata_by_period, setup_regularization_matrices
#include("lciplot.jl") # Light curve inversion main functions
#export lciplot_phase, lciplot_mjd, lciplot_vs_model_phase,lciplot_vs_model_mjd, lciplot_vs_model, lciplot_vs_model_phase_residuals
end
