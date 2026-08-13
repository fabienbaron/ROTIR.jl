module ROTIR
import OITOOLS          # for OITOOLS.set_oiplot_defaults — see src/oiplot_spheroid.jl
import OITOOLS: OIdata,
    readoifits, readoifits_multiepochs, readfits, writefits,
    set_data_filter, filter_data,
    plot_v2_residuals, plot_t3amp_residuals, plot_t3phi_residuals, plot_residuals,
    DataBlocks, data_blocks, resample_blocks, block_counts, block_weights,
    apply_block_counts, apply_block_weights, perturb_data,
    BootstrapResult, bootstrap_driver,
    visibility_ud, visibility_ldlin, visibility_ldquad, visibility_ldpow,
    visibility_Gaussian,
    cvis_to_chi2_f, cvis_to_chi2_fg
using Statistics
using LinearAlgebra
using NLopt
using Printf
using PrecompileTools
using ChainRulesCore

include("oistars.jl");
include("soft_visibility.jl");

# Convert all AbstractFloat fields of a NamedTuple to type T in one shot.
# Integer fields (surface_type, ldtype) pass through unchanged.
convert_params(::Type{T}, p::NamedTuple) where T = NamedTuple{keys(p)}(map(v -> v isa AbstractFloat ? T(v) : v, values(p)))

include("intensity.jl");     # before the chi2 paths: they all map T -> band intensity
include("geometry.jl");
include("binary_geometry.jl");
include("reflection.jl");
include("oichi2_spheroid.jl");
include("oichi2_binary.jl");
include("fused_polyft.jl");
include("shape_gradient.jl");
include("parametric_gradient.jl");
include("bootstrap.jl");
include("parametric_fit.jl");
include("orbit_fit.jl");
include("ultranest.jl");
include("rasterize.jl");
include("polyft_nfft.jl");
include("oiplot_spheroid.jl");
include("animation.jl");

# Re-export OITOOLS functions so users only need `using ROTIR`
export OIdata, readoifits, readoifits_multiepochs, readfits, writefits
export set_data_filter, filter_data
export plot_v2_residuals, plot_t3amp_residuals, plot_t3phi_residuals, plot_residuals
export DataBlocks, data_blocks, resample_blocks, block_counts, block_weights
export apply_block_counts, apply_block_weights, perturb_data
export BootstrapResult, bootstrap_driver
# Analytic component visibilities — useful for quick parametric fits (e.g. per-epoch
# binary astrometry) alongside ROTIR's tessellated surface models.
export visibility_ud, visibility_ldlin, visibility_ldquad, visibility_ldpow
export visibility_Gaussian

# Tessellation
export tessellation
export tessellation_healpix, tessellation_latlong
export nside2npix, npix2n
export tv_neighbors_healpix, tv_neighbors_healpix_visible, tv_neighbors_longlat
export upsample_map_stars, downsample_map_stars

# Stellar/binary parameters
export starparameters, binaryparameters

# Geometry: stars and binaries
export stellar_geometry
export create_star, create_star_multiepochs
export compute_separation
export oblate_const

# Geometry: two stars in one frame (binary_geometry.jl)
export binary_frame, sky_of_orbit
export create_binary_star, create_binary_geometry, create_binary_geometry_multiepochs
export projected_separation, check_binary_overlap
export occultation_weights, projected_limb_profile, limb_radius, silhouette_polygon, convex_hull_2d
export polygon_convex_clip_area
export omega_at
export roche_omega_for_volume, roche_omega_table, roche_mesh_volume, tessel_solid_angles
export finish_star

# Mutual irradiation / reflection effect (reflection.jl)
export tessel_centroids_areas, ld_bol_D0, ld_bol, crossbody_kernels
export reflection_kernels, solve_radiosity, handle_reflection

# Temperature -> band intensity (intensity.jl)
export intensity, planck_and_dT, band_of

# Geometry: Roche lobe
export update_roche_radii, get_surface_potential, synchronicity
export compute_potential_primary, compute_potential_secondary, solve_radius
export solve_R_L1, solve_R_L2, solve_R_L3, solve_lagrange_points
export brent_root, roche_polar_radius
export radius_eggleton, radius_leahy, rpole_to_fillout
export roche_volume, roche_area, roche_equivalent_radius, romberg_integrate

# Geometry: rapid rotators
export temperature_map_vonZeipel_rapid_rotator
export calc_rotspin, calc_omega

# Geometry: temperature maps
export temperature_map_vonZeipel_roche_single, compute_gravity_primary, compute_gravity_secondary
export temperature_map_vonZeipel_ellipsoid
export parametric_temperature_map, spheroid_parametric_f

# Orbits
export compute_eccentric_anomaly, compute_true_anomaly
export compute_E_NR, kepler_E, kepler_E_vec, compute_coeff, compute_xyz_rel
export binary_orbit_rel, binary_orbit_abs, binary_RV, binary_proj_plane

# OI chi2 and reconstruction
export setup_oi!, setup_polygon_ft, setup_polyflux_single, setup_polyft_single, setup_polyft_single_alt
export observables, cvis_to_obs, cvis_chi2, OI_DEFAULT_WEIGHTS, chi2s, mod360
export spheroid_chi2_f, spheroid_chi2_fg
export spheroid_chi2_allepochs_fg, spheroid_chi2_allepochs_f
export spheroid_total_variation, spheroid_crit_multiepochs_fg
export spheroid_radflat_fg, radflat_bins
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

# Rasterization (polygon -> image via Sutherland-Hodgman clipping)
export rasterize_polygon_image!, rasterize_polygon_image, rasterize_adjoint!

# NFFT (polygon -> Fourier grid via Gauss-Legendre quadrature + NFFT)
export build_gauss_samples, polyft_nfft_forward, polyft_nfft_image

# Shape gradients (joint shape + map optimization)
export rotation_matrix, dR_dinc, dR_dPA
export projected_vertices_and_derivs, shape_chi2_fg!, joint_reconstruct_oi

# Parametric gradient (Zygote-composable primitives + ChainRules rrules)
export limb_mu, mu_and_dmu
export vonzeipel_map, vonzeipel_map_and_derivs
export ld_weight, ld_and_derivs, visibility_weight
export project_geometry, interferometric_chi2, build_parametric_logπ

# Bootstrap uncertainties for parametric fits (fit_parametric needs `using Zygote`)
export epoch_blocks, resample_epochs, bootstrap_parametric, ParametricBootstrap
export fit_parametric, default_parametric_bounds, parametric_param_names
export parametric_chi2, fit_sphere_ld, fit_ellipsoid_ld
# Generic orbit fitting from OIFITS (orbit_fit.jl)
export OrbitComponent, PointSource, UniformDisk, LimbDarkenedDisk, GaussianDisk,
       EllipticalGaussian
export OrbitFitData, OrbitFitSpec, orbit_fit_data, orbit_fit_spec
export orbit_model_cvis, orbit_chi2, fit_orbit
export ORBIT_ELEMENTS

export parametric_free_indices
export fit_parametric_ultranest

# Plotting
export plot2d, plot2d_wireframe, plot2d_allepochs
export plot3d
export plot_mollweide
export draw_compass, draw_rotation_axis, draw_rotation_arrow, draw_graticules, draw_limb!
export plot_rv, plot2d_binary, add_tessel_collection!

# Animation
export binary_movie, binary_frame_maps, frames_to_movie
export sometimes_visible, never_visible, invisible_neighbors, with_invisible_neighbors, without_invisible_neighbors

# Utilities
export make_circ_spot, make_spot_move
export rl1, max_rpole
export rescale_temperature_tpole, rescale_temperature_teff

function __init__()
    # ROTIR extracts parallelism at the Julia level (threaded polygon-FT kernels in
    # fused_polyft.jl; MCMC over chains/epochs). Keep FFTW/NFFT single-threaded so their
    # internal FFTs do NOT oversubscribe against Julia's threads under `julia -t auto`
    # (NFFT.__init__ otherwise sets _use_threads = nthreads()>1). Re-enable after
    # `using ROTIR` if you drive FFTW/NFFT yourself in a single-threaded outer context:
    #     FFTW.set_num_threads(Threads.nthreads()); NFFT._use_threads[] = true
    try
        FFTW.set_num_threads(1)
        NFFT._use_threads[] = false
    catch err
        @debug "ROTIR: could not force single-threaded FFTW/NFFT" err
    end
end

include("precompile.jl")

end
