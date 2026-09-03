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
import FINUFFT
import FITSIO
import FITSIO: read_header
import Dates

include("oistars.jl");
include("surface_schema.jl");
include("surface_map_io.jl");
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
include("orbit_ties.jl");    # before orbit_fit.jl: OrbitFitSpec holds an OrbitTies
include("orbit_fit.jl");
# src/ultranest.jl and src/fit_ultranest.jl are NOT included here: they need PythonCall and
# live in ext/ROTIRUltraNestExt.jl.
include("rasterize.jl");
include("polyft_nfft.jl");
include("oiplot_spheroid_core.jl");  # mesh geometry shared by both plotting back-ends
include("animation.jl");             # pure: frame maps, value ranges, ffmpeg driver

# ── Plotting: declared here, given methods by an extension ───────────────────
#
# The drawing code itself lives in src/oiplot_spheroid.jl (matplotlib) and
# src/animation_movie.jl, loaded by ROTIRPythonPlotExt when the caller has PythonPlot. A
# Makie implementation attaches to these same names from ROTIRMakieExt.
#
# Two reasons this is not in the core load path, both borrowed from OITOOLS:
#
#  1. `using ROTIR` imported matplotlib whether or not anything was ever plotted — a cost
#     paid by every reduction script and CI run.
#  2. It is a hard blocker for a QML GUI. matplotlib probes for an interactive backend, finds
#     PySide6 in the conda environment and maps Qt 6.11 into the process; QML.jl needs its own
#     Qt to be the only one, and the collision surfaces as
#     `libQt6DBus.so: undefined symbol: _ZN14QObjectPrivateC2E16QtPrivate_6_10_2`.
#
# PythonCall is a WEAK dependency too, and for a reason measured rather than assumed: loading
# it invalidates the plot-construction code OITOOLS precompiles for its live canvas, taking
# `build_canvas` from 338 ms to 2477 ms. That 1.2 s landed on every ROTIR GUI start, sampler
# or no sampler. UltraNest is therefore in ext/ROTIRUltraNestExt.jl, which PythonCall
# triggers; `:nautilus` is the pure-Julia sampler for anyone who wants one without Python.
function plot2d end
function plot2d_binary end
function plot2d_wireframe end
function plot2d_allepochs end
function plot3d end
function plot_mollweide end
function plot_rv end
function draw_compass end
function draw_rotation_axis end
function draw_rotation_arrow end
function draw_graticules end
function draw_limb! end
function add_tessel_collection! end
function binary_movie end
# Not exported, but reached as `ROTIR.name` by demos and tests, so they must be declared in
# this module rather than left local to the extension.
function set_oiplot_defaults end
function plot2Dquad end
function set_tick_spacing end

# ── The Makie drawing layer (ROTIRMakieExt) ─────────────────────────────────
#
# Declared here for the same reason as the matplotlib stubs, and additionally because a GUI
# extension will call them: one extension cannot reach into another, so both attach to
# functions this module owns. `using Makie` — which every Makie backend brings — gives them
# methods. The `_makie` suffix exists because both back-ends can be loaded at once; it goes
# away when matplotlib does.
function plot2d_makie end
function plot3d_makie end
function plot2d_wireframe_makie end
function plot2d_binary_makie end
function plot2d_allepochs_makie end
function plot_rv_makie end
function plot_mollweide_makie end
function mollweide_xy end
function mollweide_grid end
function add_tessel_collection_makie! end
function draw_limb_makie! end
function draw_compass_makie! end
function draw_graticules_makie! end
function draw_rotation_axis_makie! end
function draw_rotation_arrow_makie! end
function plot3d_binary_makie end
# Not exported: internal to the plotting layers, but declared here so the Makie
# extension extends it and the GUI extension can reach it through ROTIR.
function _padded_cmap end
# The Nautilus sampler, implemented in ext/ROTIRNautilusExt.jl. Declared here so that
# `method = :nautilus` can be named — and refused with a message that says WHICH package is
# missing — in a session that never loads it.
function _fit_nautilus end
# NUTS, implemented in ext/ROTIRHMCExt.jl. Declared here so `method = :hmc` can be named — and
# refused with a message naming the packages — in a session that never loads them.
function _fit_hmc end
# UltraNest, implemented in ext/ROTIRUltraNestExt.jl. Declared here for the same reason as the
# other two: so `method = :ultranest` is a name this package knows, and the failure without
# PythonCall is a sentence rather than a MethodError on an underscored internal.
function _fit_ultranest end
function fit_parametric_ultranest end
# Non-reversible parallel tempering, implemented in ext/ROTIRPigeonsExt.jl. The sampler for a
# MULTIMODAL posterior, which is what a real ROTIR χ² is; declared here for the same reason as
# the other three.
function _fit_pigeons end

"""
    hmc_available() -> Bool

Whether `method = :hmc` will work here, i.e. whether `using AdvancedHMC, LogDensityProblems,
Zygote` has loaded ROTIRHMCExt.
"""
hmc_available() = !isempty(methods(_fit_hmc))

"""
    nautilus_available() -> Bool

Whether `method = :nautilus` will work in this session, i.e. whether `using Nautilus` has
loaded ROTIRNautilusExt. Checked rather than assumed because the failure otherwise appears as
a `MethodError` on an underscore-prefixed internal, which says nothing about what to install.
"""
nautilus_available() = !isempty(methods(_fit_nautilus))

"""
    ultranest_available() -> Bool

Whether `method = :ultranest` will work here, i.e. whether `using PythonCall` has loaded
ROTIRUltraNestExt. The Python `ultranest` package itself is checked separately, when the
sampler runs — this only says whether the bridge to Python exists at all.
"""
ultranest_available() = !isempty(methods(_fit_ultranest))

"""
    pigeons_available() -> Bool

Whether `method = :pigeons` will work here, i.e. whether
`using Pigeons, Distributions, LogDensityProblems, ADTypes, Zygote` has loaded
ROTIRPigeonsExt. The four besides Pigeons are Pigeons' own dependencies, so nothing extra is
installed by asking for them.
"""
pigeons_available() = !isempty(methods(_fit_pigeons))
# The GUI's entry point. Declared here so `gui()` is a name `using ROTIR` already knows about
# and the error for calling it without the stack is Julia's "no method", naming the argument
# types, rather than "not defined" — which reads as a missing feature instead of a missing
# `using GLMakie, QMLMakie, QML`.
function gui end
function star_mesh end
function scene3d end
function relative_orbit_track end
function tessel_polygons end
function map_colors end
function sky_axis_max end
function style_sky_axis! end
function rgba end

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
export sobel_gradient_healpix

export upsample_map_stars, downsample_map_stars

# Stellar/binary parameters
export starparameters, binaryparameters

# What each surface_type requires — the single declaration the validator and any
# form-generating front-end both read (src/surface_schema.jl).
export ParamSpec, SurfaceSpec, SURFACE_TYPES, SURFACE_TYPE_ORDER,
       surface_spec, surface_params, default_star_params, validate_star_params,
       ld_coefficients_used

# A surface map as a file, with the tessellation and the parameters that reproduce its χ²
# (src/surface_map_io.jl).
export save_surface_map, load_surface_map

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
export observables, cvis_to_obs, cvis_chi2, OI_DEFAULT_WEIGHTS, chi2s, chi2_breakdown, mod360
export spheroid_chi2_f, spheroid_chi2_fg
# NOTE four names were exported here without ever being defined —
# spheroid_chi2_allepochs_fg, spheroid_crit_multiepochs_fg, spheroid_total_variation and
# spheroid_l2_fg. Exporting an undefined name is legal, so `using ROTIR` was silent, but any
# code enumerating `names(ROTIR)` to build a menu hit UndefVarError. The gradient-carrying
# all-epochs criterion is `spheroid_crit_allepochs_fg`; the TV regularizers are
# `spheroid_total_variation_fg` / `_variation2_fg`, reached through `spheroid_regularization`.
export spheroid_chi2_allepochs_f
export spheroid_radflat_fg, spheroid_radialvar_fg, radflat_bins,
       spheroid_orthold_fg, orthold_direction
export spheroid_sobel_fg, spheroid_sobel2_fg, max_entropy_fg
export spheroid_harmon_bias_fg, spheroid_regularization
export image_reconstruct_oi, image_reconstruct_oi_crit, image_reconstruct_oi_chi2, image_reconstruct_oi_chi2_fg
export multires_reconstruct_oi
export cvis_to_v2, poly_to_cvis, poly_to_flux, cvis_to_t3

# Binary forward model
export binary_phase_shift, binary_cvis, binary_observables, binary_chi2_f, orbit_to_rotir_offset

# Soft visibility
export sigmoid, dsigmoid, soft_visibility

# Fused two-pass polyft (matrix-free forward/adjoint)
export compute_polyflux_and_cvis!, compute_adjoint_cvis!, compute_adjoint_vertices!
export precompute_k2_inv_im, fused_spheroid_chi2_fg, fused_cvis, POLYFT_BACKEND
export turbo_available

# Rasterization (polygon -> image via Sutherland-Hodgman clipping)
export rasterize_polygon_image!, rasterize_polygon_image, rasterize_adjoint!

# NFFT (polygon -> Fourier grid via Gauss-Legendre quadrature + NFFT)
export build_gauss_samples, polyft_nfft_forward, polyft_nfft_image, polyft_cvis_nufft, quadrature_for

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
export OrbitTies, compile_ties, apply_ties, resolve_params

export parametric_free_indices
export fit_parametric_ultranest, ultranest_available, pigeons_available

# Plotting
export plot2d, plot2d_wireframe, plot2d_allepochs
export plot3d
export plot_mollweide
export draw_compass, draw_rotation_axis, draw_rotation_arrow, draw_graticules, draw_limb!
export plot_rv, plot2d_binary, add_tessel_collection!
# Makie counterparts (ROTIRMakieExt). Suffixed while both back-ends coexist.
export plot2d_makie, plot3d_makie, plot2d_wireframe_makie
export plot2d_binary_makie, plot2d_allepochs_makie, plot_rv_makie
export plot_mollweide_makie, mollweide_xy, mollweide_grid
export add_tessel_collection_makie!, draw_limb_makie!, draw_compass_makie!,
       draw_graticules_makie!, draw_rotation_axis_makie!, draw_rotation_arrow_makie!,
       plot3d_binary_makie, star_mesh, scene3d, relative_orbit_track, gui,
       nautilus_available, hmc_available,
       tessel_polygons, map_colors, sky_axis_max, style_sky_axis!

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
