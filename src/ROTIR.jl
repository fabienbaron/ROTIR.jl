module ROTIR
using OITOOLS
using Statistics
include("geometry.jl");
include("oichi2_spheroid.jl");
include("oiplot_spheroid.jl");
include("oistars.jl");
include("geometry_rochelobe.jl");

export starparameters
export binaryparameters
export update_star
export update_binary
#export pix2vec_nest, pix2ang_nes, nside2npix, vec2ang

export readoifits_multiepochs
export oblate_const
export latlong_rochelobe
export healpix_round_star, tv_neighbours_healpix
export latlong_ellipsoid_star, tv_neighbors_longlat
export create_geometry

export calc_tempmap_vZ
export compute_teff_vonzeipel

export setup_polygon_ft
export observables, chi2s
export spheroid_chi2_f,spheroid_chi2_fg,spheroid_chi2_fg_alt,spheroid_chi2_allepochs_fg,spheroid_chi2_allepochs_f,spheroid_total_variation,spheroid_crit_multiepochs_fg,spheroid_l2_fg,spheroid_harmon_bias_fg,spheroid_regularization,spheroid_oi_reconstruct
export plot2d_temperature,plot2d_temperature_allepochs
export plot2d_intensity,plot2d_intensity_allepochs
export sometimes_visible, never_visible
export mollplot_temperature_longlat, mollplot_temperature_healpix
export plot3d_temperature
export plot2d_temperature_savefig
export plot3d_temperature_binary
export nside2npix, npix2n # Healpix
export lci_linear_inversion_frame,lci_reconstruct_mutitemporal,rescale_temperature

include("lci.jl") # Light curve inversion main functions
export lci_reconstruct, read_lci_relative,read_lci_absolute,read_lci_absolute_mjd,write_lci,modelflux_lci,modelflux_lci_rel,setup_lci, chi2_lci_f, chi2_lci_bias_fg, setup_stellar_parameters_lci, split_lcidata_by_period,setup_regularization_matrices
include("lciplot.jl") # Light curve inversion main functions
export lciplot_phase, lciplot_mjd, lciplot_vs_model_phase,lciplot_vs_model_mjd, lciplot_vs_model
end
