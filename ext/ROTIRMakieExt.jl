# Plotting in Makie, without Qt and without Python.
#
# Loaded when the caller has Makie — which any backend brings: `using CairoMakie` gets vector
# figures with no GPU and no display, `using GLMakie` gets a window, and a GUI gets it through
# QMLMakie. Choosing the backend is the caller's job and this file never names one.
#
# Everything here builds a Figure and hands it back, so none of it needs Qt. That is what will
# let a script — or a PackageCompiler build, which cannot contain PythonPlot — plot at all.
# An eventual GUI extension keeps only what is genuinely interactive: plots created once and
# then driven by Observables, because inserting into a live scene under QMLMakie allocates GPU
# buffers with no GL context bound.
module ROTIRMakieExt

using ROTIR
using Makie
using Statistics, LinearAlgebra, Printf

import ROTIR: plot2d_makie, plot3d_makie, plot2d_wireframe_makie,
              plot2d_binary_makie, plot2d_allepochs_makie, plot_rv_makie,
              plot_mollweide_makie, mollweide_xy, mollweide_grid,
              add_tessel_collection_makie!, draw_limb_makie!, draw_compass_makie!,
              draw_graticules_makie!, draw_rotation_axis_makie!,
              draw_rotation_arrow_makie!, plot3d_binary_makie, star_mesh, scene3d,
              relative_orbit_track, tessel_polygons, map_colors, sky_axis_max,
              style_sky_axis!, _padded_cmap

# Backend-agnostic mesh geometry (src/oiplot_spheroid_core.jl), shared with the matplotlib
# layer so the two cannot disagree about what they are drawing. Not exported by the parent.
using ROTIR: graticule_segments, _mesh_rotation, _spin_axis, _polar_radius,
             _mesh_surface_field, _interp_surface, _mesh_body_curve, _visible_segments,
             _latlong_dims, _map_range, _body_radius,
             convex_hull_2d, npix2nside, ang2pix_nest

include(joinpath(pkgdir(ROTIR), "src", "oiplot_spheroid_makie.jl"))

end # module
