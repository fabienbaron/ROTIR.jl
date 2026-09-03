# matplotlib plotting for ROTIR.
#
# Loaded automatically when the caller has PythonPlot. Everything here used to sit in the core
# load path, which meant `using ROTIR` imported matplotlib whether or not anything was ever
# plotted — and, more seriously, made a QML GUI impossible (see the note beside the stub
# declarations in src/ROTIR.jl).
#
# `import`, not `using`: the definitions inside oiplot_spheroid.jl are written as plain
# `function plot2d(...)`, and importing the names makes those EXTEND ROTIR's functions rather
# than defining new ones local to this module. That is what lets both source files move here
# unchanged, and it is the same arrangement OITOOLSPythonPlotExt uses.
module ROTIRPythonPlotExt

using ROTIR
using PythonPlot, PythonCall, LaTeXStrings, Statistics, LinearAlgebra, Printf
# `set_oiplot_defaults` delegates to OITOOLS' so ROTIR and OITOOLS figures match in a
# paper. An extension may name any of its parent's dependencies, and OITOOLS is one.
import OITOOLS

import ROTIR: plot2d, plot2d_binary, plot2d_wireframe, plot2d_allepochs, plot3d,
              plot_mollweide, plot_rv, plot2Dquad, set_tick_spacing,
              draw_compass, draw_rotation_axis, draw_rotation_arrow, draw_graticules,
              draw_limb!, add_tessel_collection!, set_oiplot_defaults, rgba,
              binary_movie

# The backend-agnostic mesh geometry these draw from (src/oiplot_spheroid_core.jl). Not
# exported by the parent — they are shared implementation detail, so they are named here.
using ROTIR: _polar_radius, _mesh_rotation, _spin_axis, _mesh_surface_field,
             _interp_surface, _mesh_body_curve, _visible_segments, _latlong_dims,
             _binary_value_range, _binary_axis_max, graticule_segments,
             rot_vertex, npix2nside, ang2pix_nest, longlat_ang2pix

include(joinpath(pkgdir(ROTIR), "src", "oiplot_spheroid.jl"))
include(joinpath(pkgdir(ROTIR), "src", "animation_movie.jl"))

end # module
