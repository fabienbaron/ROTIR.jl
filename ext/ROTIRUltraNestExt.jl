# The UltraNest sampler for ROTIR's parametric fits.
#
#     using ROTIR, PythonCall
#     fit_sphere_ld(data, tessels; method = :ultranest)
#
# Triggered by PythonCall rather than by an `UltraNest` package, because UltraNest is a PYTHON
# package: there is no Julia module to name, and PythonCall is the thing whose presence
# actually decides whether this can work.
#
# It is an extension at all — rather than core, which it was — because `using PythonCall`
# costs 1.2 s on every GUI start by invalidating OITOOLS' precompiled plot pipeline. See the
# header of src/fit_ultranest.jl for the measurement.
module ROTIRUltraNestExt

using ROTIR
using PythonCall
using Printf, Statistics

# Both are declared in the core package, so `method = :ultranest` and
# `fit_parametric_ultranest` can be NAMED — and refused with a message saying which import is
# missing — in a session that never loads PythonCall. These files provide the methods.
import ROTIR: _fit_ultranest, fit_parametric_ultranest

# Everything these two need from the core package. They were core files until this extension
# existed, so they use the unexported internals directly.
using ROTIR: build_parametric_logπ, default_parametric_bounds, parametric_free_indices

include(joinpath(pkgdir(ROTIR), "src", "fit_ultranest.jl"))
include(joinpath(pkgdir(ROTIR), "src", "ultranest.jl"))

end # module
