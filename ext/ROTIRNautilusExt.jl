# The Nautilus.jl sampler for ROTIR's parametric fits.
#
#     using ROTIR, Nautilus
#     fit_sphere_ld(data, tessels; method = :nautilus)
#
# Pure Julia, unlike `:ultranest`, so this is the nested sampler a session with no conda
# environment — and a PackageCompiler build, which cannot contain PythonCall — can still use.
#
# Structured like ROTIRZygoteExt: the implementation lives in src/fit_nautilus.jl and knows
# nothing about backend dispatch, and the two lines below are the whole connection.
module ROTIRNautilusExt

using ROTIR
using Nautilus
using Printf, Statistics, Random

# `_fit_nautilus` is declared in the core package (so `method = :nautilus` can name it before
# Nautilus is loaded and give a useful error if it never is); this file provides the method.
import ROTIR: _fit_nautilus

include(joinpath(pkgdir(ROTIR), "src", "fit_nautilus.jl"))

end # module
