# Non-reversible parallel tempering for ROTIR's parametric fits, through Pigeons.jl.
#
#     using ROTIR, Pigeons, Distributions, LogDensityProblems, ADTypes, Zygote
#     fit_sphere_ld(data, tessels; method = :pigeons)
#
# Five triggers rather than one, because the backend uses all of them directly: Pigeons for the
# sampler, LogDensityProblems for the target interface, Distributions and ADTypes for the
# reference distribution and the explorer's AD kind, Zygote for the gradient. All four of the
# others are Pigeons' own dependencies, so `using Pigeons` installs nothing extra and the
# import line is the only cost.
#
# Why this sampler exists beside three others: it is the one that handles a MULTIMODAL
# posterior, which is what ROTIR's χ² is on a real star. See the header of src/fit_pigeons.jl.
module ROTIRPigeonsExt

using ROTIR
using Pigeons
using Distributions, LogDensityProblems, ADTypes, Zygote
using Printf, Statistics, LinearAlgebra, Random

# `_fit_pigeons` is declared in the core package (so `method = :pigeons` can name it before
# Pigeons is loaded and give a useful error if it never is); this file provides the method.
import ROTIR: _fit_pigeons

# The box transform is shared with the NUTS backend rather than copied: the two must agree
# about what the unconstrained coordinates MEAN, or a posterior from one is not comparable
# with a posterior from the other. It lives in src/fit_hmc.jl, which belongs to a sibling
# extension — so the parent declares it and both extensions extend it.
using ROTIR: _box_transform, build_parametric_logπ, default_parametric_bounds,
             parametric_free_indices

include(joinpath(pkgdir(ROTIR), "src", "fit_pigeons.jl"))

end # module
