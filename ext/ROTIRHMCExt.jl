# NUTS for the parametric model.
#
#     using ROTIR, Zygote, LogDensityProblems, AdvancedHMC
#     fit_sphere_ld(...; method = :hmc)     # or the GUI's Model tab
#
# The sampler that uses the analytic gradient ROTIR already carries. Structured like
# ext/ROTIRNautilusExt.jl: the implementation is in src/fit_hmc.jl and knows nothing about
# backend dispatch.
#
# Zygote is a trigger as well as AdvancedHMC: the gradient comes from `Zygote.withgradient`,
# and without it this would load and then fail on the first draw.
module ROTIRHMCExt

using ROTIR
using AdvancedHMC, LogDensityProblems, Zygote
using LinearAlgebra, Statistics, Printf

import ROTIR: _fit_hmc
using ROTIR: build_parametric_logπ, default_parametric_bounds, parametric_free_indices,
             _box_transform

include(joinpath(pkgdir(ROTIR), "src", "fit_hmc.jl"))

end # module
