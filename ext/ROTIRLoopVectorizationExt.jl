# The vectorised polygon-FT kernel, switched on by `using LoopVectorization`.
#
#     using ROTIR, LoopVectorization
#     ROTIR.POLYFT_BACKEND[] = :turbo
#
# An extension rather than a dependency because of the LOAD cost, not the run cost: the kernel
# is 17x the scalar reference, but loading LoopVectorization invalidates OITOOLS' precompiled
# plot pipeline and adds 1.8 s to every GUI start. With `:nufft` the default and faster than
# `:turbo` anyway, that is a bill for a cross-check. See src/turbo_polyft.jl.
module ROTIRLoopVectorizationExt

using ROTIR
using LoopVectorization

# Declared in the core package so `:turbo` can be NAMED — and refused with a message saying
# which import is missing — in a session that never loads LoopVectorization.
import ROTIR: _cvis_turbo!

include(joinpath(pkgdir(ROTIR), "src", "turbo_polyft.jl"))

end # module
