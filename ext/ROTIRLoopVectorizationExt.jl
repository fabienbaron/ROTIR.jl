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
import ROTIR: _cvis_turbo!, _adj_cvis_turbo!, _adj_vertices_turbo!

include(joinpath(pkgdir(ROTIR), "src", "turbo_polyft.jl"))

# ANNOUNCE OURSELVES rather than let `turbo_available` go looking. It used to answer by asking
# the runtime for `methods(_cvis_turbo!)`, a reflection call in the innermost visibility
# kernel; setting a flag here is both cheaper and correct at exactly the right moment, since
# `__init__` runs when this extension loads — including when the GUI loads LoopVectorization
# on demand in the middle of a session.
function __init__()
    ROTIR.TURBO_OK[] = true
    return nothing
end

end # module
