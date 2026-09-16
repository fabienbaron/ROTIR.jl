# The FINUFFT type-3 kernel, switched on by `using FINUFFT`.
#
#     using ROTIR, FINUFFT
#     ROTIR.POLYFT_BACKEND[] = :nufft
#
# AN EXTENSION RATHER THAN A DEPENDENCY, and unlike LoopVectorization the reason is both load
# cost and SIZE. MEASURED: `using FINUFFT` is 464 ms on every `using ROTIR`, whether or not
# anything transforms; and it depends on `cufinufft_jll`, which pulls
# CUDA_Runtime_jll -> CUDA_Driver_jll and puts 455 MB of NVIDIA driver into an application
# bundle — 14 % of it — for a GPU NUFFT the application never calls. That artifact is not
# marked lazy, so `include_lazy_artifacts = false` cannot exclude it, and deleting it
# post-build kills startup inside `CUDA_Driver_jll.__init__`.
#
# WHY `:nufft` IS KEPT AT ALL, given `:t3` is faster at every mesh level and precision and
# more accurate on both test datasets. It is an INDEPENDENT implementation of the same idea by
# a mature library, and `src/type3_nufft.jl` is new: four real bugs turned up in it on the day
# it was written — a missing grid floor, a Gauss order the table does not have, an adjoint
# reading uninitialised memory, and a span calibration sitting on a boundary — and three of
# the four were found by exercising geometries its own benchmarks never produced. The exact
# kernels are the better reference where they are affordable, but at HEALPix 6 they cost
# 19.7 ms forward and 130 ms for the scalar adjoint, so this is the only cheap cross-check
# left at the meshes that matter most.
module ROTIRFINUFFTExt

using ROTIR
using FINUFFT

# Declared in the core package so `:nufft` can be NAMED — and refused with a message saying
# which import is missing — in a session that never loads FINUFFT.
import ROTIR: _finufft2d3

"""
    _finufft2d3(xs, ys, fs, iflag, tol, sk, tk)

The 2-D type-3 transform itself. A one-line shim on purpose: every choice ROTIR makes about
the quadrature, the coordinate convention and the tolerance stays in `src/polyft_nfft.jl`, so
what lives behind the extension boundary is only the C call.
"""
_finufft2d3(xs, ys, fs, iflag, tol, sk, tk) =
    FINUFFT.nufft2d3(xs, ys, fs, iflag, tol, sk, tk)

# ANNOUNCE OURSELVES rather than let `finufft_available` go looking, for the reason
# ROTIRLoopVectorizationExt does the same: `__init__` runs when this extension loads, which is
# exactly the moment the answer changes — including a mid-session load from the GUI.
function __init__()
    ROTIR.FINUFFT_OK[] = true
    return nothing
end

end # module
