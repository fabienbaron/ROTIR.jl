# polyft_nfft.jl — Fast polyft*x via NFFT at Gauss-Legendre quadrature points
#
# Instead of evaluating the continuous polygon Fourier transform via the
# closed-form edge formula at each (kx, ky) grid point, approximate each
# polygon integral by a few Gauss-Legendre quadrature points and fold the
# resulting non-uniform exponential sum into a single adjoint NFFT.
#
#     F[k] = sum_p x_weighted[p] * integral_polygon_p exp(-2*pi*i * k . r) dA
#          ~ sum_p sum_s  x_weighted[p] * J_s * w_s * exp(-2*pi*i * k . r_s)
#
# where r_s are Gauss points inside quad p, J_s is the Jacobian of the
# bilinear map [-1,1]^2 -> quad, and w_s is the Gauss weight.
# This sum is exactly an adjoint NFFT.
#
# Matches the frequency layout of `reshape(polyft*x_weighted, nx/2+1, nx)` so
# it can serve as a drop-in replacement for the dense `polyft*x` product when
# computing images on a regular grid.
#
# Adapted from planet_deconv for use with ROTIR's coordinate conventions
# (proj_west, proj_north in mas).
#
# This file only adds NEW functions -- it does not modify existing polyft code.

using NFFT
using FFTW

# Gauss-Legendre nodes and weights on [-1,1] for orders 2..10.
# Accuracy note: for oscillatory integrands exp(-2*pi*i * k . r), 2-point is
# only accurate when the phase change across the quad is < 1 rad.
# At n_hp=4: quads ~3 pixels wide => phase span ~10 rad at Nyquist => need >= 4.
# At n_hp=5: ~1.5 pixels => order 3 is usually sufficient.
# At n_hp=6: ~0.75 pixels => order 2 is sufficient.
const _GL_NODES = Dict{Int,Vector{Float64}}(
    2 => [-0.5773502691896257,  0.5773502691896257],
    3 => [-0.7745966692414834,  0.0,                 0.7745966692414834],
    4 => [-0.8611363115940526, -0.3399810435848563,  0.3399810435848563,  0.8611363115940526],
    5 => [-0.9061798459386640, -0.5384693101056831,  0.0,                 0.5384693101056831, 0.9061798459386640],
    6 => [-0.9324695142031521, -0.6612093864662645, -0.2386191860831969,  0.2386191860831969, 0.6612093864662645, 0.9324695142031521],
    8 => [-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
           0.1834346424956498,  0.5255324099163290,  0.7966664774136267,  0.9602898564975363],
   10 => [-0.9739065285171717, -0.8650633666889845, -0.6794095682990244, -0.4333953941292472, -0.1488743389816312,
           0.1488743389816312,  0.4333953941292472,  0.6794095682990244,  0.8650633666889845,  0.9739065285171717],
)
const _GL_WEIGHTS = Dict{Int,Vector{Float64}}(
    2 => [1.0,                 1.0],
    3 => [0.5555555555555556,  0.8888888888888888,  0.5555555555555556],
    4 => [0.3478548451374538,  0.6521451548625461,  0.6521451548625461,  0.3478548451374538],
    5 => [0.2369268850561891,  0.4786286704993665,  0.5688888888888889,  0.4786286704993665, 0.2369268850561891],
    6 => [0.1713244923791704,  0.3607615730481386,  0.4679139345726910,  0.4679139345726910, 0.3607615730481386, 0.1713244923791704],
    8 => [0.1012285362903763,  0.2223810344533745,  0.3137066458778873,  0.3626837833783620,
          0.3626837833783620,  0.3137066458778873,  0.2223810344533745,  0.1012285362903763],
   10 => [0.0666713443086881,  0.1494513491505806,  0.2190863625159820,  0.2692667193099963,  0.2955242247147529,
          0.2955242247147529,  0.2692667193099963,  0.2190863625159820,  0.1494513491505806,  0.0666713443086881],
)

"""
    build_gauss_samples(proj_west, proj_north, x_weighted; ngauss=4, nsub=1,
                        T=float(real(eltype(proj_west))))

Build Gauss-Legendre quadrature samples for each quadrilateral by subdividing
the reference square `[-1,1]^2` into `nsub x nsub` equal sub-squares and
applying an `ngauss x ngauss` tensor-product Gauss-Legendre rule in each.

Total samples per quad = `(ngauss * nsub)^2`.

The weight already includes the (signed) Jacobian of the bilinear
parameterization, the Gauss weight, and the `1/nsub^2` sub-square area factor.

Returns `(xs, ys, fs)` with length `(ngauss*nsub)^2 * Npix`.
"""
function build_gauss_samples(proj_west::AbstractMatrix, proj_north::AbstractMatrix,
                              x_weighted::AbstractVector;
                              ngauss::Int=4, nsub::Int=1,
                              T::Type = float(real(eltype(proj_west))))
    @assert haskey(_GL_NODES, ngauss) "ngauss must be one of $(sort(collect(keys(_GL_NODES))))"
    @assert nsub >= 1 "nsub must be >= 1"
    nodes_f64 = _GL_NODES[ngauss]
    w_f64     = _GL_WEIGHTS[ngauss]

    # Build effective node/weight list over the subdivided [-1,1] axis.
    nw = ngauss * nsub
    enodes = Vector{T}(undef, nw)
    eweights = Vector{T}(undef, nw)
    invnsub = T(1) / T(nsub)
    for j in 1:nsub
        cj = T(-1) + T(2j - 1) * invnsub
        for i in 1:ngauss
            idx = (j - 1) * ngauss + i
            enodes[idx]   = cj + T(nodes_f64[i]) * invnsub
            eweights[idx] = T(w_f64[i]) * invnsub
        end
    end

    Npix = size(proj_west, 1)
    npq  = nw * nw
    Ns   = npq * Npix
    xs   = Vector{T}(undef, Ns)
    ys   = Vector{T}(undef, Ns)
    fs   = Vector{Complex{T}}(undef, Ns)

    @inbounds for p in 1:Npix
        v1x = T(proj_west[p,1]);  v1y = T(proj_north[p,1])
        v2x = T(proj_west[p,2]);  v2y = T(proj_north[p,2])
        v3x = T(proj_west[p,3]);  v3y = T(proj_north[p,3])
        v4x = T(proj_west[p,4]);  v4y = T(proj_north[p,4])
        xw  = T(x_weighted[p])
        base = (p-1) * npq
        k = 0
        for i_eta in 1:nw
            eta  = enodes[i_eta]
            w_eta = eweights[i_eta]
            for i_xi in 1:nw
                xi  = enodes[i_xi]
                w_xi = eweights[i_xi]
                # Bilinear shape functions on [-1,1]^2
                N1 = (1-xi)*(1-eta)*T(0.25)
                N2 = (1+xi)*(1-eta)*T(0.25)
                N3 = (1+xi)*(1+eta)*T(0.25)
                N4 = (1-xi)*(1+eta)*T(0.25)
                x  = N1*v1x + N2*v2x + N3*v3x + N4*v4x
                y  = N1*v1y + N2*v2y + N3*v3y + N4*v4y
                # Jacobian of bilinear map (signed)
                dx_xi = T(0.25) * ((1-eta)*(v2x-v1x) + (1+eta)*(v3x-v4x))
                dy_xi = T(0.25) * ((1-eta)*(v2y-v1y) + (1+eta)*(v3y-v4y))
                dx_eta = T(0.25) * ((1-xi)*(v4x-v1x) + (1+xi)*(v3x-v2x))
                dy_eta = T(0.25) * ((1-xi)*(v4y-v1y) + (1+xi)*(v3y-v2y))
                J   = dx_xi*dy_eta - dx_eta*dy_xi
                k  += 1
                s   = base + k
                xs[s] = x; ys[s] = y
                fs[s] = Complex{T}(xw * J * w_xi * w_eta)
            end
        end
    end
    return xs, ys, fs
end

"""
    polyft_nfft_forward(proj_west, proj_north, x_weighted, pixsize, nx;
                        ngauss=4, nsub=1, T=float(real(eltype(proj_west))))

Compute `F[k] ~ sum_p x_weighted[p] * integral_polygon_p exp(-2*pi*i * k . r) dA`
via an adjoint NFFT at Gauss-Legendre samples inside each quadrilateral.

Returns a `(nx/2+1, nx)` complex array whose layout matches
`reshape(polyft*x_weighted, nx/2+1, nx)`, so that the caller can feed it
directly to `irfft(F, nx)` to obtain the real-space forward model.

Convention notes:
- NFFT.jl's adjoint computes `sum_j f_j exp(+2*pi*i * k . x_j)` in natural
  order `k in {-N/2, ..., N/2-1}`.
- With negated positions `pos = (-y/L, -x/L)`, the adjoint directly produces
  the standard DFT `sum f_j exp(-2*pi*i * k . r_j)`.
- `ifftshift` converts natural NFFT order to standard FFT order, then
  extracting rows `1:nx/2+1` gives the rfft layout.

# Arguments
- `proj_west, proj_north`: `(Npix, 4)` polygon vertex coordinates (mas)
- `x_weighted`: `(Npix,)` intensity weights
- `pixsize`: mas per pixel
- `nx`: grid size (pixels per side)
- `ngauss`: Gauss-Legendre order per axis (2..10)
- `nsub`: number of sub-squares per axis for subdivision
- `fftflags`: FFTW planner flags for the internal NFFT FFT. Default `FFTW.MEASURE`
  (benchmarks once per grid size, cached in FFTW wisdom → faster transforms when the
  NFFT is evaluated repeatedly). Pass `FFTW.ESTIMATE` for a single one-off call.
"""
function polyft_nfft_forward(proj_west, proj_north, x_weighted, pixsize::Real,
                              nx::Integer; ngauss::Int=4, nsub::Int=1,
                              T::Type = float(real(eltype(proj_west))), fftflags=FFTW.MEASURE)
    xs, ys, fs = build_gauss_samples(proj_west, proj_north, x_weighted;
                                      ngauss=ngauss, nsub=nsub, T=T)
    Ns = length(xs)
    L  = T(nx * pixsize)

    pos = Matrix{T}(undef, 2, Ns)
    @inbounds for s in 1:Ns
        pos[1, s] = -ys[s] / L
        pos[2, s] = -xs[s] / L
    end

    p_plan   = plan_nfft(pos, (nx, nx); fftflags=fftflags)
    fhat_nat = adjoint(p_plan) * fs

    # rfft-layout extraction: ifftshift converts natural → FFT order,
    # then rows 1:nh give the non-negative half (rfft layout).
    nh = nx ÷ 2 + 1
    return ifftshift(fhat_nat)[1:nh, :]
end

"""
    polyft_nfft_image(proj_west, proj_north, x_weighted, pixsize, nx;
                      ngauss=4, nsub=1, T=float(real(eltype(proj_west))))

Convenience function: compute the NFFT-based Fourier coefficients and
immediately inverse-FFT them to produce a real-space `(nx, nx)` image.

# Arguments
Same as `polyft_nfft_forward`.
"""
function polyft_nfft_image(proj_west, proj_north, x_weighted, pixsize::Real,
                            nx::Integer; ngauss::Int=4, nsub::Int=1,
                            T::Type = float(real(eltype(proj_west))), fftflags=FFTW.MEASURE)
    F = polyft_nfft_forward(proj_west, proj_north, x_weighted, pixsize, nx;
                             ngauss=ngauss, nsub=nsub, T=T, fftflags=fftflags)
    return fftshift(irfft(F, nx))
end


# ── Type-3: the polygon FT at SCATTERED uv, with no grid ────────────────────
#
# The functions above answer the image-deconvolution question: the polygon FT on a REGULAR
# Fourier grid, via an adjoint (type-1) NFFT. Interferometry asks a different one — the
# transform at 58 616 scattered uv points — and that is a type-3 transform, non-uniform
# sources to non-uniform targets.
#
# NFFT.jl does not have a fast one: `AbstractNFFTs` declares the two-node-set signature but
# NFFT.jl 0.14 implements only `NNDFTPlan`, the direct O(N·M) sum, which is what we are trying
# to avoid. FINUFFT's `nufft2d3` is a real type-3 and is what this uses.
#
# WHY THIS AND NOT THE RASTERISER. Rasterising the polygons onto a grid and transforming that
# was measured at 5–74× as well, and it reuses `rasterize_polygon_image` which is already
# here — but its answer disagrees with the exact polygon FT by ~5e-3 and, measured at
# nx = 128, 256 and 512, that disagreement does NOT shrink. It is not the transform: an exact
# DFT of the same image gives the same 4.97e-3, and the rasteriser conserves flux to 1.6e-4.
# It is the MODEL — a pixel-integrated image is not a polygon, and its edges are gone. The
# quadrature route has no such floor, because refining `ngauss`/`nsub` refines the integral
# itself rather than a picture of it: MEASURED at 6.8e-7 (HEALPix 3) and 2.5e-9 (4 and 5) with
# `ngauss = 4`, the latter being FINUFFT's tolerance rather than the quadrature's error.

"""
    quadrature_for(proj_west, proj_north, kx, ky; target_rad=1.0) -> (ngauss, nsub)

Pick a quadrature that resolves the OSCILLATION, not one that hopes to.

The integrand is `exp(-2πi k·r)` over a quad, so what decides the rule is the phase span
across the widest quad, `2π·|k|max·diag`. A fixed rule is therefore accurate at one mesh and
not at another — measured with `ngauss = 4, nsub = 1`, the error is 2.9e-4 at HEALPix 2 and
2.5e-9 at HEALPix 4, because coarsening the mesh makes the quads bigger and the phase across
one of them larger.

Subdivision rather than a higher order: the sub-square span falls as `1/nsub` and the error as
`1/nsub⁴`, which is the cheaper way to buy accuracy once the span is more than a couple of
radians (`ngauss` beyond ~6 buys little on an oscillatory integrand).

`target_rad = 2.5` is CALIBRATED, not chosen: on polaris the spans run 16.5 rad at HEALPix 1
down to 1.3 at HEALPix 5, and the `nsub` this returns for each is the smallest that reaches
the NUFFT's own tolerance — 2.1e-9 to 2.5e-9 at every level, where a fixed `nsub = 1` ranges
from 2.0e-2 to 2.5e-9. Uniform accuracy across the mesh is the property that made this
backend safe to default to.

NOTE THAT `:nufft` IS NO LONGER THE DEFAULT — `:t3` is, and it reads its rule from
[`quadrature_for_type3`](@ref) instead, which raises the Gauss ORDER where this one
subdivides. This rule is unchanged so that `:nufft` keeps behaving exactly as it did and
remains a usable independent cross-check.
"""
function quadrature_for(proj_west::AbstractMatrix, proj_north::AbstractMatrix,
                        kx::AbstractVector, ky::AbstractVector; target_rad::Real = 2.5)
    kmax = zero(float(eltype(kx)))
    @inbounds for i in eachindex(kx)
        k2 = kx[i]^2 + ky[i]^2
        k2 > kmax && (kmax = k2)
    end
    kmax = sqrt(kmax)
    dmax = zero(float(eltype(proj_west)))
    @inbounds for q in axes(proj_west, 1)
        d1 = hypot(proj_west[q,3] - proj_west[q,1], proj_north[q,3] - proj_north[q,1])
        d2 = hypot(proj_west[q,4] - proj_west[q,2], proj_north[q,4] - proj_north[q,2])
        dmax = max(dmax, d1, d2)
    end
    span = 2π * kmax * dmax                      # radians across the widest quad
    nsub = clamp(ceil(Int, span / target_rad), 1, 8)
    return 4, nsub
end

"""
    polyft_cvis_nufft(proj_west, proj_north, xw, kx, ky; ngauss=4, nsub=1, tol=1e-9)

Complex visibilities at SCATTERED uv points, by Gauss-Legendre quadrature over each quad
folded into one type-3 NUFFT.

Returns the same unnormalised sum as [`compute_polyflux_and_cvis!`](@ref) — the caller divides
by the flux — and agrees with it to `tol` for `ngauss ≥ 4`.

Its cost is set by the SOURCE COUNT (`nvis · (ngauss·nsub)²`) and the target count, not by the
product of the two, so it overtakes the direct kernel as the mesh is refined: MEASURED against
`:turbo` at 8.9× (HEALPix 3), 20.9× (4) and 77.6× (5) on polaris.

`ngauss` and `nsub` default to whatever [`quadrature_for`](@ref) reads off the geometry, which
is what makes the accuracy uniform across meshes: 2.1e-9 to 2.5e-9 from HEALPix 1 to 5 on
polaris, against 2.0e-2 to 2.5e-9 for a fixed 1-subdivision rule. Pass them explicitly only to
override that.

# Precision

`T` is the float type the transform RUNS in, and it defaults to `Float64` even when the mesh
is `Float32`. That pin is DELIBERATE and predates this note — it is here for accuracy — and
the numbers below are the measurement that backs it, recorded because the file carried the
choice without the reason and it read as an accidental promotion.

Single precision is available (`T = Float32`, or `T = nufft_work_type(proj_west, kx)` to
follow the inputs) and it is correct; it is simply not worth taking by default. MEASURED on one
lam And epoch, Float32 inputs throughout:

| | HEALPix 3 | HEALPix 5 |
|---|---|---|
| `T = Float64` | 1.720 ms, 1982 KiB | 2.279 ms, 3440 KiB |
| `T = Float32` | 1.755 ms, 1270 KiB | 2.150 ms, 1999 KiB |

So single buys 35-42 % of the memory and, at the mesh sizes in use, nothing at all in time —
because the type-3's cost here is split between `setpts`, which bins and sorts the points and
is index-bound rather than arithmetic-bound (~0.5 ms, unchanged by the word size), and an
`exec` whose FFT grid is set by the space-bandwidth product rather than by the precision.

What it costs is accuracy, against the exact closed-form kernel: 5.8e-10 in double against
1.9e-6 in single on lam And and 1.0e-5 on polaris. That is the Float32 type-3's own floor and
not slack — `nufft_tol(Float32)` is already at FINUFFT's single-precision limit, and asking
double for the same 1e-6 gives 7.2e-8. It also puts the TRANSFORM above the quadrature's
6.8e-7 as the limiting term, which inverts the error budget the adaptive rule was built for.

So the pin stands. `Float32` is the right choice only where memory is the binding constraint
rather than the last four digits.

What WAS broken, separately from the default, is that single precision could not be reached at
all: the targets were written `collect(Float64, 2π .* kx)`, and `2π` is a `Float64`, so a
`Float32` request promoted anyway. That is fixed — the constant is typed and the work goes
through a barrier — which is what makes `T = Float32` mean something when it is asked for.
"""
function polyft_cvis_nufft(proj_west::AbstractMatrix, proj_north::AbstractMatrix,
                           xw::AbstractVector, kx::AbstractVector, ky::AbstractVector;
                           ngauss::Union{Nothing,Int} = nothing,
                           nsub::Union{Nothing,Int} = nothing,
                           T::Type = Float64,
                           tol::Real = nufft_tol(T))
    # ADAPTIVE by default. A fixed rule is accurate at one mesh and not at another, and the
    # caller has no way to know which — `quadrature_for` reads it off the geometry.
    ag, as = quadrature_for(proj_west, proj_north, kx, ky)
    ng = ngauss === nothing ? ag : ngauss
    ns = nsub   === nothing ? as : nsub
    # THROUGH A BARRIER, because `T::Type = …` in a keyword list without a `where T` infers
    # `DataType` rather than the type itself, and everything built from it — `T(2π)`, the
    # sample arrays — comes back `Any`. Same fix as `finish_star`/`_finish_star`.
    return _polyft_cvis_nufft(T, proj_west, proj_north, xw, kx, ky, ng, ns, tol)
end

function _polyft_cvis_nufft(::Type{T}, proj_west, proj_north, xw, kx, ky,
                            ngauss::Int, nsub::Int, tol::Real) where {T}
    xs, ys, fs = build_gauss_samples(proj_west, proj_north, xw;
                                     ngauss = ngauss, nsub = nsub, T = T)
    # TYPED, and this is what makes the Float32 path possible at all. `2π .* kx` on a Float32 `kx`
    # returns Float64 — `2π` is a Float64 — so writing the targets the obvious way silently
    # promoted the transform and every later array with it. Built element by element in `T`.
    twopi = T(2π)
    sk = Vector{T}(undef, length(kx)); tk = Vector{T}(undef, length(ky))
    @inbounds for i in eachindex(kx)
        sk[i] = twopi * T(kx[i])
        tk[i] = twopi * T(ky[i])
    end
    finufft_available() ||
        error("the :nufft kernel needs FINUFFT loaded: add `using FINUFFT` to this session. " *
              "It is a weak dependency because loading it costs 464 ms on every `using ROTIR` " *
              "and drags 455 MB of CUDA driver into an application bundle, and `:t3` — the " *
              "default — is both faster and more accurate. `:nufft` is kept as the " *
              "independent cross-check, not as the working kernel.")
    # FINUFFT type 3 computes Σ_j c_j exp(iσ(s_k x_j + t_k y_j)); σ = -1 with the 2π folded
    # into the targets gives exp(-2πi k·r), which is the polygon FT's convention.
    #
    # `invokelatest` for the same reason `_cvis_forward!` needs it on `:turbo`: the GUI loads
    # FINUFFT ON DEMAND, in the middle of a session, from a callback the Qt event loop
    # dispatches — and methods added after the calling frame's world age are invisible to a
    # direct call. The cost is one dynamic dispatch per visibility computation, not per element.
    return Base.invokelatest(_finufft2d3, xs, ys, fs, -1, T(tol), sk, tk)
end

# Provided by ext/ROTIRFINUFFTExt.jl. Declared here so `:nufft` is a name this package knows
# and the failure without FINUFFT is a sentence rather than a MethodError on an underscored
# internal — the same arrangement as `_cvis_turbo!` in src/fused_polyft.jl.
function _finufft2d3 end

"""
    FINUFFT_OK

Set to `true` by `ROTIRFINUFFTExt.__init__` when that extension loads. See
[`finufft_available`](@ref).
"""
const FINUFFT_OK = Ref(false)

"""
    finufft_available() -> Bool

Whether the `:nufft` kernel will work here, i.e. whether `using FINUFFT` has loaded
ROTIRFINUFFTExt.

A CACHED FLAG set by the extension, not a `methods` lookup — same reasoning as
[`turbo_available`](@ref), and the same staleness trap avoided: the GUI can load FINUFFT
mid-session, so anything cached at ROTIR's own load time would be wrong afterwards.
"""
finufft_available() = FINUFFT_OK[]

"""
    nufft_work_type(proj_west, kx) -> Type

The float type the type-3 runs in: `Float32` when BOTH the mesh and the uv coordinates are
`Float32`, and `Float64` otherwise.

FINUFFT supports exactly these two (`FINUFFT.finufftReal`), so anything else — `Float16`,
`BigFloat` — resolves to `Float64` rather than failing inside the C call.
"""
function nufft_work_type(proj_west, kx)
    T = promote_type(float(real(eltype(proj_west))), float(real(eltype(kx))))
    return T === Float32 ? Float32 : Float64
end

"""
    nufft_tol(::Type{T}) -> Real

The tolerance to ask FINUFFT for when working in `T`.

`1e-9` in double and `1e-6` in single, and the single-precision figure is a FLOOR rather than a
preference: FINUFFT's single-precision kernels cannot deliver better than about `1e-6`, and
asking for less simply widens the spreading kernel for an accuracy the arithmetic cannot hold.

It is also not the binding term. The QUADRATURE's own error is 6.8e-7 at HEALPix 3 (see
[`quadrature_for`](@ref)), so at that mesh the two are the same size and the transform is not
what limits the answer.
"""
nufft_tol(::Type{Float32}) = 1.0e-6
nufft_tol(::Type{T}) where {T} = 1.0e-9

"""
    polyft_cvis_nufft_f64(proj_west, proj_north, xw, kx, ky; kwargs...)

[`polyft_cvis_nufft`](@ref) pinned to `Float64` — the reference the single-precision path is
measured against, in accuracy and in speed.

A PINNED ENTRY POINT rather than a frozen copy of the old body. The question it exists to
answer is "what does dropping to Float32 cost, for the same algorithm", and a copy would
answer it against whatever the algorithm used to be as soon as either changed. The exact
reference for ACCURACY is neither of these: it is the closed-form kernel (`:scalar`/`:turbo`),
which evaluates the polygon transform with no quadrature and no tolerance at all.
"""
polyft_cvis_nufft_f64(proj_west::AbstractMatrix, proj_north::AbstractMatrix,
                      xw::AbstractVector, kx::AbstractVector, ky::AbstractVector;
                      kwargs...) =
    polyft_cvis_nufft(proj_west, proj_north, xw, kx, ky; T = Float64, kwargs...)
