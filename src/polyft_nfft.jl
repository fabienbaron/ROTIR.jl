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
from 2.0e-2 to 2.5e-9. Uniform accuracy across the mesh is the property that makes the backend
safe to default to.
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
"""
function polyft_cvis_nufft(proj_west::AbstractMatrix, proj_north::AbstractMatrix,
                           xw::AbstractVector, kx::AbstractVector, ky::AbstractVector;
                           ngauss::Union{Nothing,Int} = nothing,
                           nsub::Union{Nothing,Int} = nothing, tol::Real = 1e-9)
    # ADAPTIVE by default. A fixed rule is accurate at one mesh and not at another, and the
    # caller has no way to know which — `quadrature_for` reads it off the geometry.
    ag, as = quadrature_for(proj_west, proj_north, kx, ky)
    ng = ngauss === nothing ? ag : ngauss
    ns = nsub   === nothing ? as : nsub
    xs, ys, fs = build_gauss_samples(proj_west, proj_north, xw;
                                     ngauss = ng, nsub = ns, T = Float64)
    # FINUFFT type 3 computes Σ_j c_j exp(iσ(s_k x_j + t_k y_j)); σ = -1 with the 2π folded
    # into the targets gives exp(-2πi k·r), which is the polygon FT's convention.
    return FINUFFT.nufft2d3(xs, ys, ComplexF64.(fs), -1, Float64(tol),
                            collect(Float64, 2π .* kx), collect(Float64, 2π .* ky))
end
