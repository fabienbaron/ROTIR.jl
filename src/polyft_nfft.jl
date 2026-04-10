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
    build_gauss_samples(proj_west, proj_north, x_weighted; ngauss=4, nsub=1, T=Float32)

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
                              ngauss::Int=4, nsub::Int=1, T::Type=Float32)
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
                        ngauss=4, nsub=1, T=Float32)

Compute `F[k] ~ sum_p x_weighted[p] * integral_polygon_p exp(-2*pi*i * k . r) dA`
via an adjoint NFFT at Gauss-Legendre samples inside each quadrilateral.

Returns a `(nx/2+1, nx)` complex array whose layout matches
`reshape(polyft*x_weighted, nx/2+1, nx)`, so that the caller can feed it
directly to `irfft(F, nx)` to obtain the real-space forward model.

# Arguments
- `proj_west, proj_north`: `(Npix, 4)` polygon vertex coordinates (mas)
- `x_weighted`: `(Npix,)` intensity weights
- `pixsize`: mas per pixel
- `nx`: grid size (pixels per side)
- `ngauss`: Gauss-Legendre order per axis (2..10)
- `nsub`: number of sub-squares per axis for subdivision
"""
function polyft_nfft_forward(proj_west, proj_north, x_weighted, pixsize::Real,
                              nx::Integer; ngauss::Int=4, nsub::Int=1,
                              T::Type=Float32)
    xs, ys, fs = build_gauss_samples(proj_west, proj_north, x_weighted;
                                      ngauss=ngauss, nsub=nsub, T=T)
    Ns = length(xs)
    L  = T(nx * pixsize)

    pos = Matrix{T}(undef, 2, Ns)
    @inbounds for s in 1:Ns
        pos[1, s] = ys[s] / L
        pos[2, s] = xs[s] / L
    end

    p_plan   = plan_nfft(pos, (nx, nx))
    fhat_nat = adjoint(p_plan) * fs

    # rfft-layout extraction
    nh = nx ÷ 2 + 1
    return fftshift(fhat_nat[nh:-1:1, :], 2)
end

"""
    polyft_nfft_image(proj_west, proj_north, x_weighted, pixsize, nx;
                      ngauss=4, nsub=1, T=Float32)

Convenience function: compute the NFFT-based Fourier coefficients and
immediately inverse-FFT them to produce a real-space `(nx, nx)` image.

# Arguments
Same as `polyft_nfft_forward`.
"""
function polyft_nfft_image(proj_west, proj_north, x_weighted, pixsize::Real,
                            nx::Integer; ngauss::Int=4, nsub::Int=1,
                            T::Type=Float32)
    F = polyft_nfft_forward(proj_west, proj_north, x_weighted, pixsize, nx;
                             ngauss=ngauss, nsub=nsub, T=T)
    return fftshift(irfft(F, nx))
end
