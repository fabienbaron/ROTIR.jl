# Fused two-pass polygon Fourier transform for interferometric observables.
# Eliminates the dense (Nk × Npix) polyft matrix by computing visibilities on-the-fly.
#
# Pass 1 (forward):  Accumulates complex visibilities F[k] = Σ_p polyft[k,p] * xw[p]
# Pass 2 (adjoint):  Accumulates ∂χ²/∂xw[p] and optionally ∂χ²/∂proj_west, ∂χ²/∂proj_north
#
# Adapted from planet_deconv's shape_gradient.jl for interferometric (sparse UV) use.

"""
    compute_polyflux_and_cvis!(F, polyflux, kx, ky, k2_inv_im, proj_west, proj_north, xw)

Forward pass: compute complex visibilities F[k] and pixel areas polyflux[p].
- `F`: output, length nuv — complex visibilities at each UV point
- `polyflux`: output, length npix — projected pixel areas (shoelace formula)
- `kx, ky`: UV frequencies (length nuv), pre-scaled to radians
- `k2_inv_im`: -im/(2π(kx²+ky²)) for each UV point (length nuv)
- `proj_west, proj_north`: projected quad vertices (npix × 4)
- `xw`: weighted pixel values (npix) = x .* vis_weights
"""
@views function compute_polyflux_and_cvis!(F::Vector{Complex{T}}, polyflux::Vector{T},
    kx::Vector{T}, ky::Vector{T}, k2_inv_im::Vector{Complex{T}},
    proj_west::AbstractMatrix{T}, proj_north::AbstractMatrix{T}, xw::Vector{T}) where T

    nuv = length(kx)
    npix = size(proj_west, 1)

    # Shoelace area per pixel (independent of xw — set for EVERY pixel so a zero-weight
    # pixel never leaves polyflux[p] uninitialized).
    @inbounds for p in 1:npix
        polyflux[p] = T(0.5) * (
            proj_west[p,1]*proj_north[p,2] - proj_west[p,2]*proj_north[p,1] +
            proj_west[p,2]*proj_north[p,3] - proj_west[p,3]*proj_north[p,2] +
            proj_west[p,3]*proj_north[p,4] - proj_west[p,4]*proj_north[p,3] +
            proj_west[p,4]*proj_north[p,1] - proj_west[p,1]*proj_north[p,4])
    end

    _cvis_forward!(F, kx, ky, k2_inv_im, proj_west, proj_north, xw)
    return nothing
end

"""
    fused_cvis(x, star, data; intensity_model=:linear, band=nothing) -> Vector{Complex}

Normalised complex visibilities WITHOUT building `polyft`.

The matrix-free counterpart of [`poly_to_cvis`](@ref), and identical to it in what it returns:
the same weighting (soft visibility × limb darkening), the same flux normalisation. The
difference is only that `poly_to_cvis` reads `star.polyft`, which `setup_oi!` must have built,
and this computes the sum directly.

Use it wherever the geometry changes per evaluation — a parametric fit — and `poly_to_cvis`
wherever it is fixed and the map varies, which is imaging.
"""
function fused_cvis(x, star, data; intensity_model::Symbol = :linear, band = nothing)
    T = eltype(star.proj_west)
    indx = star.index_quads_visible
    I = intensity_model === :linear ? x : intensity(x, intensity_model, band)
    xw = T.(I[indx] .* star.vis_weights[indx] .* star.ldmap[indx])
    pjx = star.proj_west[indx, :]; pjy = star.proj_north[indx, :]
    kx = T.(data.uv[1, :]) * T(-π / (180 * 3600000))
    ky = T.(data.uv[2, :]) * T( π / (180 * 3600000))
    k2 = precompute_k2_inv_im(kx, ky)
    F = Vector{Complex{T}}(undef, length(kx))
    pf = zeros(T, length(indx))
    mpjx = Matrix(pjx); mpjy = Matrix(pjy)
    if POLYFT_BACKEND[] === :nufft
        # `polyflux` is the shoelace area and is needed for the flux normalisation whichever
        # backend runs, so it is computed here rather than inside the visibility kernel.
        @inbounds for q in 1:length(indx)
            pf[q] = T(0.5) * (mpjx[q,1]*mpjy[q,2] - mpjx[q,2]*mpjy[q,1] +
                              mpjx[q,2]*mpjy[q,3] - mpjx[q,3]*mpjy[q,2] +
                              mpjx[q,3]*mpjy[q,4] - mpjx[q,4]*mpjy[q,3] +
                              mpjx[q,4]*mpjy[q,1] - mpjx[q,1]*mpjy[q,4])
        end
        Fn = polyft_cvis_nufft(mpjx, mpjy, xw, kx, ky)
        return Complex{T}.(Fn) ./ dot(pf, xw)
    end
    compute_polyflux_and_cvis!(F, pf, kx, ky, k2, mpjx, mpjy, xw)
    return F ./ dot(pf, xw)
end

"""
    POLYFT_BACKEND

Which forward kernel computes the visibilities: `:nufft` (default), `:turbo`, or `:scalar`.

All three compute the SAME quantity and are asserted against each other in
`test/test_fused_polyft.jl`. They differ in how, and therefore in how the cost scales:

| backend   | method                                   | HEALPix 3 | 4 | 5 |
|-----------|------------------------------------------|-----------|---|---|
| `:scalar` | the reference, plain Julia               | 119 ms | 425 ms | 1472 ms |
| `:turbo`  | the same sum, SIMD transcendentals       | 6.9 ms | 21.8 ms | 86 ms |
| `:nufft`  | Gauss quadrature + type-3 NUFFT          | 2.6 ms | 3.8 ms | 4.1 ms |

`:scalar` and `:turbo` are EXACT — the closed-form polygon transform, no parameters.
`:nufft` is a quadrature, accurate to 6.8e-7 at HEALPix 3 and 2.5e-9 above it, and its cost
barely moves with the mesh because it is set by the source and target counts rather than their
product. It is the default because those numbers say so, and because a fit is thousands of
these.

`:turbo` is what to fall back on if a mesh is coarse enough to worry about the quadrature (see
`polyft_cvis_nufft`: a 2-point rule is only good to 1.6e-3 at HEALPix 3), and `:scalar` is the
definition the other two are tested against — worth being able to select without rebuilding
when a fit looks wrong.
"""
const POLYFT_BACKEND = Ref(:nufft)

# `_cvis_turbo!` is provided by ext/ROTIRLoopVectorizationExt.jl. Declared here so that
# `:turbo` is a name this package knows and the failure without LoopVectorization is a
# sentence rather than a MethodError on an underscored internal.
function _cvis_turbo! end

"""
    turbo_available() -> Bool

Whether `POLYFT_BACKEND[] = :turbo` will work here, i.e. whether `using LoopVectorization` has
loaded ROTIRLoopVectorizationExt.
"""
turbo_available() = !isempty(methods(_cvis_turbo!))

function _cvis_forward!(F, kx, ky, k2, pjx, pjy, xw)
    if POLYFT_BACKEND[] === :turbo
        turbo_available() ||
            error("POLYFT_BACKEND[] = :turbo needs LoopVectorization loaded: add " *
                  "`using LoopVectorization` to this session. It is a weak dependency " *
                  "because loading it costs 1.8 s of GUI startup by invalidating OITOOLS' " *
                  "precompiled plot pipeline, and `:nufft` — the default — is faster anyway.")
        return _cvis_turbo!(F, kx, ky, k2, pjx, pjy, xw)
    end
    return _cvis_scalar!(F, kx, ky, k2, pjx, pjy, xw)
end

# `_cvis_turbo!` lives in src/turbo_polyft.jl, provided by ext/ROTIRLoopVectorizationExt.jl.
# See that file for why it is not here: loading LoopVectorization costs 1.8 s of GUI startup,
# and since `:nufft` became the default it buys a cross-check rather than the speed.

"""
    _cvis_scalar!(F, kx, ky, k2_inv_im, proj_west, proj_north, xw)

THE REFERENCE. Plain Julia, no vectorisation package, threaded over the UV index.

Kept — and tested against — because the `@turbo` kernel below rewrites the same expression
into a form LoopVectorization can vectorise (real accumulators, `ifelse` instead of a branch,
`sin`/`cos` instead of `sinc`/`cis`), and a rewrite of a numerically delicate expression needs
something to be a rewrite OF. `test/test_fused_polyft.jl` asserts the two agree across
precisions, mesh levels and surface types.

Threading over k rather than p is what makes this one race-free without accumulators: each k
owns its own `F[k]`.
"""
@views function _cvis_scalar!(F::Vector{Complex{T}}, kx::Vector{T}, ky::Vector{T},
    k2_inv_im::Vector{Complex{T}}, proj_west::AbstractMatrix{T},
    proj_north::AbstractMatrix{T}, xw::Vector{T}) where T
    nuv = length(kx); npix = size(proj_west, 1)
    Threads.@threads for k in 1:nuv
        kxk = kx[k]; kyk = ky[k]; k2k = k2_inv_im[k]
        acc = zero(Complex{T})
        @inbounds for p in 1:npix
            xw_p = xw[p]
            xw_p == zero(T) && continue
            for e in 1:4
                j1 = e; j2 = mod1(e+1, 4)
                dx = proj_west[p,j2] - proj_west[p,j1]
                dy = proj_north[p,j2] - proj_north[p,j1]
                cx = proj_west[p,j2] + proj_west[p,j1]
                cy = proj_north[p,j2] + proj_north[p,j1]
                kdd = kxk*dx + kyk*dy
                kdc = kxk*cx + kyk*cy
                cr  = kyk*dx - kxk*dy
                acc += k2k * sinc(kdd) * cis(-T(π)*kdc) * cr * xw_p
            end
        end
        F[k] = acc
    end
    return F
end

"""
    compute_adjoint_cvis!(grad_xw, adj, kx, ky, k2_inv_im, proj_west, proj_north, polyflux)

Adjoint pass: compute gradient of chi2 w.r.t. weighted pixel values.
- `grad_xw`: output, length npix — ∂χ²/∂(xw[p])
- `adj`: input, length nuv — adjoint signal in complex visibility space
- Other arguments same as forward pass.
"""
@views function compute_adjoint_cvis!(grad_xw::Vector{T},
    adj::Vector{Complex{T}}, kx::Vector{T}, ky::Vector{T},
    k2_inv_im::Vector{Complex{T}},
    proj_west::AbstractMatrix{T}, proj_north::AbstractMatrix{T}, polyflux::Vector{T}) where T

    nuv = length(kx)
    npix = size(proj_west, 1)
    grad_xw .= zero(T)

    Threads.@threads for p in 1:npix          # each p writes grad_xw[p] only — thread-safe
        acc = zero(Complex{T})

        for e in 1:4
            j1 = e
            j2 = mod1(e+1, 4)
            dx = proj_west[p,j2] - proj_west[p,j1]
            dy = proj_north[p,j2] - proj_north[p,j1]
            cx = proj_west[p,j2] + proj_west[p,j1]
            cy = proj_north[p,j2] + proj_north[p,j1]

            @inbounds for k in 1:nuv
                kdd = kx[k]*dx + ky[k]*dy
                kdc = kx[k]*cx + ky[k]*cy
                cr_k = ky[k]*dx - kx[k]*dy
                s = sinc(kdd)
                phase = cis(-T(π) * kdc)
                acc += adj[k] * k2_inv_im[k] * s * phase * cr_k
            end
        end

        grad_xw[p] = real(acc)
    end
    return nothing
end

"""
    compute_adjoint_vertices!(grad_proj_west, grad_proj_north, adj, kx, ky, k2_inv_im,
                               proj_west, proj_north, xw)

Adjoint pass for vertex positions: compute ∂χ²/∂proj_west and ∂χ²/∂proj_north.
Used for shape gradient computation.
- `grad_proj_west, grad_proj_north`: output (npix × 4) — vertex position gradients
"""
@views function compute_adjoint_vertices!(grad_proj_west::Matrix{T}, grad_proj_north::Matrix{T},
    adj::Vector{Complex{T}}, kx::Vector{T}, ky::Vector{T},
    k2_inv_im::Vector{Complex{T}},
    proj_west::AbstractMatrix{T}, proj_north::AbstractMatrix{T}, xw::Vector{T}, polyflux::Vector{T}) where T

    nuv = length(kx)
    npix = size(proj_west, 1)
    grad_proj_west .= zero(T)
    grad_proj_north .= zero(T)

    Threads.@threads for p in 1:npix          # each p writes its own grad_proj rows — thread-safe
      @inbounds begin
        xw_p = xw[p]
        xw_p == zero(T) && continue

        # DC (flux) gradient — shoelace formula derivatives
        # polyflux = 0.5 * Σ (x[j]*y[j+1] - x[j+1]*y[j])
        # ∂polyflux/∂x[j] = 0.5*(y[j+1] - y[j-1])
        # ∂polyflux/∂y[j] = 0.5*(x[j-1] - x[j+1])
        for j in 1:4
            jp = mod1(j+1, 4)
            jm = mod1(j-1, 4)
            # Note: DC contribution would be from adj[DC] * xw_p * ∂polyflux/∂vertex
            # but for interferometry the DC term is not in the UV data (no zero-spacing)
        end

        # Non-DC gradient: derivative of FT w.r.t. vertex positions
        for e in 1:4
            j1 = e
            j2 = mod1(e+1, 4)
            dx = proj_west[p,j2] - proj_west[p,j1]
            dy = proj_north[p,j2] - proj_north[p,j1]
            cx = proj_west[p,j2] + proj_west[p,j1]
            cy = proj_north[p,j2] + proj_north[p,j1]

            acc_j1x = zero(T)
            acc_j1y = zero(T)
            acc_j2x = zero(T)
            acc_j2y = zero(T)

            for k in 1:nuv
                kdd = kx[k]*dx + ky[k]*dy
                kdc = kx[k]*cx + ky[k]*cy
                cr_k = ky[k]*dx - kx[k]*dy
                s = sinc(kdd)
                phase = cis(-T(π) * kdc)

                # Derivative of sinc(kdd): dsinc = (cos(π*kdd) - sinc(kdd))/kdd
                pikdd = T(π) * kdd
                if abs(kdd) < T(1e-12)
                    ds = zero(T)  # dsinc(0) = 0
                else
                    ds = (cos(pikdd) - s) / kdd
                end

                # E = k2_inv_im * sinc(kdd) * cis(-π*kdc) * cr
                # ∂E/∂(dx) = k2_inv_im * [ds*kx * phase * cr + s * phase * ky] * xw_p
                # ∂E/∂(cx) = k2_inv_im * s * (-iπ*kx) * phase * cr * xw_p
                # And dx = x[j2]-x[j1], cx = x[j2]+x[j1], so:
                # ∂E/∂x[j1] = -∂E/∂(dx) + ∂E/∂(cx), ∂E/∂x[j2] = ∂E/∂(dx) + ∂E/∂(cx)

                base = k2_inv_im[k] * phase * xw_p
                ac = adj[k]

                # Contribution from sinc derivative (via kdd)
                dE_ddx = base * (ds * kx[k] * cr_k + s * ky[k])
                dE_ddy = base * (ds * ky[k] * cr_k - s * kx[k])

                # Contribution from phase derivative (via kdc)
                dE_dcx = base * s * cr_k * (-im * T(π) * kx[k])
                dE_dcy = base * s * cr_k * (-im * T(π) * ky[k])

                # Chain rule: dx = x2-x1, cx = x2+x1
                contrib_j1x = real(ac * (-dE_ddx + dE_dcx))
                contrib_j2x = real(ac * ( dE_ddx + dE_dcx))
                contrib_j1y = real(ac * (-dE_ddy + dE_dcy))
                contrib_j2y = real(ac * ( dE_ddy + dE_dcy))

                acc_j1x += contrib_j1x
                acc_j1y += contrib_j1y
                acc_j2x += contrib_j2x
                acc_j2y += contrib_j2y
            end

            grad_proj_west[p, j1] += acc_j1x
            grad_proj_west[p, j2] += acc_j2x
            grad_proj_north[p, j1] += acc_j1y
            grad_proj_north[p, j2] += acc_j2y
        end
      end   # close @inbounds begin
    end
    return nothing
end

"""
    precompute_k2_inv_im(kx, ky)
Precompute -im / (2π * (kx² + ky²)) for each UV point.
"""
function precompute_k2_inv_im(kx::Vector{T}, ky::Vector{T}) where T
    nuv = length(kx)
    k2_inv_im = Vector{Complex{T}}(undef, nuv)
    @inbounds for k in 1:nuv
        k2 = kx[k]^2 + ky[k]^2
        if k2 > T(1e-30)
            k2_inv_im[k] = Complex{T}(zero(T), -T(1/(2π))) / k2
        else
            k2_inv_im[k] = zero(Complex{T})
        end
    end
    return k2_inv_im
end

"""
    fused_spheroid_chi2_fg(x, g, star, data; verbose=true)

Compute chi2 and gradient using the fused two-pass approach (no polyft matrix).
Drop-in replacement for spheroid_chi2_fg.
"""
@views function fused_spheroid_chi2_fg(x, g, star, data; verbose::Bool = true)
    npix = star.npix
    T = eltype(x)
    indx = star.index_quads_visible
    w = star.vis_weights[indx] .* star.ldmap[indx]  # soft visibility × limb darkening
    xw = x[indx] .* w
    pjx = star.proj_west[indx, :]
    pjy = star.proj_north[indx, :]
    nvis = length(indx)

    # Precompute UV frequencies
    kx = data.uv[1,:] * T(-π / (180*3600000))
    ky = data.uv[2,:] * T( π / (180*3600000))
    nuv = length(kx)
    k2_inv_im = precompute_k2_inv_im(kx, ky)

    # Forward pass: compute complex visibilities
    F = Vector{Complex{T}}(undef, nuv)
    polyflux_local = Vector{T}(undef, nvis)
    compute_polyflux_and_cvis!(F, polyflux_local, kx, ky, k2_inv_im, pjx, pjy, xw)

    flux = dot(polyflux_local, xw)
    cvis_model = F / flux

    # Compute observables
    v2_model = abs2.(cvis_model[data.indx_v2])
    t3_model = cvis_model[data.indx_t3_1] .* cvis_model[data.indx_t3_2] .* cvis_model[data.indx_t3_3]
    t3amp_model = abs.(t3_model)
    t3phi_model = angle.(t3_model) * T(180/π)

    # Chi2
    chi2_v2 = sum(abs2, (v2_model - data.v2) ./ data.v2_err)
    chi2_t3amp = sum(abs2, (t3amp_model - data.t3amp) ./ data.t3amp_err)
    chi2_t3phi = sum(abs2, mod360(t3phi_model - data.t3phi) ./ data.t3phi_err)

    # Adjoint signal in complex visibility space: ∂χ²/∂cvis
    adj_cvis = zeros(Complex{T}, nuv)

    # V2 contribution: ∂χ²/∂cvis[k] = 4 * (v2_model - v2_data)/σ² * conj(cvis[k])
    for i in eachindex(data.indx_v2)
        k = data.indx_v2[i]
        adj_cvis[k] += 4 * (v2_model[i] - data.v2[i]) / data.v2_err[i]^2 * conj(cvis_model[k])
    end

    # T3amp contribution
    t3amp_res = 2 * (t3amp_model - data.t3amp) ./ data.t3amp_err.^2
    for i in eachindex(data.indx_t3_1)
        k1 = data.indx_t3_1[i]; k2 = data.indx_t3_2[i]; k3 = data.indx_t3_3[i]
        c1 = cvis_model[k1]; c2 = cvis_model[k2]; c3 = cvis_model[k3]
        a1 = abs(c1); a2 = abs(c2); a3 = abs(c3)
        adj_cvis[k1] += t3amp_res[i] * conj(c1)/a1 * a2 * a3
        adj_cvis[k2] += t3amp_res[i] * conj(c2)/a2 * a1 * a3
        adj_cvis[k3] += t3amp_res[i] * conj(c3)/a3 * a1 * a2
    end

    # T3phi contribution
    t3phi_res = mod360(t3phi_model - data.t3phi) ./ data.t3phi_err.^2
    for i in eachindex(data.indx_t3_1)
        k1 = data.indx_t3_1[i]; k2 = data.indx_t3_2[i]; k3 = data.indx_t3_3[i]
        c1 = cvis_model[k1]; c2 = cvis_model[k2]; c3 = cvis_model[k3]
        t3i = t3_model[i]
        factor = t3phi_res[i] / abs2(t3i) * conj(t3i)
        adj_cvis[k1] -= T(360/π) * im * factor * c2 * c3
        adj_cvis[k2] -= T(360/π) * im * factor * c1 * c3
        adj_cvis[k3] -= T(360/π) * im * factor * c1 * c2
    end

    # Scale adjoint for flux normalization: adj_F = adj_cvis / flux
    adj_F = adj_cvis / flux

    # Adjoint pass: compute gradient w.r.t. xw
    grad_xw = Vector{T}(undef, nvis)
    compute_adjoint_cvis!(grad_xw, adj_F, kx, ky, k2_inv_im, pjx, pjy, polyflux_local)

    # Flux correction: ∂χ²/∂xw also has a term from ∂flux/∂xw
    flux_adj = -dot(xw, grad_xw) / flux
    grad_xw .+= flux_adj * polyflux_local

    # Chain rule through soft visibility: ∂χ²/∂x = w .* ∂χ²/∂xw
    g[indx] = w .* grad_xw

    if verbose
        printstyled(@sprintf("V2: %.4f ", chi2_v2/data.nv2), color=:red)
        printstyled(@sprintf("T3A: %.4f ", chi2_t3amp/data.nt3amp), color=:blue)
        printstyled(@sprintf("T3P: %.4f ", chi2_t3phi/data.nt3phi), color=:green)
        printstyled(@sprintf("Flux: %.4f\n", flux), color=:normal)
    end
    return chi2_v2 + chi2_t3amp + chi2_t3phi
end
