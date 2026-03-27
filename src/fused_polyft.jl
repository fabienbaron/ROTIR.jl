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
    F .= zero(Complex{T})

    @inbounds for p in 1:npix
        xw_p = xw[p]
        xw_p == zero(T) && continue

        # Shoelace area for this pixel
        pf = T(0.5) * (
            proj_west[p,1]*proj_north[p,2] - proj_west[p,2]*proj_north[p,1] +
            proj_west[p,2]*proj_north[p,3] - proj_west[p,3]*proj_north[p,2] +
            proj_west[p,3]*proj_north[p,4] - proj_west[p,4]*proj_north[p,3] +
            proj_west[p,4]*proj_north[p,1] - proj_west[p,1]*proj_north[p,4])
        polyflux[p] = pf

        # Accumulate FT contribution from 4 edges of the quad
        for e in 1:4
            j1 = e
            j2 = mod1(e+1, 4)
            dx = proj_west[p,j2] - proj_west[p,j1]
            dy = proj_north[p,j2] - proj_north[p,j1]
            cx = proj_west[p,j2] + proj_west[p,j1]
            cy = proj_north[p,j2] + proj_north[p,j1]

            for k in 1:nuv
                kdd = kx[k]*dx + ky[k]*dy        # k·d (dot product with edge diff)
                kdc = kx[k]*cx + ky[k]*cy        # k·c (dot product with edge center)
                cr_k = ky[k]*dx - kx[k]*dy       # perpendicular component
                s = sinc(kdd)                     # sinc(k·d)
                phase = cis(-T(π) * kdc)          # exp(-iπ k·c)
                F[k] += k2_inv_im[k] * s * phase * cr_k * xw_p
            end
        end
    end
    return nothing
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

    @inbounds for p in 1:npix
        acc = zero(Complex{T})

        for e in 1:4
            j1 = e
            j2 = mod1(e+1, 4)
            dx = proj_west[p,j2] - proj_west[p,j1]
            dy = proj_north[p,j2] - proj_north[p,j1]
            cx = proj_west[p,j2] + proj_west[p,j1]
            cy = proj_north[p,j2] + proj_north[p,j1]

            for k in 1:nuv
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

    @inbounds for p in 1:npix
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
    w = star.vis_weights[indx]
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
        println("V2: ", chi2_v2/data.nv2, "\tT3A: ", chi2_t3amp/data.nt3amp,
                "\tT3P: ", chi2_t3phi/data.nt3phi, "\tFlux: ", flux)
    end
    return chi2_v2 + chi2_t3amp + chi2_t3phi
end
