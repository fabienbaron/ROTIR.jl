# fused_polyft_jvp.jl — FORWARD mode for the exact polygon Fourier transform.
#
# Every leaf derivative in ROTIR is hand-coded REVERSE mode (`src/parametric_gradient.jl:8`),
# which is what a gradient needs. Variational inference needs the other direction too: MGVI and
# GeoVI (VarInf.jl) require both `left_sqrt_metric` (a VJP, J'w) and `right_sqrt_metric` (a JVP,
# J v). The VJP is already available — `compute_adjoint_cvis!` and `compute_adjoint_vertices!`
# take an arbitrary complex adjoint with nothing about the chi2 residual baked in. The JVP had
# nothing to reuse: grepping `frule|jvp|pushforward|ForwardDiff` over src/ ext/ test/ finds only
# `kepler_E`'s rule. This file supplies it.
#
# THE TRANSPOSE, LITERALLY. `compute_adjoint_vertices!` runs
#     for p (facets, threaded) / for e (edges) / for k (uv)  ->  accumulate into grad_proj[p,:]
# contracting E's derivatives against `adj[k]` and reducing over k. The JVP is the same nest with
# the reduction swapped: contract against the vertex and weight TANGENTS and reduce over p into
# F_dot[k]. Threading therefore moves to the OUTER k loop, because every p contributes to every
# F_dot[k] and threading over p would race; threading over k has each task own one output.
#
# The algebra is copied term for term from the adjoint so the two cannot drift:
#     E(p,e,k) = k2_inv_im[k] * sinc(kdd) * cis(-pi*kdc) * cr_k        (times xw[p])
#     dE/d(dx) = base * (ds*kx*cr + s*ky)      dE/d(cx) = base * s * cr * (-i*pi*kx)
#     dE/d(dy) = base * (ds*ky*cr - s*kx)      dE/d(cy) = base * s * cr * (-i*pi*ky)
# with base = k2_inv_im * phase * xw[p], ds = (cos(pi*kdd) - sinc(kdd))/kdd, and the edge
# parameterisation dx = x[j2]-x[j1], cx = x[j2]+x[j1] giving the tangents
#     ddx = xdot[j2]-xdot[j1],  dcx = xdot[j2]+xdot[j1].
#
# VERIFIED three ways in `test/test_polyft_jvp.jl`: against central finite differences of the
# forward transform; against the ADJOINT IDENTITY <J v, w> = <v, J' w> using the existing,
# independently FD-checked VJP (the strongest check, since it ties the two directions together);
# and for linearity in the tangent.

"""
    compute_jvp_cvis!(F_dot, xw_dot, proj_west_dot, proj_north_dot,
                      kx, ky, k2_inv_im, proj_west, proj_north, xw)

Forward-mode (JVP) companion to [`compute_adjoint_cvis!`](@ref) and
[`compute_adjoint_vertices!`](@ref): pushes the tangents `(xw_dot, proj_west_dot,
proj_north_dot)` through the exact polygon Fourier transform, writing the directional derivative
of the unnormalised complex visibility into `F_dot`.

`F_dot[k] = Σ_p [ xw_dot[p]·E(p,k) + xw[p]·(∂E/∂proj)·proj_dot ]`, where `E` is the same
closed-form polygon kernel the forward pass evaluates. Pass `xw_dot` as zeros to differentiate
with respect to geometry alone, or `proj_*_dot` as zeros for the map alone.

This is what VarInf.jl's `right_sqrt_metric` needs. It carries no `:turbo` branch yet; the
adjoints' `_adj_*_turbo!` kernels are the template if this ever becomes hot.
"""
@views function compute_jvp_cvis!(F_dot::Vector{Complex{T}},
    xw_dot::Vector{T}, proj_west_dot::AbstractMatrix{T}, proj_north_dot::AbstractMatrix{T},
    kx::Vector{T}, ky::Vector{T}, k2_inv_im::Vector{Complex{T}},
    proj_west::AbstractMatrix{T}, proj_north::AbstractMatrix{T}, xw::Vector{T}) where T

    nuv  = length(kx)
    npix = size(proj_west, 1)
    @assert length(F_dot) == nuv
    @assert length(xw_dot) == npix == length(xw)

    Threads.@threads for k in 1:nuv        # each task writes F_dot[k] only — thread-safe
      @inbounds begin
        kxk = kx[k]; kyk = ky[k]; k2i = k2_inv_im[k]
        acc = zero(Complex{T})
        for p in 1:npix
            xw_p = xw[p]; xwd_p = xw_dot[p]
            (xw_p == zero(T) && xwd_p == zero(T)) && continue
            for e in 1:4
                j1 = e
                j2 = mod1(e + 1, 4)
                dx = proj_west[p,j2]  - proj_west[p,j1]
                dy = proj_north[p,j2] - proj_north[p,j1]
                cx = proj_west[p,j2]  + proj_west[p,j1]
                cy = proj_north[p,j2] + proj_north[p,j1]

                kdd  = kxk*dx + kyk*dy
                kdc  = kxk*cx + kyk*cy
                cr_k = kyk*dx - kxk*dy
                s     = sinc(kdd)
                phase = cis(-T(π) * kdc)

                # value: the same kernel the forward pass sums, WITHOUT xw
                E = k2i * s * phase * cr_k
                acc += xwd_p * E                       # tangent in the surface weights

                xw_p == zero(T) && continue
                pikdd = T(π) * kdd
                ds = abs(kdd) < T(1e-12) ? zero(T) : (cos(pikdd) - s) / kdd
                base = k2i * phase * xw_p

                dE_ddx = base * (ds * kxk * cr_k + s * kyk)
                dE_ddy = base * (ds * kyk * cr_k - s * kxk)
                dE_dcx = base * s * cr_k * (-im * T(π) * kxk)
                dE_dcy = base * s * cr_k * (-im * T(π) * kyk)

                # dx = x[j2]-x[j1] and cx = x[j2]+x[j1], so the tangents combine as:
                ddx = proj_west_dot[p,j2]  - proj_west_dot[p,j1]
                dcx = proj_west_dot[p,j2]  + proj_west_dot[p,j1]
                ddy = proj_north_dot[p,j2] - proj_north_dot[p,j1]
                dcy = proj_north_dot[p,j2] + proj_north_dot[p,j1]

                acc += dE_ddx*ddx + dE_ddy*ddy + dE_dcx*dcx + dE_dcy*dcy
            end
        end
        F_dot[k] = acc
      end
    end
    return F_dot
end
