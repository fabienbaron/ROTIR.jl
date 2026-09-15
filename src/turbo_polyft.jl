# The `@turbo` forward kernel, in an extension because of what LoopVectorization costs to LOAD.
#
# It is 17x the scalar reference and was the default until the type-3 NUFFT backend arrived.
# That changed the trade: `:nufft` is another 1.5-14x faster again, so `:turbo` is now the
# EXACT cross-check rather than the workhorse — and MEASURED, `using LoopVectorization`
# invalidates the plot-construction code OITOOLS precompiles for its live canvas, taking one
# `build_canvas` from 341 ms to 2685 ms and the GUI's time-to-window from 7.2 s to 9.0 s.
#
# Paying 1.8 s of every session for a cross-check is the wrong way round, and it is not
# recoverable by precompiling: `_oitools_gui()` is `nothing` throughout ROTIRGUIExt's
# workload, so OITOOLS' half of that canvas can only ever be precompiled by OITOOLS. The only
# lever is not to invalidate it, which means not loading LoopVectorization unless asked.
#
# `using LoopVectorization` brings this back and `POLYFT_BACKEND[] = :turbo` then works.

"""
    _cvis_turbo!(F, kx, ky, k2_inv_im, proj_west, proj_north, xw)

The same sum as [`_cvis_scalar!`](@ref), vectorised over the UV index.

THREE REWRITES are what LoopVectorization needs, and each changes nothing mathematically:

  * **Real accumulators.** `@turbo` does not vectorise complex arithmetic, so the real and
    imaginary parts are accumulated separately and assembled at the end.
    `k2_inv_im[k] = -i·c` is purely imaginary by construction, so only `c = -imag(k2_inv_im)`
    is carried — which also inherits the reference's zeroing of degenerate baselines
    (`precompute_k2_inv_im` returns 0 there rather than `Inf`).
  * **`sin(πa)/(πa)` with an `ifelse`** rather than `sinc`, because a branch does not
    vectorise and `ifelse` compiles to a blend.
  * **Loop order.** The reference threads over k and walks p inside; this walks p outside and
    vectorises over k, which is the contiguous axis. That inverts the threading — several
    threads now touch every `F[k]` — so each chunk of tessels gets its own accumulator and
    they are summed at the end. The chunks are formed explicitly: `Threads.threadid()` is NOT
    an index into `1:nthreads()` on Julia 1.12, which has an interactive pool as well, and
    indexing a per-thread array by it is a `BoundsError` waiting for the scheduler.
"""
function _cvis_turbo!(F::Vector{Complex{T}}, kx::Vector{T}, ky::Vector{T},
    k2_inv_im::Vector{Complex{T}}, proj_west::AbstractMatrix{T},
    proj_north::AbstractMatrix{T}, xw::Vector{T}) where T

    nuv = length(kx); npix = size(proj_west, 1)
    c = Vector{T}(undef, nuv)
    @inbounds @simd for k in 1:nuv
        c[k] = -imag(k2_inv_im[k])
    end
    # ONE chunk per ~64 tessels, capped at the thread count. The per-chunk accumulator is
    # `nuv` long whatever the tessel count is, so splitting 24 tessels eight ways costs eight
    # full-length buffers to do almost nothing: MEASURED at 41 ms of overhead against 0.3 ms
    # of work, which made the vectorised kernel SLOWER than the reference below level 3.
    nt = clamp(cld(npix, 64), 1, Threads.nthreads())
    # One allocation, and each chunk gets a CONTIGUOUS (2, nuv) slice of it — eight separate
    # arrays made the merge jump between eight 469 KiB buffers per k.
    partial = zeros(T, 2, nuv, nt)
    bounds = [round(Int, npix * (i - 1) / nt) + 1 : round(Int, npix * i / nt) for i in 1:nt]
    Threads.@sync for i in 1:nt
        Threads.@spawn _cvis_turbo_chunk!(view(partial, :, :, i), kx, ky, c,
                                          proj_west, proj_north, xw, bounds[i])
    end
    # Merged chunk-outer / k-inner, so each pass is a contiguous sweep.
    @inbounds for i in 2:nt
        @turbo for k in 1:nuv
            partial[1,k,1] += partial[1,k,i]
            partial[2,k,1] += partial[2,k,i]
        end
    end
    @inbounds @simd for k in 1:nuv
        F[k] = Complex{T}(partial[1,k,1], partial[2,k,1])
    end
    return F
end

function _cvis_turbo_chunk!(acc::AbstractMatrix{T}, kx::Vector{T}, ky::Vector{T}, c::Vector{T},
    proj_west::AbstractMatrix{T}, proj_north::AbstractMatrix{T},
    xw::Vector{T}, ps) where T

    nuv = length(kx)
    @inbounds for p in ps
        xwp = xw[p]
        xwp == zero(T) && continue          # the reference skips these; so does this
        x1 = proj_west[p,1];  x2 = proj_west[p,2];  x3 = proj_west[p,3];  x4 = proj_west[p,4]
        y1 = proj_north[p,1]; y2 = proj_north[p,2]; y3 = proj_north[p,3]; y4 = proj_north[p,4]
        d1x = x2-x1; d1y = y2-y1; c1x = x2+x1; c1y = y2+y1
        d2x = x3-x2; d2y = y3-y2; c2x = x3+x2; c2y = y3+y2
        d3x = x4-x3; d3y = y4-y3; c3x = x4+x3; c3y = y4+y3
        d4x = x1-x4; d4y = y1-y4; c4x = x1+x4; c4y = y1+y4
        @turbo for k in 1:nuv
            kxk = kx[k]; kyk = ky[k]
            a1 = kxk*d1x + kyk*d1y; b1 = -T(π)*(kxk*c1x + kyk*c1y)
            a2 = kxk*d2x + kyk*d2y; b2 = -T(π)*(kxk*c2x + kyk*c2y)
            a3 = kxk*d3x + kyk*d3y; b3 = -T(π)*(kxk*c3x + kyk*c3y)
            a4 = kxk*d4x + kyk*d4y; b4 = -T(π)*(kxk*c4x + kyk*c4y)
            p1 = T(π)*a1; p2 = T(π)*a2; p3 = T(π)*a3; p4 = T(π)*a4
            s1 = ifelse(a1 == zero(T), one(T), sin(p1)/p1)
            s2 = ifelse(a2 == zero(T), one(T), sin(p2)/p2)
            s3 = ifelse(a3 == zero(T), one(T), sin(p3)/p3)
            s4 = ifelse(a4 == zero(T), one(T), sin(p4)/p4)
            w1 = s1*(kyk*d1x - kxk*d1y); w2 = s2*(kyk*d2x - kxk*d2y)
            w3 = s3*(kyk*d3x - kxk*d3y); w4 = s4*(kyk*d4x - kxk*d4y)
            re = w1*cos(b1) + w2*cos(b2) + w3*cos(b3) + w4*cos(b4)
            im = w1*sin(b1) + w2*sin(b2) + w3*sin(b3) + w4*sin(b4)
            f = c[k] * xwp
            # −i·c·(re + i·im) = c·im − i·c·re
            acc[1,k] +=  f*im
            acc[2,k] += -f*re
        end
    end
    return acc
end

"""
    _adj_cvis_turbo!(grad_xw, adj, kx, ky, k2_inv_im, proj_west, proj_north)

The same sum as the scalar branch of [`compute_adjoint_cvis!`](@ref), vectorised over the UV
index.

THE COMPLEX ARITHMETIC IS EXPANDED, because `@turbo` does not vectorise `Complex`. The scalar
version accumulates `adj·k2·s·cis(-πkdc)·cr` and takes the real part at the end, so only the
real part need ever be formed. With `k2_inv_im[k] = -i·c` purely imaginary by construction and
`adj[k] = ar + i·ai`,

    adj·k2 = c·ai − i·c·ar
    real[(c·ai − i·c·ar)·(cos b + i·sin b)] = c·(ai·cos b + ar·sin b)

so the whole contribution is `s·cr·c·(ai·cos b + ar·sin b)` — real throughout.

THREADING IS SIMPLER THAN THE FORWARD'S. There each `F[k]` is touched by every tessel, which
forced explicit chunking and per-chunk accumulators; here each `p` writes `grad_xw[p]` and
nothing else, so `@threads` over `p` is already safe and the vectorised axis is `k` inside it.
"""
function _adj_cvis_turbo!(grad_xw::Vector{T}, adj::Vector{Complex{T}},
    kx::Vector{T}, ky::Vector{T}, k2_inv_im::Vector{Complex{T}},
    proj_west::AbstractMatrix{T}, proj_north::AbstractMatrix{T}) where T

    nuv = length(kx); npix = size(proj_west, 1)
    # The three real vectors the kernel needs, formed once rather than per tessel.
    c  = Vector{T}(undef, nuv); ar = Vector{T}(undef, nuv); ai = Vector{T}(undef, nuv)
    @inbounds @simd for k in 1:nuv
        c[k]  = -imag(k2_inv_im[k])
        ar[k] = real(adj[k]); ai[k] = imag(adj[k])
    end
    Threads.@threads for p in 1:npix
        @inbounds begin
            x1 = proj_west[p,1];  x2 = proj_west[p,2];  x3 = proj_west[p,3];  x4 = proj_west[p,4]
            y1 = proj_north[p,1]; y2 = proj_north[p,2]; y3 = proj_north[p,3]; y4 = proj_north[p,4]
            d1x = x2-x1; d1y = y2-y1; c1x = x2+x1; c1y = y2+y1
            d2x = x3-x2; d2y = y3-y2; c2x = x3+x2; c2y = y3+y2
            d3x = x4-x3; d3y = y4-y3; c3x = x4+x3; c3y = y4+y3
            d4x = x1-x4; d4y = y1-y4; c4x = x1+x4; c4y = y1+y4
            g = zero(T)
            @turbo for k in 1:nuv
                kxk = kx[k]; kyk = ky[k]; ck = c[k]; ark = ar[k]; aik = ai[k]
                a1 = kxk*d1x + kyk*d1y; b1 = -T(π)*(kxk*c1x + kyk*c1y)
                a2 = kxk*d2x + kyk*d2y; b2 = -T(π)*(kxk*c2x + kyk*c2y)
                a3 = kxk*d3x + kyk*d3y; b3 = -T(π)*(kxk*c3x + kyk*c3y)
                a4 = kxk*d4x + kyk*d4y; b4 = -T(π)*(kxk*c4x + kyk*c4y)
                p1 = T(π)*a1; p2 = T(π)*a2; p3 = T(π)*a3; p4 = T(π)*a4
                s1 = ifelse(a1 == zero(T), one(T), sin(p1)/p1)
                s2 = ifelse(a2 == zero(T), one(T), sin(p2)/p2)
                s3 = ifelse(a3 == zero(T), one(T), sin(p3)/p3)
                s4 = ifelse(a4 == zero(T), one(T), sin(p4)/p4)
                r1 = s1*(kyk*d1x - kxk*d1y); r2 = s2*(kyk*d2x - kxk*d2y)
                r3 = s3*(kyk*d3x - kxk*d3y); r4 = s4*(kyk*d4x - kxk*d4y)
                g += ck * (r1*(aik*cos(b1) + ark*sin(b1)) +
                           r2*(aik*cos(b2) + ark*sin(b2)) +
                           r3*(aik*cos(b3) + ark*sin(b3)) +
                           r4*(aik*cos(b4) + ark*sin(b4)))
            end
            grad_xw[p] = g
        end
    end
    return nothing
end

"""
    _adj_vertices_turbo!(gpw, gpn, adj, kx, ky, k2_inv_im, proj_west, proj_north, xw)

The same sums as the scalar branch of [`compute_adjoint_vertices!`](@ref), vectorised over the
UV index. THE LARGEST of the three kernels a gradient evaluation runs — 49% at nside 3 rising
to 57% at nside 5 — and therefore the one worth vectorising most.

THE ALGEBRA, expanded to real arithmetic because `@turbo` will not take `Complex`. Writing
`base = k2·cis(b)·xw` with `k2 = -i·c` purely imaginary and `b = -π·kdc`,

    base = xw·c·(sin b − i·cos b)  ≡  P + i·Q,   P = xw·c·sin b,  Q = −xw·c·cos b

The scalar version forms four derivative terms and takes `real(adj·…)` of each. Two of them
multiply `base` by a REAL factor and two by a purely IMAGINARY one:

    dE_ddx = base·R1,          R1 = ds·kx·cr + s·ky
    dE_ddy = base·R2,          R2 = ds·ky·cr − s·kx
    dE_dcx = base·(−i·M1),     M1 = π·kx·s·cr
    dE_dcy = base·(−i·M2),     M2 = π·ky·s·cr

so with `adj = ar + i·ai`, and using `dx = x₂−x₁`, `cx = x₂+x₁` for the chain rule,

    real[adj·(∓base·R + base·(−i·M))] = ar·(∓P·R + Q·M) − ai·(∓Q·R − P·M)

which is the form below. `ds`, the derivative of `sinc`, is `(cos(πkdd) − s)/kdd` and zero at
`kdd = 0`; the scalar version branches there and this uses `ifelse`, which blends instead.
"""
function _adj_vertices_turbo!(grad_proj_west::Matrix{T}, grad_proj_north::Matrix{T},
    adj::Vector{Complex{T}}, kx::Vector{T}, ky::Vector{T},
    k2_inv_im::Vector{Complex{T}}, proj_west::AbstractMatrix{T},
    proj_north::AbstractMatrix{T}, xw::Vector{T}) where T

    nuv = length(kx); npix = size(proj_west, 1)
    grad_proj_west .= zero(T); grad_proj_north .= zero(T)
    c  = Vector{T}(undef, nuv); ar = Vector{T}(undef, nuv); ai = Vector{T}(undef, nuv)
    @inbounds @simd for k in 1:nuv
        c[k]  = -imag(k2_inv_im[k])
        ar[k] = real(adj[k]); ai[k] = imag(adj[k])
    end
    Threads.@threads for p in 1:npix
        @inbounds begin
            xwp = xw[p]
            if xwp != zero(T)
                for e in 1:4
                    j1 = e; j2 = ifelse(e == 4, 1, e + 1)
                    dx = proj_west[p,j2]  - proj_west[p,j1]
                    dy = proj_north[p,j2] - proj_north[p,j1]
                    cx = proj_west[p,j2]  + proj_west[p,j1]
                    cy = proj_north[p,j2] + proj_north[p,j1]
                    a1x = zero(T); a1y = zero(T); a2x = zero(T); a2y = zero(T)
                    @turbo for k in 1:nuv
                        kxk = kx[k]; kyk = ky[k]
                        kdd = kxk*dx + kyk*dy
                        kdc = kxk*cx + kyk*cy
                        cr  = kyk*dx - kxk*dy
                        pk  = T(π)*kdd
                        s   = ifelse(kdd == zero(T), one(T), sin(pk)/pk)
                        ds  = ifelse(kdd == zero(T), zero(T), (cos(pk) - s)/kdd)
                        b   = -T(π)*kdc
                        f   = xwp * c[k]
                        P   =  f*sin(b); Q = -f*cos(b)
                        R1  = ds*kxk*cr + s*kyk
                        R2  = ds*kyk*cr - s*kxk
                        M1  = T(π)*kxk*s*cr
                        M2  = T(π)*kyk*s*cr
                        ark = ar[k]; aik = ai[k]
                        # j1 takes −R, j2 takes +R; both take +M.
                        a1x += ark*(-P*R1 + Q*M1) - aik*(-Q*R1 - P*M1)
                        a2x += ark*( P*R1 + Q*M1) - aik*( Q*R1 - P*M1)
                        a1y += ark*(-P*R2 + Q*M2) - aik*(-Q*R2 - P*M2)
                        a2y += ark*( P*R2 + Q*M2) - aik*( Q*R2 - P*M2)
                    end
                    grad_proj_west[p, j1]  += a1x
                    grad_proj_west[p, j2]  += a2x
                    grad_proj_north[p, j1] += a1y
                    grad_proj_north[p, j2] += a2y
                end
            end
        end
    end
    return nothing
end
