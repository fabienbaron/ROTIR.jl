# type3_nufft.jl — a 2-D type-3 NUFFT written for ROTIR's polygon transform and nothing else.
#
# WHAT MAKES THIS PROBLEM SPECIAL, measured rather than assumed:
#
#   * The space-bandwidth product is 12.1 x 12.7 CYCLES and is MESH-INDEPENDENT — it is fixed
#     by the star's angular size (6.4 mas on lam And) and the uv coverage (+-1 cycle/mas), not
#     by HEALPix level. So the intermediate grid is ~72x72 and lives in L2 whatever the mesh.
#   * Because the grid is that small, the binning and sorting a general type-3 does so that
#     threads can own disjoint grid tiles is solving a problem we do not have: each thread
#     takes a private grid and they are summed at the end. MEASURED as ~0.5 ms of FINUFFT's
#     1.5 ms at HEALPix 3.
#   * Sources are Gauss-Legendre nodes inside the projected quads. They are GENERATED, not
#     given, so they are never materialised — generate a node and spread it immediately.
#   * Strengths are REAL and the star is centred, so both type-3 phase factors vanish and the
#     spread grid is real until the FFT.
#   * The uv targets are FIXED for a dataset, so every target-side quantity is plan-time work.
#
# Accuracy and cost against the alternatives, on one lam And epoch (1648 uv, Float64):
#
#   HEALPix   :turbo fwd   FINUFFT fwd   this fwd | :turbo adj   FINUFFT adj   this adj
#      3        0.64 ms      1.75 ms     0.35 ms |   0.28 ms      2.98 ms     0.59 ms
#      5        4.70         2.24        0.56    |   4.13         3.55        0.92
#      6       19.69         5.92        2.01    |  17.32         8.16        3.06
#
# at 2.4e-11 forward and 8.1e-11 adjoint, against FINUFFT's 5.8e-10 at its 1e-9 tolerance.
# The adjoint is the reason this exists: `compute_adjoint_cvis!` has no quadrature branch, so
# a gradient fit or a sampler could not use the mesh-independent route at all.
#
# NO WEAK DEPENDENCIES. Everything here is FFTW plus base Julia — deliberately, because the
# point is to be selectable as a default, and a default must not need LoopVectorization.

using FFTW

# ── the kernel ───────────────────────────────────────────────────────────────────────────
#
# FINUFFT's "exponential of semicircle", which is what makes a compact kernel cheap: no
# special function, one sqrt and one exp.
@inline function es_kernel(z::T, β::T) where {T}
    a = one(T) - z * z
    return a <= zero(T) ? zero(T) : exp(β * (sqrt(a) - one(T)))
end

"""
    es_kernel_ft(t, β) -> Real

The kernel's Fourier transform, `∫_{-1}^{1} es_kernel(z,β) cos(t z) dz`, which is what the
deconvolution steps need.

SUBSTITUTED `z = sin θ`, so the integrand becomes `exp(β(cos θ − 1)) cos(t sin θ) cos θ` on
`[−π/2, π/2]`. The ES kernel has an infinite derivative at `|z| = 1` and a trapezoid rule in
`z` converges at first order there; in `θ` the integrand is analytic and the same rule is
spectrally accurate — 200 points reach machine precision and this is plan-time work anyway.
"""
function es_kernel_ft(t::T, β::T; n::Int = 200) where {T}
    s = zero(T); h = T(π) / n
    for i in 0:n
        θ = -T(π)/2 + i*h
        wq = (i == 0 || i == n) ? T(0.5) : one(T)
        s += wq * exp(β*(cos(θ) - one(T))) * cos(t*sin(θ)) * cos(θ)
    end
    return s * h
end

# PARAMETERISED ON THE FFT PLAN TYPE. Held as `Any` it was the one runtime dispatch JET found
# in either direction — `p.fft! * S` returning `Any` — and it sits once per transform rather
# than in a loop, so the cost is small; the reason to fix it is that a clean report is what
# makes the next audit readable.
"Fill one dimension's target stencil: base index into `tbase[d, :]`, weights into `tw`."
function _t3_target_stencil!(tbase::Matrix{Int32}, tw::Matrix{T}, d::Int, kv::Vector{T},
                   dk::T, nf::Int, w::Int, β::T) where {T}
    @inbounds for j in eachindex(kv)
        g = kv[j] / dk + T(nf)/2           # 0-based grid coordinate
        i0 = ceil(Int, g - T(w)/2)         # first 0-based index of the stencil
        tbase[d, j] = i0
        for q in 1:w
            tw[q, j] = es_kernel(T(2)*(g - (i0 + q - 1))/w, β)
        end
    end
    return nothing
end

"""
    t3_psi_degree(::Type) -> Int

Degree of the `1/ψ̂` Chebyshev fit: 24, which reaches 1.3e-13 — comfortably under the
transform's own 6.9e-10, and the same in both working precisions because the fit and its
evaluation are both in `Float64` whatever `T` is (see the `psicoef` field).
"""
t3_psi_degree(::Type{T}) where {T} = 24

struct Type3Plan{T,FP,W}
    nf::Int                    # fine grid points per dimension
    w::Int                     # kernel width in grid cells
    β::T
    P::T                       # spatial period the grid covers
    dk::T                      # k-grid spacing = 1/P
    # source side: the stage-B correction 1/ψ̂, as an EVEN polynomial rather than a table
    # DELIBERATELY Float64 EVEN WHEN T IS Float32, and this is the one place in the plan that
    # is not parameterised on T. `1/ψ̂` is a smooth scalar prefactor evaluated twice per source
    # against a 121-point stencil, so computing it in double is nearly free — while in single
    # it is the DOMINANT error: Clenshaw's b-values grow like deg² and the fit floors at 1.4e-5
    # against a transform that is otherwise good to 2.7e-6.
    psicoef::Vector{Float64}   # Chebyshev coefficients in x = 2(r/rmax)^2 - 1
    invrmax::Float64           # 1/rmax, rmax = P/2
    # target side — fixed geometry, so all of this is plan-time
    tbase::Matrix{Int32}       # (2, J): base grid index in x and y
    twx::Matrix{T}             # (w, J)
    twy::Matrix{T}             # (w, J)
    phihat::Vector{T}          # (nf,) deconvolution of the SPATIAL spreading kernel
    kcoef::Matrix{T}           # (deg+1, w) Horner coefficients, one column per stencil point
    deg::Int
    fft!::FP
    # REAL, not complex. The strengths are real and the centring makes the prephase trivial,
    # so nothing imaginary enters until the FFT — spreading into a complex grid would move
    # twice the memory to add zeros. MEASURED at HEALPix 3: 0.842 -> 0.775 ms, and 4.64 ->
    # 4.34 at HEALPix 6 — worth having and smaller than it looks, because the spread is
    # latency-bound on the scattered accumulate rather than bandwidth-bound.
    grids::Vector{Array{T,2}}
    scratch::Array{Complex{T},2}
    centred::Array{Complex{T},2}
    realbuf::Array{T,2}
end

"""
    plan_type3(A, kx, ky; T=Float64, w=7, σ=2.0) -> Plan

`A` is the source half-extent (the projected stellar radius, in the same units as the mesh);
`kx`, `ky` the fixed target frequencies.
"""
function plan_type3(A::Real, kx::AbstractVector, ky::AbstractVector;
                       T::Type = Float64, w::Int = 7, σ::Real = 2.0, γ::Real = 2.0,
                       tabN::Int = 0, psideg::Int = 0)
    # `tabN` is accepted and ignored: the 1/ψ̂ table it sized is now a polynomial fit.
    return _plan_type3(T, T(A), collect(T, kx), collect(T, ky), w, T(σ), T(γ), tabN, psideg)
end

# TWO OVERSAMPLINGS, and they are not the same one. `γ` pads the spatial PERIOD beyond the
# source extent, which is what oversamples the k-grid that the targets are interpolated from;
# `σ` sizes the spatial grid against the k-range that has to be reached, which is what controls
# the spreading error. Entangling them — taking the period as a fixed multiple and letting the
# grid follow — capped the accuracy at 1e-7 no matter how wide the kernel got, because the
# k-side was stuck at 1.6x.
function _plan_type3(::Type{T}, A::T, kx::Vector{T}, ky::Vector{T}, w::Int, σ::T, γ::T,
               tabN::Int, psideg_in::Int) where {T}
    # β MUST FOLLOW THE OVERSAMPLING. Hardcoded at 2.30w — FINUFFT's value FOR σ = 2 — a wider
    # grid bought nothing: w = 5 went from 5.365e-5 at σ = 2 to 5.325e-5 at σ = 4 while costing
    # 70 % more time. The standard relation is β = πw(1 − 1/2σ), which gives 2.36w at σ = 2
    # (hence the agreement) and 2.75w at σ = 4.
    β = T(π) * w * (one(T) - one(T)/(2σ))
    # ONE β FOR TWO KERNELS, and that is a real limitation rather than a simplification. The
    # spatial spreading kernel φ is oversampled by σ; the k-side interpolation kernel ψ is
    # oversampled by γ. They share β here, which is only correct when σ == γ — and the
    # measurement says so loudly: with β = πw(1−1/2σ), w = 11 gives 2.1e-11 at σ = γ = 2 but
    # 2.0e-7 at σ = 4, γ = 2. Trading grid size for kernel width — which SHOULD be available
    # here, since the FFT is 1.3 % of runtime — needs the two kernels separated first, each
    # with its own β, polynomial set and deconvolution.
    σ == γ || @warn "σ != γ: the shared kernel parameter is only calibrated for σ == γ" σ γ
    kmax = max(maximum(abs, kx), maximum(abs, ky))
    # The period must hold the sources AND the kernel's overhang, and the grid must reach kmax
    # plus the target kernel's half-width. Solved by taking a period with margin and then
    # sizing the grid; both are tiny, so slack costs nothing.
    P = T(2) * A * γ
    # A FLOOR FOR THE KERNEL OVERHANG, and it is not optional. `nf` sized from the k-reach
    # alone is fine for a well-resolved star but not for a geometry with a small
    # space-bandwidth product, and the no-wrap invariant the spread loop relies on then fails:
    # sources reach `nf/2 ± A/dx = nf/2 ± nf/2γ`, the stencil adds `w/2`, so staying inside
    # [0, nf) needs
    #
    #     nf·(1/2 − 1/2γ) ≥ w/2   ⟹   nf ≥ w / (1 − 1/γ)
    #
    # which at γ = 2 is nf ≥ 2w. Without it a short-baseline dataset (or a bootstrap resample,
    # or a binary component) gave nf = 16 against w = 11 and the spread indexed a 16x16 grid
    # at 257 — a BoundsError, found by the test suite and not by any of the benchmarks, whose
    # geometries all landed at nf = 72.
    nfmin = ceil(Int, w / (1 - 1/γ))
    nf = nextprod((2,3,5), max(ceil(Int, σ * (2*kmax*P) + w), nfmin))
    isodd(nf) && (nf += 1)
    # ASSERTED, not assumed: the spread loop has no bounds check by design, so the invariant it
    # rests on is verified once here instead.
    let dxc = P/nf, reach = nf/2 + A/dxc + w/2
        reach <= nf && A/dxc + w/2 <= nf/2 ||
            error("type-3 grid too small for the kernel: nf = $nf, w = $w, γ = $γ")
    end
    dk = one(T) / P
    dx = P / nf

    # spatial spreading kernel deconvolution, per mode, in FFT order
    phihat = Vector{T}(undef, nf)
    for m in 0:nf-1
        mm = m <= nf÷2 ? m : m - nf            # signed mode
        phihat[m+1] = (w*dx/2) * es_kernel_ft(T(π)*mm*w*dx/P, β)
    end

    # THE STAGE-B CORRECTION AS A POLYNOMIAL. This was a 131073-entry table, and it was the
    # single largest cost in the whole plan: 315 ms of the 338 ms build, because every entry is
    # one `khat` quadrature at 2.4 us. It also put a 1 MB random-access lookup in the hot loop,
    # twice per source.
    #
    # EVEN IN r, which is what makes it cheap: `khat` is even in its argument, so 1/ψ̂ is a
    # function of r² alone. Fitting in `x = 2(r/rmax)² − 1` therefore needs half the degree of
    # a fit in r, and evaluation needs no `sqrt` — one multiply for r², one for the affine map,
    # then Horner.
    #
    # CHEBYSHEV BASIS, NOT MONOMIAL, and the reason is measured. A monomial fit at Chebyshev
    # nodes reaches 4.1e-13 in double but plateaus at 4.5e-5 in single from degree 14 upwards —
    # and that plateau is the EVALUATION, not the fit: the monomial coefficients at this degree
    # are large and alternating, so Horner cancels catastrophically at Float32 precision, and
    # no amount of extra degree or a better solve helps. Chebyshev coefficients decay instead,
    # and Clenshaw costs two FMAs a term rather than Horner's one.
    # SOLVED IN Float64 WHATEVER `T` IS, then narrowed. This is plan-time work, and a
    # Vandermonde at this degree has a condition number around 1e8 — solving it in Float32 gave
    # a fit error of 8.2e-3, three orders of magnitude worse than the transform it feeds, while
    # the same solve in double and narrowed to Float32 lands at the single-precision floor.
    # Extra precision in the PLAN is not the promotion worth avoiding; extra precision in the
    # hot loop is.
    rmax = P/2
    psideg = psideg_in > 0 ? psideg_in : t3_psi_degree(T)
    psicoef = Vector{Float64}(undef, psideg+1)
    let m = psideg+1, Td = Float64
        θs = Td[Td(π)*(2i-1)/(2m) for i in 1:m]          # so x_i = cos(θ_i)
        f  = Vector{Td}(undef, m)
        wd = Td(w); dkd = Td(dk); βd = Td(β); rmd = Td(rmax)
        for i in 1:m
            r = rmd * sqrt((cos(θs[i]) + one(Td))/2)
            f[i] = one(Td) / ((wd*dkd/2) * es_kernel_ft(Td(π)*r*wd*dkd, βd))
        end
        # The discrete Chebyshev transform, straight from the definition — m is at most a few
        # dozen, so there is nothing to gain from an FFT and the direct form has no conventions
        # to get wrong.
        for k in 0:psideg
            acc = zero(Td)
            for i in 1:m
                acc += f[i] * cos(k*θs[i])
            end
            psicoef[k+1] = acc * (k == 0 ? one(Td)/m : Td(2)/m)
        end
    end

    # target stencils — fixed geometry, computed once
    J = length(kx)
    tbase = Matrix{Int32}(undef, 2, J)
    twx = Matrix{T}(undef, w, J); twy = Matrix{T}(undef, w, J)
    # MIRRORS THE SOURCE SPREAD EXACTLY — same `ceil(g - w/2)` base, same normalised argument
    # `2(g − i)/w`. Getting these two out of step is a factor-of-one error, which is to say the
    # answer is uncorrelated with the truth rather than merely inaccurate.
    #
    # TWO EXPLICIT LOOPS, not one loop over `((1, kx[j], twx), (2, ky[j], twy))`. That tuple is
    # heterogeneous — an Int, a scalar and a Matrix — so iterating it boxed on every target and
    # cost 129306 allocations in the plan build for no reason.
    _t3_target_stencil!(tbase, twx, 1, kx, dk, nf, w, β)
    _t3_target_stencil!(tbase, twy, 2, ky, dk, nf, w, β)

    # THE KERNEL AS A POLYNOMIAL, which is the difference between a toy and a spreader.
    # `es` costs a sqrt and an exp, and a w-wide stencil in 2-D needs 2w of them per source —
    # 1.3e6 `exp` calls per evaluation at HEALPix 3, which measured as the dominant cost.
    #
    # For a FIXED stencil position q the kernel value is a smooth function of one variable: the
    # fractional offset x = g − i0 − (w−1)/2, which lies in [−1/2, 1/2) by construction. So fit
    # one polynomial per q, once, and the hot loop becomes `deg` fused multiply-adds across a
    # w-vector — which LLVM vectorises, unlike `exp`.
    deg = w + 3
    kcoef = Matrix{T}(undef, deg+1, w)
    let m = deg+1
        xs = T[cos(T(π)*(2i-1)/(2m))/2 for i in 1:m]        # Chebyshev nodes on [-1/2, 1/2]
        V = [xs[i]^(j-1) for i in 1:m, j in 1:m]
        F = factorize(V)
        for q in 1:w
            rhs = [es_kernel(T(2)*(xs[i] + (w-1)/T(2) - (q-1))/w, β) for i in 1:m]
            kcoef[:, q] = F \ rhs
        end
    end

    grids = [zeros(T, nf, nf) for _ in 1:Threads.nthreads()]
    scratch = zeros(Complex{T}, nf, nf)
    # `ESTIMATE`, NOT `MEASURE`, and this is measured too. `MEASURE` benchmarks candidate
    # algorithms at plan time: 66 ms at nf = 22, 92 at 72, 134 at 120 — paid once per distinct
    # grid size per session, and it showed up as a 113 ms first call that no amount of
    # precompilation could remove, because it is runtime search rather than compilation.
    # What it buys here is NOTHING: the resulting transform is 1.01x / 0.99x / 1.00x the
    # `ESTIMATE` one at those three sizes, because at this scale FFTW's codelets are already
    # optimal and there is nothing to search. The FFT is 1.3 % of the transform anyway.
    fp = plan_fft!(scratch; flags = FFTW.ESTIMATE)
    return Type3Plan{T,typeof(fp),w}(nf, w, β, P, dk, psicoef, one(T)/rmax, tbase, twx, twy, phihat,
                              kcoef, deg, fp, grids, scratch, zeros(Complex{T}, nf, nf),
                              zeros(T, nf, nf))
end

@inline function _t3_psihat_inv(p::Type3Plan{T,FP,W}, r::T) where {T,FP,W}
    u = Float64(r) * p.invrmax
    x = muladd(2.0, u*u, -1.0)              # the even variable; no sqrt at run time
    c = p.psicoef
    @inbounds begin
        b1 = 0.0; b2 = 0.0; tx = 2.0*x
        for k in length(c):-1:2
            b1, b2 = muladd(tx, b1, c[k] - b2), b1
        end
        return T(muladd(x, b1, c[1] - b2))
    end
end

"""
    type3_psihat_fit_error(p) -> (max relative error, at r)

How well the polynomial reproduces `1/ψ̂` over the whole period, measured against the
quadrature it was fitted to. A diagnostic, not part of any transform: the end-to-end accuracy
against the exact closed-form kernel is the real test, and this says whether the fit is the
term that limits it.
"""
function type3_psihat_fit_error(p::Type3Plan{T,FP,W}; n::Int = 4001) where {T,FP,W}
    # THE REFERENCE IS IN Float64, always. Computed in `T` it reported a hard 4.5e-5 floor for
    # every degree and both bases — which was the Float32 `khat` quadrature's own error, not
    # the fit's. A diagnostic whose reference is no better than the thing it measures cannot
    # tell you when to stop.
    worst = 0.0; at = 0.0
    rmax = 1.0/p.invrmax; wd = Float64(p.w); dkd = Float64(p.dk); βd = Float64(p.β)
    for i in 0:n
        r = -rmax + 2*rmax*i/n
        exact = 1.0 / ((wd*dkd/2) * es_kernel_ft(π*r*wd*dkd, βd))
        e = abs(Float64(_t3_psihat_inv(p, T(r))) - exact) / abs(exact)
        e > worst && (worst = e; at = r)
    end
    return worst, at
end


# ── execution ────────────────────────────────────────────────────────────────────────────
#
# THE SCALING IS DERIVED, not fitted. With h_l = c_l/ψ̂(r_l) the corrected strengths:
#
#   G[m]  = Σ_l h_l φ(x_m − r_l)                             (spread)
#   Ĥ(k)  = dx · Σ_m G[m] e^{−2πi k x_m} / φ̂(k)              (FFT + deconvolve)
#   F(k)  = dk · Σ_i Ĥ[i] ψ(k − k_i)                         (gather)
#
# because ∫Ĥ(k')ψ(k−k')dk' is the transform of h(r)·ψ̂(r) = f(r), which is F. The two grid
# factors multiply to dx·dk = 1/nf — the FFT's own normalisation, which is the sanity check
# that the constants are right.
@inline function _t3_es_stencil!(kv, g::T, i0::Int, w::Int, β::T) where {T}
    @inbounds for q in 1:w
        kv[q] = es_kernel(T(2)*(g - (i0 + q - 1))/w, β)
    end
end

"""
    _t3_stencil(C, ::Val{DEG}, x, ::Val{W}) -> NTuple{W,T}

The stencil as an NTuple — what an `SVector{W,T}` is after SROA. No heap array, so the values
stay in registers from the Horner recurrence straight into the spread.

MEASURED against the heap-array + `@turbo` form it replaces: 5.2 % faster at HEALPix 3, 8 % at
4-5 and 35.2 % at 6 — while GIVING UP `@turbo`, which was itself worth 10-14 %. An
`SVector{W,T}` compiles to the same thing, so it would measure the same; the raw tuple is used
only to avoid the dependency.

WHAT STATICARRAYS WOULD ADD BEYOND THIS, AND WHY IT IS NOT TAKEN. The natural next step looks
like making `C` an `SMatrix{DEG+1,W}` and the whole thing one matvec against the power vector
`[1, x, …, x^DEG]`, using StaticArrays' unrolled matmul. Two reasons not to: 165 elements is
past the size where `SMatrix` generates good code, and the matvec form MEASURED 3.15x SLOWER
than Horner (50.4 ns against 16.0) because the power vector costs more than the parallel
contraction saves — Horner already has W independent chains, so there is no ILP to recover. It
is also the monomial basis, which is the conditioning trap that floored the ψ̂ fit at 4.5e-5.
"""
@inline function _t3_stencil(C::Matrix{T}, ::Val{DEG}, x::T, ::Val{W}) where {T,DEG,W}
    return ntuple(Val(W)) do q
        @inbounds begin
            acc = C[DEG+1, q]
            for k in DEG:-1:1
                acc = muladd(acc, x, C[k, q])
            end
            acc
        end
    end
end

"The same w kernel values, by Horner on the fractional offset — no `exp`, and vectorisable."
@inline function _t3_stencil!(kv, C::Matrix{T}, ::Val{DEG}, x::T, ::Val{W}) where {T,DEG,W}
    # `W` AND `DEG` AS TYPE PARAMETERS, so the trip counts are compile-time constants. A
    # profile of the runtime-bound version put 22 % of the whole transform in `range.iterate`
    # and `promotion.==` — the loop machinery itself, not the arithmetic. Making them static
    # lets LLVM unroll both loops: worth 8.7 % at HEALPix 3 rising to 30.7 % at 6, and the
    # output is bit-identical. This is the one place a StaticArrays-style specialisation pays;
    # a library matvec (the shape is `kv = Cᵀ·v(x)`, 165 FMA) would not, because the cost was
    # never the arithmetic.
    # PLAIN `@simd`, NOT `@turbo`. The shape is right for vectorising — w independent
    # polynomials sharing one `x` — but LoopVectorization's cost model divides by zero on the
    # nested Horner recurrence (`DivideError` out of `evaluate_cost_tile!`), so this stays on
    # LLVM's own vectoriser, which handles it. The tiles below are where `@turbo` pays.
    # ACCUMULATOR IN A REGISTER. Written the other way round — `k` outer, `q` inner, with
    # `kv[q] = muladd(kv[q], x, C[k,q])` — this was 62 % of the whole transform, and a profile
    # showed why: `kv` is a heap array, so every one of the `deg` steps loads it, does one FMA
    # and stores it back. The loop was store-forwarding bound, not arithmetic bound. With `q`
    # outer the recurrence lives in a register and there is one store per stencil point instead
    # of `deg` of them.
    @inbounds for q in 1:W
        acc = C[DEG+1, q]
        for k in DEG:-1:1
            acc = muladd(acc, x, C[k, q])
        end
        kv[q] = acc
    end
    return kv
end

"Generic type 3 on an explicit point set — the reference the fused version is checked against."
function type3_points!(out::Vector{Complex{T}}, p::Type3Plan{T,FP,W},
             xs::Vector{T}, ys::Vector{T}, cs::Vector{T}) where {T,FP,W}
    nf = p.nf; w = p.w; β = p.β; dx = p.P/nf
    for G in p.grids; fill!(G, zero(T)); end
    N = length(xs)
    Threads.@threads :static for tid in 1:Threads.nthreads()
        G = p.grids[tid]
        # `t3!` is the GENERIC point-set reference — validated against a direct sum, not a hot
        # path — so it keeps the plain `es` evaluation and heap stencil buffers. The quad
        # routines below use polynomial stencils in tuples instead.
        kvx = Vector{T}(undef, w); kvy = Vector{T}(undef, w)
        hnf = T(nf)/2; hw = T(w)/2
        @inbounds for l in tid:Threads.nthreads():N
            x = xs[l]; y = ys[l]
            f = cs[l] * _t3_psihat_inv(p, x) * _t3_psihat_inv(p, y)
            gx = x/dx + hnf; gy = y/dx + hnf
            i0 = ceil(Int, gx - hw); j0 = ceil(Int, gy - hw)
            _t3_es_stencil!(kvx, gx, i0, w, β); _t3_es_stencil!(kvy, gy, j0, w, β)
            for qy in 1:w
                jj = mod(j0 + qy - 1, nf) + 1
                fy = f * kvy[qy]
                for qx in 1:w
                    ii = mod(i0 + qx - 1, nf) + 1
                    G[ii, jj] += fy * kvx[qx]
                end
            end
        end
    end
    return _t3_finish!(out, p, dx)
end

@inline _t3_oob(i::Int, nf::Int) = (i < 0) | (i >= nf)

"""
    type3_quads!(out, p, proj_west, proj_north, xw, enodes, eweights)

The FUSED path: the quadrature nodes are generated and spread in the same loop, so the
`(xs, ys, fs)` arrays — 60k-107k entries, ~1.5 MB — are never built.

`enodes`/`eweights` are the subdivided Gauss-Legendre rule on `[-1,1]`, exactly what
`ROTIR.build_gauss_samples` builds internally.

The threading is over QUADS, not nodes, which costs nothing to arrange because a quad's nodes
are generated together — and each thread owns a private `nf x nf` grid. That is only affordable
because the grid is tiny: at nf = 64 it is 64 KB a thread, so sixteen of them fit in L2. It is
also what lets the whole of FINUFFT's `setpts` be skipped.
"""
function type3_quads!(out::Vector{Complex{T}}, p::Type3Plan{T,FP,W},
                   pw::AbstractMatrix{T}, pn::AbstractMatrix{T}, xw::AbstractVector{T},
                   enodes::Vector{T}, eweights::Vector{T}) where {T,FP,W}
    nf = p.nf; w = p.w; β = p.β; dx = p.P/nf
    for G in p.grids; fill!(G, zero(T)); end
    npix = size(pw, 1); nw = length(enodes)
    quarter = T(0.25)
    Threads.@threads :static for tid in 1:Threads.nthreads()
        G = p.grids[tid]
        kc = p.kcoef; deg = p.deg; halfw = T(w-1)/2
        vW = Val(W); vD = Val(W+3)
        # CONTIGUOUS BLOCKS per thread, not `tid:Threads.nthreads():npix`. A strided slice makes every
        # thread walk the whole of `pw`/`pn` at stride nthreads, so all of them touch every
        # cache line of the geometry. Worth 16.8 % at HEALPix 3, where each thread gets only a
        # few dozen quads, and nothing measurable above — free either way.
        chunk = cld(npix, Threads.nthreads()); lo = (tid-1)*chunk + 1; hi = min(tid*chunk, npix)
        # TYPED AND HOISTED. `nf/2` and `w/2` are Int/Int, i.e. Float64, and adding one to a
        # Float32 grid coordinate promotes the whole rest of the line — the same trap that had
        # `2π .* kx` promoting the FINUFFT targets.
        hnf = T(nf)/2; hw = T(w)/2
        @inbounds for q in lo:hi
            v1x = pw[q,1]; v2x = pw[q,2]; v3x = pw[q,3]; v4x = pw[q,4]
            v1y = pn[q,1]; v2y = pn[q,2]; v3y = pn[q,3]; v4y = pn[q,4]
            xwq = xw[q]
            xwq == zero(T) && continue
            for ie in 1:nw
                eta = enodes[ie]; we = eweights[ie]
                for ix in 1:nw
                    xi = enodes[ix]; wx = eweights[ix]
                    N1 = (1-xi)*(1-eta)*quarter; N2 = (1+xi)*(1-eta)*quarter
                    N3 = (1+xi)*(1+eta)*quarter; N4 = (1-xi)*(1+eta)*quarter
                    x = N1*v1x + N2*v2x + N3*v3x + N4*v4x
                    y = N1*v1y + N2*v2y + N3*v3y + N4*v4y
                    dxi_x = quarter*((1-eta)*(v2x-v1x) + (1+eta)*(v3x-v4x))
                    dxi_y = quarter*((1-eta)*(v2y-v1y) + (1+eta)*(v3y-v4y))
                    det_x = quarter*((1-xi)*(v4x-v1x) + (1+xi)*(v3x-v2x))
                    det_y = quarter*((1-xi)*(v4y-v1y) + (1+xi)*(v3y-v2y))
                    Jac = dxi_x*det_y - det_x*dxi_y
                    f = xwq * Jac * wx * we * _t3_psihat_inv(p, x) * _t3_psihat_inv(p, y)
                    gx = x/dx + hnf; gy = y/dx + hnf
                    i0 = ceil(Int, gx - hw); j0 = ceil(Int, gy - hw)
                    # NO WRAP-AROUND TEST, because the period carries enough margin that one
                    # cannot happen: sources satisfy |r| <= A, the grid covers 2*A*γ with γ = 2,
                    # so the stencil reaches at most nf/2 ± (A/dx + w/2) — comfortably inside.
                    # The check is a plan-time invariant instead of an inner-loop branch.
                    tvx = _t3_stencil(kc, vD, gx - i0 - halfw, vW)
                    tvy = _t3_stencil(kc, vD, gy - j0 - halfw, vW)
                    # THE WHOLE w x w TILE under one `@turbo`. The grid is REAL, which is what
                    # makes this possible at all: `@turbo` cannot vectorise `Complex`, and it
                    # is the same restriction that forced the hand-expanded real accumulators
                    # in ROTIR's own `_adj_cvis_turbo!`.
                    # No `@turbo`: LoopVectorization cannot index a tuple. Both trip counts
                    # are static, so LLVM unrolls and the tuple reads become register selects.
                    @inbounds for qy in 1:W
                        fy = f * tvy[qy]
                        base = (j0 + qy - 1)*nf + i0
                        for qx in 1:W
                            G[base + qx] += fy * tvx[qx]
                        end
                    end
                end
            end
        end
    end
    return _t3_finish!(out, p, dx)
end

"The half of the pipeline after the grid is filled: reduce, FFT, deconvolve, gather."
function _t3_finish!(out::Vector{Complex{T}}, p::Type3Plan{T,FP,W}, dx::T) where {T,FP,W}
    nf = p.nf; w = p.w
    S = p.scratch
    R = p.grids[1]
    for t in 2:Threads.nthreads(); R .+= p.grids[t]; end
    @inbounds for i in eachindex(R); S[i] = Complex{T}(R[i], zero(T)); end
    p.fft! * S
    C = p.centred
    @inbounds for ly in 0:nf-1, lx in 0:nf-1
        sgn = T(ifelse(iseven(lx + ly), 1, -1))
        C[mod(lx + nf÷2, nf) + 1, mod(ly + nf÷2, nf) + 1] =
            S[lx+1, ly+1] * sgn * dx * dx / (p.phihat[lx+1] * p.phihat[ly+1])
    end
    dk = p.dk
    @inbounds Threads.@threads :static for j in eachindex(out)
        acc = zero(Complex{T})
        i0 = p.tbase[1, j]; j0 = p.tbase[2, j]
        for qy in 1:w
            jj = j0 + qy - 1
            _t3_oob(jj, nf) && continue
            wy = p.twy[qy, j]
            for qx in 1:w
                ii = i0 + qx - 1
                _t3_oob(ii, nf) && continue
                acc += C[ii+1, jj+1] * (p.twx[qx, j] * wy)
            end
        end
        out[j] = acc * dk * dk
    end
    return out
end


# ── the adjoint ──────────────────────────────────────────────────────────────────────────
#
# WHAT IS BEING ADJOINTED. The forward is linear in the per-quad weights:
#
#     F[k] = Σ_p xw[p] · A[k,p],    A[k,p] = Σ_{l∈p} J_l w_l e^{−2πi k·r_l}
#
# and what `shape_chi2_fg!` and `interferometric_chi2`'s rrule need is
#
#     grad_xw[p] = Re( Σ_k adj[k] · A[k,p] )   =   Re( (Aᵀ adj)[p] )
#
# — the TRANSPOSE, not the conjugate transpose, which is what ROTIR's own
# `compute_adjoint_cvis!` computes (it takes `real(acc)` of an unconjugated sum).
#
# THE IMPLEMENTATION IS THE PIPELINE RUN BACKWARDS, step by step, because the transpose of a
# composition is the composition of the transposes in reverse. The gather at the targets
# becomes a spread onto the k-grid; the diagonal weighting is its own transpose; and the FFT
# is SYMMETRIC (the matrix e^{−2πi lm/nf} is symmetric in l and m), so its transpose is itself
# — not the inverse FFT, which is the mistake this arrangement exists to avoid. The centring
# permutation is its own inverse for even `nf`.
#
# The expensive end is now a GATHER over the quadrature nodes rather than a scattered
# accumulate, so it needs no private grids and no reduction: every quad reads the same grid
# and writes only its own entry.
function type3_quads_adj!(grad::Vector{T}, p::Type3Plan{T,FP,W},
                       pw::AbstractMatrix{T}, pn::AbstractMatrix{T},
                       adj::Vector{Complex{T}},
                       enodes::Vector{T}, eweights::Vector{T}) where {T,FP,W}
    nf = p.nf; w = p.w; dx = p.P/nf; dk = p.dk
    C = p.centred; fill!(C, zero(Complex{T}))
    # (1) transpose of the target gather: spread `adj` onto the centred k-grid.
    @inbounds for j in eachindex(adj)
        i0 = p.tbase[1, j]; j0 = p.tbase[2, j]
        a = adj[j] * dk * dk
        for qy in 1:w
            jj = j0 + qy - 1
            _t3_oob(jj, nf) && continue
            ay = a * p.twy[qy, j]
            for qx in 1:w
                ii = i0 + qx - 1
                _t3_oob(ii, nf) && continue
                C[ii+1, jj+1] += ay * p.twx[qx, j]
            end
        end
    end
    # (2) un-centre, deconvolve, undo the centred-grid phase, then the SAME forward FFT.
    S = p.scratch
    @inbounds for ly in 0:nf-1, lx in 0:nf-1
        sgn = T(ifelse(iseven(lx + ly), 1, -1))
        S[lx+1, ly+1] = C[mod(lx + nf÷2, nf) + 1, mod(ly + nf÷2, nf) + 1] *
                        sgn * dx * dx / (p.phihat[lx+1] * p.phihat[ly+1])
    end
    p.fft! * S
    # (3) transpose of the fused spread: gather at each node, contract into its quad.
    # ONLY THE REAL PART IS EVER USED: every step after this gather multiplies by real kernel
    # weights, a real Jacobian and real quadrature weights, so Re(Σ z·a) = Σ Re(z)·a. Taking it
    # here halves the traffic AND gives `@turbo` a real array it can vectorise.
    Sr = p.realbuf
    @inbounds for idx in eachindex(Sr); Sr[idx] = real(S[idx]); end
    npix = size(pw, 1); nw = length(enodes); quarter = T(0.25)
    fill!(grad, zero(T))
    Threads.@threads :static for tid in 1:Threads.nthreads()
        kc = p.kcoef; deg = p.deg; halfw = T(w-1)/2
        vW = Val(W); vD = Val(W+3)
        # CONTIGUOUS BLOCKS per thread, not `tid:Threads.nthreads():npix`. A strided slice makes every
        # thread walk the whole of `pw`/`pn` at stride nthreads, so all of them touch every
        # cache line of the geometry. Worth 16.8 % at HEALPix 3, where each thread gets only a
        # few dozen quads, and nothing measurable above — free either way.
        chunk = cld(npix, Threads.nthreads()); lo = (tid-1)*chunk + 1; hi = min(tid*chunk, npix)
        # TYPED AND HOISTED. `nf/2` and `w/2` are Int/Int, i.e. Float64, and adding one to a
        # Float32 grid coordinate promotes the whole rest of the line — the same trap that had
        # `2π .* kx` promoting the FINUFFT targets.
        hnf = T(nf)/2; hw = T(w)/2
        @inbounds for q in lo:hi
            v1x = pw[q,1]; v2x = pw[q,2]; v3x = pw[q,3]; v4x = pw[q,4]
            v1y = pn[q,1]; v2y = pn[q,2]; v3y = pn[q,3]; v4y = pn[q,4]
            acc = zero(T)
            for ie in 1:nw
                eta = enodes[ie]; we = eweights[ie]
                for ix in 1:nw
                    xi = enodes[ix]; wx = eweights[ix]
                    N1 = (1-xi)*(1-eta)*quarter; N2 = (1+xi)*(1-eta)*quarter
                    N3 = (1+xi)*(1+eta)*quarter; N4 = (1-xi)*(1+eta)*quarter
                    x = N1*v1x + N2*v2x + N3*v3x + N4*v4x
                    y = N1*v1y + N2*v2y + N3*v3y + N4*v4y
                    dxi_x = quarter*((1-eta)*(v2x-v1x) + (1+eta)*(v3x-v4x))
                    dxi_y = quarter*((1-eta)*(v2y-v1y) + (1+eta)*(v3y-v4y))
                    det_x = quarter*((1-xi)*(v4x-v1x) + (1+xi)*(v3x-v2x))
                    det_y = quarter*((1-xi)*(v4y-v1y) + (1+xi)*(v3y-v2y))
                    Jac = dxi_x*det_y - det_x*dxi_y
                    gx = x/dx + hnf; gy = y/dx + hnf
                    i0 = ceil(Int, gx - hw); j0 = ceil(Int, gy - hw)
                    tvx = _t3_stencil(kc, vD, gx - i0 - halfw, vW)
                    tvy = _t3_stencil(kc, vD, gy - j0 - halfw, vW)
                    # TUPLES HERE TOO. This gather read `kvx`/`kvy` — the heap arrays the
                    # forward used to fill — after the stencils moved into tuples, i.e. it read
                    # `undef` memory and produced garbage while still timing fast. It also
                    # removes the last `@turbo` in the file, and with it the LoopVectorization
                    # dependency: this route is meant to be usable as a default, and a default
                    # must not need a weak dependency.
                    sr = zero(T)
                    @inbounds for qy in 1:W
                        base = (j0 + qy - 1)*nf + i0
                        t = zero(T)
                        for qx in 1:W
                            t = muladd(Sr[base + qx], tvx[qx], t)
                        end
                        sr = muladd(t, tvy[qy], sr)
                    end
                    acc += sr * Jac * wx * we * _t3_psihat_inv(p, x) * _t3_psihat_inv(p, y)
                end
            end
            grad[q] = acc
        end
    end
    return grad
end

# ── the ROTIR-facing layer ───────────────────────────────────────────────────────────────

"""
    quadrature_for_type3(proj_west, proj_north, kx, ky) -> (ngauss, nsub)

The Gauss-Legendre rule to integrate each quad with, for THIS transform.

SEPARATE FROM [`quadrature_for`](@ref) on purpose, and the difference is the point. That rule
prefers SUBDIVISION: it holds `ngauss = 4` and raises `nsub` until the phase span across a
sub-cell is under 2.5 rad. This one raises the ORDER instead and subdivides only when the span
outruns the highest rule available.

WHY, MEASURED. Gauss-Legendre is spectrally accurate on an analytic integrand while
subdivision is only algebraic, so order buys accuracy far more cheaply. At HEALPix 3 on lam
And, `ngauss = 6, nsub = 1` matches `ngauss = 4, nsub = 3` to the same 4.5e-11 with 15192
nodes against 60768 — four times fewer, and 2.35x less time. `quadrature_for` could not see
this because it was calibrated against FINUFFT, where every choice floored at its own 1e-9
tolerance and the difference was invisible.

The order is read off the phase span and CALIBRATED, not derived: the classical Gauss bound
`(span/2)^2n/(2n)!` is far too pessimistic here, because the integrand carries a smooth
Jacobian and the error that matters is relative to the whole transform. Measured requirements
for ~5e-11 on lam And: span 5.33 -> 6 points, 1.40 -> 4, 0.70 -> 3.

AND THE MARGIN IS NOT OPTIONAL. The first calibration, `ceil(2.5 + span/2)`, reproduced those
three points exactly and still sat on a boundary: polaris's span of 5.0 fell to 5 points where
the requirement is 6, and the forward error went from 2.5e-11 on lam And to 5.8e-9 there —
worse than FINUFFT's own 2.4e-9 on the same data. Accuracy claims for this route are
SPAN-dependent, so the rule carries a point of slack.
"""
function quadrature_for_type3(proj_west::AbstractMatrix, proj_north::AbstractMatrix,
                              kx::AbstractVector, ky::AbstractVector)
    k2m = zero(float(eltype(kx)))
    @inbounds for i in eachindex(kx)
        v = kx[i]^2 + ky[i]^2
        v > k2m && (k2m = v)
    end
    dmax = zero(float(eltype(proj_west)))
    @inbounds for q in axes(proj_west, 1)
        d1 = hypot(proj_west[q,3] - proj_west[q,1], proj_north[q,3] - proj_north[q,1])
        d2 = hypot(proj_west[q,4] - proj_west[q,2], proj_north[q,4] - proj_north[q,2])
        dmax = max(dmax, d1, d2)
    end
    span = 2π * sqrt(k2m) * dmax
    # SNAPPED UP TO AN ORDER THE TABLE HAS. `_GL_NODES` carries 2,3,4,5,6,8,10 — there is no
    # 7-point or 9-point rule — so the raw formula could return a width `t3_gauss_rule` then
    # refuses. Every geometry I benchmarked happened to land on 3-6 and none of them caught
    # it; the precompile workload's synthetic uv did, immediately. Rounded UP, never down,
    # because down would silently cost accuracy.
    want = clamp(ceil(Int, 3 + span/2), 3, T3_GAUSS_ORDERS[end])
    ng = T3_GAUSS_ORDERS[findfirst(>=(want), T3_GAUSS_ORDERS)]
    ns = max(1, ceil(Int, span / (2 * T3_GAUSS_ORDERS[end])))
    return ng, ns
end

"""
    T3_GAUSS_ORDERS

The Gauss-Legendre orders `_GL_NODES` actually tabulates, sorted — the set
[`quadrature_for_type3`](@ref) may choose from. There is no 7- or 9-point rule, so a formula
that computes a width has to be snapped onto this list rather than used directly.
"""
const T3_GAUSS_ORDERS = (2, 3, 4, 5, 6, 8, 10)

"""
    t3_gauss_rule(ngauss, nsub, T) -> (nodes, weights)

The subdivided Gauss-Legendre rule on `[-1,1]`, as two length-`ngauss·nsub` vectors.

The same rule [`build_gauss_samples`](@ref) forms internally — but returned as the 1-D rule
rather than applied, because this transform generates its nodes inside the spread loop and
never materialises them.
"""
function t3_gauss_rule(ngauss::Int, nsub::Int, ::Type{T} = Float64) where {T}
    haskey(_GL_NODES, ngauss) || error("ngauss must be one of $(sort(collect(keys(_GL_NODES))))")
    n = _GL_NODES[ngauss]; wq = _GL_WEIGHTS[ngauss]
    nw = ngauss * nsub
    en = Vector{T}(undef, nw); ew = Vector{T}(undef, nw)
    inv = one(T) / T(nsub)
    for j in 1:nsub, i in 1:ngauss
        cj = -one(T) + T(2j - 1) * inv
        en[(j-1)*ngauss + i] = cj + T(n[i]) * inv
        ew[(j-1)*ngauss + i] = T(wq[i]) * inv
    end
    return en, ew
end

"""
    TYPE3_PLANS

Cache of [`Type3Plan`](@ref)s, because building one costs 0.43 ms and a fit calls the transform
thousands of times with the SAME uv targets.

The source half-extent is quantised UPWARD to quarter-octaves before it becomes part of the
key. That is not a rounding convenience, it is what makes the cache correct as the star's size
moves during a fit: MEASURED, a plan built at radius 3.2 holds to 6.0e-10 when the star shrinks
to 2.4 but degrades to 7.0e-8 at 4.0 and 1.5e-3 at 5.6, because the k-side oversampling
`γ_eff = P/2A` falls below 2. Quantising up means a growing star lands on a new, larger plan
instead of silently using one that no longer covers it.

Keyed on a hash of the target coordinates rather than their identity: `fused_cvis_parts`
rebuilds `kx`/`ky` on every call, so `objectid` would miss every time.
"""
const TYPE3_PLANS = Dict{Any,Any}()
const TYPE3_PLANS_LOCK = ReentrantLock()
const TYPE3_PLANS_MAX = 32

"Quantise the source half-extent upward to a quarter-octave, so a growing star gets a new plan."
t3_extent_bin(A::Real) = exp2(ceil(log2(float(A)) * 4) / 4)

"""
    type3_plan_for(proj_west, proj_north, kx, ky; w=11, σ=2.0, γ=2.0, T=…) -> Type3Plan

The cached plan for this geometry, built on first use. See [`TYPE3_PLANS`](@ref).
"""
function type3_plan_for(proj_west::AbstractMatrix, proj_north::AbstractMatrix,
                        kx::AbstractVector, ky::AbstractVector;
                        w::Int = 11, σ::Real = 2.0, γ::Real = 2.0,
                        T::Type = Float64)
    A = zero(float(real(eltype(proj_west))))
    @inbounds for i in eachindex(proj_west); A = max(A, abs(proj_west[i])); end
    @inbounds for i in eachindex(proj_north); A = max(A, abs(proj_north[i])); end
    # `Float64` WHATEVER THE MESH IS, for the same measured reason `polyft_cvis_nufft` pins it:
    # single precision buys ~10 % (0.179 ms against 0.197 at HEALPix 3) and costs four orders
    # of accuracy (2.1e-6 against 2.5e-11). Pass `T = Float32` to take that trade deliberately.
    Aq = t3_extent_bin(A)
    key = (T, w, Float64(σ), Float64(γ), Aq, length(kx), hash(kx), hash(ky))
    lock(TYPE3_PLANS_LOCK) do
        p = get(TYPE3_PLANS, key, nothing)
        p === nothing || return p
        # A CAP, not an LRU: the number of distinct (dataset, epoch, extent) triples in a
        # session is small, and a plan is ~1.2 MB. Clearing wholesale is simpler than tracking
        # use order and costs one rebuild.
        length(TYPE3_PLANS) >= TYPE3_PLANS_MAX && empty!(TYPE3_PLANS)
        q = plan_type3(Aq, kx, ky; T = T, w = w, σ = σ, γ = γ)
        TYPE3_PLANS[key] = q
        return q
    end
end

"""
    type3_cvis(proj_west, proj_north, xw, kx, ky; ngauss=nothing, nsub=nothing, w=11)

Complex visibilities at scattered uv points by the type-3 route — the drop-in counterpart of
[`polyft_cvis_nufft`](@ref), returning the same unnormalised sum for the caller to divide by
the flux.

`ngauss`/`nsub` default to [`quadrature_for_type3`](@ref), which reads them off the geometry.
"""
function type3_cvis(proj_west::AbstractMatrix, proj_north::AbstractMatrix,
                    xw::AbstractVector, kx::AbstractVector, ky::AbstractVector;
                    ngauss::Union{Nothing,Int} = nothing,
                    nsub::Union{Nothing,Int} = nothing, w::Int = 11)
    ag, as = quadrature_for_type3(proj_west, proj_north, kx, ky)
    ng = ngauss === nothing ? ag : ngauss
    ns = nsub === nothing ? as : nsub
    p = type3_plan_for(proj_west, proj_north, kx, ky; w = w)
    T = _t3_plan_type(p)
    en, ew = t3_gauss_rule(ng, ns, T)
    out = Vector{Complex{T}}(undef, length(kx))
    return type3_quads!(out, p, _t3_as(T, proj_west), _t3_as(T, proj_north),
                        _t3_as(T, xw), en, ew)
end

"""
    type3_cvis_adj!(grad_xw, proj_west, proj_north, adj, kx, ky; ngauss=nothing, nsub=nothing, w=11)

`grad_xw[p] = Re(Σ_k adj[k] ∂F[k]/∂xw[p])` by the type-3 route — the adjoint that
[`compute_adjoint_cvis!`](@ref) computes for the exact operator, computed here for the
quadrature one. They agree to 8.1e-11.
"""
function type3_cvis_adj!(grad_xw::AbstractVector, proj_west::AbstractMatrix,
                         proj_north::AbstractMatrix, adj::AbstractVector,
                         kx::AbstractVector, ky::AbstractVector;
                         ngauss::Union{Nothing,Int} = nothing,
                         nsub::Union{Nothing,Int} = nothing, w::Int = 11)
    ag, as = quadrature_for_type3(proj_west, proj_north, kx, ky)
    ng = ngauss === nothing ? ag : ngauss
    ns = nsub === nothing ? as : nsub
    p = type3_plan_for(proj_west, proj_north, kx, ky; w = w)
    T = _t3_plan_type(p)
    en, ew = t3_gauss_rule(ng, ns, T)
    g = grad_xw isa Vector{T} ? grad_xw : Vector{T}(undef, length(grad_xw))
    type3_quads_adj!(g, p, _t3_as(T, proj_west), _t3_as(T, proj_north),
                     _t3_as(Complex{T}, adj), en, ew)
    g === grad_xw || copyto!(grad_xw, g)
    return grad_xw
end

_t3_plan_type(::Type3Plan{T}) where {T} = T
_t3_as(::Type{S}, A::AbstractArray{S}) where {S} = A isa Array ? A : Array(A)
_t3_as(::Type{S}, A::AbstractArray) where {S} = convert(Array{S}, A)
