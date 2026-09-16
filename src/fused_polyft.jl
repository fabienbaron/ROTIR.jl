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
fused_cvis(x, star, data; intensity_model::Symbol = :linear, band = nothing) =
    ((F, flux) = fused_cvis_parts(x, star, data; intensity_model = intensity_model,
                                  band = band); F ./ flux)

"""
    fused_cvis_parts(x, star, data; intensity_model, band, extra_weights) -> (F, flux)

The two HALVES of the matrix-free visibility computation: the unnormalised polygon transform
and the total flux it should be divided by.

[`fused_cvis`](@ref) is just their quotient. They are exposed separately because a BINARY
cannot use the quotient: two components combine as `(F1 + F2·phase) / (flux1 + flux2)`, which
needs each star's transform and flux before either is normalised. Without this,
[`binary_cvis`](@ref) had no matrix-free route at all and every binary χ² had to build the
dense `nuv × npix` matrix through `setup_oi!` — twice, once per component.

`extra_weights` multiplies the per-tessel weights after limb darkening and soft visibility,
which is where mutual occultation enters.
"""
function fused_cvis_parts(x, star, data; intensity_model::Symbol = :linear, band = nothing,
                          extra_weights = nothing)
    T = eltype(star.proj_west)
    indx = star.index_quads_visible
    I = intensity_model === :linear ? x : intensity(x, intensity_model, band)
    xw = T.(I[indx] .* star.vis_weights[indx] .* star.ldmap[indx])
    extra_weights === nothing || (xw = xw .* T.(extra_weights[indx]))
    pjx = star.proj_west[indx, :]; pjy = star.proj_north[indx, :]
    kx = T.(data.uv[1, :]) * T(-π / (180 * 3600000))
    ky = T.(data.uv[2, :]) * T( π / (180 * 3600000))
    k2 = precompute_k2_inv_im(kx, ky)
    F = Vector{Complex{T}}(undef, length(kx))
    pf = zeros(T, length(indx))
    mpjx = Matrix(pjx); mpjy = Matrix(pjy)
    bk = polyft_backend()
    if bk === :nufft || bk === :t3
        # `polyflux` is the shoelace area and is needed for the flux normalisation whichever
        # backend runs, so it is computed here rather than inside the visibility kernel.
        @inbounds for q in 1:length(indx)
            pf[q] = T(0.5) * (mpjx[q,1]*mpjy[q,2] - mpjx[q,2]*mpjy[q,1] +
                              mpjx[q,2]*mpjy[q,3] - mpjx[q,3]*mpjy[q,2] +
                              mpjx[q,3]*mpjy[q,4] - mpjx[q,4]*mpjy[q,3] +
                              mpjx[q,4]*mpjy[q,1] - mpjx[q,1]*mpjy[q,4])
        end
        Fn = bk === :t3 ? type3_cvis(mpjx, mpjy, xw, kx, ky) :
                          polyft_cvis_nufft(mpjx, mpjy, xw, kx, ky)
        return Complex{T}.(Fn), dot(pf, xw)
    end
    compute_polyflux_and_cvis!(F, pf, kx, ky, k2, mpjx, mpjy, xw)
    return F, dot(pf, xw)
end

"""
    POLYFT_BACKEND

The PROCESS-WIDE forward kernel: `:t3` (default), `:nufft`, `:turbo`, or `:scalar`. A task can
override it for its own duration with [`with_polyft_backend`](@ref); read the answer that
applies here and now with [`polyft_backend`](@ref) rather than this Ref.

All three compute the SAME quantity and are asserted against each other in
`test/test_fused_polyft.jl`. They differ in how, and therefore in how the cost scales:

| backend   | method                                   | HEALPix 3 | 4 | 5 |
|-----------|------------------------------------------|-----------|---|---|
| `:scalar` | the reference, plain Julia               | 119 ms | 425 ms | 1472 ms |
| `:turbo`  | the same sum, SIMD transcendentals       | 6.9 ms | 21.8 ms | 86 ms |
| `:nufft`  | Gauss quadrature + type-3 NUFFT          | 2.6 ms | 3.8 ms | 4.1 ms |
| `:t3`     | the same quadrature, ROTIR's own type 3  | 0.20 ms | 0.28 ms | 0.62 ms |

`:scalar` and `:turbo` are EXACT — the closed-form polygon transform, no parameters.
Both quadrature routes are accurate to a tolerance rather than exactly, and their cost barely
moves with the mesh because it is set by the source and target counts rather than their
product.

`:t3` IS THE DEFAULT, and it is the default over `:nufft` on measurement rather than novelty:
it is faster at every mesh level and in both precisions (0.20 ms against 1.7 at HEALPix 3, 2.3
against 6.0 at HEALPix 6), MORE accurate on both test datasets (2.5e-11 against 5.8e-10 on
lam And, 1.1e-10 against 2.4e-9 on polaris), needs no weak dependency, and is the only
quadrature route with an adjoint. `:nufft` remains as the independent cross-check — it is a
different implementation of the same idea through FINUFFT, which is what makes it worth
keeping.

`:t3` is the same quadrature through `src/type3_nufft.jl` instead of FINUFFT. It is faster than
every other backend in double at every mesh level, more accurate than `:nufft` (2.5e-11 against
5.8e-10), needs no weak dependency — and is the ONLY route with an adjoint, so it is the only
one a gradient fit or a sampler can use. In SINGLE precision the picture differs: `:turbo`
gains 3.5-5.6x from Float32 where `:t3` gains only 10-25 %, so `:turbo` wins the Float32
adjoint at HEALPix 3-4 and `:t3` only from 5.

`:turbo` is what to fall back on if a mesh is coarse enough to worry about the quadrature (see
`polyft_cvis_nufft`: a 2-point rule is only good to 1.6e-3 at HEALPix 3), and `:scalar` is the
definition the other two are tested against — worth being able to select without rebuilding
when a fit looks wrong.
"""
const POLYFT_BACKEND = Ref(:t3)

"""
    POLYFT_BACKEND_SCOPE

A per-TASK override of [`POLYFT_BACKEND`](@ref); `nothing` — the default — means "no override,
use the process-wide setting". Set it with [`with_polyft_backend`](@ref) and read the answer
with [`polyft_backend`](@ref); nothing should read either this or `POLYFT_BACKEND` directly.

WHY A SCOPE AND NOT JUST THE REF. The right kernel is a property of the CODE PATH, not of the
caller's taste, and the two paths disagree:

  * `fused_cvis` / `fused_cvis_parts` — one χ² for a table, a derivative-free trial point —
    has a `:nufft` branch, and `:nufft` is nearly mesh-independent there (2.6/3.8/4.1 ms at
    HEALPix 3/4/5).
  * `interferometric_chi2` and `shape_chi2_fg!` go through `_cvis_forward!` and the two
    adjoints, which have NO `:nufft` branch — they run `:scalar` unless `:turbo` is selected.
    That is where a gradient fit and a sampler live, and `:turbo` takes the three kernels a
    gradient evaluation runs from 11.4 ms to 0.7 ms at HEALPix 3.

So a fit wants to select `:turbo` for its own duration WITHOUT changing what the rest of the
process computes: a GUI holds its worker on one thread while the event loop keeps answering on
another, and flipping a global Ref under a running optimiser changes the objective function
mid-line-search — value and gradient from different kernels. A scope is task-local and is
inherited by tasks started inside it (`Threads.@threads` inside the kernels included), which is
exactly the lifetime wanted.

The one place it does NOT reach is another PROCESS: a distributed Pigeons run would see the
process default on its workers. Neither does the Ref, so this is not a regression.
"""
const POLYFT_BACKEND_SCOPE = Base.ScopedValues.ScopedValue{Union{Nothing, Symbol}}(nothing)

"""
    polyft_backend() -> Symbol

Which forward kernel to run here and now: the innermost [`with_polyft_backend`](@ref) scope if
there is one, and [`POLYFT_BACKEND`](@ref) otherwise. Every kernel dispatch reads this.

The read is one scope lookup per KERNEL INVOCATION — once per epoch per evaluation, not per
element or per uv point — against milliseconds of kernel, so it does not register.
"""
function polyft_backend()
    s = POLYFT_BACKEND_SCOPE[]
    return s === nothing ? POLYFT_BACKEND[] : s
end

"""
    with_polyft_backend(f, kind::Symbol)

Run `f()` with `kind` as the forward kernel, for this task and any task it starts, then restore
whatever was in force. See [`POLYFT_BACKEND_SCOPE`](@ref) for why a fit does this rather than
assigning to [`POLYFT_BACKEND`](@ref).

    with_polyft_backend(:turbo) do
        fit_parametric(data, tess, tepochs, base; free = names)
    end
"""
with_polyft_backend(f, kind::Symbol) =
    Base.ScopedValues.with(f, POLYFT_BACKEND_SCOPE => kind)

# `_cvis_turbo!` is provided by ext/ROTIRLoopVectorizationExt.jl. Declared here so that
# `:turbo` is a name this package knows and the failure without LoopVectorization is a
# sentence rather than a MethodError on an underscored internal.
function _cvis_turbo! end

# The two ADJOINT kernels, same arrangement: declared here so `:turbo` is a name this package
# knows, defined in ext/ROTIRLoopVectorizationExt.jl.
#
# WHY THEY MATTER MORE THAN THE FORWARD. MEASURED on one lam And epoch, the three kernels a
# gradient evaluation runs:
#
#     nside    forward   adj_cvis   adj_verts   vertex share
#       3      1.5 ms     3.6 ms      4.9 ms        49%
#       4      5.2 ms     9.3 ms     18.1 ms        55%
#       5     18.0 ms    32.7 ms     67.3 ms        57%
#
# So the forward is 15-19% of a gradient and the two adjoints are the rest. `:turbo` on the
# forward alone was capped at that; this is where a gradient fit's time actually is.
function _adj_cvis_turbo! end
function _adj_vertices_turbo! end

"""
    TURBO_OK

Set to `true` by `ROTIRLoopVectorizationExt.__init__` when that extension loads. See
[`turbo_available`](@ref).
"""
const TURBO_OK = Ref(false)

"""
    turbo_available() -> Bool

Whether the `:turbo` kernel will work here, i.e. whether `using LoopVectorization` has loaded
ROTIRLoopVectorizationExt.

A CACHED FLAG, not `!isempty(methods(_cvis_turbo!))`. That spelling asks the runtime's method
table on every call, which is a reflection call returning `Any` — JET flags it as the only
runtime dispatch left anywhere in the differentiable log-posterior, reported from inside
`_cvis_forward!`, the innermost visibility kernel.

To be accurate about what that cost: the call sits inside `if polyft_backend() === :turbo`,
and the default backend is `:t3`, so the default path never reaches it. JET sees it because
it analyses both branches. The fix is worth making for the `:turbo` path, which the GUI does
offer, and because a clean report is what makes the next audit readable — not because the
default was paying for it.

The extension SETS the flag rather than this function asking, which also removes a staleness
trap: the GUI loads LoopVectorization on demand mid-session, so anything cached at ROTIR's own
load time would have been wrong for the rest of that session.
"""
turbo_available() = TURBO_OK[]

function _cvis_forward!(F, kx, ky, k2, pjx, pjy, xw)
    if polyft_backend() === :turbo
        turbo_available() ||
            error("the :turbo kernel needs LoopVectorization loaded: add " *
                  "`using LoopVectorization` to this session. It is a weak dependency " *
                  "because loading it costs 1.8 s of GUI startup by invalidating OITOOLS' " *
                  "precompiled plot pipeline, and `:t3` — the default — is faster anyway.")
        # `invokelatest`, and it is NOT optional. `_cvis_turbo!` gets its methods when
        # LoopVectorization is loaded, and the GUI loads it ON DEMAND — in the middle of a
        # session, from a callback the Qt event loop dispatches. Methods added after the
        # calling frame's world age are invisible to it, so a direct call lands on the stub
        # and throws a MethodError that the χ² path catches and reports as "χ² failed": the
        # backend appears to switch while every number comes back empty. The GUI's event loop
        # makes this permanent rather than transient, since every later callback still runs
        # in the world age fixed when `QML.exec()` was entered.
        #
        # The cost is one dynamic dispatch per VISIBILITY COMPUTATION — not per element — so
        # against milliseconds of kernel it does not register.
        return Base.invokelatest(_cvis_turbo!, F, kx, ky, k2, pjx, pjy, xw)
    end
    return _cvis_scalar!(F, kx, ky, k2, pjx, pjy, xw)
end

# `_cvis_turbo!` lives in src/turbo_polyft.jl, provided by ext/ROTIRLoopVectorizationExt.jl.
# See that file for why it is not here: loading LoopVectorization costs 1.8 s of GUI startup,
# and since a quadrature route became the default it buys a cross-check rather than the speed.

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

    bk = polyft_backend()
    if bk === :turbo && turbo_available()
        # `invokelatest` for the same reason the forward needs it: the GUI can load
        # LoopVectorization mid-session, and methods added after this frame's world age are
        # invisible to a direct call.
        return Base.invokelatest(_adj_cvis_turbo!, grad_xw, adj, kx, ky, k2_inv_im,
                                 proj_west, proj_north)
    end
    # NO `:t3` BRANCH HERE, deliberately, even though `type3_cvis_adj!` exists and is tested.
    # It is the adjoint of the QUADRATURE operator, and it is only meaningful paired with the
    # quadrature FORWARD — but every caller of this function (`shape_chi2_fg!`, the
    # `interferometric_chi2` rrule, `fused_spheroid_chi2_fg`) computes its forward with
    # `compute_polyflux_and_cvis!`, which is the EXACT closed form. Wiring it in would pair an
    # exact value with an approximate gradient of a different operator. It is exported for a
    # future gradient route built on `fused_cvis`, which is where it would be consistent.
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

    if polyft_backend() === :turbo && turbo_available()
        return Base.invokelatest(_adj_vertices_turbo!, grad_proj_west, grad_proj_north, adj,
                                 kx, ky, k2_inv_im, proj_west, proj_north, xw)
    end
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
