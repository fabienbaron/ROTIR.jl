# parametric_gradient.jl
# ---------------------------------------------------------------------------
# Zygote-composable parametric forward model for rapid-rotator interferometry.
#
# Design: expose the expensive / non-AD-friendly pieces as pure PRIMITIVES, each
# with a hand-coded ChainRulesCore.rrule (reusing the existing FD-validated
# adjoints). Zygote then composes them and performs the parameter fan-in
# accumulation automatically. All leaf derivatives are HAND-CODED (no ForwardDiff).
#
# The only complex arithmetic (visibilities) is hidden inside the single REAL
# primitive `interferometric_chi2(xw, proj_west, proj_north, …) -> chi2::Real`, so
# no complex-cotangent conventions are ever exposed to the AD.
#
# Everything is type-generic and preserves the input eltype (Float32 by default).
# Only the FD test harness (demos/test_gradients.jl) runs in Float64.
# ---------------------------------------------------------------------------

using ChainRulesCore

# The intensity model (`intensity`, `planck_and_dT`, and its rrule) lives in
# src/intensity.jl — it is included earlier in ROTIR.jl so the χ² paths in
# oichi2_spheroid.jl / oichi2_binary.jl can share it. Step 1 of the parametric
# chain calls it unchanged.

# ===========================================================================
# von Zeipel temperature map + hand-coded parameter derivatives  ─ Step 3
# ===========================================================================
# star_map[p] = tpole · R[p]^β,  R = g_θ/g_pole,  with (rapid rotator):
#   r_θ = rpole·f(fev·sinθ),  ω = fev·√(8GM/27rpole³),  g_pole = GM/rpole²
#   g_r = −GM/r_θ² + r_θ(ω sinθ)²,  g_tt = ω² r_θ sinθ cosθ,  g_θ = √(g_r²+g_tt²)
# θ (colatitude) is intrinsic (param-independent); sinθ, cosθ are precomputed.
"""
    vonzeipel_map_and_derivs(rpole, fev, β, tpole, sinθ, cosθ; GM=1) ->
        (x, dx_drpole, dx_dfev, dx_dβ, dx_dtpole)

Per-tessel von Zeipel temperature map and its hand-coded analytic derivatives.
"""
function vonzeipel_map_and_derivs(rpole::T, fev::T, β::T, tpole::T,
                                  sinθ::AbstractVector{T}, cosθ::AbstractVector{T};
                                  GM::T = one(T)) where {T}
    n = length(sinθ)
    x       = Vector{T}(undef, n)
    dx_drp  = Vector{T}(undef, n)
    dx_dfev = Vector{T}(undef, n)
    dx_dβ   = Vector{T}(undef, n)
    dx_dtp  = Vector{T}(undef, n)

    ωc      = sqrt(T(8) * GM / (T(27) * rpole^3))    # ω = fev·ωc,  ωc ∝ rpole^(-3/2)
    ω       = fev * ωc
    dω_dfev = ωc
    dω_drp  = -T(1.5) * ω / rpole
    g_pole  = GM / rpole^2
    dlgpole_drp = -T(2) / rpole                      # ∂ln g_pole/∂rpole
    ω2 = ω * ω

    @inbounds for i in 1:n
        s = sinθ[i]; c = cosθ[i]
        a = fev * s
        f, fp = f_rapid_rot_and_deriv(a)             # guarded (pole → (1,0))
        rt = rpole * f
        drt_drp  = f
        drt_dfev = rpole * fp * s

        g_r  = -GM / (rt * rt) + rt * (ω * s)^2
        g_t  = ω2 * rt * s * c
        gθ2  = g_r * g_r + g_t * g_t
        gθ   = sqrt(gθ2)
        R    = gθ / g_pole
        xi   = tpole * R^β
        x[i] = xi

        dx_dtp[i] = xi / tpole
        dx_dβ[i]  = xi * log(R)

        # ∂g_r/∂q = 2GM/r_θ³·∂r_θ + sinθ²(ω²·∂r_θ + 2 r_θ ω ∂ω)
        # ∂g_t/∂q = (2 ω ∂ω r_θ + ω² ∂r_θ) sinθ cosθ
        dgr_drp = T(2)*GM/(rt^3)*drt_drp + s*s*(ω2*drt_drp + T(2)*rt*ω*dω_drp)
        dgt_drp = (T(2)*ω*dω_drp*rt + ω2*drt_drp) * s * c
        dlgθ_drp = (g_r*dgr_drp + g_t*dgt_drp) / gθ2
        dx_drp[i] = xi * β * (dlgθ_drp - dlgpole_drp)

        dgr_dfev = T(2)*GM/(rt^3)*drt_dfev + s*s*(ω2*drt_dfev + T(2)*rt*ω*dω_dfev)
        dgt_dfev = (T(2)*ω*dω_dfev*rt + ω2*drt_dfev) * s * c
        dlgθ_dfev = (g_r*dgr_dfev + g_t*dgt_dfev) / gθ2
        dx_dfev[i] = xi * β * dlgθ_dfev              # ∂ln g_pole/∂fev = 0
    end
    return x, dx_drp, dx_dfev, dx_dβ, dx_dtp
end

"""
    vonzeipel_map(rpole, fev, β, tpole, sinθ, cosθ; GM=1) -> x

Forward von Zeipel temperature map (Zygote primitive; rrule uses the hand-coded
derivatives above).
"""
vonzeipel_map(rpole, fev, β, tpole, sinθ, cosθ; GM = one(eltype(sinθ))) =
    first(vonzeipel_map_and_derivs(rpole, fev, β, tpole, sinθ, cosθ; GM = GM))

function ChainRulesCore.rrule(::typeof(vonzeipel_map), rpole, fev, β, tpole,
                              sinθ, cosθ; GM = one(eltype(sinθ)))
    x, dx_drp, dx_dfev, dx_dβ, dx_dtp =
        vonzeipel_map_and_derivs(rpole, fev, β, tpole, sinθ, cosθ; GM = GM)
    function vonzeipel_pullback(x̄raw)
        x̄ = unthunk(x̄raw)
        return (NoTangent(),
                dot(x̄, dx_drp), dot(x̄, dx_dfev), dot(x̄, dx_dβ), dot(x̄, dx_dtp),
                NoTangent(), NoTangent())
    end
    return x, vonzeipel_pullback
end

# ===========================================================================
# Limb-darkening weight + hand-coded derivatives  ─ Step 3
# ===========================================================================
# μ from nz via the shared `mu_and_dmu` (geometry.jl); ld from `compute_ldmap`'s laws.
"""
    ld_and_derivs(nz, ldtype, ld1, ld2) -> (ld, dld_dnz, dld_dld1, dld_dld2)

Per-tessel limb-darkening map and hand-coded derivatives w.r.t. nz (through μ),
ld1 and ld2. μ = max(nz,0) (shared `mu_and_dmu`); consistent with the forward
`compute_ldmap` via `limb_mu`.
"""
function ld_and_derivs(nz::AbstractVector{T}, ldtype::Integer, ld1::T, ld2::T,
                       ld3::T = zero(T), ld4::T = zero(T)) where {T}
    n = length(nz)
    ld       = Vector{T}(undef, n)
    dld_dnz  = Vector{T}(undef, n)
    dld_dld1 = Vector{T}(undef, n)
    dld_dld2 = Vector{T}(undef, n)
    @inbounds for i in 1:n
        μ, dμ = mu_and_dmu(nz[i])
        if ldtype == 1                       # linear: 1 − ld1(1−μ)
            ld[i]       = one(T) - ld1*(one(T) - μ)
            dld_dμ      = ld1
            dld_dld1[i] = -(one(T) - μ)
            dld_dld2[i] = zero(T)
        elseif ldtype == 2                   # quadratic: 1 − ld1(1−μ) − ld2(1−μ)²
            m1          = one(T) - μ
            ld[i]       = one(T) - ld1*m1 - ld2*m1*m1
            dld_dμ      = ld1 + T(2)*ld2*m1
            dld_dld1[i] = -m1
            dld_dld2[i] = -m1*m1
        elseif ldtype == 4
            # Claret (2000): 1 − a₁(1−μ^½) − a₂(1−μ) − a₃(1−μ^{3/2}) − a₄(1−μ²).
            #
            # All FOUR terms are in the forward map, so this is the same law `compute_ldmap`
            # applies — a gradient path evaluating a two-term approximation of it would be
            # optimising a different model from the one being reported.
            #
            # Derivatives are carried for a₁ and a₂ only, because `build_parametric_logπ`'s θ
            # has room for two LD coefficients and widening it would change a vector every
            # caller of `fit_parametric` passes. a₃ and a₄ are held at their model values,
            # exactly as `ld2` is under the power law.
            sq          = sqrt(μ)
            m32         = μ * sq                       # μ^{3/2}
            ld[i]       = one(T) - ld1*(one(T) - sq) - ld2*(one(T) - μ) -
                          ld3*(one(T) - m32) - ld4*(one(T) - μ*μ)
            # d/dμ = a₁/(2√μ) + a₂ + (3/2)a₃√μ + 2a₄μ. The 1/√μ diverges at the limb, where
            # `dμ` is zero anyway and the product is dropped below.
            dld_dμ      = μ > zero(T) ?
                          ld1/(T(2)*sq) + ld2 + T(1.5)*ld3*sq + T(2)*ld4*μ :
                          ld2 + T(2)*ld4*μ
            dld_dld1[i] = -(one(T) - sq)
            dld_dld2[i] = -(one(T) - μ)
        else                                 # ldtype 3, Hestroffer μ^ld1
            ldi         = μ > zero(T) ? μ^ld1 : zero(T)
            ld[i]       = ldi
            dld_dμ      = μ > zero(T) ? ld1 * μ^(ld1 - one(T)) : zero(T)
            dld_dld1[i] = μ > zero(T) ? ldi * log(μ) : zero(T)
            dld_dld2[i] = zero(T)
        end
        # μ-derivative flows to nz only where dμ ≠ 0 (dμ = 0 on the backside/limb,
        # which also avoids 0·Inf from dld_dμ at μ = 0).
        dld_dnz[i] = iszero(dμ) ? zero(T) : dld_dμ * dμ
    end
    return ld, dld_dnz, dld_dld1, dld_dld2
end

"""
    ld_weight(nz, ldtype, ld1, ld2, ld3 = 0, ld4 = 0) -> ld

Forward LD map (Zygote primitive).
"""
ld_weight(nz, ldtype, ld1, ld2, ld3 = zero(ld1), ld4 = zero(ld1)) =
    first(ld_and_derivs(nz, ldtype, ld1, ld2, ld3, ld4))

function ChainRulesCore.rrule(::typeof(ld_weight), nz, ldtype, ld1, ld2,
                              ld3 = zero(ld1), ld4 = zero(ld1))
    ld, dld_dnz, dld_dld1, dld_dld2 = ld_and_derivs(nz, ldtype, ld1, ld2, ld3, ld4)
    function ld_pullback(l̄draw)
        l̄d = unthunk(l̄draw)
        n̄z  = l̄d .* dld_dnz
        l̄d1 = dot(l̄d, dld_dld1)
        l̄d2 = dot(l̄d, dld_dld2)
        # One tangent per argument, `ld3`/`ld4` included: they are held fixed here (see
        # `ld_and_derivs`), so their tangent is exactly zero rather than absent — a pullback
        # returning too few tangents is a silent shape error inside Zygote.
        return (NoTangent(), n̄z, NoTangent(), l̄d1, l̄d2, ZeroTangent(), ZeroTangent())
    end
    return ld, ld_pullback
end

# ===========================================================================
# Soft visibility weight (sigmoid) primitive
# ===========================================================================
"""
    visibility_weight(nz, κ) -> vw

Soft visibility σ(κ·nz) (Zygote primitive; reuses `sigmoid`/`dsigmoid`).
"""
visibility_weight(nz::AbstractVector{T}, κ) where {T} = sigmoid.(T(κ) .* nz)

function ChainRulesCore.rrule(::typeof(visibility_weight), nz::AbstractVector{T}, κ) where {T}
    κT = T(κ)
    vw = sigmoid.(κT .* nz)
    function vis_pullback(v̄wraw)
        v̄w = unthunk(v̄wraw)
        n̄z = v̄w .* dsigmoid.(κT .* nz) .* κT
        return (NoTangent(), n̄z, NoTangent())
    end
    return vw, vis_pullback
end

# ===========================================================================
# Geometry projection primitive  ─ Step 4
# ===========================================================================
# Reuses projected_vertices_and_derivs (shape_gradient.jl) for BOTH the forward
# projected vertices/normals and the reverse contraction. θ order [rpole, ω(fev),
# inc, PA] matches the rapid-rotator layout there.
"""
    project_geometry(rpole, fev, inc, PA, tessels, t, base_params) ->
        (proj_west, proj_north, nz)

Projected quad vertices (all npix) and line-of-sight normal component, as a
function of the geometry parameters (Zygote primitive).
"""
function project_geometry(rpole, fev, inc, PA, tessels, t, base_params)
    sp = merge(base_params, (rpole = rpole, frac_escapevel = fev,
                             inclination = inc, position_angle = PA))
    pw, pn, _, _, nz, _ = projected_vertices_and_derivs(tessels, sp, t; nparams = 4)
    return pw, pn, nz
end

function ChainRulesCore.rrule(::typeof(project_geometry), rpole, fev, inc, PA,
                              tessels, t, base_params)
    sp = merge(base_params, (rpole = rpole, frac_escapevel = fev,
                             inclination = inc, position_angle = PA))
    pw, pn, dpw, dpn, nz, dnz = projected_vertices_and_derivs(tessels, sp, t; nparams = 4)
    _active(z) = !(z === nothing || z isa ChainRulesCore.AbstractZero)
    function geom_pullback(Δ)
        p̄w = unthunk(Δ[1]); p̄n = unthunk(Δ[2]); n̄z = unthunk(Δ[3])
        T = eltype(dnz)
        usew = _active(p̄w); usen = _active(p̄n); usez = _active(n̄z)
        g = ntuple(4) do j
            s = zero(T)
            if usew
                @inbounds for p in axes(dpw, 1), v in 1:4
                    s += p̄w[p, v] * dpw[p, v, j]
                end
            end
            if usen
                @inbounds for p in axes(dpn, 1), v in 1:4
                    s += p̄n[p, v] * dpn[p, v, j]
                end
            end
            if usez
                @inbounds for p in axes(dnz, 1)
                    s += n̄z[p] * dnz[p, j]
                end
            end
            s
        end
        return (NoTangent(), g[1], g[2], g[3], g[4],
                NoTangent(), NoTangent(), NoTangent())
    end
    return (pw, pn, nz), geom_pullback
end

# ===========================================================================
# Interferometric χ² primitive (real; hides all complex visibility math)  ─ Step 4
# ===========================================================================
# Forward: F = polyFT(xw, proj); flux = Σ polyflux·xw; cvis = F/flux; χ²(V²,T3amp,T3phi).
# rrule pullback replicates the FD-validated adjoint of shape_chi2_fg!:
#   x̄w  via compute_adjoint_cvis!  + flux-normalization correction
#   p̄roj via compute_adjoint_vertices! + shoelace flux correction
"""
    interferometric_chi2(xw, proj_west, proj_north, kx, ky, k2_inv_im, data) -> chi2

Single-epoch interferometric χ² (V² + T3amp + T3phi) from LD/visibility-weighted
per-tessel values `xw` and projected vertices. Real-valued Zygote primitive.
`kx, ky, k2_inv_im` are the pre-scaled UV frequencies (see `precompute_k2_inv_im`).
"""
function interferometric_chi2(xw::AbstractVector{T}, proj_west::AbstractMatrix{T},
                              proj_north::AbstractMatrix{T}, kx::Vector{T}, ky::Vector{T},
                              k2_inv_im::Vector{Complex{T}}, data) where {T}
    nuv = length(kx); npix = length(xw)
    F  = Vector{Complex{T}}(undef, nuv)
    pf = zeros(T, npix)   # zeros: compute_polyflux_and_cvis! skips xw==0 pixels, leaving
                          # them here at 0 (correct: their flux contribution pf·xw is 0)
    compute_polyflux_and_cvis!(F, pf, kx, ky, k2_inv_im, proj_west, proj_north, xw)
    flux = dot(pf, xw)
    cvis = F ./ flux
    v2   = abs2.(cvis[data.indx_v2])
    t3   = cvis[data.indx_t3_1] .* cvis[data.indx_t3_2] .* cvis[data.indx_t3_3]
    t3a  = abs.(t3)
    t3p  = angle.(t3) .* T(180/π)
    return sum(abs2, (v2 .- data.v2) ./ data.v2_err) +
           sum(abs2, (t3a .- data.t3amp) ./ data.t3amp_err) +
           sum(abs2, mod360(t3p .- data.t3phi) ./ data.t3phi_err)
end

function ChainRulesCore.rrule(::typeof(interferometric_chi2),
                              xw::AbstractVector{T}, proj_west::AbstractMatrix{T},
                              proj_north::AbstractMatrix{T}, kx::Vector{T}, ky::Vector{T},
                              k2_inv_im::Vector{Complex{T}}, data) where {T}
    nuv = length(kx); npix = length(xw)
    F  = Vector{Complex{T}}(undef, nuv)
    pf = zeros(T, npix)   # zeros: compute_polyflux_and_cvis! skips xw==0 pixels, leaving
                          # them here at 0 (correct: their flux contribution pf·xw is 0)
    compute_polyflux_and_cvis!(F, pf, kx, ky, k2_inv_im, proj_west, proj_north, xw)
    flux = dot(pf, xw)
    cvis = F ./ flux
    v2model = abs2.(cvis[data.indx_v2])
    t3model = cvis[data.indx_t3_1] .* cvis[data.indx_t3_2] .* cvis[data.indx_t3_3]
    t3amod  = abs.(t3model)
    t3pmod  = angle.(t3model) .* T(180/π)
    chi2 = sum(abs2, (v2model .- data.v2) ./ data.v2_err) +
           sum(abs2, (t3amod .- data.t3amp) ./ data.t3amp_err) +
           sum(abs2, mod360(t3pmod .- data.t3phi) ./ data.t3phi_err)

    function chi2_pullback(c̄raw)
        c̄ = unthunk(c̄raw)                         # real scalar (usually 1)
        # ∂χ²/∂cvis  (same construction as shape_chi2_fg!)
        adj_cvis = zeros(Complex{T}, nuv)
        @inbounds for i in eachindex(data.indx_v2)
            k = data.indx_v2[i]
            adj_cvis[k] += 4*(v2model[i]-data.v2[i])/data.v2_err[i]^2 * conj(cvis[k])
        end
        t3amp_res = 2 .* (t3amod .- data.t3amp) ./ data.t3amp_err.^2
        @inbounds for i in eachindex(data.indx_t3_1)
            k1=data.indx_t3_1[i]; k2=data.indx_t3_2[i]; k3=data.indx_t3_3[i]
            c1=cvis[k1]; c2=cvis[k2]; c3=cvis[k3]
            a1=abs(c1); a2=abs(c2); a3=abs(c3)
            adj_cvis[k1] += t3amp_res[i]*conj(c1)/a1*a2*a3
            adj_cvis[k2] += t3amp_res[i]*conj(c2)/a2*a1*a3
            adj_cvis[k3] += t3amp_res[i]*conj(c3)/a3*a1*a2
        end
        t3phi_res = mod360(t3pmod .- data.t3phi) ./ data.t3phi_err.^2
        @inbounds for i in eachindex(data.indx_t3_1)
            k1=data.indx_t3_1[i]; k2=data.indx_t3_2[i]; k3=data.indx_t3_3[i]
            c1=cvis[k1]; c2=cvis[k2]; c3=cvis[k3]; t3i=t3model[i]
            factor = t3phi_res[i]/abs2(t3i)*conj(t3i)
            adj_cvis[k1] -= T(360/π)*im*factor*c2*c3
            adj_cvis[k2] -= T(360/π)*im*factor*c1*c3
            adj_cvis[k3] -= T(360/π)*im*factor*c1*c2
        end
        adj_F = adj_cvis ./ flux

        grad_xw = Vector{T}(undef, npix)
        compute_adjoint_cvis!(grad_xw, adj_F, kx, ky, k2_inv_im, proj_west, proj_north, pf)
        flux_adj = -dot(xw, grad_xw) / flux
        grad_xw .+= flux_adj .* pf

        gpw = Matrix{T}(undef, npix, 4); gpn = Matrix{T}(undef, npix, 4)
        compute_adjoint_vertices!(gpw, gpn, adj_F, kx, ky, k2_inv_im,
                                  proj_west, proj_north, xw, pf)
        # ∂flux/∂proj shoelace correction (shape_gradient.jl:417-425)
        @inbounds for p in 1:npix
            fa = flux_adj * xw[p]
            for j in 1:4
                jp = mod1(j+1, 4); jm = mod1(j-1, 4)
                gpw[p, j] += fa * T(0.5) * (proj_north[p, jp] - proj_north[p, jm])
                gpn[p, j] += fa * T(0.5) * (proj_west[p, jm]  - proj_west[p, jp])
            end
        end

        return (NoTangent(), c̄ .* grad_xw, c̄ .* gpw, c̄ .* gpn,
                NoTangent(), NoTangent(), NoTangent(), NoTangent())
    end
    return chi2, chi2_pullback
end

# ===========================================================================
# Composed parametric log-likelihood  ─ Step 6
# ===========================================================================
"""
    build_parametric_logπ(data_epochs, tessels, tepochs, base_params;
                          intensity_model=:linear, band=nothing, κ=50, GM=1,
                          tpole_free=false, gravity_law=nothing,
                          logprior=nothing) -> logπ(θ)

Return a Zygote-differentiable closure `θ -> -0.5·χ²(θ) + logprior(θ)` for the
rapid-rotator parametric model. θ = [rpole, ω, inc, PA, β, ld1, ld2]
(+ tpole if `tpole_free`). All epoch/data constants are captured once.

`gravity_law` selects the gravity-darkening law — `:vonzeipel` or `:elr`, see
[`gravity_law_spec`](@ref). The default `nothing` takes it from `base_params.gravity_law`,
so the law is a property of the MODEL and a caller who set it there need not repeat it here;
pass it explicitly only to override. β is θ[5] under either law, and holding it at 1/4 to
recover Espinosa Lara & Rieutord's published exponent is done by leaving `beta` out of the
fit's `free` set, not by changing the law.

Compute the gradient with `Zygote.gradient(logπ, θ)` (load Zygote yourself).
"""
function build_parametric_logπ(data_epochs, tessels, tepochs, base_params;
                               intensity_model::Symbol = :linear, band = nothing,
                               κ = 50, GM = 1, tpole_free::Bool = false,
                               gravity_law = nothing, logprior = nothing)
    T = eltype(tessels.unit_xyz)
    colat = tessels.unit_spherical[:, 5, 2]
    sinθ = T.(sin.(colat)); cosθ = T.(cos.(colat))
    ldtype = base_params.ldtype
    # THE LAW AS A `Val`, resolved once out here. Both laws have the same signature and the
    # same rrule shape, so all the closure needs is which primitive to call — and as a type
    # parameter rather than a captured Symbol, so the call stays statically dispatched and
    # Zygote sees one concrete primitive per closure rather than a branch.
    lawv = Val(gravity_law_name(gravity_law === nothing ? base_params : gravity_law))
    # Claret's third and fourth coefficients, from the model rather than from θ: the θ vector
    # carries two LD coefficients and widening it would change every caller of
    # `fit_parametric`. They are still part of the FORWARD law, so the gradient path evaluates
    # the same limb darkening `compute_ldmap` does.
    ld3_base = T(hasproperty(base_params, :ld3) ? base_params.ld3 : 0)
    ld4_base = T(hasproperty(base_params, :ld4) ? base_params.ld4 : 0)
    tpole_base = T(base_params.tpole)
    κT = T(κ); GMT = T(GM)
    nepochs = length(data_epochs)
    kxs = Vector{Vector{T}}(undef, nepochs)
    kys = Vector{Vector{T}}(undef, nepochs)
    k2s = Vector{Vector{Complex{T}}}(undef, nepochs)
    ts  = T.(tepochs)
    for ep in 1:nepochs
        d = data_epochs[ep]
        kx = T.(d.uv[1, :] .* T(-π/(180*3600000)))
        ky = T.(d.uv[2, :] .* T( π/(180*3600000)))
        kxs[ep] = kx; kys[ep] = ky; k2s[ep] = precompute_k2_inv_im(kx, ky)
    end

    return function logπ(θ)
        R = eltype(θ)
        rpole = θ[1]; fev = θ[2]; inc = θ[3]; PA = θ[4]; β = θ[5]; ld1 = θ[6]; ld2 = θ[7]
        tpole = tpole_free ? θ[8] : R(tpole_base)
        x = gravity_map(lawv, rpole, fev, β, tpole, sinθ, cosθ; GM = GMT)
        Imap = intensity(x, intensity_model, band)
        chi2 = sum(1:nepochs) do ep
            pw, pn, nz = project_geometry(rpole, fev, inc, PA, tessels, ts[ep], base_params)
            ld = ld_weight(nz, ldtype, ld1, ld2, ld3_base, ld4_base)
            vw = visibility_weight(nz, κT)
            xw = Imap .* vw .* ld
            interferometric_chi2(xw, pw, pn, kxs[ep], kys[ep], k2s[ep], data_epochs[ep])
        end
        val = -R(0.5) * chi2
        return logprior === nothing ? val : val + logprior(θ)
    end
end

# ===========================================================================
# The SPHERE's log-posterior  ─ NUTS on a limb-darkened uniform disc
# ===========================================================================
# The parametric log-posterior above IS the rapid rotator: its θ starts `[rpole,
# frac_escapevel, …]` and it calls `vonzeipel_map` directly. A sphere needs its own, and
# crucially a SHORTER θ, because most of the rapid rotator's parameters are not identifiable
# on a sphere at all:
#
#   * `inclination` and `position_angle` — a uniform limb-darkened sphere looks the same from
#     every direction. MEASURED on six lam And epochs at HEALPix 3: swinging the inclination
#     55 degrees moves V² by 4.8e-4 in relative terms and a 137-degree position-angle swing by
#     4.6e-4, against 1.7e-2 for a 2 % change in radius — 36x smaller, and it is HEALPix
#     FACETING rather than physics: the 60->120 degree reflection maps the mesh onto itself and
#     gives 3.9e-7. Sampling these would be sampling discretisation noise over a flat prior.
#   * `tpole` — a uniform map is a constant, and `cvis = F/flux` divides a constant out. It is
#     a pure scale on the visibilities and carries no information in them.
#   * `beta` — gravity darkening on a sphere has nothing to act on.
#
# What is left is what `fit_sphere_ld` has always fitted, and what an interferometer actually
# measures for a single star: the ANGULAR RADIUS and the LIMB DARKENING. So
#
#     θ = [radius, ld…]      with as many LD coefficients as the law reads
#
# and the map is `ones`: not an approximation, but the statement that the temperature cancels.

"""
    project_sphere_geometry(radius, tessels, t, base_params) -> (pw, pn, nz)

Projected vertices and face normals of a SPHERE of angular radius `radius`, differentiable in
`radius`.

`projected_vertices_and_derivs` already carries the sphere's own θ layout —
`[radius, inclination, position_angle]`, `nparams = 3` — so only the wrapper and its pullback
are new. The orientation comes from `base_params` and is NOT differentiated: it is flat to
within the mesh's faceting (see the note above), and a gradient through a flat direction is
numerical noise the sampler would chase.
"""
function project_sphere_geometry(radius, tessels, t, base_params)
    sp = merge(base_params, (radius = radius,))
    pw, pn, _, _, nz, _ = projected_vertices_and_derivs(tessels, sp, t; nparams = 3)
    return pw, pn, nz
end

function ChainRulesCore.rrule(::typeof(project_sphere_geometry), radius,
                              tessels, t, base_params)
    sp = merge(base_params, (radius = radius,))
    pw, pn, dpw, dpn, nz, dnz = projected_vertices_and_derivs(tessels, sp, t; nparams = 3)
    _active(z) = !(z === nothing || z isa ChainRulesCore.AbstractZero)
    function sphere_geom_pullback(Δ)
        p̄w = unthunk(Δ[1]); p̄n = unthunk(Δ[2]); n̄z = unthunk(Δ[3])
        T = eltype(dnz)
        s = zero(T)
        # Column 1 of the sphere layout is `radius`; columns 2 and 3 are the orientation, which
        # is held fixed here and so contributes nothing.
        if _active(p̄w)
            @inbounds for p in axes(dpw, 1), v in 1:4
                s += p̄w[p, v] * dpw[p, v, 1]
            end
        end
        if _active(p̄n)
            @inbounds for p in axes(dpn, 1), v in 1:4
                s += p̄n[p, v] * dpn[p, v, 1]
            end
        end
        if _active(n̄z)
            @inbounds for p in axes(dnz, 1)
                s += n̄z[p] * dnz[p, 1]
            end
        end
        return (NoTangent(), s, NoTangent(), NoTangent(), NoTangent())
    end
    return (pw, pn, nz), sphere_geom_pullback
end

"""
    sphere_param_names(ldtype) -> Vector{String}
    default_sphere_bounds(ldtype) -> (lb, ub)

The sphere's θ layout: `radius` followed by exactly the limb-darkening coefficients the law
reads (`ld_coefficients_used`). A coefficient the law ignores is perfectly unconstrained, so
it is not in the vector — the same rule the panel applies when it greys those fields.

`radius` is bounded strictly away from zero: the visibilities are normalised by the total flux,
so a zero radius is 0/0 and the objective is NaN exactly at the bound.
"""
# The limb-darkening coefficients the sphere's θ can carry: those the law reads AND that have
# a derivative. `ld_and_derivs` differentiates `ld1` and `ld2`; the rrule returns
# `ZeroTangent()` for Claret's `ld3`/`ld4`, so they stay in the forward law and out of θ.
_sphere_nld(ldtype::Integer) = min(length(ld_coefficients_used(ldtype)), 2)

sphere_param_names(ldtype::Integer) =
    ["radius"; [String(s) for s in ld_coefficients_used(ldtype)[1:_sphere_nld(ldtype)]]]

function default_sphere_bounds(ldtype::Integer)
    nld = _sphere_nld(ldtype)
    lb = Float64[1e-3]; append!(lb, fill(-1.0, nld))
    # The schema's LD bounds: `ld1` spans -1..2, because Hestroffer's exponent runs past 1;
    # the remaining coefficients are -1..1.
    ub = Float64[Inf]
    nld >= 1 && push!(ub, 2.0)
    append!(ub, fill(1.0, max(nld - 1, 0)))
    return lb, ub
end

"""
    sphere_free_indices(free, ldtype) -> Vector{Int}

Which entries of the sphere's θ a fit may move, from a list of NAMES (`"radius"`, `"ld1"`, …).

`nothing` means all of them. A name the law does not use — or that has no derivative, which is
`ld3`/`ld4` — is not in θ at all, so asking for it is an error rather than a silent no-op: a
parameter you believe is being sampled and is not is the worst of the three outcomes.
"""
function sphere_free_indices(free, ldtype::Integer)
    names = sphere_param_names(ldtype)
    free === nothing && return collect(1:length(names))
    idx = Int[]
    for f in free
        s = String(f)
        i = findfirst(==(s), names)
        i === nothing && error("sphere_free_indices: `$(s)` is not a sphere parameter; " *
                               "θ is $(names) for ldtype $(ldtype)")
        push!(idx, i)
    end
    return sort!(unique(idx))
end

"""
    build_sphere_logπ(data_epochs, tessels, tepochs, base_params; kwargs...) -> logπ(θ)

The log-posterior of a limb-darkened SPHERE, differentiable by Zygote.

`θ = [radius, ld…]`; see [`sphere_param_names`](@ref) for why it is that short. The
temperature map is `ones` because a uniform map cancels in `cvis = F/flux`, so this returns
the same value for every `tpole` — which is the honest statement, not a shortcut.
"""
function build_sphere_logπ(data_epochs, tessels, tepochs, base_params;
                           κ = 50, logprior = nothing)
    T = eltype(tessels.unit_xyz)
    ldtype = base_params.ldtype
    nld = _sphere_nld(ldtype)
    ld3_base = T(hasproperty(base_params, :ld3) ? base_params.ld3 : 0)
    ld4_base = T(hasproperty(base_params, :ld4) ? base_params.ld4 : 0)
    ld1_base = T(hasproperty(base_params, :ld1) ? base_params.ld1 : 0)
    ld2_base = T(hasproperty(base_params, :ld2) ? base_params.ld2 : 0)
    κT = T(κ)
    nepochs = length(data_epochs)
    kxs = Vector{Vector{T}}(undef, nepochs)
    kys = Vector{Vector{T}}(undef, nepochs)
    k2s = Vector{Vector{Complex{T}}}(undef, nepochs)
    ts  = T.(tepochs)
    for ep in 1:nepochs
        d = data_epochs[ep]
        kx = T.(d.uv[1, :] .* T(-π/(180*3600000)))
        ky = T.(d.uv[2, :] .* T( π/(180*3600000)))
        kxs[ep] = kx; kys[ep] = ky; k2s[ep] = precompute_k2_inv_im(kx, ky)
    end

    return function logπ(θ)
        R = eltype(θ)
        radius = θ[1]
        # ld1 AND ld2 ONLY, and `nld` is capped to match. Claret's law reads four
        # coefficients, but `ld_and_derivs` computes `dld/dld1` and `dld/dld2` and the rrule
        # returns `ZeroTangent()` for the other two — so `ld3` and `ld4` are part of the
        # FORWARD law and have no derivative. Putting them in θ made them flat directions,
        # and that passes a finite-difference check silently because FD finds the same zero.
        # `build_parametric_logπ` takes them from the model for exactly this reason.
        ld1 = nld >= 1 ? θ[2] : R(ld1_base)
        ld2 = nld >= 2 ? θ[3] : R(ld2_base)
        chi2 = sum(1:nepochs) do ep
            pw, pn, nz = project_sphere_geometry(radius, tessels, ts[ep], base_params)
            ld = ld_weight(nz, ldtype, ld1, ld2, ld3_base, ld4_base)
            vw = visibility_weight(nz, κT)
            # UNIFORM map, hence no `Imap` factor: a constant divides out of the normalised
            # visibility, so `tpole` is not a parameter of this problem at all.
            xw = vw .* ld
            interferometric_chi2(xw, pw, pn, kxs[ep], kys[ep], k2s[ep], data_epochs[ep])
        end
        val = -R(0.5) * chi2
        return logprior === nothing ? val : val + logprior(θ)
    end
end

# ===========================================================================
# The ELLIPSOID's log-posterior  ─ NUTS on a triaxial von Zeipel star
# ===========================================================================
# The third of three. Where the sphere's θ is short because almost nothing about a sphere is
# identifiable, an ellipsoid's is the longest of the three, because both halves of the model
# respond to it: the projected GEOMETRY depends on the three radii and the orientation, and the
# temperature MAP depends on the three radii, on β and — only under a non-linear intensity law
# — on `tpole`.
#
#     θ = [rx, ry, rz, inc, PA, β, ld1, ld2]      (+ tpole when `tpole_free`)
#
# A FIXED layout, like the rapid rotator's and unlike the sphere's, because at eight entries an
# ldtype-dependent length buys nothing: `ellipsoid_free_indices` refuses `ld2` when the law does
# not read it, which is where that protection belongs (the panel's `_set_state!` refuses the
# same thing for the same reason).
#
# `tpole` IS A PURE SCALE under `:linear`. The map is `tpole·f_i` with `f` independent of it,
# every step to the visibility is linear, and `cvis = F/flux` divides it out — so freeing it
# there is a flat direction. Under `:planck` the map's CONTRAST changes with `tpole` and it
# becomes identifiable, so `tpole_free` is allowed only there, and refused otherwise rather
# than silently sampled.

"""
    project_ellipsoid_geometry(rx, ry, rz, inc, PA, tessels, t, base_params) -> (pw, pn, nz)

Projected vertices and face normals of a triaxial ellipsoid, differentiable in all five.

`projected_vertices_and_derivs` already carries the ellipsoid's θ layout —
`[rx, ry, rz, inc, PA]`, `nparams = 5` — so this is the wrapper and its pullback only.
"""
function project_ellipsoid_geometry(rx, ry, rz, inc, PA, tessels, t, base_params)
    sp = merge(base_params, (radius_x = rx, radius_y = ry, radius_z = rz,
                             inclination = inc, position_angle = PA))
    pw, pn, _, _, nz, _ = projected_vertices_and_derivs(tessels, sp, t; nparams = 5)
    return pw, pn, nz
end

function ChainRulesCore.rrule(::typeof(project_ellipsoid_geometry), rx, ry, rz, inc, PA,
                              tessels, t, base_params)
    sp = merge(base_params, (radius_x = rx, radius_y = ry, radius_z = rz,
                             inclination = inc, position_angle = PA))
    pw, pn, dpw, dpn, nz, dnz = projected_vertices_and_derivs(tessels, sp, t; nparams = 5)
    _active(z) = !(z === nothing || z isa ChainRulesCore.AbstractZero)
    function ellipsoid_geom_pullback(Δ)
        p̄w = unthunk(Δ[1]); p̄n = unthunk(Δ[2]); n̄z = unthunk(Δ[3])
        T = eltype(dnz)
        usew = _active(p̄w); usen = _active(p̄n); usez = _active(n̄z)
        g = ntuple(5) do j
            s = zero(T)
            if usew
                @inbounds for p in axes(dpw, 1), v in 1:4
                    s += p̄w[p, v] * dpw[p, v, j]
                end
            end
            if usen
                @inbounds for p in axes(dpn, 1), v in 1:4
                    s += p̄n[p, v] * dpn[p, v, j]
                end
            end
            if usez
                @inbounds for p in axes(dnz, 1)
                    s += n̄z[p] * dnz[p, j]
                end
            end
            s
        end
        return (NoTangent(), g[1], g[2], g[3], g[4], g[5],
                NoTangent(), NoTangent(), NoTangent())
    end
    return (pw, pn, nz), ellipsoid_geom_pullback
end

"""
    ellipsoid_map(rx, ry, rz, β, tpole, tessels, base_params) -> Vector

The ellipsoid's von Zeipel temperature map, differentiable in all five arguments.

The forward value is `temperature_map_vonZeipel_ellipsoid`'s, and the pullback is the closed
form in [`temperature_map_vonZeipel_ellipsoid_derivs`](@ref) — which is why that function
exists. The map does NOT depend on the orientation (`r_i` is a norm), so `inc` and `PA` are
not arguments here at all rather than arguments with a zero derivative.
"""
function ellipsoid_map(rx, ry, rz, β, tpole, tessels, base_params)
    sp = merge(base_params, (radius_x = rx, radius_y = ry, radius_z = rz,
                             beta = β, tpole = tpole))
    Tmap, = temperature_map_vonZeipel_ellipsoid_derivs(sp, tessels)
    return Tmap
end

function ChainRulesCore.rrule(::typeof(ellipsoid_map), rx, ry, rz, β, tpole,
                              tessels, base_params)
    sp = merge(base_params, (radius_x = rx, radius_y = ry, radius_z = rz,
                             beta = β, tpole = tpole))
    Tmap, drx, dry, drz, dβ, dtp = temperature_map_vonZeipel_ellipsoid_derivs(sp, tessels)
    function ellipsoid_map_pullback(T̄)
        t̄ = unthunk(T̄)
        return (NoTangent(), dot(t̄, drx), dot(t̄, dry), dot(t̄, drz),
                dot(t̄, dβ), dot(t̄, dtp), NoTangent(), NoTangent())
    end
    return Tmap, ellipsoid_map_pullback
end

"""
    ellipsoid_param_names(; tpole_free = false) -> Vector{String}
    default_ellipsoid_bounds(; tpole_free = false) -> (lb, ub)
    ellipsoid_free_indices(free, ldtype; tpole_free = false) -> Vector{Int}

The ellipsoid's θ layout, its box, and the mapping from names to positions.

The radii are bounded strictly away from zero for the reason
`default_parametric_bounds` gives: the visibilities are normalised by the total flux, so a
zero radius is 0/0 and the objective is NaN exactly at the bound.

`ellipsoid_free_indices` REFUSES a limb-darkening coefficient the current law does not read —
`ld_and_derivs` differentiates only `ld1` and `ld2`, and a coefficient nothing reads is
perfectly unconstrained, so freeing it would add a direction the posterior is flat along.
"""
ellipsoid_param_names(; tpole_free::Bool = false) =
    tpole_free ?
    ["radius_x", "radius_y", "radius_z", "inclination", "position_angle",
     "beta", "ld1", "ld2", "tpole"] :
    ["radius_x", "radius_y", "radius_z", "inclination", "position_angle",
     "beta", "ld1", "ld2"]

function default_ellipsoid_bounds(; tpole_free::Bool = false)
    lb = [1e-3, 1e-3, 1e-3,   0.0, -180.0, 0.0,  -1.0, -1.0]
    ub = [Inf,  Inf,  Inf,  180.0,  180.0, 1.0,   2.0,  1.0]
    if tpole_free
        push!(lb, 0.0); push!(ub, Inf)
    end
    return lb, ub
end

function ellipsoid_free_indices(free, ldtype::Integer; tpole_free::Bool = false)
    names = ellipsoid_param_names(; tpole_free = tpole_free)
    free === nothing && return collect(1:length(names))
    used = ld_coefficients_used(ldtype)
    idx = Int[]
    for f in free
        s = String(f)
        i = findfirst(==(s), names)
        i === nothing && error("ellipsoid_free_indices: `$(s)` is not an ellipsoid " *
                               "parameter; θ is $(names)")
        (s in ("ld1", "ld2") && !(Symbol(s) in used)) &&
            error("ellipsoid_free_indices: `$(s)` is not read by limb-darkening law " *
                  "$(ldtype), so it is perfectly unconstrained; do not free it")
        push!(idx, i)
    end
    return sort!(unique(idx))
end

"""
    build_ellipsoid_logπ(data_epochs, tessels, tepochs, base_params; kwargs...) -> logπ(θ)

The log-posterior of a triaxial von Zeipel ellipsoid, differentiable by Zygote.

`θ = [rx, ry, rz, inc, PA, β, ld1, ld2]`, plus `tpole` when `tpole_free`. Both halves of the
model are differentiated: the geometry through
[`project_ellipsoid_geometry`](@ref) and the temperature map through
[`ellipsoid_map`](@ref), whose pullbacks wrap the analytic derivatives rather than asking
Zygote to trace the mesh construction.

`tpole_free = true` requires `intensity_model = :planck`: under `:linear` the temperature is a
pure multiplicative scale that `cvis = F/flux` divides out, and sampling it would be sampling a
flat direction.
"""
function build_ellipsoid_logπ(data_epochs, tessels, tepochs, base_params;
                              intensity_model::Symbol = :linear, band = nothing,
                              κ = 50, tpole_free::Bool = false, logprior = nothing)
    (tpole_free && intensity_model !== :planck) &&
        error("build_ellipsoid_logπ: `tpole` is a pure scale under intensity_model = " *
              ":$(intensity_model) and divides out of the normalised visibility — it is only " *
              "identifiable under :planck. Either pass intensity_model = :planck or leave " *
              "tpole_free = false.")
    T = eltype(tessels.unit_xyz)
    ldtype = base_params.ldtype
    ld3_base = T(hasproperty(base_params, :ld3) ? base_params.ld3 : 0)
    ld4_base = T(hasproperty(base_params, :ld4) ? base_params.ld4 : 0)
    tpole_base = T(base_params.tpole)
    κT = T(κ)
    nepochs = length(data_epochs)
    kxs = Vector{Vector{T}}(undef, nepochs)
    kys = Vector{Vector{T}}(undef, nepochs)
    k2s = Vector{Vector{Complex{T}}}(undef, nepochs)
    ts  = T.(tepochs)
    for ep in 1:nepochs
        d = data_epochs[ep]
        kx = T.(d.uv[1, :] .* T(-π/(180*3600000)))
        ky = T.(d.uv[2, :] .* T( π/(180*3600000)))
        kxs[ep] = kx; kys[ep] = ky; k2s[ep] = precompute_k2_inv_im(kx, ky)
    end

    return function logπ(θ)
        R = eltype(θ)
        rx = θ[1]; ry = θ[2]; rz = θ[3]; inc = θ[4]; PA = θ[5]
        β  = θ[6]; ld1 = θ[7]; ld2 = θ[8]
        tpole = tpole_free ? θ[9] : R(tpole_base)
        x = ellipsoid_map(rx, ry, rz, β, tpole, tessels, base_params)
        Imap = intensity(x, intensity_model, band)
        chi2 = sum(1:nepochs) do ep
            pw, pn, nz = project_ellipsoid_geometry(rx, ry, rz, inc, PA, tessels,
                                                    ts[ep], base_params)
            ld = ld_weight(nz, ldtype, ld1, ld2, ld3_base, ld4_base)
            vw = visibility_weight(nz, κT)
            xw = Imap .* vw .* ld
            interferometric_chi2(xw, pw, pn, kxs[ep], kys[ep], k2s[ep], data_epochs[ep])
        end
        val = -R(0.5) * chi2
        return logprior === nothing ? val : val + logprior(θ)
    end
end
