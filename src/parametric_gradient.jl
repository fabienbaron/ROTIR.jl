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
function ld_and_derivs(nz::AbstractVector{T}, ldtype::Integer, ld1::T, ld2::T) where {T}
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
    ld_weight(nz, ldtype, ld1, ld2) -> ld

Forward LD map (Zygote primitive).
"""
ld_weight(nz, ldtype, ld1, ld2) = first(ld_and_derivs(nz, ldtype, ld1, ld2))

function ChainRulesCore.rrule(::typeof(ld_weight), nz, ldtype, ld1, ld2)
    ld, dld_dnz, dld_dld1, dld_dld2 = ld_and_derivs(nz, ldtype, ld1, ld2)
    function ld_pullback(l̄draw)
        l̄d = unthunk(l̄draw)
        n̄z  = l̄d .* dld_dnz
        l̄d1 = dot(l̄d, dld_dld1)
        l̄d2 = dot(l̄d, dld_dld2)
        return (NoTangent(), n̄z, NoTangent(), l̄d1, l̄d2)
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
                          tpole_free=false, logprior=nothing) -> logπ(θ)

Return a Zygote-differentiable closure `θ -> -0.5·χ²(θ) + logprior(θ)` for the
rapid-rotator parametric model. θ = [rpole, ω, inc, PA, β, ld1, ld2]
(+ tpole if `tpole_free`). All epoch/data constants are captured once.

Compute the gradient with `Zygote.gradient(logπ, θ)` (load Zygote yourself).
"""
function build_parametric_logπ(data_epochs, tessels, tepochs, base_params;
                               intensity_model::Symbol = :linear, band = nothing,
                               κ = 50, GM = 1, tpole_free::Bool = false,
                               logprior = nothing)
    T = eltype(tessels.unit_xyz)
    colat = tessels.unit_spherical[:, 5, 2]
    sinθ = T.(sin.(colat)); cosθ = T.(cos.(colat))
    ldtype = base_params.ldtype
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
        x = vonzeipel_map(rpole, fev, β, tpole, sinθ, cosθ; GM = GMT)
        Imap = intensity(x, intensity_model, band)
        chi2 = sum(1:nepochs) do ep
            pw, pn, nz = project_geometry(rpole, fev, inc, PA, tessels, ts[ep], base_params)
            ld = ld_weight(nz, ldtype, ld1, ld2)
            vw = visibility_weight(nz, κT)
            xw = Imap .* vw .* ld
            interferometric_chi2(xw, pw, pn, kxs[ep], kys[ep], k2s[ep], data_epochs[ep])
        end
        val = -R(0.5) * chi2
        return logprior === nothing ? val : val + logprior(θ)
    end
end
