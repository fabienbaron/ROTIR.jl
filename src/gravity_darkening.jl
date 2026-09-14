# Gravity darkening beyond von Zeipel: the Espinosa Lara & Rieutord (2011) law.
#
# von Zeipel (1924) gives T_eff ∝ g_eff^β, which follows from assuming the star is BAROTROPIC —
# and barotropy is incompatible with radiative equilibrium in a rotating star (Eddington 1925).
# The law is therefore only valid for slow rotation, and interferometry of fast rotators has
# repeatedly found it OVERESTIMATES the pole-to-equator temperature contrast. The usual patch is
# to fit β and find it below 1/4, which hides a model error inside a free parameter.
#
# ESPINOSA LARA & RIEUTORD (2011), A&A 533, A43 — `papers/aa17252-11.pdf` — drop barotropy and
# assume instead that the flux is ANTI-PARALLEL to effective gravity, `F = -f(r,θ) g_eff`, which
# holds to better than half a degree even at 90% of break-up (their Fig. 1). Imposing
# `∇·F = 0` on a Roche surface then fixes the latitudinal flux profile with NO free parameter,
# and their eq. (31) is
#
#     T_eff = (L/4πσGM)^{1/4} · √(tan ϑ / tan θ) · g_eff^{1/4}
#
# where ϑ(r̃,θ) is an auxiliary angle defined implicitly by their eq. (24),
#
#     cos ϑ + ln tan(ϑ/2) = (1/3) ω² r̃³ cos³θ + cos θ + ln tan(θ/2)                   (24)
#
# `√(tan ϑ / tan θ)` is the whole difference from von Zeipel: it is the extra latitudinal
# structure barotropy misses. As ω → 0, eq. (24) gives ϑ = θ, the factor becomes 1, and eq. (31)
# collapses to eq. (33) — von Zeipel with β = 1/4 exactly.
#
# THE EXPONENT IS KEPT FREE HERE. Eq. (31) fixes it at 1/4 because it assumes a grey radiative
# atmosphere; a real atmosphere, a convective envelope, or a limb-darkening law absorbing part
# of the profile all move it. Writing the flux form of eq. (31) as `T ∝ (F_ω g_eff)^β` recovers
# eq. (31) at β = 1/4 and leaves β available, so von Zeipel and ELR can be compared at the SAME
# number of free parameters — which is what an evidence ratio needs.
#
# WHAT `q` IS, AND WHY IT SIMPLIFIES. Everything ω- and radius-dependent in eq. (24) enters
# through the single product `q = ω² r̃³`, where ω = Ω/Ω_k is the rotation rate in units of the
# Keplerian rate at the equator and r̃ = r/R_e. In ROTIR's parameters that is
#
#     ω  = fev · √(8/27) · f(fev)^{3/2},   r̃ = f(fev·sinθ) / f(fev)
#     ⇒  q = (8/27) · fev² · f(fev·sinθ)³
#
# with `f` the Roche shape factor `f_rapid_rot_and_deriv`. **`f(fev)` cancels**, so `q` — and
# hence the entire ELR correction — depends only on `fev` and θ, never on `rpole`. The polar
# radius enters the map exactly where it does for von Zeipel, through `g_eff/g_pole`.
#
# Two checks the implementation is held to, both independent of the code that computes it:
#
#   * ω → 0 must reproduce `vonzeipel_map` to round-off (eq. 33).
#   * the pole-to-equator ratio must equal ELR eq. (32) in closed form,
#
#         T_e/T_p = √(2/(2+ω²)) (1-ω²)^{1/12} exp(-(4/3) ω²/(2+ω²)³)
#
#     which follows from eqs. (27), (28) and r̃_p = 2/(2+ω²) — the last being eq. (30) at θ = 0,
#     and itself a check on the ω mapping above since it must equal 1/f(fev).

"The Roche rotation rate ω = Ω/Ω_k that ROTIR's `frac_escapevel` corresponds to."
elr_omega(fev::T) where {T} =
    fev * sqrt(T(8) / T(27)) * f_rapid_rot_and_deriv(fev)[1]^T(1.5)

"""
    elr_q_and_deriv(fev, sinθ) -> (q, dq/dfev)

`q = ω² r̃³ = (8/27) fev² f(fev sinθ)³`, the only combination of rotation rate and radius that
ELR eq. (24) depends on, and its derivative.

`f(fev)` cancels out of `ω² r̃³` algebraically, which is why `rpole` does not appear: the ELR
correction is a function of the rotation fraction and the colatitude alone.
"""
function elr_q_and_deriv(fev::T, s::T) where {T}
    c = T(8) / T(27)
    f, fp = f_rapid_rot_and_deriv(fev * s)
    q = c * fev * fev * f^3
    dq = c * fev * f * f * (T(2) * f + T(3) * fev * fp * s)
    return q, dq
end

# ELR eq. (24)'s left-hand side and its ϑ-derivative.
#
#     G(ϑ) = cos ϑ + ln tan(ϑ/2),      G'(ϑ) = cos²ϑ / sin ϑ
#
# G is strictly increasing on (0, π/2], mapping it onto (-∞, 0]. Its derivative VANISHES at
# π/2, which is why the equator gets a closed form instead of a Newton solve.
_elr_G(ϑ::T) where {T} = cos(ϑ) + log(tan(ϑ / 2))
_elr_dG(ϑ::T) where {T} = cos(ϑ)^2 / sin(ϑ)

"""
    elr_flux_factor(q, θ) -> F_ω

The ELR latitudinal flux factor `F_ω = tan²ϑ / tan²θ`, with ϑ from eq. (24).

Three regimes, and the two outer ones are EXACT limits rather than approximations — both were
derived by expanding eq. (24) and both appear in the paper as eqs. (27) and (28):

  * θ → 0:    `F_ω → exp((2/3) q)`            (both sides of eq. 24 go as `1 + ln(angle/2)`)
  * θ → π/2:  `F_ω → (1 - q)^{-2/3}`          (both sides go as `-angle³/3`, giving
                                               δ = ε(1-q)^{1/3} for θ = π/2 - ε, ϑ = π/2 - δ)
  * between:  Newton on eq. (24), from the von Zeipel guess ϑ = θ.

A southern colatitude is folded to its northern mirror first: the star is symmetric about its
equator, so `F_ω(π - θ) = F_ω(θ)`.

The limits are not merely convenient. `tan θ` diverges at the equator and `ln tan(θ/2)`
diverges at the pole, so the ratio is 0/0 and ∞/∞ there respectively; and near π/2 `G'(ϑ) → 0`,
so Newton would take arbitrarily large steps. Away from the ends the solve converges in three
or four iterations because G is monotone and smooth.
"""
function elr_flux_factor(q::T, θ::T) where {T}
    # Non-rotating: ϑ = θ identically, and this is the von Zeipel limit (eq. 33).
    q <= zero(T) && return one(T)
    # THE SOUTHERN HEMISPHERE. Eq. (24) is written for θ ∈ (0, π/2]; a rigidly rotating star is
    # symmetric about its equator, so `F_ω(π - θ) = F_ω(θ)` and the colatitude folds. Without
    # this the `θ >= π/2 - ε` guard below catches EVERY southern colatitude and hands back the
    # equatorial closed form — measured 0.6 % too hot at θ = 0.2 rad on the wrong side, and
    # invisible to a finite-difference check because the analytic derivative makes the same
    # mistake self-consistently. The test that catches it is `T(θ) == T(π - θ)`.
    θ = θ > T(π) / 2 ? T(π) - θ : θ
    ε = T(1e-4)
    θ <= ε        && return exp(T(2) / T(3) * q)
    θ >= T(π)/2 - ε && return (max(one(T) - q, eps(T)))^(-T(2) / T(3))
    rhs = q * cos(θ)^3 / T(3) + cos(θ) + log(tan(θ / 2))
    ϑ = θ
    for _ in 1:40
        r = _elr_G(ϑ) - rhs
        d = _elr_dG(ϑ)
        d <= eps(T) && break
        step = r / d
        # Clamped into the open interval: ϑ must stay in (0, π/2), where G is defined and
        # monotone. An unclamped Newton step can leave it on the first iteration at high q.
        ϑn = clamp(ϑ - step, T(1e-9), T(π) / 2 - T(1e-9))
        abs(ϑn - ϑ) <= eps(T) * max(one(T), ϑ) && (ϑ = ϑn; break)
        ϑ = ϑn
    end
    return (tan(ϑ) / tan(θ))^2
end

"""
    elr_flux_factor_and_dq(q, θ) -> (F_ω, ∂F_ω/∂q)

`elr_flux_factor` with its derivative, by implicit differentiation of eq. (24).

Differentiating `G(ϑ) = RHS(θ, q)` at fixed θ gives `∂ϑ/∂q = (cos³θ / 3) · sin ϑ / cos²ϑ`, and
`F_ω = tan²ϑ/tan²θ` then gives

    ∂F_ω/∂q = (2/3) · cos³θ · sin²ϑ / (cos⁵ϑ · tan²θ)

The two limiting regimes are differentiated in closed form instead, which also keeps the
derivative continuous with the value: `(2/3) exp((2/3)q)` at the pole and
`(2/3)(1-q)^{-5/3}` at the equator.
"""
function elr_flux_factor_and_dq(q::T, θ::T) where {T}
    q <= zero(T) && return one(T), T(2) / T(3)          # dF/dq at q = 0 is 2/3 either way
    # THE SOUTHERN HEMISPHERE. Eq. (24) is written for θ ∈ (0, π/2]; a rigidly rotating star is
    # symmetric about its equator, so `F_ω(π - θ) = F_ω(θ)` and the colatitude folds. Without
    # this the `θ >= π/2 - ε` guard below catches EVERY southern colatitude and hands back the
    # equatorial closed form — measured 0.6 % too hot at θ = 0.2 rad on the wrong side, and
    # invisible to a finite-difference check because the analytic derivative makes the same
    # mistake self-consistently. The test that catches it is `T(θ) == T(π - θ)`.
    θ = θ > T(π) / 2 ? T(π) - θ : θ
    ε = T(1e-4)
    if θ <= ε
        F = exp(T(2) / T(3) * q)
        return F, T(2) / T(3) * F
    elseif θ >= T(π)/2 - ε
        u = max(one(T) - q, eps(T))
        F = u^(-T(2) / T(3))
        return F, T(2) / T(3) * u^(-T(5) / T(3))
    end
    rhs = q * cos(θ)^3 / T(3) + cos(θ) + log(tan(θ / 2))
    ϑ = θ
    for _ in 1:40
        r = _elr_G(ϑ) - rhs
        d = _elr_dG(ϑ)
        d <= eps(T) && break
        ϑn = clamp(ϑ - r / d, T(1e-9), T(π) / 2 - T(1e-9))
        abs(ϑn - ϑ) <= eps(T) * max(one(T), ϑ) && (ϑ = ϑn; break)
        ϑ = ϑn
    end
    tϑ = tan(ϑ); tθ = tan(θ)
    F = (tϑ / tθ)^2
    dF = T(2) / T(3) * cos(θ)^3 * sin(ϑ)^2 / (cos(ϑ)^5 * tθ^2)
    return F, dF
end

"""
    elr_map(rpole, fev, β, tpole, sinθ, cosθ; GM=1) -> Vector

The Espinosa Lara & Rieutord temperature map, normalised to `tpole` at the pole.

    T(θ) = tpole · [ (F_ω(θ)/F_ω(pole)) · (g_eff(θ)/g_pole) ]^β

Same signature as [`vonzeipel_map`](@ref) so the two are interchangeable wherever a rapid
rotator's map is built, and identical to it when `fev → 0`. At β = 1/4 this is ELR eq. (31);
β is free so the two laws can be compared with the same parameter count.
"""
function elr_map(rpole::T, fev::T, β::T, tpole::T,
                 sinθ::AbstractVector{T}, cosθ::AbstractVector{T};
                 GM::T = one(T)) where {T}
    x, = elr_map_and_derivs(rpole, fev, β, tpole, sinθ, cosθ; GM = GM)
    return x
end

"""
    elr_map_and_derivs(rpole, fev, β, tpole, sinθ, cosθ; GM=1)
        -> (x, dx_drpole, dx_dfev, dx_dβ, dx_dtpole)

[`elr_map`](@ref) and its analytic derivatives, in the shape
`vonzeipel_map_and_derivs` returns them, so the gradient and sampling paths need no new
plumbing.

The gravity half is von Zeipel's and is differentiated exactly as there. The new half is the
flux factor, which depends on `fev` only — through `q = (8/27) fev² f(fev sinθ)³` — so

    ∂ln(F_ω(θ)/F_ω(pole))/∂fev = (∂F_θ/∂q · ∂q_θ/∂fev)/F_θ - (∂F_p/∂q · ∂q_p/∂fev)/F_p

and `rpole`, `β` and `tpole` see the factor only through the value they already multiply.
"""
function elr_map_and_derivs(rpole::T, fev::T, β::T, tpole::T,
                            sinθ::AbstractVector{T}, cosθ::AbstractVector{T};
                            GM::T = one(T)) where {T}
    n = length(sinθ)
    x       = Vector{T}(undef, n)
    dx_drp  = Vector{T}(undef, n)
    dx_dfev = Vector{T}(undef, n)
    dx_dβ   = Vector{T}(undef, n)
    dx_dtp  = Vector{T}(undef, n)

    ωc      = sqrt(T(8) * GM / (T(27) * rpole^3))
    ω       = fev * ωc
    dω_dfev = ωc
    dω_drp  = -T(1.5) * ω / rpole
    g_pole  = GM / rpole^2
    dlgpole_drp = -T(2) / rpole
    ω2 = ω * ω

    # The POLE's flux factor, the normalisation. `q` there is `(8/27) fev²` (f(0) = 1).
    qp, dqp = elr_q_and_deriv(fev, zero(T))
    Fp, dFp = elr_flux_factor_and_dq(qp, zero(T))
    dlnFp_dfev = dFp * dqp / Fp

    @inbounds for i in 1:n
        s = sinθ[i]; c = cosθ[i]
        θi = atan(s, c)                      # colatitude, from the precomputed sin/cos
        a = fev * s
        f, fp = f_rapid_rot_and_deriv(a)
        rt = rpole * f
        drt_drp  = f
        drt_dfev = rpole * fp * s

        g_r  = -GM / (rt * rt) + rt * (ω * s)^2
        g_t  = ω2 * rt * s * c
        gθ2  = g_r * g_r + g_t * g_t
        gθ   = sqrt(gθ2)
        Rg   = gθ / g_pole

        qi, dqi = elr_q_and_deriv(fev, s)
        Fi, dFi = elr_flux_factor_and_dq(qi, θi)
        Rf = Fi / Fp                          # the ELR correction, 1 at the pole by design

        R    = Rf * Rg
        xi   = tpole * R^β
        x[i] = xi

        dx_dtp[i] = xi / tpole
        dx_dβ[i]  = xi * log(R)

        # The gravity half, exactly as in `vonzeipel_map_and_derivs`.
        dgr_drp = T(2)*GM/(rt^3)*drt_drp + s*s*(ω2*drt_drp + T(2)*rt*ω*dω_drp)
        dgt_drp = (T(2)*ω*dω_drp*rt + ω2*drt_drp) * s * c
        dlgθ_drp = (g_r*dgr_drp + g_t*dgt_drp) / gθ2
        dx_drp[i] = xi * β * (dlgθ_drp - dlgpole_drp)

        dgr_dfev = T(2)*GM/(rt^3)*drt_dfev + s*s*(ω2*drt_dfev + T(2)*rt*ω*dω_dfev)
        dgt_dfev = (T(2)*ω*dω_dfev*rt + ω2*drt_dfev) * s * c
        dlgθ_dfev = (g_r*dgr_dfev + g_t*dgt_dfev) / gθ2
        # ...plus the flux factor, which is what makes this law different from von Zeipel.
        dlnRf_dfev = dFi * dqi / Fi - dlnFp_dfev
        dx_dfev[i] = xi * β * (dlgθ_dfev + dlnRf_dfev)
    end
    return x, dx_drp, dx_dfev, dx_dβ, dx_dtp
end

"""
    elr_temperature_ratio(ω) -> T_eq / T_pole

ELR eq. (32): the equator-to-pole effective-temperature ratio at β = 1/4, in closed form.

    T_e/T_p = √(2/(2+ω²)) · (1-ω²)^{1/12} · exp(-(4/3) ω²/(2+ω²)³)

Kept because it is an INDEPENDENT check on [`elr_map`](@ref): it follows from eqs. (27), (28)
and `r̃_p = 2/(2+ω²)` without touching the Newton solve or the map assembly, so agreement
between the two tests the whole chain.
"""
elr_temperature_ratio(ω::T) where {T} =
    sqrt(T(2) / (T(2) + ω^2)) * (one(T) - ω^2)^(one(T)/T(12)) *
    exp(-T(4) / T(3) * ω^2 / (T(2) + ω^2)^3)

# =========================================================================================
# The law registry: which gravity-darkening law a model uses
# =========================================================================================
# TWO laws, and the choice is physical rather than a matter of taste. von Zeipel is derived
# for a barotropic star and is the right thing for a SLOW rotator, where the assumption costs
# nothing and the law is one line. ELR drops barotropy and is the right thing for a FAST one,
# where von Zeipel is known to overestimate the pole-to-equator contrast.
#
# β IS FREE IN BOTH. ELR eq. (31) pins it at 1/4 by assuming a grey radiative atmosphere, and
# ROTIR keeps it free anyway: a real atmosphere, a convective envelope, or a limb-darkening
# law absorbing part of the latitudinal profile all move it, and the two laws then have the
# SAME parameter count, which is what makes an evidence ratio between them mean something.
# Pinning it at 1/4 to recover the published law exactly is the ordinary free/fixed control —
# leave `beta` out of `free`, or set its state to "fixed" in the GUI.
#
# WHY THERE IS NO THIRD ENTRY. Two more laws are in `papers/` and neither is a drop-in:
#
#   * Espinosa Lara & Rieutord (2012), A&A 547, A32 — doi:10.1051/0004-6361/201219942 — the
#     same construction for a BINARY, where `g_eff` comes from the Roche potential of two
#     bodies. It would attach to `surface_type = 3`, not here.
#   * Zorec et al. (2017), A&A 606, A32 — doi:10.1051/0004-6361/201730818 — ELR generalised to
#     surface DIFFERENTIAL rotation `Ω(θ) = Ω_0[1 + α Υ(θ)]`. It replaces the surface geometry
#     as well as the flux profile: a differentially rotating surface has no rotational
#     potential, so the shape follows from `g_eff·ds = 0` rather than from ROTIR's Roche
#     factor `f_rapid_rot_and_deriv`, and their eq. (29) for ϑ carries an integral
#     `∫_θ^ϑ sin⁻¹x cos²x (dΥ/dx) dx` inside the implicit equation. Adding it means a new
#     radius solve and a new derivative chain, not a new map.
#     Note it is also what `B_rot` in the `surface_type = 2` schema is reserved for.

"""
    GravityLawSpec

One gravity-darkening law: what to call it, what to show, and where it comes from.

`code` is the integer the schema and the GUI carry it as, the way `ldtype` is carried.
"""
struct GravityLawSpec
    name::Symbol
    code::Int
    label::String
    short::String
    reference::String
    doc::String
end

"""
The implemented gravity-darkening laws, in `code` order. See [`gravity_law_spec`](@ref).
"""
const GRAVITY_LAWS = (
    GravityLawSpec(:vonzeipel, 1, "von Zeipel", "vZ",
        "von Zeipel (1924), MNRAS 84, 665, doi:10.1093/mnras/84.9.665",
        "T_eff = tpole·(g_eff/g_pole)^β. Assumes barotropy, which is incompatible with " *
        "radiative equilibrium in a rotating star, so it is a slow-rotation law; it " *
        "overestimates the pole-to-equator contrast for a fast rotator."),
    GravityLawSpec(:elr, 2, "Espinosa Lara-Rieutord", "ELR",
        "Espinosa Lara & Rieutord (2011), A&A 533, A43, doi:10.1051/0004-6361/201117252",
        "T_eff = tpole·[(F_ω(θ)/F_ω(0))·(g_eff/g_pole)]^β, with the latitudinal flux " *
        "factor F_ω = tan²ϑ/tan²θ from their eq. (24). Assumes the flux is anti-parallel " *
        "to effective gravity instead of barotropy, which holds to better than half a " *
        "degree even near break-up. Reduces to von Zeipel as the rotation goes to zero."),
)

"""
    gravity_law_spec(x) -> GravityLawSpec

The law named by `x`, which may be its `Symbol` (`:vonzeipel`, `:elr`), its string, its
integer `code`, a `GravityLawSpec` itself, or a `star_params`-like object carrying a
`gravity_law` field. An object without that field is von Zeipel, which is what every model
written before the laws became selectable means.
"""
gravity_law_spec(s::GravityLawSpec) = s
function gravity_law_spec(name::Symbol)
    i = findfirst(l -> l.name === name, GRAVITY_LAWS)
    i === nothing && throw(ArgumentError(
        "unknown gravity_law `:$(name)`; implemented: " *
        join((":$(l.name)" for l in GRAVITY_LAWS), ", ")))
    return GRAVITY_LAWS[i]
end
gravity_law_spec(name::AbstractString) = gravity_law_spec(Symbol(name))
function gravity_law_spec(code::Integer)
    i = findfirst(l -> l.code == code, GRAVITY_LAWS)
    i === nothing && throw(ArgumentError(
        "unknown gravity_law code $(code); implemented: " *
        join(("$(l.code) = $(l.name)" for l in GRAVITY_LAWS), ", ")))
    return GRAVITY_LAWS[i]
end
gravity_law_spec(p) = gravity_law_spec(hasproperty(p, :gravity_law) ? p.gravity_law :
                                       :vonzeipel)

"`gravity_law_spec(x).code`, the integer form the schema and the GUI carry."
gravity_law_code(x) = gravity_law_spec(x).code
"`gravity_law_spec(x).name`, the symbol the fit functions take."
gravity_law_name(x) = gravity_law_spec(x).name
"The `code => label` pairs a form needs to offer the laws."
gravity_law_choices() = [l.code => l.label for l in GRAVITY_LAWS]

"""
    gravity_map(law, rpole, fev, β, tpole, sinθ, cosθ; GM=1) -> Vector
    gravity_map_and_derivs(law, rpole, fev, β, tpole, sinθ, cosθ; GM=1)

The rapid rotator's temperature map under `law`, and the same with its analytic derivatives
with respect to `(rpole, fev, β, tpole)`.

`law` is a `Val` of the law's name, so the choice is a type parameter and a closure that
captures it stays type-stable: resolve it once with `Val(gravity_law_name(x))` outside the
hot loop and pass it in. [`vonzeipel_map`](@ref) and [`elr_map`](@ref) are the two
implementations and are interchangeable — same signature, same return shape.
"""
gravity_map(::Val{:vonzeipel}, args...; kwargs...) = vonzeipel_map(args...; kwargs...)
gravity_map(::Val{:elr}, args...; kwargs...) = elr_map(args...; kwargs...)
gravity_map_and_derivs(::Val{:vonzeipel}, args...; kwargs...) =
    vonzeipel_map_and_derivs(args...; kwargs...)
gravity_map_and_derivs(::Val{:elr}, args...; kwargs...) =
    elr_map_and_derivs(args...; kwargs...)

# The Zygote primitive for the ELR map, identical in shape to `vonzeipel_map`'s: the forward
# value and the four analytic derivatives come out of one call, and the pullback is four dot
# products. `sinθ`/`cosθ` are geometry, not parameters, hence `NoTangent`.
function ChainRulesCore.rrule(::typeof(elr_map), rpole, fev, β, tpole,
                              sinθ, cosθ; GM = one(eltype(sinθ)))
    x, dx_drp, dx_dfev, dx_dβ, dx_dtp =
        elr_map_and_derivs(rpole, fev, β, tpole, sinθ, cosθ; GM = GM)
    function elr_map_pullback(x̄)
        v = unthunk(x̄)
        return (NoTangent(), dot(v, dx_drp), dot(v, dx_dfev), dot(v, dx_dβ),
                dot(v, dx_dtp), NoTangent(), NoTangent())
    end
    return x, elr_map_pullback
end
