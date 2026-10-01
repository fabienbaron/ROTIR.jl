# stellar_physics.jl
# ---------------------------------------------------------------------------
# Physical units for the rapid rotator: absolute mass, radius and surface gravity.
#
# WHY THIS FILE EXISTS. ROTIR's geometry is angular throughout — `rpole` is in mas —
# and its gravity-darkening maps only ever use the RATIO `g_θ/g_pole`, so they run
# with `GM = 1` and never need a mass (see `vonzeipel_map_and_derivs`). That is
# correct and it is also a dead end: a model atmosphere is indexed by (Teff, logg),
# and `logg` is an ABSOLUTE gravity in cgs. Without it the emergent intensity cannot
# be predicted, so limb darkening has to be FITTED — and `ld1` then absorbs the error
# in whichever gravity-darkening law was assumed. That degeneracy (corr(beta, ld1)
# = −0.63 at fev ≈ 0.8, measured on six lambda And epochs) is what drove a beta Cas
# NUTS fit to a negative limb-darkening coefficient. See src/gravity_darkening.jl and
# demos/gravity_law_comparison.jl.
#
# WHAT IS DERIVED, AND WHY NOT MASS AS A FREE PARAMETER. For a Roche rotator
#
#     Ω = fev · Ω_crit,    Ω_crit = √(8GM / 27 R_p³),    Ω = 2π / P_rot
#
# so (rpole, d, fev, mass, P_rot) are FIVE quantities with ONE relation — four free.
# ROTIR already carries `rpole`, `fev` and `rotation_period`, so adding the distance
# alone closes the system and mass comes out DERIVED:
#
#     M = 27 R_p³ Ω² / (8 G fev²)
#
# This costs one new parameter instead of two, and `d` has a tight Gaia prior, so the
# posterior barely widens. It also repairs an existing inconsistency for free:
# `rotate_star` spins the surface at 2πt/rotation_period while `update_radii_rapidrot`
# and the temperature maps take the rotation rate from `frac_escapevel`, and NOTHING
# tied the two. They are not independent for a Roche surface; now they are tied.
#
# Checked against beta Cas (rpole = 0.849 mas, d = 16.8 pc, P = 0.893 d, fev = 0.92):
# R_p = 3.07 R☉, M = 1.93 M☉, logg_pole = 3.75 — all three match an F2III-IV star.
#
# WHY THE GRAVITY IS RECOMPUTED HERE rather than returned from the map functions.
# `vonzeipel_map_and_derivs` computes `g_θ` internally but discards everything except
# `g_θ/g_pole`, and it does so with `GM = 1` and `rpole` in MAS. Its `dx_drpole` is
# therefore a derivative with respect to an ANGULAR radius. Feeding it physical units
# to extract an absolute gravity would silently change what that derivative means and
# break every rrule downstream of it. The gravity kernel below is the same algebra in
# SI, kept deliberately separate, and `roche_gravity_and_derivs` is the single place
# it lives so the two cannot drift apart in form.
# ---------------------------------------------------------------------------

# IAU 2015 Resolution B3 nominal values, plus the IAU definition of the astronomical
# unit. These are DEFINED constants, not measurements, so they carry no uncertainty.
const _AU_M      = 1.495978707e11      # m, exact by IAU definition
const _RSUN_M    = 6.957e8             # m, IAU 2015 nominal solar radius
const _GMSUN_SI  = 1.3271244e20        # m³/s², IAU 2015 nominal solar mass parameter
const _DAY_S     = 86400.0             # s
const _LN10      = log(10.0)

# One mas at one parsec is 10⁻³ AU — that is what a parsec IS. So a polar radius in
# mas times a distance in pc is a length, with this single factor.
const _MASPC_TO_M = 1e-3 * _AU_M

"Polar radius in metres, from `rpole` in mas and `d` in pc."
@inline polar_radius_m(rpole_mas::T, d_pc::T) where {T} =
    rpole_mas * d_pc * T(_MASPC_TO_M)

"Polar radius in solar radii — the number a reader recognises."
@inline polar_radius_rsun(rpole_mas::T, d_pc::T) where {T} =
    polar_radius_m(rpole_mas, d_pc) / T(_RSUN_M)

"Angular rotation rate in rad/s, from a period in days."
@inline angular_rate(period_days::T) where {T} = T(2π) / (period_days * T(_DAY_S))

"""
    derive_mass(rpole_mas, d_pc, fev, period_days) -> M / M☉

The mass a Roche rotator must have for its shape (`fev`) and its spin (`period_days`)
to be consistent at the given polar radius.

    Ω = 2π/P,    M = 27 R_p³ Ω² / (8 G fev²)

from `Ω = fev·√(8GM/27R_p³)`, i.e. ELR/Roche's definition of the critical rate. The
result is in SOLAR masses, obtained by dividing by the nominal `GM☉` rather than by
multiplying out `G` and `M☉` separately, which avoids carrying `G`'s uncertainty.

Returns `Inf` as `fev → 0`: a spherical star spinning at a finite rate needs infinite
mass to make that rate a vanishing fraction of critical. Callers that allow `fev = 0`
must not ask for a mass.
"""
@inline function derive_mass(rpole_mas::T, d_pc::T, fev::T, period_days::T) where {T}
    Rp = polar_radius_m(rpole_mas, d_pc)
    Ω  = angular_rate(period_days)
    return T(27) * Rp^3 * Ω^2 / (T(8) * T(_GMSUN_SI) * fev^2)
end

"""
    derive_mass_and_dlog(rpole_mas, d_pc, fev, period_days)
        -> (M, ∂lnM/∂rpole, ∂lnM/∂d, ∂lnM/∂fev, ∂lnM/∂P)

[`derive_mass`](@ref) with its LOGARITHMIC derivatives, which are what the chain rule
wants downstream and which are exact one-liners here:

    ln M = const + 3 ln R_p + 2 ln Ω − 2 ln fev,   R_p ∝ rpole·d,   Ω ∝ 1/P

so every derivative is a bare power. No cancellation, nothing to guard.
"""
@inline function derive_mass_and_dlog(rpole_mas::T, d_pc::T, fev::T,
                                      period_days::T) where {T}
    M = derive_mass(rpole_mas, d_pc, fev, period_days)
    return M, T(3)/rpole_mas, T(3)/d_pc, -T(2)/fev, -T(2)/period_days
end

# ---------------------------------------------------------------------------
# The Roche effective gravity, in SI
# ---------------------------------------------------------------------------
# Same algebra as `vonzeipel_map_and_derivs` (src/parametric_gradient.jl:38) and
# `elr_map_and_derivs` (src/gravity_darkening.jl:213), with two differences that matter:
#
#   * GM and r are PHYSICAL (m³/s², m), so `g` comes out in m/s² and can be logged.
#   * Ω is taken straight from the period, NOT rebuilt as fev·√(8GM/27R³). The two are
#     equal by construction once the mass is derived, and using the period directly
#     makes ∂Ω/∂rpole = ∂Ω/∂d = ∂Ω/∂fev = 0 — which removes three terms from every
#     derivative below and one sqrt from the inner loop.

"""
    roche_gravity_and_derivs(Rp, GM, Ω, fev, sinθ, cosθ,
                             dRp_dq, dGM_dq, dΩ_dq) -> (g, ∂g/∂q)

Effective gravity magnitude at one colatitude on a Roche surface, in the same units as
`GM`/`Rp`, and its derivative with respect to one parameter `q`.

    r_θ = R_p · f(fev sinθ)
    g_r = −GM/r_θ² + r_θ (Ω sinθ)²          radial
    g_θ = Ω² r_θ sinθ cosθ                  latitudinal
    g   = √(g_r² + g_θ²)

`q` enters only through the three sensitivities passed in — `dRp_dq`, `dGM_dq`, `dΩ_dq`
— plus, for `q = fev`, the shape factor's own `f'`, which is why `fev` is passed
separately from its sensitivity. `f` is [`f_rapid_rot_and_deriv`](@ref), which is
series-guarded near zero (the textbook form cancels catastrophically there).
"""
@inline function roche_gravity_and_derivs(Rp::T, GM::T, Ω::T, fev::T,
                                          s::T, c::T,
                                          dRp_dq::T, dGM_dq::T, dΩ_dq::T,
                                          dfev_dq::T) where {T}
    f, fp = f_rapid_rot_and_deriv(fev * s)
    rt    = Rp * f
    # r_θ depends on q through R_p and, when q is fev, through f(fev·sinθ).
    drt_dq = dRp_dq * f + Rp * fp * s * dfev_dq

    Ω2   = Ω * Ω
    g_r  = -GM / (rt * rt) + rt * (Ω * s)^2
    g_t  = Ω2 * rt * s * c
    g2   = g_r * g_r + g_t * g_t
    g    = sqrt(g2)

    dgr_dq = T(2)*GM/(rt^3)*drt_dq - dGM_dq/(rt*rt) +
             s*s*(Ω2*drt_dq + T(2)*rt*Ω*dΩ_dq)
    dgt_dq = (T(2)*Ω*dΩ_dq*rt + Ω2*drt_dq) * s * c
    return g, (g_r*dgr_dq + g_t*dgt_dq) / g
end

"""
    logg_map(rpole_mas, d_pc, fev, period_days, sinθ, cosθ) -> Vector

Per-tessel `log₁₀(g_eff)` in **cgs** for a rapid rotator, with the mass derived from
the spin (see [`derive_mass`](@ref)).

This is the quantity a model atmosphere is indexed by, alongside `Teff`. It varies by
0.3–1.0 dex from pole to equator on a fast rotator — which is exactly why a grid with
a single scalar `logg` (as `src/di.jl`'s `modelGrid` has) cannot serve this surface.
"""
function logg_map(rpole_mas::T, d_pc::T, fev::T, period_days::T,
                  sinθ::AbstractVector{T}, cosθ::AbstractVector{T}) where {T}
    x, = logg_map_and_derivs(rpole_mas, d_pc, fev, period_days, sinθ, cosθ)
    return x
end

"""
    logg_map_and_derivs(rpole_mas, d_pc, fev, period_days, sinθ, cosθ)
        -> (logg, dlogg_drpole, dlogg_dd, dlogg_dfev, dlogg_dP)

[`logg_map`](@ref) and its analytic derivatives, in the shape the other
`*_map_and_derivs` functions in this package return them, so the gradient path needs
no new plumbing.

The chain is short because the mass's log-derivatives are bare powers
([`derive_mass_and_dlog`](@ref)) and because Ω comes straight from the period:

    ∂logg/∂q = (1/ln10) · (1/g) · ∂g/∂q

`cosθ` enters only through the latitudinal term `g_θ ∝ sinθ cosθ`, which vanishes at
both the pole and the equator — so `logg` is stationary there, as it must be.
"""
function logg_map_and_derivs(rpole_mas::T, d_pc::T, fev::T, period_days::T,
                             sinθ::AbstractVector{T},
                             cosθ::AbstractVector{T}) where {T}
    n = length(sinθ)
    lg      = Vector{T}(undef, n)
    dlg_drp = Vector{T}(undef, n)
    dlg_dd  = Vector{T}(undef, n)
    dlg_dfev= Vector{T}(undef, n)
    dlg_dP  = Vector{T}(undef, n)

    Rp = polar_radius_m(rpole_mas, d_pc)
    Ω  = angular_rate(period_days)
    M, dlnM_drp, dlnM_dd, dlnM_dfev, dlnM_dP =
        derive_mass_and_dlog(rpole_mas, d_pc, fev, period_days)
    GM = M * T(_GMSUN_SI)

    # R_p ∝ rpole·d, so its sensitivities are just R_p/rpole and R_p/d. Ω depends on
    # the period alone. GM's sensitivities come from the mass's log-derivatives.
    dRp_drp = Rp / rpole_mas;  dRp_dd = Rp / d_pc
    dGM_drp = GM * dlnM_drp;   dGM_dd = GM * dlnM_dd
    dGM_dfev= GM * dlnM_dfev;  dGM_dP = GM * dlnM_dP
    dΩ_dP   = -Ω / period_days
    z = zero(T); o = one(T)
    invln10 = one(T) / T(_LN10)

    @inbounds for i in 1:n
        s = sinθ[i]; c = cosθ[i]
        # One call per parameter. The value `g` is identical across the four; taking it
        # from the first is what keeps this exact rather than four nearly-equal values.
        g, dg_drp  = roche_gravity_and_derivs(Rp, GM, Ω, fev, s, c, dRp_drp, dGM_drp, z, z)
        _, dg_dd   = roche_gravity_and_derivs(Rp, GM, Ω, fev, s, c, dRp_dd,  dGM_dd,  z, z)
        _, dg_dfev = roche_gravity_and_derivs(Rp, GM, Ω, fev, s, c, z,       dGM_dfev,z, o)
        _, dg_dP   = roche_gravity_and_derivs(Rp, GM, Ω, fev, s, c, z,       dGM_dP,  dΩ_dP, z)

        # g is in m/s²; logg is quoted in cgs, hence the ×100.
        lg[i]       = log10(g * T(100))
        k           = invln10 / g
        dlg_drp[i]  = k * dg_drp
        dlg_dd[i]   = k * dg_dd
        dlg_dfev[i] = k * dg_dfev
        dlg_dP[i]   = k * dg_dP
    end
    return lg, dlg_drp, dlg_dd, dlg_dfev, dlg_dP
end

# ---------------------------------------------------------------------------
# Reported quantities
# ---------------------------------------------------------------------------

"Polar surface gravity, log₁₀(cgs). The normalisation the maps divide by."
@inline function logg_pole(rpole_mas::T, d_pc::T, fev::T, period_days::T) where {T}
    Rp = polar_radius_m(rpole_mas, d_pc)
    GM = derive_mass(rpole_mas, d_pc, fev, period_days) * T(_GMSUN_SI)
    return log10(GM / (Rp * Rp) * T(100))
end

"""
    equatorial_velocity(rpole_mas, d_pc, fev, period_days) -> km/s

True equatorial rotation velocity, `Ω R_eq`, with `R_eq = R_p·f(fev)` from the Roche
shape factor. Not `Ω R_p` — the equator of a fast rotator is up to 1.5x the polar
radius, and that factor is the whole point of the Roche surface.
"""
@inline function equatorial_velocity(rpole_mas::T, d_pc::T, fev::T,
                                     period_days::T) where {T}
    Rp = polar_radius_m(rpole_mas, d_pc)
    Req = Rp * f_rapid_rot_and_deriv(fev)[1]
    return angular_rate(period_days) * Req / T(1000)   # m/s → km/s
end

# NOT named `vsini`: that is the single most common variable name in a Doppler-imaging
# script, and Julia refuses assignment to an imported binding ("cannot assign a value to
# imported variable"). Exporting `vsini` would break every script with a `vsini = 25.0`
# line in it. The recognisable name survives as a FIELD of `derived_quantities`.
"Projected equatorial velocity, `v_eq sin i`, in km/s. `inclination` in degrees."
@inline projected_veq(rpole_mas::T, d_pc::T, fev::T, period_days::T,
                      inclination_deg::T) where {T} =
    equatorial_velocity(rpole_mas, d_pc, fev, period_days) *
    sin(inclination_deg * T(π) / T(180))

"True whether `p` carries the distance that makes physical units available."
has_physical_scale(p) = hasproperty(p, :d) && hasproperty(p, :rpole) &&
                        hasproperty(p, :frac_escapevel) &&
                        hasproperty(p, :rotation_period)

"""
    derived_quantities(star_params) -> NamedTuple

Everything physical that `star_params` implies but does not store: `mass` [M☉],
`rpole_rsun`, `req_rsun`, `logg_pole`, `veq` and `vsini` [km/s], `omega` [rad/s].

Meant for reporting alongside a fit — a posterior on `(rpole, d, fev, P)` is a
posterior on the mass too, and this is what turns one into the other. Returns
`nothing` when the parameters carry no distance, rather than inventing one.
"""
function derived_quantities(p)
    has_physical_scale(p) || return nothing
    T = float(typeof(p.rpole))
    rp = T(p.rpole); d = T(p.d); fev = T(p.frac_escapevel)
    P = T(p.rotation_period)
    inc = hasproperty(p, :inclination) ? T(p.inclination) : T(90)
    return (mass       = derive_mass(rp, d, fev, P),
            rpole_rsun = polar_radius_rsun(rp, d),
            req_rsun   = polar_radius_rsun(rp, d) * f_rapid_rot_and_deriv(fev)[1],
            logg_pole  = logg_pole(rp, d, fev, P),
            veq        = equatorial_velocity(rp, d, fev, P),
            vsini      = projected_veq(rp, d, fev, P, inc),
            omega      = angular_rate(P))
end

# ---------------------------------------------------------------------------
# The Zygote seam
# ---------------------------------------------------------------------------
# WHY THIS IS REQUIRED, not an optimisation. `logg_map` fills its output with `setindex!`,
# and Zygote refuses to differentiate through array mutation — `Zygote.gradient` on it fails
# outright with "Mutating arrays is not supported". Without this rrule the distance cannot
# reach a gradient-based fit at all, so `d` could be added to the parameter vector and would
# silently contribute nothing.
#
# Same shape as `vonzeipel_map`'s rrule (src/parametric_gradient.jl:99): one call produces the
# value and all four derivatives, and the pullback is four dot products. `sinθ`/`cosθ` are
# geometry, fixed by the tessellation, hence `NoTangent`.

function ChainRulesCore.rrule(::typeof(logg_map), rpole, d, fev, period_days, sinθ, cosθ)
    x, dx_drp, dx_dd, dx_dfev, dx_dP =
        logg_map_and_derivs(rpole, d, fev, period_days, sinθ, cosθ)
    function logg_map_pullback(x̄raw)
        x̄ = unthunk(x̄raw)
        return (NoTangent(), dot(x̄, dx_drp), dot(x̄, dx_dd), dot(x̄, dx_dfev), dot(x̄, dx_dP),
                NoTangent(), NoTangent())
    end
    return x, logg_map_pullback
end

"""
    limb_mu_vec(nz) -> μ

`max.(nz, 0)` as a named function, for the μ a provider is indexed by.

Written as a broadcast rather than a hand-written primitive because Zygote differentiates
`max` correctly on its own — the subgradient at `nz = 0` is what [`mu_and_dmu`](@ref) already
picks (zero on the back side), so the AD and the forward paths agree without a rule.
"""
limb_mu_vec(nz) = max.(nz, zero(eltype(nz)))

"""
    derived_summary_text(star_params) -> String

One line of the PHYSICAL quantities a `star_params` implies but does not store — mass, polar
and equatorial radius, polar `logg`, equatorial and projected rotation velocity. Empty when
the model carries no distance, or is not a rapid rotator.

Lives in the core rather than in the GUI extension even though the GUI is its only caller:
it is a pure function of the parameters, so putting it here is what makes it testable without
GLMakie and QML loaded. `shell_derived_summary` is the thin session-reading wrapper.

WHY THIS EARNS ITS PLACE IN THE PANEL. The four numbers a user types — `rpole` in mas,
`frac_escapevel`, `rotation_period`, `d` — imply a mass through
M = 27R_p³Ω²/(8G·fev²), and nothing on screen said what it was. A fit that wandered to
rpole 0.98 / fev 0.59 implied a 7.1 Msun star for beta Cas and looked unremarkable in the
form; it would have been obvious here at a glance. `validate_star_params` flags an
implausible mass, but only once it is already absurd.
"""
function derived_summary_text(p)
    Int(get(p, :surface_type, -1)) == 2 || return ""
    q = derived_quantities(p)
    q === nothing && return ""
    return Printf.@sprintf("M = %.3f M⊙   ·   R_p = %.3f, R_eq = %.3f R⊙   ·   logg_p = %.3f   ·   v_eq = %.1f, vsini = %.1f km/s",
                           q.mass, q.rpole_rsun, q.req_rsun, q.logg_pole, q.veq, q.vsini)
end
