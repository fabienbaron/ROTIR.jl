# intensity.jl
# ---------------------------------------------------------------------------
# Temperature → band intensity.
#
# ROTIR's forward model historically used the temperature map *itself* as the
# surface brightness (`model = :linear`). That is the Rayleigh–Jeans limit and it
# is a good approximation only in H and K. As soon as shorter wavelengths are in
# play (CHARA now reaches R, and V is on the horizon) the Wien curvature matters:
# for a Spica-like pair (25300 K / 20585 K) the linear proxy misstates the
# component flux ratio by −3.2 % in K, −4.1 % in H, −5.5 % in J, −11.0 % in R and
# −13.1 % in V. Hence `model = :planck`.
#
# The Planck form here is NON-DIMENSIONAL: the λ⁵ prefactor and every physical
# constant are dropped, because visibilities are flux-normalized and only the
# *shape* of I(T) across the surface matters. This stays valid for a two-component
# binary as long as both components are evaluated at the same λ (they are — one
# `band` per OIFITS wavelength block), since the dropped prefactor then cancels
# from the flux ratio too.
#
# This file is included EARLY in ROTIR.jl so that both the χ² paths
# (oichi2_spheroid.jl, oichi2_binary.jl) and the AD path (parametric_gradient.jl)
# can use it.
# ---------------------------------------------------------------------------

const _PLANCK_C2 = 1.438776877e-2   # hc/k  [m·K]

# Non-dimensional band intensity. The global λ⁵ (and any constant) prefactor is
# dropped because it cancels in the flux-normalized visibilities; only the SHAPE
# of I(T) across the disk matters. dB/dT > 0 (the Wien curvature that breaks the
# tpole degeneracy). λ (`band`) in metres, T in Kelvin.
@inline function planck_and_dT(Tk::T, λ::T) where {T}
    x   = T(_PLANCK_C2) / (λ * Tk)
    em1 = expm1(x)                       # e^x − 1  (stable for small x)
    ex  = em1 + one(T)                   # e^x
    B     = one(T) / em1
    dB_dT = ex * T(_PLANCK_C2) / (em1 * em1 * λ * Tk * Tk)
    return B, dB_dT
end

"""
    intensity(x, model, band) -> I

Map the temperature map `x` (K) to a per-tessel intensity. `model = :linear`
(default, `I = x`, the historical linear proxy) or `:planck` (`band` = wavelength
in metres). Global scale is irrelevant (visibilities are normalized).
"""
function intensity(x::AbstractVector{T}, model::Symbol, band) where {T}
    if model === :linear
        return copy(x)
    elseif model === :planck
        band === nothing && error("intensity: model = :planck requires `band` (wavelength in metres)")
        λ = T(band)
        I = similar(x)
        @inbounds for i in eachindex(x)
            I[i], _ = planck_and_dT(x[i], λ)
        end
        return I
    else
        error("intensity: unknown model $(model) (use :linear or :planck)")
    end
end

function _intensity_dT(x::AbstractVector{T}, model::Symbol, band) where {T}
    if model === :linear
        return ones(T, length(x))
    else
        λ = T(band); dI = similar(x)
        @inbounds for i in eachindex(x)
            _, dI[i] = planck_and_dT(x[i], λ)
        end
        return dI
    end
end

function ChainRulesCore.rrule(::typeof(intensity), x::AbstractVector, model::Symbol, band)
    I = intensity(x, model, band)
    dI = _intensity_dT(x, model, band)
    function intensity_pullback(Ī)
        x̄ = unthunk(Ī) .* dI
        return (NoTangent(), x̄, NoTangent(), NoTangent())
    end
    return I, intensity_pullback
end

"""
    band_of(data) -> λ

Mean effective wavelength (metres) of an `OIdata` block, for use as the `band`
argument of [`intensity`](@ref). `readoifits` already splits the data by
wavelength bin, so within one block λ is nearly constant; this is the
representative value for that block.
"""
band_of(data) = sum(data.uv_lam) / length(data.uv_lam)
