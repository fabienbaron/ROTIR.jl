# intensity_provider.jl
# ---------------------------------------------------------------------------
# Where the emergent intensity comes from: I(Teff, logg, μ, λ).
#
# WHY AN INTERFACE. `src/intensity.jl` offers two intensity models, `:linear` (the
# Rayleigh-Jeans proxy) and `:planck`. Neither knows about μ, so limb darkening is
# applied separately as `ldmap` with FITTED coefficients — and a fitted `ld1` absorbs
# the error in whichever gravity-darkening law was assumed, because gravity darkening
# and limb darkening both remove flux from the limb and an interferometer measures only
# their sum. That is how a beta Cas NUTS fit reached a NEGATIVE ld1. A model-atmosphere
# intensity carries the μ dependence itself, so there is nothing left to fit.
#
# But no single atmosphere code spans this repo's targets. Korg.jl is pure Julia and
# excellent, and `interpolate_marcs` stops at Teff = 8000 K: beta Cas (7208 K) fits,
# Vega (~9500 K), Regulus, beta Lyr and Spica (25300/20585 K) do not. PMOIRED has the
# same ceiling — its ATLAS9 fetch defaults to TeffMax = 12000. So the provider is an
# INTERFACE with swappable backends, and the hot stars are served by tabulated Kurucz /
# TLUSTY grids rather than by a synthesis code ROTIR has to drive.
#
# WHO OWNS μ. Exactly one of the provider and `ldmap` may apply limb darkening, or it is
# counted twice — `src/di.jl`'s `setup_di` defaults to `ld = true` on top of a SPECTRUM
# intensity and does precisely that. `owns_mu(provider)` declares which, and
# `check_provider_consistency` refuses the combination that double-counts rather than
# leaving it to be noticed in a fit residual.
#
# WHY THE INTERPOLATION IS HAND-ROLLED. The derivative of a multilinear interpolant is
# piecewise CONSTANT in each coordinate. `di.jl` built its `dRdT` by evaluating
# `Interpolations.gradient` at the grid nodes and then interpolating THAT linearly, which
# is piecewise linear — so its derivative disagreed with its own interpolant inside every
# cell. The 16-corner form below returns the value and all four exact partials from one
# cell lookup, which is both cheaper and actually consistent. It also avoids a dependency
# whose interpolants carry no ChainRules rrule.
# ---------------------------------------------------------------------------

"""
    IntensityProvider

Source of the per-tessel emergent intensity `I(Teff, logg, μ, λ)`.

Implement for a new backend:

| function | meaning |
|---|---|
| `provider_intensity(p, Teff, logg, μ, λ)` | the intensity, elementwise over the vectors |
| `provider_intensity_and_derivs(p, Teff, logg, μ, λ)` | it plus `∂I/∂Teff`, `∂I/∂logg`, `∂I/∂μ`, `∂I/∂λ` |
| `owns_mu(p)` | `true` if `I` already contains limb darkening |
| `needs_logg(p)` | `true` if `logg` is actually read (so a distance is required) |
| `provider_support(p)` | `(Teff, logg, λ)` ranges, for refusing silent extrapolation |

`∂I/∂λ` carries the gradient with respect to the velocity field: a velocity-resolved model
evaluates each tessel at its rest-frame wavelength `λ_obs/(1 + v_los/c)`.
"""
abstract type IntensityProvider end

"""
    owns_mu(provider) -> Bool

Whether the provider's intensity already contains limb darkening (then `ldtype` must be 0).
"""
owns_mu(::IntensityProvider) = false
"""
    needs_logg(provider) -> Bool

Whether the provider reads `logg`, in which case the model needs a distance `d`.
"""
needs_logg(::IntensityProvider) = false
"""
    provider_support(provider) -> (Teff = (lo,hi), logg = (lo,hi), λ = (lo,hi))

The parameter ranges the provider covers.
"""
provider_support(::IntensityProvider) =
    (Teff = (0.0, Inf), logg = (-Inf, Inf), λ = (0.0, Inf))

"""
    check_provider_consistency(provider, star_params) -> Vector{String}

Problems with using `provider` with `star_params`; empty if none.

  * `owns_mu(provider)` with `ldtype ≠ 0`: limb darkening counted twice.
  * `needs_logg(provider)` without [`has_physical_scale`](@ref): `logg` undefined.
"""
function check_provider_consistency(p::IntensityProvider, star_params)
    msgs = String[]
    if owns_mu(p) && hasproperty(star_params, :ldtype) && Int(star_params.ldtype) != 0
        push!(msgs, "$(nameof(typeof(p))) already carries the μ dependence, but " *
                    "ldtype = $(star_params.ldtype) applies a limb-darkening law on top " *
                    "of it — limb darkening would be counted twice. Set ldtype = 0.")
    end
    if needs_logg(p) && !has_physical_scale(star_params)
        push!(msgs, "$(nameof(typeof(p))) is indexed by logg, which needs a distance `d` " *
                    "to derive a mass from (rpole, fev, rotation_period); " *
                    "see src/stellar_physics.jl")
    end
    return msgs
end

# ===========================================================================
# Planck: the existing behaviour, behind the new interface
# ===========================================================================

"""
    PlanckProvider()

Non-dimensional Planck intensity, the same as `intensity(T, :planck, …)`.

No μ or `logg` dependence (`owns_mu` and `needs_logg` are false); limb darkening comes from
the `ldtype` law.
"""
struct PlanckProvider <: IntensityProvider end

function provider_intensity(::PlanckProvider, Teff::AbstractVector{T},
                            logg, μ, λ) where {T}
    I = similar(Teff)
    @inbounds for i in eachindex(Teff)
        I[i], _ = planck_and_dT(Teff[i], T(_λ_at(λ, i)))
    end
    return I
end

function provider_intensity_and_derivs(::PlanckProvider, Teff::AbstractVector{T},
                                       logg, μ, λ) where {T}
    n = length(Teff)
    I  = similar(Teff); dT = similar(Teff)
    dλ = zeros(T, n)
    @inbounds for i in 1:n
        λi = T(_λ_at(λ, i))
        I[i], dT[i] = planck_and_dT(Teff[i], λi)
        # ∂B/∂λ of the non-dimensional 1/(e^{c₂/λT} − 1). Carried because a
        # velocity-resolved model shifts λ per tessel, and a zero here would silently
        # make the velocity field invisible to the gradient.
        x   = T(_PLANCK_C2) / (λi * Teff[i])
        em1 = expm1(x)
        dλ[i] = (em1 + one(T)) * x / (em1 * em1 * λi)
    end
    return I, dT, zeros(T, n), zeros(T, n), dλ
end

"""
    RayleighJeansProvider()

Rayleigh–Jeans limit of [`PlanckProvider`](@ref): `I = λT/c₂`, linear in temperature, with
the same normalisation. Equal to `intensity(T, :linear, …)` up to a constant. `owns_mu` is
false.
"""
struct RayleighJeansProvider <: IntensityProvider end

function provider_intensity(::RayleighJeansProvider, Teff::AbstractVector{T},
                            logg, μ, λ) where {T}
    I = similar(Teff)
    @inbounds for i in eachindex(Teff)
        I[i] = T(_λ_at(λ, i)) * Teff[i] / T(_PLANCK_C2)
    end
    return I
end

function provider_intensity_and_derivs(p::RayleighJeansProvider, Teff::AbstractVector{T},
                                       logg, μ, λ) where {T}
    n = length(Teff)
    I = provider_intensity(p, Teff, logg, μ, λ)
    dT = T[T(_λ_at(λ, i)) / T(_PLANCK_C2) for i in 1:n]
    dλ = T[Teff[i] / T(_PLANCK_C2) for i in 1:n]
    return I, dT, zeros(T, n), zeros(T, n), dλ
end

# λ may be one scalar for the whole block or one value per tessel (the Doppler-shifted
# case). Branching here keeps every backend's inner loop identical.
@inline _λ_at(λ::Number, ::Integer) = λ
@inline _λ_at(λ::AbstractVector, i::Integer) = λ[i]

# ===========================================================================
# A rectilinear 4-D grid, with exact partials
# ===========================================================================

"""
    RectGrid4{T}

Values on a rectilinear `(Teff, logg, μ, λ)` grid, multilinearly interpolated.

Axes are strictly increasing and μ starts above 0; `values` is `(nTeff, nlogg, nμ, nλ)`;
`λ` is in metres.
[`TabulatedProvider`](@ref) wraps one as an [`IntensityProvider`](@ref).
"""
struct RectGrid4{T<:AbstractFloat}
    Teff::Vector{T}
    logg::Vector{T}
    μ::Vector{T}
    λ::Vector{T}
    values::Array{T,4}
    function RectGrid4(Teff::Vector{T}, logg::Vector{T}, μ::Vector{T}, λ::Vector{T},
                       values::Array{T,4}) where {T}
        size(values) == (length(Teff), length(logg), length(μ), length(λ)) ||
            throw(DimensionMismatch(
                "RectGrid4: values is $(size(values)) but the axes are " *
                "$((length(Teff), length(logg), length(μ), length(λ)))"))
        for (nm, ax) in (("Teff", Teff), ("logg", logg), ("μ", μ), ("λ", λ))
            length(ax) >= 2 ||
                throw(ArgumentError("RectGrid4: the $nm axis needs at least 2 nodes"))
            issorted(ax) && allunique(ax) ||
                throw(ArgumentError("RectGrid4: the $nm axis must be sorted and strictly " *
                                    "increasing"))
        end
        # A μ = 0 node is a trap, not a convenience: `src/di.jl:260-262` filled it with
        # zeros because SPECTRUM was not run there, so every interpolation in the last
        # cell blended toward zero instead of toward the true limb intensity. Refuse it.
        first(μ) <= zero(T) &&
            throw(ArgumentError("RectGrid4: the μ axis starts at $(first(μ)); a μ ≤ 0 " *
                                "node makes the outermost cell interpolate toward a " *
                                "value no atmosphere produced. Start at ~1e-3 " *
                                "(TLUSTY's own first ray is 0.001)."))
        new{T}(Teff, logg, μ, λ, values)
    end
end

Base.eltype(::RectGrid4{T}) where {T} = T
Base.size(g::RectGrid4) = size(g.values)
# The per-dimension form too, so `size(g, 4)` means the λ count rather than a MethodError.
Base.size(g::RectGrid4, d::Integer) = size(g.values, d)
Base.ndims(::RectGrid4) = 4

# Cell index and normalised offset, clamped to the grid. Clamping rather than throwing is
# deliberate: a sampler WILL step outside the box, and an exception there kills the whole
# chain. `provider_support` is how a caller checks intent up front, and
# `TabulatedProvider`'s `strict` flag is how a validation run demands refusal instead.
@inline function _cell(ax::Vector{T}, x::T) where {T}
    n = length(ax)
    i = searchsortedlast(ax, x)
    i = i < 1 ? 1 : (i > n - 1 ? n - 1 : i)
    h = ax[i+1] - ax[i]
    t = (x - ax[i]) / h
    return i, t, one(T) / h
end

"""
    interp4_and_grad(g, T, lg, m, l) -> (I, ∂I/∂T, ∂I/∂lg, ∂I/∂m, ∂I/∂l)

Multilinear interpolation on `g` with all four exact partials, from one cell lookup.

The value is `Σ₁₆ wᵢ vᵢ` over the cell's corners with `wᵢ` a product of four 1-D weights;
each partial replaces that coordinate's weight factor by `∓1/h`. Because this
differentiates the interpolant that is actually evaluated, the derivative is consistent
with the value everywhere — unlike interpolating a table of nodal derivatives, which is
what `di.jl` did and which disagrees with its own interpolant inside every cell.
"""
@inline function interp4_and_grad(g::RectGrid4{T}, Tk::T, lg::T, m::T, l::T) where {T}
    iT, tT, hT = _cell(g.Teff, Tk)
    ig, tg, hg = _cell(g.logg, lg)
    im, tm, hm = _cell(g.μ,    m)
    il, tl, hl = _cell(g.λ,    l)

    v = g.values
    I = zero(T); dT = zero(T); dg = zero(T); dm = zero(T); dl = zero(T)
    @inbounds for a in 0:1, b in 0:1, c in 0:1, d in 0:1
        val = v[iT+a, ig+b, im+c, il+d]
        # 1-D weight per axis, and the sign its derivative carries.
        wT = a == 0 ? one(T) - tT : tT;  sT = a == 0 ? -hT : hT
        wg = b == 0 ? one(T) - tg : tg;  sg = b == 0 ? -hg : hg
        wm = c == 0 ? one(T) - tm : tm;  sm = c == 0 ? -hm : hm
        wl = d == 0 ? one(T) - tl : tl;  sl = d == 0 ? -hl : hl
        I  += val *  wT * wg * wm * wl
        dT += val *  sT * wg * wm * wl
        dg += val *  wT * sg * wm * wl
        dm += val *  wT * wg * sm * wl
        dl += val *  wT * wg * wm * sl
    end
    return I, dT, dg, dm, dl
end

"""
    TabulatedProvider(grid; strict = false, name = "tabulated")

A precomputed `(Teff, logg, μ, λ)` intensity grid as an [`IntensityProvider`](@ref).

Used for Kurucz ATLAS9, Korg/MARCS and any grid saved with [`save_intensity_grid`](@ref).
`owns_mu` and `needs_logg` are true, so it requires `ldtype = 0` and a distance.

Out-of-range `(Teff, logg, λ)` are extrapolated linearly from the edge cell; `strict = true`
raises an error instead.
"""
struct TabulatedProvider{T} <: IntensityProvider
    grid::RectGrid4{T}
    strict::Bool
    name::String
end
TabulatedProvider(g::RectGrid4; strict::Bool = false, name::AbstractString = "tabulated") =
    TabulatedProvider(g, strict, String(name))

owns_mu(::TabulatedProvider)   = true
needs_logg(::TabulatedProvider) = true
provider_support(p::TabulatedProvider) =
    (Teff = (first(p.grid.Teff), last(p.grid.Teff)),
     logg = (first(p.grid.logg), last(p.grid.logg)),
     λ    = (first(p.grid.λ),    last(p.grid.λ)))

function _assert_in_support(p::TabulatedProvider{T}, Teff, logg, λ) where {T}
    p.strict || return nothing
    s = provider_support(p)
    for (nm, v, (lo, hi)) in (("Teff", Teff, s.Teff), ("logg", logg, s.logg))
        mn, mx = extrema(v)
        (mn < lo || mx > hi) && throw(ArgumentError(
            "$(p.name): $nm spans [$mn, $mx], outside the grid's [$lo, $hi]. " *
            "Extrapolating a model atmosphere is not meaningful; widen the grid or " *
            "constrain the fit."))
    end
    λmn, λmx = λ isa Number ? (λ, λ) : extrema(λ)
    (λmn < s.λ[1] || λmx > s.λ[2]) && throw(ArgumentError(
        "$(p.name): λ spans [$λmn, $λmx], outside the grid's [$(s.λ[1]), $(s.λ[2])]"))
    return nothing
end

# THE OUTPUT ELTYPE FOLLOWS THE INPUTS, NOT THE GRID. ROTIR meshes default to Float32 while
# `load_intensity_grid` defaults to Float64, so returning the grid's eltype would silently
# promote every Float32 model to Float64 on the most common combination there is — and then
# the whole downstream visibility computation with it. `intensity` in src/intensity.jl sets
# the precedent by returning `similar(x)`, and `convert_params` narrows parameters to the mesh
# type for the same reason.
#
# The interpolation itself still runs in the GRID's precision (the query is widened, which is
# lossless when the mesh is narrower), so only the final result is narrowed. Nothing is lost
# beyond that one conversion.
_out_eltype(Teff, logg, μ) = float(promote_type(eltype(Teff), eltype(logg), eltype(μ)))

function provider_intensity(p::TabulatedProvider{T}, Teff::AbstractVector,
                            logg::AbstractVector, μ::AbstractVector, λ) where {T}
    _assert_in_support(p, Teff, logg, λ)
    S = _out_eltype(Teff, logg, μ)
    I = Vector{S}(undef, length(Teff))
    @inbounds for i in eachindex(Teff)
        v, _, _, _, _ = interp4_and_grad(p.grid, T(Teff[i]), T(logg[i]), T(μ[i]),
                                         T(_λ_at(λ, i)))
        I[i] = S(v)
    end
    return I
end

function provider_intensity_and_derivs(p::TabulatedProvider{T}, Teff::AbstractVector,
                                       logg::AbstractVector, μ::AbstractVector,
                                       λ) where {T}
    _assert_in_support(p, Teff, logg, λ)
    n = length(Teff); S = _out_eltype(Teff, logg, μ)
    I = Vector{S}(undef, n); dT = Vector{S}(undef, n); dg = Vector{S}(undef, n)
    dm = Vector{S}(undef, n); dl = Vector{S}(undef, n)
    @inbounds for i in 1:n
        v, a, b, c, d = interp4_and_grad(p.grid, T(Teff[i]), T(logg[i]), T(μ[i]),
                                         T(_λ_at(λ, i)))
        I[i] = S(v); dT[i] = S(a); dg[i] = S(b); dm[i] = S(c); dl[i] = S(d)
    end
    return I, dT, dg, dm, dl
end

# ===========================================================================
# The Zygote seam
# ===========================================================================
# Every partial is DIAGONAL — output element i depends only on input element i — so the
# pullback is four elementwise products, exactly like `intensity`'s rrule in
# src/intensity.jl. The difference is that there are now four differentiable inputs
# instead of one, which is what carries the distance (through logg), the geometry
# (through μ) and the velocity field (through λ) into the gradient.

"""
    provider_map(provider, Teff, logg, μ, λ) -> I

Zygote primitive for a provider's intensity. Differentiable in `Teff`, `logg`, `μ` and
`λ`; `provider` itself is a constant.
"""
provider_map(p::IntensityProvider, Teff, logg, μ, λ) =
    provider_intensity(p, Teff, logg, μ, λ)

function ChainRulesCore.rrule(::typeof(provider_map), p::IntensityProvider,
                              Teff, logg, μ, λ)
    I, dT, dg, dm, dl = provider_intensity_and_derivs(p, Teff, logg, μ, λ)
    function provider_map_pullback(Ībar)
        Ī = unthunk(Ībar)
        # λ may be a scalar shared by every tessel, in which case its cotangent is the
        # sum rather than an elementwise product.
        λ̄ = λ isa Number ? dot(Ī, dl) : Ī .* dl
        return (NoTangent(), NoTangent(), Ī .* dT, Ī .* dg, Ī .* dm, λ̄)
    end
    return I, provider_map_pullback
end

# ===========================================================================
# Synthetic grids, for tests
# ===========================================================================

"""
    analytic_test_grid(; kwargs...) -> (RectGrid4, f)

A `RectGrid4` sampled from a smooth closed-form `f(Teff, logg, μ, λ)`, together with `f`.

For testing the interpolation and its partials:

    f = (Teff/6000)³ · 10^(−0.4(logg−4)) · (0.35 + 0.65μ^0.8) · (λ/2.2e-6)^(−1.7)

The interpolant equals `f` at the nodes; between nodes the error is second order in the
cell size.
"""
function analytic_test_grid(; nT = 9, ng = 7, nm = 11, nl = 6,
                            Teff = (4000.0, 9000.0), logg = (2.5, 4.5),
                            μ = (1e-3, 1.0), λ = (1.5e-6, 2.5e-6))
    f(t, g, m, l) = (t/6000)^3 * 10^(-0.4*(g-4)) * (0.35 + 0.65*m^0.8) *
                    (l/2.2e-6)^(-1.7)
    Ta = collect(range(Teff...; length = nT))
    ga = collect(range(logg...; length = ng))
    ma = collect(range(μ...;    length = nm))
    la = collect(range(λ...;    length = nl))
    v  = [f(t, g, m, l) for t in Ta, g in ga, m in ma, l in la]
    return RectGrid4(Ta, ga, ma, la, v), f
end

# ===========================================================================
# Persisting a grid
# ===========================================================================
# FITS, via FITSIO, which is ALREADY a hard dependency. `src/di.jl`'s reintegration plan
# reached for JLD2 for the same job; a 4-D array plus four axis vectors needs nothing that
# a primary HDU and four image extensions cannot hold, and the file is then readable by
# anything in astronomy rather than only by Julia. `src/surface_map_io.jl` and
# `src/surface_geometry_io.jl` are the existing precedents for FITS I/O here.
#
# The λ axis is stored in METRES, matching `band`/`band_of` and the rest of ROTIR. A
# builder working in Ångström (Korg does) converts on the way in, once, rather than leaving
# two conventions live in the same file.

const _GRID_FITS_VERSION = 1

"""
    save_intensity_grid(path, grid; comment = "") -> path

Write a [`RectGrid4`](@ref) to a FITS file: the 4-D intensity as the primary HDU, the four
axes as image extensions `TEFF`, `LOGG`, `MU` and `LAMBDA`. Reading it back needs no Korg.
"""
function save_intensity_grid(path::AbstractString, g::RectGrid4; comment::AbstractString = "")
    hdr = FITSIO.FITSHeader(
        ["GRIDVER", "NTEFF", "NLOGG", "NMU", "NLAMBDA", "LAMUNIT", "COMMENT"],
        Any[_GRID_FITS_VERSION, length(g.Teff), length(g.logg), length(g.μ), length(g.λ),
            "m", String(comment)],
        ["intensity-grid format version",
         "Teff axis length", "logg axis length", "mu axis length", "lambda axis length",
         "lambda axis unit (metres)", "provenance"])
    FITSIO.FITS(path, "w") do f
        FITSIO.write(f, g.values; header = hdr)
        for (k, (nm, ax)) in enumerate((("TEFF", g.Teff), ("LOGG", g.logg),
                                        ("MU", g.μ), ("LAMBDA", g.λ)))
            FITSIO.write(f, collect(ax))
            # EXTNAME has to be set with `write_key` AFTER the HDU exists. Passing it inside
            # a `FITSHeader` to `write` is silently DROPPED — it is a structural keyword that
            # CFITSIO owns — and the file then reads back with no EXTNAME at all, so every
            # axis lookup fails. `f[k+1]` is the HDU just written (1 is the primary).
            FITSIO.write_key(f[k+1], "EXTNAME", nm)
        end
    end
    return path
end

"""
    load_intensity_grid(path; T = Float64) -> RectGrid4{T}

Read a grid written by [`save_intensity_grid`](@ref). Axes are found by `EXTNAME`.
"""
function load_intensity_grid(path::AbstractString; T::Type = Float64)
    FITSIO.FITS(path, "r") do f
        v = T.(read(f[1]))
        ax = Dict{String,Vector{T}}()
        for i in 2:length(f)
            h = FITSIO.read_header(f[i])
            nm = haskey(h, "EXTNAME") ? String(h["EXTNAME"]) : ""
            nm in ("TEFF", "LOGG", "MU", "LAMBDA") && (ax[nm] = T.(vec(FITSIO.read(f[i]))))
        end
        for nm in ("TEFF", "LOGG", "MU", "LAMBDA")
            haskey(ax, nm) || error("load_intensity_grid: $path has no $nm extension")
        end
        # The primary header records the axis lengths independently, so a file whose
        # extensions were reordered or truncated fails here rather than producing a grid
        # whose axes are silently attached to the wrong dimensions.
        pri = FITSIO.read_header(f[1])
        for (kw, nm) in (("NTEFF", "TEFF"), ("NLOGG", "LOGG"), ("NMU", "MU"),
                         ("NLAMBDA", "LAMBDA"))
            haskey(pri, kw) || continue
            pri[kw] == length(ax[nm]) || error(
                "load_intensity_grid: $path declares $kw = $(pri[kw]) but its $nm " *
                "extension has $(length(ax[nm])) entries")
        end
        return RectGrid4(ax["TEFF"], ax["LOGG"], ax["MU"], ax["LAMBDA"], v)
    end
end

# ===========================================================================
# Korg, as a grid BUILDER
# ===========================================================================
# Declared here, implemented in ext/ROTIRKorgExt.jl — the same stub-plus-extension pattern
# as `fit_parametric` (src/bootstrap.jl:320 / ext/ROTIRZygoteExt.jl). Korg pulls in a
# MARCS atmosphere artifact and a large linelist, so it is a weak dependency: `using Korg`
# activates this.
#
# Korg is NOT a new provider type. It produces a `RectGrid4`, which `TabulatedProvider`
# already evaluates — so there is exactly one evaluation path, one rrule and one set of
# tests, and the synthesis code is needed only when a grid is first built.

"""
    build_korg_grid(; Teff, logg, μ, λ, kwargs...) -> RectGrid4

Synthesize a `(Teff, logg, μ, λ)` specific-intensity grid with Korg.jl.

Requires `using Korg`. MARCS covers Teff 2800–8000 K and logg −0.5 to 5.5; hotter stars
need a Kurucz grid ([`read_kurucz_intensity`](@ref)) or a grid loaded with
[`load_intensity_grid`](@ref).
"""
function build_korg_grid end

"""
    korg_provider(; kwargs...) -> TabulatedProvider

[`build_korg_grid`](@ref) wrapped as a [`TabulatedProvider`](@ref); `strict` is passed to it
and the other keywords to `build_korg_grid`. Requires `using Korg`.
"""
function korg_provider end

# ===========================================================================
# One provider from a form's worth of settings
# ===========================================================================
# The four backends above are constructed in four different ways — a bare struct, a FITS read,
# an ATLAS9 parse, a synthesis run — and a front-end that offers the CHOICE needs one function
# that turns the choice plus its settings into a provider. The settings themselves, their
# labels, units, defaults and which backend reads which, are declared once in
# `INTENSITY_PARAMS` (src/surface_schema.jl); this is the consumer of that declaration.
#
# It belongs here rather than in the GUI because nothing about it is a GUI concern: a demo
# script that wants "the H-band Kurucz grid for this pack" asks for it the same way, and a
# second copy of the mapping inside the extension would be free to drift from the schema.

"""
    build_intensity_provider(cfg, path) -> IntensityProvider

Build the provider selected by `cfg[:intensity_source]` (codes in
[`INTENSITY_PARAMS`](@ref): 1 Planck, 2 grid file, 3 Kurucz pack, 4 Korg, 5 Rayleigh–Jeans).

`cfg` maps the numeric fields of `INTENSITY_PARAMS` to values (a `Dict` or a NamedTuple).
`path` is `grid_path`, read only for a custom grid (`grid_choice = 0`). Raises on a bad
configuration; [`validate_intensity_source`](@ref) checks one without building it.

Cost: 1 and 5 are immediate, 2 reads a FITS file, 3 parses the 67.5 MB Kurucz pack
(downloaded on first use), 4 runs `grid_nodes²` Korg syntheses (seconds to minutes). Nothing
is cached; save a built grid with [`save_intensity_grid`](@ref) and load it as source 2.
"""
function build_intensity_provider(cfg, path::AbstractString = "")
    src = round(Int, float(_cfg(cfg, :intensity_source, INTENSITY_PARAMS[1].default)))
    λlo, λhi = intensity_lambda_range(cfg)
    src == 1 && return PlanckProvider()
    src == 5 && return RayleighJeansProvider()
    if src == 2
        # A SHIPPED grid is found by name inside the package; only a custom one uses `path`.
        f = intensity_grid_path(cfg, path)
        isfile(f) || throw(ArgumentError("build_intensity_provider: no such grid file: $f"))
        return TabulatedProvider(load_intensity_grid(f); name = basename(f))
    end
    if src == 3
        # THE PACK IS FETCHED, NOT ASKED FOR: one file serves every star, so its location is
        # ROTIR's business — see `kurucz_pack_path`.
        pack = kurucz_pack_path(; download = true)
        g = read_kurucz_intensity(pack; λrange = (λlo, λhi),
                                  Teff = (float(_cfg(cfg, :grid_teff_min, 3500.0)),
                                          float(_cfg(cfg, :grid_teff_max, 10000.0))),
                                  logg = (float(_cfg(cfg, :grid_logg_min, 2.0)),
                                          float(_cfg(cfg, :grid_logg_max, 4.5))))
        return TabulatedProvider(g; name = "kurucz:" * basename(pack))
    end
    src == 4 || throw(ArgumentError("build_intensity_provider: unknown intensity_source $src"))
    # THE WEAK DEPENDENCY, reported as a missing package rather than as a MethodError on a
    # stub. `build_korg_grid` is declared with no methods (see above), so without `using Korg`
    # the call fails with "no method matching build_korg_grid(; …)" — which reads as a wrong
    # call rather than as a package that is not loaded.
    isempty(methods(build_korg_grid)) &&
        throw(ArgumentError("Korg synthesis needs `using Korg`: it is a weak dependency, " *
                            "so ROTIR carries only the stub until it is loaded"))
    n = round(Int, float(_cfg(cfg, :grid_nodes, 6.0)))
    g = build_korg_grid(; Teff = (float(_cfg(cfg, :grid_teff_min, 3500.0)),
                                  float(_cfg(cfg, :grid_teff_max, 7500.0)), n),
                          logg = (float(_cfg(cfg, :grid_logg_min, 2.0)),
                                  float(_cfg(cfg, :grid_logg_max, 4.5)), n),
                          # Korg takes Angstrom; `intensity_lambda_range` gave metres.
                          λ = (λlo * 1e10, λhi * 1e10))
    return TabulatedProvider(g; name = "korg")
end

# ===========================================================================
# The Kurucz ATLAS9 pack: one file for every star, so ROTIR keeps it
# ===========================================================================
const KURUCZ_PACK_URL = "http://kurucz.harvard.edu/grids/gridp00/ip00k2.pck19"
const KURUCZ_PACK_ENV = "ROTIR_KURUCZ_PACK"

"""
    intensity_cache_dir() -> String

Per-user cache for downloaded intensity data: `scratchspaces/ROTIR/intensity` under the first
Julia depot.
"""
intensity_cache_dir() = joinpath(first(DEPOT_PATH), "scratchspaces", "ROTIR", "intensity")

"""
    kurucz_pack_path(; download = false) -> String

Path of the solar-metallicity ATLAS9 intensity pack `ip00k2.pck19` (67.5 MB, Teff
3500–50000 K): `ENV["$(KURUCZ_PACK_ENV)"]` if set, otherwise in [`intensity_cache_dir`](@ref).
With `download = true` a missing pack is downloaded from $(KURUCZ_PACK_URL); otherwise the
path is returned whether or not the file exists.
"""
function kurucz_pack_path(; download::Bool = false)
    forced = get(ENV, KURUCZ_PACK_ENV, "")
    isempty(forced) || return abspath(expanduser(forced))
    dest = joinpath(intensity_cache_dir(), "ip00k2.pck19")
    (isfile(dest) || !download) && return dest
    mkpath(dirname(dest))
    tmp = dest * ".part"
    @info "Downloading the Kurucz ATLAS9 intensity pack (67.5 MB) once" url = KURUCZ_PACK_URL dest
    # Downloads is a stdlib that ROTIR does not otherwise depend on; reach it the way Base's own
    # deprecated `download` does, by package id, rather than adding a dependency every
    # environment that develops ROTIR would have to re-resolve for.
    D = Base.require(Base.PkgId(Base.UUID("f43a241f-c20a-4ad4-852c-f6b1247861c6"), "Downloads"))
    Base.invokelatest(D.download, KURUCZ_PACK_URL, tmp)
    mv(tmp, dest; force = true)                 # atomic: a half-written pack is never "found"
    return dest
end

# A `Dict`, a NamedTuple or anything else `get`/`getproperty` reaches — the GUI holds its
# settings in a `Dict{Symbol,Float64}` and a script would pass a NamedTuple.
_cfg(cfg::AbstractDict, k::Symbol, d) = get(cfg, k, d)
_cfg(cfg, k::Symbol, d) = hasproperty(cfg, k) ? getproperty(cfg, k) : d
