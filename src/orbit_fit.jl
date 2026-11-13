# =======================================================================================
# Generic orbit fitting from OIFITS
# =======================================================================================
# Fit the relative orbit of a two-component system directly to interferometric
# observables, given nothing but `readoifits` output and a starting guess.
#
# The pieces existed already but were scattered: `src/orbits.jl` has the Kepler solver and
# Thiele-Innes projection, OITOOLS supplies the analytic component visibilities, and the
# per-epoch precompute that makes this fast lived inside a bespoke beta Lyr script. This
# file is the join.
#
# WHY IT IS FAST. Two things, both of which matter more than they look:
#
#   * One Kepler solve per distinct TIME, not per uv point. An OIFITS night has thousands
#     of uv samples spread over tens of exposures, so the separation only needs computing
#     at the unique timestamps and gathering back. On beta Lyr that is 1613 uv points
#     against 50 distinct times — 32x fewer Kepler solves — and it takes the likelihood to
#     ~0.05 ms for 7 nights.
#   * Analytic component visibilities. A uniform disk is a Bessel function; there is no
#     mesh, no polygon Fourier transform. Fitting six orbital elements does not need a
#     resolved surface, and nested sampling over 10 parameters is only practical at this
#     cost per likelihood.
#
# CONVENTIONS. Component 1 sits at the origin and component 2 is displaced by the orbit,
# following `orbit_to_rotir_offset` (see src/oichi2_binary.jl):
#
#     cvis = (V1 + f * V2 * exp(-2pi*i*(u*Δra + v*Δdec))) / (1 + f)
#
# with `f` the flux of component 2 relative to component 1, and Δra/Δdec in mas. East is
# -West, matching ROTIR's internal (West, North, toward-observer) frame.

# ---------------------------------------------------------------------------------------
# Component models
# ---------------------------------------------------------------------------------------
"""
    OrbitComponent

A resolved (or unresolved) component of a binary, in the *analytic* visibility sense: it
knows how to turn its own parameters into complex visibilities on a uv grid, and which of
those parameters may be fitted.

Implementations: [`PointSource`](@ref), [`UniformDisk`](@ref), [`LimbDarkenedDisk`](@ref),
[`GaussianDisk`](@ref), [`EllipticalGaussian`](@ref).
"""
abstract type OrbitComponent end

"Unresolved point source. No free parameters."
struct PointSource <: OrbitComponent end

"Uniform disk of angular `diameter` (mas)."
struct UniformDisk <: OrbitComponent
    diameter::Float64
end
UniformDisk(; diameter) = UniformDisk(diameter)

"""
Limb-darkened disk of angular `diameter` (mas). `law` selects the profile and therefore how
many coefficients are used: `:linear` (ld1), `:quadratic` (ld1, ld2), `:power` (ld1).
"""
struct LimbDarkenedDisk <: OrbitComponent
    diameter::Float64
    law::Symbol
    ld1::Float64
    ld2::Float64
end
LimbDarkenedDisk(; diameter, law = :linear, ld1 = 0.3, ld2 = 0.0) =
    LimbDarkenedDisk(diameter, law, ld1, ld2)

"Circular Gaussian of full width at half maximum `fwhm` (mas)."
struct GaussianDisk <: OrbitComponent
    fwhm::Float64
end
GaussianDisk(; fwhm) = GaussianDisk(fwhm)

"""
Elliptical Gaussian: `fwhm` (mas) along the major axis, axis `ratio` in (0,1]. `pa`
(degrees) is the angle of the minor axis, measured from +RA toward +Dec; the major axis is at
position angle `-pa` (North through East).

`pa` is taken modulo 180: the profile is symmetric under a half turn, so a [0,360) range
would hold two exact copies of every solution.
"""
struct EllipticalGaussian <: OrbitComponent
    fwhm::Float64
    ratio::Float64
    pa::Float64
end
EllipticalGaussian(; fwhm, ratio = 1.0, pa = 0.0) = EllipticalGaussian(fwhm, ratio, pa)

# --- introspection: names, values, default bounds -------------------------------------
component_param_names(::PointSource)        = Symbol[]
component_param_names(::UniformDisk)        = [:diameter]
component_param_names(c::LimbDarkenedDisk)  =
    c.law === :quadratic ? [:diameter, :ld1, :ld2] : [:diameter, :ld1]
component_param_names(::GaussianDisk)       = [:fwhm]
component_param_names(::EllipticalGaussian) = [:fwhm, :ratio, :pa]

component_param_values(::PointSource)          = Float64[]
component_param_values(c::UniformDisk)         = [c.diameter]
component_param_values(c::LimbDarkenedDisk)    =
    c.law === :quadratic ? [c.diameter, c.ld1, c.ld2] : [c.diameter, c.ld1]
component_param_values(c::GaussianDisk)        = [c.fwhm]
component_param_values(c::EllipticalGaussian)  = [c.fwhm, c.ratio, c.pa]

"Default (lo, hi) box for a component's own parameters. Override per fit via `bounds`."
component_param_bounds(::PointSource)          = (Float64[], Float64[])
component_param_bounds(::UniformDisk)          = ([0.01], [20.0])
function component_param_bounds(c::LimbDarkenedDisk)
    c.law === :quadratic ? ([0.01, 0.0, -0.5], [20.0, 1.0, 1.0]) : ([0.01, 0.0], [20.0, 1.0])
end
component_param_bounds(::GaussianDisk)         = ([0.01], [20.0])
# pa upper bound is 180, not 360 — see the EllipticalGaussian docstring.
component_param_bounds(::EllipticalGaussian)   = ([0.01, 0.05, 0.0], [20.0, 1.0, 180.0])

# --- the visibilities ------------------------------------------------------------------
# `uv` is 2 x n (metres/wavelength), `ρ` its radius, precomputed once per epoch.
# OITOOLS 0.15 REMOVED the standalone `visibility_*(param, uv)` functions — commit 6f83ff4
# dropped `include("vis_functions.jl")` and their export block together, and
# OITOOLS/docs/src/api/modeling.md records the removal in 0.13.2. The replacement is the model
# API: build a model from a parameter dict, then evaluate it.
#
# EQUIVALENCE WAS MEASURED, not assumed, over 400 random uv points against the old
# implementations recovered from `git show 6f83ff4^:src/vis_functions.jl`:
#
#     visibility_ud([D], uv)           -> "c,ud"      => D                  max|Δ| = 2.2e-16
#     visibility_ldlin([D,u], uv)      -> "c,ldlin"   => D, "c,u"     => u  max|Δ| = 4.4e-16
#     visibility_ldquad([D,u,w], uv)   -> "c,ldquad"  => D, …              max|Δ| = 8.9e-16
#     visibility_ldpow([D,α], uv)      -> "c,ldpow"   => D, "c,alpha" => α  max|Δ| = 6.7e-16
#
# `test/test_orbit_fit.jl` pins these numerically so a future OITOOLS change cannot move them
# quietly. Every component also needs its flux, `"c,f" => 1.0`.
#
# `real.()`, not the raw result: `model_to_vis` returns `ComplexF64` because a model in general
# carries an offset, but these components are CENTRED, so the imaginary part is exactly zero
# (checked: `all(iszero, imag.(v))`). The other `component_vis` methods return real vectors and
# the phase factor in `_orbit_cvis` is what makes the sum complex — keeping this real means the
# promotion happens in one place rather than per component.
_vis(d::Dict{String,Any}, uv) = real.(model_to_vis(dict_to_model(d, String[]), Float64[], uv))

component_vis(::PointSource, p, uv, ρ) = ones(eltype(ρ), length(ρ))
component_vis(::UniformDisk, p, uv, ρ) = _vis(Dict{String,Any}("c,ud" => p[1],
                                                               "c,f"  => 1.0), uv)

function component_vis(c::LimbDarkenedDisk, p, uv, ρ)
    c.law === :linear    && return _vis(Dict{String,Any}("c,ldlin" => p[1], "c,u" => p[2],
                                                         "c,f" => 1.0), uv)
    c.law === :quadratic && return _vis(Dict{String,Any}("c,ldquad" => p[1], "c,u" => p[2],
                                                         "c,w" => p[3], "c,f" => 1.0), uv)
    c.law === :power     && return _vis(Dict{String,Any}("c,ldpow" => p[1], "c,alpha" => p[2],
                                                         "c,f" => 1.0), uv)
    error("LimbDarkenedDisk: unknown law $(c.law) (use :linear, :quadratic or :power)")
end

# NOT the model API, and not OITOOLS at all. A `GaussianDisk` is circular, and the elliptical
# method below already computes exactly that at `ratio = 1, pa = 0` — verified against the old
# `visibility_Gaussian([F, 0, 0], uv)` at FWHM 0.5/3/12 mas to max|Δ| = 0.0/1.1e-16/3.1e-16.
#
# Using it avoids the one genuinely dangerous part of this migration: the model API measures
# position angle from North while the old function measured from the u axis, so the mapping is
# `pa = 90 − φ` and writing the obvious `pa = φ` yields a silently misoriented ellipse (0.024
# on the same test points — wrong, but not an error). A circular Gaussian has no orientation to
# get wrong, so the trap simply does not arise here.
component_vis(c::GaussianDisk, p, uv, ρ) = component_vis(EllipticalGaussian(p[1], 1.0, 0.0),
                                                         (p[1], 1.0, 0.0), uv, ρ)

function component_vis(::EllipticalGaussian, p, uv, ρ)
    fwhm, ratio, pa = p[1], p[2], p[3]
    φ  = deg2rad(pa)
    up =  uv[1, :] .* cos(φ) .+ uv[2, :] .* sin(φ)
    vp = -uv[1, :] .* sin(φ) .+ uv[2, :] .* cos(φ)
    ρe = sqrt.((up .* ratio).^2 .+ vp.^2)
    σ  = fwhm / MAS_PER_RAD_ORBFIT / (2 * sqrt(2 * log(2)))
    return exp.(-2 * (pi * σ)^2 .* ρe.^2)
end

const MAS_PER_RAD_ORBFIT = 2.0626480624709636e8

# ---------------------------------------------------------------------------------------
# Orbital elements
# ---------------------------------------------------------------------------------------
"Orbital elements, in the order the parameter vector uses."
# `domega` (apsidal motion, deg/day) is APPENDED rather than inserted next to the other
# angles, so every existing index keeps its meaning — T0 stays at 7, and only `_unpack`
# moves. It is off by default (0.0) and not in the default free list: it matters over a
# long baseline and is unmeasurable over a short one.
const ORBIT_ELEMENTS = (:a, :i, :Omega, :omega, :e, :P, :T0, :dP, :domega)

const ORBIT_ELEMENT_BOUNDS = Dict(
    :a     => (0.01, 50.0),      # mas
    :i     => (0.0, 180.0),      # deg; >90 is retrograde
    :Omega => (0.0, 180.0),      # deg; see the degeneracy note in `fit_orbit`
    :omega => (0.0, 360.0),      # deg; undefined at e = 0
    :e     => (0.0, 0.95),
    :P     => (0.01, 1.0e5),     # days
    :T0    => (-Inf, Inf),       # filled from P unless given
    :dP    => (-1.0e-4, 1.0e-4), # d/d
    # deg/day. Spica's apsidal motion is ~2.6 deg/yr = 7.1e-3 deg/d (U ~ 105-110 yr),
    # so this box is generous by roughly an order of magnitude either way.
    :domega => (-0.05, 0.05),
)

# ---------------------------------------------------------------------------------------
# Per-epoch precompute
# ---------------------------------------------------------------------------------------
"""
Everything about the data that does not change during a fit: uv coordinates, their radii,
and — the important part — the *distinct* observation times with an index mapping each uv
point back to one of them. That is what collapses the Kepler solves.
"""
struct OrbitFitData{T}
    data::Vector                 # the OIdata blocks
    uv::Vector{Matrix{T}}
    rho::Vector{Vector{T}}
    tsrt::Vector{Vector{Float64}}  # distinct times per epoch (JD)
    tidx::Vector{Vector{Int}}      # uv point -> index into tsrt
    nv2::Int
    nt3amp::Int
    nt3phi::Int
    ntot::Int
end

"""
    orbit_fit_data(data; T = Float64) -> OrbitFitData

Precompute the per-epoch quantities a fit needs. `data` is a single `OIdata` or a vector of
them (one per night/epoch), as returned by `readoifits`.

Times are reduced to their distinct values: an OIFITS night typically has thousands of uv
points but only tens of exposures, and the orbit only has to be solved once per exposure.
"""
function orbit_fit_data(data; T::Type = Float64)
    blocks = data isa AbstractVector ? collect(data) : [data]
    uv   = Vector{Matrix{T}}(undef, length(blocks))
    rho  = Vector{Vector{T}}(undef, length(blocks))
    tsrt = Vector{Vector{Float64}}(undef, length(blocks))
    tidx = Vector{Vector{Int}}(undef, length(blocks))
    for (k, d) in enumerate(blocks)
        uv[k]  = T.(d.uv)
        rho[k] = sqrt.(uv[k][1, :].^2 .+ uv[k][2, :].^2)
        t      = uv_times(d)
        u      = unique(t)
        sort!(u)
        pos    = Dict(v => j for (j, v) in enumerate(u))
        tsrt[k] = u
        tidx[k] = [pos[x] for x in t]
    end
    OrbitFitData{T}(blocks, uv, rho, tsrt, tidx,
                    sum(d.nv2 for d in blocks), sum(d.nt3amp for d in blocks),
                    sum(d.nt3phi for d in blocks),
                    sum(d.nv2 + d.nt3amp + d.nt3phi for d in blocks))
end

# ---------------------------------------------------------------------------------------
# Parameter vector
# ---------------------------------------------------------------------------------------
"""
    OrbitFitSpec

The layout of a fit: which quantities exist, their starting values, their bounds, and
whether each is free, fixed or tied. The parameter vector is

    [a, i, Omega, omega, e, P, T0, dP, domega,
     <comp1 params>...,  <comp2 params>...,  f]

with component parameter names prefixed `c1_` / `c2_` so `free = [:c2_pa]` is unambiguous.
The orbital block is `ORBIT_ELEMENTS`; index off `length(ORBIT_ELEMENTS)` rather than a
literal, since it has grown once already (`domega` was appended for apsidal motion).

Every parameter is exactly one of:

* **free** — listed in `free`, sampled within `lo`/`hi`;
* **fixed** — absent from both, holding its starting value;
* **tied** — listed in `ties`, computed from the others by [`compile_ties`](@ref).

Build one with [`orbit_fit_spec`](@ref), and turn a free-parameter vector into a full one
with [`resolve_params`](@ref) — which is the only supported way to do so, because it is
also what applies the ties.
"""
struct OrbitFitSpec{TT}
    names::Vector{Symbol}
    values::Vector{Float64}
    lo::Vector{Float64}
    hi::Vector{Float64}
    free::Vector{Int}
    comp1::OrbitComponent
    comp2::OrbitComponent
    n1::Int          # number of comp1 parameters
    n2::Int
    # Radial velocities. `rvinst` lists the table's instruments, the first being the
    # reference (no offset); `nrv` is the length of the RV block that follows `f`:
    # K1, K2, gamma, one dgamma_<instrument> per further instrument, jit_a, jit_b. Both are
    # empty/0 for a fit without velocities.
    rvinst::Vector{String}
    nrv::Int
    # true: the block starts with plx, q (the amplitudes follow from the orbit); false: K1, K2.
    rvphys::Bool
    # Parameters computed from the others (src/orbit_ties.jl). The type parameter matters:
    # declared as plain `::OrbitTies` this field is abstract, every `spec.ties.fn` call goes
    # through dynamic dispatch, and resolving a one-line tie costs 5x the time and an extra
    # allocation per likelihood evaluation.
    ties::TT
end

"""
    orbit_fit_spec(comp1, comp2; elements, flux_ratio = 1.0, free = nothing,
                   bounds = Dict(), ties = Dict(), rv = nothing,
                   rv_params = NamedTuple()) -> OrbitFitSpec

Assemble the parameter layout. `elements` is a NamedTuple of orbital elements (any of
$(ORBIT_ELEMENTS); missing ones take sensible defaults, `e = 0`, `omega = 90`, `dP = 0`).

`free` is a vector of parameter names to fit; the default frees the four elements a
two-component interferometric orbit actually constrains — `a, i, Omega, T0` — plus the flux
ratio and every component parameter. `bounds` overrides individual boxes, e.g.
`Dict(:a => (0.5, 1.5))`.

With `rv` (an [`RVData`](@ref)) the vector gains an RV block after `f`:

| name | meaning | default start | box |
|---|---|---|---|
| `plx` | parallax, mas | from [`rv_start`](@ref) and the elements | [1e-4, 1e4] |
| `q` | mass ratio M2/M1 | idem | [0.01, 100] |
| `gamma` | systemic velocity, km/s | from `rv_start` | [−1000, 1000] |
| `dgamma_<instrument>` | offset of each instrument after the first, km/s | 0 | [−50, 50] |
| `jit_a`, `jit_b` | jitter of each star, km/s, added in quadrature to σ | 0 | [0, 50] |

The semi-amplitudes are not parameters: with the interferometry fitting `a`, `i`, `P` and `e`,
they follow from `plx` and `q` through [`rv_amplitudes`](@ref). `rv_parametrization =
:amplitudes` fits `K1`, `K2` (km/s, box [0, 1000]) in their place instead.

`rv_params` overrides starting values by name, e.g. `(plx = 12.99, jit_b = 2.0)`. The default
`free` then also includes `plx, q, gamma` (or `K1, K2, gamma`) and the offsets; the jitters are
fitted only when listed. Velocities break the astrometric (Omega, omega) → (Omega+180°, omega+180°) mirror,
so `Omega`'s box becomes [0, 360) unless `bounds` sets it.
"""
function orbit_fit_spec(comp1::OrbitComponent, comp2::OrbitComponent;
                        elements, flux_ratio::Real = 1.0,
                        free = nothing, bounds = Dict{Symbol,Tuple{Float64,Float64}}(),
                        ties = Dict{Symbol,String}(), rv = nothing,
                        rv_params = NamedTuple(), rv_parametrization::Symbol = :physical)
    rv_parametrization in (:physical, :amplitudes) ||
        error("orbit_fit_spec: rv_parametrization is :physical or :amplitudes")
    rvphys = rv_parametrization === :physical
    el = merge((a = 1.0, i = 90.0, Omega = 0.0, omega = 90.0, e = 0.0,
                P = 1.0, T0 = 0.0, dP = 0.0, domega = 0.0), elements)
    haskey(elements, :P) || error("orbit_fit_spec: `elements` must include the period P")

    names  = collect(ORBIT_ELEMENTS)
    values = Float64[getfield(el, k) for k in ORBIT_ELEMENTS]
    lo     = Float64[ORBIT_ELEMENT_BOUNDS[k][1] for k in ORBIT_ELEMENTS]
    hi     = Float64[ORBIT_ELEMENT_BOUNDS[k][2] for k in ORBIT_ELEMENTS]
    # T0 defaults to one full period centred on the supplied value. One period, not more:
    # combined with Omega restricted to [0,180) this is exactly one fundamental domain of
    # the (Omega, T0) -> (Omega+180, T0+P/2) degeneracy that holds at e = 0.
    lo[7], hi[7] = el.T0 - el.P/2, el.T0 + el.P/2

    for (pref, c) in ((:c1_, comp1), (:c2_, comp2))
        cn, cv = component_param_names(c), component_param_values(c)
        cl, ch = component_param_bounds(c)
        append!(names,  Symbol.(pref, cn))
        append!(values, cv); append!(lo, cl); append!(hi, ch)
    end
    push!(names, :f); push!(values, Float64(flux_ratio)); push!(lo, 0.0); push!(hi, 100.0)

    rvinst = String[]
    if rv !== nothing
        rvinst = rv_instruments(rv)
        offn = rv_offset_name.(rvinst[2:end])
        allunique(offn) || error("orbit_fit_spec: instrument names $(rvinst) do not give " *
                                 "distinct parameter names")
        st = rv_start(rv)
        if rvphys
            # parallax and mass ratio that give the data's amplitudes at these elements
            p0, q0 = rv_physical(el.a, el.i, el.P, el.e, st.K1, st.K2)
            (isfinite(p0) && p0 > 0) || (p0 = 1.0)
            (isfinite(q0) && q0 > 0) || (q0 = 1.0)
            rvn = Symbol[:plx, :q, :gamma, offn..., :jit_a, :jit_b]
            rvv = Float64[p0, clamp(q0, 0.01, 100.0), st.gamma, zeros(length(offn))..., 0.0, 0.0]
            rvl = Float64[1e-4, 0.01, -1000.0, fill(-50.0, length(offn))..., 0.0, 0.0]
            rvh = Float64[1e4, 100.0, 1000.0, fill(50.0, length(offn))..., 50.0, 50.0]
        else
            rvn = Symbol[:K1, :K2, :gamma, offn..., :jit_a, :jit_b]
            rvv = Float64[st.K1, st.K2, st.gamma, zeros(length(offn))..., 0.0, 0.0]
            rvl = Float64[0.0, 0.0, -1000.0, fill(-50.0, length(offn))..., 0.0, 0.0]
            rvh = Float64[1000.0, 1000.0, 1000.0, fill(50.0, length(offn))..., 50.0, 50.0]
        end
        for (k, v) in pairs(rv_params)
            j = findfirst(==(Symbol(k)), rvn)
            j === nothing && error("orbit_fit_spec: unknown RV parameter $k in `rv_params`; " *
                                   "available: $(join(string.(rvn), ", "))")
            rvv[j] = Float64(v)
        end
        append!(names, rvn); append!(values, rvv); append!(lo, rvl); append!(hi, rvh)
        hi[3] = 360.0                                # Omega: the mirror is broken
    elseif !isempty(rv_params)
        error("orbit_fit_spec: `rv_params` given without `rv`")
    end
    nrv = length(names) - (length(ORBIT_ELEMENTS) + length(component_param_names(comp1)) +
                           length(component_param_names(comp2)) + 1)

    for (k, (l, h)) in bounds
        j = findfirst(==(k), names)
        j === nothing && error("orbit_fit_spec: unknown parameter $k in `bounds`")
        lo[j], hi[j] = l, h
    end

    n1 = length(component_param_names(comp1))
    n2 = length(component_param_names(comp2))
    if free === nothing
        # e and omega are excluded by default: at e = 0 omega is undefined, and e itself is
        # poorly constrained by a short baseline. P belongs to an ephemeris, not to a few
        # nights of visibilities. Free them explicitly if you mean to.
        deflt = Symbol[:a, :i, :Omega, :T0, :f]
        append!(deflt, Symbol.(:c1_, component_param_names(comp1)))
        append!(deflt, Symbol.(:c2_, component_param_names(comp2)))
        rv === nothing || append!(deflt, [(rvphys ? (:plx, :q) : (:K1, :K2))..., :gamma,
                                          rv_offset_name.(rvinst[2:end])...])
        free = deflt
    end
    compiled_ties = compile_ties(names, ties)
    # A tie overwrites its target on every evaluation, so a parameter that is both free and
    # tied would have its sampled value silently discarded — the fit would report an
    # uncertainty for a quantity the likelihood never saw. Refuse rather than mislead.
    tied_and_free = intersect(Set(compiled_ties.names), Set(free))
    isempty(tied_and_free) ||
        error("orbit_fit_spec: $(join(string.(collect(tied_and_free)), ", ")) " *
              "listed as both free and tied; a tie overwrites the sampled value")

    idx = Int[]
    avail = join(string.(names), ", ")
    for k in free
        j = findfirst(==(k), names)
        j === nothing && error("orbit_fit_spec: unknown free parameter $k; available: $avail")
        push!(idx, j)
    end
    for j in idx
        lo[j] <= values[j] <= hi[j] ||
            error("orbit_fit_spec: start value $(values[j]) for $(names[j]) is outside " *
                  "its bounds [$(lo[j]), $(hi[j])]")
    end
    OrbitFitSpec(names, values, lo, hi, sort!(idx), comp1, comp2, n1, n2, rvinst, nrv,
                 rvphys, compiled_ties)
end

# ---------------------------------------------------------------------------------------
# Forward model and chi2
# ---------------------------------------------------------------------------------------
"Split `θ` into (orbital elements NamedTuple, comp1 params, comp2 params, flux ratio)."
@inline function _unpack(spec::OrbitFitSpec, θ)
    a, i, Om, om, e, P, T0, dP, dom =
        θ[1], θ[2], θ[3], θ[4], θ[5], θ[6], θ[7], θ[8], θ[9]
    nel = length(ORBIT_ELEMENTS)
    p1 = view(θ, (nel + 1):(nel + spec.n1))
    p2 = view(θ, (nel + 1 + spec.n1):(nel + spec.n1 + spec.n2))
    f  = θ[nel + spec.n1 + spec.n2 + 1]
    return (a = a, i = i, Omega = Om, omega = om, e = e, P = P, T0 = T0, dP = dP,
            domega = dom), p1, p2, f
end

# The elements as `binary_RV` / `orbit_to_rotir_offset` read them.
_orbit_bp(el) = (i = el.i, Ω = el.Omega, ω = el.omega, P = el.P, a = el.a, e = el.e,
                 T0 = el.T0, q = 1.0, dP = el.dP, dω = el.domega)

"`(K1, K2, gamma, offsets, jitter)` from a resolved parameter vector with an RV block."
function _rv_block(spec::OrbitFitSpec, θ)
    j = length(ORBIT_ELEMENTS) + spec.n1 + spec.n2 + 1      # index of f
    noff = length(spec.rvinst) - 1
    offsets = Dict{String,eltype(θ)}(spec.rvinst[m + 1] => θ[j + 3 + m] for m in 1:noff)
    jitter  = (a = θ[j + 4 + noff], b = θ[j + 5 + noff])
    # ORBIT_ELEMENTS order: a, i, Omega, omega, e, P, …
    K1, K2 = spec.rvphys ? rv_amplitudes(θ[1], θ[2], θ[6], θ[5], θ[j + 1], θ[j + 2]) :
                           (θ[j + 1], θ[j + 2])
    return K1, K2, θ[j + 3], offsets, jitter
end

function _rv_check(spec::OrbitFitSpec, rv)
    spec.nrv > 0 || error("RV data given, but the spec has no RV block: build it with " *
                          "orbit_fit_spec(...; rv = rv)")
    rv_instruments(rv) == spec.rvinst ||
        error("the RV table's instruments $(rv_instruments(rv)) differ from the spec's " *
              "$(spec.rvinst)")
    return nothing
end

"""
    orbit_model_cvis(spec, θ, fd, k) -> Vector{Complex}

Complex visibilities for epoch `k`. Component 1 sits at the origin, component 2 is
displaced by the orbit; `f` is the flux of 2 relative to 1.
"""
function orbit_model_cvis(spec::OrbitFitSpec, θ, fd::OrbitFitData, k::Int)
    # Public entry point, so it resolves ties itself; `orbit_chi2` resolves them once for
    # all epochs and calls `_orbit_cvis` directly. `apply_ties` returns θ untouched when
    # there are none, so an untied model pays nothing here.
    _orbit_cvis(spec, apply_ties(spec.ties, θ), fd, k)
end

"""
    resolve_params(spec, p) -> θ

Build the full parameter vector from the free subset `p`: start from the fixed values,
scatter the free ones, then apply the ties. This is the ONLY place a full parameter vector
is constructed, which is the point — when the free-to-full mapping is written out at each
call site instead, it is possible (and did happen) to build one that skips the ties, so a
fit reports tied parameters still holding their starting values while quoting a χ² computed
from the correct ones.

Follows the scheme OITOOLS uses in `resolvers.jl`: free values, fixed values and derived
expressions resolved together in one pass rather than layered as separate copies.
"""
function resolve_params(spec::OrbitFitSpec, p)
    θ = copy(spec.values)
    @inbounds for (j, i) in enumerate(spec.free); θ[i] = p[j]; end
    isempty(spec.ties) && return θ
    vals = spec.ties.fn(θ)
    @inbounds for (j, i) in enumerate(spec.ties.targets); θ[i] = vals[j]; end
    return θ
end

"Epoch-`k` visibilities for a parameter vector whose ties are ALREADY resolved."
function _orbit_cvis(spec::OrbitFitSpec, θ, fd::OrbitFitData, k::Int)
    el, p1, p2, f = _unpack(spec, θ)
    bp = _orbit_bp(el)
    # One Kepler solve per DISTINCT time, then gather back to the uv points.
    ows, decs = orbit_to_rotir_offset(bp, fd.tsrt[k])
    ras = -ows                                   # East = -West
    uv, ρ, ti = fd.uv[k], fd.rho[k], fd.tidx[k]
    v1 = component_vis(spec.comp1, p1, uv, ρ)
    v2 = component_vis(spec.comp2, p2, uv, ρ)
    ph = cis.((-2pi / MAS_PER_RAD_ORBFIT) .*
              (view(uv, 1, :) .* ras[ti] .+ view(uv, 2, :) .* decs[ti]))
    return (v1 .+ f .* v2 .* ph) ./ (1 + f)
end

"""
    orbit_chi2(spec, θ, fd; weights = OI_DEFAULT_WEIGHTS, split = false, rv = nothing)

Total over all epochs, or `(χ²_V2, χ²_T3amp, χ²_T3phi)` with `split = true`.

With `rv` (the spec must have been built with the same table) the velocities are added:
the total then includes [`rv_nll2`](@ref), which equals χ²_RV when the jitters are zero, and
`split = true` returns a fourth entry, χ²_RV with jitter-inflated errors.
"""
orbit_chi2(spec::OrbitFitSpec, θ, fd::OrbitFitData;
           weights = OI_DEFAULT_WEIGHTS, split::Bool = false, rv = nothing) =
    _orbit_chi2(spec, apply_ties(spec.ties, θ), fd; weights = weights, split = split, rv = rv)

"χ² for a parameter vector whose ties are ALREADY resolved (see `resolve_params`)."
function _orbit_chi2(spec::OrbitFitSpec, θr, fd::OrbitFitData;
                     weights = OI_DEFAULT_WEIGHTS, split::Bool = false, rv = nothing)
    c = zeros(Float64, 3)
    # All-zero weights skip the visibilities entirely: an RV-only objective then costs the
    # Kepler solves for the velocities and nothing else.
    oi = any(!iszero, (weights[1], weights[2], weights[3]))
    for k in (oi ? eachindex(fd.data) : 1:0)
        cvis = _orbit_cvis(spec, θr, fd, k)
        d = fd.data[k]
        v2m, t3am, t3pm = cvis_to_obs(cvis, d)
        c[1] += sum(abs2, (v2m  .- d.v2)    ./ d.v2_err)
        c[2] += sum(abs2, (t3am .- d.t3amp) ./ d.t3amp_err)
        c[3] += sum(abs2, mod360(t3pm .- d.t3phi) ./ d.t3phi_err)
    end
    c .*= (weights[1], weights[2], weights[3])
    if rv !== nothing
        _rv_check(spec, rv)
        el, _, _, _ = _unpack(spec, θr)
        K1, K2, γ, off, jit = _rv_block(spec, θr)
        bp = _orbit_bp(el)
        split && return (c[1], c[2], c[3],
                         rv_chi2(bp, rv; K1 = K1, K2 = K2, gamma = γ, offsets = off, jitter = jit))
        return sum(c) + rv_nll2(bp, rv; K1 = K1, K2 = K2, gamma = γ, offsets = off, jitter = jit)
    end
    return split ? (c[1], c[2], c[3]) : sum(c)
end

"""
    orbit_residuals(spec, θ, fd; rv = nothing) -> Vector

Normalised residuals `(model − data)/σ` for a full parameter vector `θ` (ties are applied
here): V², T3amp and T3φ (wrapped to ±180°), epoch by epoch, then with `rv` the velocity
residuals `(model − rv)/√(σ² + s²)`. `sum(abs2, r)` equals the χ² parts of
[`orbit_chi2`](@ref) with default weights.
"""
function orbit_residuals(spec::OrbitFitSpec, θ, fd::OrbitFitData; rv = nothing)
    θr = apply_ties(spec.ties, θ)
    out = Float64[]
    for k in eachindex(fd.data)
        d = fd.data[k]
        v2m, t3am, t3pm = cvis_to_obs(_orbit_cvis(spec, θr, fd, k), d)
        append!(out, (v2m .- d.v2) ./ d.v2_err)
        append!(out, (t3am .- d.t3amp) ./ d.t3amp_err)
        append!(out, mod360(t3pm .- d.t3phi) ./ d.t3phi_err)
    end
    if rv !== nothing
        _rv_check(spec, rv)
        el, _, _, _ = _unpack(spec, θr)
        K1, K2, γ, off, jit = _rv_block(spec, θr)
        append!(out, .-rv_residuals(_orbit_bp(el), rv; K1 = K1, K2 = K2, gamma = γ,
                                    offsets = off, jitter = jit))
    end
    return out
end

# ---------------------------------------------------------------------------------------
# The fitter
# ---------------------------------------------------------------------------------------
"""
    orbit_periodic(spec::OrbitFitSpec) -> Vector{Int}

The free parameters whose box is one full period of an angle, as indices into `spec.free`:
`Omega` and `omega` over 360°, and a component's `pa` over 180°, the half turn under which
an elliptical profile is unchanged. `Omega` is periodic only with radial velocities, whose
box is [0, 360); without them it spans [0, 180), one half of the astrometric mirror, and is
not.
"""
function orbit_periodic(spec::OrbitFitSpec)
    out = Int[]
    for (j, k) in enumerate(spec.free)
        n = spec.names[k]
        period = n in (:Omega, :omega) ? 360.0 : endswith(String(n), "_pa") ? 180.0 : 0.0
        period > 0 && isapprox(spec.hi[k] - spec.lo[k], period; atol = 1e-9) && push!(out, j)
    end
    return out
end

"""
    fit_orbit(data, comp1, comp2; elements, kwargs...) -> NamedTuple

Fit the relative orbit of a two-component system to interferometric data.

`data` is a single `OIdata` or a vector of them (one per epoch), straight from
`readoifits`. `comp1` and `comp2` are [`OrbitComponent`](@ref)s — component 1 sits at the
origin, component 2 is displaced by the orbit.

# Keywords
* `elements` — NamedTuple of starting orbital elements. `P` is required; the rest default
  to `a = 1, i = 90, Omega = 0, omega = 90, e = 0, T0 = 0, dP = 0`.
* `flux_ratio = 1.0` — flux of component 2 relative to component 1.
* `free` — parameter names to fit. Default: `a, i, Omega, T0, f` plus every component
  parameter. Component parameters are prefixed `c1_` / `c2_`.
* `ties` — parameters computed from the others instead of fitted, as
  `Dict(:c2_pa => "-Omega")`. See [`compile_ties`](@ref). A parameter cannot be both free
  and tied.
* `bounds` — `Dict(:name => (lo, hi))` overriding individual boxes.
* `method` — `:neldermead` (NLopt, fast, returns a point estimate), or `:ultranest` or
  `:nautilus` (nested sampling, returns a posterior and log Z). The samplers treat a free
  angle whose box is one full period as periodic ([`orbit_periodic`](@ref)), so a mode
  across the box edge is not split, and quote its median and interval around the circle.
* `model` — `:analytic` (default) or `:tessellated`.
* `weights` — observable weights, default `OI_DEFAULT_WEIGHTS` (V², T3amp, T3φ on).
* `rv` — radial velocities ([`RVData`](@ref)), fitted together with the interferometry; see
  [`orbit_fit_spec`](@ref) for the parameters this adds (K1, K2, gamma, per-instrument
  offsets, jitters). `rv_params` sets their starting values; `rv_parametrization` chooses
  `plx, q` (default) or `K1, K2`.
* `maxeval` (`:neldermead`), `min_num_live_points`, `use_stepsampler` (`:ultranest`),
  `n_live`, `n_eff` (`:nautilus`) — passed to the chosen fitter.

# Returns
A NamedTuple with `params` (full vector), `elements`, `chi2`, `chi2_split`, `chi2_red`,
`spec`, `data`, `rv`, and — for `:ultranest` / `:nautilus` — `posterior`, `logz`, `q16`,
`q84`. With velocities, `chi2_split` has a fourth entry (χ²_RV), `chi2_red` divides by the
interferometric plus RV point count, and `rvpar` holds the RV block by name. The fit
minimises χ²_OI + [`rv_nll2`](@ref), which is the plain χ² when the jitters are fixed at 0.

# Degeneracies worth knowing
At `e = 0` the argument of periastron `omega` is undefined and is not fitted by default.
The residual `(Omega, T0) -> (Omega + 180°, T0 + P/2)` degeneracy is broken by restricting
`Omega` to `[0, 180)` while `T0` spans one full period — together exactly one fundamental
domain. Radial velocities break the (Omega, omega) mirror, and with `rv` Omega spans
`[0, 360)`. An elliptical Gaussian's `pa` is likewise capped at 180°, since the profile is
symmetric under a half turn.

`P` is not free by default: for most systems an eclipse or spectroscopic ephemeris pins it
far better than a few nights of visibilities can, and freeing it mostly adds an
unconstrained direction. Free it explicitly as a consistency check.

# Example
```julia
data = readoifits("betlyr.oifits")[1, 1]
res = fit_orbit(data,
                UniformDisk(diameter = 0.57),
                EllipticalGaussian(fwhm = 0.90, ratio = 0.61, pa = 73.7);
                elements = (P = 12.9414, e = 0.0, a = 0.865,
                            i = 93.5, Omega = 73.7, T0 = 2454283.043),
                flux_ratio = 0.81,
                method = :neldermead)
res.elements.a, res.chi2_red
```
"""
function fit_orbit(data, comp1::OrbitComponent, comp2::OrbitComponent;
                   elements, flux_ratio::Real = 1.0, free = nothing,
                   bounds = Dict{Symbol,Tuple{Float64,Float64}}(),
                   ties = Dict{Symbol,String}(),
                   rv = nothing, rv_params = NamedTuple(),
                   rv_parametrization::Symbol = :physical,
                   method::Symbol = :neldermead, model::Symbol = :analytic,
                   weights = OI_DEFAULT_WEIGHTS, maxeval::Int = 20_000,
                   threaded::Bool = true,
                   min_num_live_points::Int = 400, use_stepsampler::Bool = true,
                   n_live::Int = 1000, n_eff::Real = 10_000,
                   verbose::Bool = true)
    model === :analytic || model === :tessellated ||
        error("fit_orbit: model must be :analytic or :tessellated (got $model)")
    if model === :tessellated
        error("""
              fit_orbit: model = :tessellated is not implemented yet.

              The analytic path fits components as uniform/limb-darkened disks and
              Gaussians, which is what an orbit fit needs. A tessellated fit would carry
              ROTIR `stellar_geometry` meshes through `binary_cvis`, giving Roche shapes,
              gravity darkening and irradiation — but it needs its own parameterisation
              (rpole, tpole, fillout, ...) rather than the component library here, and it
              is orders of magnitude slower per likelihood.

              For a tessellated binary today, drive `create_binary_geometry` and
              `binary_chi2_f` directly — see demos/spica_binary_roche.jl.""")
    end

    fd   = orbit_fit_data(data)
    spec = orbit_fit_spec(comp1, comp2; elements = elements, flux_ratio = flux_ratio,
                          free = free, bounds = bounds, ties = ties,
                          rv = rv, rv_params = rv_params,
                          rv_parametrization = rv_parametrization)
    idx  = spec.free
    θ0   = copy(spec.values)
    nrvp = rv === nothing ? 0 : length(rv)
    npts = fd.ntot + nrvp

    chi2_of = let spec = spec, fd = fd, weights = weights, rv = rv
        p -> _orbit_chi2(spec, resolve_params(spec, p), fd; weights = weights, rv = rv)
    end

    if verbose
        @printf("fit_orbit: %d epochs, %d points (V² %d, T3amp %d, T3φ %d, RV %d)\n",
                length(fd.data), npts, fd.nv2, fd.nt3amp, fd.nt3phi, nrvp)
        @printf("           %d free of %d parameters: %s\n", length(idx), length(spec.names),
                join(string.(spec.names[idx]), ", "))
        @printf("           starting χ²/n = %.3f\n", chi2_of(θ0[idx]) / npts)
    end

    local best, extra
    if method === :neldermead
        opt = NLopt.Opt(:LN_NELDERMEAD, length(idx))
        NLopt.lower_bounds!(opt, spec.lo[idx]); NLopt.upper_bounds!(opt, spec.hi[idx])
        NLopt.maxeval!(opt, maxeval); NLopt.xtol_rel!(opt, 1e-8)
        NLopt.min_objective!(opt, (p, g) -> chi2_of(p))
        (_, pbest, _) = NLopt.optimize(opt, θ0[idx])
        best = pbest; extra = NamedTuple()
    elseif method in (:ultranest, :nautilus)
        # `threaded` nests the sampler's concurrency inside the χ²'s own threading, which is
        # safe since `:t3`'s scratch moved into a pool (see the note on `fit_sphere_ld`). This
        # objective mostly does not reach `:t3` at all — `model = :analytic` goes through
        # OITOOLS' `dict_to_model`/`model_to_vis` — so that path was checked separately:
        # 800 concurrent `component_vis` calls, all bit-identical to the serial answer.
        per = orbit_periodic(spec)
        verbose && !isempty(per) &&
            println("           periodic: ", join(string.(spec.names[idx[per]]), ", "))
        r = _fit_nested(method, chi2_of, string.(spec.names[idx]), spec.lo[idx], spec.hi[idx];
                        min_num_live_points = min_num_live_points, threaded = threaded,
                        use_stepsampler = use_stepsampler, n_live = n_live, n_eff = n_eff,
                        periodic = per, verb = verbose)
        best  = r.median
        extra = (posterior = r.samples, logz = r.logz, logzerr = r.logzerr,
                 q16 = r.q16, q84 = r.q84)
    else
        error("fit_orbit: method must be :neldermead, :ultranest or :nautilus " *
              "(got $method)")
    end

    θbest = resolve_params(spec, best)
    cs    = _orbit_chi2(spec, θbest, fd; weights = weights, split = true, rv = rv)
    el    = NamedTuple{ORBIT_ELEMENTS}(ntuple(j -> θbest[j], length(ORBIT_ELEMENTS)))
    verbose && @printf("           final    χ²/n = %.3f  (V² %.2f, T3amp %.2f, T3φ %.2f%s)\n",
                       sum(cs)/npts, cs[1]/max(fd.nv2,1), cs[2]/max(fd.nt3amp,1),
                       cs[3]/max(fd.nt3phi,1),
                       rv === nothing ? "" : @sprintf(", RV %.2f", cs[4]/max(nrvp,1)))
    rvpar = NamedTuple()
    if rv !== nothing
        j = length(spec.names) - spec.nrv
        nt = NamedTuple{Tuple(spec.names[j+1:end])}(Tuple(θbest[j+1:end]))
        if spec.rvphys
            # what plx and q mean: the amplitudes, the distance and the dynamical masses
            K1, K2 = rv_amplitudes(el.a, el.i, el.P, el.e, nt.plx, nt.q)
            dm = dynamical_masses(el.a, el.P, nt.plx, nt.q)
            nt = merge(nt, (; K1, K2, M = dm.M, M1 = dm.M1, M2 = dm.M2, d = dm.d))
        end
        rvpar = (rvpar = nt,)
    end
    return merge((params = θbest, names = spec.names, elements = el,
                  chi2 = sum(cs), chi2_split = cs, chi2_red = sum(cs)/npts,
                  spec = spec, data = fd, rv = rv), rvpar, extra)
end
