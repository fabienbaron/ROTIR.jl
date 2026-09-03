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
Elliptical Gaussian: `fwhm` (mas) along the major axis, axis `ratio` in (0,1], major axis
at position angle `pa` (degrees).

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
component_vis(::PointSource, p, uv, ρ) = ones(eltype(ρ), length(ρ))
component_vis(::UniformDisk, p, uv, ρ) = visibility_ud([p[1]], uv)

function component_vis(c::LimbDarkenedDisk, p, uv, ρ)
    c.law === :linear    && return visibility_ldlin([p[1], p[2]], uv)
    c.law === :quadratic && return visibility_ldquad([p[1], p[2], p[3]], uv)
    c.law === :power     && return visibility_ldpow([p[1], p[2]], uv)
    error("LimbDarkenedDisk: unknown law $(c.law) (use :linear, :quadratic or :power)")
end

# OITOOLS' `visibility_Gaussian` is ELLIPTICAL: (FWHM, inclination_deg, PA_deg). A circular
# Gaussian is the i = 0 case, where the aspect factor cos(i) is 1 and PA drops out.
component_vis(::GaussianDisk, p, uv, ρ) = visibility_Gaussian([p[1], 0.0, 0.0], uv)

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
        t      = _uv_times(d)
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

"""
Observation time (JD) for every uv point of an `OIdata` block.

The uv table is NOT a plain [V2; T3_1; T3_2; T3_3] concatenation — `OIdata` carries index
arrays that scatter each observable into a shared `nuv`-length uv list, and a single uv
point can be shared between observables. Times must be scattered the same way, or they are
silently assigned to the wrong baselines.
"""
function _uv_times(d)
    t = zeros(Float64, d.nuv)
    t[d.indx_v2]   .= d.v2_mjd
    t[d.indx_t3_1] .= d.t3_mjd
    t[d.indx_t3_2] .= d.t3_mjd
    t[d.indx_t3_3] .= d.t3_mjd
    return t .+ 2400000.5
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
    # Parameters computed from the others (src/orbit_ties.jl). The type parameter matters:
    # declared as plain `::OrbitTies` this field is abstract, every `spec.ties.fn` call goes
    # through dynamic dispatch, and resolving a one-line tie costs 5x the time and an extra
    # allocation per likelihood evaluation.
    ties::TT
end

"""
    orbit_fit_spec(comp1, comp2; elements, flux_ratio = 1.0, free = nothing,
                   bounds = Dict()) -> OrbitFitSpec

Assemble the parameter layout. `elements` is a NamedTuple of orbital elements (any of
$(ORBIT_ELEMENTS); missing ones take sensible defaults, `e = 0`, `omega = 90`, `dP = 0`).

`free` is a vector of parameter names to fit; the default frees the four elements a
two-component interferometric orbit actually constrains — `a, i, Omega, T0` — plus the flux
ratio and every component parameter. `bounds` overrides individual boxes, e.g.
`Dict(:a => (0.5, 1.5))`.
"""
function orbit_fit_spec(comp1::OrbitComponent, comp2::OrbitComponent;
                        elements, flux_ratio::Real = 1.0,
                        free = nothing, bounds = Dict{Symbol,Tuple{Float64,Float64}}(),
                        ties = Dict{Symbol,String}())
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
    OrbitFitSpec(names, values, lo, hi, sort!(idx), comp1, comp2, n1, n2, compiled_ties)
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
    f  = θ[end]
    return (a = a, i = i, Omega = Om, omega = om, e = e, P = P, T0 = T0, dP = dP,
            domega = dom), p1, p2, f
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
    bp = (i = el.i, Ω = el.Omega, ω = el.omega, P = el.P, a = el.a, e = el.e,
          T0 = el.T0, q = 1.0, dP = el.dP, dω = el.domega)
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
    orbit_chi2(spec, θ, fd; weights = OI_DEFAULT_WEIGHTS, split = false)

Total χ² over all epochs, or `(χ²_V2, χ²_T3amp, χ²_T3phi)` with `split = true`.
"""
orbit_chi2(spec::OrbitFitSpec, θ, fd::OrbitFitData;
           weights = OI_DEFAULT_WEIGHTS, split::Bool = false) =
    _orbit_chi2(spec, apply_ties(spec.ties, θ), fd; weights = weights, split = split)

"χ² for a parameter vector whose ties are ALREADY resolved (see `resolve_params`)."
function _orbit_chi2(spec::OrbitFitSpec, θr, fd::OrbitFitData;
                     weights = OI_DEFAULT_WEIGHTS, split::Bool = false)
    c = zeros(Float64, 3)
    for k in eachindex(fd.data)
        cvis = _orbit_cvis(spec, θr, fd, k)
        d = fd.data[k]
        v2m, t3am, t3pm = cvis_to_obs(cvis, d)
        c[1] += sum(abs2, (v2m  .- d.v2)    ./ d.v2_err)
        c[2] += sum(abs2, (t3am .- d.t3amp) ./ d.t3amp_err)
        c[3] += sum(abs2, mod360(t3pm .- d.t3phi) ./ d.t3phi_err)
    end
    c .*= (weights[1], weights[2], weights[3])
    return split ? (c[1], c[2], c[3]) : sum(c)
end

# ---------------------------------------------------------------------------------------
# The fitter
# ---------------------------------------------------------------------------------------
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
* `method` — `:neldermead` (NLopt, fast, returns a point estimate) or `:ultranest`
  (nested sampling, returns a posterior and log Z).
* `model` — `:analytic` (default) or `:tessellated`.
* `weights` — observable weights, default `OI_DEFAULT_WEIGHTS` (V², T3amp, T3φ on).
* `maxeval`, `min_num_live_points`, `use_stepsampler` — passed to the chosen fitter.

# Returns
A NamedTuple with `params` (full vector), `elements`, `chi2`, `chi2_split`, `chi2_red`,
`spec`, `data`, and — for `:ultranest` — `posterior`, `logz`, `q16`, `q84`.

# Degeneracies worth knowing
At `e = 0` the argument of periastron `omega` is undefined and is not fitted by default.
The residual `(Omega, T0) -> (Omega + 180°, T0 + P/2)` degeneracy is broken by restricting
`Omega` to `[0, 180)` while `T0` spans one full period — together exactly one fundamental
domain. An elliptical Gaussian's `pa` is likewise capped at 180°, since the profile is
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
                   method::Symbol = :neldermead, model::Symbol = :analytic,
                   weights = OI_DEFAULT_WEIGHTS, maxeval::Int = 20_000,
                   min_num_live_points::Int = 400, use_stepsampler::Bool = true,
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
                          free = free, bounds = bounds, ties = ties)
    idx  = spec.free
    θ0   = copy(spec.values)

    chi2_of = let spec = spec, fd = fd, weights = weights
        p -> _orbit_chi2(spec, resolve_params(spec, p), fd; weights = weights)
    end

    if verbose
        @printf("fit_orbit: %d epochs, %d points (V² %d, T3amp %d, T3φ %d)\n",
                length(fd.data), fd.ntot, fd.nv2, fd.nt3amp, fd.nt3phi)
        @printf("           %d free of %d parameters: %s\n", length(idx), length(spec.names),
                join(string.(spec.names[idx]), ", "))
        @printf("           starting χ²/n = %.3f\n", chi2_of(θ0[idx]) / fd.ntot)
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
        r = _fit_nested(method, chi2_of, string.(spec.names[idx]), spec.lo[idx], spec.hi[idx];
                        min_num_live_points = min_num_live_points,
                        use_stepsampler = use_stepsampler, verb = verbose)
        best  = r.median
        extra = (posterior = r.samples, logz = r.logz, logzerr = r.logzerr,
                 q16 = r.q16, q84 = r.q84)
    else
        error("fit_orbit: method must be :neldermead, :ultranest or :nautilus " *
              "(got $method)")
    end

    θbest = resolve_params(spec, best)
    cs    = _orbit_chi2(spec, θbest, fd; weights = weights, split = true)
    el    = NamedTuple{ORBIT_ELEMENTS}(ntuple(j -> θbest[j], length(ORBIT_ELEMENTS)))
    verbose && @printf("           final    χ²/n = %.3f  (V² %.2f, T3amp %.2f, T3φ %.2f)\n",
                       sum(cs)/fd.ntot, cs[1]/max(fd.nv2,1), cs[2]/max(fd.nt3amp,1),
                       cs[3]/max(fd.nt3phi,1))
    return merge((params = θbest, names = spec.names, elements = el,
                  chi2 = sum(cs), chi2_split = cs, chi2_red = sum(cs)/fd.ntot,
                  spec = spec, data = fd), extra)
end
