# What a `star_params` NamedTuple has to contain, per surface type.
#
# `star_params` is an untyped NamedTuple built by hand in every demo script. A misspelt or
# missing field does not fail there: it fails several calls deep, inside `compute_radii` or
# `compute_ldmap`, as a `type NamedTuple has no field radius_x` with no indication of which
# surface type wanted it or what it meant. Worse, some omissions do not raise at all —
# `compute_ldmap` has branches for `ldtype` 1, 2 and 3 and falls off the end returning
# `nothing` for anything else, so an `ldtype = 0` produces a `MethodError` on `nothing` much
# later, in the visibility model.
#
# This table is the single declaration of the answer. It is ground truth taken from the code
# that reads the fields — `compute_radii` (src/geometry.jl:201), `compute_ldmap`
# (src/geometry.jl:250), `rotate_star` (src/geometry.jl:229), the `temperature_map_vonZeipel_*`
# family and, for Roche, `compute_separation` / `update_roche_radii` — not from the docs, which
# disagree with it (docs/src/guides/surfaces.md spells the orbital elements in ASCII; the code
# reads `Ω`, `ω`, `dω`, and ASCII names fail at runtime).
#
# Two consumers: `validate_star_params` turns a wrong NamedTuple into one readable message
# before any geometry is built, and a GUI generates its parameter form — labels, units,
# defaults, bounds, which fields even apply — from the same table rather than hardcoding a
# second copy of it that can drift.

"""
    ParamSpec

One field of a `star_params` NamedTuple: what it is called, what it means, and what a
plausible value looks like.

`lo`/`hi` are *plausibility* bounds for a fit or a form, not hard validity limits — the
validator warns outside them and refuses nothing. `kind` is `:float`, `:int` or `:choice`;
`choices` is non-empty only for `:choice`, mapping each allowed value to its meaning.
"""
struct ParamSpec
    name::Symbol
    label::String
    unit::String
    default::Float64
    lo::Float64
    hi::Float64
    group::Symbol      # :geometry, :thermal, :limbdark, :orientation, :orbit
    kind::Symbol       # :float, :int, :choice
    choices::Vector{Pair{Int,String}}
    doc::String
end

ParamSpec(name, label, unit, default, lo, hi, group, doc;
          kind::Symbol = :float, choices = Pair{Int,String}[]) =
    ParamSpec(name, label, unit, default, lo, hi, group, kind,
              collect(Pair{Int,String}, choices), doc)

"""
    SurfaceSpec

Everything a given `surface_type` needs, split into what must be present and what only some
code paths read.

`required` is exactly the set whose absence raises. `optional` is read on some paths and
defaulted on others — `dP`/`dω` are zero for a fixed ephemeris, `B_rot` is carried by
`starparameters` but the differential-rotation code that used it is not currently wired in.
"""
struct SurfaceSpec
    code::Int
    name::Symbol
    label::String
    required::Vector{ParamSpec}
    optional::Vector{ParamSpec}
    doc::String
end

# ── fields shared by every surface type ─────────────────────────────────────────────────
#
# `rotation_period` is required even for a sphere: `rotate_star` divides by it unconditionally,
# so a missing or zero value is a division by zero rather than "no rotation".
const _ORIENTATION = [
    ParamSpec(:inclination, "Inclination", "deg", 60.0, 0.0, 180.0, :orientation,
              "Spin-axis inclination. >90 deg views the retrograde pole."),
    ParamSpec(:position_angle, "Position angle", "deg", 0.0, -180.0, 360.0, :orientation,
              "Spin-axis PA on the sky, North through East."),
    ParamSpec(:rotation_period, "Rotation period", "d", 1.0, 1e-6, 1e5, :orientation,
              "Divides the epoch to give the rotation phase; never 0."),
]

const _THERMAL = [
    ParamSpec(:tpole, "Polar temp.", "K", 6000.0, 1000.0, 60000.0, :thermal,
              "Effective temperature at the pole, before gravity darkening."),
]

# β enters every von Zeipel map. 0.25 is the radiative value, ~0.08 convective (Lucy 1967).
const _BETA = ParamSpec(:beta, "Grav. darkening β", "", 0.25, 0.0, 0.5, :thermal,
                        "T ∝ g^β. 0.25 radiative, ~0.08 convective.")

const _LIMBDARK = [
    ParamSpec(:ldtype, "LD law", "", 3.0, 1.0, 4.0, :limbdark,
              "Which law `compute_ldmap` applies. Anything outside 1–4 silently " *
              "returns nothing and fails later in the visibility model.";
              kind = :choice,
              choices = [1 => "linear: 1 − u(1−μ)",
                         2 => "quadratic: 1 − a(1−μ) − b(1−μ)²",
                         3 => "Hestroffer power law: μ^α",
                         4 => "Claret 4-parameter"]),
    ParamSpec(:ld1, "LD coeff 1", "", 0.2, -1.0, 2.0, :limbdark,
              "u (linear), a (quadratic), α (power law) or a₁ (Claret)."),
    ParamSpec(:ld2, "LD coeff 2", "", 0.0, -1.0, 1.0, :limbdark,
              "b of the quadratic law, or a₂ of Claret's."),
    ParamSpec(:ld3, "LD coeff 3", "", 0.0, -1.0, 1.0, :limbdark,
              "a₃ of Claret's four-parameter law; read only when ldtype = 4."),
    ParamSpec(:ld4, "LD coeff 4", "", 0.0, -1.0, 1.0, :limbdark,
              "a₄ of Claret's four-parameter law; read only when ldtype = 4."),
]

# ── orbital elements, for the Roche surface ─────────────────────────────────────────────
#
# UNICODE, and that is not cosmetic: `compute_coeff`, `omega_at` and `binary_orbit_abs` read
# `bparams.Ω`, `.ω` and `.dω` by those exact names. An ASCII `Omega` is a different field and
# the NamedTuple lookup fails at runtime.
const _ORBIT = [
    ParamSpec(:P, "Orbital period", "d", 10.0, 1e-4, 1e6, :orbit,
              "Sets both the separation history and, when tidally locked, the spin."),
    ParamSpec(:a, "Semi-major", "mas", 3.0, 1e-4, 1e4, :orbit,
              "Of the RELATIVE orbit. The default is 6x the default `rpole`, which keeps " *
              "the default Roche surface well inside its lobe — at a ~ rpole the star " *
              "strains against L1 and gravity darkening spans tens of thousands of K."),
    ParamSpec(:e, "Eccentricity", "", 0.0, 0.0, 0.99, :orbit, ""),
    ParamSpec(:T0, "Periastron", "JD", 2450000.0, 0.0, 1e7, :orbit, ""),
    ParamSpec(:i, "Orbital incl.", "deg", 90.0, 0.0, 180.0, :orbit,
              "Distinct from the component's own spin `inclination`."),
    ParamSpec(:Ω, "Ω, asc. node", "deg", 0.0, -180.0, 360.0, :orbit,
              "Unicode Ω — an ASCII `Omega` is a different field and fails at runtime."),
    ParamSpec(:ω, "ω, periapsis", "deg", 0.0, -180.0, 360.0, :orbit,
              "Of the RELATIVE orbit (secondary about primary). Unicode ω."),
    ParamSpec(:q, "Mass ratio q", "", 0.5, 1e-3, 1e3, :orbit,
              "M_companion/M_self for the Roche potential: q for the primary, 1/q " *
              "for the secondary. NOT the same convention on both components."),
]

const _ORBIT_OPTIONAL = [
    ParamSpec(:dP, "Ṗ, period rate", "d/d", 0.0, -1.0, 1.0, :orbit,
              "Quadratic ephemeris. 0 for a constant period."),
    ParamSpec(:dω, "ω̇, apsidal", "deg/d", 0.0, -1.0, 1.0, :orbit,
              "Unicode dω. 0 for a fixed apsidal line."),
    ParamSpec(:d, "Distance", "pc", 100.0, 0.1, 1e6, :orbit,
              "Carried for physical-unit conversions; the geometry is in mas throughout."),
]

# `rotate_star` spins the component at `rotation_period` regardless of the orbit, so for a
# Roche surface these two defaults have to agree: a component turning 10x faster than it
# orbits is a near-break-up rotator, and the resulting von Zeipel map spans tens of thousands
# of K purely because two defaults disagreed. Tidal locking is the normal case and the only
# self-consistent default; a non-synchronous rotator is a deliberate choice the caller makes.
const _ORIENTATION_SYNC = [
    _ORIENTATION[1],
    _ORIENTATION[2],
    ParamSpec(:rotation_period, "Rotation period", "d", 10.0, 1e-6, 1e5, :orientation,
              "Tidally locked components have rotation_period == P; that is the default " *
              "here. Setting it away from P models a non-synchronous rotator."),
]

const _FILLOUT = [
    ParamSpec(:fillout_factor_primary, "Fill-out, primary", "", -1.0, -1.0, 1.5, :geometry,
              "Roche-lobe fill-out. -1 means: use `rpole` to define the equipotential."),
    ParamSpec(:fillout_factor_secondary, "Fill-out, secondary", "", -1.0, -1.0, 1.5, :geometry,
              "As above, for the secondary."),
]

"""
    SURFACE_TYPES

`surface_type` code → [`SurfaceSpec`](@ref). The codes are the integers `compute_radii`
branches on, and nothing else is implemented: a `surface_type` outside 0–3 falls through that
`if`/`elseif` chain leaving `xyz` and `r` as empty `Vector{Any}`, which then fails in
`finish_star` with an unrelated-looking error.
"""
const SURFACE_TYPES = Dict{Int,SurfaceSpec}(
    0 => SurfaceSpec(0, :sphere, "Sphere",
        vcat([ParamSpec(:radius, "Radius", "mas", 1.0, 1e-4, 1e3, :geometry,
                        "Uniform radius. Note this type uses `radius`, NOT `rpole`.")],
             _THERMAL, _LIMBDARK, _ORIENTATION),
        ParamSpec[],
        "Uniform sphere at `tpole`; no gravity darkening, so `beta` is not read."),

    1 => SurfaceSpec(1, :ellipsoid, "Triaxial ellipsoid",
        vcat([ParamSpec(:radius_x, "Radius x", "mas", 1.2, 1e-4, 1e3, :geometry,
                        "Body-frame semi-axis along x; also the normalisation of the " *
                        "ellipsoid von Zeipel map."),
              ParamSpec(:radius_y, "Radius y", "mas", 1.0, 1e-4, 1e3, :geometry, ""),
              ParamSpec(:radius_z, "Radius z", "mas", 1.0, 1e-4, 1e3, :geometry, "")],
             _THERMAL, [_BETA], _LIMBDARK, _ORIENTATION),
        ParamSpec[],
        "Triaxial ellipsoid with a von Zeipel map normalised on `radius_x`."),

    2 => SurfaceSpec(2, :rapid_rotator, "Rapid rotator",
        vcat([ParamSpec(:rpole, "Polar radius", "mas", 1.0, 1e-4, 1e3, :geometry, ""),
              ParamSpec(:frac_escapevel, "v_eq / v_crit", "", 0.5, 0.0, 0.999, :geometry,
                        "Equatorial rotation as a fraction of critical. 1 is break-up: " *
                        "the equatorial radius diverges as it is approached.")],
             _THERMAL, [_BETA], _LIMBDARK, _ORIENTATION),
        [ParamSpec(:B_rot, "Diff. rotation B", "", 0.0, -1.0, 1.0, :orientation,
                   "Carried by `starparameters` but NOT currently read: the " *
                   "differential-rotation path is not wired in (see rotate_star).")],
        "Roche-model oblate rotator with von Zeipel gravity darkening."),

    3 => SurfaceSpec(3, :roche, "Roche lobe (binary component)",
        vcat([ParamSpec(:rpole, "Polar radius", "mas", 0.5, 1e-4, 1e3, :geometry,
                        "Defines the equipotential when the matching fill-out factor is -1.")],
             _FILLOUT, _THERMAL, [_BETA], _LIMBDARK, _ORIENTATION_SYNC, _ORBIT),
        _ORBIT_OPTIONAL,
        "One component of a binary. The shape depends on the INSTANTANEOUS separation, so " *
        "the full orbit is part of the surface definition, not merely of where it is drawn."),
)

"""
    SURFACE_TYPE_ORDER

Codes in the order a selector should list them: increasing complexity, which is also
increasing `surface_type`.
"""
const SURFACE_TYPE_ORDER = [0, 1, 2, 3]

"""
    surface_spec(x) -> SurfaceSpec

Look up a surface type by code (`3`), by name (`:roche`), or from anything carrying a
`surface_type` field — a `star_params` NamedTuple or a `stellar_geometry`.
"""
surface_spec(code::Integer) = get(SURFACE_TYPES, Int(code)) do
    throw(ArgumentError("unknown surface_type $code; implemented: " *
                        join(("$c ($(SURFACE_TYPES[c].name))" for c in SURFACE_TYPE_ORDER), ", ")))
end
function surface_spec(name::Symbol)
    for c in SURFACE_TYPE_ORDER
        SURFACE_TYPES[c].name === name && return SURFACE_TYPES[c]
    end
    throw(ArgumentError("unknown surface name :$name; implemented: " *
                        join((":$(SURFACE_TYPES[c].name)" for c in SURFACE_TYPE_ORDER), ", ")))
end
surface_spec(x) = surface_spec(x.surface_type)

"""
    surface_params(x; optional=true) -> Vector{ParamSpec}

Every field the given surface type reads, required first. `optional=false` gives only the
fields whose absence raises.
"""
function surface_params(x; optional::Bool = true)
    s = surface_spec(x)
    return optional ? vcat(s.required, s.optional) : copy(s.required)
end

"""
    default_star_params(x; T=Float64, kwargs...) -> NamedTuple

A complete, valid `star_params` for a surface type, with every field at its schema default,
overridden by whatever is passed as a keyword.

`surface_type` and `ldtype` stay `Int` — `compute_radii` and `compute_ldmap` branch on them
with `==`, and a Float64 3.0 compares equal but reads as a continuous parameter to anything
that iterates the schema. Everything else is `T`, so a Float32 model stays Float32 rather than
being silently widened by the defaults.

Unknown keywords are an error, not a silent addition: passing `Omega = 45` instead of `Ω = 45`
is exactly the mistake this table exists to catch, and adding it as a new field would let it
through.
"""
function default_star_params(x; T::Type = Float64, kwargs...)
    s = surface_spec(x)
    specs = vcat(s.required, s.optional)
    known = Set(p.name for p in specs)
    unknown = setdiff(keys(kwargs), known)
    if !isempty(unknown)
        msg = "not fields of surface_type $(s.code) ($(s.name)): " *
              join(sort(collect(String.(unknown))), ", ")
        # The one near-miss worth naming, because it is the near-miss that actually happens.
        any(n in (:Omega, :omega, :domega, :dOmega, :Omega_dot) for n in unknown) &&
            (msg *= " — the orbital elements are Unicode: Ω, ω, dω")
        throw(ArgumentError(msg))
    end
    pairs_ = Pair{Symbol,Any}[:surface_type => s.code]
    for p in specs
        v = get(kwargs, p.name, p.default)
        push!(pairs_, p.name => p.kind === :float ? T(v) : Int(v))
    end
    return NamedTuple(pairs_)
end

"""
    validate_star_params(p) -> Vector{String}

Every problem with a `star_params`, as messages; empty means it will build.

Reports, in order: an unimplemented `surface_type`; missing required fields; an `ldtype`
outside the implemented laws; values outside their plausibility bounds. The bound messages are
warnings — a real star can sit outside a range this table calls plausible — but a missing
field or a bad `ldtype` will raise once the geometry is built, which is what this is for.

Deliberately does not throw. The GUI marks fields; a script decides for itself:

    issues = validate_star_params(p); isempty(issues) || error(join(issues, "\\n"))
"""
function validate_star_params(p)
    msgs = String[]
    hasproperty(p, :surface_type) ||
        return ["missing `surface_type`: cannot tell which fields are required"]
    if !haskey(SURFACE_TYPES, Int(p.surface_type))
        return ["surface_type $(p.surface_type) is not implemented (0 sphere, " *
                "1 ellipsoid, 2 rapid rotator, 3 Roche); `compute_radii` would leave the " *
                "geometry empty and fail later in `finish_star`"]
    end
    s = surface_spec(p)
    for spec in s.required
        if !hasproperty(p, spec.name)
            push!(msgs, "missing `$(spec.name)` ($(spec.label)), required by surface_type " *
                        "$(s.code) ($(s.name))")
        end
    end
    if hasproperty(p, :ldtype) && !(Int(p.ldtype) in (1, 2, 3, 4))
        push!(msgs, "ldtype = $(p.ldtype) is not one of 1 (linear), 2 (quadratic), " *
                    "3 (power law), 4 (Claret); `compute_ldmap` returns nothing for it and " *
                    "the failure surfaces much later, in the visibility model")
    end
    for spec in vcat(s.required, s.optional)
        hasproperty(p, spec.name) || continue
        spec.kind === :choice && continue
        v = getproperty(p, spec.name)
        v isa Number || continue
        if !(spec.lo <= v <= spec.hi)
            push!(msgs, "$(spec.name) = $v is outside the plausible range " *
                        "[$(spec.lo), $(spec.hi)]$(isempty(spec.unit) ? "" : " " * spec.unit)")
        end
    end
    # Fields the caller probably meant to spell in Unicode. Cheap to check, and it is the
    # single most common way a hand-written Roche parameter set fails.
    for (ascii, uni) in (("Omega", :Ω), ("omega", :ω), ("domega", :dω))
        if hasproperty(p, Symbol(ascii)) && !hasproperty(p, uni)
            push!(msgs, "`$ascii` is present but `$uni` is not — the orbit code reads the " *
                        "Unicode name and will not see this field")
        end
    end
    return msgs
end

"""
    ld_coefficients_used(ldtype) -> Tuple{Symbol,...}

Which of `ld1`/`ld2` the given law actually reads. A form greys out the rest; a fitter must
not let a parameter that is never read float free, since it is then perfectly unconstrained.
"""
ld_coefficients_used(ldtype::Integer) =
    Int(ldtype) == 4 ? (:ld1, :ld2, :ld3, :ld4) :
    Int(ldtype) == 2 ? (:ld1, :ld2) :
    Int(ldtype) in (1, 3) ? (:ld1,) : ()
