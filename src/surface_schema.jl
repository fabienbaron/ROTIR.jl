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
                        "T ∝ g^β, and under the ELR law T ∝ (F·g)^β. 0.25 radiative, " *
                        "~0.08 convective; set it fixed at 0.25 for ELR as published.")

# WHICH gravity-darkening law, for the surface types that have a choice. Carried as an integer
# code the way `ldtype` is, so a form renders it as a combo with no new machinery; the fit
# functions take either the code or the symbol. β stays free under both laws — pinning it at
# 1/4 to recover the published ELR result exactly is the ordinary free/fixed control.
const _GRAVITY_LAW = ParamSpec(:gravity_law, "Gravity law", "", 1.0, 1.0, 2.0, :thermal,
                               "1 von Zeipel, T ∝ g^β — derived for a barotropic star, so " *
                               "the slow-rotation law; it overestimates the pole-to-equator " *
                               "contrast for a fast rotator. 2 Espinosa Lara & Rieutord " *
                               "(2011), which adds their latitudinal flux factor and is the " *
                               "one to use for a fast rotator. Identical as the rotation " *
                               "goes to zero.";
                               kind = :choice,
                               choices = [1 => "von Zeipel", 2 => "Espinosa Lara-Rieutord"])

const _LIMBDARK = [
    ParamSpec(:ldtype, "LD law", "", 3.0, 0.0, 4.0, :limbdark,
              "Which law `compute_ldmap` applies: 0 none; 1 linear, 1 − u(1−μ); " *
              "2 quadratic, 1 − a(1−μ) − b(1−μ)²; 3 Hestroffer power law, μ^α; " *
              "4 Claret four-parameter. Anything outside 0–4 silently returns nothing " *
              "and fails later in the visibility model. " *
              "0 is for a MODEL-ATMOSPHERE intensity, which already carries the full μ " *
              "dependence — applying a law on top of it double-counts limb darkening. It " *
              "also removes ld1/ld2 from the fit, which is the point: gravity and limb " *
              "darkening both take flux out of the limb, an interferometer measures only " *
              "their sum, and a fitted ld1 absorbs the error in whichever " *
              "gravity-darkening law was assumed.";
              kind = :choice,
              choices = [0 => "none (from atmosphere)",
                         1 => "linear",
                         2 => "quadratic",
                         3 => "power law",
                         4 => "Claret-4"]),
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

# The distance, shared by every surface type that needs a physical scale. A factory rather
# than a constant because the two consumers file it under different form groups — for the
# Roche surface it belongs with the orbit, for a rapid rotator it IS a geometry parameter,
# since it is what converts the angular `rpole` into a length. Label, unit and bounds stay
# in one place either way.
_distance_spec(group::Symbol, doc::String; default::Float64 = 100.0) =
    ParamSpec(:d, "Distance", "pc", default, 0.1, 1e6, group, doc)

const _ORBIT_OPTIONAL = [
    ParamSpec(:dP, "Ṗ, period rate", "d/d", 0.0, -1.0, 1.0, :orbit,
              "Quadratic ephemeris. 0 for a constant period."),
    ParamSpec(:dω, "ω̇, apsidal", "deg/d", 0.0, -1.0, 1.0, :orbit,
              "Unicode dω. 0 for a fixed apsidal line."),
    _distance_spec(:orbit,
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
              # 0.95, not 0.5. A rapid rotator is worth looking at when it is visibly
              # oblate with real gravity darkening, and at 0.5 it is very nearly a sphere —
              # so "+ model" on this surface type showed something that did not look like
              # what it is for. 0.999 remains the upper bound; break-up is at 1.
              ParamSpec(:frac_escapevel, "v_eq / v_crit", "", 0.95, 0.0, 0.999, :geometry,
                        "Equatorial rotation as a fraction of critical. 1 is break-up: " *
                        "the equatorial radius diverges as it is approached.")],
             _THERMAL, [_BETA, _GRAVITY_LAW], _LIMBDARK, _ORIENTATION),
        [ParamSpec(:B_rot, "Diff. rotation B", "", 0.0, -1.0, 1.0, :orientation,
                   "Surface differential rotation: Ω(θ) = Ω₀(1 − B_rot·cos²θ), so " *
                   "`rotation_period` is the EQUATORIAL period and B_rot > 0 makes the " *
                   "poles lag. Read by `los_velocity` (src/velocity_field.jl) and so by " *
                   "every velocity-resolved observable; still NOT read by `rotate_star`, " *
                   "which spins the surface rigidly for spot phasing. Note a " *
                   "differentially rotating surface has no rotational potential, so its " *
                   "SHAPE does not follow the Roche factor — B_rot is a kinematic " *
                   "perturbation on a figure that is not self-consistent with it."),
         # Label kept short: the form's column fits `label + unit + 3 <= 20` characters and
         # elides silently past it, which `test/gui/runtests.jl:274` exists to catch.
         ParamSpec(:vgamma, "γ, systemic", "km/s", 0.0, -2000.0, 2000.0,
                   :orientation,
                   "Adds to every tessel's line-of-sight velocity, positive RECEDING. " *
                   "Shifts a line profile bodily without changing its shape, so it is what " *
                   "a spectroscopic dataset constrains and the interferometric " *
                   "observables (V², T3) are blind to. Read by `los_velocity`."),
         # 16 pc, NOT the 100 pc the Roche surface uses, and the value is COUPLED to the
         # other defaults on this surface type. Adding `d` gave `rotation_period` physical
         # meaning it did not have before — it used to only set a rotation phase — so the
         # defaults now have to be jointly possible. At 100 pc, rpole = 1 mas with
         # fev = 0.95 and P = 1 d implies a 499 Msun star, because M scales as d³; at 16 pc
         # it is 2.04 Msun, an ordinary F-type rapid rotator, and the whole default set
         # describes a real object (compare beta Cas: 0.849 mas, 16.8 pc, 0.92, 0.893 d).
         # `test_gravity_darkening.jl:68` asserts these defaults validate cleanly, which is
         # what caught the inconsistency.
         _distance_spec(:geometry,
                   "Converts the angular `rpole` into a length, which is what makes MASS " *
                   "and absolute `logg` available — see src/stellar_physics.jl. Mass is " *
                   "DERIVED, not fitted: (rpole, d, fev, rotation_period) are four " *
                   "quantities and Ω = fev·√(8GM/27R_p³) = 2π/P is one relation, so " *
                   "M = 27R_p³Ω²/(8G·fev²). That also ties `rotation_period` to " *
                   "`frac_escapevel`, which were previously free to disagree about how " *
                   "fast the star turns. Optional: a model using :linear or :planck " *
                   "intensity never needs it. The 16 pc default is chosen so that it, " *
                   "`rpole`, `frac_escapevel` and `rotation_period` together imply a " *
                   "2.04 Msun star rather than an impossible one."; default = 16.0)],
        "Roche-model oblate rotator, gravity-darkened by von Zeipel or by Espinosa Lara & "  *
        "Rieutord (2011) — see `gravity_law`."),

    3 => SurfaceSpec(3, :roche, "Roche lobe",
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
        # A `:choice` field may be named rather than numbered — `gravity_law = :elr` as well
        # as `gravity_law = 2` — because a bare integer code says nothing to a reader.
        if p.kind === :choice && (v isa Symbol || v isa AbstractString)
            p.name === :gravity_law ||
                throw(ArgumentError("$(p.name) takes a number, not a name; the choices are " *
                                    join(("$(k) = $(vv)" for (k, vv) in p.choices), ", ")))
            v = gravity_law_code(v)
        end
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
    if hasproperty(p, :ldtype) && !(Int(p.ldtype) in (0, 1, 2, 3, 4))
        push!(msgs, "ldtype = $(p.ldtype) is not one of 0 (none), 1 (linear), " *
                    "2 (quadratic), 3 (power law), 4 (Claret); `compute_ldmap` returns " *
                    "nothing for it and the failure surfaces much later, in the " *
                    "visibility model")
    end
    if hasproperty(p, :gravity_law)
        try
            gravity_law_spec(p.gravity_law)
        catch e
            push!(msgs, sprint(showerror, e))
        end
    end
    # NUMERIC SANITY over every declared field of this surface type. BLOCKING, and driving the
    # GUI's build gate is the whole point: `build_epoch_star` does
    # `isempty(validate_star_params(p)) || return nothing`, so this is the only thing between a
    # half-typed number and a degenerate mesh.
    #
    # THE CASE THIS EXISTS FOR, measured rather than imagined. `rpole = 0` on a Roche binary
    # takes the equipotential through a zero radius and the temperature map comes back
    # ENTIRELY NaN — and nothing downstream raises. `_map_range`'s `minimum` PROPAGATES NaN
    # instead of erroring, its `pmax - pmin < 1.0` widening never fires because `NaN < 1.0` is
    # false, and the NaN lands in a Makie Colorbar whose tick machinery is not total either:
    # PlotUtils reports "No strict ticks found" and Ryu then throws
    # `InexactError: convert(UInt64, …)` out of `writefixed`. That happens inside a QML
    # callback, so the window freezes with nothing printed. The range half of this loop is
    # what used to prevent it; `test/gui/runtests.jl`'s shape fuzz is what noticed when it
    # stopped. The non-finite half is new, because the range half SKIPS non-finite values (a
    # NaN fails every comparison, so `!(lo <= v <= hi)` would report it as "outside the
    # plausible range", which is true but useless) and `rpole = NaN` reaches the same crash.
    for spec in Iterators.flatten((s.required, s.optional))
        hasproperty(p, spec.name) || continue
        spec.kind === :choice && continue
        v = getproperty(p, spec.name)
        v isa Number || continue
        if v isa AbstractFloat && !isfinite(v)
            push!(msgs, "`$(spec.name)` is $v, which cannot describe a surface")
        elseif !(spec.lo <= v <= spec.hi)
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

"A field that is present, numeric, finite and strictly positive — the four things the derived
mass needs before it can be computed at all."
_finite_pos(p, f) = hasproperty(p, f) && (v = getproperty(p, f); v isa Number &&
                                          isfinite(float(v)) && float(v) > 0)

"""
    advise_star_params(p) -> Vector{String}

Things worth telling the user about `p` that do NOT make the model unbuildable.

SEPARATE FROM [`validate_star_params`](@ref) FOR A CONCRETE REASON: the GUI uses that one as
a hard GATE — `isempty(validate_star_params(p)) || return nothing` guards `epoch_chi2` and the
plot paths — so any message added there silently disables functionality rather than informing
anyone. Both checks below are observations about plausibility, not validity: the geometry
builds fine either way, and a fitter is entitled to walk through implausible values on its way
somewhere sensible.

A caller that wants everything a user should see joins the two lists; a caller deciding
whether the model can be BUILT consults `validate_star_params` alone.
"""
function advise_star_params(p)
    msgs = String[]
    # TOTAL BY CONSTRUCTION. `shell_validate_model` calls this on EVERY keystroke, from inside
    # `shell_set_param`, and the GUI test suite feeds those callbacks every shape QML can send —
    # empty strings, "NaN", huge numbers, a surface_type that is not a surface type. A throw
    # here propagates out of `shell_set_param` and the form dies on a typo. `validate_star_params`
    # is total for the same reason and guards its `surface_spec` call the same way. Both of the
    # guards below were added after an unguarded `surface_spec` call here threw out of 46 of
    # those callback tests.
    hasproperty(p, :surface_type) || return msgs
    stype = try Int(p.surface_type) catch; return msgs end
    haskey(SURFACE_TYPES, stype) || return msgs
    # LDTYPE = 0 WITHOUT AN ATMOSPHERE. The zero branch of `compute_ldmap` returns ones, which
    # is correct only when an `IntensityProvider` supplies the μ dependence instead. On its own
    # it means the star has NO limb darkening — and that does not raise anywhere, it just fits
    # a worse model. The GUI offers `ldtype = 0` in a combo box, so this message is the only
    # thing between a click and a silently wrong model. It cannot see the provider (that is a
    # fit argument, not a model field), so it always fires for ldtype = 0;
    # `check_provider_consistency` is the complement, firing when a provider IS present and
    # ldtype is NOT 0.
    ldt = hasproperty(p, :ldtype) ? (try Int(p.ldtype) catch; nothing end) : nothing
    if ldt == 0
        push!(msgs, "ldtype = 0 applies NO limb darkening. That is intended only when an " *
                    "`IntensityProvider` supplies the μ dependence instead — see " *
                    "`check_provider_consistency` and src/intensity_provider.jl. Without one, " *
                    "the model has a limb as bright as the disc centre.")
    end
    # THE DERIVED MASS. For a rapid rotator carrying a distance, (rpole, d, fev,
    # rotation_period) over-determine the mass — see src/stellar_physics.jl. An implausible
    # value is the visible symptom of `rotation_period` disagreeing with `frac_escapevel`,
    # which nothing used to catch because the two were read by different functions and never
    # compared.
    if stype == 2 && has_physical_scale(p) && _finite_pos(p, :frac_escapevel) &&
       _finite_pos(p, :rpole) && _finite_pos(p, :d) && _finite_pos(p, :rotation_period)
        m = derive_mass(float(p.rpole), float(p.d), float(p.frac_escapevel),
                        float(p.rotation_period))
        # `isfinite` first: a zero period or radius gives Inf/NaN, and a NaN fails every
        # comparison so the `!(lo <= m <= hi)` form would report it as out of range with a
        # "NaN Msun" message. Saying nothing is better than saying that.
        if isfinite(m) && !(0.05 <= m <= 200.0)
            push!(msgs, "the mass implied by rpole = $(p.rpole) mas, d = $(p.d) pc, " *
                        "frac_escapevel = $(p.frac_escapevel) and rotation_period = " *
                        "$(p.rotation_period) d is $(round(m, sigdigits=3)) Msun, outside " *
                        "[0.05, 200]; for a Roche rotator these four are not independent " *
                        "(M = 27R_p³Ω²/(8G·fev²)), so this usually means rotation_period " *
                        "and frac_escapevel disagree about how fast the star turns")
        end
    end
    # NO RANGE PASS HERE. It lives in `validate_star_params`, which is where it has always
    # lived and where it has to stay: the GUI gates the build on that list, and a value outside
    # a `ParamSpec`'s declared bounds is how a degenerate mesh gets built (see the comment on
    # that loop). Moving it here to stop the DERIVED-MASS message below from disabling
    # `epoch_chi2` took the geometry guard with it, and `rpole = 0` went back to freezing the
    # window. Only the mass check needed to move; it is the one that fires on a model that
    # builds perfectly well.
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
