# The Orbit perspective: a binary's elements, a fit, and a picture of the orbit.
#
# TWO FRAMES, KEPT APART, which is the thing this tab gets wrong if it is written casually.
#
#   * The **orbital elements** — a, i, Ω, ω, e, P, T0, dP, dω — describe where the secondary
#     sits relative to the primary. They are the same whatever the stars are made of.
#   * The **star model** says what sits at those two positions: either ANALYTIC 2-D profiles
#     (uniform disk, limb-darkened disk, Gaussian) or ROTIR's own TESSELLATED 3-D components,
#     which is what makes this ROTIR rather than a generic astrometry tool — Roche shapes from
#     the instantaneous separation, gravity darkening, mutual irradiation, occultation.
#
# A parameter belongs to exactly one of the two, and the panel is laid out that way: elements
# above, star model below. `q`, `rpole` and `tpole` are star-model quantities and never
# elements — the relative astrometric orbit does not constrain a mass ratio at all.
#
# The 3-D path is fully wired for RENDERING (`_render_epochs` builds real
# `create_binary_geometry` surfaces). For FITTING, `fit_orbit(model = :tessellated)` still
# throws by design; selecting it here passes the choice straight through, so the refusal comes
# from `fit_orbit` itself with its own explanation rather than from a hidden guard here.

"""
    OrbitEntry

Everything the Orbit tab holds: the two components, the elements, what is free, and which
proximity effects the rendering includes.

`elements` covers `ORBIT_ELEMENTS` plus `f` (the flux ratio) and the component parameters,
which is exactly the parameter vector `orbit_fit_spec` assembles — kept as one flat dictionary
so the form, the bounds and the ties are all keyed the same way.
"""
mutable struct OrbitEntry
    name::String
    kind1::Symbol                       # :uniform, :limbdark, :gaussian, :ellgauss, :point
    kind2::Symbol
    params::Dict{Symbol,Float64}        # elements + :f + c1_*/c2_*
    free::Set{Symbol}
    bounds::Dict{Symbol,Tuple{Float64,Float64}}
    ties::Dict{Symbol,String}
    # The star-model frame: what kind of star sits at each orbital position, and the
    # parameters that kind needs. Separate from the elements above by construction.
    model::Symbol                       # :analytic (2-D profiles) or :tessellated (3-D)
    render::Bool                        # draw the star surfaces at each epoch
    roche::Bool                         # recompute the Roche shape from D(t) per epoch
    irradiation::Bool
    occultation::Bool
    # Tessellated rendering parameters, which the analytic fit never sees.
    rpole1::Float64
    rpole2::Float64
    tpole1::Float64
    tpole2::Float64
    q::Float64
end

"The component kinds the tab offers, and what each contributes to the parameter vector."
const ORBIT_COMPONENT_KINDS = [
    ("uniform",  "uniform disk",        [:diameter]),
    ("limbdark", "limb-darkened disk",  [:diameter, :ld1]),
    ("gaussian", "circular Gaussian",   [:fwhm]),
    ("ellgauss", "elliptical Gaussian", [:fwhm, :ratio, :pa]),
    ("point",    "point source",        Symbol[]),
]

_kind_params(k::Symbol) =
    something(findfirst(x -> Symbol(x[1]) == k, ORBIT_COMPONENT_KINDS), 0) == 0 ? Symbol[] :
    ORBIT_COMPONENT_KINDS[findfirst(x -> Symbol(x[1]) == k, ORBIT_COMPONENT_KINDS)][3]

"Build the `OrbitComponent` a kind and its parameter values describe."
function _orbit_component(kind::Symbol, p::AbstractDict, prefix::Symbol)
    g(n, d) = get(p, Symbol("$(prefix)_$(n)"), d)
    kind === :uniform  && return UniformDisk(diameter = g(:diameter, 1.0))
    kind === :limbdark && return LimbDarkenedDisk(diameter = g(:diameter, 1.0),
                                                  law = :linear, ld1 = g(:ld1, 0.3))
    kind === :gaussian && return GaussianDisk(fwhm = g(:fwhm, 1.0))
    kind === :ellgauss && return EllipticalGaussian(fwhm = g(:fwhm, 1.0),
                                                    ratio = g(:ratio, 1.0), pa = g(:pa, 0.0))
    return PointSource()
end

"Sensible starting values for a component's parameters."
function _seed_component!(p::Dict{Symbol,Float64}, kind::Symbol, prefix::Symbol)
    for n in _kind_params(kind)
        k = Symbol("$(prefix)_$(n)")
        haskey(p, k) || (p[k] = n === :diameter ? 0.5 : n === :fwhm ? 0.5 :
                                n === :ratio ? 1.0 : n === :ld1 ? 0.3 : 0.0)
    end
    return p
end

"""
    _analytic_radius(kind, params, prefix) -> mas

The outer radius the ANALYTIC profile of one component draws at.

The same extent `_profile_rings` uses, which is half the diameter for a disc and the FWHM
itself for a Gaussian — so a tessellated radius taken from here renders the component at the
size its 2-D profile has.
"""
_analytic_radius(kind::Symbol, p::AbstractDict, pre::Symbol) =
    kind in (:uniform, :limbdark) ? Float64(get(p, Symbol("$(pre)_diameter"), 0.5)) / 2 :
    kind in (:gaussian, :ellgauss) ? Float64(get(p, Symbol("$(pre)_fwhm"), 0.5)) :
    0.02

"""
    default_orbit() -> OrbitEntry

A two-disk orbit at neutral values. Deliberately not a copy of β Lyr or Spica: a default that
looks like a real system invites its numbers being left in place by accident.
"""
function default_orbit()
    p = Dict{Symbol,Float64}(:a => 1.0, :i => 90.0, :Omega => 0.0, :omega => 90.0,
                             :e => 0.0, :P => 10.0, :T0 => 2450000.0, :dP => 0.0,
                             :domega => 0.0, :f => 1.0)
    _seed_component!(p, :uniform, :c1); _seed_component!(p, :uniform, :c2)
    b = Dict{Symbol,Tuple{Float64,Float64}}()
    for (k, v) in ORBIT_ELEMENT_BOUNDS; b[k] = v; end
    b[:f] = (0.0, 10.0)
    for pre in (:c1, :c2), n in (:diameter, :fwhm)
        b[Symbol("$(pre)_$(n)")] = (0.01, 20.0)
    end
    b[:c1_ld1] = (0.0, 1.0); b[:c2_ld1] = (0.0, 1.0)
    b[:c1_ratio] = (0.05, 1.0); b[:c2_ratio] = (0.05, 1.0)
    b[:c1_pa] = (-180.0, 180.0); b[:c2_pa] = (-180.0, 180.0)
    # The tessellated polar radii are HALF THE DEFAULT DIAMETERS, so the two star models draw
    # the same star. They were 0.2 and 0.12 against two 0.5 mas discs, so switching from
    # "analytic 2-D profiles" to the 3-D surfaces changed the size of both components — which
    # reads as one of the two views having a factor wrong.
    return OrbitEntry("orbit", :uniform, :uniform, p, Set([:a, :i, :Omega, :T0]), b,
                      Dict{Symbol,String}(), :analytic, false, true, false, false,
                      _analytic_radius(:uniform, p, :c1), _analytic_radius(:uniform, p, :c2),
                      20000.0, 15000.0, 0.5)
end

"""
    orbit_param_names(o) -> Vector{Symbol}

The parameter rows the form shows, in the order they appear.

The nine ELEMENTS always; the flux ratio and the per-component profile parameters
(`c1_diameter`, `c2_fwhm`, …) only under the ANALYTIC star model, because that is the only
model that reads them. The tessellated path sizes its components from `rpole1`/`rpole2` and
takes their relative brightness from the temperatures, so `c1_diameter` shown beside a 3-D
binary is a number nothing uses — it sat in the same flat table as the elements and was drawn
whatever the star model was.

The values are not deleted, only unlisted: switching back to analytic brings them back as they
were, which is what makes the two star models comparable at all.
"""
function orbit_param_names(o::OrbitEntry; analytic::Bool = o.model === :analytic)
    ns = Symbol[collect(ORBIT_ELEMENTS)...]
    analytic || return ns
    push!(ns, :f)
    append!(ns, Symbol("c1_$(n)") for n in _kind_params(o.kind1))
    append!(ns, Symbol("c2_$(n)") for n in _kind_params(o.kind2))
    return ns
end

"Units and one-line meanings, for the form. Empty where the quantity is a ratio."
const ORBIT_PARAM_INFO = Dict(
    :a => ("mas", "semi-major axis of the RELATIVE orbit"),
    :i => ("deg", "inclination; >90 is retrograde"),
    :Omega => ("deg", "longitude of the ascending node"),
    :omega => ("deg", "argument of periapsis; undefined at e = 0"),
    :e => ("", "eccentricity"),
    :P => ("d", "orbital period"),
    :T0 => ("JD", "time of periastron passage"),
    :dP => ("d/d", "period derivative; 0 for a constant period"),
    :domega => ("deg/d", "apsidal motion"),
    :f => ("", "flux ratio, secondary / primary"))

_orbit_info(n::Symbol) = get(ORBIT_PARAM_INFO, n,
    (occursin("pa", String(n)) ? "deg" : occursin("ratio", String(n)) ? "" : "mas",
     "component parameter"))

"""
    orbit_elements_nt(o) -> NamedTuple

The nine orbital elements as `orbit_fit_spec` wants them, ties applied.
"""
function orbit_elements_nt(o::OrbitEntry)
    v = apply_orbit_ties(o)
    return NamedTuple(k => get(v, k, 0.0) for k in ORBIT_ELEMENTS)
end

"""
    apply_orbit_ties(o) -> Dict

Parameter values with every tie expression evaluated, to a fixed point — the same treatment
model ties get, and for the same reason: `Omega` and a component position angle are routinely
tied to each other, and the order they are declared in should not decide the answer.
"""
function apply_orbit_ties(o::OrbitEntry)
    vals = copy(o.params)
    isempty(o.ties) && return vals
    for _ in 1:(length(o.ties) + 1)
        changed = false
        for (name, expr) in o.ties
            v = eval_tie(expr, vals)
            v === nothing && continue
            if !haskey(vals, name) || vals[name] != v
                vals[name] = v; changed = true
            end
        end
        changed || return vals
    end
    @warn "orbit tie expressions did not settle; a cycle is likely"
    return vals
end

"""
    _orbit_binary(sh) -> ModelEntry or nothing

The Model tab's binary, when there is one: the PRIMARY entry, whose `companion` is the
secondary. `nothing` when there is no model, or the model is a single star.

One description of the system, shared between the two tabs. The Orbit tab says where the
components are; the Model tab says what they are made of.
"""
function _orbit_binary(sh)
    m = current_model(sh.session)
    (m === nothing || m.companion === nothing) && return nothing
    return m
end

"""
    orbit_bparams(o) -> NamedTuple

The `binaryparameters` shape the geometry and the plotting take, from the fitted elements.

`q` is NOT an orbital element — an astrometric orbit of the RELATIVE motion does not constrain
the mass ratio — so it comes from the rendering block, where it is a stated assumption rather
than something the fit produced.
"""
function orbit_bparams(o::OrbitEntry; binary = nothing)
    v = apply_orbit_ties(o)
    # BOTH components are oriented by the orbit, not by themselves: `create_binary_geometry`
    # takes the spin axis from the orbital inclination and node, so a per-component
    # `inclination` in the Model tab would be overwritten here anyway. Everything else —
    # radius, temperature, the limb-darkening law and its coefficients, β, the rotation period
    # — comes from the component when there is one.
    orient(sp) = (180.0 - v[:i], v[:Omega] - 180.0)
    function comp(m, rpole_fallback, tpole_fallback)
        m === nothing && return starparameters(rpole_fallback, tpole_fallback, 0.0, 3, 0.2,
                                               0.0, 0.25, 0.0, 180.0 - v[:i],
                                               v[:Omega] - 180.0, 0.0, v[:P])
        p = m.params
        g(k, d) = Float64(get(p, k, d))
        # `rpole` for a Roche or rapid-rotator component, `radius` for a sphere: the schema
        # gives them different names and a binary can legitimately hold either.
        r = haskey(p, :rpole) ? g(:rpole, rpole_fallback) : g(:radius, rpole_fallback)
        return starparameters(r, g(:tpole, tpole_fallback), g(:frac_escapevel, 0.0),
                              round(Int, g(:ldtype, 3.0)), g(:ld1, 0.2), g(:ld2, 0.0),
                              g(:beta, 0.25), g(:B_rot, 0.0),
                              180.0 - v[:i], v[:Omega] - 180.0, 0.0,
                              g(:rotation_period, v[:P]))
    end
    m1 = binary === nothing ? nothing : binary
    m2 = binary === nothing ? nothing : binary.companion
    s1 = comp(m1, o.rpole1, o.tpole1)
    s2 = comp(m2, o.rpole2, o.tpole2)
    # `q` follows the components when they are supplying everything else — the primary's Roche
    # potential carries it — and otherwise stays the rendering assumption it has always been.
    q = m1 === nothing ? o.q : Float64(get(m1.params, :q, o.q))
    return binaryparameters(s1, s2, 100.0, v[:i], v[:Omega], v[:omega], v[:P], v[:a],
                            v[:e], v[:T0], q, [1.0, 1.0], v[:dP], v[:domega])
end

# ── save and load ────────────────────────────────────────────────────────────────────────

"""
    save_orbit(o, path) -> path

Write the orbit as TOML: the elements, the bounds, the ties, the component kinds and the
rendering block.

TOML rather than a Julia snippet, because this file is DATA — it is edited by hand between
sessions, diffed against a published solution, and read back by a script that is not this GUI.
The command log is where the runnable version lives.
"""
function save_orbit(o::OrbitEntry, path::AbstractString)
    d = Dict{String,Any}(
        "name" => o.name,
        "components" => Dict("c1" => String(o.kind1), "c2" => String(o.kind2)),
        "parameters" => Dict(String(k) => v for (k, v) in o.params),
        "free" => sort(String.(collect(o.free))),
        "bounds" => Dict(String(k) => [v[1], v[2]] for (k, v) in o.bounds),
        "ties" => Dict(String(k) => v for (k, v) in o.ties),
        "star_model" => Dict("kind" => String(o.model)),
        "rendering" => Dict("render" => o.render, "roche" => o.roche,
                            "irradiation" => o.irradiation, "occultation" => o.occultation,
                            "rpole1" => o.rpole1, "rpole2" => o.rpole2,
                            "tpole1" => o.tpole1, "tpole2" => o.tpole2, "q" => o.q))
    mkpath(dirname(abspath(path)))
    open(io -> TOML.print(io, d), path, "w")
    return path
end

"""
    load_orbit(path) -> OrbitEntry

Read one back. Missing keys keep the default, so a file written by an older version — or one
trimmed by hand to the few elements that matter — still loads.
"""
function load_orbit(path::AbstractString)
    d = TOML.parsefile(String(path))
    o = default_orbit()
    o.name = get(d, "name", o.name)
    if haskey(d, "components")
        o.kind1 = Symbol(get(d["components"], "c1", String(o.kind1)))
        o.kind2 = Symbol(get(d["components"], "c2", String(o.kind2)))
    end
    _seed_component!(o.params, o.kind1, :c1); _seed_component!(o.params, o.kind2, :c2)
    for (k, v) in get(d, "parameters", Dict()); v isa Real && (o.params[Symbol(k)] = Float64(v)); end
    haskey(d, "free") && (o.free = Set(Symbol.(d["free"])))
    for (k, v) in get(d, "bounds", Dict())
        v isa AbstractVector && length(v) == 2 &&
            (o.bounds[Symbol(k)] = (Float64(v[1]), Float64(v[2])))
    end
    for (k, v) in get(d, "ties", Dict()); v isa AbstractString && (o.ties[Symbol(k)] = v); end
    haskey(d, "star_model") &&
        (o.model = Symbol(get(d["star_model"], "kind", String(o.model))))
    r = get(d, "rendering", Dict())
    # The tessellated radii default to the ANALYTIC sizes, so an orbit file that says nothing
    # about rendering still draws the same star under both star models. Both shipped orbits
    # carried `default_orbit`'s old 0.2 and 0.12 — never revised when the files were written —
    # and Spica's analytic primary is 0.894 mas across, so its 3-D surfaces came out at less
    # than half the size of its 2-D profiles. A file that DOES give a radius still wins: a
    # component can legitimately render at a size its analytic profile does not have.
    o.rpole1 = _analytic_radius(o.kind1, o.params, :c1)
    o.rpole2 = _analytic_radius(o.kind2, o.params, :c2)
    o.render      = get(r, "render", o.render)
    o.roche       = get(r, "roche", o.roche)
    o.irradiation = get(r, "irradiation", o.irradiation)
    o.occultation = get(r, "occultation", o.occultation)
    for f in (:rpole1, :rpole2, :tpole1, :tpole2, :q)
        haskey(r, String(f)) && setfield!(o, f, Float64(r[String(f)]))
    end
    return o
end

# ── the canvas ───────────────────────────────────────────────────────────────────────────

"""
    OrbitCanvas

The orbit on the sky: the relative track, the secondary's position at every observed epoch,
and — optionally — the two star surfaces rendered at each of those epochs.

The rendered surfaces are ONE `poly!` holding every tessel of every epoch, not one plot per
epoch. That is the plot-once rule: a plot inserted after the window exists allocates GPU
buffers with no context bound. It is also what makes the render toggle instant — turning it off
assigns an empty vector rather than tearing plots down.
"""
struct OrbitCanvas
    figure::Any
    axis::Any
    track::Makie.Observable{Vector{Makie.Point2f}}
    trackplot::Any
    marks::Makie.Observable{Vector{Makie.Point2f}}
    marksplot::Any
    labelpos::Makie.Observable{Vector{Makie.Point2f}}
    labeltext::Makie.Observable{Vector{String}}
    labelplot::Any
    polys::Makie.Observable{Vector{Vector{Makie.Point2f}}}
    polycolors::Makie.Observable{Vector{Makie.RGBAf}}
    polyplot::Any
    primary::Makie.Observable{Vector{Makie.Point2f}}
    primaryplot::Any
    # The synthetic TIME CURSOR: one position on the track, drawn filled so it cannot be read
    # as one more observed epoch. Empty unless the Times row is ticked. A plot of its own,
    # created here like everything else on this canvas — nothing may be inserted after the
    # window exists (see the note at the top of livecanvas.jl).
    cursor::Makie.Observable{Vector{Makie.Point2f}}
    cursorplot::Any
    # A VECTOR of labels, not a String, even though there is only ever one. `text!` with a
    # positions vector treats a lone String as one text block and then demands exactly one
    # position, so an empty cursor (0 positions, 1 block) fails to resolve — CairoMakie says
    # "Text blocks and positions have different lengths: 1 != 0" and the offscreen save dies.
    # Both observables must move together, which is what the epoch labels already do.
    cursorlabel::Makie.Observable{Vector{String}}
    cursorlabelplot::Any
    message::Makie.Observable{String}
    messageplot::Any
    colormap::Makie.Observable{Any}
    cbarlimits::Makie.Observable{Tuple{Float32,Float32}}
    cbarlabel::Makie.Observable{String}
    colorbar::Any
    # The full extent `show_orbit!` framed. A zoom is bounded against it and a right-click
    # returns to it — the same role `homespan` plays on `SkyCanvas`. A `Ref` because the
    # struct is immutable and this is rewritten on every redraw.
    homespan::Base.RefValue{Float64}
end

function build_orbit_canvas(fig)
    ax = Makie.Axis(fig[1, 1]; xlabel = "x ← E (mas)", ylabel = "y → N (mas)",
                    aspect = Makie.DataAspect())
    ax.xgridvisible = false; ax.ygridvisible = false

    # The rendered surfaces go in FIRST, so the track and the markers draw over them.
    polys  = Makie.Observable(Vector{Makie.Point2f}[])
    pcols  = Makie.Observable(Makie.RGBAf[])
    pplot  = Makie.poly!(ax, polys; color = pcols, strokecolor = pcols, strokewidth = 0.3)

    track  = Makie.Observable(Makie.Point2f[])
    tplot  = Makie.lines!(ax, track; color = (:grey35, 0.9), linewidth = 1.2)
    marks  = Makie.Observable(Makie.Point2f[])
    mplot  = Makie.scatter!(ax, marks; color = (:black, 0.0), strokecolor = :black,
                            strokewidth = 1.2, markersize = 9)
    lpos   = Makie.Observable(Makie.Point2f[])
    ltext  = Makie.Observable(String[])
    lplot  = Makie.text!(ax, lpos; text = ltext, fontsize = 10, color = (:black, 0.75),
                         align = (:left, :bottom), offset = (4, 4))
    prim   = Makie.Observable(Makie.Point2f[])
    pmplot = Makie.scatter!(ax, prim; color = :black, marker = :cross, markersize = 12)

    # FILLED and coloured, against the epoch marks' hollow black rings: the epoch marks say
    # "we observed here" and the cursor says "this is the time on screen". Two things that
    # mean different things must not look alike.
    curs   = Makie.Observable(Makie.Point2f[])
    cplot  = Makie.scatter!(ax, curs; color = (:dodgerblue, 0.95), strokecolor = :black,
                            strokewidth = 0.8, markersize = 11)
    ctext  = Makie.Observable(String[])
    ctplot = Makie.text!(ax, curs; text = ctext, fontsize = 10, color = (:dodgerblue4, 0.95),
                         align = (:left, :top), offset = (5, -5))

    msg  = Makie.Observable("")
    msgp = Makie.text!(ax, [Makie.Point2f(0, 0)]; text = msg, fontsize = 15,
                       color = Makie.RGBAf(0.45, 0.5, 0.55, 1), align = (:center, :center))
    cmap = Makie.Observable{Any}(_padded_cmap("gist_heat"))
    lims = Makie.Observable((0.0f0, 1.0f0))
    # An Observable label: the 3-D components are coloured by temperature and the analytic
    # ones by relative brightness, and a bar that says "T (K)" over a brightness ramp is
    # worse than no bar at all.
    clab = Makie.Observable("T (K)")
    cbar = Makie.Colorbar(fig[1, 2]; colormap = cmap, colorrange = lims, label = clab)

    c = OrbitCanvas(fig, ax, track, tplot, marks, mplot, lpos, ltext, lplot,
                    polys, pcols, pplot, prim, pmplot, curs, cplot, ctext, ctplot,
                    msg, msgp, cmap, lims, clab, cbar,
                    Base.RefValue(0.0))
    idle!(c, "no orbit yet")
    return c
end

"""
    zoom_step!(c::OrbitCanvas, steps; at = nothing) -> c
    reset_zoom!(c::OrbitCanvas) -> c

Wheel zoom for the orbit view, BOUNDED — the same treatment the sky and Mollweide views get.

**Unbounded zoom here could take the machine down.** QMLMakie delivers 120 units per wheel
detent, so Makie's own `ScrollZoom` applies `(1 - speed)^120` for a single notch; on the canvases
that already route through [`zoom_step!`](@ref) that is divided by `WHEEL_DETENT` first, but this
one still had Makie's interaction and so had neither the division nor the span limits. One notch
was enough to drive the limits to something the renderer could not cope with.

Bounds are the shared `ZOOM_MIN_SPAN`/`ZOOM_MAX_SPAN` fractions of the framed extent, so they
mean the same thing for a 1 mas orbit and a 200 mas one.
"""
zoom_step!(c::OrbitCanvas, steps::Real; at = nothing) =
    (zoom_step!(c.axis, c.homespan[], steps; at = at); c)

reset_zoom!(c::OrbitCanvas) = begin
    a = c.homespan[] / 2
    a > 0 || return c
    Makie.xlims!(c.axis, a, -a); Makie.ylims!(c.axis, -a, a)   # East to the left, as framed
    return c
end

_install_zoom!(c::OrbitCanvas) =
    _install_zoom!(c.axis, s -> zoom_step!(c, s), () -> reset_zoom!(c))

_recolor!(::OrbitCanvas) = nothing

"""
    show_orbit!(canvas, o, tepochs; render, ...)

Draw the orbit, the observed epochs on it, and — when `render` is on — the two surfaces at each
of those epochs.

Rendering is OFF by default and that is not only cost: an orbit with six epochs drawn as six
pairs of Roche surfaces is a picture of the geometry, not of the astrometry, and the thing being
looked at here is usually whether the elements put the secondary where the data say it is.
"""
function show_orbit!(c::OrbitCanvas, o::OrbitEntry, tepochs::AbstractVector;
                     nside_exp::Int = 3, T::DataType = Float32, binary = nothing,
                     cursor::Union{Nothing,Real} = nothing)
    bp = orbit_bparams(o; binary = binary)
    busy!(c)
    c.track[]   = _orbit_track_2d(bp; t0 = _track_epoch(bp, tepochs))
    c.primary[] = [Makie.Point2f(0, 0)]

    pts = Makie.Point2f[]; lab = String[]
    for (k, t) in enumerate(tepochs)
        ow, on = orbit_to_rotir_offset(bp, t)
        push!(pts, Makie.Point2f(-ow, on)); push!(lab, string(k))
    end
    c.marks[] = pts; c.labelpos[] = pts; c.labeltext[] = lab

    # THE SYNTHETIC TIME CURSOR — visualization only. `tepochs` still drives the numbered
    # marks, because those mean "we observed here" and a chosen time is not an observation.
    # What a cursor changes is which SINGLE time the surfaces are drawn at, and that is what
    # makes the pair walkable frame by frame: the only way to watch a Roche pair's shape
    # change around an eccentric orbit, since the shape follows D(t).
    if cursor === nothing
        c.cursor[] = Makie.Point2f[]; c.cursorlabel[] = String[]
    else
        cw, cn = orbit_to_rotir_offset(bp, Float64(cursor))
        # POSITION FIRST, then the label: the two are one plot's inputs and a moment with a
        # label and no position to put it at is what the length mismatch above is.
        c.cursor[] = [Makie.Point2f(-cw, cn)]
        c.cursorlabel[] = [Printf.@sprintf("JD %.3f", Float64(cursor))]
    end
    # What the RENDERING sweeps over: every observed epoch, or the one time on the cursor.
    rendert = cursor === nothing ? collect(Float64, tepochs) : [Float64(cursor)]

    if o.render && !isempty(rendert)
        # WHICH stars, from the star model: tessellated surfaces with real Roche shapes and
        # gravity darkening, or the analytic profiles the fit uses. Drawing the 3-D surfaces
        # under an analytic fit would be a picture of a model that is not being fitted.
        polys, cols, lo, hi, lab = o.model === :tessellated ?
            (_render_epochs(o, bp, rendert; nside_exp = nside_exp, T = T)..., "T (K)") :
            (_render_analytic(o, bp, rendert)..., "relative brightness")
        c.polys[] = polys; c.polycolors[] = cols
        c.cbarlimits[] = (Float32(lo), Float32(hi))
        c.cbarlabel[] = lab
        _colorbar_visible!(c, true)
    else
        c.polys[] = Vector{Makie.Point2f}[]; c.polycolors[] = Makie.RGBAf[]
        _colorbar_visible!(c, false)
    end

    # Frame the whole orbit plus whatever the stars add, with a margin.
    r = isempty(c.track[]) ? 1.0 :
        maximum(max(abs(p[1]), abs(p[2])) for p in c.track[])
    r = max(r, o.rpole1 + o.rpole2) * 1.15
    # Re-framed only when the extent actually changed — a different orbit, or stars of a
    # different size. A redraw that leaves the geometry alone keeps the user's zoom, as the
    # sky view does (`_frame!` in src/gui/livecanvas.jl).
    if !(c.homespan[] > 0 && isapprox(c.homespan[], 2r; rtol = 1e-9))
        Makie.xlims!(c.axis, r, -r); Makie.ylims!(c.axis, -r, r)   # East to the left
    end
    c.homespan[] = 2r
    return c
end

"""
    _render_analytic(o, bp, tepochs) -> (polys, colors, lo, hi)

Both components as their ANALYTIC profiles, at each epoch.

The 2-D counterpart of [`_render_epochs`](@ref), and the model the fit actually uses: a
uniform disk is a filled circle, a limb-darkened disk the same circle shaded by its law, a
Gaussian a set of rings following `exp(-4 ln2 (r/FWHM)²)`. Drawn as rings rather than as a
flat blob because the profile is the whole difference between the kinds — a uniform disk and
a Gaussian of the same width are the same picture otherwise, and choosing between them is
what the component selector is for.

Brightness is RELATIVE: the primary peaks at 1 and the secondary at the flux ratio `f`, which
is what `f` means in the fit and the only photometric statement an astrometric orbit makes.
"""
function _render_analytic(o::OrbitEntry, bp, tepochs)
    p = apply_orbit_ties(o)
    f = max(0.0, get(p, :f, 1.0))
    polys = Vector{Makie.Point2f}[]; vals = Float64[]
    for t in tepochs
        ow, on = orbit_to_rotir_offset(bp, t)
        for (kind, pre, cx, cy, scale) in ((o.kind1, :c1, 0.0, 0.0, 1.0),
                                           (o.kind2, :c2, -ow, on, f))
            pl, vl = _profile_rings(kind, p, pre, cx, cy)
            append!(polys, pl); append!(vals, vl .* scale)
        end
    end
    isempty(vals) && return (polys, Makie.RGBAf[], 0.0, 1.0)
    lo, hi = _map_range(vals)
    return (polys, map_colors(vals, _padded_cmap("gist_heat"), (lo, hi)), lo, hi)
end

# One component as filled rings: the polygons and their relative brightness.
#
# `nring` is 10 because the rings are the profile — fewer and a limb-darkened disk reads as a
# uniform one — and the count is fixed rather than adaptive since these are drawn once per
# epoch and never animated.
function _profile_rings(kind::Symbol, p::AbstractDict, pre::Symbol,
                        cx::Real, cy::Real; nring::Int = 10, nang::Int = 72)
    g(n, d) = get(p, Symbol("$(pre)_$(n)"), d)
    # Outer extent, and the brightness at a given SEMI-distance along the major axis.
    if kind === :uniform
        R = g(:diameter, 1.0) / 2;                       I = r -> 1.0
        ratio, pa = 1.0, 0.0
    elseif kind === :limbdark
        R = g(:diameter, 1.0) / 2; u = g(:ld1, 0.3)
        I = r -> (μ = sqrt(max(0.0, 1 - (r / R)^2)); 1 - u * (1 - μ))
        ratio, pa = 1.0, 0.0
    elseif kind === :gaussian
        w = g(:fwhm, 1.0); R = w
        I = r -> exp(-4 * log(2) * (r / w)^2)
        ratio, pa = 1.0, 0.0
    elseif kind === :ellgauss
        w = g(:fwhm, 1.0); R = w
        I = r -> exp(-4 * log(2) * (r / w)^2)
        ratio, pa = clamp(g(:ratio, 1.0), 0.01, 1.0), g(:pa, 0.0)
    else                                   # a point source has no extent to draw
        R = 0.02; I = r -> 1.0
        ratio, pa = 1.0, 0.0
    end
    R > 0 || return (Vector{Makie.Point2f}[], Float64[])

    # Position angle is measured East of North, and the plot's x is East (drawn to the left by
    # the reversed axis), so the major-axis direction is (sin pa, cos pa).
    sp, cp = sind(pa), cosd(pa)
    pt(a, θ) = Makie.Point2f(cx + a * (cos(θ) * sp + ratio * sin(θ) * cp),
                             cy + a * (cos(θ) * cp - ratio * sin(θ) * sp))
    # TWO angular ranges, and the difference is the bug that put a slice out of every disc.
    #
    # `θs` is open — 0 to 2π-Δ — which is what a FILLED polygon wants, since `poly!` closes
    # the last point back to the first itself. An ANNULUS is drawn as one polygon with a slit:
    # the outer arc, then the inner arc reversed. There the arcs have to CLOSE, so both need
    # the 2π endpoint; without it the outer arc stopped at 2π-Δ and the polygon jumped
    # straight inwards, leaving a Δ-wide wedge (5 degrees at `nang = 72`) missing from all
    # nine annuli of every component. Only the central disc, the one case that takes `θs`,
    # came out whole.
    θs = range(0, 2π, length = nang + 1)[1:nang]
    θfull = range(0, 2π, length = nang + 1)

    polys = Vector{Makie.Point2f}[]; vals = Float64[]
    edges = range(0, R, length = nring + 1)
    for k in 1:nring
        rin, rout = edges[k], edges[k + 1]
        ring = if rin == 0
            [pt(rout, θ) for θ in θs]
        else
            # Outer arc, then the inner one reversed: one closed polygon with a hole cut out
            # of it, which is what `poly!` draws without needing a second primitive.
            vcat([pt(rout, θ) for θ in θfull], [pt(rin, θ) for θ in reverse(θfull)])
        end
        push!(polys, ring)
        push!(vals, I((rin + rout) / 2))
    end
    return (polys, vals)
end

"""
    _track_epoch(bp, tepochs) -> t0

Which periastron passage the drawn track belongs to.

It matters only when the apsidal line MOVES: `omega_at` advances ω by `dω·(t − T0)`, so an
ellipse sampled from `T0` carries the orbit's orientation THEN, not at the epochs being marked
on it. Spica's elements are anchored in 2007 and its ω̇ is 0.0071 deg/day, so its 2015 epochs
sit 21 degrees of apsidal rotation away from the `T0` ellipse — off the curve, which is the
visible symptom. The track is drawn over the period containing the MIDDLE of the data instead,
which is the orientation the marks are mostly computed at.

With `dω = 0`, the overwhelmingly common case, this changes nothing at all: every periastron
passage gives the same ellipse.
"""
function _track_epoch(bp, tepochs)
    isempty(tepochs) && return bp.T0
    bp.P > 0 || return bp.T0
    return bp.T0 + round((sum(tepochs) / length(tepochs) - bp.T0) / bp.P) * bp.P
end

# The relative orbit over one period, in PLOT coordinates (x = East to the left, which the
# reversed axis then draws). `orbit_to_rotir_offset` gives (West, North); x = -West.
function _orbit_track_2d(bp; npoints::Int = 360, t0 = bp.T0)
    return [begin
                ow, on = orbit_to_rotir_offset(bp, t0 + f * bp.P)
                Makie.Point2f(-ow, on)
            end for f in range(0, 1, length = npoints)]
end

"""
    _render_epochs(o, bp, tepochs; nside_exp, T) -> (polys, colors, lo, hi)

Both surfaces at every epoch, as one flat polygon list.

This is where the proximity switches act, and only here — `fit_orbit(model = :analytic)` never
sees them:

  * **Roche per epoch**: `create_binary_geometry` recomputes each component's shape from the
    INSTANTANEOUS separation `D(t)`, so on an eccentric orbit the stars breathe. Off, every
    epoch is drawn at the periastron shape, which is what a fixed-shape model assumes.
  * **Mutual irradiation**: `handle_reflection` solves the radiosity between the two surfaces,
    so the facing hemispheres are heated. It is iterative, and of the three switches it has
    by far the highest impact on performance.
  * **Occultation**: `occultation_weights` zeroes the tessels one star hides behind the other,
    which only does anything near conjunction.
"""
function _render_epochs(o::OrbitEntry, bp, tepochs; nside_exp::Int = 3, T::DataType = Float32)
    tess1 = tessellation_healpix(nside_exp; T = T)
    tess2 = tessellation_healpix(nside_exp; T = T)
    base = (surface_type = 3, ldtype = 3, ld1 = 0.2, ld2 = 0.0, beta = 0.25, d = 100.0,
            fillout_factor_primary = -1, fillout_factor_secondary = -1,
            inclination = 180.0 - bp.i, position_angle = bp.Ω - 180.0,
            rotation_period = bp.P, i = bp.i, Ω = bp.Ω, ω = bp.ω, P = bp.P, a = bp.a,
            e = bp.e, T0 = bp.T0, dP = bp.dP, dω = bp.dω)
    p1 = merge(base, (rpole = o.rpole1, tpole = o.tpole1, q = o.q))
    p2 = merge(base, (rpole = o.rpole2, tpole = o.tpole2, q = 1 / o.q))

    polys = Vector{Makie.Point2f}[]; vals = Float64[]
    # Frozen shape: build once at periastron and reuse the geometry's epoch argument only for
    # the orientation. That is what "Roche off" means — a fixed shape, not a missing one.
    tfix = o.roche ? nothing : bp.T0
    for t in tepochs
        te = tfix === nothing ? t : t
        s1, s2 = create_binary_geometry(tess1, p1, tess2, p2, bp,
                                        o.roche ? t : bp.T0)
        # The offsets have to follow the ACTUAL epoch even when the shape is frozen, or the
        # stars are drawn in the wrong place.
        ow, on = orbit_to_rotir_offset(bp, t)
        m1 = Float64.(parametric_temperature_map(p1, s1))
        m2 = Float64.(parametric_temperature_map(p2, s2; secondary = true))
        if o.irradiation
            try
                m1, m2 = handle_reflection(s1, m1, s2, m2)[1:2]
            catch err
                @warn "irradiation skipped at t = $t: $(sprint(showerror, err))"
            end
        end
        w1, w2 = o.occultation ? _occultation_or_ones(s1, s2) :
                 (ones(length(s1.index_quads_visible)), ones(length(s2.index_quads_visible)))
        for (st, mp, w, dx, dy) in ((s1, m1, w1, 0.0, 0.0), (s2, m2, w2, ow, on))
            append!(polys, tessel_polygons(st; offset_west = dx, offset_north = dy))
            append!(vals, Float64.(mp[st.index_quads_visible]) .* w)
        end
    end
    isempty(vals) && return (polys, Makie.RGBAf[], 0.0, 1.0)
    # Occulted tessels come out at zero and would drag the colour scale to black, so the range
    # is taken over what is actually VISIBLE.
    vis = filter(>(0), vals)
    lo, hi = isempty(vis) ? (0.0, 1.0) : _map_range(vis)
    return (polys, map_colors(vals, _padded_cmap("gist_heat"), (lo, hi)), lo, hi)
end

# `occultation_weights` needs both meshes to have their offsets set, which
# `create_binary_geometry` does. A failure here is not worth taking the picture down for.
function _occultation_or_ones(s1, s2)
    try
        w1, w2 = occultation_weights(s1, s2)
        return (Float64.(w1[s1.index_quads_visible]), Float64.(w2[s2.index_quads_visible]))
    catch err
        @warn "occultation skipped: $(sprint(showerror, err))"
        return (ones(length(s1.index_quads_visible)), ones(length(s2.index_quads_visible)))
    end
end

# ── the callbacks ────────────────────────────────────────────────────────────────────────

"""
    shell_orbit_params() -> String

The orbit form, one row per parameter:

    name\tunit\tvalue\tstate\tlo\thi\ttie\tdoc
"""
function shell_orbit_params()
    sh = _sh()
    o = sh.orbit
    rows = String[]
    for n in orbit_param_names(o)
        u, doc = _orbit_info(n)
        v  = get(o.params, n, 0.0)
        lo, hi = get(o.bounds, n, (-Inf, Inf))
        tie = get(o.ties, n, "")
        st = !isempty(tie) ? "tied" : (n in o.free ? "free" : "fixed")
        push!(rows, join((n, u, repr(v), st, repr(lo), repr(hi), tie, doc), "\t"))
    end
    return join(rows, "\n")
end

function shell_set_orbit_param(name, value)
    sh = _sh()
    n = _qmlname(name)
    n === nothing && return NO_QML_NAME
    v = _qmlreal(value, NaN)
    isnan(v) && return "not a number: $(value)"
    sh.orbit.params[n] = v
    refresh_orbit!(sh)
    return ""
end

function shell_set_orbit_state(name, state)
    sh = _sh(); o = sh.orbit
    n = _qmlname(name)
    n === nothing && return NO_QML_NAME
    s = _qmlstr(state)
    if s == "free";      push!(o.free, n); delete!(o.ties, n)
    elseif s == "fixed"; delete!(o.free, n); delete!(o.ties, n)
    elseif s == "tied";  delete!(o.free, n); haskey(o.ties, n) || (o.ties[n] = "")
    else return "unknown state: $s"
    end
    return ""
end

function shell_set_orbit_bound(name, lo, hi)
    sh = _sh()
    n = _qmlname(name)
    n === nothing && return NO_QML_NAME
    l = _qmlreal(lo, NaN); h = _qmlreal(hi, NaN)
    (isnan(l) || isnan(h)) && return "bounds must be numbers"
    l < h || return "lower bound must be below upper"
    sh.orbit.bounds[n] = (l, h)
    return ""
end

function shell_set_orbit_tie(name, expr)
    sh = _sh(); o = sh.orbit
    n = _qmlname(name)
    n === nothing && return NO_QML_NAME
    e = _qmlstr(expr)
    o.ties[n] = e; delete!(o.free, n)
    v = eval_tie(e, o.params)
    v === nothing && return isempty(strip(e)) ? "" : "does not evaluate yet"
    o.params[n] = v
    refresh_orbit!(sh)
    return Printf.@sprintf("= %.6g", v)
end

"""
    shell_orbit_star_models() -> String

`key\tlabel\tdoc` per star model — the choice of WHAT sits at the two orbital positions,
which is a different question from where they are.
"""
shell_orbit_star_models() = join((
    "analytic\tanalytic 2-D profiles\t" *
    "uniform / limb-darkened disks and Gaussians displaced by the orbit — what an " *
    "astrometric orbit fit needs, and what fit_orbit implements",
    "tessellated\ttessellated 3-D components\t" *
    "ROTIR surfaces: Roche shape from the instantaneous separation, gravity darkening, " *
    "irradiation, occultation — rendered here; fit_orbit does not fit them yet"), "\n")

"`analytic` or `tessellated` — which star model the tab is on."
shell_orbit_star_model() = String(_sh().orbit.model)

"""
    shell_set_orbit_star_model(kind) -> String

Switch the star model. The ELEMENTS are untouched: they describe the orbit, not the stars, so
swapping a pair of disks for a pair of Roche lobes must not move the secondary.
"""
function shell_set_orbit_star_model(kind)
    sh = _sh(); o = sh.orbit
    k = Symbol(_qmlstr(kind))
    k in (:analytic, :tessellated) || return "star model must be analytic or tessellated"
    # The 3-D components ARE the Model tab's binary — its two parameter sets, its limb
    # darkening, its β — rather than a second thin description that can disagree with it.
    # Without one there is nothing to render, so this refuses rather than quietly falling back
    # to a hardcoded pair of Roche spheres.
    if k === :tessellated && _orbit_binary(sh) === nothing
        msg = "define a binary on the Model tab first: add a model, tick secondary, and its " *
              "two components become the 3-D stars here"
        # To the CONSOLE as well as the status bar: it is a refusal with an instruction in it,
        # and the status bar is overwritten by the next thing that happens.
        console!(sh, msg)
        return msg
    end
    o.model = k
    # A 3-D component only exists on screen once it is drawn, so choosing it turns rendering
    # on; the switch back does not turn it off again, since by then it is the user's setting.
    k === :tessellated && (o.render = true)
    refresh_orbit!(sh)
    # To the console, not only the status line: this decides what a fit will be handed, and a
    # choice with that reach belongs in the transcript beside the fit it changes.
    console!(sh, "star model: $(k)")
    return "star model: $(k)"
end

"`key\tlabel\tparams` per component kind."
shell_orbit_component_kinds() =
    join(("$(k)\t$(l)\t$(join(String.(ps), ','))" for (k, l, ps) in ORBIT_COMPONENT_KINDS), "\n")

"`kind1\tkind2` — which profile each component currently is."
shell_orbit_components() = (o = _sh().orbit; "$(o.kind1)\t$(o.kind2)")

"""
    shell_set_orbit_component(which, kind) -> String

Change one component's profile. Its parameters are seeded if they are new and LEFT ALONE
otherwise, so switching a disk to a Gaussian and back does not lose the diameter.
"""
function shell_set_orbit_component(which, kind)
    sh = _sh(); o = sh.orbit
    k = Symbol(_qmlstr(kind))
    k in (Symbol(x[1]) for x in ORBIT_COMPONENT_KINDS) || return "unknown kind $(kind)"
    w = _qmlstr(which)
    if w == "1"; o.kind1 = k; _seed_component!(o.params, k, :c1)
    elseif w == "2"; o.kind2 = k; _seed_component!(o.params, k, :c2)
    else return "component must be 1 or 2"
    end
    refresh_orbit!(sh)
    return "component $(w): $(k)"
end

"""
    shell_orbit_options() -> String

`name\tlabel\ton\tdoc` for the rendering switches.
"""
function shell_orbit_options()
    o = _sh().orbit
    # "show stars" applies to BOTH star models — an analytic pair of disks is as drawable as a
    # pair of Roche lobes, and the question "where are the two components at this epoch and do
    # they overlap" is the same question either way. The other three are surface physics and
    # exist only for the tessellated components; offered under the analytic model they would
    # be three ticks that do nothing.
    rows = Tuple{String,String,Bool,String}[
        ("render", "show stars", o.render,
         "draw both components at each observed epoch; off shows the orbit alone, which is " *
         "what the astrometry actually constrains")]
    if o.model === :tessellated
        append!(rows, [
            ("roche", "recompute Roche shape per epoch", o.roche,
             "the shape follows the instantaneous separation D(t); off freezes it at " *
             "periastron"),
            ("irradiation", "mutual irradiation", o.irradiation,
             "solves the radiosity between the surfaces; high impact on performance"),
            ("occultation", "eclipses / occultation", o.occultation,
             "hides the tessels one star is behind; only acts near conjunction")])
    end
    return join(("$(n)\t$(l)\t$(v ? "1" : "0")\t$(d)" for (n, l, v, d) in rows), "\n")
end

function shell_set_orbit_option(name, on)
    sh = _sh(); o = sh.orbit
    v = _qmlstr(on) == "1"
    n = _qmlstr(name)
    n == "render"      ? (o.render = v)      :
    n == "roche"       ? (o.roche = v)       :
    n == "irradiation" ? (o.irradiation = v) :
    n == "occultation" ? (o.occultation = v) :
    return "unknown option $(name)"
    refresh_orbit!(sh)
    return ""
end

"`rpole1\trpole2\ttpole1\ttpole2\tq` — the rendering block, which the analytic fit never sees."
shell_orbit_render_params() = (o = _sh().orbit;
    Printf.@sprintf("%.6g\t%.6g\t%.6g\t%.6g\t%.6g", o.rpole1, o.rpole2, o.tpole1, o.tpole2, o.q))

function shell_set_orbit_render_param(name, value)
    sh = _sh()
    v = _qmlreal(value, NaN)
    isnan(v) && return "not a number: $(value)"
    f = _qmlname(name)
    f === nothing && return NO_QML_NAME
    f in (:rpole1, :rpole2, :tpole1, :tpole2, :q) || return "unknown field $(f)"
    setfield!(sh.orbit, f, v)
    refresh_orbit!(sh)
    return ""
end

"""
    refresh_orbit!(sh) -> String

Redraw the orbit from the current elements and the current dataset's epochs.

With no dataset the track is still drawn — an orbit is a model and can be looked at before any
data are loaded — but with no epoch markers, because there are no epochs.
"""
function refresh_orbit!(sh::ShellState)
    sh.orbitcanvas === nothing && return sh.status
    d = current_dataset(sh.session)
    # JULIAN DATE, not MJD. The elements are quoted in JD — Spica's `T0` is 2454189.4 — and
    # every other consumer of an epoch time converts: `companion_offsets` adds 2400000.5, and
    # `fit_orbit` does it for the whole uv table in `_uv_times`. This did not, so the epoch
    # marks were evaluated 2.4 million days before the elements describe.
    #
    # The visible symptom was the marks not falling on the drawn track, and that is worth
    # spelling out because it is not the obvious consequence. A wrong TIME alone would only
    # put a mark at the wrong PHASE of the right ellipse, still on the curve. What takes it
    # off the curve is apsidal motion: `omega_at` advances ω by `dω·(t − T0)`, so an offset of
    # −2.4e6 days rotated the apsidal line by ~96 degrees for Spica's ω̇ of 0.0071 deg/day
    # (mod 360) while `_orbit_track_2d`, which samples from `T0`, drew the unrotated ellipse.
    # Two different ellipses, which is exactly what it looked like.
    tep = d === nothing ? Float64[] : d.mjd .+ 2_400_000.5
    try
        show_orbit!(sh.orbitcanvas, sh.orbit, tep; binary = _orbit_binary(sh),
                    nside_exp = sh.nside_exp[], T = sh.precision[],
                    cursor = cursor_time(sh))
    catch err
        idle!(sh.orbitcanvas, "could not draw the orbit — check the elements")
        console!(sh, "orbit draw failed: $(sprint(showerror, err))")
    end
    return sh.status
end

# ── the synthetic time axis ──────────────────────────────────────────────────────────────
#
# VISUALIZATION ONLY, and that is the whole design constraint. These times never become epochs:
# the χ² is computed at the observations' times and nowhere else, the numbered marks on the
# orbit keep meaning "we observed here", and nothing here touches the dataset or the model. What
# a cursor decides is which single time the two surfaces are DRAWN at — which is what makes the
# orbit walkable frame by frame, and the only way to watch a Roche pair's shape change around an
# eccentric orbit, since the shape follows D(t).

"""
    cursor_times(sh) -> StepRangeLen or nothing

The times the cursor can sit on, `start:step:stop` in JULIAN DATES.

`nothing` when the Times row is unticked. Defaults, when the fields have not been set, are one
PERIOD FROM PERIASTRON in 48 steps — the same span `_orbit_track_2d` samples, so the cursor
walks exactly the ellipse that is drawn. Taken from the elements each time rather than stored,
so editing `P` or `T0` moves the range with them.
"""
function cursor_times(sh::ShellState)
    t = sh.times
    t.on || return nothing
    v = apply_orbit_ties(sh.orbit)
    P = Float64(get(v, :P, 10.0))
    T0 = Float64(get(v, :T0, 2.45e6))
    a = isfinite(t.start) ? t.start : T0
    b = isfinite(t.stop)  ? t.stop  : T0 + (P > 0 ? P : 1.0)
    st = isfinite(t.step) && t.step > 0 ? t.step : (b - a) / 48
    # A range with a non-positive span or step is not a range: one frame, at the start.
    (st > 0 && b > a) || return range(a, a; length = 1)
    return range(a, b; step = st)
end

"The time on the cursor now, or `nothing` when the Times row is unticked."
function cursor_time(sh::ShellState)
    ts = cursor_times(sh)
    ts === nothing && return nothing
    return ts[clamp(sh.times.index, 1, length(ts))]
end

"""
    shell_times() -> String

`on\tstart\tstop\tstep\tindex\tcount\ttime` — everything the Times row draws, including the
RESOLVED start/stop/step, so the fields show the defaults that are in force rather than blanks
the reader has to know the meaning of.
"""
function shell_times()
    sh = _sh()
    ts = cursor_times(sh)
    # Resolved even when the row is unticked: ticking it should reveal the numbers it will use,
    # not empty boxes.
    was = sh.times.on
    sh.times.on = true
    tr = cursor_times(sh)
    sh.times.on = was
    n = length(tr)
    i = clamp(sh.times.index, 1, n)
    return join((was ? "1" : "0", _fmt_param(first(tr)), _fmt_param(last(tr)),
                 _fmt_param(n > 1 ? step(tr) : 0.0), i, n,
                 _fmt_param(ts === nothing ? first(tr) : ts[i])), "\t")
end

"""
    shell_set_times(on, start, stop, step) -> String

Tick the Times row and set its range. Empty or unparseable fields keep the resolved default.
"""
function shell_set_times(on, start, stop, step)
    sh = _sh()
    sh.times.on = _qmlstr(on) == "1"
    a = _qmlreal(start, NaN); b = _qmlreal(stop, NaN); st = _qmlreal(step, NaN)
    sh.times.start = a
    sh.times.stop  = b
    sh.times.step  = st
    ts = cursor_times(sh)
    ts === nothing && (refresh_orbit!(sh); return "times off — the marks are the observations")
    sh.times.index = clamp(sh.times.index, 1, length(ts))
    refresh_orbit!(sh)
    return Printf.@sprintf("times: %d frames, JD %.3f to %.3f", length(ts), first(ts), last(ts))
end

"""
    shell_set_time_index(i) -> String
    shell_step_time(delta) -> String

Move the cursor: to a frame, or by a number of frames. Both clamp rather than wrap — an orbit
is periodic but a TIME RANGE is not, and silently jumping from the last frame back to the first
hides that the end has been reached.
"""
function shell_set_time_index(i)
    sh = _sh()
    ts = cursor_times(sh)
    ts === nothing && return "times off"
    sh.times.index = clamp(round(Int, _qmlreal(i, 1.0)), 1, length(ts))
    refresh_orbit!(sh)
    return Printf.@sprintf("JD %.3f  (frame %d of %d)", ts[sh.times.index],
                           sh.times.index, length(ts))
end

function shell_step_time(delta)
    sh = _sh()
    ts = cursor_times(sh)
    ts === nothing && return "times off"
    return shell_set_time_index(clamp(sh.times.index + round(Int, _qmlreal(delta, 0.0)),
                                      1, length(ts)))
end

"""
    shell_fit_orbit(method, maxeval) -> String

Fit the free parameters against the loaded epochs, on a worker.

The star model chosen on the tab is passed straight through as `fit_orbit`'s `model`. Only
`:analytic` is implemented; `:tessellated` throws there, with `fit_orbit`'s own explanation of
what a tessellated fit would need. Passing it rather than intercepting it is deliberate — one
refusal, from the function that owns the limitation, instead of a second copy of it here that
can go stale.

On the analytic path the proximity switches change the RENDERING and not this χ²: analytic
profiles have no surface to shade.
"""
function shell_fit_orbit(method, maxeval)
    sh = _sh(); o = sh.orbit
    d = current_dataset(sh.session)
    d === nothing && return "no dataset"
    isempty(o.free) && return "nothing is free — mark at least one parameter free"
    meth = Symbol(_qmlstr(method))
    # No `:ultranest`: the GUI is Python-free by construction, so the method that would need
    # PythonCall is not offered here at all. `fit_orbit(...; method = :ultranest)` still works
    # from a script that loads it.
    meth in (:neldermead, :nautilus) ||
        return "orbit fitting takes :neldermead or :nautilus"
    meth === :nautilus && !nautilus_available() &&
        return "add `using Nautilus` to this session"

    c1 = _orbit_component(o.kind1, apply_orbit_ties(o), :c1)
    c2 = _orbit_component(o.kind2, apply_orbit_ties(o), :c2)
    el = orbit_elements_nt(o)
    # FILTERED to the parameters this component pair actually has. The session keeps bounds for
    # every kind so that switching a Gaussian to a disk and back does not lose them, but
    # `orbit_fit_spec` rejects a bound for a parameter that is not in the vector — and its
    # error names the parameter without saying that the component kind is why.
    # `analytic = true` unconditionally: a fit IS analytic — `fit_orbit(model = :tessellated)`
    # throws by design — so the vector it assembles has the component parameters in it whatever
    # the tab is currently DRAWING. Without this, selecting the 3-D star model would silently
    # drop a free `c1_diameter` from the fit instead of letting `fit_orbit` refuse with its own
    # explanation.
    known = Set(orbit_param_names(o; analytic = true))
    fr = sort([n for n in o.free if n in known])
    bd = Dict(k => v for (k, v) in o.bounds if k in known)
    ti = Dict(k => v for (k, v) in o.ties if k in known && !isempty(strip(v)))
    fl = get(o.params, :f, 1.0)
    data = d.data
    it = max(10, Int(maxeval))
    console!(sh, "fit_orbit(data, $(o.kind1), $(o.kind2); free = $(fr), " *
                 "method = :$(meth), model = :$(o.model))"; kind = :cmd)
    log!(sh.session,
         "res = fit_orbit(data, " * _component_source(o.kind1, o.params, :c1) * ", " *
         _component_source(o.kind2, o.params, :c2) * ";\n" *
         "                elements = " * _namedtuple_literal(el) * ",\n" *
         "                flux_ratio = $(fl), free = $(_literal(fr)),\n" *
         "                bounds = $(_literal(bd)), ties = $(_literal(ti)),\n" *
         "                method = :$(meth), model = :$(o.model), maxeval = $(it))";
         note = "orbit fit", binding = "res")
    mdl = o.model
    return start_job!(sh, :orbit, function (stop)
        fit_orbit(data, c1, c2; elements = el, flux_ratio = fl, free = fr, bounds = bd,
                  ties = ti, method = meth, model = mdl, maxeval = it, verbose = true)
    end)
end

# The component constructor as source, for the exported script.
function _component_source(kind::Symbol, p::AbstractDict, pre::Symbol)
    g(n, d) = get(p, Symbol("$(pre)_$(n)"), d)
    kind === :uniform  && return "UniformDisk(diameter = $(g(:diameter, 1.0)))"
    kind === :limbdark && return "LimbDarkenedDisk(diameter = $(g(:diameter, 1.0)), " *
                                 "law = :linear, ld1 = $(g(:ld1, 0.3)))"
    kind === :gaussian && return "GaussianDisk(fwhm = $(g(:fwhm, 1.0)))"
    kind === :ellgauss && return "EllipticalGaussian(fwhm = $(g(:fwhm, 1.0)), " *
                                 "ratio = $(g(:ratio, 1.0)), pa = $(g(:pa, 0.0)))"
    return "PointSource()"
end

"The last orbit fit, as `name\tvalue\terror` rows."
shell_orbit_result() = _sh().lastorbit

"""
    orbit_dir() -> String

The folder orbits are loaded from and saved to, created on first use.

Per USER rather than per project, beside the appearance settings and for the same reason: an
orbit is a description of a SYSTEM, not of a dataset, and the β Lyr elements are the same ones
whichever night's data is open. It also has to be writable, which a package directory under a
system-wide depot is not.

Seeded, the first time it is created, with the orbits that ship in `demos/orbits` — so the
picker opens on something rather than on an empty folder. Only on creation: a file deleted
from here stays deleted.
"""
function orbit_dir()
    base = Sys.iswindows() ? get(ENV, "APPDATA", homedir()) :
           get(ENV, "XDG_CONFIG_HOME", joinpath(homedir(), ".config"))
    dir = joinpath(base, "rotir", "orbits")
    isdir(dir) && return dir
    try
        mkpath(dir)
        shipped = ROTIR.resource("demos", "orbits")
        if shipped !== nothing
            for f in readdir(shipped)
                endswith(f, ".toml") || continue
                cp(joinpath(shipped, f), joinpath(dir, f); force = false)
            end
        end
    catch err
        # A read-only home is not a reason to fail to open a picker.
        @warn "could not create the orbit folder" dir exception = err
        return pwd()
    end
    return dir
end

"""
    image_dir() -> String

Where saved views go, created on first use.

Beside the orbits, and per user for the same reason: a PNG of a fit is a thing you keep, and the
working directory a GUI happens to have been launched from is not somewhere anyone looks for it
later. It was the working directory until now, so "Save view" scattered `rotir_star3d-2.png`
through whatever tree the launcher was run in.
"""
function image_dir()
    base = Sys.iswindows() ? get(ENV, "APPDATA", homedir()) :
           get(ENV, "XDG_CONFIG_HOME", joinpath(homedir(), ".config"))
    dir = joinpath(base, "rotir", "images")
    isdir(dir) && return dir
    try
        mkpath(dir)
    catch err
        @warn "could not create the image folder" dir exception = err
        return pwd()
    end
    return dir
end

"Put a bare filename in [`orbit_dir`](@ref); leave anything with a directory alone."
_in_orbit_dir(p::AbstractString) = isempty(dirname(p)) ? joinpath(orbit_dir(), p) : p

function shell_save_orbit(path)
    sh = _sh()
    p = String(path); startswith(p, "file://") && (p = p[8:end])
    isempty(splitext(p)[2]) && (p *= ".toml")
    # A bare name from the Save button lands in the orbit folder rather than in whatever
    # directory Julia happens to have been started from.
    p = _in_orbit_dir(p)
    try
        p = unique_path(p)
        save_orbit(sh.orbit, p)
        sh.status = "wrote $(p)"
    catch err
        sh.status = "could not write $(p): $(sprint(showerror, err))"
    end
    console!(sh, sh.status)
    return sh.status
end

function shell_load_orbit(path)
    sh = _sh()
    p = String(path); startswith(p, "file://") && (p = p[8:end])
    try
        sh.orbit = load_orbit(p)
        sh.status = "loaded $(basename(p))"
        console!(sh, sh.status)
        refresh_orbit!(sh)
    catch err
        sh.status = "could not read $(basename(p)): $(sprint(showerror, err))"
        console!(sh, sh.status)
    end
    return sh.status
end
