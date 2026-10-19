# Everything QML calls into.
#
# One rule shapes this file: **nothing crosses the QML boundary except strings and numbers.**
# No Observables, no JuliaPropertyMap, no Julia objects. A table crosses as `\n`-separated rows
# of `\t`-separated fields; a boolean crosses as "1"/"0". That is not a limitation of QML.jl —
# it is what makes the shell testable as text, with no display and no Qt, which is how every
# callback below is exercised.
#
# QML holds no state of its own. Every edit goes to Julia and the table is re-read, so what is
# on screen is always what a script would see. That is what lets the command log claim it
# reproduces the session.
#
# The other rule: **GL work only on the GUI thread.** Long-running calls go to a worker via
# `GuiJob`; the worker computes and nothing else. Every canvas update happens in `job_finish!`,
# back on the GUI thread.

"""
    TimeCursor

The Orbit tab's synthetic time axis: `start:step:stop`, and which of those times is current.

**Visualization only.** These times are not epochs and never become any: the χ² is computed at
the observations' times and nowhere else, and the epoch marks on the orbit keep meaning "we
observed here". What this drives is the picture — where the cursor sits on the track and, when
rendering is on, which single time the two surfaces are drawn at. That is what makes the orbit
walkable frame by frame, which is the only way to watch a Roche pair's shape change around an
eccentric orbit, since the shape follows `D(t)`.

Times are JULIAN DATES, like `T0` and like everything `orbit_to_rotir_offset` reads — the one
place that passed MJD is exactly where the epoch marks went wrong (see plans/gui_todo.md #4).

Lives on `ShellState` rather than on `OrbitEntry`: an `OrbitEntry` is saved and loaded as TOML
and describes a SYSTEM, while a cursor describes what is on screen right now.
"""
mutable struct TimeCursor
    on::Bool
    start::Float64
    stop::Float64
    step::Float64
    index::Int         # 1-based, into `start:step:stop`
end

"A cursor covering one period from periastron, which is what the track already samples."
TimeCursor() = TimeCursor(false, NaN, NaN, NaN, 1)

"""
    ShellState

One mutable bundle holding the session, the canvases and the console, reached through the
global `SHELL[]` because Qt callbacks take no arguments.

**Constructed by KEYWORD, and every field but `session` has a default.** It was positional, and
with forty fields that made adding one a four-place edit: `gui()`, the headless test harness,
and — the one that is easy to miss, because it lives in another file and only runs at
precompile time — the `@compile_workload` in `ext/ROTIRGUIExt.jl`. Adding `times` broke exactly
that third site, and the failure is a `MethodError` printing all thirty-nine argument types on
one line, which says nothing about which field is missing.

Three comments in this struct used to end "appended rather than placed with the other canvases
because every positional construction would otherwise have to be renumbered". That reason is
gone: a new field can go where it belongs, with a default, and no call site changes at all.

The defaults are the values every caller was spelling out identically anyway — and they had
already drifted, with the workload passing `nothing` for the orbit canvas and the harness `""`
for a status the window computes.
"""
Base.@kwdef mutable struct ShellState
    session::Session
    star::Any           = nothing   # StarCanvas
    moll::Any           = nothing   # MollCanvas
    # The Imaging tab's own pair. Separate figures, not a second view of the two above: a Makie
    # Figure belongs to one MakieArea at a time, and a reconstruction must not overwrite the
    # plot being read on another tab.
    imsky::Any          = nothing   # SkyCanvas
    immoll::Any         = nothing   # MollCanvas
    msky::Any           = nothing   # SkyCanvas — the Model tab's own, and its DEFAULT view
    # The Explore tab's observable plot: an OITOOLS `LiveCanvas`, reached through its GUI
    # extension rather than reimplemented. See `build_obs_canvas`.
    obs::Any            = nothing
    obskind::Base.RefValue{Symbol} = Ref(:v2)
    obscolor::Base.RefValue{Symbol} = Ref(:baseline)
    # What the data are compared with: :data, :model, :image, :residmodel or :residimage.
    obscompare::Base.RefValue{Symbol} = Ref(:model)
    # One panel per baseline (triplet, station) against wavelength, as OITOOLS draws it.
    obspanels::Base.RefValue{Bool} = Ref(false)
    # The description of each point on the observable plot, index-aligned with what is drawn
    # (`update_canvas!` returns it), and the last one a click identified.
    obsinfo::Vector{String} = String[]
    picktext::String = ""
    # What the surface views COLOUR BY: the temperature map itself, or the emergent intensity —
    # limb darkening multiplied in, and optionally a Planck conversion at the observing
    # wavelength. Two different physical quantities, and the one you want depends on whether
    # you are reading the model or reading what the interferometer sees.
    intensity::Base.RefValue{Bool} = Ref(true)
    intensity_model::Base.RefValue{Symbol} = Ref(:linear)
    band::Base.RefValue{Float64} = Ref(0.0)   # metres; 0 means "take it from the data"
    # Which decorations the sky views draw. One set for every panel, because they are the same
    # picture seen on different tabs and three different sets of annotations for one star would
    # be three things to keep in your head rather than one.
    decor::Dict{Symbol,Bool} = Dict(:limb => true, :compass => true,
                                    :graticules => false, :spin => false,
                                    :plotmesh => false)
    # Graticule settings. Spacing in DEGREES, which is how anyone reading a map thinks about
    # it; `graticule_segments` counts lines, so the two are converted where they meet.
    gratlat::Base.RefValue{Float64} = Ref(30.0)
    gratlon::Base.RefValue{Float64} = Ref(30.0)
    gratcolor::Base.RefValue{String} = Ref("black")
    # The colour range of every surface view: `nothing` autoscales each panel to its own map,
    # a `(lo, hi)` pair pins them all to the same scale so two maps can be compared by eye.
    # In the units of whatever is being coloured (K or arbitrary intensity), so it reverts to
    # autoscale when that quantity changes.
    colorrange::Base.RefValue{Any} = Ref{Any}(nothing)
    # Which observables enter the χ², as `(V², T3amp, T3φ)` weights of 1 or 0 — one set for
    # the Model tab's fits and one for the Imaging tab, because a common sequence is to fit the
    # shape on V² and T3amp and then image with the closure phases as well.
    fitweights::Base.RefValue{NTuple{3,Float64}} = Ref((1.0, 1.0, 1.0))
    imgweights::Base.RefValue{NTuple{3,Float64}} = Ref((1.0, 1.0, 1.0))
    # The mesh the model is evaluated on, and the float type it is built in. Both belong to the
    # MODEL rather than to a plot: they decide what is being fitted, not how it is drawn.
    tessel::Base.RefValue{Symbol} = Ref(:healpix)   # :longlat is not wired yet
    nside_exp::Base.RefValue{Int} = Ref(3)   # log2(nside); 3 for a result, 2 to try something
    precision::Base.RefValue{DataType} = Ref{DataType}(Float32)
    # WHO PICKS THE FORWARD KERNEL. `true` — the default — means the CODE PATH picks: `:t3`
    # for the interactive χ² and the derivative-free engines, which reach its branch through
    # `fused_cvis`, and `:turbo` for the gradient path and the samplers, which go through
    # `interferometric_chi2` / `shape_chi2_fg!` and have no `:nufft` branch to reach. `false`
    # means the user forced one from the panel and every path honours it, which is what makes
    # "run the reference kernel and see if the answer moves" possible.
    kernel_auto::Base.RefValue{Bool} = Ref(true)
    # The Orbit perspective. One orbit at a time, like the model — a session fitting two
    # unrelated binaries at once is not a thing anyone does, and a selector for it would be
    # another way to have the wrong one selected.
    orbit::Any                  = default_orbit()   # OrbitEntry
    orbitcanvas::Any            = nothing   # OrbitCanvas
    lastorbit::String           = ""        # the fit result table
    console::Vector{String} = String[]
    status::String = ""
    job::Any            = nothing   # GuiJob or nothing
    jobkind::Symbol = :none
    # WHICH ENGINE the running fit is using. `jobkind` says what kind of job it is
    # (:fit / :reconstruct / :orbit); this says which optimiser, which is what decides whether
    # Stop can do anything at all — see `STOP_AWARE`.
    jobmethod::Symbol = :none
    lastfit::String     = ""        # rendered result table of the last fit
    imaging_defaults::Dict{Symbol,Any} = Dict{Symbol,Any}()
    # Cache for `epoch_chi2`, keyed on what it actually depends on. Recomputing it is seconds
    # of work — the geometry at every epoch plus `setup_oi!`, i.e. the polygon FT against every
    # uv point — and before this it ran on every tab switch.
    chi2cache::Base.RefValue{Any} = Ref{Any}(nothing)
    chi2key::Base.RefValue{Any} = Ref{Any}(nothing)
    # The posterior panel: the canvas, and which two parameters it is showing.
    post::Any                   = nothing   # PostCanvas
    postx::Base.RefValue{Int} = Ref(1)      # marginal / pair x
    posty::Base.RefValue{Int} = Ref(2)      # pair y
    # Whether the observable plot's y axis is logarithmic.
    obslog::Base.RefValue{Bool} = Ref(false)
    # The Imaging tab's own 3-D scene. A reconstruction is a map like any other and is looked
    # at the same three ways; it needs a scene of its own because a Figure belongs to one
    # MakieArea at a time — pointing two tabs at one figure leaves one blank and draws the
    # other detached across the panel beside it.
    imstar::Any                 = nothing   # StarCanvas
    # The map the Imaging views are currently showing, as `(star, x, title)` — or `nothing`
    # before the first reconstruction. Kept because those views are drawn ONCE, when a run
    # finishes, so a later change to how a map is drawn (the intensity tick, a decoration, the
    # graticule) had nothing to redraw from and silently did nothing on this tab.
    lastmap::Base.RefValue{Any} = Ref{Any}(nothing)
    # The Orbit tab's synthetic time cursor.
    times::TimeCursor = TimeCursor()
end

const SHELL = Ref{Any}(nothing)

"""
Whether Qt's event loop is turning, so a redraw asked for from QML can actually complete.

**A Makie redraw before the first frame deadlocks the window, silently.** Qt's render thread
sits in `MakieViewport::setup_buffer` waiting on `ForeignThreadManager::begin_julia` for the
Julia lock; a callback that redraws during `Component.onCompleted` holds that lock through a
full Makie pass, the render thread never gets in, and the window stays BLACK for ever. Nothing
throws, so stdout, the console and the QML error channel are all empty — the only evidence is
a backtrace showing the render thread blocked on a mutex.

Set from [`shell_job_poll`](@ref), which is driven by a QML `Timer` and therefore cannot run
until the loop is live. Cleared by `gui` before `loadqml`, so a second window in one session
starts from false rather than inheriting the first one's.
"""
const GUI_LIVE = Ref(false)

_sh() = SHELL[]::ShellState

# ── console ─────────────────────────────────────────────────────────────────────────────

"Append a line to the console pane. `kind = :cmd` marks a replayable call."
function console!(sh::ShellState, line::AbstractString; kind::Symbol = :info)
    push!(sh.console, kind === :cmd ? "» " * String(line) : String(line))
    length(sh.console) > 4000 && deleteat!(sh.console, 1:500)
    return line
end

"""
    shell_refresh() -> String

Redraw every canvas from the current session, and return the status line.

Called from QML's `refreshAll()`, which runs once the window is up and again on every tab
change. Without it a session populated BEFORE `gui()` — the launcher's command-line
arguments, say — showed empty plots until something was edited, because the canvases are only
filled by the actions that change them.
"""
function shell_refresh()
    sh = _sh()
    refresh_data_tab!(sh)
    refresh_model_tab!(sh)
    refresh_posterior!(sh)
    refresh_orbit!(sh)
    return sh.status
end

shell_console() = join(_sh().console, "\n")
shell_ready() = "1"
shell_status() = _sh().status

"""
    shell_ui_scale(scale, dpi) -> ""

QML reports what it measured, once, from `Component.onCompleted`. See src/gui/scaling.jl for
why detection lives on that side.
"""
function shell_ui_scale(scale, dpi)
    set_ui_scale!(Float64(scale), Float64(dpi))
    # OITOOLS' live canvas sizes its type from ITS OWN copy of the scale, so the value has to
    # be pushed there too or the observable plot comes out at a different size from everything
    # beside it.
    O = _oitools_gui()
    O === nothing || O.set_ui_scale!(Float64(scale), Float64(dpi))
    return ""
end

# ── datasets ────────────────────────────────────────────────────────────────────────────

"""
    shell_open_many(paths, add = "0") -> String

Load several OIFITS at once. `paths` is newline-separated, as the picker builds it.

One call, not a loop over `shell_open`: `readoifits_multiepochs` reads the whole set together,
and the epoch times are only meaningful relative to the earliest of them — opening six files
one at a time would recompute that origin six times and log six lines for one action.
"""
function shell_open_many(paths, add = "0", dosplit = "1")
    sh = _sh()
    ps = String[strip(replace(p, r"^file://" => "")) for p in split(String(paths), '\n')
                if !isempty(strip(p))]
    isempty(ps) && return sh.status
    length(ps) == 1 && return shell_open(ps[1], add, dosplit)
    append = _qmlstr(add) == "1" && current_dataset(sh.session) !== nothing
    try
        e = load_dataset!(sh.session, ps; add = append, split = _qmlstr(dosplit) == "1")
        _log_dataset!(sh, e)
        console!(sh, "$(append ? "add" : "load") $(length(ps)) files — $(e.name) now has " *
                     "$(length(e.data)) epoch(s), $(sum(e.nv2)) V², $(sum(e.nt3)) T3")
        sh.status = "$(e.name): $(length(e.data)) epoch(s)"
        sh.band[] <= 0 && (sh.band[] = _data_band(sh))
        refresh_data_tab!(sh)
    catch err
        sh.status = "could not read $(length(ps)) files: $(sprint(showerror, err))"
        console!(sh, sh.status)
    end
    return sh.status
end

"""
    shell_remove_epoch(i) -> String

Drop epoch `i` from the current dataset.

Dropping the LAST epoch removes the dataset itself: a dataset with no data is not a state the
rest of the shell has any use for, and leaving an empty row in the selector is worse than
leaving nothing.

The epoch times are recomputed afterwards, because they are measured from the earliest
remaining epoch — removing the first one shifts every other.
"""
function shell_remove_epoch(i)
    sh = _sh()
    d = current_dataset(sh.session)
    d === nothing && return sh.status
    n = Int(i)
    n in eachindex(d.data) || return sh.status
    for f in (:data, :mjd, :tepochs, :nv2, :nt3); deleteat!(getfield(d, f), n); end
    if isempty(d.data)
        return shell_close_dataset()
    end
    d.tepochs = d.mjd .- minimum(d.mjd)
    sh.session.current_epoch = clamp(sh.session.current_epoch, 1, length(d.data))
    sh.chi2key[] = nothing                       # the answer depended on the removed epoch
    _log_dataset!(sh, d)
    console!(sh, "removed epoch $n; $(length(d.data)) left")
    sh.status = "$(d.name): $(length(d.data)) epoch(s)"
    refresh_data_tab!(sh)
    return sh.status
end

"""
    shell_close_dataset() -> String

Remove the current dataset entirely, and select the one before it.
"""
function shell_close_dataset()
    sh = _sh()
    d = current_dataset(sh.session)
    d === nothing && return sh.status
    deleteat!(sh.session.datasets, sh.session.current_dataset)
    sh.session.current_dataset = min(sh.session.current_dataset, length(sh.session.datasets))
    sh.session.current_epoch = 1
    sh.chi2key[] = nothing
    # The log entry for the closed dataset goes with it: an exported script that loaded a file
    # and then never used it is a script that does not reproduce the session.
    filter!(e -> e.binding != "data", sh.session.log)
    nd = current_dataset(sh.session)
    nd === nothing || _log_dataset!(sh, nd)
    console!(sh, "closed $(d.name)")
    sh.status = nd === nothing ? "no dataset loaded" :
                "$(nd.name): $(length(nd.data)) epoch(s)"
    refresh_data_tab!(sh)
    return sh.status
end

"""
    shell_open(path, add = "0", dosplit = "1") -> String

Load an OIFITS and show its first epoch.

`add = "1"` appends it to the CURRENT dataset rather than starting a new one — which is what
opening a second night means for a target like λ And, whose six nights live in six files.
`dosplit = "0"` takes the file as ONE epoch instead of cutting it on gaps in its V² timestamps.

Returns a status line.
"""
function shell_open(path, add = "0", dosplit = "1")
    sh = _sh()
    p = String(path)
    startswith(p, "file://") && (p = p[8:end])
    append = _qmlstr(add) == "1" && current_dataset(sh.session) !== nothing
    try
        e = load_dataset!(sh.session, p; add = append, split = _qmlstr(dosplit) == "1")
        _log_dataset!(sh, e)
        console!(sh, "$(append ? "add" : "load"): $(basename(p)) — $(e.name) now has " *
                     "$(length(e.data)) epoch(s), $(sum(e.nv2)) V², $(sum(e.nt3)) T3")
        sh.status = "$(e.name): $(length(e.data)) epoch(s)"
        sh.band[] <= 0 && (sh.band[] = _data_band(sh))
        refresh_data_tab!(sh)
    catch err
        sh.status = "could not read $(basename(p)): $(sprint(showerror, err))"
        console!(sh, sh.status)
    end
    return sh.status
end

# The log records the whole dataset as ONE call each time it changes, rather than one line per
# click. That is what the equivalent script looks like — nobody writes six `readoifits` calls
# and then sorts them — and it is what makes the exported script reproduce the epoch ORDER,
# which the rotation phase depends on.
function _log_dataset!(sh::ShellState, e::DatasetEntry)
    isempty(sh.session.log) || (sh.session.log[end].binding == "data" &&
                                 pop!(sh.session.log))
    code = length(e.paths) == 1 ?
        "data = split_epochs(readoifits($(_literal(e.paths[1])); warn = false, " *
        "verbose = false)[1,1])[1]" :
        "files = " * _literal(e.paths) * "\n" *
        "data  = collect(readoifits_multiepochs(files; warn = false, verbose = false)[1, :])"
    code *= "\ntepochs = [d.mean_mjd for d in data]; tepochs .-= minimum(tepochs)"
    log!(sh.session, code; note = "dataset $(e.name)", binding = "data")
    return e
end

"""
    shell_datasets() -> String

`name\tnepochs\tnV2\tnT3` per dataset, one per line. Empty when nothing is loaded.
"""
shell_datasets() = join(("$(d.name)\t$(length(d.data))\t$(sum(d.nv2))\t$(sum(d.nt3))"
                         for d in _sh().session.datasets), "\n")

"Index of the current dataset, 1-based; 0 when none is loaded."
shell_current_dataset() = string(_sh().session.current_dataset)

"""
    shell_current_epoch() -> String

Which epoch the session is on, 1-based.

The table has to READ this rather than remember its own selection. `load_dataset!` with
`add = true` moves the session to the epoch it just appended, and a table that only restored
its previous row then highlighted epoch 1 while every plot beside it drew epoch 2 — the two
disagreeing silently, which is worse than either being wrong.
"""
shell_current_epoch() = string(_sh().session.current_epoch)

"""
    shell_epochs() -> String

`index\tMJD\tdays\tv2\tt3amp\tt3phi\thave_chi2` per epoch of the current dataset.

Both times are shown: `MJD` identifies the night, `days` is the time since the first epoch and
is what the rotation phase is actually computed from.

The three observables are SEPARATE columns, because that is the question the panel is for: a
model can fit V² well and closure phase badly, and a single combined χ² hides exactly that. Each
cell carries the reduced χ² for that observable when a model is set, and the point COUNT
otherwise — `have_chi2` says which, so the header can be labelled honestly rather than showing
"χ²ᵣ" above a column of counts.

`chi2r` is empty rather than 0 when there is no model to compare against: an empty cell reads
as "not computed", a zero reads as a perfect fit.
"""
function shell_epochs()
    sh = _sh()
    d = current_dataset(sh.session)
    d === nothing && return ""
    chi2 = epoch_chi2(sh)
    rows = String[]
    for i in eachindex(d.data)
        if chi2 === nothing
            v2, ta, tp = string(d.nv2[i]), string(d.data[i].nt3amp), string(d.nt3[i])
        else
            b = chi2[i]
            f(x) = isfinite(x) ? Printf.@sprintf("%.3g", x) : "—"
            v2, ta, tp = f(b.v2r), f(b.t3ampr), f(b.t3phir)
        end
        push!(rows, Printf.@sprintf("%d\t%.4f\t%.3f\t%s\t%s\t%s\t%s", i, d.mjd[i],
                                    d.tepochs[i], v2, ta, tp, chi2 === nothing ? "0" : "1"))
    end
    return join(rows, "\n")
end

"""
    shell_select_dataset(i) -> String

Make dataset `i` current. One-based, and out of range is a no-op rather than an error: the
index comes from a list view whose model can be a frame behind.
"""
function shell_select_dataset(i)
    sh = _sh()
    n = Int(i)
    n in eachindex(sh.session.datasets) || return sh.status
    sh.session.current_dataset = n
    sh.session.current_epoch = 1
    refresh_data_tab!(sh)
    return sh.status
end

function shell_select_epoch(i)
    sh = _sh()
    d = current_dataset(sh.session)
    d === nothing && return sh.status
    n = clamp(Int(i), 1, length(d.data))
    sh.session.current_epoch = n
    refresh_data_tab!(sh)
    return sh.status
end

# ── models ──────────────────────────────────────────────────────────────────────────────

"""
    shell_surface_types() -> String

`code\tname\tlabel\tdoc` for every implemented surface type, from the schema — so the selector
cannot offer a type `compute_radii` does not branch on — plus BINARY.

Binary is not a surface type and has no schema entry: it is two components, each with a surface
type of its own. It belongs in this list all the same, because this list answers "what am I
building" and a binary is one of the answers. `BINARY_CODE` is negative so it cannot collide
with a real `surface_type`, and `shell_add_model` recognises it.
"""
const BINARY_CODE = -1

shell_surface_types() =
    join(vcat(["$(c)\t$(SURFACE_TYPES[c].name)\t$(SURFACE_TYPES[c].label)\t" *
               replace(SURFACE_TYPES[c].doc, "\n" => " ")
               for c in SURFACE_TYPE_ORDER],
              ["$(BINARY_CODE)\tbinary\tRoche Binary\t" *
               "two Roche components sharing one orbit, which is taken from the Orbit tab. " *
               "Each component keeps its own surface type, radius, temperature and limb " *
               "darkening; change either with the selectors that appear. The χ² is that of " *
               "the PAIR against the data"]), "\n")

"""
    shell_add_model(surface_type) -> String

A new model at the schema defaults. Status line back.
"""
function shell_add_model(surface_type)
    sh = _sh()
    # ONE model at a time. A session with several was offering a selector for something nobody
    # was using, and the second model then silently decided what the χ² column and the
    # reconstruction were about. Replacing is what "+ model" means here.
    empty!(sh.session.models)
    sh.session.current_model = 0
    sh.chi2key[] = nothing
    st = Int(surface_type)
    # BINARY: a Roche primary and a Roche companion. Roche because that surface type exists
    # precisely for a star in a binary — its shape follows the instantaneous separation — and
    # because a binary of two spheres is the case where the orbit is all that matters, which
    # the Orbit tab already fits analytically and faster.
    if st == BINARY_CODE
        m = add_model!(sh.session, 3)
        shell_set_binary("1", 3)
        # BY THE ORBIT, unlike the `secondary` tick. Ticking `secondary` on an arbitrary model
        # defaults to a fixed offset because one displacement is what a snapshot constrains
        # (see `DEFAULT_COMPANION_OFFSET`); but a Roche shape is computed FROM the
        # instantaneous separation, so this entry already has an orbit in it, and placing the
        # pair by a hand-typed offset while shaping it from an orbit says two things at once.
        m.companion.place = :orbit
        # AND THE ORBIT TAB TAKES ITS ELEMENTS FROM THE MODEL, not the other way round.
        #
        # This direction is the whole difference between a usable default and a broken one.
        # The Roche schema's defaults are SELF-CONSISTENT — rpole 0.5 mas against a = 3 mas at
        # q = 0.5 puts each component at 38% of its Roche lobe — while the Orbit tab's default
        # is a generic two-disk orbit at a = 1 mas that knows nothing about them. Seeding the
        # other way let the first refresh overwrite a = 3 with a = 1, which is rpole/a = 0.5
        # against a lobe of 0.44: **114% of the Roche lobe**. Past L1 the radii diverge toward
        # the Lagrange point rather than failing — measured max |r| = 2.09 mas for rpole = 0.5
        # — so the two components merged into one blob, and the degenerate geometry reached
        # Makie as `(NaN, NaN)` axis limits, whose exception escapes `QML.julia_call` and
        # freezes the window.
        #
        # After this the Orbit tab IS the owner and `_sync_orbit_from_tab!` is correct; it just
        # must not be the one to decide what a brand-new Roche pair looks like.
        _seed_orbit_from_model!(sh, m)
        console!(sh, "model $(m.name): binary, two roche components")
        sh.status = "model $(m.name) (binary)"
        refresh_model_tab!(sh)
        return sh.status
    end
    m = add_model!(sh.session, st)
    log!(sh.session, "$(m.name) = " * _namedtuple_literal(star_params(m));
         note = "model $(m.name)", binding = m.name)
    console!(sh, "model $(m.name): surface_type $(m.surface_type) at schema defaults")
    sh.status = "model $(m.name)"
    refresh_model_tab!(sh)
    return sh.status
end

# FOUR fields: the numeric `surface_type` stays, because it is what every branch in the
# geometry code compares against and what the command log has to write, and the LABEL is added
# beside it for the header to show. A reader comparing models wants "Rapid rotator", not the
# integer `compute_radii` dispatches on.
shell_models() =
    join(("$(m.name)\t$(m.surface_type)\t$(length(m.free))\t$(_fit_type_label(m.surface_type))"
          for m in _sh().session.models), "\n")

"""
    shell_clear_model() -> String

Drop the current model. The surface views go idle and the epoch table falls back to point
counts, which is the honest state when there is nothing to compare the data against.
"""
function shell_clear_model()
    sh = _sh()
    isempty(sh.session.models) && return sh.status
    nm = sh.session.models[end].name
    empty!(sh.session.models)
    sh.session.current_model = 0
    sh.chi2key[] = nothing
    sh.lastfit = ""
    filter!(e -> !startswith(e.note, "model "), sh.session.log)
    console!(sh, "cleared model $(nm)")
    sh.status = "no model"
    refresh_both!(sh)
    return sh.status
end

function shell_select_model(i)
    sh = _sh()
    n = Int(i)
    n in eachindex(sh.session.models) || return sh.status
    sh.session.current_model = n
    refresh_model_tab!(sh)
    return sh.status
end

# The ORBIT group is dropped once there is a companion: it belongs to the PAIR, and a copy in
# each component frame is two editable copies of one physical orbit that can disagree. See
# `shell_binary_orbit_params`, which shows it once, in the Positions frame.
"""
    shell_params() -> String

The parameter form, one row per field:

    name\tlabel\tunit\tvalue\tstate\tlo\thi\ttie\tgroup\tkind\tchoices\tdoc

`state` is `free` / `fixed` / `tied` — three-way, because a tick box cannot express
"expression". `choices` is `value=meaning` pairs joined by `|`, empty for a continuous field,
and `kind` says which of the two the form should render.

Rows a given `ldtype` does not use are still listed, with `state` forced to `fixed`: hiding
them would make the form jump around as the law is changed, and a value that is never read
must not be left floating in a fit.

`ldtype = 0` IS THE EXCEPTION, and it is a different case rather than a harder version of the
same one. There is no law at all then, so `ld1`–`ld4` are not coefficients the law happens not
to read — they are coefficients of nothing, and the μ dependence comes from an intensity
provider instead. Those four rows are therefore dropped.

`INTENSITY_PARAMS` LEADS THE TABLE, under its own `group`, for every `ldtype`: converting a
temperature to an emitted intensity is a layer on top of any temperature map, independent of
whether a μ law is in force. See `_params_table`.
"""
function shell_params()
    sh = _sh()
    m = current_model(sh.session)
    m === nothing && return ""
    # SYNC FIRST, so the rows show the value that will actually be drawn rather than the one
    # the next build is about to overwrite.
    _sync_orbit_from_tab!(sh, m)
    return _params_table(m; orbit = m.companion === nothing, frozen = _derived_names(m))
end

"""
    _fmt_param(v) -> String

One parameter value, for a text field a person reads and edits.

`repr` of a Float64 gives every digit it takes to round-trip the binary value — a fitted
radius came back as `1.6007181518514104`, seventeen characters in a box that fits about ten,
so the field scrolled and showed only the tail. The digits were not information either: the
mesh and the polygon transform are evaluated in the model's own precision, Float32 by default,
which carries about 7 significant decimal digits. Everything past that is the printing of a
binary fraction, not a measurement.

Eight significant digits, with `%g` so a whole number stays whole (`5000`, not `5000.0000`)
and a small one does not become an exponent. Wide enough to round-trip Float32 exactly, short
enough to fit.
"""
function _fmt_param(v)
    x = Float64(v)
    isfinite(x) || return string(x)
    s = Printf.@sprintf("%.8g", x)
    # `%g` can still produce an exponent for genuinely large or tiny values (T0 is a JD near
    # 2.45e6, and a bound can be 1e-4); that is correct and short, so it is left alone.
    return s
end

"""
    _params_table(m) -> String

The form rows for ONE model. Shared by the primary and, for a binary, the companion — so the
two forms cannot drift apart, and a column added here appears in both.

`orbit = false` leaves the orbit group out, which is what both component frames of a BINARY
ask for: those elements describe the pair, not either star, and are shown once in the Positions
frame instead (`shell_binary_orbit_params`).

`frozen` names rows that are DERIVED rather than editable — the orientation a component takes
from the orbit that places it (see `_derived_names`). They are still listed, with their value
and a doc saying where it comes from, because hiding them would make the form jump around and
would hide the answer to "what inclination is this actually being drawn at".
"""
function _params_table(m; orbit::Bool = true, frozen = Set{Symbol}())
    specs = surface_params(m.surface_type)
    ldt = round(Int, get(m.params, :ldtype, 3.0))
    used = ld_coefficients_used(ldt)
    rows = String[]
    # ONE row builder, so the provider fields and the surface fields cannot disagree
    # about the column order or the escaping.
    row(name, label, unit, value, state, lo, hi, tie, group, kind, ch, doc, inert) =
        push!(rows, join((name, label, unit, value, state, lo, hi, tie, group, kind, ch,
                          replace(doc, "\n" => " ", "\t" => " "),
                          inert ? "1" : "0"), "\t"))
    # ── the intensity provider FIRST, for every law ─────────────────────────────────────
    #
    # TEMPERATURE → INTENSITY IS A LAYER ON TOP OF ANY TEMPERATURE MAP, not a property of the
    # limb-darkening choice: whatever the surface's temperatures are, something has to say what
    # brightness each one emits at the observing wavelength, and that question is the same
    # whether or not a μ law is in force. So the provider rows lead the form and are shown under
    # every `ldtype`, rather than appearing at the bottom only when the law is switched off.
    #
    # APPENDED TO THE SAME TABLE rather than published as a table of its own. The three
    # existing tables (`shell_params`, `shell_position_params`, `shell_binary_orbit_params`)
    # each cost a ListModel, an append loop and a `pcomp` branch in all four setters; these
    # rows are per-component exactly as the surface fields are, so they travel with them and
    # the companion's form gets them for nothing. `pgroup` is what sections them under their
    # own heading, which is the mechanism the form already uses for geometry and thermal.
    if any(ps -> ps.name === :ldtype, specs)
        isrc = round(Int, get(m.params, :intensity_source, INTENSITY_PARAMS[1].default))
        iused = intensity_params_used(isrc, m.params)
        for ps in INTENSITY_PARAMS
            # ONLY WHAT THE SOURCE READS, plus the selector. Unlike a limb-darkening
            # coefficient — which still belongs to the model under a law that ignores it — a
            # grid's Teff window means nothing at all to Planck or Rayleigh–Jeans, so greying
            # it would be showing a setting of something that is not there. The rows come and
            # go with the source, and with the grid choice (a path only for a custom grid).
            (ps.name === :intensity_source || ps.name in iused) || continue
            v = ps.kind === :path ? m.grid_path :
                _fmt_param(get(m.params, ps.name, ps.default))
            ch = join(("$(k)=$(vv)" for (k, vv) in ps.choices), "|")
            # No state, no tie and no bounds columns worth having: none of these is a
            # coordinate a fit can move (`:choice`, `:path` and `:setting` all say so), and
            # `_set_state!` refuses one if the form ever offers it.
            row(ps.name, ps.label, ps.unit, v, "fixed",
                _fmt_param(ps.lo), _fmt_param(ps.hi), "", ps.group, ps.kind, ch, ps.doc, false)
        end
    end
    for ps in specs
        (orbit || ps.group !== :orbit) || continue
        # LD COEFFICIENTS DISAPPEAR AT ldtype = 0, rather than being greyed like the ones a
        # law simply does not read. The difference is real: under laws 1-4 the unused
        # coefficients still BELONG to the model — showing them greyed is how the form says
        # which ones the law in force uses, and switching the law brings them back with their
        # values. At ldtype = 0 there is no law at all, so there is nothing left for a
        # coefficient to be a coefficient OF: the μ dependence comes from the intensity
        # provider, whose own fields lead the form. `ld_coefficients_used(0) === ()`
        # is the test, which is why it is spelled as `isempty(used)` and not as `ldt == 0`.
        (ps.name in (:ld1, :ld2, :ld3, :ld4) && isempty(used)) && continue
        v  = get(m.params, ps.name, ps.default)
        lo, hi = get(m.bounds, ps.name, (ps.lo, ps.hi))
        tie = get(m.ties, ps.name, "")
        inert = (ps.group === :limbdark && ps.name in (:ld1, :ld2, :ld3, :ld4) &&
                 !(ps.name in used)) || ps.name in frozen
        state = !isempty(tie) ? "tied" :
                (inert ? "fixed" : (ps.name in m.free ? "free" : "fixed"))
        ch = join(("$(k)=$(vv)" for (k, vv) in ps.choices), "|")
        # WHY A ROW IS GREYED, said here. The panel used to carry the one literal it had a
        # use for ("not used by the current LD law") and applied it to every inert row — so a
        # frozen orientation, which is inert for a completely different reason, claimed to be
        # an unused limb-darkening coefficient. The side that decides is the side that
        # explains.
        doc = ps.name in frozen ?
              "from the orbit, which is what places the secondary: " *
              "inclination = 180 - i, position angle = Ω - 180. " * ps.doc :
              inert ? "not used by the current limb-darkening law. " * ps.doc : ps.doc
        row(ps.name, ps.label, ps.unit, _fmt_param(v), state, _fmt_param(lo), _fmt_param(hi),
            tie, ps.group, ps.kind, ch, doc, inert)
    end
    return join(rows, "\n")
end

# Where a new secondary starts, in `(West, North, toward-observer)` mas.
#
# NOT the origin: two components at the same place are one picture of one star, the χ² of a
# single object with twice the flux, and nothing on screen to say the second exists. A visible
# separation is the state you can then move, and 10 mas is wide for a resolved binary — which
# is the point, since it is easier to bring two stars together than to discover one hiding
# inside the other.
#
# `:offset` rather than `:orbit` to begin with, for the same reason: it is self-contained. An
# orbit needs elements set up on another tab before it places anything, and a binary that
# draws nothing until you visit a second tab looks broken.
const DEFAULT_COMPANION_OFFSET = (10.0, 10.0, 0.0)

"""
    shell_binary() -> String

`"1"` when the current model is a binary, `"0"` otherwise, `""` with no model.
"""
function shell_binary()
    m = current_model(_sh().session)
    m === nothing && return ""
    return m.companion === nothing ? "0" : "1"
end

"""
    shell_companion_surface_type() -> String

The companion's surface type as a code, or `""` when there is no companion.

The panel needs this to SHOW what the companion is, not merely to set it: the combo beside
"Secondary" had no way to read the current value and so sat blank, even on a binary whose
companion had a perfectly good type. `shell_companion_type` is the setter; this is the missing
getter.
"""
function shell_companion_surface_type()
    m = current_model(_sh().session)
    (m === nothing || m.companion === nothing) && return ""
    # The FIELD, not `params[:surface_type]` — the constructor deliberately excludes that key
    # from the parameter dict (`if k !== :surface_type`), so reading it there silently returned
    # the fallback and the combo reported a Roche companion whatever it actually was.
    return string(m.companion.surface_type)
end

"""
    shell_set_binary(on, surface_type) -> String

Give the current model a companion, or take it away.

The companion is a full model of its own — its own surface type, parameters, bounds and ties —
because the two stars of a binary are rarely the same thing: Spica is two limb-darkened
spheroids of very different size, β Lyr is a star and a disc-shrouded companion. What they
SHARE is the orbit, and that lives on the Orbit tab rather than being copied into either.

Defaults to a Roche component, since that is the surface type whose whole reason for existing
is that it sits in a binary.
"""
function shell_set_binary(on, surface_type = 3)
    sh = _sh()
    m = current_model(sh.session)
    m === nothing && return "no model"
    if _qmlstr(on) == "1"
        if m.companion === nothing
            st = round(Int, Float64(surface_type))
            st in SURFACE_TYPE_ORDER || return "unknown surface type $(st)"
            spec = surface_spec(st)
            pr = default_star_params(st)
            m.companion = ModelEntry(
                m.name * "_B", st,
                Dict{Symbol,Float64}(k => Float64(v) for (k, v) in pairs(pr)
                                     if k !== :surface_type),
                Set{Symbol}(),
                Dict{Symbol,Tuple{Float64,Float64}}(ps.name => (ps.lo, ps.hi)
                                                    for ps in surface_params(st)),
                Dict{Symbol,String}(), true, nothing, :offset, DEFAULT_COMPANION_OFFSET,
                # The companion starts with no grid of its own. Not inherited from the
                # primary: the two stars of a binary are rarely the same kind of thing — one
                # grid rarely spans both — which is the whole reason the companion is a
                # component with its own parameters rather than a set of numbers on the
                # primary.
                "")
            # The orbit is the PAIR's, so the companion starts from the primary's copy of it
            # rather than from the schema defaults — and with `q` INVERTED, which the Roche
            # potential requires and which a schema default cannot express.
            _sync_shared_orbit!(m)
            console!(sh, "binary: added a $(spec.name) companion")
        end
    else
        m.companion === nothing || console!(sh, "binary: companion removed")
        m.companion = nothing
    end
    sh.chi2key[] = nothing
    refresh_both!(sh)
    return m.companion === nothing ? "single star" : "binary"
end

# The companion's form is the SAME table as the primary's, read from the companion instead —
# so the QML that draws one row draws either, and a column added to one appears in both.
_companion(sh) = (m = current_model(sh.session); m === nothing ? nothing : m.companion)

"""
    shell_binary_placement() -> String

`mode\tx\ty\tz` — how the secondary is placed, and the offset if it is a fixed one.
Empty when the model is not a binary.
"""
function shell_binary_placement()
    c = _companion(_sh())
    c === nothing && return ""
    return join((c.place, _fmt_param(c.offset[1]), _fmt_param(c.offset[2]),
                 _fmt_param(c.offset[3])), "\t")
end

"""
    shell_set_binary_placement(mode, x, y, z) -> String

Place the secondary by the ORBIT or at a fixed offset in mas.

`:orbit` takes the elements from the Orbit tab and gives a different separation at every epoch.
`:offset` is one displacement throughout — what a snapshot supports, where there is no orbital
motion in the data to constrain elements with. The offset is `(West, North, toward-observer)`;
`z` does not enter the interferometric phase, which sees only the sky projection, but it
decides which component is in front for the 3-D view and for occultation.
"""
function shell_set_binary_placement(mode)
    sh = _sh()
    c = _companion(sh)
    c === nothing && return "not a binary"
    md = _qmlname(mode)
    md === nothing && return "placement must be orbit or offset"
    md in (:orbit, :offset) || return "placement must be orbit or offset"
    # THE MODE ONLY. The offset itself belongs to the position parameters and is written by
    # `shell_set_position_param`; this used to take x, y and z as well, and since the panel
    # stopped passing them once they became parameters, their defaults of "0" were silently
    # overwriting the separation. Switching to the orbit and back therefore moved the secondary
    # onto the primary, where a binary looks exactly like one star.
    #
    # Keeping it here is also wrong in principle: a placement is a choice between two sources,
    # and switching away from one must not destroy the other's numbers — you have to be able to
    # compare them.
    c.place = md
    sh.chi2key[] = nothing
    refresh_both!(sh)
    return md === :offset ?
        Printf.@sprintf("secondary at (%.4g, %.4g) mas", c.offset[1], c.offset[2]) :
        "secondary placed by the orbit"
end

"""
    shell_params2() -> String

[`shell_params`](@ref) for the companion; empty when the model is not a binary.
"""
function shell_params2()
    sh = _sh()
    m = current_model(sh.session)
    (m === nothing || m.companion === nothing) && return ""
    _sync_orbit_from_tab!(sh, m)
    return _params_table(m.companion; orbit = false, frozen = _derived_names(m))
end

"""
    shell_set_param2(name, value) -> String

[`shell_set_param`](@ref) for the companion.
"""
function shell_set_param2(name, value)
    sh = _sh()
    c = _companion(sh)
    c === nothing && return "not a binary"
    n = _qmlname(name)
    n === nothing && return NO_QML_NAME
    n in _derived_names(current_model(sh.session)) && return ORBIT_ORIENTS
    if _set_text_param!(c, n, value)
        refresh_both!(sh)
        return ""
    end
    v = _qmlreal(value, NaN)
    isnan(v) && return "not a number: $(value)"
    c.params[n] = v
    n === :ldtype && _reset_unused_ld!(c)
    n === :intensity_source && _note_intensity_source(sh, v)
    sh.chi2key[] = nothing
    refresh_both!(sh)
    return ""
end

"""
    shell_set_param_state2(name, state) -> String
    shell_set_bound2(name, lo, hi) -> String
    shell_set_tie2(name, expr) -> String

The companion's state, bounds and ties — [`shell_set_param_state`](@ref) and friends, on the
secondary.

The secondary is a component in its own right, not a set of numbers attached to the primary: it
has its own surface type, its own free set and its own bounds, because the two stars of a
binary are rarely the same kind of thing. So it gets the same four setters and the same form.
"""
function shell_set_param_state2(name, state)
    sh = _sh(); c = _companion(sh)
    c === nothing && return "not a binary"
    n = _qmlname(name)
    n === nothing && return NO_QML_NAME
    n in _derived_names(current_model(sh.session)) && return ORBIT_ORIENTS
    return _set_state!(sh, c, n, _qmlstr(state))
end

function shell_set_bound2(name, lo, hi)
    c = _companion(_sh())
    c === nothing && return "not a binary"
    l = _qmlreal(lo, NaN); h = _qmlreal(hi, NaN)
    (isnan(l) || isnan(h)) && return "bounds must be numbers"
    l < h || return "lower bound must be below upper"
    n = _qmlname(name)
    n === nothing && return NO_QML_NAME
    c.bounds[n] = (l, h)
    return ""
end

function shell_set_tie2(name, expr)
    sh = _sh(); c = _companion(sh)
    c === nothing && return "not a binary"
    n = _qmlname(name)
    n === nothing && return NO_QML_NAME
    e = strip(_qmlstr(expr))
    isempty(e) ? delete!(c.ties, n) : (c.ties[n] = String(e))
    sh.chi2key[] = nothing
    refresh_both!(sh)
    return ""
end

# ── the secondary's POSITION, as parameters rather than as three text fields ─────────────
#
# Free/fixed/tied and bounds, like every other parameter, because that is what they are: for a
# snapshot binary the separation IS the model, and `x` and `y` are the two numbers a fit would
# move. They are kept apart from the surface parameters — their own frame, their own table —
# because they describe where a component is rather than what it is made of, and mixing the two
# is the same frame confusion the Orbit tab is laid out to avoid.
#
# Stored on the companion under names that cannot collide with a schema field, so they share
# `free`, `bounds` and `ties` with the surface parameters and need no second state machine.
const POSITION_PARAMS = ((:pos_x, "x (W)", 1), (:pos_y, "y (N)", 2), (:pos_z, "z (obs)", 3))

# How far a fit may move the secondary, in mas, unless the bounds are widened by hand.
#
# ±100 mas rather than ±1000: a long-baseline interferometer's field of view is tens of
# milliarcseconds, and beyond it the phase term oscillates faster than the uv sampling can
# follow — so a far-away companion contributes nothing and the fit is free to park it there.
# MEASURED on iota Peg with ±1000: Nelder-Mead shrank the primary to 0.008 mas and put the
# secondary at 958 mas, which fits 224 points better than the starting guess by making one
# component disappear. A degenerate minimum reached legitimately, and a bound is the only
# thing that rules it out.
const POSITION_BOUND = (-100.0, 100.0)

"""
    shell_position_params() -> String

The secondary's position, in the same 13 columns [`shell_params`](@ref) uses — so the form
draws it with the same delegate, and a column added there appears here too.

Empty when the model is not a binary. Rows are marked INERT when the orbit places the
secondary, since then these numbers are not what decides where it is.
"""
function shell_position_params()
    c = _companion(_sh())
    c === nothing && return ""
    inert = c.place !== :offset
    rows = String[]
    for (n, label, i) in POSITION_PARAMS
        v = c.offset[i]
        lo, hi = get(c.bounds, n, POSITION_BOUND)
        tie = get(c.ties, n, "")
        state = !isempty(tie) ? "tied" : (inert ? "fixed" : (n in c.free ? "free" : "fixed"))
        doc = n === :pos_z ?
              "toward the observer. Not an interferometric observable — the phase sees only " *
              "the sky projection — but it decides which component is in front" :
              "offset from the primary, which is the origin"
        push!(rows, join((n, label, "mas", _fmt_param(v), state, _fmt_param(lo),
                          _fmt_param(hi), tie, "position", "number", "", doc,
                          inert ? "1" : "0"), "\t"))
    end
    return join(rows, "\n")
end

"Set one position component, in mas."
function shell_set_position_param(name, value)
    sh = _sh(); c = _companion(sh)
    c === nothing && return "not a binary"
    v = _qmlreal(value, NaN)
    isnan(v) && return "not a number: $(value)"
    n = _qmlname(name)
    n === nothing && return NO_QML_NAME
    i = findfirst(q -> q[1] === n, POSITION_PARAMS)
    i === nothing && return "unknown position parameter $(name)"
    o = collect(c.offset); o[POSITION_PARAMS[i][3]] = v
    c.offset = (o[1], o[2], o[3])
    sh.chi2key[] = nothing
    refresh_both!(sh)
    return ""
end

"Free, fix or tie one position component."
function shell_set_position_state(name, state)
    sh = _sh(); c = _companion(sh)
    c === nothing && return "not a binary"
    c.place === :offset ||
        return "the orbit places the secondary; switch to a fixed offset to free these"
    n = _qmlname(name)
    n === nothing && return NO_QML_NAME
    return _set_state!(sh, c, n, _qmlstr(state))
end

"Bounds for one position component."
shell_set_position_bound(name, lo, hi) = shell_set_bound2(name, lo, hi)

# ── the ORBIT, once, for the PAIR ────────────────────────────────────────────────────────
#
# A Roche component's SHAPE depends on the instantaneous separation, so the orbital elements
# are part of its surface definition and the schema lists them on surface type 3. In a binary
# that puts them on BOTH components, and the Primary and Secondary frames each showed an
# editable copy of one physical orbit — two sets of numbers, side by side, free to disagree.
# They are shown once here instead, in the frame that already answers "where is the secondary".
#
# WHO OWNS WHAT is the only part that needs thought:
#
#   * The ELEMENTS (P, a, e, T0, i, Ω, ω, Ṗ, ω̇) are the Orbit tab's whenever the ORBIT is what
#     places the secondary — that is what "secondary by the orbit" means — so they are shown
#     read-only here and copied into both components, which is what makes the SHAPE agree with
#     the position. Under a fixed offset nothing else owns them: a Roche shape still needs a
#     separation, and these rows are then the only place to give it one, so they are editable.
#   * `q` and `d` are never elements. The Orbit tab's own header says why: a relative
#     astrometric orbit does not constrain a mass ratio. They stay here, editable, always.
#
# `q` is also the one quantity whose NUMBER differs between the components. The Roche potential
# reads M_companion/M_self, which is `q` on the primary and `1/q` on the secondary
# (`update_roche_radii`'s own note, and `demos/spica_binary_roche.jl` writes `1.0/q_binary` by
# hand). The row is in the PRIMARY's convention — M₂/M₁, as published solutions quote it — and
# the secondary's copy is written inverted. It was not: a default binary had q = 0.5 on both
# components, which is a different system on each side of the same model.
const ORBIT_PARAM_NAMES = Set{Symbol}(ps.name for ps in surface_params(3)
                                      if ps.group === :orbit)

# Schema name -> the Orbit tab's name for the same element. `q` and `d` are absent on purpose:
# they are not elements. The tab spells three of these ASCII (`Omega`, `omega`, `domega`) while
# the geometry code reads the Unicode fields, so the mapping cannot be an identity.
const ORBIT_ELEMENT_OF = Dict{Symbol,Symbol}(
    :P => :P, :a => :a, :e => :e, :T0 => :T0, :i => :i,
    :Ω => :Omega, :ω => :omega, :dP => :dP, :dω => :domega)

"Both components of a binary, primary first; just the primary when there is no companion."
_components(m) = m.companion === nothing ? (m,) : (m, m.companion)

"The shared orbit rows to show, in schema order, from whichever component carries them."
function _orbit_specs(m)
    out = similar(surface_params(3), 0)
    seen = Set{Symbol}()
    for c in _components(m), ps in surface_params(c.surface_type)
        (ps.group === :orbit && !(ps.name in seen)) || continue
        push!(seen, ps.name); push!(out, ps)
    end
    return out
end

"Which component carries the state, bounds and tie of a shared orbit parameter."
_orbit_owner(m, n::Symbol) =
    haskey(m.params, n) ? m :
    (m.companion !== nothing && haskey(m.companion.params, n)) ? m.companion : nothing

"One shared orbit parameter, in the primary's convention."
function _orbit_value(m, n::Symbol)
    haskey(m.params, n) && return m.params[n]
    c = m.companion
    (c !== nothing && haskey(c.params, n)) || return NaN
    return n === :q ? 1 / c.params[n] : c.params[n]
end

"Write one shared orbit parameter into every component that carries it, each in its own sense."
function _put_orbit!(m, n::Symbol, v::Real)
    haskey(m.params, n) && (m.params[n] = Float64(v))
    c = m.companion
    if c !== nothing && haskey(c.params, n)
        c.params[n] = n === :q ? 1 / Float64(v) : Float64(v)
    end
    return m
end

"""
    _sync_shared_orbit!(m) -> m

Push the primary's resolved orbital elements onto the companion, `q` inverted.

Called before the geometry is built rather than only on edit, so it also covers the two paths
that do not go through the panel: a TIED element, whose value is not in `params` at all until
`apply_model_ties` runs, and a FIT, which writes trial values straight into the owner. Without
it a free shared `a` moved the primary's shape and left the secondary's at its old separation —
the same orbit disagreeing with itself, one component at a time.

The companion's own free set and ties for these names are dropped: its copy is DERIVED, and a
derived value that is independently free is a second coordinate for one physical quantity.
"""
function _sync_shared_orbit!(m)
    c = m.companion
    c === nothing && return m
    v = apply_model_ties(m)
    for n in ORBIT_PARAM_NAMES
        (haskey(v, n) && haskey(c.params, n)) || continue
        c.params[n] = n === :q ? 1 / v[n] : v[n]
        delete!(c.free, n); delete!(c.ties, n)
    end
    return m
end

"""
    _sync_orbit_from_tab!(sh, m) -> m

Copy the Orbit tab's elements into a binary whose secondary is placed BY that orbit.

With `place === :orbit` the separation at each epoch comes from the Orbit tab, while the Roche
SHAPE is computed from the component's own elements — so two orbits decided the two halves of
one picture, and the star could be shaped by a 10-day period while being moved on a 4-day one.
One orbit: the tab's. `q` and `d` are left alone, being star-model quantities rather than
elements (`orbit_bparams` reads `q` off the primary for exactly that reason).
"""
function _sync_orbit_from_tab!(sh::ShellState, m)
    c = m.companion
    (c === nothing || c.place !== :orbit) && return m
    v = apply_orbit_ties(sh.orbit)
    for (n, e) in ORBIT_ELEMENT_OF
        haskey(v, e) || continue
        for comp in _components(m)
            haskey(comp.params, n) && (comp.params[n] = Float64(v[e]))
        end
    end
    # AND THE ORIENTATION. Both components are oriented by the ORBIT, not by themselves:
    # `create_binary_geometry` takes the spin axis from the orbital inclination and node, and
    # `orbit_bparams` has applied exactly this conversion since the Orbit tab was written —
    #
    #     inclination = 180 - i,   position_angle = Ω - 180
    #
    # while the Model tab passed each component's OWN `inclination`/`position_angle` to
    # `create_star_multiepochs`. So the Orbit plot and the Model tab oriented the same pair
    # differently; `graticule_segments` records the same trap from the other side, measuring
    # the two answers 102.6 degrees apart on a Spica-like binary.
    #
    # Derived, therefore, and not independent: the panel shows these two rows inert (see
    # `_derived_names`) and `binary_fit_names` drops them, because a value this overwrites on
    # the next build must not also be a coordinate a fit can move.
    inc = 180.0 - Float64(get(v, :i, 90.0))
    pa  = Float64(get(v, :Omega, 0.0)) - 180.0
    for comp in _components(m), (n, val) in ((:inclination, inc), (:position_angle, pa))
        haskey(comp.params, n) || continue
        comp.params[n] = val
        delete!(comp.free, n); delete!(comp.ties, n)
    end
    return m
end

"""
    ORBIT_ORIENTED

The component parameters an ORBIT owns, when the orbit is what places the secondary.

Not the orbital elements — those are `ORBIT_ELEMENT_OF`. These are the two fields of each
COMPONENT that `create_binary_geometry` derives from the orbit and would overwrite anyway.
"""
const ORBIT_ORIENTED = (:inclination, :position_angle)

"""
    _seed_orbit_from_model!(sh, m) -> sh

Copy a model's orbital elements INTO the Orbit tab — the reverse of
[`_sync_orbit_from_tab!`](@ref), and used once, when a Roche binary is created.

The two sets of defaults are independent and only one can win. A component's come from the
surface schema and are chosen to sit comfortably inside the Roche lobe; the Orbit tab's
describe a generic pair of discs. For a model that is being built right now, the model's are
the ones that were asked for.
"""
function _seed_orbit_from_model!(sh::ShellState, m)
    for (n, e) in ORBIT_ELEMENT_OF
        haskey(m.params, n) || continue
        sh.orbit.params[e] = Float64(m.params[n])
        # A tie on an element would immediately overwrite what was just seeded.
        delete!(sh.orbit.ties, e)
    end
    # `q` is not an element (a relative astrometric orbit does not constrain a mass ratio), but
    # the tab keeps one for rendering when there is no model; keep it in step anyway.
    haskey(m.params, :q) && (sh.orbit.q = Float64(m.params[:q]))
    refresh_orbit!(sh)
    return sh
end

"""
    _derived_names(m) -> Set{Symbol}

The parameter names a component frame must show INERT, because something else owns them.

Empty for a single star and for a binary placed by a fixed offset — there, a component's
orientation is its own. Non-empty only while the ORBIT places the secondary, which is exactly
when `_sync_orbit_from_tab!` writes these from the elements.
"""
_derived_names(m) =
    (m.companion !== nothing && m.companion.place === :orbit) ?
    Set{Symbol}(ORBIT_ORIENTED) : Set{Symbol}()

"""
    shell_binary_orbit_params() -> String

The PAIR's orbit, in the same 13 columns [`shell_params`](@ref) uses — so the Positions frame
draws it with the same delegate as everything else.

Empty unless the model is a binary whose components carry orbital elements at all: two spheres
have no orbit to shape them, and for them this frame is the offset and nothing more. Element
rows are marked INERT while the orbit places the secondary, because then the Orbit tab owns
them.
"""
function shell_binary_orbit_params()
    sh = _sh()
    m = current_model(sh.session)
    (m === nothing || m.companion === nothing) && return ""
    _sync_orbit_from_tab!(sh, m)
    _sync_shared_orbit!(m)
    byorbit = m.companion.place === :orbit
    rows = String[]
    for ps in _orbit_specs(m)
        own = _orbit_owner(m, ps.name)
        own === nothing && continue
        inert = byorbit && haskey(ORBIT_ELEMENT_OF, ps.name)
        lo, hi = get(own.bounds, ps.name, (ps.lo, ps.hi))
        tie = get(own.ties, ps.name, "")
        state = !isempty(tie) ? "tied" :
                (inert ? "fixed" : (ps.name in own.free ? "free" : "fixed"))
        doc = inert ? "from the Orbit tab, which is what places the secondary. " * ps.doc :
              ps.doc
        push!(rows, join((ps.name, ps.label, ps.unit, _fmt_param(_orbit_value(m, ps.name)),
                          state, _fmt_param(lo), _fmt_param(hi), tie, "orbit", ps.kind,
                          join(("$(k)=$(vv)" for (k, vv) in ps.choices), "|"),
                          replace(doc, "\n" => " ", "\t" => " "),
                          inert ? "1" : "0"), "\t"))
    end
    return join(rows, "\n")
end

"The message an orientation edit gets back while the orbit is what decides it."
const ORBIT_ORIENTS = "the orbit orients both components (inclination = 180 - i, " *
                      "position angle = Ω - 180); change i or Ω, or switch to a fixed offset"

"The message an element edit gets back while the Orbit tab is the one that owns it."
const ORBIT_TAB_OWNS = "the orbit places the secondary; set its elements on the Orbit tab, " *
                       "or switch to a fixed offset"

"Set one shared orbital parameter, in every component that carries it."
function shell_set_binary_orbit_param(name, value)
    sh = _sh()
    m = current_model(sh.session)
    (m === nothing || m.companion === nothing) && return "not a binary"
    n = _qmlname(name)
    n === nothing && return NO_QML_NAME
    n in ORBIT_PARAM_NAMES || return "unknown orbital parameter $(n)"
    haskey(ORBIT_ELEMENT_OF, n) && m.companion.place === :orbit && return ORBIT_TAB_OWNS
    v = _qmlreal(value, NaN)
    isnan(v) && return "not a number: $(value)"
    # The secondary's copy is `1/q`, so zero is not a value this can carry.
    n === :q && v == 0 && return "mass ratio cannot be zero"
    _put_orbit!(m, n, v)
    # AND BACK TO THE ORBIT TAB, so the two tabs cannot describe different systems. The sync
    # already ran the other way — `_sync_orbit_from_tab!` pushes the tab's elements into both
    # components whenever the orbit is what places the secondary — but under a FIXED OFFSET
    # nothing owned these, so an element edited here left the Orbit tab drawing the old one.
    # One orbit, edited from either end.
    #
    # `q` and `d` do not propagate: they are not elements. The Orbit tab's `q` is a rendering
    # assumption used only when there is no model to take it from, and `orbit_bparams` already
    # prefers the model's.
    if haskey(ORBIT_ELEMENT_OF, n)
        sh.orbit.params[ORBIT_ELEMENT_OF[n]] = v
        refresh_orbit!(sh)
    end
    sh.chi2key[] = nothing
    refresh_both!(sh)
    return ""
end

"Free, fix or tie one shared orbital parameter."
function shell_set_binary_orbit_state(name, state)
    sh = _sh()
    m = current_model(sh.session)
    (m === nothing || m.companion === nothing) && return "not a binary"
    n = _qmlname(name)
    n === nothing && return NO_QML_NAME
    own = _orbit_owner(m, n)
    own === nothing && return "unknown orbital parameter $(n)"
    haskey(ORBIT_ELEMENT_OF, n) && m.companion.place === :orbit && return ORBIT_TAB_OWNS
    return _set_state!(sh, own, n, _qmlstr(state))
end

"Bounds for one shared orbital parameter."
function shell_set_binary_orbit_bound(name, lo, hi)
    m = current_model(_sh().session)
    (m === nothing || m.companion === nothing) && return "not a binary"
    n = _qmlname(name)
    n === nothing && return NO_QML_NAME
    own = _orbit_owner(m, n)
    own === nothing && return "unknown orbital parameter $(n)"
    l = _qmlreal(lo, NaN); h = _qmlreal(hi, NaN)
    (isnan(l) || isnan(h)) && return "bounds must be numbers"
    l < h || return "lower bound must be below upper"
    own.bounds[n] = (l, h)
    return ""
end

"Tie one shared orbital parameter to an expression."
function shell_set_binary_orbit_tie(name, expr)
    sh = _sh()
    m = current_model(sh.session)
    (m === nothing || m.companion === nothing) && return "not a binary"
    n = _qmlname(name)
    n === nothing && return NO_QML_NAME
    own = _orbit_owner(m, n)
    own === nothing && return "unknown orbital parameter $(n)"
    haskey(ORBIT_ELEMENT_OF, n) && m.companion.place === :orbit && return ORBIT_TAB_OWNS
    e = strip(_qmlstr(expr))
    isempty(e) ? delete!(own.ties, n) : (own.ties[n] = String(e))
    sh.chi2key[] = nothing
    refresh_both!(sh)
    return ""
end

"""
    shell_companion_type(surface_type) -> String

Change the companion's surface type, rebuilding it at that type's defaults.
"""
function shell_companion_type(surface_type)
    sh = _sh()
    m = current_model(sh.session)
    (m === nothing || m.companion === nothing) && return "not a binary"
    m.companion = nothing
    return shell_set_binary("1", surface_type)
end

"""
    _reset_unused_ld!(m) -> m

Zero the limb-darkening coefficients the current law does not read, and release them.

`compute_ldmap` reads `ld1` for the linear and power laws, `ld1`+`ld2` for the quadratic and
all four for Claret. A value left in `ld2` after switching from quadratic to linear is not
harmless: it is invisible in the fit, it reappears the moment the law is switched back, and if
it was FREE it stays in the parameter vector as a coordinate that changes nothing — which is
how an optimiser spends its evaluations on a direction the χ² is flat along.

So the unused ones are set to zero, dropped from `free`, and their ties removed. The panel
greys them, which is the visible half of the same rule.
"""
function _reset_unused_ld!(m)
    used = ld_coefficients_used(round(Int, get(m.params, :ldtype, 3.0)))
    for n in (:ld1, :ld2, :ld3, :ld4)
        n in used && continue
        haskey(m.params, n) && (m.params[n] = 0.0)
        delete!(m.free, n)
        delete!(m.ties, n)
    end
    return m
end

"""
    _qmlreal(x, default) -> Float64

A number arriving from QML, whatever shape QML chose to send it in.

**QML does not always send text.** A `TextField` hands over a string, but a `SpinBox` hands
over its integer `value` and a slider a real — so a callback that assumed text and wrote
`tryparse(Float64, String(x))` threw the moment it was wired to a spin box:

    shell_set_graticule(lat_deg::Int32, lon_deg::Int32, colour::QML.QStringAllocated)
    The type `String` exists, but no method is defined for this combination of argument types

and the exception escaped through `QML.julia_call` and froze the window — the same failure
mode as a `nothing` name (see [`_qmlname`](@ref)). `nothing` is handled here too, for the same
reason it is there: a recycled delegate can call with no value at all.
"""
_qmlreal(x, default::Real) =
    x === nothing ? Float64(default) :
    x isa Real    ? Float64(x) :
    something(tryparse(Float64, String(x)), Float64(default))

"""
    _qmlstr(x, default = "") -> String

Text arriving from QML, or `default` when QML sent nothing. See [`_qmlreal`](@ref).

TOTAL, on every shape QML sends, because the alternative is a frozen window. `String(x)` alone
has no method for `Int32`, `Float64` or `Nothing` — a `CheckBox` bound through a numeric
property, a `Slider`, and a recycled delegate respectively — and each of those is a
`MethodError` that escapes `QML.julia_call`. A number becomes its decimal text, which no flag
comparison will match, so the caller REFUSES instead of freezing.
"""
_qmlstr(x, default::AbstractString = "") =
    x === nothing         ? String(default) :
    x isa AbstractString  ? String(x) :
    x isa Real            ? string(x) :
    String(x)                     # QML's own string type, which has a `String` constructor

"""
    _qmlname(x) -> Symbol or nothing

A parameter name arriving from QML, or `nothing` when QML had none to give.

**QML can hand us `undefined`, and it arrives here as `nothing`.** A ListView recycles its
delegates, so a handler such as `onEditingFinished` — which fires on FOCUS LOSS — can run after
its row's model context is gone; the roles then read `undefined`. `String(nothing)` throws, the
exception escapes through `QML.julia_call`, and the whole window stops responding with nothing
on screen and nothing in the console to say why.

The QML side captures the row's identity at construction so this should not arise, but a
callback that cannot be reached from a dead row is worth guaranteeing on this side too: every
entry point that takes a name from QML goes through here and REFUSES rather than throws.
"""
_qmlname(x) = x === nothing ? nothing : Symbol(_qmlstr(x))

"The refusal every name-taking entry point returns when QML had no name to give."
const NO_QML_NAME = "no parameter named by the form (the row went away before the edit landed)"

"""
    _set_text_param!(m, n, value) -> Bool

Take a TEXT-valued form field, and say whether this was one.

Exactly one field is text rather than a number — `grid_path`, the intensity grid file — and it
arrives through the same `shell_set_param` as everything else, because the form has one value
column and the schema is what says how to render it. The TEST IS THE KIND, not the name: a
second `:path` field would travel this path with nothing to change here.

Stripped, because a path pasted from a terminal or dropped from a file manager carries
whitespace, and a trailing space makes `isfile` false for a file that is plainly there. A
`file://` prefix is removed for the same reason — that is what a drag-and-drop hands over, and
it is what the rest of this file strips at every other path entry point.
"""
function _set_text_param!(m, n::Symbol, value)
    ps = intensity_spec(n)
    (ps === nothing || ps.kind !== :path) && return false
    p = strip(_qmlstr(value))
    startswith(p, "file://") && (p = p[8:end])
    m.grid_path = String(p)
    return true
end

"""
    shell_set_param(name, value) -> String

Set one parameter. Non-numeric text is refused rather than silently zeroing the field, which
is what a half-typed number would otherwise do on every keystroke.
"""
function shell_set_param(name, value)
    sh = _sh()
    m = current_model(sh.session)
    m === nothing && return "no model"
    n = _qmlname(name)
    n === nothing && return NO_QML_NAME
    n in _derived_names(m) && return ORBIT_ORIENTS
    if _set_text_param!(m, n, value)
        refresh_model_tab!(sh)
        return ""
    end
    v = _qmlreal(value, NaN)
    isnan(v) && return "not a number: $(value)"
    m.params[n] = v
    n === :ldtype && _reset_unused_ld!(m)
    n === :intensity_source && _note_intensity_source(sh, v)
    refresh_model_tab!(sh)
    return ""
end

"""
    _note_intensity_source(sh, src)

Say, in the console, WHAT CONSUMES the intensity source that was just chosen.

This is here because the honest answer is "not everything yet", and a selector whose effect a
user has to discover by comparing pictures is worse than one that says so. `ROTIR.build_intensity_provider`
turns these settings into a provider, and `build_parametric_logπ` is the one χ² that takes one
(`provider = …`, src/parametric_gradient.jl) — the panel's own surface build and its χ² readout
go through `compute_ldmap`, which has no provider argument, so at `ldtype = 0` they draw and
score a surface with no limb darkening whichever backend is selected here.

Once per CHANGE, not once per refresh: `refresh_model_tab!` runs on every keystroke in the
form, so a note emitted from the table builder would fill the console.
"""
function _note_intensity_source(sh::ShellState, src)
    i = round(Int, Float64(src))
    j = findfirst(c -> c.first == i, INTENSITY_PARAMS[1].choices)
    name = j === nothing ? "$(i)" : INTENSITY_PARAMS[1].choices[j].second
    console!(sh, "intensity source: $(name)")
    i in (1, 5) && return nothing                 # analytic: nothing to fetch or synthesize
    if i == 3
        pack = ROTIR.kurucz_pack_path()
        console!(sh, isfile(pack) ? "  Kurucz pack: $(pack)" :
                     "  the Kurucz pack (67.5 MB) is downloaded once, on first use, to $(pack)")
    end
    console!(sh, "  note: this is recorded and checked here, and `build_intensity_provider` " *
                 "builds it for a script or for the gradient path's `provider =`; the " *
                 "panel's own surface build and χ² still go through `compute_ldmap`, which " *
                 "takes no provider, so what is DRAWN at ldtype = 0 has no limb darkening")
    return nothing
end

"""
    shell_set_param_state(name, state) -> String

`free`, `fixed`, or `tied` (with the expression set separately by `shell_set_tie`).
"""
function shell_set_param_state(name, state)
    sh = _sh()
    m = current_model(sh.session)
    m === nothing && return "no model"
    n = _qmlname(name)
    n === nothing && return NO_QML_NAME
    n in _derived_names(m) && return ORBIT_ORIENTS
    return _set_state!(sh, m, n, _qmlstr(state))
end

"""
    _set_state!(sh, m, n, s) -> String

Free, fix or tie one parameter of ONE component. Shared by the primary, the secondary and the
position parameters, so all three obey the same refusals rather than three copies of them.
"""
function _set_state!(sh::ShellState, m, n::Symbol, s::AbstractString)
    # A DISCRETE field has no free/fixed/tied to make: `ldtype` selects which limb-darkening
    # law `compute_ldmap` applies, and an optimiser cannot walk 1→2→3→4 as a continuous
    # coordinate — every backend would sample fractional laws that do not exist. The panel
    # hides the selector for these, and this refuses it, so the rule lives with the schema
    # rather than only in the QML that happens to draw it.
    specs = surface_params(m.surface_type)
    i = findfirst(ps -> ps.name === n, specs)
    if i !== nothing && specs[i].kind === :choice && s != "fixed"
        return "$(specs[i].label) is a discrete choice, not a value that can be $(s)"
    end
    # NOR ON A PROVIDER SETTING. `INTENSITY_PARAMS` describes where the intensity comes from
    # (a backend, a file, the window to read out of it), and none of it is a coordinate: a
    # grid's Teff axis is a property of the FILE, so freeing it would be a direction the χ² is
    # not merely flat along but undefined along. The form hides the selector for these kinds;
    # this is the same rule where the schema can enforce it rather than the QML that draws it.
    let ips = intensity_spec(n)
        if ips !== nothing && s != "fixed"
            return "$(ips.label) configures where the intensity comes from; it is not a " *
                   "coordinate and cannot be $(s)"
        end
    end
    # A limb-darkening coefficient the current law does not read is not a coordinate: freeing
    # it would add a direction the χ² is exactly flat along, and tying it would compute a value
    # nothing reads. The panel greys these; this is the same rule where it can be enforced.
    if n in (:ld1, :ld2, :ld3, :ld4) && s != "fixed" &&
       !(n in ld_coefficients_used(round(Int, get(m.params, :ldtype, 3.0))))
        return "$(n) is not used by the current limb-darkening law"
    end
    if s == "free"
        push!(m.free, n); delete!(m.ties, n)
    elseif s == "fixed"
        delete!(m.free, n); delete!(m.ties, n)
    elseif s == "tied"
        delete!(m.free, n)
        haskey(m.ties, n) || (m.ties[n] = "")
    else
        return "unknown state: $s"
    end
    refresh_model_tab!(sh)
    return ""
end

function shell_set_bound(name, lo, hi)
    sh = _sh()
    m = current_model(sh.session)
    m === nothing && return "no model"
    l = _qmlreal(lo, NaN); h = _qmlreal(hi, NaN)
    (isnan(l) || isnan(h)) && return "bounds must be numbers"
    l < h || return "lower bound must be below upper"
    n = _qmlname(name)
    n === nothing && return NO_QML_NAME
    m.bounds[n] = (l, h)
    return ""
end

"""
    shell_set_tie(name, expr) -> String

Set a tie expression and report whether it currently evaluates.

Evaluation failure is reported, not raised: a half-typed expression is the normal state of a
text field somebody is editing.
"""
function shell_set_tie(name, expr)
    sh = _sh()
    m = current_model(sh.session)
    m === nothing && return "no model"
    n = _qmlname(name)
    n === nothing && return NO_QML_NAME
    e = _qmlstr(expr)
    m.ties[n] = e
    delete!(m.free, n)
    v = eval_tie(e, m.params)
    v === nothing && return isempty(strip(e)) ? "" : "does not evaluate yet"
    m.params[n] = v
    refresh_model_tab!(sh)
    return Printf.@sprintf("= %.6g", v)
end

"""
    shell_derived_summary() -> String

[`derived_summary_text`](@ref) for the model the panel is showing. Empty string when there is
nothing to say, which is what makes the QML label hide itself.
"""
function shell_derived_summary()
    sh = _sh()
    m = current_model(sh.session)
    m === nothing && return ""
    return try
        ROTIR.derived_summary_text(star_params(m))
    catch
        # A half-edited model (a cleared field, a period of 0) must not take the panel down.
        ""
    end
end

"""
    shell_validate_model() -> String

Every problem the schema can see with the current model, one per line. Empty means it builds.
"""
function shell_validate_model()
    sh = _sh()
    m = current_model(sh.session)
    m === nothing && return "no model"
    # BOTH lists: the user should see plausibility advice too. Only `validate_star_params`
    # gates whether the model can be built (see `advise_star_params` for why they are split).
    p = star_params(m)
    # AND THE PROVIDER, which `validate_star_params` deliberately knows nothing about: none of
    # `INTENSITY_PARAMS` is a `star_params` field, and that function is the build GATE — a
    # wrong grid path must not stop the panel drawing the surface, it must say so in the same
    # place every other problem is reported. Both components, because each carries its own.
    prov = vcat(_intensity_messages(m, ""),
                m.companion === nothing ? String[] :
                _intensity_messages(m.companion, "secondary: "))
    return join(vcat(validate_star_params(p), advise_star_params(p), prov), "\n")
end

"[`validate_intensity_source`](@ref) for one component, prefixed so a binary says which."
_intensity_messages(m, prefix::AbstractString) =
    [prefix * msg for msg in validate_intensity_source(
        round(Int, get(m.params, :ldtype, 3.0)), m.params, m.grid_path)]

# ── the plots ───────────────────────────────────────────────────────────────────────────

"""
    shell_tessellation() -> String

`kind\tnside_exp\tnpix\tprecision` — what the model is evaluated on.
"""
function shell_tessellation()
    sh = _sh()
    n = sh.nside_exp[]
    return "$(sh.tessel[])\t$(n)\t$(12 * (2^n)^2)\t$(sh.precision[])"
end

"""
    shell_set_tessellation(kind, nside_exp, precision) -> String

Choose the mesh the model is evaluated on, and the float type it is built in.

`nside_exp` is log₂(nside), which is what `tessellation_healpix` takes — and it is clamped to
2–6. Level 3 up is what a result should be quoted at; level 2 (192 tessels) is too coarse to
carry a spot but converges in seconds, which is what makes it the level to check a setup on
before committing to one. Above 6 the polygon FT dominates everything else the GUI does.

`:longlat` is accepted by the geometry but is NOT wired here: `create_star` builds it, and the
regularizers and the shape gradients are all written against the HEALPix neighbour structure,
so offering it would mean a mesh half the panel does not work on.

`precision` is `Float32` by default. It is what `tessellation_healpix` itself defaults to, it
halves the memory the polygon FT needs, and the visibilities are not measured to seven digits
— but the analytic shape fit needs Float64, and switches to it on its own rather than making
that the default for everything.
"""
function shell_set_tessellation(kind, nside_exp, precision)
    sh = _sh()
    k = Symbol(_qmlstr(kind))
    k === :healpix || return "only :healpix is wired; :longlat has no regularizers or gradients"
    sh.tessel[] = k
    sh.nside_exp[] = clamp(Int(nside_exp), 2, 6)
    p = _qmlstr(precision)
    sh.precision[] = p == "Float64" ? Float64 : Float32
    sh.chi2key[] = nothing                    # the mesh changed, so the χ² did
    refresh_both!(sh)
    return "mesh: $(k) n=$(sh.nside_exp[]) ($(12 * (2^sh.nside_exp[])^2) tessels, $(sh.precision[]))"
end

"The tessellation the model is evaluated on, built from the current settings."
model_tessellation(sh::ShellState) =
    tessellation_healpix(sh.nside_exp[]; T = sh.precision[])

"""
    shell_polyft_backends() -> String

`key\tlabel\tdoc` per forward-kernel backend.

THREE, and the ordering is the measured one. The rasterise-plus-FFT route is deliberately not
among them: it is fast, but it disagrees with the exact transform by ~5e-3 and that does not
shrink from nx = 128 to 512 — a pixel-integrated image is not a polygon, and an exact DFT of
the same image reproduces the same 4.97e-3, so it is the model rather than the transform.
"""
shell_polyft_backends() = join((
    "auto\tAutomatic (per fit)\t" *
    "let each code path pick: the type-3 kernel for this panel's χ² and for Nelder-Mead, " *
    "BOBYQA and Nautilus, which can reach it; the exact vectorised kernel for VMLMB, NUTS " *
    "and Pigeons, which cannot — their gradient goes through kernels that have no quadrature " *
    "branch, and it is 16x on the three of them. The first such fit LOADS LoopVectorization, " *
    "which pauses the window for a few seconds, once per session" *
    (ROTIR.turbo_available() ? " (already loaded)" : ""),
    "t3\tQuadrature + ROTIR type 3\t" *
    "the same Gauss quadrature as nufft, through ROTIR's own type-3 transform instead of " *
    "FINUFFT: 0.20 ms at HEALPix 3 and 2.2 ms at 6 against 1.7 and 6.0, and more accurate " *
    "(1.1e-10 against 2.4e-9 on polaris). It is also the only quadrature route with an " *
    "ADJOINT, so unlike nufft a gradient fit can use it — but the SHAPE gradient also needs " *
    "the vertex adjoint, which has no type-3 form and would fall back to the plain reference, " *
    "so leave gradient fits and samplers on turbo",
    "nufft\tQuadrature + NUFFT (fastest)\t" *
    "Gauss-Legendre over each tessel folded into one type-3 NUFFT; 2.6 ms at HEALPix 3 and " *
    "4.1 ms at 5, where the exact kernels take 7 ms and 86 ms — the cost barely moves with " *
    "the mesh. Accurate to 6.8e-7 at HEALPix 3 and 2.5e-9 above it",
    "turbo\tExact, vectorised\t" *
    "the closed-form polygon transform with SIMD transcendentals — no approximation " *
    "parameters at all, 17x the reference at HEALPix 3. Selecting it LOADS " *
    "LoopVectorization, which pauses the window for a few seconds, once per session" *
    (ROTIR.turbo_available() ? " (already loaded)" : ""),
    "scalar\tExact, reference\t" *
    "the same arithmetic in plain Julia; the definition the other two are tested against, " *
    "and what to fall back on if a fit looks wrong"), "\n")

"""
    shell_polyft_backend() -> String

Which forward kernel the panel is SET to — `"auto"`, or the one the user forced.

Not the same question as `ROTIR.polyft_backend()`, which answers what is running right now: a
fit scopes its own kernel for its duration, so during one of them the two differ by design.
"""
shell_polyft_backend() =
    _sh().kernel_auto[] ? "auto" : String(ROTIR.POLYFT_BACKEND[])

"""
    TURBO_LOAD_FAILED

Whether loading LoopVectorization has already been tried in this session and failed. A
PROCESS-level flag rather than a `ShellState` field: it records a property of the installation,
which a new shell does not change.
"""
const TURBO_LOAD_FAILED = Ref(false)

"""
    _ensure_turbo!(sh) -> String

Load LoopVectorization if it is not loaded, so the `:turbo` kernel has methods. Returns `""`
on success and a sentence on failure, which the caller reports and falls back from.

`:turbo` needs LoopVectorization, which ROTIR does NOT load: measured, loading it invalidates
OITOOLS' precompiled canvas code and takes one `build_canvas` from 341 ms to 2685 ms — 1.8 s
onto every GUI start, whether or not this session ever fits anything. So it is loaded HERE,
when something asks for it, and never otherwise.

SYNCHRONOUSLY, ON THE GUI THREAD, which is why this is a separate function and why every
caller of it sits BEFORE `start_job!` rather than inside a worker closure: package loading
takes locks that the render thread also wants, and a `using` on a worker can deadlock against
it. A few seconds of frozen window is the honest cost, and the panel says so before the click.
It is only ever paid once — the next caller finds the extension loaded.
"""
function _ensure_turbo!(sh::ShellState)
    return _ensure_weakdep!(sh, :LoopVectorization, ROTIR.turbo_available, TURBO_LOAD_FAILED,
                            "the turbo kernel",
                            "A bundle cannot carry it: `create_app` builds a multiversioned " *
                            "sysimage while VectorizationBase specialises to the build " *
                            "machine, and the link aborts — see app/Project.toml",
                            "ROTIRLoopVectorizationExt")
end

"""
    _ensure_finufft!(sh) -> String

The same, for FINUFFT and the `:nufft` kernel.

FINUFFT became a weak dependency for two measured reasons: 464 ms on every `using ROTIR`
whether or not anything transforms, and 455 MB of CUDA driver in an application bundle (14 %
of it) by way of `cufinufft_jll`, for a GPU NUFFT nothing calls. `:nufft` is kept as the
independent cross-check against `:t3`, which is the default and both faster and more accurate.
"""
function _ensure_finufft!(sh::ShellState)
    return _ensure_weakdep!(sh, :FINUFFT, ROTIR.finufft_available, FINUFFT_LOAD_FAILED,
                            "the nufft kernel",
                            "A bundle deliberately ships without it: it would add 455 MB of " *
                            "CUDA driver for a GPU transform nothing calls — see " *
                            "ext/ROTIRFINUFFTExt.jl",
                            "ROTIRFINUFFTExt")
end

"""
    FINUFFT_LOAD_FAILED

As [`TURBO_LOAD_FAILED`](@ref), for FINUFFT.
"""
const FINUFFT_LOAD_FAILED = Ref(false)

"""
    _ensure_weakdep!(sh, pkg, available, failed, what, bundle_note, extname) -> String

Load a weak dependency on demand so a kernel that needs it can be selected mid-session.
Returns `""` on success and a sentence on failure, which the caller reports and falls back from.

ONCE PER SESSION, INCLUDING THE FAILURE, which is what `failed` is for: a bundle ships without
these packages by design, and without the flag every fit would retry the load and print the
same sentence again.

SYNCHRONOUSLY, ON THE GUI THREAD, which is why every caller sits BEFORE `start_job!` rather
than inside a worker closure: package loading takes locks the render thread also wants, and a
`using` on a worker can deadlock against it. A few seconds of frozen window is the honest
cost, and the panel says so before the click.
"""
function _ensure_weakdep!(sh::ShellState, pkg::Symbol, available::Function,
                          failed::Base.RefValue{Bool}, what::AbstractString,
                          bundle_note::AbstractString, extname::AbstractString)
    available() && return ""
    failed[] && return "$(what) is not available in this installation ($(pkg) did not load)"
    # ASKED BEFORE IT IS TRIED. `identify_package` resolves against the active project — the
    # same thing the `using` below would consult — so a bundle answers `nothing` here and gets
    # a sentence that says what is wrong instead of a caught `ArgumentError` from `require`.
    if Base.identify_package(String(pkg)) === nothing
        failed[] = true
        return "this installation does not ship $(pkg), so $(what) is unavailable. " *
               bundle_note
    end
    console!(sh, "loading $(pkg) for $(what) (one-time)…")
    try
        # Into `Main`, whose load path is the ACTIVE PROJECT — `bin`, which carries these as
        # direct deps for exactly this. An extension module resolves only the parent's deps
        # plus its own triggers, so a `using` in here would not find a weakdep of ROTIR that
        # is not one of ROTIRGUIExt's triggers. The binding is not what is wanted anyway:
        # loading the package is what activates the extension and gives the stubs their
        # methods.
        Core.eval(Main, Expr(:using, Expr(:., pkg)))
        # The documented way to force extension activation. Not needed in practice — loading
        # the trigger activates the extension synchronously — but harmless, and the check
        # below is what actually reports a failure.
        Base.retry_load_extensions()
    catch e
        failed[] = true
        return "$(what) needs $(pkg), which failed to load: " * sprint(showerror, e)
    end
    if !available()
        failed[] = true
        return "$(pkg) loaded but $(extname) did not"
    end
    return ""
end

"""
    shell_set_polyft_backend(kind) -> String

Choose the forward kernel: `"auto"` to let each code path pick (see `ShellState.kernel_auto`),
or one of the three to force it everywhere. Invalidates the χ² cache, since the number it holds
came from the kernel that was selected when it was computed — and the three do not agree to
better than the quadrature's 6.8e-7 at HEALPix 3.
"""
function shell_set_polyft_backend(kind)
    sh = _sh()
    k = Symbol(_qmlstr(kind))
    k in (:auto, :t3, :nufft, :turbo, :scalar) ||
        return "backend must be auto, t3, nufft, turbo or scalar"
    # AUTO IS NOT A KERNEL, so it does not go into `POLYFT_BACKEND`, which the library reads
    # and the test suite assigns to. It sets the process default back to the one the forward
    # paths want and lets each fit scope its own — and it deliberately does NOT load
    # LoopVectorization here: that cost belongs to the first fit that needs it, not to a click.
    if k === :auto
        sh.kernel_auto[] = true
        ROTIR.POLYFT_BACKEND[] = :t3
        sh.chi2key[] = nothing
        console!(sh, "polyft backend: auto (t3 here, turbo for gradient fits and samplers)")
        refresh_both!(sh)
        return "polyft backend: auto"
    end
    if k === :turbo
        err = _ensure_turbo!(sh)
        isempty(err) || return err
    end
    if k === :nufft
        # FINUFFT is a weak dependency now, so `:nufft` is loaded on demand exactly as
        # `:turbo` is. `:t3` — the default — needs nothing.
        err = _ensure_finufft!(sh)
        isempty(err) || return err
    end
    sh.kernel_auto[] = false
    ROTIR.POLYFT_BACKEND[] = k
    sh.chi2key[] = nothing
    console!(sh, "polyft backend: $(k)")
    refresh_both!(sh)
    return "polyft backend: $(k)"
end

"""
    shell_set_zoom_step(x) -> String

How far one wheel detent zooms, as a factor. Zero restores Makie's own 1.11.

Applied to the three sky axes, which are the ones that zoom: the Mollweide is a whole-surface
projection with nothing to zoom into, and the χ² panel is a bar chart.
"""
function shell_set_zoom_step(x)
    sh = _sh()
    v = set_zoom_step!(Float64(x))
    # The 2-D handlers read the step per event, but `Axis3` keeps Makie's `ScrollZoom`, whose
    # speed is fixed when it is constructed — so those two are re-registered.
    for c in (sh.star, sh.imstar)
        c === nothing || _install_zoom!(c.scene)
    end
    return v > 0 ? Printf.@sprintf("wheel zoom: %.2fx per detent", v) :
                   "wheel zoom: default"
end

"""
    shell_version() -> String

Which code is actually running.

The first thing to establish about any bug report, and the settings panel is where a user
already comes to look at what the WINDOW is doing rather than at their data. `pkgversion` of
the parent package, not of this extension — an extension reports the version of whatever
built it, which for a dev checkout is not a number anybody can act on.
"""
function shell_version()
    v = try
        string(pkgversion(ROTIR))
    catch
        "unknown"
    end
    return "ROTIR v" * v * "  ·  Julia " * string(VERSION)
end

"The decorations the sky views can draw, in the order the ticks appear."
# `:spin` is ONE tick driving both the axis and the arrow: they are two halves of the same
# annotation — the line says where the pole is, the arrow says which way it turns — and either
# alone reads as a drawing mistake.
# `:graticules` is LAST because it is the only one with settings of its own: ticking it opens
# a spacing and a colour beside it, and a control that grows the row should grow it at the end
# rather than pushing everything after it sideways.
const DECORATIONS = [(:limb, "limb"), (:compass, "N/E"), (:spin, "spin axis + arrow"),
                     (:plotmesh, "mesh"), (:graticules, "graticules")]

"""
    shell_decorations() -> String

`name\tlabel\ton` per decoration, for the tick boxes.
"""
shell_decorations() = join(("$(k)\t$(l)\t$(_sh().decor[k] ? "1" : "0")"
                            for (k, l) in DECORATIONS), "\n")

"""
    shell_graticule() -> String

`lat_step_deg\tlon_step_deg\tcolour` — the graticule settings behind the tick.
"""
shell_graticule() = (sh = _sh();
    Printf.@sprintf("%.4g\t%.4g\t%s", sh.gratlat[], sh.gratlon[], sh.gratcolor[]))

"The colours the graticule dropdown offers, dark and light so it reads on either end of a map."
const GRATICULE_COLORS = ["black", "white", "grey", "red", "cyan", "yellow"]

shell_graticule_colors() = join(GRATICULE_COLORS, "\n")

"""
    shell_set_graticule(lat_deg, lon_deg, colour) -> String

Spacing (in degrees) and colour of the graticule.

Steps are clamped to 5–90 degrees: below 5 the lines merge into a grey wash at any resolution
this GUI draws, and above 90 there is at most one parallel to see.
"""
function shell_set_graticule(lat_deg, lon_deg, colour)
    sh = _sh()
    la = _qmlreal(lat_deg, sh.gratlat[])
    lo = _qmlreal(lon_deg, sh.gratlon[])
    sh.gratlat[] = clamp(la, 5.0, 90.0)
    sh.gratlon[] = clamp(lo, 5.0, 90.0)
    c = _qmlstr(colour)
    c in GRATICULE_COLORS && (sh.gratcolor[] = c)
    refresh_both!(sh)
    refresh_image_tab!(sh)
    return ""
end

"""
    shell_set_decoration(name, on) -> String

Turn one decoration on or off, everywhere.
"""
function shell_set_decoration(name, on)
    sh = _sh()
    k = _qmlname(name)
    k === nothing && return NO_QML_NAME
    haskey(sh.decor, k) || return "unknown decoration $(name)"
    sh.decor[k] = _qmlstr(on) == "1"
    refresh_both!(sh)
    refresh_image_tab!(sh)
    return ""
end

"""
    shell_surface_field() -> String

`intensity\tmodel\tband_um` — what the surface views are colouring by.
"""
function shell_surface_field()
    sh = _sh()
    # The wavelength REPORTED is the one that would be used, not the sentinel: showing "0" in
    # a field labelled λ (µm) is showing the caller's shorthand back to them instead of the
    # number the Planck conversion is about to use.
    b = sh.band[] > 0 ? sh.band[] : _data_band(sh)
    return Printf.@sprintf("%s\t%s\t%.4f", sh.intensity[] ? "1" : "0", sh.intensity_model[],
                           b * 1e6)
end

"""
    shell_set_surface_field(intensity, model, band_um) -> String

Colour the surface views by TEMPERATURE or by emergent INTENSITY.

Intensity is what the data actually constrain: the limb-darkening map multiplied in, and — with
`model = "planck"` — the temperature converted to a surface brightness at the observing
wavelength, which is strongly non-linear across a gravity-darkened star. Temperature is what
the model IS. Neither is the right default for both questions, hence the tick.

`band_um` is in micrometres for the form's sake; 0 means "take the mean wavelength from the
current dataset", which is almost always what is wanted and is one fewer number to look up.
"""
function shell_set_surface_field(intensity, model, band_um)
    sh = _sh()
    newint = _qmlstr(intensity) == "1"
    m = Symbol(_qmlstr(model))
    m in (:linear, :planck) || return "unknown intensity model $(model)"
    # A pinned range is in the units of the old quantity, and would be meaningless on the new.
    (newint != sh.intensity[] || m !== sh.intensity_model[]) && (sh.colorrange[] = nothing)
    sh.intensity[] = newint
    sh.intensity_model[] = m
    b = _qmlreal(band_um, 0.0)
    sh.band[] = b > 0 ? b * 1e-6 : _data_band(sh)
    # Fall back to the data's own mean wavelength whenever the field is empty, zero or
    # unparseable — the band the observations were taken in is the only defensible default for
    # a conversion whose whole point is that it depends on the band.
    if m === :planck && sh.band[] <= 0
        return "Planck needs a wavelength: load a dataset or type one in µm"
    end
    refresh_data_tab!(sh)
    refresh_model_tab!(sh)
    refresh_image_tab!(sh)
    return sh.intensity[] ? "colouring by intensity ($(m))" : "colouring by temperature"
end

# The dataset's own mean wavelength, in metres. Better than any default: the Planck conversion
# is being asked for precisely because the answer depends on the band the data were taken in.
function _data_band(sh::ShellState)
    d = current_dataset(sh.session)
    d === nothing && return 0.0
    try
        return Float64(Statistics.mean(d.data[1].uv_lam))
    catch
        return 0.0
    end
end

"""
    surface_values(sh, tmap, star; visible_only) -> Vector{Float64}

The quantity the surface views draw, honouring the intensity tick.

`visible_only` picks between the two indexing conventions the drawing layer uses: the 2-D
views take one value per VISIBLE tessel, the 3-D and Mollweide views one per pixel. Getting
that wrong scrambles the map with no error, which is why it is one function rather than two
call sites.
"""
function surface_values(sh::ShellState, tmap, star; visible_only::Bool)
    sh.intensity[] || return Float64.(visible_only ? tmap[star.index_quads_visible] : tmap)
    b = sh.band[] > 0 ? sh.band[] : nothing
    Imap = sh.intensity_model[] === :linear ? tmap :
           b === nothing ? tmap : ROTIR.intensity(Float64.(tmap), :planck, b)
    visible_only || return Float64.(Imap)
    # LIMB DARKENING ONLY ON THE VIEW THAT HAS A LIMB.
    #
    # `ldmap` is a viewing-geometry quantity: it is zero wherever the surface faces away, on
    # exactly half the tessels (measured: 384 non-zero of 768). Multiplying a WHOLE-SURFACE
    # map by it therefore blanked the far side — the Mollweide and the 3-D view came out with
    # half the star black, which read as a map that had never been populated.
    #
    # The orthographic view is the one where limb darkening belongs, and it is the one the
    # intensity tick was added for: it shows what the interferometer sees, which is the
    # brightness of the visible hemisphere, darkened towards the limb. A Mollweide is a map OF
    # THE SURFACE, where "how foreshortened is this patch from here" is not a property of the
    # surface at all — and the 3-D view can be rotated, so a per-tessel factor baked for one
    # viewpoint would be wrong from every other.
    return Float64.(Imap[star.index_quads_visible]) .*
           Float64.(star.ldmap[star.index_quads_visible])
end

"""
    visible_values(sh, allv, star) -> Vector{Float64}

The orthographic view's values, from the whole-surface ones already computed.

The two differ by exactly the limb-darkening factor, so the expensive half — the `intensity`
evaluation, a Planck integral per tessel when that law is chosen, on every keystroke in the
parameter form — is done once by [`surface_values`](@ref) and this applies the cheap half.

It WAS a plain slice, `allv[star.index_quads_visible]`, and that was correct only while limb
darkening went onto both conventions. It does not any more (see `surface_values`), so a slice
silently drew the orthographic view without a limb — flat across the disk, exactly as the
temperature view — and the intensity tick stopped doing anything on the two tabs that took
this route.
"""
visible_values(sh::ShellState, allv::AbstractVector, star) =
    sh.intensity[] ?
    Float64.(allv[star.index_quads_visible]) .*
        Float64.(star.ldmap[star.index_quads_visible]) :
    Float64.(allv[star.index_quads_visible])

"""
    whole_surface_label(sh) -> String

The colour-bar label for a view that has no limb — the 3-D scene and the Mollweide.

NOT the same rule as the orthographic view's, and the difference is the honest statement of
what the intensity tick does. The tick means two things at once: switch the quantity from a
temperature to an emergent intensity, and multiply limb darkening in. Only the first survives
on a whole-surface view, and under the `linear` law the first is the identity — `surface_values`
returns `tmap` untouched — so a "linear intensity" on these views IS the temperature map, and
labelling it `I (arb.)` would claim a quantity that is not there.

Under `planck` it is a real surface brightness at the observing wavelength, which is a
different quantity from T and is labelled as one. These two colour bars read `T (K)` as a
CONSTANT before this existed, so a Planck map was drawn under a temperature label.
"""
whole_surface_label(sh::ShellState) =
    (sh.intensity[] && sh.intensity_model[] === :planck) ? "I (arb.)" : "T (K)"


"""
    build_epoch_star(sh) -> (star, values) or nothing

The current model evaluated at the current epoch's time, ready to draw.

Returns `nothing` rather than raising when there is no model or the parameters do not build:
the preview is redrawn on every keystroke in the form, and a transiently invalid parameter set
must leave the last good picture on screen rather than take the window down.
"""
function build_epoch_star(sh::ShellState)
    m = current_model(sh.session)
    m === nothing && return nothing
    d = current_dataset(sh.session)
    t = d === nothing ? 0.0 : d.tepochs[clamp(sh.session.current_epoch, 1, length(d.tepochs))]
    p = star_params(m)
    isempty(validate_star_params(p)) || return nothing
    try
        tess = model_tessellation(sh)
        star = create_star(tess, p, t; secondary = m.secondary)
        tmap = parametric_temperature_map(p, star; secondary = m.secondary)
        c = m.companion
        c === nothing && return (star, tmap)
        # A BINARY: the companion at the same epoch, and where the orbit puts it. The offset
        # is in the sky frame the polygons are drawn in, which is what `orbit_to_rotir_offset`
        # returns — and it needs JD, while the epoch times here are relative days.
        p2 = star_params(c)
        isempty(validate_star_params(p2)) || return (star, tmap)
        star2 = create_star(tess, p2, t; secondary = true)
        tmap2 = parametric_temperature_map(p2, star2; secondary = true)
        off = try
            if c.place === :offset
                (c.offset[1], c.offset[2])
            elseif d === nothing
                (0.0, 0.0)
            else
                i = clamp(sh.session.current_epoch, 1, length(d.mjd))
                only(companion_offsets(sh, m, d.mjd[i:i]))
            end
        catch
            (0.0, 0.0)
        end
        return (star, tmap, star2, tmap2, off)
    catch err
        console!(sh, "preview failed: $(sprint(showerror, err))")
        return nothing
    end
end

"""
    check_mesh!(sh, mesh, colors, what) -> Bool

Report a mesh whose colour vector does not match it, before the renderer chokes on it.

QMLMakie's render callback is `catch; Core.println("exception in render")` — no binding, so
the exception is discarded and a broken plot says only that something went wrong, on a thread
whose stdout is the launching terminal. Every render failure found in this GUI has had to be
reconstructed by elimination for that reason.

This checks the one invariant that is cheap to test and known to kill the render: a mesh needs
exactly one colour per vertex. A zero-length colour vector fails to build a render object at
all ("Failed to resolve gl_renderobject"), and a short one is worse — it fails at draw time,
after the plot exists. `update_state_before_display!` does NOT catch either: it was tried here
and returned cleanly on a deliberately mismatched pair, which is why this checks the invariant
directly rather than running the renderer's preparation step.
"""
function check_mesh!(sh::ShellState, mesh, colors, what::AbstractString)
    mesh isa Makie.GeometryBasics.Mesh || return true
    n = length(Makie.GeometryBasics.coordinates(mesh))
    m = length(colors)
    n == m && return true
    console!(sh, "$(what) will not render: $(n) vertices but $(m) colours")
    return false
end

# The decoration flags as keywords for `show_map!`. `:spin` expands into the two the drawing
# layer actually takes.
"""
    refresh_model_tab!(sh)

Redraw the 3-D preview and the Mollweide from the current model. GUI thread only.
"""
function _decor(sh::ShellState)
    d = sh.decor
    # DEGREES STRAIGHT THROUGH. This used to convert to line COUNTS
    # (`nlat = round(180/gratlat)`), and the conversion could not express what was asked:
    # `graticule_segments` spaces `nlat` parallels at `180/(nlat + 1)`, so a requested 75
    # degrees became 2 parallels 60 degrees apart. It now takes the spacing itself and places
    # parallels at multiples of it from the equator — 75 gives -75, 0, +75.
    return (limb = d[:limb], compass = d[:compass], graticules = d[:graticules],
            plotmesh = d[:plotmesh],
            rotation_axis = d[:spin], rotation_arrow = d[:spin],
            graticule_dlat = sh.gratlat[], graticule_dlon = sh.gratlon[],
            graticule_color = Symbol(sh.gratcolor[]),
            colorrange = sh.colorrange[])
end

function refresh_model_tab!(sh::ShellState; got = build_epoch_star(sh))
    m = current_model(sh.session)
    if got === nothing
        msg = current_model(sh.session) === nothing ?
              "no model — choose a surface type and press + model" :
              "the current parameters do not build; see the message under the form"
        idle!(sh.msky, msg); idle!(sh.star, msg); idle!(sh.moll, msg)
        return sh
    end
    star, tmap = got[1], got[2]
    # A binary carries two more entries and the sky offset between them.
    bin2 = length(got) >= 5 ? (got[3], got[4], got[5]) : nothing
    # All three, every time. They are three views of one map, and refreshing only the visible
    # one would mean a stale picture appearing the moment the view is switched — which is
    # worse than the cost, because switching is meant to be instant.
    lbl = sh.intensity[] ? "I (arb.)" : "T (K)"
    sh.msky.cbarlabel[] = lbl
    # The 3-D scene and the Mollweide have no limb, so they label by their own rule.
    sh.star.cbarlabel[] = whole_surface_label(sh)
    sh.moll.cbarlabel[] = whole_surface_label(sh)
    # ONE evaluation, three views. This called `surface_values` three times — twice with
    # exactly the same arguments — and it runs on every keystroke in the parameter form. With
    # the intensity tick on and Planck selected that is a Planck evaluation over every tessel,
    # done twice for nothing. `visible_values` takes the orthographic view's values off the
    # full ones instead, which is the cheap half of the difference between them.
    allv = surface_values(sh, tmap, star; visible_only = false)
    visv = visible_values(sh, allv, star)
    if bin2 === nothing
        show_map!(sh.msky, star, visv; star_params = star_params(m), _decor(sh)...)
    else
        # Both components, one colour range. The 3-D and Mollweide views stay on the PRIMARY:
        # a Mollweide is a map of one surface, and two of them in one frame is two pictures.
        star2, tmap2, off = bin2
        v2 = surface_values(sh, tmap2, star2; visible_only = true)
        # The COMPANION's own parameters too: every decoration is drawn for both components
        # now, and `_polar_radius` reads `rpole`/`radius` — the primary's would draw the
        # secondary's spin axis at the primary's length.
        show_binary_map!(sh.msky, star, visv, star2, v2, off;
                         star_params = star_params(m),
                         star_params2 = star_params(m.companion), _decor(sh)...)
    end
    if bin2 === nothing
        show_star3d!(sh.star, star, allv; colorrange = sh.colorrange[])
    else
        # BOTH components in 3-D too, at the same separation the sky view draws — the z term
        # included, which the sky view drops. That is the whole reason to look at the pair in
        # three dimensions: an eclipse reads as one star in front of the other only once the
        # line-of-sight offset is there and the camera is tilted off the observer's axis.
        star2, tmap2, off = bin2
        allv2 = surface_values(sh, tmap2, star2; visible_only = false)
        c = m.companion
        zoff = c === nothing ? 0.0 : (c.place === :offset ? c.offset[3] : 0.0)
        # The relative orbit track only when an orbit is what places them. Drawing a closed
        # ellipse through a pair positioned by hand would be a claim the data do not make.
        track = if c !== nothing && c.place === :orbit
            try
                relative_orbit_track(orbit_bparams(sh.orbit; binary = m))
            catch
                Makie.Point3f[]
            end
        else
            Makie.Point3f[]
        end
        show_binary3d!(sh.star, star, allv, star2, allv2,
                       (Float64(off[1]), Float64(off[2]), zoff); track = track,
                       colorrange = sh.colorrange[])
    end
    show_mollweide!(sh.moll, allv, star; colorrange = sh.colorrange[])
    # Say which panel is wrong and how, rather than leaving "exception in render" as the only
    # evidence. Cheap: two length comparisons against work already done.
    check_mesh!(sh, sh.moll.mesh[], sh.moll.colors[], "the Mollweide")
    check_mesh!(sh, sh.star.mesh[], sh.star.colors[], "the 3-D view")
    return sh
end

"""
    refresh_data_tab!(sh)

Redraw the Explore tab: its observable plot, for the current epoch. `got` is accepted and
ignored, so the callers that share one geometry build between the tabs need not know which
tabs use it.
"""
refresh_data_tab!(sh::ShellState; got = nothing) = (refresh_obs!(sh); sh)

"""
    companion_offsets(sh, m, mjd) -> Vector{Tuple{Float64,Float64}}

Where the secondary sits, in `(West, North)` mas, at each of `mjd`.

The two placements answer the same question with different amounts of model. `:orbit` asks the
Orbit tab, which gives a different separation at every epoch and needs elements the data have
to constrain. `:offset` is one displacement, the same throughout — the honest model for a
snapshot, where there is no orbital motion to see and seven elements would be fitted to a
separation and a position angle.

Returned per epoch either way, so nothing downstream has to know which is in force.
"""
function companion_offsets(sh::ShellState, m, mjd::AbstractVector)
    c = m.companion
    c === nothing && return fill((0.0, 0.0), length(mjd))
    if c.place === :offset
        return fill((c.offset[1], c.offset[2]), length(mjd))
    end
    bp = orbit_bparams(sh.orbit; binary = m)
    # The elements are in JD and datasets carry MJD; getting that wrong puts the companion in
    # the wrong place with no error anywhere.
    return [orbit_to_rotir_offset(bp, t + 2_400_000.5) for t in mjd]
end

"""
    _binary_model_state(sh, d, m, c, p, p2, key) -> NamedTuple or nothing

The binary half of [`model_state`](@ref): both components, their separation at each epoch, and
the χ² of the PAIR against the data.

Built the way `demos/spica_binary_roche.jl` builds it, and for the same reason it does it that
way: each component is created on its own with `create_star_multiepochs`, its own
`secondary` convention, and the two are brought together only at the visibility stage by
`binary_phase_shift`. `create_binary_geometry` exists and sets `center_offsets` instead, but
the offset the transform needs is per-EPOCH and in the sky frame the polygon FT uses, which is
exactly what `orbit_to_rotir_offset` returns.

The orbit comes from the Orbit tab — one description of the system, shared, rather than a
second copy that can disagree with it.
"""
function _binary_model_state(sh::ShellState, d, m, c, p, p2, key)
    try
        tess = model_tessellation(sh)
        offs = companion_offsets(sh, m, d.mjd)
        stars1 = create_star_multiepochs(tess, p,  d.tepochs; secondary = false)
        stars2 = create_star_multiepochs(tess, p2, d.tepochs; secondary = true)
        x1 = parametric_temperature_map(p,  stars1[1]; secondary = false)
        x2 = parametric_temperature_map(p2, stars2[1]; secondary = true)
        chi2 = map(eachindex(d.data)) do i
            ox, oy = offs[i]
            ph = binary_phase_shift(d.data[i].uv, ox, oy)
            v2m, t3am, t3pm = binary_observables(x1, stars1[i], x2, stars2[i], d.data[i], ph)
            dat = d.data[i]
            cv2 = sum(abs2, (v2m   .- dat.v2)    ./ dat.v2_err)
            ca  = sum(abs2, (t3am  .- dat.t3amp) ./ dat.t3amp_err)
            cp  = sum(abs2, ROTIR.mod360(t3pm .- dat.t3phi) ./ dat.t3phi_err)
            n   = dat.nv2 + dat.nt3amp + dat.nt3phi
            (v2 = cv2, t3amp = ca, t3phi = cp, total = cv2 + ca + cp,
             nv2 = dat.nv2, nt3amp = dat.nt3amp, nt3phi = dat.nt3phi, ndata = n,
             v2r = cv2 / max(dat.nv2, 1), t3ampr = ca / max(dat.nt3amp, 1),
             t3phir = cp / max(dat.nt3phi, 1), totalr = (cv2 + ca + cp) / max(n, 1))
        end
        out = (stars = stars1, x = x1, chi2 = chi2,
               stars2 = stars2, x2 = x2, offsets = offs)
        sh.chi2key[] = key; sh.chi2cache[] = out
        return out
    catch err
        console!(sh, "binary χ² failed: $(sprint(showerror, err))")
        sh.chi2key[] = nothing
        return nothing
    end
end

"""
    model_state(sh) -> NamedTuple or nothing

The current model built against the current dataset, CACHED: the geometry at every epoch with
`setup_oi!` already run, the map those parameters generate, and the per-epoch
[`chi2_breakdown`](@ref).

`nothing` — not a vector of zeros — when either the dataset or the model is missing, so a panel
with nothing to compare against can say so rather than showing a fit that was never computed.

ONE cache for all three, because they come from one computation. `create_star_multiepochs`
plus `setup_oi!` is the polygon Fourier transform against every uv point, and both the χ²
column and the model prediction drawn over the data need exactly it. Computing it twice was
what made switching from V² to closure phase cost 130 ms of arithmetic on six λ And epochs
when the answer had not changed.

The key is the dataset identity, its epoch count, the resolved parameters AND the mesh —
level and element type. The mesh was missing from it, so changing the HEALPix level left the
χ² column showing numbers from the previous mesh until something else invalidated the entry.
"""
function model_state(sh::ShellState)
    d = current_dataset(sh.session)
    m = current_model(sh.session)
    (d === nothing || m === nothing) && return nothing
    # ONE orbit for the pair, resolved before anything is built from it: the Orbit tab's when
    # it is what places the secondary, the primary's otherwise, and the companion's copy
    # derived from it either way.
    _sync_orbit_from_tab!(sh, m)
    _sync_shared_orbit!(m)
    p = star_params(m)
    isempty(validate_star_params(p)) || return nothing
    c = m.companion
    p2 = c === nothing ? nothing : star_params(c)
    if c !== nothing && !isempty(validate_star_params(p2))
        return nothing
    end
    # The ORBIT is part of the key for a binary: the separation at each epoch comes from the
    # Orbit tab, so editing an element there changes this χ² and must invalidate the cache.
    key = (objectid(d), length(d.data), m.secondary, p, p2,
           c === nothing ? nothing : c.place,
           c === nothing ? nothing : c.offset,
           c === nothing || c.place !== :orbit ? nothing : apply_orbit_ties(sh.orbit),
           c === nothing || c.place !== :orbit ? nothing : sh.orbit.q,
           sh.nside_exp[], sh.precision[], sh.tessel[])
    sh.chi2key[] === key && return sh.chi2cache[]
    c === nothing || return _binary_model_state(sh, d, m, c, p, p2, key)
    try
        tess  = model_tessellation(sh)
        stars = create_star_multiepochs(tess, p, d.tepochs; secondary = m.secondary)
        # NO `setup_oi!`. This is the panel's χ² table and the model prediction drawn over the
        # data — one evaluation per parameter change, on a geometry that changed with it. The
        # dense matrix was built for it and read once: MEASURED at 1423 ms for one refresh on
        # polaris, which is the delay after "+ model" that started this. `observables` sees an
        # empty `polyft` and takes the matrix-free route.
        x = parametric_temperature_map(p, stars[1]; secondary = m.secondary)
        out = (stars = stars, x = x,
               chi2 = [chi2_breakdown(x, stars[i], d.data[i]) for i in eachindex(d.data)])
        sh.chi2key[] = key; sh.chi2cache[] = out
        return out
    catch err
        console!(sh, "χ² failed: $(sprint(showerror, err))")
        sh.chi2key[] = nothing
        return nothing
    end
end

"""
    epoch_chi2(sh) -> Vector{NamedTuple} or nothing

Per-epoch [`chi2_breakdown`](@ref) of the current model against the current dataset.

The χ² half of [`model_state`](@ref), which is where the caching lives. Without it this ran on
every tab switch, and it is the single most expensive thing the GUI does.
"""
epoch_chi2(sh::ShellState) = (s = model_state(sh); s === nothing ? nothing : s.chi2)

# ── the observable plot ──────────────────────────────────────────────────────────────────
#
# V² against baseline, closure phase, uv coverage: exactly the plots the OITOOLS GUI draws, and
# literally its plots — `build_canvas` and `update_canvas!` come from OITOOLSGUIExt, reached at
# RUNTIME through `Base.get_extension`. Re-implementing them here would mean two versions of
# the same figure that could disagree about which array feeds an observable or what colour a
# baseline is, which is the failure the shared `PlotData` layer exists to prevent.
#
# Extensions cannot import each other, but this one and OITOOLSGUIExt are triggered by the same
# three packages, so by the time any of this runs both exist.
#
# What ROTIR adds is the second series: the current model's prediction, over the data.

"The OITOOLS GUI extension, or `nothing` if it is not loaded."
_oitools_gui() = Base.get_extension(OITOOLS, :OITOOLSGUIExt)

"""
    build_obs_canvas(fig) -> canvas

An OITOOLS `LiveCanvas`, or `nothing` without OITOOLSGUIExt. The model's prediction and the
residuals are drawn by the canvas itself (`update_canvas!`'s `overlay` and `residual`), the
same way OITOOLS' Explore draws them.
"""
function build_obs_canvas(fig)
    O = _oitools_gui()
    O === nothing && return nothing
    ax = Makie.Axis(fig[1, 1])
    c  = O.build_canvas(fig, ax)
    # The canvas carries a hidden heatmap for OITOOLS' reconstructed image, on cell edges
    # −0.5…1.5. A figure first displayed with a log y axis — "Save plot…" of a log V² — runs
    # log10 over those edges and fails to render at all. ROTIR never shows an image here, so
    # the layer is moved to positive coordinates where a log axis can transform it.
    c.imagey[] = Float32[0.5, 1.0, 1.5]
    # Right-click returns the view to the data. `build_canvas` gives this axis Makie's own
    # `:limitreset`, which resets to `ax.limits[]` — and `update_canvas!` calls `autolimits!`,
    # leaving that automatic, so the reset had nothing to go back to and did nothing visible.
    # `autolimits!` is the honest "everything" here, since unlike a sky view there is no framed
    # extent to return to: the plot is whatever the data span.
    haskey(Makie.interactions(ax), :limitreset) &&
        Makie.deregister_interaction!(ax, :limitreset)
    Makie.register_interaction!(ax, :limitreset) do event::Makie.MouseEvent, axis
        if event.type === Makie.MouseEventTypes.rightclick
            Makie.autolimits!(axis)
            return Makie.Consume(true)
        end
        return Makie.Consume(false)
    end
    return c
end

# Which PLOTS each availability flag enables.
#
# The two vocabularies are not the same and mapping one onto the other by name was wrong in
# both directions. `observable_availability` answers "does the file have this table?" and
# returns `v2, t3amp, t3phi, cvis, flux, diffvis`; `OBS_SPECS` — what `update_canvas!` can
# actually draw — is keyed `v2, t3amp, t3amp_max, t3phi, t3phi_max, visamp, visphi, flux,
# diffphi, diffvisamp`.
#
# The `_max` variants are deliberately NOT listed. They are the same observable against a
# different x, and OITOOLS does not offer them as views either — it uses them for residual
# plots. A view per combination doubles the button row to say one thing, so the choice of
# baseline is a TICK beside the plot and `shell_set_obs_view` appends the suffix.
#
# So `cvis` named a view that does not exist (it is the flag for "there are visamp/visphi
# data", and `kind = :cvis` raises "unknown plot"), which is the empty panel polaris offered;
# and `visamp`, `visphi`, `diffphi` and `diffvisamp` could never be reached at all, because no
# flag is called that. The `_max` variants plot a closure quantity against the LONGEST baseline
# of its triangle rather than the geometric mean — a different question about the same data,
# and the one to ask when resolution is what matters.
const OBS_KINDS_FOR_FLAG = Dict(
    "v2"      => ["v2"],
    "t3amp"   => ["t3amp"],
    "t3phi"   => ["t3phi"],
    "cvis"    => ["visamp", "visphi"],
    "flux"    => ["flux"],
    "diffvis" => ["diffphi", "diffvisamp"],
)

"""
    shell_obs_kinds() -> String

Which observable views the current dataset can show, as `key\tavailable` per line.

Delegated to OITOOLS' own availability check, so ROTIR cannot offer a panel for a table the
file does not have.
"""
function shell_obs_kinds()
    sh = _sh()
    d = current_dataset(sh.session)
    d === nothing && return ""
    always = ["uv"]                    # geometry: present whenever there is data at all
    O = _oitools_gui()
    O === nothing && return join(("$(k)\t1" for k in vcat(always, ["v2","t3amp","t3phi"])), "\n")
    flags = try
        # OITOOLS returns "v2=1,t3amp=1,t3phi=1,cvis=0,..." — comma-separated `key=value`, NOT
        # the tab-separated rows the rest of this bridge uses. Parsing it here rather than in
        # QML keeps one format crossing the boundary; assuming its shape was why the Explore
        # tab silently fell back to the sky map, with no error anywhere.
        O.observable_flags_string(reshape(d.data, 1, length(d.data)))
    catch
        ""
    end
    rows = ["$(k)\t1" for k in always]
    for part in split(flags, ',')
        kv = split(strip(part), '=')
        length(kv) == 2 || continue
        k = String(strip(kv[1]))
        avail = strip(kv[2])
        for kind in get(OBS_KINDS_FOR_FLAG, k, String[])
            push!(rows, "$(kind)\t$(avail)")
        end
    end
    # A file whose tables confuse the check should still be explorable.
    length(rows) == length(always) &&
        append!(rows, ("$(k)\t1" for k in ("v2", "t3amp", "t3phi")))
    return join(rows, "\n")
end

"""
    shell_set_obs_view(kind, color, compare, logy, maxbaseline = "0", panels = "0") -> String

Redraw the observable plot, with the options OITOOLS' Explore offers: which quantity, what to
colour by, what to compare the data with, a logarithmic y axis, closure quantities against the
longest baseline, and one panel per baseline (or triplet, or station).

`compare` is `"data"`, `"model"`, `"image"`, `"residmodel"` or `"residimage"` — the data alone,
the model's or the current reconstruction's observables over them, or `(model − data)/σ`.
`"1"` and `"0"` are read as `"model"` and `"data"`.
"""
function shell_set_obs_view(kind, color, compare, logy, maxbaseline = "0", panels = "0")
    sh = _sh()
    k = _qmlstr(kind)
    # A closure quantity can be plotted against the geometric mean of its triangle's baselines
    # or against the LONGEST leg — the resolution the closure actually probes. Same data, so it
    # is a modifier on the view rather than a view of its own.
    if _qmlstr(maxbaseline) == "1" && k in ("t3amp", "t3phi")
        k *= "_max"
    end
    c = _qmlstr(compare)
    c = c == "1" ? "model" : c == "0" ? "" : c
    sh.obskind[]    = Symbol(k)
    sh.obscolor[]   = Symbol(_qmlstr(color))
    sh.obscompare[] = c in ("model", "image", "residmodel", "residimage") ? Symbol(c) : :data
    sh.obslog[]     = _qmlstr(logy) == "1"
    sh.obspanels[]  = _qmlstr(panels) == "1"
    return refresh_obs!(sh)
end

"""
    shell_obs_grouping_noun(kind) -> String

What one panel of the per-group view holds for `kind` — `"baseline"`, `"triplet"` or
`"station"` — or `""` when the view has no groups (uv coverage). From OITOOLS, so the tick is
labelled the way OITOOLS labels it.
"""
function shell_obs_grouping_noun(kind)
    k = Symbol(_qmlstr(kind))
    O = _oitools_gui()
    (O === nothing || !haskey(O.OBS_SPECS, k)) && return ""
    return O.grouping_noun(k)
end

"""
    shell_obs_compare_available() -> String

`"model=1,image=0"`: whether the current model, and the current reconstruction, can be compared
with the data — what the Compare choices are enabled by.
"""
function shell_obs_compare_available()
    sh = _sh()
    d = current_dataset(sh.session)
    m = current_model(sh.session)
    okm = d !== nothing && m !== nothing && model_state(sh) !== nothing
    oki = d !== nothing && m !== nothing && m.companion === nothing &&
          current_image(sh.session) !== nothing
    return "model=$(okm ? 1 : 0),image=$(oki ? 1 : 0)"
end

"The kinds a log y axis applies to: the positive ones, as in OITOOLS."
const OBS_LOG_KINDS = (:v2, :t3amp, :t3amp_max, :visamp, :flux, :diffvisamp)

"The kinds ROTIR predicts, and so can overplot or turn into residuals."
const OBS_PREDICTED = (:v2, :t3amp, :t3amp_max, :t3phi, :t3phi_max)

"""
    _obs_logscale(sh) -> Bool

Whether the observable plot's y axis is logarithmic for the view it is currently showing.

The tick is remembered across views, but it applies only to positive quantities: a closure
phase takes both signs, uv coverage is a geometry, and a residual is centred on zero. Decided
here rather than only in QML, because the callback is also reachable from a script and from the
tests, and a log axis on a phase is not a rescaling — it is half the points thrown away.
"""
_obs_logscale(sh::ShellState) = sh.obslog[] && sh.obskind[] in OBS_LOG_KINDS &&
                                !(sh.obscompare[] in (:residmodel, :residimage))

"""
    _predicted_obs(sh, dat, epoch, source) -> (; v2, t3amp, t3phi) or an error String

The observables the current model (`source = :model`) or the current reconstruction
(`:image`) predicts for one epoch, in the data's own point order.

The model comes from [`model_state`](@ref), so it is the same evaluation the χ² column is
computed from, both components of a binary included. The reconstruction is evaluated on the
current model's geometry at the map's own HEALPix level, with the intensity law the Imaging
tab reconstructs with — the same three things the reconstruction itself used.
"""
function _predicted_obs(sh::ShellState, dat, epoch::Int, source::Symbol)
    m = current_model(sh.session)
    m === nothing && return "no model to compare with"
    if source === :model
        st = model_state(sh)
        st === nothing && return "the model's parameters do not build a star"
        v2, t3a, t3p = if hasproperty(st, :stars2)
            ox, oy = st.offsets[epoch]
            binary_observables(st.x, st.stars[epoch], st.x2, st.stars2[epoch], dat,
                               binary_phase_shift(dat.uv, ox, oy))
        else
            observables(st.x, st.stars[epoch], dat)
        end
        return (; v2, t3amp = t3a, t3phi = t3p)
    end
    im = current_image(sh.session)
    im === nothing && return "no reconstruction yet — run one in Imaging first"
    m.companion === nothing ||
        return "a reconstruction is one surface and the model is a binary"
    ii = _imaging_intensity(sh, m)
    isempty(ii.err) || return ii.err
    d = current_dataset(sh.session)
    star = create_star(tessellation_healpix(im.nside_exp), star_params(m), d.tepochs[epoch];
                       secondary = m.secondary)
    v2, t3a, t3p = observables(im.x, star, dat; intensity_model = ii.model, band = ii.band)
    return (; v2, t3amp = t3a, t3phi = t3p)
end

"""
    _obs_compare(sh, dat, epoch) -> (; overlay, residual) or an error String

What the Compare choice adds to the plot, as the two vectors OITOOLS' `update_canvas!` takes:
the prediction to draw over the data, or `(model − data)/σ` to draw instead of them. Both
`nothing` for the data alone, and for a view ROTIR does not predict (uv, visibilities, flux).
"""
function _obs_compare(sh::ShellState, dat, epoch::Int)
    none = (; overlay = nothing, residual = nothing)
    cmp = sh.obscompare[]
    kind = sh.obskind[]
    (cmp === :data || !(kind in OBS_PREDICTED)) && return none
    src = cmp in (:model, :residmodel) ? :model : :image
    # Nothing to compare WITH is not an error, only the data alone: the choice stays set, so a
    # model added later appears over the data without the box being touched again.
    current_model(sh.session) === nothing && return none
    src === :image && current_image(sh.session) === nothing && return none
    p = try
        _predicted_obs(sh, dat, epoch, src)
    catch err
        "$(src) observables failed: $(sprint(showerror, err))"
    end
    p isa String && return p
    base = kind in (:t3amp_max,) ? :t3amp : kind in (:t3phi_max,) ? :t3phi : kind
    y = Float64.(getproperty(p, base))
    if cmp in (:residmodel, :residimage)
        obs = Float64.(getfield(dat, base)); err = Float64.(getfield(dat, Symbol(base, :_err)))
        r = base === :t3phi ? ROTIR.mod360(y .- obs) ./ err : (y .- obs) ./ err
        return (; overlay = nothing, residual = r)
    end
    # On a log axis OITOOLS drops the data's non-positive and unmeasured points before drawing,
    # and an overlay that no longer matches the data point for point is not drawn at all — so
    # the same points are dropped from it here.
    if _obs_logscale(sh)
        obs = Float64.(getfield(dat, base))
        OITOOLS.hide_unmeasured!(obs, Float64.(getfield(dat, Symbol(base, :_err))))
        y = y[isfinite.(obs) .& (obs .> 0)]
    end
    return (; overlay = y, residual = nothing)
end

"""
    _draw_obs!(sh, canvas, dat, epoch) -> String

Draw one epoch into an OITOOLS canvas the way the Explore tab is set: the per-group panels, or
the single plot with the Compare choice on it. Returns a status line; throws what the canvas
throws, for the caller to report.
"""
function _draw_obs!(sh::ShellState, c, dat, epoch::Int)
    O = _oitools_gui()
    kind = sh.obskind[]
    if sh.obspanels[] && haskey(O.OBS_SPECS, kind)
        c === sh.obs && (sh.obsinfo = String[]; sh.picktext = "")
        O.show_panels!(c, true)
        n = O.update_panels!(c, dat, kind)
        return "$(kind): $(n) $(O.grouping_noun(kind))$(n == 1 ? "" : "s"), against wavelength"
    end
    O.show_panels!(c, false)
    cmp = _obs_compare(sh, dat, epoch)
    msg = ""
    if cmp isa String
        msg = cmp
        cmp = (; overlay = nothing, residual = nothing)
    end
    info = O.update_canvas!(c, dat, kind; color = sh.obscolor[], logscale = _obs_logscale(sh),
                            overlay = cmp.overlay, residual = cmp.residual)
    c === sh.obs && (sh.obsinfo = info; sh.picktext = "")
    return msg
end

"""
    refresh_obs!(sh) -> String

Push the current epoch's data into the observable plot, with whatever the Explore tab compares
it with.

ONE epoch, not the concatenation: the epochs are separate observations of a rotating star, and
a V² curve stacked across six nights is six curves drawn on top of each other. The epoch is the
one the Explore tab's table selects.
"""
function refresh_obs!(sh::ShellState)
    sh.obs === nothing && return sh.status
    d = current_dataset(sh.session)
    d === nothing && return sh.status
    i = clamp(sh.session.current_epoch, 1, length(d.data))
    msg = try
        _draw_obs!(sh, sh.obs, d.data[i], i)
    catch err
        console!(sh, "could not draw $(sh.obskind[]): $(sprint(showerror, err))")
        return sh.status
    end
    isempty(msg) || console!(sh, msg)
    # `update_canvas!` empties the filter preview; the points the Filter panel would remove
    # from this epoch and this observable go back on.
    sh.obspanels[] || _draw_filter_preview!(sh)
    return isempty(msg) ? sh.status : msg
end

# ── view controls ────────────────────────────────────────────────────────────────────────

"""
    shell_colormaps() -> String

The colormap button labels, one per line.
"""
shell_colormaps() = colormap_names()

"""
    shell_set_colormap(name) -> String

Switch every surface canvas to one colormap, and re-colour what is already drawn.

All of them together rather than per-panel: the same map is on screen in three views at once,
and three different ramps for one temperature field is a way to misread it, not a feature.
"""
function shell_set_colormap(name)
    sh = _sh()
    ok = false
    for c in (sh.msky, sh.star, sh.moll, sh.imsky, sh.immoll)
        ok |= set_colormap!(c, _qmlstr(name))
    end
    # The Mollweide computes its colours from the map, so it has to be redrawn rather than
    # re-tinted; the sky canvases keep their values and `set_colormap!` remaps them in place.
    ok && refresh_model_tab!(sh)
    ok && refresh_data_tab!(sh)
    return ok ? "colormap: $(name)" : "unknown colormap $(name)"
end

"""
    shell_set_colorrange(lo, hi) -> String

Pin every surface view to the colour range `[lo, hi]`, or return them to autoscaling when
either bound is empty or not a number. Bounds given in the wrong order are swapped.
"""
function shell_set_colorrange(lo, hi)
    sh = _sh()
    a = _qmlreal(lo, NaN); b = _qmlreal(hi, NaN)
    if !(isfinite(a) && isfinite(b)) || a == b
        sh.colorrange[] = nothing
        msg = "colour range: auto"
    else
        a > b && ((a, b) = (b, a))
        sh.colorrange[] = (Float32(a), Float32(b))
        msg = Printf.@sprintf("colour range: %.4g – %.4g", a, b)
    end
    refresh_model_tab!(sh)
    refresh_data_tab!(sh)
    refresh_image_tab!(sh)
    return msg
end

"""
    shell_colorrange(tab, view) -> String

`"auto\tlo\thi"` for the canvas `view` (0 orthographic, 1 Mollweide, 2 3-D) of `tab`
(`"model"` or `"image"`): whether the range is autoscaled, and the range that canvas is
drawn with — the autoscaled one included, so the fields can show it. Empty when that canvas
has not been built.
"""
function shell_colorrange(tab, view)
    sh = _sh()
    v = round(Int, _qmlreal(view, 0))
    cs = _qmlstr(tab) == "image" ? (sh.imsky, sh.immoll, sh.imstar) :
                                   (sh.msky, sh.moll, sh.star)
    c = cs[clamp(v + 1, 1, 3)]
    (c === nothing || !hasproperty(c, :cbarlimits)) && return ""
    lo, hi = c.cbarlimits[]
    return Printf.@sprintf("%d\t%.6g\t%.6g", sh.colorrange[] === nothing ? 1 : 0, lo, hi)
end

"""
    shell_reset_zoom() -> String

Every sky view back to the framing it chose for the current star.

Kept as a callback although nothing in the UI calls it any more — the gesture is a RIGHT-CLICK
on the plot (see [`install_interactions!`](@ref)), as in the OITOOLS GUI. A test, or a script
driving the shell, still wants a way to ask for it by name.
"""
function shell_reset_zoom()
    sh = _sh()
    for c in (sh.msky, sh.imsky); reset_zoom!(c); end
    return ""
end

"""
    install_interactions!(sh)

Mouse behaviour that belongs to the PLOT rather than to a button.

Two things, both on every sky canvas, and both read from Makie's own event stream because
neither has a QML equivalent:

  * **Right-click resets the zoom**, as in the OITOOLS GUI. A button for this was the first
    arrangement and it was wrong: the gesture belongs where the eye already is, and a control
    that is only useful while looking at the plot should not be at the top of the panel.
    Only a right-click that did not MOVE counts — a right-drag is Makie's rectangle zoom, and
    resetting at the end of one would undo what the user just did.
  * **The zoom is bounded** on every limit change, which catches the scroll wheel and the
    drag-zoom alike. `Consume(false)` throughout: swallowing these events would take pan and
    zoom with them.

Called once, from `gui()`, after the canvases exist and before the window does.
"""
function install_interactions!(sh::ShellState)
    for c in (sh.msky, sh.imsky)
        # ONE handler owns the wheel and the right-click, and it enforces the bounds BEFORE it
        # moves the view — see `_install_zoom!`. What was here instead was Makie's own
        # `ScrollZoom` plus a `finallimits` listener that pulled the view back afterwards, and
        # the two fought:
        #
        #   * QMLMakie delivers 120 units per wheel click, so Makie's `0.9^120` zoomed by
        #     ~270000x per notch and every click ended pinned against a bound;
        #   * the listener re-entered itself, because the clamp set x and then y and each
        #     write fired it again — 4216 nested calls, ending in "Can't set x limits to the
        #     same value 0.0", swallowed by QMLMakie's render callback;
        #   * and Makie's `:limitreset` returned to `ax.limits[]`, the rectangle the clamp had
        #     just written, so right-click went back to the last zoom rather than to the star.
        #
        # Refusing an out-of-range step is also why no listener is needed to catch the
        # rectangle drag-zoom: that gesture is bounded by what can be drawn on screen, where
        # the wheel was not bounded by anything.
        _install_zoom!(c)
    end
    # The other two view types get it as well: a Mollweide is an `Axis` like the sky and takes
    # the same treatment, and the 3-D scene is an `Axis3` whose zoom cannot be reimplemented
    # here but receives the same oversized event — so it gets a speed that cancels the detent.
    for c in (sh.moll, sh.immoll)
        c === nothing || _install_zoom!(c)
    end
    for c in (sh.star, sh.imstar)
        c === nothing || _install_zoom!(c.scene)
    end
    # The ORBIT view too, and it is the one that most needed it: without a limiter it kept
    # Makie's own `ScrollZoom`, which receives QMLMakie's 120-unit event undivided — a single
    # wheel notch was enough to drive the limits somewhere the renderer could not survive.
    sh.orbitcanvas === nothing || _install_zoom!(sh.orbitcanvas)
    return sh
end

# ── background jobs ─────────────────────────────────────────────────────────────────────

"""
    GuiJob

A long-running Julia call, off the GUI thread, with its output readable while it runs.

Copied from OITOOLS' `GuiJob` because the mechanism is a property of the environment, not of
either package. The split matters: the worker computes and NOTHING else. Every canvas update
stays on the GUI thread, in `job_finish!`, because a GL call from a worker crashes.

The stdout capture needs BOTH `Base.stdout` assignment and a `dup2` on descriptor 1 — the
engines print through both paths — into a `Filesystem.File` rather than the `IOStream` `mktemp`
returns: an IOStream buffers, and a long optimiser trace then arrives in one lump at the end
instead of filling the console as it runs. `redirect_stdout()` onto a pipe HANGS under
`QML.exec()`, where Qt owns the thread and libuv's loop is never pumped.
"""
mutable struct GuiJob
    task::Task
    kind::Symbol
    path::String
    io::Base.Filesystem.File
    saved::Any
    savedfd::Cint
    stop::Base.RefValue{Bool}
    started::Float64
    # MID-RUN PROGRESS, written by the WORKER and read by the GUI thread.
    #
    # One slot holding the latest report rather than a queue: the GUI polls at 200 ms and an
    # engine can report far faster than that, so anything but the most recent is already stale
    # by the time it could be drawn. The worker overwrites; the GUI thread reads whatever is
    # there. A NamedTuple is immutable and the assignment is a single pointer store, so a read
    # racing a write sees the old report or the new one, never half of each — which is all the
    # synchronisation a dropped frame needs.
    #
    # NOTHING HERE TOUCHES GL. The worker only stores; every draw happens in `shell_job_poll`,
    # on the GUI thread, for the same reason `job_succeeded!` does.
    progress::Base.RefValue{Any}
    drawn::Base.RefValue{Int}          # the report number already on screen
end

function start_job!(sh::ShellState, kind::Symbol, f; progress = Ref{Any}(nothing))
    sh.job === nothing || return "a job is already running"
    path, tmp = mktemp(); close(tmp)
    io = Base.Filesystem.open(path, Base.Filesystem.JL_O_WRONLY |
                                    Base.Filesystem.JL_O_CREAT |
                                    Base.Filesystem.JL_O_TRUNC, 0o644)
    saved   = Base.stdout
    savedfd = Sys.isunix() ? ccall(:dup, Cint, (Cint,), 1) : Cint(-1)
    flush(stdout)
    Base.stdout = IOContext(io, :color => false)
    savedfd >= 0 && ccall(:dup2, Cint, (Cint, Cint),
                          Base.cconvert(Cint, Base.fd(io)), Cint(1))
    stop = Ref(false)
    task = Threads.@spawn f(stop)
    sh.job = GuiJob(task, kind, path, io, saved, savedfd, stop, time(), progress, Ref(0))
    sh.jobkind = kind
    return ""
end

"Whatever the running job has printed so far."
function job_output(sh::ShellState)
    j = sh.job
    j === nothing && return ""
    return try; read(j.path, String); catch; ""; end
end

"""
    job_finish!(sh) -> (; ok, result, err, output)

Put stdout back, collect what the worker produced, and delete the job. **GUI thread only.**
"""
function job_finish!(sh::ShellState)
    j = sh.job
    j === nothing && return (; ok = false, result = nothing, err = nothing, output = "")
    Base.stdout = j.saved
    if j.savedfd >= 0
        ccall(:dup2, Cint, (Cint, Cint), j.savedfd, Cint(1))
        ccall(:close, Cint, (Cint,), j.savedfd)
    end
    try; close(j.io); catch; end
    out = try; read(j.path, String); catch; ""; end
    rm(j.path; force = true)
    res, err = nothing, nothing
    try
        res = fetch(j.task)
    catch e
        err = e isa TaskFailedException ? e.task.exception : e
    end
    sh.job = nothing
    return (; ok = err === nothing, result = res, err, output = out)
end

shell_job_running() = (sh = SHELL[]; sh !== nothing && sh.job !== nothing) ? "1" : "0"

# WHICH ENGINES ACTUALLY CHECK THE STOP FLAG. `start_job!` hands its worker a `stop::Ref`, and
# the only thing that reads it is the OBJECTIVE closure `_run_fit` is given — so the two NLopt
# searches notice, and VMLMB, NUTS, Pigeons and the shape path never see the Ref at all:
# `_run_gradient_fit`, `_run_hmc_fit`, `_run_pigeons_fit` and `_run_shape_fit` are called
# without it. MEASURED in the GUI: Stop pressed during a VMLMB fit left the elapsed counter
# climbing and the run went to its full budget.
#
# Even where it IS read it is a soft stop rather than an abort: the objective starts returning
# 1e30, so NLopt keeps searching a now-constant function and burns the rest of its budget
# quickly instead of returning at once.
#
# Until each runner threads the Ref into its optimiser this says what it can and cannot do,
# rather than answering "stopping…" to every engine and leaving the user watching a counter.
const STOP_AWARE = (:neldermead, :bobyqa)

# WHICH ENGINES REACH THE QUADRATURE KERNEL AND WHICH DO NOT. This is the whole reason a fit
# selects its own kernel rather than inheriting the panel's.
#
# `:nufft` exists in ONE place — `fused_cvis_parts`. So:
#
#   * Nelder-Mead, BOBYQA, Nautilus and every binary fit evaluate the objective closure in
#     `shell_fit_start`, which calls `parametric_chi2` → `fused_cvis` → that branch. They are
#     already on the fastest kernel there is (2.6/3.8/4.1 ms at HEALPix 3/4/5) and want
#     nothing done to them.
#   * `:t3` — ROTIR's own type 3 — is faster than `:nufft` everywhere and more accurate, and
#     it is the only quadrature route with an adjoint. It is still NOT what a gradient fit
#     wants: the shape gradient also calls `compute_adjoint_vertices!`, which computes
#     ∂/∂vertex rather than a transpose and therefore has no type-3 form, so selecting `:t3`
#     drops that kernel to the plain reference — MEASURED 3.07/12.0/48.7 ms at HEALPix 3/4/5
#     against `:turbo`'s 0.41/1.57/6.37. Hence `:turbo` below, not `:t3`.
#   * VMLMB (`fit_parametric` → `build_parametric_logπ`), NUTS and Pigeons all evaluate
#     `interferometric_chi2`, and the shape route evaluates `shape_chi2_fg!`; both go through
#     `_cvis_forward!` and the two adjoint kernels, NONE of which has a `:nufft` branch. Left
#     alone they run the plain-Julia reference. MEASURED on one lam And epoch, the three
#     kernels one gradient evaluation runs: 11.4 ms → 0.7 ms at nside 3, 32.3 → 2.2 at 4,
#     124.9 → 9.0 at 5. A 20-draw NUTS fit made 2240 gradient evaluations, so that is ~24 s
#     against a one-time 1.8 s load.
#
# `:gradient` covers both gradient routes — `_run_gradient_fit` and `_run_shape_fit` — since
# `shell_fit_start` picks between them on `gradient_fit_kind` after this point.
const KERNEL_TURBO_METHODS = (:gradient, :hmc, :pigeons)

"""
    _fit_kernel(sh, meth) -> Union{Symbol,Nothing}

Which forward kernel this engine should run under, or `nothing` to leave the process default
alone. See [`KERNEL_TURBO_METHODS`](@ref).

`:scalar` — the reference — is returned rather than `nothing` wherever the console tells the
user the fit runs the reference kernel, so that the claim is true BY CONSTRUCTION. `nothing`
means "no scope at all", under which the worker runs whatever `POLYFT_BACKEND` happens to be,
normally `:t3`. There is no live discrepancy today, because the three
[`KERNEL_TURBO_METHODS`](@ref) evaluate `interferometric_chi2` / `shape_chi2_fg!`, neither of
which has a `:t3` branch — but that is a property of a tuple someone may edit, not of this
function, and adding `:neldermead` to it would silently make the message false.

A FORCED selection wins. Someone who picked `:scalar` from the panel is cross-checking a
result that looked wrong, and silently running something else under their fit is exactly the
thing that would make that check worthless.
"""
function _fit_kernel(sh::ShellState, meth::Symbol)
    sh.kernel_auto[] || return nothing
    meth in KERNEL_TURBO_METHODS || return nothing
    # AND NOT IF IT IS KNOWN NOT TO WORK HERE — a bundle, where LoopVectorization is absent by
    # necessity (see `TURBO_LOAD_FAILED`). Answering without letting the load fail again keeps
    # the console clean. `:scalar`, not `nothing`: the fit then genuinely runs the reference
    # kernel rather than the process default, which is what every comment here claims.
    TURBO_LOAD_FAILED[] && return :scalar
    return :turbo
end

function shell_job_stop()
    sh = _sh()
    j = sh.job
    j === nothing && return "nothing running"
    j.stop[] = true
    sh.jobkind === :fit || return "stopping…"
    m = sh.jobmethod
    m in STOP_AWARE &&
        return "stopping $(fit_method_short(m)) — it checks between evaluations"
    return "$(fit_method_short(m)) cannot be interrupted yet; it will run to its budget"
end

"""
    shell_job_poll() -> String

Called by a QML timer every 200 ms. `running\telapsed\toutput`.

When the worker has finished, this is where the result is taken up and the canvases redrawn —
on the GUI thread, which is the only place a GL call may happen.
"""
function shell_job_poll()
    # Driven by a QML Timer, so reaching this line proves the event loop is turning and a
    # redraw can complete. See `GUI_LIVE`.
    GUI_LIVE[] = true
    sh = _sh()
    j = sh.job
    j === nothing && return "0\t0\t\t"
    if istaskdone(j.task)
        kind = j.kind
        r = job_finish!(sh)
        isempty(r.output) || console!(sh, r.output)
        if r.ok
            job_succeeded!(sh, kind, r.result)
        else
            sh.status = "$(kind) failed: $(sprint(showerror, r.err))"
            console!(sh, sh.status)
        end
        return "0\t0\t\t"
    end
    return Printf.@sprintf("1\t%.1f\t%s\t%s", time() - j.started, job_output(sh),
                           draw_progress!(sh, j))
end

"""
    draw_progress!(sh, job) -> String

Draw the running job's latest report and describe it in one line. **GUI thread only.**

This is the whole of "watch it converge": the worker leaves a map in the job's slot, and the
poll — which already runs on the GUI thread every 200 ms — puts it on the canvases. Redrawing
is skipped when the report has not changed, so a job that reports slower than the poll costs
one comparison per tick rather than a render pass.

A failed draw disables further ones instead of repeating five times a second into the console.
"""
function draw_progress!(sh::ShellState, j::GuiJob)
    p = j.progress[]
    p === nothing && return ""
    line = Printf.@sprintf("%d evaluations · f = %.6g", p.n, p.f)
    p.n == j.drawn[] && return line
    j.drawn[] = p.n
    try
        show_reconstruction!(sh, p.star, p.x;
                             title = Printf.@sprintf("running — %d evals", p.n))
    catch err
        console!(sh, "live map failed, continuing without it: $(sprint(showerror, err))")
        j.progress[] = nothing
    end
    return line
end

# ── fitting ─────────────────────────────────────────────────────────────────────────────

"""
    FIT_METHODS

What the Fit button offers. `:neldermead` and `:bobyqa` are NLopt's derivative-free local
searches; `:nautilus` and `:pigeons` sample and give posteriors rather than a point.

Nelder–Mead is listed first because it is the one that copes with a bad starting point, but it
is a local search on an easy landscape — for anything with a real degeneracy a sampler is the
honest answer and the extra cost is the price of knowing that. `:pigeons` goes further and is
the only one here that sees a MULTIMODAL posterior at all.
"""
const FIT_METHODS = [
    # key            short            what the combo shows
    ("gradient",   "VMLMB",      "VMLMB (quasi-Newton, analytic gradient)"),
    ("neldermead", "Nelder–Mead", "Nelder–Mead (NLopt, local)"),
    ("bobyqa",     "BOBYQA",     "BOBYQA (NLopt, local, quadratic model)"),
    ("hmc",        "NUTS",       "NUTS (Hamiltonian, uses the analytic gradient)"),
    ("nautilus",   "Nautilus",   "Nautilus (importance nested sampling, with evidence)"),
    ("pigeons",    "Pigeons",    "Pigeons (parallel tempering, variational — multimodal)"),
    ("ess",        "ESS",        "ESS (ensemble slice sampling, gradient-free)"),
]

"""
    fit_method_short(method) -> String

The ALGORITHM's name for a table column: `:gradient` reads "VMLMB", `:hmc` reads "NUTS".

The fits history used to print the key, so a row said "gradient" — which names ROTIR's
internal route to the optimiser rather than the optimiser, and is the one entry in the table
whose key is not the algorithm's name. An unrecognised method falls back to its key, since a
fit loaded from a newer session should still say something.
"""
function fit_method_short(method)
    k = String(method)
    i = findfirst(e -> e[1] == k, FIT_METHODS)
    return i === nothing ? k : FIT_METHODS[i][2]
end
"""
    FIT_BUDGETS

What the one budget box means for each engine, and what a sensible value of it is.

THE NUMBER IS THE ENGINE'S OWN UNIT. The box used to say "evaluations" whatever was selected,
carry 2000 whatever was selected, and be divided down inside each runner — `÷50` for VMLMB
iterations, `÷7` for NUTS draws, `log2(n) − 5` for Pigeons rounds. So one typed number meant
five different things, none of them the one on the label, and switching engine silently
changed the size of the run by two orders of magnitude in either direction. Each engine now
names its own unit and carries its own default, and the runners take the number as given.

`lo`/`hi`/`step` are what the spin box needs; `tip` is its tooltip.
"""
const FIT_BUDGETS = Dict{String,NamedTuple}(
    # NLopt counts calls to the objective, and a derivative-free simplex needs thousands.
    "neldermead" => (label = "evaluations", default = 2000, lo = 50, hi = 200_000, step = 100,
                     tip = "objective evaluations (NLopt maxeval)"),
    "bobyqa"     => (label = "evaluations", default = 2000, lo = 50, hi = 200_000, step = 100,
                     tip = "objective evaluations (NLopt maxeval)"),
    # A gradient step is worth ~50 simplex evaluations, which is why this default is small.
    "gradient"   => (label = "iterations", default = 200, lo = 10, hi = 5000, step = 10,
                     tip = "VMLMB iterations; each one costs a value and a gradient"),
    # Kept after warm-up. A quarter again are spent adapting on top of this.
    "hmc"        => (label = "draws", default = 400, lo = 50, hi = 20_000, step = 50,
                     tip = "posterior draws kept; adaptation adds 75% more on top"),
    # The effective sample size, which is what a nested run's cost actually scales with.
    "nautilus"   => (label = "eff. samples", default = 2000, lo = 200, hi = 50_000, step = 200,
                     tip = "effective posterior sample size; the live population follows it"),
    # 2^rounds scans, so this is a small number and every increment doubles the run.
    "pigeons"    => (label = "rounds", default = 8, lo = 6, hi = 12, step = 1,
                     tip = "2^rounds scans — each extra round DOUBLES the run. The " *
                           "variational reference starts tuning at round 2, so more rounds " *
                           "also buy a better-fitted reference, not just more scans."),
    # One sweep moves every walker once, and the ensemble is 2*nfree walkers, so a sweep costs
    # ~5 objective calls per walker. Chen (2025) gives O(d) autocorrelation for this family, so
    # a few thousand sweeps is the right order for the handful of parameters a model fit has.
    "ess"        => (label = "sweeps", default = 3000, lo = 200, hi = 100_000, step = 200,
                     tip = "sweeps kept after an equal burn-in. Each moves all 2*nfree " *
                           "walkers once, at roughly 5 objective calls per walker."),
)

"The budget spec for one method, falling back to the NLopt one for anything unlisted."
fit_budget(method) =
    get(FIT_BUDGETS, String(method),
        (label = "evaluations", default = 2000, lo = 50, hi = 200_000, step = 100,
         tip = "objective evaluations"))

"""
    shell_fit_budget(method) -> String

`label\tdefault\tlo\thi\tstep\ttip` for the budget box, so the panel relabels and
re-ranges it when the engine changes instead of carrying one unit across all of them.
"""
function shell_fit_budget(method)
    b = fit_budget(_qmlstr(method, "neldermead"))
    return join((b.label, b.default, b.lo, b.hi, b.step, b.tip), "\t")
end

# NO `:ultranest`. It is not gated here, it is absent: the GUI is Python-free by construction,
# and a listed method is a method something in this process can load PythonCall to satisfy.
# `:nautilus` is the same kind of sampler in pure Julia — importance nested sampling, an
# evidence, uncertainties — and agrees with UltraNest to 1.4e-6 mas on a λ And sphere fit.
# UltraNest is still there for a script: `using PythonCall` loads ROTIRUltraNestExt and
# `fit_sphere_ld(...; method = :ultranest)` works exactly as before.

"""
    shell_fit_methods() -> String

`key\tlabel` per line. A sampler is listed only when the session can actually run it —
Nautilus needs `using Nautilus`, Pigeons needs its five imports — because offering one that
then fails on click is worse than not offering it, and the schema-driven form has no place to
explain the difference. The launcher loads both, so in a normal session both are there.

UltraNest is not in the list at all; see `FIT_METHODS`.
"""
function shell_fit_methods()
    sh = SHELL[]
    m = sh === nothing ? nothing : current_model(sh.session)
    rows = String[]
    for (k, _short, v) in FIT_METHODS
        k == "nautilus" && !nautilus_available() && continue
        # Tempering runs on the same parametric log-posterior NUTS does, so it is offered in
        # the same place: the rapid rotator, with the gradient stack loaded.
        # Tempering runs on the same three log-posteriors NUTS does, and needs LESS: with the
        # slice explorer it uses no gradient at all. It was offered on the rapid rotator alone
        # only because `_fit_pigeons` was hardcoded to that surface's θ.
        k == "pigeons" && !(pigeons_available() && m !== nothing &&
                            m.surface_type in (0, 1, 2)) && continue
        # NUTS needs the analytic gradient, so it is offered exactly where the gradient path
        # is: the rapid rotator, with Zygote loaded.
        # All three non-Roche surfaces have a differentiable log-posterior now:
        # `build_sphere_logπ`, `build_ellipsoid_logπ` and `build_parametric_logπ`.
        k == "hmc" && !(hmc_available() && m !== nothing &&
                        m.surface_type in (0, 1, 2)) && continue
        # Offered only where the gradient is CONSISTENT with the objective being minimised.
        # See `gradient_fit_kind` — this is not the same set as "has analytic shape
        # derivatives", and the difference is measurable.
        k == "gradient" && (m === nothing || gradient_fit_kind(m) === :none) && continue
        # ESS is pure Julia and lives in the package, so it needs no optional import. It does
        # need N = 2d walkers (Chen 2025, arXiv:2505.02987) — fine for a model fit's handful of
        # parameters, which is exactly the regime it is offered for.
        push!(rows, "$(k)\t$(v)")
    end
    return join(rows, "\n")
end

"Whether ROTIRZygoteExt has loaded, i.e. whether `fit_parametric` exists."
fit_parametric_available() = !isempty(methods(ROTIR.fit_parametric))

"""
    _fit_ok(p) -> Bool

Whether a derivative-free objective should accept `p` at all: buildable AND inside its
plausible ranges. The range half is a soft barrier — it is what keeps a Nelder-Mead search
from wandering into `beta = 12` — and is deliberately NOT applied by the display and compute
paths, which gate on `validate_star_params` alone.
"""
_fit_ok(p) = isempty(validate_star_params(p)) && isempty(advise_star_params(p))

"""
    gradient_fit_kind(model) -> :parametric | :shape | :none

Which analytic-gradient fit, if any, is CONSISTENT with this model.

The distinction that matters is not "does ROTIR have shape derivatives" — it has them for
surface types 0, 1 and 2 (src/shape_gradient.jl), and `joint_reconstruct_oi` uses all three.
It is whether the gradient matches the objective a PARAMETRIC fit minimises, where the
temperature map is pinned by the model rather than free:

  * `:parametric` — surface type 2 with Zygote loaded. `build_parametric_logπ` differentiates
    the von Zeipel map as well as the geometry, so `beta`, `ld1`, `ld2` and `tpole` are
    fittable and the gradient is exact.
  * `:shape` — surface type 0. Its map is `tpole` everywhere and does not depend on any shape
    parameter, so holding it fixed inside `shape_chi2_fg!` is exact rather than an
    approximation.
  * `:none` — everything else. Type 1's map is normalised on `radius_x` and type 2's on
    `rpole`/`frac_escapevel`, so a shape-only gradient omits a real `∂map/∂θ` term.

That omission is not academic. MEASURED on six λ And epochs with four free shape parameters,
a shape-only gradient on type 2 converged in 10 s to χ²ᵣ = 13.43 without moving `inclination`
or `position_angle` at all, while Nelder–Mead took 220 s and reached χ²ᵣ = 12.18 — the
inconsistent gradient stalls, and a "fast" button that quietly returns a worse fit is worse
than no button. So it is not offered there, and the panel says why.
"""
function gradient_fit_kind(m)
    m === nothing && return :none
    m.surface_type == 2 && fit_parametric_available() && return :parametric
    m.surface_type == 0 && return :shape
    # THE ELLIPSOID, now over its whole θ. `shape_chi2_fg!` used to hold the temperature map
    # fixed while the radii moved, which is the inconsistency this used to refuse — measured
    # at 6.6 % against finite differences of the objective, and 13 % in `radius_y`. The map's
    # analytic derivatives now go into the gradient
    # (`temperature_map_vonZeipel_ellipsoid_derivs`), and the whole θ agrees with FD to 1e-9.
    m.surface_type == 1 && return :shape
    return :none
end

function _no_gradient_reason(m)
    m.surface_type == 3 &&
        return "no analytic gradient for a Roche surface; use :neldermead or :nautilus"
    m.surface_type == 2 &&
        return "the rapid rotator's gradient needs Zygote (its temperature map is part of " *
               "the model): add `using Zygote` to this session"
    return "surface_type $(m.surface_type)'s temperature map depends on its shape " *
           "parameters, and only the geometry has analytic derivatives — the mismatch " *
           "stalls the optimiser. Use :neldermead or :nautilus."
end

"""
    SHAPE_THETA

ROTIR parameter name -> position in `shape_chi2_fg!`'s theta vector, per surface type.

`shape_chi2_fg!` (src/shape_gradient.jl) has ANALYTIC derivatives of the projected geometry and
the face normals for surface types 0, 1 and 2 — so the gradient path is not rapid-rotator-only,
which is what the parametric path (`fit_parametric`) happens to be. What type 2 gets in
addition is a differentiable TEMPERATURE MAP (`vonzeipel_map`'s rrule), which is why `beta`,
`ld1`, `ld2` and `tpole` are fittable there and not here: with a free map the limb-darkening
coefficients are largely degenerate with it, and `shape_chi2_fg!` holds them fixed by design.
"""
const SHAPE_THETA = Dict(
    0 => Dict(:radius => 1, :inclination => 2, :position_angle => 3),
    1 => Dict(:radius_x => 1, :radius_y => 2, :radius_z => 3,
              :inclination => 4, :position_angle => 5),
    2 => Dict(:rpole => 1, :frac_escapevel => 2, :inclination => 3, :position_angle => 4))

# DERIVED from `ParametricLayout` rather than written out, so it cannot drift from the θ the
# fit actually takes. Built from the `tpole_free` layout because that is the longest legacy
# one and the GUI always passes the full vector, moving only the named entries.
# The WIDEST layout — everything the gradient path CAN fit, which is what this table is for:
# it answers "could this parameter be freed at all", and the per-fit layout (built from what
# the user actually freed) answers "which slot is it in". Narrower here and `d`/`vgamma` would
# show in the form but be refused on free, which is exactly the trap this fixes.
"""
    PARAMETRIC_THETA

ROTIR parameter name -> position in `fit_parametric`'s theta vector.

`fit_parametric` takes a FIXED 7- (or 8-) element vector,
`[rpole, omega, inc, PA, beta, ld1, ld2]` (+ `tpole`), which is the rapid rotator's parameter
set under different names — `omega` is `frac_escapevel`, `inc` is `inclination`, `PA` is
`position_angle`. The GUI's model is "whichever fields are marked free", so this table is what
lets the two meet, and a free parameter with no entry here is what makes the fast path
inapplicable.
"""
const PARAMETRIC_THETA = let L = ROTIR.parametric_layout(; tpole_free = true,
                                                           distance_free = true,
                                                           vgamma_free = true)
    Dict(L.fields[i] => L.names[i] for i in eachindex(L.names))
end

"""
    _fit_layout(snap, names) -> ParametricLayout

The parameter vector for ONE fit: the rapid rotator's base θ plus exactly the optional
entries the user freed, and with the limb-darkening coefficients dropped when `ldtype` does
not read them.

Built per fit rather than once, because the layout IS the choice of what is being fitted —
`d` is only in θ if the user freed it, and under `ldtype = 0` (a model-atmosphere intensity)
`ld1`/`ld2` leave θ altogether.
"""
_fit_layout(snap, names) =
    ROTIR.parametric_layout(; tpole_free     = true,
                              ldtype        = round(Int, get(snap.params, :ldtype, 3.0)),
                              distance_free = :d in names,
                              vgamma_free   = :vgamma in names)

# Panel parameter symbol -> its slot in `L`. Direct on the FIELD, with no detour through the
# θ alias, so a layout that omits an entry simply has no index for it.
function _theta_slot(L, n)
    j = findfirst(==(n), L.fields)
    j === nothing && error("$(n) is not a parameter of the gradient path for this model " *
                           "(θ is $(L.names)); if it is `ld1`/`ld2`, note ldtype = 0 " *
                           "removes them because the intensity provider owns limb darkening")
    return j
end

# ── fitting a binary ────────────────────────────────────────────────────────────────────
#
# A binary's parameter vector spans THREE places: the primary's free set, the secondary's, and
# the secondary's position. They are kept as `(component, name)` pairs rather than flattened
# into one dictionary because `radius` means a different number on each star — the primary and
# the secondary of Spica differ by a factor of two — and a single namespace would silently tie
# them together.
#
# Component 1 is the primary, 2 the secondary, 3 the position. Only the label differs; the
# fitter sees a plain vector either way.
_fit_component_of(n::Symbol) = n in (:pos_x, :pos_y, :pos_z) ? 3 : 2

"The `(component, name)` pairs a fit will move, primary first."
function binary_fit_names(m)
    c = m.companion
    # An orbital ELEMENT is not a coordinate while the Orbit tab is what places the secondary.
    # The panel shows those rows read-only for that reason, and `_sync_orbit_from_tab!`
    # overwrites on the next build whatever a fit had put there — so a free one would be a
    # direction the χ² only appears to respond to. They are fitted on the Orbit tab instead.
    frozen = c !== nothing && c.place === :orbit
    # The orbital ELEMENTS, and the ORIENTATION the orbit derives (`ORBIT_ORIENTED`): both are
    # rewritten by `_sync_orbit_from_tab!` on the next build, so a fit that moved them would be
    # moving a direction the χ² only appears to respond to.
    movable(n) = !(frozen && (haskey(ORBIT_ELEMENT_OF, n) || n in ORBIT_ORIENTED))
    out = Tuple{Int,Symbol}[(1, n) for n in sort(collect(m.free)) if movable(n)]
    c === nothing && return out
    for n in sort(collect(c.free))
        # A position is only a coordinate when a fixed offset is what places the secondary; if
        # the orbit does, these numbers are not read and freeing them would add flat directions.
        (_fit_component_of(n) == 3 && c.place !== :offset) && continue
        movable(n) || continue
        push!(out, (_fit_component_of(n), n))
    end
    return out
end

"How a `(component, name)` pair reads in the results table."
fit_label(comp::Int, n::Symbol) = comp == 2 ? "2:" * String(n) : String(n)

"""
    _fit_adopt!(m, label, v) -> Bool

Write one of a recorded fit's values back into model `m`, returning whether it landed.

The label is what the results table shows, which is [`fit_label`](@ref)'s output: a bare name
for the primary, `"2:name"` for a companion, and one of `POSITION_PARAMS` for the offset. So
this is `fit_label` inverted, plus the check that makes it safe — the field has to EXIST on the
target. Adopting a sphere's fit into a rapid rotator must not invent `radius` on a model whose
`surface_type` never reads it; it reports that name as skipped instead.
"""
function _fit_adopt!(m, label::Symbol, v::Real)
    str = String(label)
    if startswith(str, "2:")
        n = Symbol(str[3:end])
        m.companion === nothing && return false
        haskey(m.companion.params, n) || return false
        m.companion.params[n] = Float64(v)
        return true
    end
    n = Symbol(str)
    if haskey(m.params, n)
        m.params[n] = Float64(v)
        return true
    end
    # The companion's PLACEMENT, which is neither star's parameter but a property of the pair.
    i = findfirst(q -> q[1] === n, POSITION_PARAMS)
    if i !== nothing && m.companion !== nothing
        o = collect(m.companion.offset)
        o[POSITION_PARAMS[i][3]] = Float64(v)
        m.companion.offset = (o[1], o[2], o[3])
        return true
    end
    return false
end

"Current value, lower and upper bound of one `(component, name)` pair."
function _fit_triple(m, (comp, n))
    if comp == 1
        return (m.params[n], m.bounds[n][1], m.bounds[n][2])
    end
    c = m.companion
    comp == 2 && return (c.params[n], get(c.bounds, n, (-Inf, Inf))[1],
                         get(c.bounds, n, (-Inf, Inf))[2])
    i = findfirst(q -> q[1] === n, POSITION_PARAMS)
    lo, hi = get(c.bounds, n, POSITION_BOUND)
    return (c.offset[POSITION_PARAMS[i][3]], lo, hi)
end

"Write one trial value back into whichever component it belongs to."
function _fit_put!(m, (comp, n), v)
    if comp == 1
        m.params[n] = v
    elseif comp == 2
        m.companion.params[n] = v
    else
        i = findfirst(q -> q[1] === n, POSITION_PARAMS)
        o = collect(m.companion.offset); o[POSITION_PARAMS[i][3]] = v
        m.companion.offset = (o[1], o[2], o[3])
    end
    return m
end

"""
    FitHealth

Evaluation counters for a GUI fit objective, so a run in which *every* evaluation failed is
refused rather than reported as a successful fit.

This exists because of a measured failure, not a hypothetical one. Before `:t3` became
re-entrant, 40 000 of 40 000 likelihood calls threw; each objective's bare `catch` turned them
into the 1e30 sentinel; and the panel reported `log(Z) = -5e29` (finite), `N_eff = 1` and the
FLAT PRIOR as the posterior, after exhausting its budget because `f_live` never moves when no
live point has a finite likelihood. Nothing in the result said anything had gone wrong. The
cause is fixed, but the swallowing is not, and it will hide the next one just as well.

Written from several threads — `_fit_nested` evaluates the likelihood concurrently — so the
counters are atomics and the first exception is taken under a lock.

**No I/O happens here.** The objective runs on worker tasks; the report is made once,
afterwards, by [`_fit_health_report!`](@ref).

A `stop[]` return is deliberately NOT counted. Pressing Stop makes the objective return the
sentinel for the rest of the run, and counting those would turn "the user stopped it" into
"every evaluation failed".
"""
mutable struct FitHealth
    nevals    ::Threads.Atomic{Int}
    nsentinel ::Threads.Atomic{Int}
    nthrown   ::Threads.Atomic{Int}
    first_err ::Any
    lk        ::ReentrantLock
end
FitHealth() = FitHealth(Threads.Atomic{Int}(0), Threads.Atomic{Int}(0),
                        Threads.Atomic{Int}(0), nothing, ReentrantLock())

@inline _eval!(h::FitHealth) = (Threads.atomic_add!(h.nevals, 1); nothing)

"A refused trial point: counted, and 1e30 returned so the search keeps going."
@inline function _refused!(h::FitHealth)
    Threads.atomic_add!(h.nsentinel, 1)
    return 1e30
end

"A trial point whose evaluation THREW. The first exception is kept for the report."
function _refused!(h::FitHealth, err)
    Threads.atomic_add!(h.nsentinel, 1)
    Threads.atomic_add!(h.nthrown, 1)
    if h.first_err === nothing
        lock(h.lk) do
            h.first_err === nothing && (h.first_err = err)
        end
    end
    return 1e30
end

"""
    _fit_health_report!(h, meth) -> String

Report what a fit objective swallowed. Returns a note for the status line, or `""`.

Two behaviours, both of which the previous bare `catch` lacked:

  * if ANY evaluation threw, the first exception is printed. The worker's stdout is captured
    to the job's file and surfaced by `job_output`, so a `println` here is the channel that
    reaches the user without touching `sh.console` from a worker thread.
  * if EVERY evaluation failed, this THROWS. `job_finish!` catches it into `err` and the panel
    reports a failed job, instead of a fit whose "result" is the starting point and whose
    evidence is a finite number with no meaning.

A no-op for the gradient, HMC, Pigeons and shape engines, which never call `obj` — `nevals`
stays 0 and there is nothing to report.
"""
function _fit_health_report!(h::FitHealth, meth)
    n, ns, nt = h.nevals[], h.nsentinel[], h.nthrown[]
    n == 0 && return ""
    firstmsg() = h.first_err === nothing ? "" :
                 sprint(showerror, h.first_err; context = :limit => true)
    if nt > 0
        println("fit warning: $nt of $n evaluations of the $meth objective threw; " *
                "first was: ", firstmsg())
    end
    if ns == n
        error("every one of the $n evaluations of the $meth objective failed, so the " *
              "reported result would be the starting point and the evidence meaningless" *
              (nt > 0 ? "; first exception: " * firstmsg() :
               " — all were refused by `validate_star_params`, so check the parameter bounds"))
    end
    ns == 0 && return ""
    return Printf.@sprintf(" (%d of %d evaluations refused)", ns, n)
end

"""
    _binary_objective(snap, tess, data, mjd, tepochs) -> (θ -> χ²)

The criterion a binary fit minimises: the χ² of the PAIR, over whatever the parameter vector
happens to span.

Matrix-free, through the same `binary_observables` route the panel's χ² uses — which is what
makes this feasible at all. `fit_orbit` refuses tessellated components partly because they were
"orders of magnitude slower per likelihood", and that was true while every evaluation built an
`nuv × npix` matrix per component through `setup_oi!`. It is 13 ms now.
"""
function _binary_objective(snap, names, tess, data, mjd, tepochs, orbit, stop, health;
                           weights = (1, 1, 1))
    return function (θ)
        stop[] && return 1e30          # NOT counted: see `FitHealth`
        _eval!(health)
        for (k, nm) in enumerate(names); _fit_put!(snap, nm, θ[k]); end
        # A shared orbital element is written into the OWNER only; the companion's copy is
        # derived and has to follow before its shape is built (see `_sync_shared_orbit!`).
        _sync_shared_orbit!(snap)
        p1 = star_params(snap); p2 = star_params(snap.companion)
        (isempty(validate_star_params(p1)) && isempty(validate_star_params(p2))) ||
            return _refused!(health)
        try
            c = snap.companion
            offs = if c.place === :offset
                fill((c.offset[1], c.offset[2]), length(data))
            else
                bp = orbit_bparams(orbit; binary = snap)
                [orbit_to_rotir_offset(bp, t + 2_400_000.5) for t in mjd]
            end
            s1 = create_star_multiepochs(tess, p1, tepochs; secondary = false)
            s2 = create_star_multiepochs(tess, p2, tepochs; secondary = true)
            x1 = parametric_temperature_map(p1, s1[1]; secondary = false)
            x2 = parametric_temperature_map(p2, s2[1]; secondary = true)
            tot = 0.0
            for i in eachindex(data)
                ph = binary_phase_shift(data[i].uv, offs[i][1], offs[i][2])
                tot += binary_chi2_f(x1, s1[i], x2, s2[i], data[i], ph; weights = weights)
            end
            return isfinite(tot) ? tot : _refused!(health)
        catch err
            return _refused!(health, err)
        end
    end
end

"""
    shell_free_count() -> String

How many parameters a fit would move, as a decimal string.

The panel greys the Fit button on this rather than counting rows itself, because for a binary
the free parameters are spread over three tables — the primary, the secondary and the
positions — and a position is only free when a fixed offset is what places the secondary.
`binary_fit_names` already resolves all of that; counting in QML would be a second, quietly
different answer.
"""
function shell_free_count()
    m = current_model(_sh().session)
    m === nothing && return "0"
    return string(m.companion === nothing ? length(m.free) : length(binary_fit_names(m)))
end

"""
    _sampler_theta_error(model, sh, names, what) -> String

`""` if every free parameter is one this surface's sampler θ contains, and a sentence saying
why not otherwise.

SHARED BY NUTS AND PIGEONS, because they sample the same three log-posteriors and the check is
about the θ layout rather than the algorithm. It is a function rather than two copies for a
specific reason: the copies drifted, and Pigeons spent a while offered on the rapid rotator
alone while NUTS ran on all three.

Refused rather than dropped silently: a parameter you believe is being sampled and is not is
worse than being told you cannot sample it.
"""
function _sampler_theta_error(m, sh, names, what::AbstractString)
    if m.surface_type == 0
        # A SPHERE samples `[radius, ld…]` and nothing else, because that is all a sphere puts
        # into normalised visibilities — its orientation is flat to within the mesh's faceting
        # and its temperature is a scale the normalisation divides out (both measured; see
        # `build_sphere_logπ`).
        θn = ROTIR.sphere_param_names(round(Int, get(m.params, :ldtype, 3.0)))
        bad = filter(n -> !(String(n) in θn), names)
        isempty(bad) ||
            return "$(what) on a sphere samples " * join(θn, ", ") * " — " *
                   join(bad, ", ") * " carries no information in the visibilities"
    elseif m.surface_type == 1
        ldt = round(Int, get(m.params, :ldtype, 3.0))
        # `tpole` only under `:planck`: a uniform scale on the map divides out of the
        # normalised visibility, so under `:linear` it is a flat direction. The panel's
        # intensity law decides, which is why this is checked here and not in the layout.
        tpf = :tpole in names && sh.intensity_model[] === :planck
        # THIS TEST COMES FIRST, and it did not used to. With `tpole` free under `:linear`,
        # `tpf` is false, so `ellipsoid_free_indices` throws "`tpole` is not an ellipsoid
        # parameter" and the `catch` below returned THAT — leaving the sentence underneath
        # unreachable. The lower-level message is accurate and useless: it names the θ layout
        # rather than the reason, and does not say what to do. This one says both.
        (:tpole in names && !tpf) &&
            return "tpole is a pure scale under a linear intensity law and divides out of " *
                   "the normalised visibility — switch the intensity model to planck or fix " *
                   "tpole"
        try
            ROTIR.ellipsoid_free_indices([String(n) for n in names], ldt; tpole_free = tpf)
        catch err
            return "$(what) on an ellipsoid: " * first(split(sprint(showerror, err), '\n'))
        end
    elseif m.surface_type == 2
        bad = filter(n -> !haskey(PARAMETRIC_THETA, n), names)
        isempty(bad) || return "not parameters of the $(what) path: " * join(bad, ", ")
    else
        return "$(what) runs on a sphere, an ellipsoid or a rapid rotator; surface type " *
               "$(m.surface_type) has no differentiable log-posterior"
    end
    return ""
end

"""
    _weighted_ndata(data, w) -> Int

The number of points a χ² weighted by `w = (V², T3amp, T3φ)` sums over: an observable with
weight zero contributes no terms, so it contributes no degrees of freedom either.
"""
_weighted_ndata(data, w) = sum(d -> (w[1] > 0) * d.nv2 + (w[2] > 0) * d.nt3amp +
                                    (w[3] > 0) * d.nt3phi, data; init = 0)

const OBS_NAMES = ("V²", "T3amp", "T3φ")

"`\"\"` when every observable is in, else `\" — V², T3amp only\"` for the console."
_weights_note(w) = all(>(0), w) ? "" :
    " — " * join((OBS_NAMES[k] for k in 1:3 if w[k] > 0), ", ") * " only"

"`w` as the literal a script passes: `(1, 1, 0)`."
_weights_src(w) = "(" * join((isinteger(x) ? string(Int(x)) : string(x) for x in w), ", ") * ")"

_weights_ref(sh, tab) = String(tab) == "image" ? sh.imgweights : sh.fitweights

"""
    shell_obs_weights(tab) -> String

Which observables the `"model"` tab's fits or the `"image"` tab's reconstructions use, as
`"v2\tt3amp\tt3phi"` with each `1` or `0`.
"""
shell_obs_weights(tab) =
    join((x > 0 ? "1" : "0" for x in _weights_ref(_sh(), _qmlstr(tab))[]), "\t")

"""
    shell_set_obs_weights(tab, v2, t3amp, t3phi) -> String

Tick or untick V², T3amp and T3φ for the `"model"` tab's fits or the `"image"` tab's
reconstructions. Unticking all three is allowed here and refused when a run starts, so the
boxes can be changed in any order.
"""
function shell_set_obs_weights(tab, v2, t3amp, t3phi)
    sh = _sh()
    t = _qmlstr(tab)
    on(x) = (x === true || x == 1 || _qmlstr(x) in ("true", "1")) ? 1.0 : 0.0
    w = (on(v2), on(t3amp), on(t3phi))
    _weights_ref(sh, t)[] = w
    what = t == "image" ? "reconstructions" : "fits"
    return any(>(0), w) ? "$(what) use " * join((OBS_NAMES[k] for k in 1:3 if w[k] > 0), ", ") :
                          "$(what) use no observable — tick at least one"
end

"""
    shell_fit(method, budget) -> String

Fit the FREE parameters of the current model against every epoch of the current dataset.

Generic over `parametric_chi2` rather than routed to `fit_sphere_ld` / `fit_ellipsoid_ld`:
those take their own fixed parameter orders per shape, while the GUI's model is "whichever
fields the user marked free", which is exactly what this walks. Tied parameters are re-derived
inside the objective, so a tie holds at every trial point rather than only at the start.

`budget` is in the ENGINE's own unit — evaluations for NLopt, iterations for VMLMB, draws for
NUTS, effective samples for Nautilus, rounds for Pigeons. See [`FIT_BUDGETS`](@ref), which is
also what tells the panel how to label the box.

THE MODEL IS NOT MODIFIED. The result is recorded as a fit and shown in the panel;
[`shell_adopt_fit`](@ref) is what copies it into the model's parameters.
"""
function shell_fit(method, maxeval)
    sh = _sh()
    d = current_dataset(sh.session)
    d === nothing && return "no dataset"
    m = current_model(sh.session)
    m === nothing && return "no model"
    isbin = m.companion !== nothing
    bnames = isbin ? binary_fit_names(m) : Tuple{Int,Symbol}[]
    names = isbin ? Symbol[n for (_, n) in bnames] : sort(collect(m.free))
    isempty(names) && return isbin ?
        "nothing is free — mark a parameter free on either component, or a position" :
        "nothing is free — mark at least one parameter free"
    meth = Symbol(_qmlstr(method))
    # Named before the membership test, so asking for it gets the REASON rather than
    # "unknown method" — it is a method ROTIR has and the GUI declines, not a typo.
    meth === :ultranest &&
        return "the GUI does not offer :ultranest — it would load Python into a process that " *
               "must not have it. Use :nautilus here, or :ultranest from a script."
    meth in (:neldermead, :bobyqa, :nautilus, :gradient, :hmc, :pigeons) ||
        return "unknown method $(method)"
    meth === :nautilus && !nautilus_available() &&
        return "Nautilus needs `using Nautilus` in this session"
    if meth === :pigeons
        pigeons_available() ||
            return "Pigeons needs Pigeons, Distributions, LogDensityProblems, ADTypes and " *
                   "Zygote in this session"
        # THE SAME CHECK NUTS MAKES, through the same function. `_fit_pigeons` used to be
        # hardcoded to the rapid rotator's θ and this branch refused everything else; it
        # dispatches on `model` now, and sharing the validation is what stops the two drifting
        # apart again — which is how Pigeons came to be offered on one surface while NUTS ran
        # on three.
        err = _sampler_theta_error(m, sh, names, "Pigeons")
        isempty(err) || return err
    end
    if meth === :hmc
        hmc_available() ||
            return "NUTS needs AdvancedHMC, LogDensityProblems and Zygote in this session"
        err = _sampler_theta_error(m, sh, names, "NUTS")
        isempty(err) || return err
    end
    if meth === :gradient
        kind = gradient_fit_kind(m)
        kind === :none && return _no_gradient_reason(m)
        known = kind === :parametric ? keys(PARAMETRIC_THETA) : keys(SHAPE_THETA[m.surface_type])
        bad = filter(n -> !(n in known), names)
        isempty(bad) || return "not parameters of the gradient path for this surface type: " *
                              join(bad, ", ") * " (free: " * join(sort(collect(known)), ", ") * ")"
    end

    # A BINARY is fitted by the derivative-free methods only. The gradient path differentiates
    # the single-star parametric model, and NUTS and Pigeons both need `PARAMETRIC_THETA`,
    # which maps one star's parameters onto that model's θ — and there is no binary θ.
    # `binary_chi2_fg` is not that gradient: it differentiates the two SURFACE MAPS, which is
    # what imaging varies, not the shape and orientation parameters this panel fits.
    # Nelder-Mead, BOBYQA and Nautilus need only the criterion, which exists.
    if isbin && !(meth in (:neldermead, :bobyqa, :nautilus))
        return "$(meth) cannot fit a binary: it needs a gradient or the parametric θ, and " *
               "neither exists for two components. Use neldermead, bobyqa or nautilus"
    end

    θ0, lb, ub = if isbin
        tr = [_fit_triple(m, nm) for nm in bnames]
        (Float64[x[1] for x in tr], Float64[x[2] for x in tr], Float64[x[3] for x in tr])
    else
        ([m.params[n] for n in names], [m.bounds[n][1] for n in names],
         [m.bounds[n][2] for n in names])
    end
    # PARENTHESISED. `&&` binds tighter than `||`, so `a || b && return` parses as
    # `a || (b && return)` — and when `a` is true the `||` short-circuits and the return never
    # runs. A starting value BELOW its lower bound was therefore accepted silently, which is
    # the half of this check that matters: a sampler started outside its box transforms to
    # infinity, and an optimiser clamps without saying so.
    if any(θ0 .< lb) || any(θ0 .> ub)
        bad = findfirst(k -> θ0[k] < lb[k] || θ0[k] > ub[k], eachindex(θ0))
        return "$(names[bad]) = $(θ0[bad]) is outside its bounds " *
               "[$(lb[bad]), $(ub[bad])]"
    end

    # Snapshot the model: the worker must not read a struct the GUI thread can edit under it.
    # The COMPANION is copied too, deeply — a shallow copy would share the parameter dictionary
    # the objective writes to on every trial point, which is exactly the struct this snapshot
    # exists to keep away from the GUI thread.
    csnap = m.companion === nothing ? nothing :
        ModelEntry(m.companion.name, m.companion.surface_type, copy(m.companion.params),
                   copy(m.companion.free), copy(m.companion.bounds), copy(m.companion.ties),
                   m.companion.secondary, nothing, m.companion.place, m.companion.offset,
                   m.companion.grid_path)
    snap = ModelEntry(m.name, m.surface_type, copy(m.params), copy(m.free), copy(m.bounds),
                      copy(m.ties), m.secondary, csnap, m.place, m.offset, m.grid_path)
    orbit_snap = deepcopy(sh.orbit)
    mjd = d === nothing ? Float64[] : copy(d.mjd)
    data = d.data
    tepochs = copy(d.tepochs)
    w = sh.fitweights[]
    nd = _weighted_ndata(data, w)
    nd > 0 || return "no data to fit: tick at least one observable the dataset has"
    # THE ENGINE'S OWN UNIT, floored at 1 and not at 10: each runner clamps to its own range
    # (see `FIT_BUDGETS`), and a central floor of 10 would have turned a Pigeons budget of 8
    # ROUNDS into 10 — four times the run the panel asked for.
    it = max(1, Int(maxeval))
    unit = fit_budget(meth).label

    # THE PANEL'S MESH, not a hardcoded one. This built `tessellation_healpix(4)` — nside 16,
    # 3072 tessels — while the mesh box above it read "n = 3, 768 tessels", so every fit did
    # four times the work per likelihood that the panel said it would, and a fit that takes
    # seconds from a script took minutes here. Captured before the job starts: the worker must
    # not read a Ref the GUI thread can change under it.
    nexp = sh.nside_exp[]
    prec = sh.precision[]
    sh.jobmethod = meth
    # THE KERNEL THIS ROUTE WANTS, resolved and LOADED here — on the GUI thread, before the
    # worker exists, because `_ensure_turbo!` must not run on a worker. Captured as a value
    # and scoped over the worker below, so it is snapshotted exactly the way `nexp` and `prec`
    # are: a click on the kernel combo mid-fit must not change the objective function under a
    # running line search, which is what a bare assignment to the global Ref would do.
    kernel = _fit_kernel(sh, meth)
    if kernel === :turbo
        err = _ensure_turbo!(sh)
        if !isempty(err)
            # Not fatal. The reference kernel computes the same thing, only slower, so a
            # failed load costs time rather than the fit.
            console!(sh, err * " — this fit runs the reference kernel instead")
            # `:scalar`, not `nothing`. `nothing` means no scope, under which the worker runs
            # the process default — normally `:t3` — and the line just printed would be false.
            kernel = :scalar
        end
    end
    # EVERY PATH RUNS AT `prec` NOW. Each runner builds its mesh at the precision the box
    # asks for, and every library entry point narrows the Float64 θ the panel hands it to the
    # mesh's type once — `fit_parametric` and `_fit_hmc` both do
    # `collect(eltype(tessels.unit_xyz), θ0)` — so nothing promotes per evaluation. The shape
    # path used to be an exception, forcing Float64 and printing `prec` anyway; it honours the
    # box now, so this line is true for all six engines.
    console!(sh, "fit $(join(names, ", ")) by $(meth), $(it) $(unit), " *
                 "HEALPix level $(nexp) ($(12 * (2^nexp)^2) tessels, $(prec))" *
                 _weights_note(w); kind = :cmd)
    body = function (stop)
        tess = tessellation_healpix(nexp; T = prec)
        health = FitHealth()
        obj = isbin ?
            _binary_objective(snap, bnames, tess, data, mjd, tepochs, orbit_snap, stop, health;
                              weights = w) :
            function (θ)
            stop[] && return 1e30          # NOT counted: see `FitHealth`
            _eval!(health)
            for (k, n) in enumerate(names); snap.params[n] = θ[k]; end
            p = star_params(snap)
            isempty(validate_star_params(p)) || return _refused!(health)
            c = try
                parametric_chi2(p, tess, data, tepochs; weights = w)
            catch err
                return _refused!(health, err)
            end
            return isfinite(c) ? c : _refused!(health)
        end
        best, chi2, extra, post, nevals = if meth === :hmc
            _run_hmc_fit(snap, data, tepochs, names, θ0, lb, ub, nexp, prec, it; weights = w)
        elseif meth === :pigeons
            # Same three-way split as `:hmc` above: the θ layout is the surface's.
            if snap.surface_type == 0
                _run_sphere_pigeons_fit(snap, data, tepochs, names, lb, ub, nexp, prec, it;
                                        weights = w)
            elseif snap.surface_type == 1
                _run_ellipsoid_pigeons_fit(snap, data, tepochs, names, lb, ub, nexp, prec, it;
                                           weights = w)
            else
                _run_pigeons_fit(snap, data, tepochs, names, θ0, lb, ub, nexp, prec, it; weights = w)
            end
        elseif meth !== :gradient
            _run_fit(meth, obj, θ0, lb, ub, names, it)
        elseif gradient_fit_kind(snap) === :parametric
            _run_gradient_fit(snap, data, tepochs, names, θ0, lb, ub, nexp, prec, it; weights = w)
        else
            _run_shape_fit(snap, data, tepochs, names, θ0, lb, ub, nexp, prec, it; weights = w)
        end
        # Say what the objective swallowed, and refuse a run in which everything failed.
        # Placed here rather than inside the objective: this is the worker, whose stdout the
        # job captures, and `sh.console` must not be written from a worker thread.
        healthnote = _fit_health_report!(health, meth)
        if isbin
            for (k, nm) in enumerate(bnames); _fit_put!(snap, nm, best[k]); end
        else
            for (k, n) in enumerate(names); snap.params[n] = best[k]; end
        end
        # The WHOLE result travels back, not a summary of it. `job_succeeded!` turns this into
        # a `FitEntry` on the GUI thread; nothing is reduced here, on the worker, where the
        # draws would otherwise go out of scope and be collected.
        # LABELS for a binary, so "radius" and "2:radius" are told apart in the table and in
        # the fit history — they are different numbers on different stars.
        labels = isbin ? [Symbol(fit_label(c, n)) for (c, n) in bnames] : names
        return (; params = Dict(n => best[k] for (k, n) in enumerate(names)),
                # `binary` carries the (component, name) pairs so the GUI thread can put each
                # value back on the star it came from; without it every result would land on
                # the primary.
                binary = isbin ? bnames : nothing,
                table = _fit_table(labels, best, Dict(Symbol(fit_label(c, n)) => v
                                                      for ((c, n), v) in
                                                      zip(isbin ? bnames :
                                                          [(1, x) for x in names],
                                                          [get(extra, n, NaN) for n in names])
                                                      if isfinite(v)),
                                   ),
                names = labels, best = best, errs = extra, post = post, nevals = nevals,
                method = meth, model = snap.name, surface_type = snap.surface_type,
                chi2 = chi2, ndata = nd,
                status = Printf.@sprintf("fit done: χ²ᵣ = %.4f over %d points", chi2 / nd, nd) *
                         healthnote)
    end
    # The scope wraps the WHOLE worker, so every task the engine starts inside it — the
    # `Threads.@threads` loops in the adjoint kernels included — inherits the choice, and the
    # GUI thread answering the event loop beside it keeps the kernel the panel is set to.
    return start_job!(sh, :fit, kernel === nothing ? body :
                      (stop -> ROTIR.with_polyft_backend(() -> body(stop), kernel)))
end

"""
    _run_gradient_fit(snap, data, tepochs, names, θ0, lb, ub, nexp, prec, maxeval) -> (best, chi2, extra)

The FAST fit: `fit_parametric`, i.e. VMLMB driven by a Zygote gradient of the analytic
parametric model.

Everything else here evaluates `parametric_chi2`, which rebuilds the geometry and re-evaluates
the polygon FT per trial point — the same cost as one likelihood call in a sampler. This path
instead differentiates the model itself (src/parametric_gradient.jl provides the rrules), so a
fit that takes a derivative-free search thousands of evaluations takes a gradient method tens
of iterations. It is what `demos/rapid_rotator_betCas_param_fit.jl` uses.

Restricted to `surface_type = 2` because that is the parameter vector `fit_parametric` is
written for; `shell_fit_methods` only offers it there.
"""
function _run_gradient_fit(snap, data, tepochs, names, θ0, lb, ub, nexp, prec, maxeval;
                           weights = (1, 1, 1))
    # Full theta from the model, INCLUDING the frozen entries: `fit_parametric` wants the whole
    # vector and moves only the indices named in `free`.
    L = _fit_layout(snap, names)
    getp(k) = Float64(get(snap.params, k, 0.0))
    θfull = Float64[getp(f) for f in L.fields]
    lo = copy(L.lower); hi = copy(L.upper)
    free = String[]
    for (k, n) in enumerate(names)
        j = _theta_slot(L, n)
        push!(free, L.names[j])
        # The panel's bounds win over the defaults; the layout carries an infinite upper
        # bound on rpole and tpole, which a box-constrained search accepts but which throws
        # away whatever the user typed.
        lo[j] = lb[k]; hi[j] = ub[k]
    end
    tess = tessellation_healpix(nexp; T = prec)
    base = star_params(snap)
    θ̂, chi2r, info = ROTIR.fit_parametric(data, tess, tepochs, base;
                                          θ0 = θfull, free = free, lb = lo, ub = hi,
                                          tpole_free = true, layout = L, weights = weights,
                                          maxiter = max(10, maxeval), verb = true)
    best = [θ̂[_theta_slot(L, n)] for n in names]
    return best, chi2r * info.ndof, Dict{Symbol,Float64}(), nothing, 0
end

"""
    _run_hmc_fit(snap, data, tepochs, names, θ0, lb, ub, nexp, prec, maxeval) -> (best, chi2, extra)

NUTS over the free parametric parameters, returning the posterior MEDIAN as the point estimate
and half the 16–84 interval as the error.

The median, not the maximum: a Hamiltonian sampler is run for the posterior, and the mode of a
skewed one is not the number to quote beside a symmetric error bar. The panel's budget is the
number of DRAWS TO KEEP; adaptation is three quarters as many again on top, because that is
where the step size and the mass matrix are set and starving it is the commonest way to get
divergences.
"""

function _run_hmc_fit(snap, data, tepochs, names, θ0, lb, ub, nexp, prec, maxeval;
                      weights = (1, 1, 1))
    # A SPHERE has its own, much shorter θ, so it takes its own route. Its orientation and its
    # polar temperature carry no information in normalised visibilities — MEASURED, see
    # `build_sphere_logπ` — so what is sampled is the angular radius and the limb darkening,
    # which is exactly what an interferometer measures for a single star.
    snap.surface_type == 0 &&
        return _run_sphere_hmc_fit(snap, data, tepochs, names, lb, ub, nexp, prec, maxeval;
                                   weights = weights)
    snap.surface_type == 1 &&
        return _run_ellipsoid_hmc_fit(snap, data, tepochs, names, lb, ub, nexp, prec, maxeval;
                                      weights = weights)
    L = _fit_layout(snap, names)
    getp(k) = Float64(get(snap.params, k, 0.0))
    θfull = Float64[getp(f) for f in L.fields]
    lo = copy(L.lower); hi = copy(L.upper)
    free = String[]
    for (k, n) in enumerate(names)
        j = _theta_slot(L, n)
        push!(free, L.names[j])
        lo[j] = lb[k]; hi[j] = ub[k]
    end
    tess = tessellation_healpix(nexp; T = prec)
    base = star_params(snap)
    ndraw = clamp(Int(maxeval), 50, 20_000)
    r = ROTIR._fit_hmc(data, tess, tepochs, base; θ0 = θfull, free = free, lb = lo, ub = hi,
                       weights = weights,
                       tpole_free = true, layout = L, n_samples = ndraw,
                       n_adapt = 3ndraw ÷ 4, verb = true)
    # `_fit_hmc` returns the free parameters in the order `parametric_free_indices` sorts them,
    # which is not the order the panel listed them in.
    order = sortperm([_theta_slot(L, n) for n in names])
    best = zeros(length(names)); errs = Dict{Symbol,Float64}()
    for (slot, k) in enumerate(order)
        best[k] = r.median[slot]
        errs[names[k]] = (r.q84[slot] - r.q16[slot]) / 2
    end
    # `invperm`, not `order`: `order` maps a panel slot to a sampler column, and the columns
    # have to be permuted the other way to land in panel order.
    # `r.chi2` is the χ² AT THE MEDIAN, which `_fit_hmc` and `_fit_pigeons` now return. It used
    # to be `NaN` here, so every sampler row in the fit history showed `NaN` in the χ²ᵣ column
    # — in the one table built for comparing a sphere against a rapid rotator, which cannot be
    # done if half the rows have no number. The evidence is the better comparison where there
    # is one, but a sampler without an evidence still has a χ².
    return best, r.chi2, errs,
           _posterior(r; order = invperm(order),
                      diagnostics = Printf.@sprintf("%d draws, %d divergences",
                                                    size(r.samples, 1), r.divergences)), 0
end

"""
    _run_sphere_hmc_fit(snap, data, tepochs, names, lb, ub, nexp, prec, maxeval)

NUTS on a SPHERE: `θ = [radius, ld…]`, and nothing else.

`names` is what the panel marked free, already checked against the sphere's θ by `shell_fit` —
`inclination` and `position_angle` are flat to within the mesh's own faceting on a sphere, and
`tpole` divides out of the normalised visibility, so neither reaches this far.
"""
function _run_sphere_hmc_fit(snap, data, tepochs, names, lb, ub, nexp, prec, maxeval;
                             weights = (1, 1, 1))
    base = star_params(snap)
    θnames = ROTIR.sphere_param_names(base.ldtype)
    lo, hi = ROTIR.default_sphere_bounds(base.ldtype)
    lo = collect(Float64, lo); hi = collect(Float64, hi)
    θfull = Float64[get(snap.params, Symbol(n), 0.0) for n in θnames]
    free = String[]
    for (k, n) in enumerate(names)
        j = findfirst(==(String(n)), θnames)
        push!(free, θnames[j])
        lo[j] = lb[k]; hi[j] = ub[k]
    end
    tess = tessellation_healpix(nexp; T = prec)
    ndraw = clamp(Int(maxeval), 50, 20_000)
    r = ROTIR._fit_hmc(data, tess, tepochs, base; θ0 = θfull, free = free, lb = lo, ub = hi,
                       weights = weights,
                       model = :sphere, n_samples = ndraw, n_adapt = 3ndraw ÷ 4, verb = true)
    # Sampler columns come back in sorted-index order, not the order the panel listed them.
    order = sortperm([findfirst(==(String(n)), θnames) for n in names])
    best = zeros(length(names)); errs = Dict{Symbol,Float64}()
    for (slot, k) in enumerate(order)
        best[k] = r.median[slot]
        errs[names[k]] = (r.q84[slot] - r.q16[slot]) / 2
    end
    # `r.chi2` is the χ² AT THE MEDIAN, which `_fit_hmc` and `_fit_pigeons` now return. It used
    # to be `NaN` here, so every sampler row in the fit history showed `NaN` in the χ²ᵣ column
    # — in the one table built for comparing a sphere against a rapid rotator, which cannot be
    # done if half the rows have no number. The evidence is the better comparison where there
    # is one, but a sampler without an evidence still has a χ².
    return best, r.chi2, errs,
           _posterior(r; order = invperm(order),
                      diagnostics = Printf.@sprintf("%d draws, %d divergences",
                                                    size(r.samples, 1), r.divergences)), 0
end

"""
    _run_ellipsoid_hmc_fit(snap, data, tepochs, names, lb, ub, nexp, prec, maxeval)

NUTS on a triaxial ellipsoid: `θ = [rx, ry, rz, inc, PA, β, ld1, ld2]`, plus `tpole` when the
intensity law is Planck.

The longest of the three θ vectors, because both halves of the model respond to it — the
projected geometry through the radii and the orientation, the temperature map through the radii
and β. `tpole` joins only under `:planck`: under `:linear` the map is `tpole·f` and
`cvis = F/flux` divides the scale out, so it would be a flat direction (measured at 8.8e-16 —
machine zero — against 1.1e-2 under Planck).
"""
function _run_ellipsoid_hmc_fit(snap, data, tepochs, names, lb, ub, nexp, prec, maxeval;
                                weights = (1, 1, 1))
    sh = _sh()
    base = star_params(snap)
    tpf = :tpole in names && sh.intensity_model[] === :planck
    θnames = ROTIR.ellipsoid_param_names(; tpole_free = tpf)
    lo, hi = ROTIR.default_ellipsoid_bounds(; tpole_free = tpf)
    lo = collect(Float64, lo); hi = collect(Float64, hi)
    θfull = Float64[Float64(get(snap.params, Symbol(n), 0.0)) for n in θnames]
    free = String[]
    for (k, n) in enumerate(names)
        j = findfirst(==(String(n)), θnames)
        push!(free, θnames[j])
        lo[j] = lb[k]; hi[j] = ub[k]
        θfull[j] = clamp(θfull[j], lb[k], ub[k])
    end
    tess = tessellation_healpix(nexp; T = prec)
    ndraw = clamp(Int(maxeval), 50, 20_000)
    r = ROTIR._fit_hmc(data, tess, tepochs, base; θ0 = θfull, free = free, lb = lo, ub = hi,
                       weights = weights,
                       model = :ellipsoid, tpole_free = tpf,
                       intensity_model = sh.intensity_model[],
                       band = sh.band[] > 0 ? sh.band[] : nothing,
                       n_samples = ndraw, n_adapt = 3ndraw ÷ 4, verb = true)
    order = sortperm([findfirst(==(String(n)), θnames) for n in names])
    best = zeros(length(names)); errs = Dict{Symbol,Float64}()
    for (slot, k) in enumerate(order)
        best[k] = r.median[slot]
        errs[names[k]] = (r.q84[slot] - r.q16[slot]) / 2
    end
    # `r.chi2` is the χ² AT THE MEDIAN, which `_fit_hmc` and `_fit_pigeons` now return. It used
    # to be `NaN` here, so every sampler row in the fit history showed `NaN` in the χ²ᵣ column
    # — in the one table built for comparing a sphere against a rapid rotator, which cannot be
    # done if half the rows have no number. The evidence is the better comparison where there
    # is one, but a sampler without an evidence still has a χ².
    return best, r.chi2, errs,
           _posterior(r; order = invperm(order),
                      diagnostics = Printf.@sprintf("%d draws, %d divergences",
                                                    size(r.samples, 1), r.divergences)), 0
end

"""
    PIGEONS_GUI_OPTS

Tempering options shared by all three Pigeons runners, in ONE place.

# A VARIATIONAL REFERENCE, and why it is not optional here

Pigeons anneals FROM a reference TO the posterior. The default fixed reference is a wide
normal on the unconstrained coordinates, and for a ROTIR posterior that reference holds
essentially no mass where the target lives — the box spans radius 0.05–20 mas while the
posterior is ~6e-4 mas wide. The whole communication barrier then sits in the first
infinitesimal step away from the reference, and no number of chains bridges it: extra rungs
resolve a narrow singularity instead of crossing it.

MEASURED on α Cen A (sphere, radius + ld1, V² + T3amp, 2^8 scans):

    configuration                Λ        round trips    local barrier at β≈0
    fixed reference, 10 chains   5.482          2              50 918
    fixed reference, 30 chains   6.082          3             236 792
    variational,     10 chains   0.659         62                1.26
    stabilized,      10+10       5.162         70                  —

The variational reference collapses the reference-end barrier by a factor of 40 000 and turns
2 round trips into 62 for 1.75x the wall time — 10 s per round trip against 176 s. All four
agree on the posterior, so the fixed reference was not WRONG, only unable to demonstrate that
it had explored anything.

# Why STABILIZED rather than plain variational

`n_chains_variational > 0` keeps a fixed AND a variational reference on separate legs
(Surjanovic et al. 2022). It costs ~1.7x plain variational and got the most round trips here.
It is the default because the GUI offers Pigeons specifically FOR a multimodal posterior, and
Pigeons' own documentation warns that fitting a single Gaussian reference "can cause
catastrophic forgetting of modes" — the exact failure this method exists to avoid. A sampler
that reports beautiful round trips while having collapsed onto one tooth is worse than one
that reports zero.

# first_tuning_round

`GaussianReference` does nothing until that round, and its default is 6. With the budget's
minimum of 6 rounds the default would tune only in the last round and return a result
BIT-IDENTICAL to a fixed reference — a silent no-op, not an error. 2 leaves every round after
the cheap ones to adapt in.
"""
const PIGEONS_GUI_OPTS = (n_chains = 8, n_chains_variational = 8,
                          variational_first_round = 2, explorer = :slice)

"""
    _run_sphere_pigeons_fit(snap, data, tepochs, names, lb, ub, nexp, prec, maxeval)

Tempering on a SPHERE's θ, `[radius, ld…]`.

Pigeons was offered on the rapid rotator alone because `_fit_pigeons` was hardcoded to that
surface's θ, not because tempering needs anything the sphere cannot give — with the default
slice explorer it needs no gradient at all, which is a weaker requirement than the NUTS that
has been offered here all along. What tempering adds over NUTS on a short, near-unimodal θ is
not mixing but the EVIDENCE: comparing a sphere against an ellipsoid against a rapid rotator
by `log(Z)` needs a `log(Z)` from each, and Nautilus was the only other source (measured at
836 s against Nelder-Mead's 4.4 s on a lam And sphere).
"""
function _run_sphere_pigeons_fit(snap, data, tepochs, names, lb, ub, nexp, prec, maxeval;
                                 weights = (1, 1, 1))
    base = star_params(snap)
    θnames = ROTIR.sphere_param_names(base.ldtype)
    lo, hi = ROTIR.default_sphere_bounds(base.ldtype)
    lo = collect(Float64, lo); hi = collect(Float64, hi)
    θfull = Float64[get(snap.params, Symbol(n), 0.0) for n in θnames]
    free = String[]
    for (k, n) in enumerate(names)
        j = findfirst(==(String(n)), θnames)
        push!(free, θnames[j])
        lo[j] = lb[k]; hi[j] = ub[k]
        θfull[j] = clamp(θfull[j], lb[k], ub[k])
    end
    tess = tessellation_healpix(nexp; T = prec)
    nr = clamp(Int(maxeval), 6, 12)
    r = ROTIR._fit_pigeons(data, tess, tepochs, base; θ0 = θfull, free = free, lb = lo, ub = hi,
                           weights = weights,
                           model = :sphere, n_rounds = nr, PIGEONS_GUI_OPTS...,
                           multithreaded = Threads.nthreads() > 1, verb = true)
    return _pigeons_result(r, names, θnames, nr)
end

"""
    _run_ellipsoid_pigeons_fit(snap, data, tepochs, names, lb, ub, nexp, prec, maxeval)

Tempering on an ELLIPSOID's θ.

`tpole_free` IS NOT UNCONDITIONAL HERE, and that is the one place this cannot simply copy the
rapid-rotator runner. Under `:linear` the polar temperature is a pure scale that the flux
normalisation divides out, so `build_ellipsoid_logπ` REFUSES to free it; it carries information
only under `:planck`, where the map's contrast depends on it. The same `tpf` test as
`_run_ellipsoid_hmc_fit`, for the same reason.
"""
function _run_ellipsoid_pigeons_fit(snap, data, tepochs, names, lb, ub, nexp, prec, maxeval;
                                    weights = (1, 1, 1))
    sh = _sh()
    base = star_params(snap)
    tpf = :tpole in names && sh.intensity_model[] === :planck
    θnames = ROTIR.ellipsoid_param_names(; tpole_free = tpf)
    lo, hi = ROTIR.default_ellipsoid_bounds(; tpole_free = tpf)
    lo = collect(Float64, lo); hi = collect(Float64, hi)
    θfull = Float64[Float64(get(snap.params, Symbol(n), 0.0)) for n in θnames]
    free = String[]
    for (k, n) in enumerate(names)
        j = findfirst(==(String(n)), θnames)
        push!(free, θnames[j])
        lo[j] = lb[k]; hi[j] = ub[k]
        θfull[j] = clamp(θfull[j], lb[k], ub[k])
    end
    tess = tessellation_healpix(nexp; T = prec)
    nr = clamp(Int(maxeval), 6, 12)
    r = ROTIR._fit_pigeons(data, tess, tepochs, base; θ0 = θfull, free = free, lb = lo, ub = hi,
                           weights = weights,
                           model = :ellipsoid, tpole_free = tpf,
                           intensity_model = sh.intensity_model[],
                           band = sh.band[] > 0 ? sh.band[] : nothing,
                           n_rounds = nr, PIGEONS_GUI_OPTS...,
                           multithreaded = Threads.nthreads() > 1, verb = true)
    return _pigeons_result(r, names, θnames, nr)
end

"""
    _pigeons_result(r, names, θnames, n_rounds) -> (best, chi2, errs, post, nevals)

The shared tail of the three Pigeons runners: reorder the sampler's columns into the panel's
order, take the median and a half-interquantile width, and print the evidence.

Sampler columns come back in sorted-index order, not the order the panel listed them.
"""
function _pigeons_result(r, names, θnames, nr::Int)
    Printf.@printf("Pigeons: log(Z) = %.4f, %d round trips over 2^%d scans\n",
                   r.logz, r.round_trips, nr)
    # ROUND TRIPS ARE THE DIAGNOSTIC, so say what a low count means rather than leaving a
    # number on the console. A ladder nothing traverses has not visited the other modes, and
    # the posterior is then the single-basin answer a cheaper sampler would have given — which
    # is exactly the result that was once taken at face value on β Cas. With the variational
    # reference this should be tens; a handful means the reference is still a poor fit to the
    # target, and more ROUNDS (not chains) is the lever.
    r.round_trips < 5 && println("  ⚠ only $(r.round_trips) round trip(s): the temperature " *
                                 "ladder was barely traversed, so treat this posterior as " *
                                 "ONE basin, not the whole story. Raise `rounds`.")
    order = sortperm([findfirst(==(String(n)), θnames) for n in names])
    best = zeros(length(names)); errs = Dict{Symbol,Float64}()
    for (slot, k) in enumerate(order)
        best[k] = r.median[slot]
        errs[names[k]] = (r.q84[slot] - r.q16[slot]) / 2
    end
    return best, r.chi2, errs,
           _posterior(r; order = invperm(order),
                      diagnostics = Printf.@sprintf("%d draws, %d chains, %d round trips",
                                                    size(r.samples, 1), r.n_chains,
                                                    r.round_trips)), 0
end

"""
    _run_pigeons_fit(snap, data, tepochs, names, θ0, lb, ub, nexp, prec, maxeval) -> (best, chi2, extra)

Non-reversible parallel tempering over the free parametric parameters, returning the posterior
MEDIAN and half the 16–84 interval, exactly as the NUTS path does.

Chosen over NUTS when the posterior may be MULTIMODAL, which a ROTIR χ² routinely is — NUTS
samples the basin it starts in and reports a confident interval that says nothing about the
others. The round-trip count is what says whether the ladder was actually traversed, and it is
put in the console rather than only in the return value, because a run with no round trips has
not answered the question it was chosen for.

`maxeval` becomes `n_rounds` on a log₂ scale: the run does 2^n_rounds scans and each round
costs as much as everything before it, so the knob has to be a power rather than a count. The
explorer is the slice sampler — no gradient — since ROTIR posteriors are badly conditioned
often enough that a preconditioned MALA is the riskier default (see `_fit_pigeons`).
"""
function _run_pigeons_fit(snap, data, tepochs, names, θ0, lb, ub, nexp, prec, maxeval;
                          weights = (1, 1, 1))
    L = _fit_layout(snap, names)
    getp(k) = Float64(get(snap.params, k, 0.0))
    θfull = Float64[getp(f) for f in L.fields]
    lo = copy(L.lower); hi = copy(L.upper)
    free = String[]
    for (k, n) in enumerate(names)
        j = _theta_slot(L, n)
        push!(free, L.names[j])
        lo[j] = lb[k]; hi[j] = ub[k]
    end
    tess = tessellation_healpix(nexp; T = prec)
    base = star_params(snap)
    # 2^n_rounds scans, and the panel asks for the ROUNDS directly now. Still clamped to
    # 6..12: below six the ladder has not tuned itself and the answer is noise, above twelve
    # one run is longer than a working session.
    nr = clamp(Int(maxeval), 6, 12)
    r = ROTIR._fit_pigeons(data, tess, tepochs, base; θ0 = θfull, free = free, lb = lo, ub = hi,
                           weights = weights,
                           tpole_free = true, n_rounds = nr, PIGEONS_GUI_OPTS...,
                           multithreaded = Threads.nthreads() > 1, verb = true)
    # The rapid rotator's panel names are not its θ names — `PARAMETRIC_THETA` maps between
    # them — so this one reorders on the mapped name rather than the panel's.
    Printf.@printf("Pigeons: log(Z) = %.4f, %d round trips over 2^%d scans\n",
                   r.logz, r.round_trips, nr)
    order = sortperm([_theta_slot(L, n) for n in names])
    best = zeros(length(names)); errs = Dict{Symbol,Float64}()
    for (slot, k) in enumerate(order)
        best[k] = r.median[slot]
        errs[names[k]] = (r.q84[slot] - r.q16[slot]) / 2
    end
    # `r.chi2` is the χ² AT THE MEDIAN, which `_fit_hmc` and `_fit_pigeons` now return. It used
    # to be `NaN` here, so every sampler row in the fit history showed `NaN` in the χ²ᵣ column
    # — in the one table built for comparing a sphere against a rapid rotator, which cannot be
    # done if half the rows have no number. The evidence is the better comparison where there
    # is one, but a sampler without an evidence still has a χ².
    return best, r.chi2, errs,
           _posterior(r; order = invperm(order),
                      diagnostics = Printf.@sprintf("%d draws, %d chains, %d round trips",
                                                    size(r.samples, 1), r.n_chains,
                                                    r.round_trips)), 0
end

"""
    _run_shape_fit(snap, data, tepochs, names, θ0, lb, ub, nexp, prec, maxiter) -> (best, chi2, extra)

The analytic SHAPE fit: `shape_chi2_fg!` under VMLMB.

Hand-coded derivatives of the projected vertices and the face normals — no AD, so this needs no
Zygote — for surface types 0, 1 and 2. The temperature map is held FIXED at the model's own
map: `shape_chi2_fg!` also returns a gradient with respect to the map, and letting both move is
`joint_reconstruct_oi`, which is a reconstruction rather than a parameter fit and belongs on
the Imaging tab.
"""
function _run_shape_fit(snap, data, tepochs, names, θ0, lb, ub, nexp, prec, maxiter;
                       verb::Bool = true, weights = (1, 1, 1))
    layout = SHAPE_THETA[snap.surface_type]
    n = maximum(values(layout))
    inv = Dict(v => k for (k, v) in layout)
    # THE PANEL'S ELEMENT TYPE, not a forced Float64. `shape_chi2_fg!` is generic in a single
    # `T` shared by θ, both gradients, the map and the tessellation — so honouring `prec` is
    # just a matter of building all five at the same type, which is what this does.
    #
    # MEASURED on a sphere's diameter, one epoch of bet Cas at HEALPix 3:
    #
    #                 kernel      per call    vmlmb      result
    #     Float32     5.14 ms     1457 KiB    21 calls   1.027463
    #     Float64     6.37 ms     2854 KiB    12 calls   1.025099
    #
    # Float32 is a fifth faster per evaluation and halves the working set, and it needs more
    # iterations to satisfy the same tolerances on a noisier gradient — so the wall clock comes
    # out a wash and the answer lands 0.23 % away. Set the precision box to Float64 when the
    # diameter itself is the result; Float32 is the right default for everything else, and it
    # is what the rest of the model is evaluated in.
    T  = prec
    θ  = T[T(get(snap.params, inv[j], 0.0)) for j in 1:n]
    lo = fill(T(-Inf), n); hi = fill(T(Inf), n)
    # Only the named parameters move. VMLMB has no frozen-variable facility, so the others are
    # pinned by collapsing their bounds onto their current value — which is exact, and cheaper
    # than the scatter-matrix reduction the Zygote path uses.
    for j in 1:n; lo[j] = θ[j]; hi[j] = θ[j]; end
    for (k, nm) in enumerate(names)
        j = layout[nm]
        lo[j] = T(lb[k]); hi[j] = T(ub[k])
        θ[j]  = clamp(θ[j], T(lb[k]), T(ub[k]))
    end
    tess = tessellation_healpix(nexp; T = T)
    # The base parameters come off the panel as Float64 (`ModelEntry.params` is a
    # `Dict{Symbol,Float64}`), and `merge`ing θ into them would leave a NamedTuple with mixed
    # fields. `convert_params` narrows the floats and leaves `surface_type` and `ldtype` as the
    # `Int`s the geometry branches on.
    base = ROTIR.convert_params(T, star_params(snap))
    star = create_star(tess, base, T(tepochs[1]); secondary = snap.secondary)
    xmap = T.(parametric_temperature_map(base, star; secondary = snap.secondary))
    gθ = zeros(T, n); gx = zeros(T, length(xmap))
    teps = collect(T, tepochs)
    npts = max(_weighted_ndata(data, weights), 1)
    calls = Ref(0)
    function fg!(x, g)
        calls[] += 1
        # `parametric_map`: the temperature map is a FUNCTION of θ here, not a free field, so
        # it is recomputed from θ every call and its derivative enters the gradient. Held
        # fixed — which is what this did — the gradient is 6.6 % wrong overall for an
        # ellipsoid and 13 % wrong in `radius_y` (measured against finite differences of the
        # objective itself), and that mismatch is what stalls the optimiser. A no-op for a
        # sphere, whose uniform map depends on nothing in θ.
        c = shape_chi2_fg!(gθ, gx, xmap, collect(T, x), data, tess, base, teps;
                           parametric_map = true, weights = weights)
        g .= gθ ./ npts
        # `verb` exists for the PRECOMPILE WORKLOAD, which drives this function to compile it
        # and has nowhere for a trace to go: the lines landed in the package build log, twelve
        # of them, between Pkg's own output. During a real fit `start_job!` has redirected
        # stdout into the job's file and the GUI console streams it, which is the point.
        verb && calls[] % 10 == 1 &&
            Printf.@printf("  it %3d   χ²ᵣ = %.6f\n", calls[], c / npts)
        return c / npts
    end
    # ITERATIONS AS ASKED. This divided by 50, which was right while the panel's box said
    # "evaluations" and carried 2000 — and wrong the moment each engine started naming its own
    # unit (`FIT_BUDGETS`): a box reading "iterations 200" ran twenty. `_run_gradient_fit` was
    # converted with the budget table; this path was missed.
    θ̂ = OptimPackNextGen.vmlmb(fg!, θ; lower = lo, upper = hi, maxiter = max(10, maxiter),
                               blmvm = false, verb = false)
    # COMPUTED at `T`, REPORTED at Float64: the fit history, the parameter table and
    # `ModelEntry.params` are all Float64, and widening three numbers at the end costs nothing.
    # The same split `_fit_hmc` makes with its quantiles.
    best = Float64[θ̂[layout[nm]] for nm in names]
    c = shape_chi2_fg!(gθ, gx, xmap, θ̂, data, tess, base, teps; parametric_map = true,
                       weights = weights)
    return best, Float64(c), Dict{Symbol,Float64}(), nothing, 0
end

# One place that knows how each backend is driven, so `shell_fit` stays about the model.
function _run_fit(meth::Symbol, obj, θ0, lb, ub, names, maxeval)
    if meth in (:ultranest, :nautilus)
        # Both return the same fields, so the panel does not care which ran. The error column
        # is HALF the 16-84 interval, i.e. a symmetric 1σ-equivalent: a table cannot show an
        # asymmetric interval in one number, and the quantiles themselves are on `un`.
        #
        # THE BUDGET COMES FROM THE PANEL. `_fit_nested` defaults to `n_eff = 10_000`, which is
        # the setting that took 836 s on a two-parameter sphere fit — so a user who typed 2000
        # into a box labelled "evaluations" got a forty-minute run and no way to shorten it.
        # `n_eff` is the effective sample size, which is what the cost actually scales with;
        # `n_live` follows it, floored so the sampler still has a population to work with.
        # THREADED, which it could not be until `:t3` became re-entrant. `_fit_nautilus`
        # evaluates the likelihood from several tasks at once; `obj` above is
        # `parametric_chi2`, which runs the process default forward kernel `:t3`. That kernel
        # used to keep its mutable scratch in the globally cached `Type3Plan`, so every
        # concurrent evaluation wrote another's grids: the calls threw, `obj`'s own `catch`
        # turned each one into the 1e30 sentinel, and the run came back with `log(Z) = -5e29`,
        # `N_eff = 1` and the FLAT PRIOR as its posterior — after exhausting its budget,
        # because `f_live` never moves when no live point has a finite likelihood. A silently
        # wrong fit, reported in the panel as a successful one. The scratch now comes from
        # `t3_checkout!`, verified bit-identical to the old kernel and clean over 2960
        # concurrent calls that previously corrupted about 1 in 2000.
        un = _fit_nested(meth, obj, String.(names), lb, ub; verb = true,
                         threaded = true,
                         n_eff = clamp(maxeval, 200, 50_000),
                         n_live = clamp(maxeval ÷ 4, 100, 2000),
                         min_num_live_points = clamp(maxeval ÷ 4, 100, 2000))
        return un.median, obj(un.median),
               Dict(n => ((un.q84[k] - un.q16[k]) / 2) for (k, n) in enumerate(names)),
               _posterior(un; diagnostics = Printf.@sprintf("%d samples", size(un.samples, 1))),
               0
    end
    if meth === :ess
        # Gradient-free, affine-invariant, and a SAMPLER rather than a search: it returns the
        # same (median, error, posterior) triple the nested samplers do, so the panel does not
        # care which ran. The objective is a chi2 and already returns a large finite sentinel
        # for a refused parameter set, so the log-posterior is -chi2/2 inside the box and -Inf
        # outside it; a uniform prior on the box is what the schema's bounds already mean.
        nfree = length(θ0)
        logp = function (θ)
            all(lb[i] <= θ[i] <= ub[i] for i in eachindex(θ)) || return -Inf
            c = obj(θ)
            isfinite(c) ? -0.5 * c : -Inf
        end
        # 2d walkers is the scaling requirement, not a preference. They start in a small ball
        # about θ0, clamped inside the box so no walker begins at -Inf (the slice level would
        # then be -Inf and the walker could never move).
        nw = max(4, 2 * nfree + (isodd(2 * nfree) ? 1 : 0))
        span = [min(ub[i] - lb[i], max(abs(θ0[i]), 1.0)) for i in 1:nfree]
        x0 = Matrix{Float64}(undef, nw, nfree)
        for w in 1:nw, i in 1:nfree
            x0[w, i] = clamp(θ0[i] + 0.02 * span[i] * randn(), lb[i] + 1e-9, ub[i] - 1e-9)
        end
        nsweep = clamp(maxeval, 200, 100_000)
        ch, dg = ROTIR.ensemble_slice(logp, x0; nsweep = nsweep, burn = max(200, nsweep ÷ 2),
                                      thin = max(1, nsweep ÷ 2000), verb = false)
        F = ROTIR.flatten_chain(ch)
        med = [median(@view F[:, i]) for i in 1:nfree]
        q16 = [quantile(@view(F[:, i]), 0.16) for i in 1:nfree]
        q84 = [quantile(@view(F[:, i]), 0.84) for i in 1:nfree]
        R = maximum(ROTIR.rhat(ch)); tau = ROTIR.iat(ch)
        # R-hat and tau are REPORTED, not hidden: a sampler that has not converged gives a
        # median and an error bar that look exactly like a converged one, and this is the only
        # place the user can see the difference. tau is a lower bound when it approaches
        # nsweep/50 — see `iat`'s docstring.
        diag = Printf.@sprintf("%d samples, %d walkers; R-hat %.4f%s, tau %.0f sweeps%s",
                               size(F, 1), nw, R, R > 1.01 ? " (NOT converged)" : "",
                               tau, tau > nsweep/50 ? " (lower bound)" : "")
        return med, obj(med),
               Dict(n => ((q84[k] - q16[k]) / 2) for (k, n) in enumerate(names)),
               _posterior((samples = F, q16 = q16, q84 = q84); diagnostics = diag),
               dg.n_eval
    end
    alg = meth === :bobyqa ? :LN_BOBYQA : :LN_NELDERMEAD
    opt = NLopt.Opt(alg, length(θ0))
    NLopt.lower_bounds!(opt, lb)
    # NLopt rejects an infinite upper bound for BOBYQA's initial trust region, and the schema
    # deliberately carries generous rather than infinite ones — but a caller can still widen a
    # bound to Inf by hand, so it is clamped here rather than assumed finite.
    NLopt.upper_bounds!(opt, [isfinite(u) ? u : 1e12 for u in ub])
    NLopt.maxeval!(opt, maxeval)
    NLopt.xtol_rel!(opt, 1e-6)
    NLopt.min_objective!(opt, (θ, g) -> obj(θ))
    chi2, best, _ = NLopt.optimize(opt, copy(θ0))
    # No fourth value: an optimiser has no posterior, and inventing an empty one that the
    # panel then draws as a spike would be worse than saying there is none.
    #
    # The EVALUATION COUNT is what it has instead, and it is worth reporting: it says whether
    # the search converged or simply ran out of budget, which the χ² alone does not. NLopt
    # counts them itself, so this is the number of calls that actually happened rather than
    # the ceiling the panel asked for — a Nelder-Mead that stops on `xtol_rel` reports fewer.
    return best, chi2, Dict{Symbol,Float64}(), nothing, NLopt.numevals(opt)
end

"""
    _posterior(r; diagnostics = "", order = nothing) -> NamedTuple

The sampler-independent part of a result, in the ONE shape [`FitEntry`](@ref) stores.

`order` reorders the columns when the sampler returned the free parameters in its own order
rather than the panel's — `_fit_hmc` and `_fit_pigeons` both sort by
`parametric_free_indices`, which is not how the form lists them, and a posterior plotted
against the wrong parameter name is the kind of wrong that looks right.
"""
function _posterior(r; diagnostics::AbstractString = "", order = nothing)
    S = hasproperty(r, :samples) ? Matrix{Float64}(r.samples) : zeros(0, 0)
    q16 = hasproperty(r, :q16) ? Vector{Float64}(r.q16) : Float64[]
    q84 = hasproperty(r, :q84) ? Vector{Float64}(r.q84) : Float64[]
    if order !== nothing && !isempty(order)
        isempty(S)   || (S = S[:, order])
        isempty(q16) || (q16 = q16[order])
        isempty(q84) || (q84 = q84[order])
    end
    lz  = hasproperty(r, :logz)    ? Float64(r.logz)    : NaN
    lze = hasproperty(r, :logzerr) ? Float64(r.logzerr) : NaN
    return (samples = S, q16 = q16, q84 = q84, logz = lz, logzerr = lze,
            diagnostics = String(diagnostics))
end

# `name\tvalue\terror` — the error column empty for a local search, which does not produce one.
# Empty rather than 0: a zero uncertainty is a claim, and Nelder-Mead makes no such claim.
_fit_table(names, best, errs) =
    join((Printf.@sprintf("%s\t%.6g\t%s", n, best[k],
                          haskey(errs, n) ? Printf.@sprintf("%.3g", errs[n]) : "")
          for (k, n) in enumerate(names)), "\n")

# ── imaging ─────────────────────────────────────────────────────────────────────────────

"""
    REGULARIZER_KINDS

Every regularizer `spheroid_regularization` dispatches on, with what the panel needs to offer
it: a default weight, and the ONE extra scalar (if any) the caller has to choose.

Taken from the dispatch chain itself (src/oichi2_spheroid.jl:995-1019) rather than from the
docs, and it is longer than it looks: `mean`, `bias` and `radialvar` are implemented and were
missing from the first version of this table, which would have made them unreachable from the
GUI while `spheroid_regularization` happily accepted them from a script.

A regularizer entry is FOUR elements — `[name, weight, aux, subset]` — and only the first two
are numbers the panel can supply:

  * `aux` is a precomputed STRUCTURE for most of them: the Sobel operator for `sobel`/`sobel2`,
    the neighbour Laplacian for `tv`/`tv2`, the radial binning for `radflat`/`radialvar`, the
    degenerate direction for `orthold`. Only `bias` takes a plain number there.
  * `subset` is a PIXEL INDEX SET, not a number, and for the radial regularizers it must be
    the exact index set their binning was built from or the call errors.

So the panel chooses names, weights and the one scalar knob, and [`build_regularizers`](@ref)
constructs the rest at run time. Passing four numbers — which is what a naive `name:a:b:c`
string invites — cannot work at all.

Fields: name, default weight, extra-knob label ("" when there is none), extra default, SHORT
label, doc.

TWO descriptions, not one shortened description. The row was asking a single string to name the
maths AND say when to reach for the regularizer, and it overflowed its column doing it — the
column is what is left of the row after a tick, the name, the weight and, on three of the ten
rows, an extra knob and its label, so it is ~90 px narrower exactly where the strings were
longest. Shortening again would have cost the advice, which is the half a user actually needs.

So the SHORT label is the formula alone, sized to fit the narrow case, and the doc carries the
advice into the tooltip, which has as much room as it likes.
"""
const REGULARIZER_KINDS = [
    ("sobel",     1e1, "",      0.0, "∫|∇x| dΩ",
     "∫|∇x| dΩ — edge-preserving"),
    ("sobel2",    1e1, "",      0.0, "∫|∇x|² dΩ",
     "∫|∇x|² dΩ from the 8-neighbour gradient fit — smooth"),
    ("laplace_beltrami", 1e1, "", 0.0, "∫|∇x|² (LB)",
     "∫|∇x|² dΩ, cotangent Laplace–Beltrami — smooth, same meaning at every resolution"),
    ("tv",        1e1, "",      0.0, "‖Lx‖",
     "‖Lx‖ — curvature in L1"),
    ("tv2",       1e1, "",      0.0, "‖Lx‖²",
     "‖Lx‖² — curvature in L2"),
    ("mem",       1e1, "",      0.0, "entropy",
     "maximum entropy, per-pixel"),
    ("mean",      1e1, "",      0.0, "from the mean",
     "departure from the mean"),
    ("bias",      1e1, "B",     2.0, "harmonic bias",
     "harmonic bias; the knob B is the asymmetry"),
    ("radflat",   1e2, "nbins", 6.0, "flat in annuli",
     "flat BETWEEN annuli — single epoch only"),
    ("radialvar", 1e2, "nbins", 6.0, "annulus variance",
     "variance WITHIN annuli"),
    ("orthold",   1e2, "",      0.0, "⊥ LD direction",
     "penalises the limb-darkening-degenerate direction"),
]

"""
    shell_regularizer_kinds() -> String

`name\tweight\textra_label\textra_default\tshort\tdoc` per line. The panel shows the extra
field only where `extra_label` is non-empty, which is what keeps eight of the ten rows to two
controls; it draws `short` in the row and puts `doc` in the tooltip.
"""
shell_regularizer_kinds() =
    join(("$(n)\t$(w)\t$(el)\t$(ed)\t$(sh)\t$(doc)"
          for (n, w, el, ed, sh, doc) in REGULARIZER_KINDS), "\n")

"""
    parse_regularizers(spec) -> Vector{NamedTuple}

`"sobel2:10.0:0;bias:1.0:2.5"` → `[(name="sobel2", weight=10.0, extra=0.0), …]`.

Three fields, not four: the panel can only choose a name, a weight and at most one scalar. The
structures and index sets the engine actually wants are built by
[`build_regularizers`](@ref).

Unknown names and unparseable numbers are dropped rather than raising — the panel builds this
string from its own widgets, so anything malformed is a bug to see in the console rather than a
reason to refuse the run.
"""
function parse_regularizers(spec::AbstractString)
    known = Set(first(r) for r in REGULARIZER_KINDS)
    out = NamedTuple{(:name, :weight, :extra),Tuple{String,Float64,Float64}}[]
    for part in split(spec, ';'; keepempty = false)
        f = split(part, ':')
        length(f) >= 2 || continue
        name = String(strip(f[1]))
        name in known || continue
        w = tryparse(Float64, String(strip(f[2])))
        w === nothing && continue
        e = length(f) >= 3 ? something(tryparse(Float64, String(strip(f[3]))), 0.0) : 0.0
        push!(out, (name = name, weight = w, extra = e))
    end
    return out
end

# The exported script has to build the operators too, and `_literal` of a sparse matrix is not
# something anyone would want in a .jl file. So the log records the CONSTRUCTION rather than the
# result — which is also what a reader wants to see and edit.
function _regularizer_source(specs, n)
    isempty(specs) && return "Any[]"
    # Anything a regularizer has to CONSTRUCT first goes on its own line before the array.
    # Putting it in a trailing comment inside the literal (which the first version did) makes
    # the following comma part of the comment, and the exported script does not parse.
    pre   = String[]
    parts = String[]
    for s in specs
        if s.name in ("sobel", "sobel2")
            push!(parts, "[\"$(s.name)\", $(s.weight), sobel, 1:length(x0)]")
            "sobel = sobel_gradient_healpix($n)" in pre ||
                push!(pre, "sobel = sobel_gradient_healpix($n)")
        elseif s.name == "laplace_beltrami"
            push!(parts, "[\"laplace_beltrami\", $(s.weight), lb, 1:length(x0)]")
            "lb = laplace_beltrami_healpix($n)" in pre ||
                push!(pre, "lb = laplace_beltrami_healpix($n)")
        elseif s.name in ("tv", "tv2")
            push!(parts, "[\"$(s.name)\", $(s.weight), tvinfo, 1:length(x0)]")
            "tvinfo = tv_neighbors_healpix($n)" in pre ||
                push!(pre, "tvinfo = tv_neighbors_healpix($n)")
        elseif s.name in ("mem", "mean")
            push!(parts, "[\"$(s.name)\", $(s.weight), nothing, 1:length(x0)]")
        elseif s.name == "bias"
            push!(parts, "[\"bias\", $(s.weight), $(s.extra), 1:length(x0)]")
        elseif s.name in ("radflat", "radialvar")
            nb = max(2, round(Int, s.extra))
            push!(parts, "[\"$(s.name)\", $(s.weight), bins, bins.idx]")
            any(startswith(q, "bins = ") for q in pre) ||
                push!(pre, "bins = radflat_bins(stars[1]; nbins = $nb)")
        elseif s.name == "orthold"
            push!(parts, "[\"orthold\", $(s.weight), od, stars[1].index_quads_visible]")
            "od = orthold_direction(stars[1], x0, params)" in pre ||
                push!(pre, "od = orthold_direction(stars[1], x0, params)")
        end
    end
    body = "regs  = Any[" * join(parts, ",\n             ") * "]"
    return isempty(pre) ? body : join(pre, "\n") * "\n" * body
end

"""
    build_regularizers(specs, n, star, x0, star_params; subset = nothing) -> Vector{Any}

Turn the panel's choices into the four-element entries `spheroid_regularization` wants.

This is where the operators are actually built, and each is built ONCE per run rather than per
iteration — `sobel_gradient_healpix` and `tv_neighbors_healpix` are sparse assemblies over the
whole mesh, and the criterion is evaluated thousands of times.

Two subsets are in play and they are not the same. The mesh regularizers act on every pixel;
the radial ones act only on `star.index_quads_visible`, and their fourth element must be the
index set their binning was built from — `radflat_bins` checks and errors otherwise, which is
the right behaviour but a confusing one to meet by accident.

`orthold` additionally needs the STARTING map and the limb-darkening parameters, because the
direction it suppresses is the one a change in `ld1` would move the model along. It is skipped
with a warning when that direction is degenerate — an `ld1` that has no effect on this map —
rather than taking the run down.
"""
function build_regularizers(specs, n::Integer, star, x0, star_params; subset = nothing)
    npix = star.npix
    allpix = 1:npix
    out = Any[]
    sobel_op = nothing; tv_op = nothing; bins = nothing
    for s in specs
        if s.name in ("sobel", "sobel2")
            sobel_op === nothing && (sobel_op = sobel_gradient_healpix(Int(n)))
            push!(out, Any[s.name, s.weight, sobel_op, allpix])
        elseif s.name == "laplace_beltrami"
            push!(out, Any[s.name, s.weight, laplace_beltrami_healpix(Int(n)), allpix])
        elseif s.name in ("tv", "tv2")
            tv_op === nothing && (tv_op = tv_neighbors_healpix(Int(n)))
            push!(out, Any[s.name, s.weight, tv_op, allpix])
        elseif s.name in ("mem", "mean")
            # `aux` is unused by these two, but element 4 IS read — `x[reg[4]]` — so a
            # two-element entry raises on the subset lookup, not on the missing aux.
            push!(out, Any[s.name, s.weight, nothing, allpix])
        elseif s.name == "bias"
            push!(out, Any[s.name, s.weight, Float64(s.extra), allpix])
        elseif s.name in ("radflat", "radialvar")
            # These weight their annuli by `star.polyflux`, which `setup_oi!` fills. Said
            # plainly here because the failure otherwise is a BoundsError on an empty vector
            # several calls down, which names neither the regularizer nor the missing step.
            isempty(star.polyflux) &&
                error("$(s.name) needs the Fourier setup: call setup_oi!(data, stars) before " *
                      "building it (the reconstruction path does this for you).")
            nb = max(2, round(Int, s.extra))
            bins === nothing && (bins = radflat_bins(star; nbins = nb, subset = subset))
            push!(out, Any[s.name, s.weight, bins, bins.idx])
        elseif s.name == "orthold"
            od = try
                orthold_direction(star, x0, star_params; subset = subset)
            catch err
                @warn "orthold skipped: $(sprint(showerror, err))"
                nothing
            end
            od === nothing && continue
            push!(out, Any[s.name, s.weight, od,
                           subset === nothing ? star.index_quads_visible : subset])
        end
    end
    return out
end

"""
    _imaging_intensity(sh, m) -> (; model, band, label, err)

How the Imaging tab turns the map into intensity: the model's own "Intensity from", which is
what the Model form shows at its top. A reconstruction honours it rather than carrying a
second, private conversion that could disagree with the form.

`image_reconstruct_oi` takes an `intensity_model` symbol, so the analytic sources map onto it
directly — Rayleigh–Jeans is `:linear` (intensity linear in the map, the historical
behaviour) and Planck is `:planck` at the data's mean wavelength. The model-atmosphere sources
are refused with `err` set: they carry their own μ dependence through an `IntensityProvider`,
which the reconstruction criterion does not take.
"""
function _imaging_intensity(sh::ShellState, m)
    src = round(Int, get(m.params, :intensity_source, INTENSITY_PARAMS[1].default))
    j = findfirst(c -> c.first == src, INTENSITY_PARAMS[1].choices)
    name = j === nothing ? string(src) : INTENSITY_PARAMS[1].choices[j].second
    src == 5 && return (; model = :linear, band = nothing, label = name, err = "")
    if src == 1
        b = _data_band(sh)
        b > 0 || return (; model = :planck, band = nothing, label = name,
                         err = "Planck needs a wavelength: load a dataset first")
        return (; model = :planck, band = b,
                label = Printf.@sprintf("%s at λ = %.4f µm (mean of the data)", name, b * 1e6),
                err = "")
    end
    return (; model = :linear, band = nothing, label = name,
            err = "the Imaging tab reconstructs with Rayleigh–Jeans or Planck; '$(name)' is " *
                  "a model atmosphere carrying its own μ dependence, which " *
                  "image_reconstruct_oi does not take — choose one of the two above")
end

"""
    shell_imaging_intensity() -> String

`code\tlabel\terror` for the Imaging tab's "Intensity from" row: the model's source, how a
reconstruction will apply it, and why it would refuse (empty when it will not). `""` without a
model.
"""
function shell_imaging_intensity()
    sh = _sh()
    m = current_model(sh.session)
    m === nothing && return ""
    src = round(Int, get(m.params, :intensity_source, INTENSITY_PARAMS[1].default))
    ii = _imaging_intensity(sh, m)
    return "$(src)\t$(ii.label)\t$(ii.err)"
end

"""
    shell_reconstruct(nside_exp, regspec, maxiter) -> String

Start a surface reconstruction on a worker.

`regspec` is `name:weight:a:b` entries separated by `;`, exactly as the panel builds it. Parsed
here rather than in QML so the regularizer list has one authority.
"""
function shell_reconstruct(nside_exp, regspec, maxiter)
    sh = _sh()
    d = current_dataset(sh.session)
    d === nothing && return "no dataset"
    m = current_model(sh.session)
    m === nothing && return "no model — the reconstruction needs a geometry to sit on"
    # A BINARY is refused rather than reconstructed as its primary. Everything below builds one
    # `stars` vector from `star_params(m)`, which is the primary's parameters — so with a
    # companion in the model this fitted ONE surface against data containing two, and reported
    # a χ² for it. That is not a worse map, it is a meaningless one, and nothing on screen said
    # so. `binary_reconstruct_oi` is the two-component route and is reachable from the REPL;
    # the panel is not wired to it yet.
    m.companion === nothing ||
        return "the Imaging tab reconstructs one component: untick secondary on the Model " *
               "tab, or use binary_reconstruct_oi from the REPL, which images both"
    ii = _imaging_intensity(sh, m)
    isempty(ii.err) || return ii.err
    imodel, iband = ii.model, ii.band
    specs = parse_regularizers(String(regspec))
    # Level 2 (192 tessels) is coarse for science but it is the level a fit converges on
    # in seconds, which is what makes it the right place to check a setup before committing
    # to one. 6 remains the ceiling.
    n = clamp(Int(nside_exp), 2, 6)
    it = max(1, Int(maxiter))
    p = star_params(m)
    tepochs = copy(d.tepochs)
    data = d.data
    sec = m.secondary
    w = sh.imgweights[]
    _weighted_ndata(data, w) > 0 ||
        return "no data to image: tick at least one observable the dataset has"
    console!(sh, "image_reconstruct_oi(...; nside=2^$n, maxiter=$it, intensity: $(ii.label))" *
                 _weights_note(w); kind = :cmd)
    ikw = (imodel === :linear ? "" :
           Printf.@sprintf(", intensity_model = :%s, band = %.6e", imodel, iband)) *
          (w == (1.0, 1.0, 1.0) ? "" : ", weights = $(_weights_src(w))")
    log!(sh.session,
         "tess  = tessellation_healpix($n)\n" *
         "stars = create_star_multiepochs(tess, $(m.name), tepochs)\n" *
         # `setup_oi!` BEFORE the regularizers, and it was missing entirely: the worker does
         # it, the exported script did not, and `radflat_bins` then indexed an empty
         # `star.polyflux` — a script that parsed, read correctly, and died on the fourth line.
         # Nothing caught it until the round trip was actually RUN in the test suite.
         "setup_oi!(data, stars)\n" *
         "x0    = parametric_temperature_map($(m.name), stars[1])\n" *
         _regularizer_source(specs, n) * "\n" *
         "x = image_reconstruct_oi(x0, data, stars; regularizers = regs, maxiter = $it$(ikw))";
         note = "reconstruct on a HEALPix level-$n sphere", binding = "x")
    # The slot the worker reports into. Made HERE so the job closure can capture it, since the
    # `GuiJob` that holds it does not exist until `start_job!` returns.
    prog = Ref{Any}(nothing)
    return start_job!(sh, :image, function (stop)
        tess  = tessellation_healpix(n)
        stars = create_star_multiepochs(tess, p, tepochs; secondary = sec)
        setup_oi!(data, stars)
        x0 = Float64.(parametric_temperature_map(p, stars[1]; secondary = sec))
        # Built on the worker, because building them needs the star, which is built here. Once
        # per run, not per iteration: these are sparse assemblies over the whole mesh.
        regs = build_regularizers(specs, n, stars[1], x0, p)
        # `copy(x)` — the vector handed to the callback is VMLMB's own working array and is
        # about to be overwritten. `star` travels with the report so the GUI thread can draw
        # without rebuilding a geometry the worker already has.
        cb = (x, nev, f) -> (prog[] = (; n = nev, f = f, x = copy(x), star = stars[1]))
        # (x_start, DATA, stars) — that order, checked against src/oichi2_spheroid.jl:1032.
        x  = image_reconstruct_oi(x0, data, stars; regularizers = regs,
                                  maxiter = it, verbose = true,
                                  intensity_model = imodel, band = iband, weights = w,
                                  callback = cb, callback_every = 25)
        return (; x, n, regs, stars, data, intensity_model = imodel, band = iband, weights = w)
    end; progress = prog)
end

"""
    job_succeeded!(sh, kind, result)

Take up a finished worker's result. **GUI thread only** — this is where the GL work happens.
"""
function job_succeeded!(sh::ShellState, kind::Symbol, result)
    if kind === :image
        # The SAME conversion the reconstruction used, or the χ² reported for the map is the χ²
        # of a different model.
        b = chi2_breakdown(result.x, result.stars, result.data;
                           intensity_model = get(result, :intensity_model, :linear),
                           band = get(result, :band, nothing))
        # THE χ² THE MAP WAS FITTED TO: an observable left out of the reconstruction is left
        # out of the number reported for it. The per-epoch plot still shows all three.
        w = get(result, :weights, (1.0, 1.0, 1.0))
        chi2 = w[1] * b.v2 + w[2] * b.t3amp + w[3] * b.t3phi
        nd = _weighted_ndata(result.data, w)
        e = add_image!(sh.session, "map$(length(sh.session.images) + 1)", result.x,
                       result.n, chi2, nd, result.regs,
                       something(current_dataset(sh.session), (; name = "")).name)
        sh.status = Printf.@sprintf("%s: χ²ᵣ = %.4f over %d points", e.name, chi2 / max(nd, 1),
                                    nd) * _weights_note(w)
        console!(sh, sh.status)
        show_reconstruction!(sh, result.stars[1], result.x; title = e.name)
    elseif kind === :orbit
        # `fit_orbit` returns `params` (the resolved FULL vector), `names`, `chi2_red` and its
        # `spec`, whose `free` holds the indices it moved. A nested sampler adds `q16`/`q84`.
        r = result
        o = sh.orbit
        names = collect(r.names)
        for (k, n) in enumerate(names)
            haskey(o.params, n) && (o.params[n] = r.params[k])
        end
        idx = r.spec.free
        errs = hasproperty(r, :q16) ?
               Dict(names[i] => (r.q84[j] - r.q16[j]) / 2 for (j, i) in enumerate(idx)) :
               Dict{Symbol,Float64}()
        sh.lastorbit = join((Printf.@sprintf("%-12s %12.6g  %s", names[i], r.params[i],
                                             haskey(errs, names[i]) ?
                                             Printf.@sprintf("± %.3g", errs[names[i]]) : "")
                             for i in idx), "\n")
        sh.status = Printf.@sprintf("orbit fit: χ²ᵣ = %.4f  (V² %.2f, T3amp %.2f, T3φ %.2f)",
                                    r.chi2_red, r.chi2_split[1] / max(r.data.nv2, 1),
                                    r.chi2_split[2] / max(r.data.nt3amp, 1),
                                    r.chi2_split[3] / max(r.data.nt3phi, 1))
        console!(sh, sh.status)
        refresh_orbit!(sh)
    elseif kind === :skyexport
        r = result
        sh.status = "wrote $(r.path)  ($(join(r.dims, "×")): nx, nx, channels, epochs)"
        console!(sh, sh.status)
        if r.chk !== nothing
            w = r.chk.worst
            console!(sh, Printf.@sprintf("  OITOOLS vs ROTIR on the data: χ² added %.2e, rel |ΔV| %.2e, max |ΔV²| %.2e, closure phase %.3f° rms",
                                         w.chi2dist, w.rel_dV, w.max_dV2, w.rms_dT3phi_deg))
            w.chi2dist < 0.01 ||
                console!(sh, "  the image adds χ² — use band-limited rendering, a smaller " *
                             "pixel or a larger field")
        end
    elseif kind === :fit
        sh.lastfit = result.table
        sh.status = result.status
        console!(sh, sh.status)
        e = _record_fit!(sh, result)
        # THE MODEL IS NOT TOUCHED. A finished fit used to write its parameters straight into
        # the form, which meant the starting point was gone the moment the fit returned: there
        # was no way to compare the two, to re-run from where you began, or to decline a result
        # that had walked to a bound. The numbers live on the fit — the row in the history and
        # the "Optimized parameters" frame — and `shell_adopt_fit` is what moves them into the
        # model, when and if the user asks.
        current_model(sh.session) === nothing ||
            console!(sh, "$(e.name): the model is unchanged — press Adopt to take these values")
        # The tables still refresh: the history gained a row and the parameter frame follows
        # the selection, even though no parameter moved.
        refresh_model_tab!(sh)
        refresh_data_tab!(sh)
        refresh_posterior!(sh)
    end
    return sh
end

"""
    show_reconstruction!(sh, star, x; title)

Put a reconstructed map into all three Imaging views, and into the Model tab's 3-D scene.

The same three views the Model tab has, with the same decorations and the same colormap: a
reconstruction is a map on a star, and there is no reason to look at it differently from the
map the model generates — which is also the comparison being made most of the time.
"""
function show_reconstruction!(sh::ShellState, star, x; title::AbstractString = "")
    m = current_model(sh.session)
    sp = m === nothing ? nothing : star_params(m)
    sh.lastmap[] = (star = star, x = Vector{Float64}(x), title = String(title))
    # Through `surface_values`, exactly as the Model tab does. This drew `x` raw, so the
    # intensity tick — and with it the limb darkening, which is the whole visible difference
    # between a temperature map and what the interferometer sees — did nothing on this tab.
    # A reconstructed map is a temperature map like any other; there is no reason for the two
    # tabs to disagree about how one is displayed.
    allv = surface_values(sh, x, star; visible_only = false)
    sh.imsky.cbarlabel[] = sh.intensity[] ? "I (arb.)" : "T (K)"
    sh.immoll.cbarlabel[] = whole_surface_label(sh)
    sh.imstar === nothing || (sh.imstar.cbarlabel[] = whole_surface_label(sh))
    show_map!(sh.imsky, star, visible_values(sh, allv, star);
              title = title, star_params = sp, _decor(sh)...)
    show_mollweide!(sh.immoll, allv, star; colorrange = sh.colorrange[])
    sh.imstar === nothing || show_star3d!(sh.imstar, star, allv; colorrange = sh.colorrange[])
    # And the Model tab's own 3-D scene, which is where the recovered map is compared with the
    # parametric one it started from.
    show_star3d!(sh.star, star, allv; colorrange = sh.colorrange[])
    return sh
end

"""
    refresh_image_tab!(sh)

Redraw the Imaging views from the last reconstruction, if there is one.

Unlike the other two tabs there is nothing to recompute: the map came from a run that has
finished, and only the way it is DRAWN can have changed. Without this the tab was frozen at
whatever the run left behind — every view option was live on the Model tab and dead here.
"""
function refresh_image_tab!(sh::ShellState)
    lm = sh.lastmap[]
    if lm === nothing
        msg = "no reconstruction yet — set up the regularizers and press Reconstruct"
        idle!(sh.imsky, msg); idle!(sh.immoll, msg)
        sh.imstar === nothing || idle!(sh.imstar, msg)
        return sh
    end
    show_reconstruction!(sh, lm.star, lm.x; title = lm.title)
    return sh
end

"""
    refresh_both!(sh)

Refresh the Data and Model tabs from ONE geometry build.

Six callbacks refresh both — a parameter edit, an epoch change, a mesh change — and each of the
two used to call `build_epoch_star` for itself, so `create_star` and
`parametric_temperature_map` ran twice for the same parameters at the same epoch. They take a
`got` keyword for exactly this reason; nothing else needs to know.
"""
function refresh_both!(sh::ShellState)
    got = build_epoch_star(sh)
    refresh_data_tab!(sh)
    refresh_model_tab!(sh; got = got)
    return sh
end

"""
    refresh_posterior!(sh)

Redraw the posterior panel from the selected fit. GUI thread only.
"""
function refresh_posterior!(sh::ShellState)
    sh.post === nothing && return sh
    f = current_fit(sh.session)
    if f === nothing
        idle!(sh.post, "no fit yet — run one on this tab")
        return sh
    end
    if isempty(f.samples)
        idle!(sh.post, "$(f.name): $(f.method) is an optimiser and has no posterior")
        return sh
    end
    nfree = size(f.samples, 2)
    i = clamp(sh.postx[], 1, nfree)
    j = clamp(sh.posty[], 1, nfree)
    sh.postx[] = i; sh.posty[] = j
    show_posterior!(sh.post, f.samples, f.names, i, j, f.best[i], f.q16[i], f.q84[i])
    return sh
end

"""
    shell_set_posterior_pair(i, j) -> String

Which parameter the marginal shows and which pair the scatter shows, both 1-based.
"""
function shell_set_posterior_pair(i, j)
    sh = _sh()
    sh.postx[] = max(1, Int(i)); sh.posty[] = max(1, Int(j))
    refresh_posterior!(sh)
    return sh.status
end

shell_posterior_pair() = (sh = _sh(); "$(sh.postx[])\t$(sh.posty[])")

"""
    _record_fit!(sh, result) -> FitEntry

Turn a finished fit into a [`FitEntry`](@ref) the session keeps.

Runs on the GUI thread, from `job_succeeded!`, because that is where the worker's return value
first exists in a place that outlives the job.
"""
function _record_fit!(sh::ShellState, result)
    names = collect(Symbol, result.names)
    best  = collect(Float64, result.best)
    errs  = [Float64(get(result.errs, n, NaN)) for n in names]
    post  = result.post
    q16 = post === nothing || isempty(post.q16) ? fill(NaN, length(names)) : post.q16
    q84 = post === nothing || isempty(post.q84) ? fill(NaN, length(names)) : post.q84
    S   = post === nothing ? zeros(0, 0) : post.samples
    e = FitEntry("$(result.method)_$(length(sh.session.fits) + 1)", result.method,
                 result.model, result.surface_type, names, best, errs, q16, q84, S,
                 post === nothing ? NaN : post.logz,
                 post === nothing ? NaN : post.logzerr,
                 Float64(result.chi2), Int(result.ndata),
                 post === nothing ? "" : post.diagnostics,
                 hasproperty(result, :nevals) ? Int(result.nevals) : 0)
    add_fit!(sh.session, e)
    isempty(S) ||
        console!(sh, "posterior kept: $(size(S, 1)) draws x $(size(S, 2)) parameters" *
                     (isempty(e.diagnostics) ? "" : "  ($(e.diagnostics))"))
    isfinite(e.logz) &&
        console!(sh, Printf.@sprintf("log(Z) = %.4f ± %.4f", e.logz, e.logzerr))
    return e
end

"""
    _work_done(f) -> String

How much work a fit did, for the one column that reports it.

A sampler answers in DRAWS and a local optimiser in criterion EVALUATIONS. They are not the
same quantity and there is no converting between them, but they answer the same question — how
hard did this run work — and a table with a column for each would leave one of them empty on
every row. So the column holds whichever the method produced, and the row already says which
method that was.

The evaluation count is NLopt's own, not the budget the panel asked for: a Nelder-Mead that
stops on `xtol_rel` reports fewer than its ceiling, and that difference is exactly the
difference between converging and running out of budget.
"""
_work_done(f) =
    size(f.samples, 1) > 0 ? string(size(f.samples, 1)) :
    f.nevals > 0           ? string(f.nevals) : "—"

"""
    shell_fits() -> String

`name\tmethod\tmodel\tsurface_type\tchi2r\tlogz\tlogzerr\tnsamples\tdiagnostics` per fit.

THE MODEL COMPARISON TABLE, which is what the evidence is for. Δχ² cannot compare a sphere
with a rapid rotator when χ²ᵣ ≫ 1 — more parameters always fit better — and log(Z) is the
number that can, because the prior volume a model spends is already in it. The samplers were
computing it all along and the GUI was dropping it on the floor.
"""
function shell_fits()
    sh = _sh()
    rows = String[]
    for f in sh.session.fits
        push!(rows, join((f.name, fit_method_short(f.method), f.model,
                          _fit_type_label(f.surface_type),
                          f.ndata > 0 ? Printf.@sprintf("%.4f", f.chi2 / f.ndata) : "—",
                          isfinite(f.logz) ? Printf.@sprintf("%.3f", f.logz) : "—",
                          isfinite(f.logzerr) ? Printf.@sprintf("%.3f", f.logzerr) : "—",
                          _work_done(f), f.diagnostics), "\t"))
    end
    return join(rows, "\n")
end

"""
    _fit_type_label(code) -> String

The surface type's NAME for the fits table, not its integer code.

The column used to print `2`, which is the number `compute_radii` branches on and means
nothing to a reader comparing a sphere against a rapid rotator — which is the entire purpose
of the table. An unknown code still prints as a number rather than raising: a fit recorded by
a newer ROTIR and reloaded here should show something.
"""
_fit_type_label(code::Integer) =
    haskey(SURFACE_TYPES, Int(code)) ? surface_spec(Int(code)).label : string(code)

"Which fit the posterior panel is showing, 1-based; `0` when there is none."
shell_current_fit() = string(_sh().session.current_fit)

"""
    shell_fit_has_posterior() -> "1" | "0"

Whether the SELECTED fit kept draws, so the panel can disable a posterior button rather than
offer a view that would come up empty. An optimiser reports a point and no distribution.
"""
function shell_fit_has_posterior()
    f = current_fit(_sh().session)
    return (f !== nothing && size(f.samples, 1) > 0) ? "1" : "0"
end

"""
    shell_adopt_fit() -> String

Copy the SELECTED fit's parameters into the current model, and report what moved.

A fit already writes its own result back when it finishes, so this is for the OTHER direction
in time: pick an earlier row out of the history and put those numbers back in the form, after
a later fit has overwritten them or after the form has been edited by hand.

Only parameters the current model actually has are written. A fit of a sphere selected while a
rapid rotator is loaded shares `radius` with nothing, and silently inventing fields is how a
model ends up with parameters its `surface_type` never reads.
"""
function shell_adopt_fit()
    sh = _sh()
    f = current_fit(sh.session)
    f === nothing && return "no fit selected"
    m = current_model(sh.session)
    m === nothing && return "no model to adopt into"
    took = String[]; skipped = String[]
    for (k, n) in enumerate(f.names)
        isfinite(f.best[k]) || continue
        # `_fit_put!` resolves a component-qualified label ("2:rpole") onto the right star of a
        # binary, and a bare name onto the primary — the same routing the fit itself uses when
        # it writes its result back, so adopting cannot disagree with fitting.
        if _fit_adopt!(m, n, f.best[k])
            push!(took, String(n))
        else
            push!(skipped, String(n))
        end
    end
    isempty(took) && return "$(f.name): none of its parameters belong to $(m.name)"
    # THE EXPORTED SCRIPT has to end up at the same model. `shell_add_model` logs the
    # parameters as they were when the model was created, and a fit does not move them any
    # more — so adopting is the point at which the script's model has to be rewritten, or the
    # replay would reproduce the starting guess and not the answer.
    log!(sh.session, "$(m.name) = " * _namedtuple_literal(star_params(m));
         note = "adopt $(f.name)", binding = m.name)
    refresh_model_tab!(sh)
    refresh_data_tab!(sh)
    msg = "adopted $(length(took)) parameter$(length(took) == 1 ? "" : "s") from $(f.name)"
    isempty(skipped) || (msg *= "; $(join(skipped, ", ")) not in this model")
    sh.status = msg
    console!(sh, msg)
    return msg
end

function shell_select_fit(i)
    sh = _sh()
    n = Int(i)
    n in eachindex(sh.session.fits) || return sh.status
    sh.session.current_fit = n
    refresh_posterior!(sh)
    # REPORT THE FIT THAT WAS SELECTED. This used to hand back `sh.status`, which is whatever
    # the last thing to run put there — so clicking an older row in the history answered with
    # the NEWEST fit's χ², next to a parameter table that had just switched to the older one.
    f = sh.session.fits[n]
    sh.status = f.ndata > 0 ?
        Printf.@sprintf("%s (%s, %s): χ²ᵣ = %.4f over %d points%s", f.name, f.model,
                        _fit_type_label(f.surface_type), f.chi2 / f.ndata, f.ndata,
                        isfinite(f.logz) ? Printf.@sprintf(", log(Z) = %.3f", f.logz) : "") :
        "$(f.name) ($(f.model))"
    return sh.status
end

"""
    shell_fit_params() -> String

`name\tvalue\terr\tq16\tq84` for the selected fit — the parameter list the posterior panel
plots, and the source of its parameter selectors.
"""
function shell_fit_params()
    f = current_fit(_sh().session)
    f === nothing && return ""
    fmt(v) = isfinite(v) ? Printf.@sprintf("%.6g", v) : "—"
    return join((join((String(f.names[k]), fmt(f.best[k]), fmt(f.err[k]),
                       fmt(f.q16[k]), fmt(f.q84[k])), "\t")
                 for k in eachindex(f.names)), "\n")
end

"""
    shell_imaging_context() -> String

`dataset\tepochs\tmodel\tsurface_type\tmesh` — what a reconstruction started now would
actually run on.

Shown at the top of the Imaging panel because those choices are made on OTHER tabs. A
reconstruction sits on the current model's geometry and fits the current dataset's epochs, and
before this the panel gave no sign of which either was — a run started after switching datasets
looked identical to one started before.
"""
function shell_imaging_context()
    sh = _sh()
    d = current_dataset(sh.session)
    m = current_model(sh.session)
    ds = d === nothing ? "—" : d.name
    ne = d === nothing ? 0 : length(d.data)
    mn = m === nothing ? "—" : m.name
    st = m === nothing ? "—" : string(surface_spec(m.surface_type).name)
    n  = sh.nside_exp[]
    return "$(ds)\t$(ne)\t$(mn)\t$(st)\tHEALPix n=$(n) ($(12 * (2^n)^2) tessels, $(sh.precision[]))"
end

# One row per reconstruction for the Imaging tab's Maps table: name, HEALPix level, tessels,
# χ²ᵣ to four significant figures (a fixed four DECIMALS made 1234.5678 wider than its column
# and ran it into the dataset name), number of data, dataset, and 1 for the current map — the
# one "Save model map" and the sky export use.
function shell_images()
    s = _sh().session
    return join((Printf.@sprintf("%s\t%d\t%d\t%s\t%d\t%s\t%d", e.name, e.nside_exp,
                                 12 * 4^e.nside_exp,
                                 e.ndata > 0 ? Printf.@sprintf("%.4g", e.chi2 / e.ndata) : "—",
                                 e.ndata, e.source, k == s.current_image ? 1 : 0)
                 for (k, e) in enumerate(s.images)), "\n")
end

shell_last_fit() = _sh().lastfit

# ── appearance settings ──────────────────────────────────────────────────────────────────
#
# Saved per USER, not per project: a scale that suits this monitor suits it for every dataset,
# and a settings file inside a repository is one more thing to gitignore. The path follows the
# platform convention, beside OITOOLS' own rather than inside it — the two GUIs have different
# panels and sharing one file would mean each silently dropping the other's keys.

"Where the appearance defaults live."
function gui_settings_file()
    base = Sys.iswindows() ? get(ENV, "APPDATA", homedir()) :
           get(ENV, "XDG_CONFIG_HOME", joinpath(homedir(), ".config"))
    return joinpath(base, "rotir", "gui.toml")
end

"""
    shell_save_settings(payload) -> String

Write the settings panel's values, and return where they went ("" on failure).

`payload` is `key\tvalue` lines, the same tab-separated shape as the rest of this bridge.
Values are stored as they ARRIVE: QML owns what the keys mean, so a key this version does not
understand is round-tripped rather than dropped, and an older ROTIR opening a newer file keeps
the settings it cannot use.
"""
function shell_save_settings(payload)
    sh = SHELL[]
    path = gui_settings_file()
    try
        d = isfile(path) ? TOML.parsefile(path) : Dict{String,Any}()
        for line in split(String(payload), '\n')
            isempty(line) && continue
            f = split(line, '\t')
            length(f) == 2 || continue
            v = tryparse(Float64, f[2])
            d[f[1]] = v === nothing ? String(f[2]) : v
        end
        mkpath(dirname(path))
        open(io -> TOML.print(io, d), path, "w")
        sh === nothing || console!(sh, "saved appearance defaults to " * path)
        return path
    catch err
        sh === nothing || console!(sh, "could not save settings: $(sprint(showerror, err))")
        return ""
    end
end

"""
    shell_settings_path() -> String

Where the appearance config is kept, whether or not anything has been saved there.

The panel shows it on a line of its own, always present and saying what it is: empty until a
button was pressed, it grew the panel the first time anything was saved, and a bare truncated
path beside the buttons reads as a stray fragment rather than as an answer.
"""
shell_settings_path() = gui_settings_file()

"""
    shell_reset_settings() -> String

Delete the saved appearance defaults, and say where they were.

The panel's reset has to remove the FILE and not merely the values on screen: settings are
applied at startup, so a window restored to its built-in look while the file still says
otherwise would come back tweaked at the next launch — which is precisely the confusion the
button exists to end.

Returns the path removed, or `""` when there was nothing saved, which is also the answer when
the settings have never been touched.
"""
function shell_reset_settings()
    sh = SHELL[]
    path = gui_settings_file()
    isfile(path) || return ""
    try
        rm(path)
        sh === nothing || console!(sh, "removed saved appearance defaults: " * path)
        return path
    catch err
        sh === nothing || console!(sh, "could not remove settings: " * sprint(showerror, err);
                                   kind = :err)
        return ""
    end
end

"""
    shell_load_settings() -> String

The saved defaults as `key\tvalue` lines; empty when there are none.

Called from QML before the first plot is drawn. A missing or unreadable file is NOT an error —
the built-in defaults are a perfectly good answer, and a corrupt settings file must not stop the
window from opening.
"""
function shell_load_settings()
    path = gui_settings_file()
    isfile(path) || return ""
    try
        d = TOML.parsefile(path)
        return join(("$(k)\t$(v)" for (k, v) in sort(collect(d), by = first)), "\n")
    catch err
        sh = SHELL[]
        sh === nothing || console!(sh, "ignoring an unreadable settings file at $path")
        return ""
    end
end

"""
    shell_set_plot_scale(x) / shell_set_marker_size(x) -> String

Turn the live plot typography by hand. Zero restores the value computed from the screen.
"""
function shell_set_plot_scale(x)
    v = set_plot_scale!(_qmlreal(x, 0.0))
    sh = SHELL[]
    # `GUI_LIVE`, not just `sh !== nothing`: `applySavedSettings` calls this during
    # `Component.onCompleted`, before the first frame, and redrawing there deadlocks the window
    # black. The value is still applied; the startup path redraws through `refreshAll()` a
    # moment later.
    if sh !== nothing && GUI_LIVE[]
        refresh_both!(sh)
    end
    return v == 0 ? "plot scale: from the screen" : Printf.@sprintf("plot scale: %.2f", v)
end

"""
    shell_plot_scale() -> String

`in_force\tplot_scale_user\tmarker_size_user\tzoom_step`, tab separated.

The panel shows the FIRST in its box and stores the others: `plot_scale_user` may be 0 for
"computed from the screen", and a box reading 0 would say nothing about what the plots are
actually doing. Four fields rather than two because the panel fills the marker and zoom boxes
from the same call, and asking three times would let them disagree.
"""
shell_plot_scale() = join((live_plot_scale(), PLOT_SCALE_USER[], marker_size_user(),
                           zoom_per_detent()), '\t')

"""
    marker_size_user() -> Float64

The data-point size the panel last set, or 0 for "each plot's own default".

Held by OITOOLS rather than here: the Explore tab's plot IS `build_obs_canvas`, which comes from
`OITOOLSGUIExt`, so the number that governs it is that package's `MARKER_SIZE_USER`. A second
copy here would be the one this panel showed and the other the one the plot obeyed.
"""
function marker_size_user()
    O = _oitools_gui()
    O === nothing && return 0.0
    return try
        Float64(O.MARKER_SIZE_USER[])
    catch err
        @debug "OITOOLS has no MARKER_SIZE_USER; the marker box will read auto" err
        0.0
    end
end

"""
    shell_set_marker_size(x) -> String

Set the data-point size in the Explore tab's plot. Zero restores each plot's own default.

Forwarded to OITOOLS, which owns that canvas — see [`marker_size_user`](@ref).
"""
function shell_set_marker_size(x)
    v = _qmlreal(x, 0.0)
    O = _oitools_gui()
    O === nothing && return "no plot to size: the OITOOLS GUI extension is not loaded"
    got = try
        Float64(O.set_marker_size!(v))
    catch err
        return "could not set the marker size: " * sprint(showerror, err)
    end
    # NO REDRAW HERE. This is called from `applySavedSettings` during
    # `Component.onCompleted`, i.e. while Qt's render thread is inside
    # `MakieViewport::setup_buffer` waiting on `begin_julia` for the Julia lock. Redrawing from
    # this thread holds that lock through a full Makie pass, the render thread never gets in,
    # and the window stays BLACK for ever — with no error anywhere, because nothing threw.
    # The panel's own handler redraws after the frame; the startup path is followed by
    # `refreshAll()` regardless.
    return got == 0 ? "marker size: each plot's own" : Printf.@sprintf("marker size: %.0f px", got)
end

# ── export ──────────────────────────────────────────────────────────────────────────────

function shell_export(path)
    sh = _sh()
    p = String(path)
    startswith(p, "file://") && (p = p[8:end])
    p = unique_path(p)
    try
        export_script(sh.session, p)
        sh.status = "wrote $(p)"
    catch err
        sh.status = "could not write $(p): $(sprint(showerror, err))"
    end
    console!(sh, sh.status)
    return sh.status
end

shell_script() = export_script(_sh().session)

"""
    unique_path(path) -> String

`path` if nothing is there, otherwise `name-2.ext`, `name-3.ext`, … until one is free.

EVERY save action in this GUI writes a FIXED name — `rotir_map.fits`, `rotir_orbit.toml`,
`rotirgui_session.jl`, `rotir_sky.png` — because the picker is a file CHOOSER and inventing a
save dialog for a one-click action is a different widget. That is a reasonable trade until the
second save, which silently destroyed the first: reconstruct, look, change a regularizer,
reconstruct, save — and the map you wanted to compare against is gone with no message.

Numbering rather than a timestamp: the files are meant to be found again by eye, and
`rotir_map-2.fits` sorts and reads better than an epoch second.
"""
function unique_path(path::AbstractString)
    p = String(path)
    ispath(p) || return p
    base, ext = splitext(p)
    for k in 2:9999
        cand = "$(base)-$(k)$(ext)"
        ispath(cand) || return cand
    end
    return p
end

# ── sky images for OITOOLS ───────────────────────────────────────────────────────────────
#
# The export dialog's two calls. Everything about WHAT a sky cube is — the layout OITOOLS uses,
# band-limited vs area rendering, the normalisation, the self-check — lives in
# src/sky_export.jl, where it is tested without a window; these only gather the session's
# state into a call and report back.

# What the dialog would export FROM: the model, or one of the Imaging tab's reconstructions.
function _sky_source(sh::ShellState, which::AbstractString)
    m = current_model(sh.session)
    m === nothing && return nothing, "no model — the sky images need a geometry"
    # A BINARY is drawn whole, both components in one image (src/sky_export.jl); only the
    # MODEL can be one, since a reconstruction here is a single surface.
    (m.companion !== nothing && which != "model") &&
        return nothing, "a reconstruction is one component; export the binary model instead"
    ii = _imaging_intensity(sh, m)
    isempty(ii.err) || return nothing, replace(ii.err, "the Imaging tab reconstructs" => "sky export renders")
    prm = star_params(m)
    if which == "model"
        nexp = sh.nside_exp[]
        tess = tessellation_healpix(nexp; T = Float64)
        d = current_dataset(sh.session)
        star = create_star(tess, prm, d === nothing ? 0.0 : d.tepochs[1]; secondary = m.secondary)
        x = Float64.(parametric_temperature_map(prm, star; secondary = m.secondary))
        c = m.companion
        comp2 = nothing
        if c !== nothing
            p2 = star_params(c)
            msgs = validate_star_params(p2)
            isempty(msgs) || return nothing, "secondary: " * first(msgs)
            star2 = create_star(tess, p2, d === nothing ? 0.0 : d.tepochs[1]; secondary = true)
            comp2 = (; params = p2, x = Float64.(parametric_temperature_map(p2, star2; secondary = true)),
                     tess, secondary = true)
        end
        return (; x, tess, prm, ii, label = c === nothing ? "model ($(m.name))" :
                                             "binary model ($(m.name) + secondary)",
                secondary = m.secondary, comp2), ""
    end
    k = findfirst(e -> e.name == which, sh.session.images)
    k === nothing && return nothing, "no reconstruction named $(which)"
    e = sh.session.images[k]
    return (; x = Float64.(e.x), tess = tessellation_healpix(e.nside_exp; T = Float64), prm, ii,
            label = "reconstruction $(e.name)", secondary = m.secondary, comp2 = nothing), ""
end

# The components to draw at these MJDs: the primary at the origin and, for a binary, the
# companion where the model puts it at each epoch (`companion_offsets`: a fixed offset or the
# Orbit tab's orbit).
function _sky_components(sh::ShellState, src, mjd)
    comps = Any[(; params = src.prm, x = src.x, tess = src.tess, secondary = src.secondary)]
    if src.comp2 !== nothing
        offs = companion_offsets(sh, current_model(sh.session), collect(Float64, mjd))
        push!(comps, merge(src.comp2, (; offsets = offs)))
    end
    return comps
end

"""
    shell_sky_export_defaults() -> String

What the export dialog pre-fills, one `key\tvalue` per line: the sources (`model` and every
reconstruction, the current one first), the dataset's epochs (MJD) and channels (µm), the
suggested pixel size (mas), image size and why, the intensity conversion, and `error` when
nothing can be exported.
"""
function shell_sky_export_defaults()
    sh = _sh()
    m = current_model(sh.session)
    d = current_dataset(sh.session)
    out = String[]
    put(k, v) = push!(out, "$(k)\t$(v)")
    srcs = ["model"; [e.name for e in sh.session.images]]
    im = current_image(sh.session)
    im === nothing || (srcs = [im.name; filter(!=(im.name), srcs)])
    put("sources", join(srcs, "|"))
    src, err = _sky_source(sh, first(srcs))
    put("error", err)
    if d !== nothing
        put("epochs", join((Printf.@sprintf("%.4f", t) for t in d.mjd), "|"))
        λs, _ = _data_channels(d)
        put("channels", join((Printf.@sprintf("%.5f", λ * 1e6) for λ in λs), "|"))
    else
        put("epochs", ""); put("channels", "")
    end
    if src !== nothing
        tep = d === nothing ? [0.0] : d.tepochs
        mj  = d === nothing ? [0.0] : d.mjd
        g = sky_grid(_sky_components(sh, src, mj), tep; data = d === nothing ? nothing : d.data)
        put("pixsize", Printf.@sprintf("%.5g", g.pixsize))
        put("nx", string(g.nx))
        put("reason", g.reason * Printf.@sprintf("; %.1f MB per plane", g.mb_per_plane))
        put("intensity", src.ii.label)
    end
    return join(out, "\n")
end

"""
    _data_channels(d) -> (λ, Δλ)

The dataset's CHANNELS, in metres, with their bandwidths — one per spectral channel, not one per
distinct wavelength value. An instrument's channel wavelengths drift from night to night (MIRC-X:
channel 1 is 1.515382 µm on one night and 1.515206 µm on the next), so the raw distinct values
are one per channel PER NIGHT. Values closer than a quarter of the channel width are the same
channel; its plane is rendered at their mean and its `EFF_BAND` is their median bandwidth, which
is the matching tolerance OITOOLS uses (half `EFF_BAND`).
"""
function _data_channels(d)
    lam = Float64.(reduce(vcat, (vec(x.uv_lam) for x in d.data)))
    dlam = Float64.(reduce(vcat, (vec(x.uv_dlam) for x in d.data)))
    o = sortperm(lam); lam = lam[o]; dlam = dlam[o]
    isempty(lam) && return Float64[], Float64[]
    w = filter(>(0), dlam)
    tol = isempty(w) ? 1e-3 * Statistics.median(lam) : 0.25 * Statistics.median(w)
    groups = [[1]]
    for i in 2:length(lam)
        lam[i] - lam[groups[end][1]] > tol ? push!(groups, [i]) : push!(groups[end], i)
    end
    return [Statistics.mean(lam[g]) for g in groups],
           [(b = filter(>(0), dlam[g]); isempty(b) ? 0.0 : Statistics.median(b)) for g in groups]
end

"""
    shell_sky_export(spec) -> String

Render and write a sky cube from the dialog's `key=value` lines:

  * `path` — output FITS file (made unique, `.fits` added);
  * `source` — `model` or a reconstruction's name;
  * `epochs` — `data` (every epoch), `data:1,3` (those, 1-based), `mjd:59000.1,59001.2`, or
    `phases:N` (N rotation phases from the first epoch);
  * `channels` — `data`, `mean`, or `um:1.55,1.60` (micrometres);
  * `pixsize` (mas), `nx`, `rendering` (`bandlimited`|`area`), `normalize` (`plane`|`relative`),
    `check` (`1` to compare OITOOLS's visibilities from the cube with ROTIR's on the data).

Runs on the worker; the result and the check are reported in the console.
"""
function shell_sky_export(spec)
    sh = _sh()
    kv = Dict{String,String}()
    for line in split(String(spec), "\n")
        i = findfirst('=', line); i === nothing && continue
        kv[strip(line[1:i-1])] = strip(line[i+1:end])
    end
    src, err = _sky_source(sh, get(kv, "source", "model"))
    src === nothing && return err
    d = current_dataset(sh.session)
    # EPOCHS: ROTIR epoch times (days from the dataset's first MJD, the session's convention)
    # and the MJD recorded beside them; `dataidx` is which data epochs the check can use.
    em = get(kv, "epochs", "data")
    mjd0 = d === nothing ? 0.0 : minimum(d.mjd)
    tep, mjd, dataidx = if startswith(em, "data")
        d === nothing && return "no dataset: choose MJDs or rotation phases instead"
        idx = em == "data" ? collect(eachindex(d.mjd)) :
              [parse(Int, strip(t)) for t in split(em[6:end], ',') if !isempty(strip(t))]
        all(in(eachindex(d.mjd)), idx) || return "epochs: indices must be 1..$(length(d.mjd))"
        d.tepochs[idx], d.mjd[idx], idx
    elseif startswith(em, "mjd:")
        ms = [parse(Float64, strip(t)) for t in split(em[5:end], ',') if !isempty(strip(t))]
        isempty(ms) && return "epochs: no MJD given"
        ms .- mjd0, ms, Int[]
    elseif startswith(em, "phases:")
        n = parse(Int, strip(em[8:end])); n >= 1 || return "epochs: at least one phase"
        P = Float64(get(src.prm, :rotation_period, 1.0))
        t1 = d === nothing ? 0.0 : d.tepochs[1]
        t = t1 .+ P .* (0:n-1) ./ n
        collect(t), collect(t .+ mjd0), Int[]
    else
        return "epochs: unknown choice $(em)"
    end
    cm = get(kv, "channels", "mean")
    bands = nothing
    λs = if cm == "data"
        d === nothing && return "no dataset: give the wavelengths"
        λc, bands = _data_channels(d)
        λc
    elseif cm == "mean"
        b = _data_band(sh); b > 0 || return "no dataset: give the wavelengths"
        [b]
    elseif startswith(cm, "um:")
        [1e-6 * parse(Float64, strip(t)) for t in split(cm[4:end], ',') if !isempty(strip(t))]
    else
        return "channels: unknown choice $(cm)"
    end
    isempty(λs) && return "channels: none"
    pix = tryparse(Float64, get(kv, "pixsize", "")); nx = tryparse(Int, get(kv, "nx", ""))
    (pix === nothing || pix <= 0) && return "pixel size must be a positive number (mas)"
    (nx === nothing || nx < 8) && return "image size must be an integer ≥ 8"
    isodd(nx) && (nx += 1)
    rendering = Symbol(get(kv, "rendering", "bandlimited"))
    normalize = Symbol(get(kv, "normalize", "plane"))
    check = get(kv, "check", "1") == "1" && !isempty(dataidx)
    p = String(get(kv, "path", "rotir_sky.fits")); startswith(p, "file://") && (p = p[8:end])
    endswith(lowercase(p), ".fits") || (p *= ".fits")
    p = unique_path(p)
    imodel = src.ii.model
    comps = _sky_components(sh, src, mjd)
    # The size is the caller's to choose, but say so before minutes of rendering: a wide binary
    # at Nyquist/4 is thousands of pixels on a side.
    mb = nx^2 * 4 * length(λs) * length(tep) / 2^20
    mb > 4096 && return Printf.@sprintf("that cube would be %.1f GB — use fewer epochs/channels, a coarser pixel or a smaller field", mb / 1024)
    console!(sh, Printf.@sprintf("sky export: %s, %d epoch(s) × %d channel(s), %d×%d px at %.4g mas, %s, %s",
                                 src.label, length(tep), length(λs), nx, nx, pix, rendering,
                                 src.ii.label); kind = :cmd)
    log!(sh.session,
         "cube, meta = render_sky_cube(components; tepochs = $(repr(tep)), " *
         "wavelengths = $(repr(λs)), pixsize = $(pix), nx = $(nx), rendering = :$(rendering), " *
         "intensity_model = :$(imodel), normalize = :$(normalize), mjd = $(repr(mjd)))\n" *
         "write_sky_cube($(repr(p)), cube, meta)";
         note = "sky images for OITOOLS", binding = "cube")
    data = d === nothing ? nothing : d.data
    return start_job!(sh, :skyexport, function (stop)
        cube, meta = render_sky_cube(comps; tepochs = tep, wavelengths = λs,
                                     pixsize = pix, nx = nx, rendering = rendering,
                                     intensity_model = imodel, normalize = normalize, mjd = mjd,
                                     bandwidths = bands)
        write_sky_cube(p, cube, meta; comment = src.label)
        chk = check ? check_sky_cube(cube, meta, comps, data[dataidx]) : nothing
        return (; path = p, dims = size(cube), chk)
    end)
end

# ── the surface map as a file ───────────────────────────────────────────────────────────
#
# The command log exports the SCRIPT that would rebuild a session. This exports the RESULT: the
# per-tessel map, with the tessellation level and every parameter beside it, so the χ² can be
# recomputed from the file alone. The two answer different questions — "how did you get this"
# and "what is this" — and a reconstruction that took twenty minutes should not have to be
# re-run to be looked at again.

"""
    shell_save_map(path) -> String

Write the current map to FITS.

WHICH map: the reconstructed image if there is one, otherwise the current model's parametric
temperature map. A reconstruction is what a user means by "the map" once they have one, and
until then the model's own map is the only thing on screen.

Everything needed to reproduce the χ² travels with it — the HEALPix level, every resolved
parameter (ties applied, since a tie is a way of writing a value and not a value of its own),
and the epoch times of the dataset it was fitted against.
"""
function shell_save_map(path)
    sh = _sh()
    p = String(path); startswith(p, "file://") && (p = p[8:end])
    endswith(lowercase(p), ".fits") || (p *= ".fits")
    p = unique_path(p)
    m = current_model(sh.session)
    im = current_image(sh.session)
    if m === nothing && im === nothing
        sh.status = "nothing to save — add a model or run a reconstruction first"
        console!(sh, sh.status); return sh.status
    end
    d = current_dataset(sh.session)
    try
        # The image carries the level it was reconstructed at, which is not necessarily the
        # one the Model tab is showing now; saving the panel's level against the image's map
        # is how a file becomes unloadable.
        nexp = im === nothing ? sh.nside_exp[] : im.nside_exp
        prm  = m === nothing ? NamedTuple() : star_params(m)
        x = if im !== nothing
            im.x
        else
            star = create_star(tessellation_healpix(nexp; T = sh.precision[]), prm,
                               d === nothing ? 0.0 : d.tepochs[1]; secondary = m.secondary)
            Float64.(parametric_temperature_map(prm, star; secondary = m.secondary))
        end
        save_surface_map(p, x, prm; nside_exp = nexp,
                         tepochs = d === nothing ? nothing : d.tepochs,
                         mjd     = d === nothing ? nothing : d.mjd,
                         secondary = m === nothing ? false : m.secondary,
                         field = :temperature,
                         chi2 = im === nothing ? NaN : im.chi2,
                         ndata = im === nothing ? 0 : im.ndata,
                         comment = im === nothing ? "parametric map" : im.name)
        log!(sh.session, "save_surface_map($(repr(p)), x, params; nside_exp = $(nexp))";
             note = "save map")
        sh.status = "wrote $(p) — $(length(x)) tessels at level $(nexp)"
    catch err
        sh.status = "could not write $(p): $(sprint(showerror, err))"
    end
    console!(sh, sh.status)
    return sh.status
end

"""
    shell_load_map(path) -> String

Read a map back, as a model and an image.

The parameters become a model entry, so the form shows what produced the map and a fit can be
started from it; the values become an image entry, so the views draw the map itself rather
than the parametric map its parameters would generate. Both, because a reconstruction is
exactly the case where the two differ.
"""
function shell_load_map(path)
    sh = _sh()
    p = String(path); startswith(p, "file://") && (p = p[8:end])
    local mp
    try
        mp = load_surface_map(p)
    catch err
        sh.status = "could not read $(p): $(sprint(showerror, err))"
        console!(sh, sh.status); return sh.status
    end
    mp.tessellation === :healpix ||
        return (sh.status = "$(basename(p)) is on a $(mp.tessellation) tessellation, " *
                            "which the GUI cannot show yet"; console!(sh, sh.status); sh.status)
    if haskey(mp.params, :surface_type)
        st = Int(mp.params.surface_type)
        # ONE model at a time, the same rule "+ model" follows: a second model left behind
        # would silently decide what the χ² column and the next reconstruction were about.
        empty!(sh.session.models); sh.session.current_model = 0; sh.chi2key[] = nothing
        m = add_model!(sh.session, st; name = "loaded_" * basename(p), secondary = mp.secondary)
        for (k, v) in pairs(mp.params)
            k === :surface_type && continue
            haskey(m.params, k) && (m.params[k] = Float64(v))
        end
    end
    sh.nside_exp[] = mp.nside_exp
    add_image!(sh.session, basename(p), mp.x, mp.nside_exp, mp.chi2, mp.ndata, Any[], p)
    log!(sh.session, "m = load_surface_map($(repr(p)))"; note = "load map", binding = "m")
    st = "loaded $(basename(p)) — $(length(mp.x)) tessels at level $(mp.nside_exp)" *
         (isfinite(mp.chi2) ? Printf.@sprintf(", χ² = %.4g", mp.chi2) : "")
    console!(sh, st)
    refresh_both!(sh)
    sh.status = st
    return sh.status
end

# ── the model GEOMETRY as a file ────────────────────────────────────────────────────────
#
# The third of three things this panel can write, and the three are deliberately different:
#
#   * the COMMAND LOG is the script that would rebuild the session — "how did you get this";
#   * the MAP is what the surface is emitting — "what is this";
#   * the GEOMETRY is the surface itself — "what shape was it".
#
# A map file already carries the parameters, so in principle it could rebuild the mesh. In
# practice it would rebuild it with whatever `create_star` does the day it is read: the Roche
# equipotential is a root solve, the rapid rotator's radii come from a cubic, the visibility
# clip is a sigmoid with a κ. Writing the mesh is what turns "the geometry comes out the same"
# from an assumption into something that can be checked — and `shell_load_geometry` checks it.

"""
    shell_save_geometry(path) -> String

Write the current model's GEOMETRY to FITS: the mesh, the parameters that built it, the
surface type, the tessellation, and — for a binary — the companion and how it is placed.

NOT the map values. `shell_save_map` has those, and the two are separate files because they
have separate lifetimes: one geometry carries every map ever fitted on it.
"""
function shell_save_geometry(path)
    sh = _sh()
    p = String(path); startswith(p, "file://") && (p = p[8:end])
    endswith(lowercase(p), ".fits") || (p *= ".fits")
    p = unique_path(p)
    m = current_model(sh.session)
    if m === nothing
        sh.status = "no model — add one before saving a geometry"
        console!(sh, sh.status); return sh.status
    end
    got = build_epoch_star(sh)
    if got === nothing
        sh.status = "the model does not build — fix the parameters first"
        console!(sh, sh.status); return sh.status
    end
    d = current_dataset(sh.session)
    c = m.companion
    try
        # The star EXACTLY as it is drawn, at the epoch on screen — `build_epoch_star` is what
        # the canvases are fed, so the file and the picture cannot disagree.
        star2 = length(got) >= 5 ? got[3] : nothing
        save_star_geometry(p, got[1], star_params(m);
                           nside_exp = sh.nside_exp[],
                           tessellation = sh.tessel[],
                           star2 = star2,
                           params2 = star2 === nothing ? nothing : star_params(c),
                           place = c === nothing ? :offset : c.place,
                           offset = c === nothing ? (0.0, 0.0, 0.0) : c.offset,
                           tepochs = d === nothing ? nothing : d.tepochs,
                           mjd     = d === nothing ? nothing : d.mjd,
                           secondary = m.secondary,
                           comment = m.name)
        log!(sh.session, "save_star_geometry($(repr(p)), star, params; " *
                         "nside_exp = $(sh.nside_exp[]))"; note = "save geometry")
        sh.status = "wrote $(p) — $(got[1].npix) tessels" *
                    (star2 === nothing ? "" : " + $(star2.npix) for the secondary")
    catch err
        sh.status = "could not write $(p): $(sprint(showerror, err))"
    end
    console!(sh, sh.status)
    return sh.status
end

"""
    shell_load_geometry(path) -> String

Read a geometry back as the current model, and say how closely it rebuilds.

The model is rebuilt from the stored PARAMETERS — that is what the rest of the panel edits and
what every χ² is computed from — and the stored MESH is then compared against the rebuild. The
largest vertex disagreement goes to the console in mas.

Zero for a surface that is a multiplication — a sphere, an ellipsoid. Float32 round-off, of
order 1e-7 of the radius, for one whose radii come from a root solve: a Roche component
rebuilds to 6e-8 mas on a 0.7 mas polar radius, which is the mesh precision and not a
disagreement. Anything larger than that names, in one number, a geometry the current code no
longer reproduces — which is exactly what a parameter-only file cannot tell you.
"""
function shell_load_geometry(path)
    sh = _sh()
    p = String(path); startswith(p, "file://") && (p = p[8:end])
    local g
    try
        g = load_star_geometry(p)
    catch err
        sh.status = "could not read $(p): $(sprint(showerror, err))"
        console!(sh, sh.status); return sh.status
    end
    g.tessellation === :healpix ||
        return (sh.status = "$(basename(p)) is on a $(g.tessellation) tessellation, " *
                            "which the GUI cannot show yet"; console!(sh, sh.status); sh.status)
    haskey(g.params, :surface_type) ||
        return (sh.status = "$(basename(p)) carries no parameters — nothing to rebuild from";
                console!(sh, sh.status); sh.status)

    # ONE model at a time, the same rule "+ model" and "Load map" follow.
    empty!(sh.session.models); sh.session.current_model = 0; sh.chi2key[] = nothing
    m = add_model!(sh.session, Int(g.params.surface_type);
                   name = "loaded_" * splitext(basename(p))[1], secondary = g.secondary)
    for (k, v) in pairs(g.params)
        k === :surface_type && continue
        haskey(m.params, k) && (m.params[k] = Float64(v))
    end
    if g.params2 !== nothing && haskey(g.params2, :surface_type)
        shell_set_binary("1", Int(g.params2.surface_type))
        c = m.companion
        for (k, v) in pairs(g.params2)
            k === :surface_type && continue
            haskey(c.params, k) && (c.params[k] = Float64(v))
        end
        c.place = g.place
        c.offset = g.offset
    end
    sh.nside_exp[] = g.nside_exp
    sh.precision[] = eltype(g.star.vertices_xyz)

    # THE CHECK the file exists for. Same epoch, same level, same precision — so any
    # difference is the geometry code, not the inputs.
    dev = try
        r = create_star(tessellation_healpix(g.nside_exp; T = sh.precision[]),
                        star_params(m), sh.precision[](g.t); secondary = m.secondary)
        size(r.vertices_xyz) == size(g.star.vertices_xyz) ?
            maximum(abs.(Float64.(r.vertices_xyz) .- Float64.(g.star.vertices_xyz))) : NaN
    catch err
        console!(sh, "could not rebuild the geometry to compare: $(sprint(showerror, err))")
        NaN
    end
    # `@sprintf` takes a LITERAL format string — a concatenation is an ArgumentError at
    # macro-expansion time, i.e. when the extension is precompiled, not when this runs.
    console!(sh, isnan(dev) ? "geometry: the rebuild could not be compared" :
                 dev == 0   ? "geometry: the rebuild is identical to the stored mesh" :
                 Printf.@sprintf("geometry: the rebuild differs from the stored mesh by at most %.3g mas", dev))
    log!(sh.session, "g = load_star_geometry($(repr(p)))"; note = "load geometry", binding = "g")
    refresh_both!(sh)
    sh.status = "loaded $(basename(p)) — $(g.star.npix) tessels at level $(g.nside_exp)" *
                (g.star2 === nothing ? "" : ", binary")
    console!(sh, sh.status)
    return sh.status
end

# ── Qt conflict check ───────────────────────────────────────────────────────────────────

"""
    check_qt_conflict()

Warn if matplotlib has already mapped a Qt into this process.

ROTIR keeps matplotlib behind `ROTIRPythonPlotExt` precisely so that `using ROTIR` does not do
this, but a session that loaded PythonPlot for its own reasons will have. The symptom is not a
clean error — it is a loader failure deep in Qt:

    libQt6DBus.so: undefined symbol: _ZN14QObjectPrivateC2E16QtPrivate_6_10_2

which is what a second, incompatible Qt in the same address space looks like. Naming it here
turns twenty minutes of confusion into one line.
"""
function check_qt_conflict()
    loaded = String[]
    for (_, m) in Base.loaded_modules
        n = String(nameof(m))
        n in ("PythonPlot", "PyPlot") && push!(loaded, n)
    end
    isempty(loaded) && return nothing
    @warn """
        $(join(loaded, " and ")) is loaded in this session.

        matplotlib probes for an interactive backend and can map a second Qt into this
        process, which QML cannot survive — the failure appears as an undefined Qt symbol
        rather than as anything mentioning matplotlib.

        Start the GUI in a session that has not imported it, or set MPLBACKEND=Agg before
        the import.
        """
    return nothing
end
