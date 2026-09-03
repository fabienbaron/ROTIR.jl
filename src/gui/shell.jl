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
    ShellState

One mutable bundle holding the session, the canvases and the console, reached through the
global `SHELL[]` because Qt callbacks take no arguments.
"""
mutable struct ShellState
    session::Session
    sky::Any            # SkyCanvas
    star::Any           # StarCanvas
    moll::Any           # MollCanvas
    chi2::Any           # Chi2Canvas
    # The Imaging tab's own pair. Separate figures, not a second view of the two above: a Makie
    # Figure belongs to one MakieArea at a time, and a reconstruction must not overwrite the
    # plot being read on another tab.
    imsky::Any          # SkyCanvas
    immoll::Any         # MollCanvas
    msky::Any           # SkyCanvas — the Model tab's own, and its DEFAULT view
    # The Data tab's observable plot: an OITOOLS `LiveCanvas`, reached through its GUI
    # extension rather than reimplemented. See `build_obs_canvas`.
    obs::Any
    obsmodel::Any       # Observable{Vector{Point2f}} — the model prediction drawn over it
    obskind::Base.RefValue{Symbol}
    obscolor::Base.RefValue{Symbol}
    obsoverlay::Base.RefValue{Bool}
    # What the surface views COLOUR BY: the temperature map itself, or the emergent intensity —
    # limb darkening multiplied in, and optionally a Planck conversion at the observing
    # wavelength. Two different physical quantities, and the one you want depends on whether
    # you are reading the model or reading what the interferometer sees.
    intensity::Base.RefValue{Bool}
    intensity_model::Base.RefValue{Symbol}
    band::Base.RefValue{Float64}          # metres; 0 means "take it from the data"
    # Which decorations the sky views draw. One set for every panel, because they are the same
    # picture seen on different tabs and three different sets of annotations for one star would
    # be three things to keep in your head rather than one.
    decor::Dict{Symbol,Bool}
    # Graticule settings. Spacing in DEGREES, which is how anyone reading a map thinks about
    # it; `graticule_segments` counts lines, so the two are converted where they meet.
    gratlat::Base.RefValue{Float64}
    gratlon::Base.RefValue{Float64}
    gratcolor::Base.RefValue{String}
    # The mesh the model is evaluated on, and the float type it is built in. Both belong to the
    # MODEL rather than to a plot: they decide what is being fitted, not how it is drawn.
    tessel::Base.RefValue{Symbol}         # :healpix (:longlat is not wired yet)
    nside_exp::Base.RefValue{Int}         # log2(nside); 3 for a result, 2 to try something
    precision::Base.RefValue{DataType}
    # The Orbit perspective. One orbit at a time, like the model — a session fitting two
    # unrelated binaries at once is not a thing anyone does, and a selector for it would be
    # another way to have the wrong one selected.
    orbit::Any                            # OrbitEntry
    orbitcanvas::Any                      # OrbitCanvas
    lastorbit::String                     # the fit result table
    console::Vector{String}
    status::String
    job::Any            # GuiJob or nothing
    jobkind::Symbol
    lastfit::String     # rendered result table of the last fit
    imaging_defaults::Dict{Symbol,Any}
    # Cache for `epoch_chi2`, keyed on what it actually depends on. Recomputing it is seconds
    # of work — the geometry at every epoch plus `setup_oi!`, i.e. the polygon FT against every
    # uv point — and before this it ran on every tab switch.
    chi2cache::Base.RefValue{Any}
    chi2key::Base.RefValue{Any}
    # The posterior panel: the canvas, and which two parameters it is showing. Appended rather
    # than placed with the other canvases because every positional construction of this struct
    # — the tests, the precompile workload — would otherwise have to be renumbered.
    post::Any                             # PostCanvas
    postx::Base.RefValue{Int}             # marginal / pair x
    posty::Base.RefValue{Int}             # pair y
    # Whether the observable plot's y axis is logarithmic. Appended for the same reason as
    # `post`: the field order is what every positional construction of this struct spells out.
    obslog::Base.RefValue{Bool}
    # The Imaging tab's own 3-D scene. A reconstruction is a map like any other and is looked
    # at the same three ways; it needs a scene of its own because a Figure belongs to one
    # MakieArea at a time — pointing two tabs at one figure leaves one blank and draws the
    # other detached across the panel beside it.
    imstar::Any                           # StarCanvas
end

const SHELL = Ref{Any}(nothing)

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
    append = String(add) == "1" && current_dataset(sh.session) !== nothing
    try
        e = load_dataset!(sh.session, ps; add = append, split = String(dosplit) == "1")
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
    append = String(add) == "1" && current_dataset(sh.session) !== nothing
    try
        e = load_dataset!(sh.session, p; add = append, split = String(dosplit) == "1")
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
        "data = split_epochs(readoifits($(_literal(e.paths[1])))[1,1])[1]" :
        "files = " * _literal(e.paths) * "\n" *
        "data  = collect(readoifits_multiepochs(files; warn = false, verbose = false)[1, :])"
    code *= "\ntepochs = [d.mean_mjd for d in data]; tepochs .-= minimum(tepochs)"
    log!(sh.session, code; note = "dataset $(e.name)", binding = "data")
    return e
end

"""
    shell_select_dataset(i) -> String

Make dataset `i` current. One-based, and out of range is a no-op rather than an error: the
index comes from a list view whose model can be a frame behind.
"""

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
cannot offer a type `compute_radii` does not branch on.
"""
shell_surface_types() = join(("$(c)\t$(SURFACE_TYPES[c].name)\t$(SURFACE_TYPES[c].label)\t" *
                              replace(SURFACE_TYPES[c].doc, "\n" => " ")
                              for c in SURFACE_TYPE_ORDER), "\n")

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
    m = add_model!(sh.session, Int(surface_type))
    log!(sh.session, "$(m.name) = " * _namedtuple_literal(star_params(m));
         note = "model $(m.name)", binding = m.name)
    console!(sh, "model $(m.name): surface_type $(m.surface_type) at schema defaults")
    sh.status = "model $(m.name)"
    refresh_model_tab!(sh)
    return sh.status
end

shell_models() = join(("$(m.name)\t$(m.surface_type)\t$(length(m.free))"
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
"""
function shell_params()
    sh = _sh()
    m = current_model(sh.session)
    m === nothing && return ""
    used = ld_coefficients_used(round(Int, get(m.params, :ldtype, 3.0)))
    rows = String[]
    for ps in surface_params(m.surface_type)
        v  = get(m.params, ps.name, ps.default)
        lo, hi = get(m.bounds, ps.name, (ps.lo, ps.hi))
        tie = get(m.ties, ps.name, "")
        inert = ps.group === :limbdark && ps.name in (:ld1, :ld2, :ld3, :ld4) &&
                !(ps.name in used)
        state = !isempty(tie) ? "tied" :
                (inert ? "fixed" : (ps.name in m.free ? "free" : "fixed"))
        ch = join(("$(k)=$(vv)" for (k, vv) in ps.choices), "|")
        push!(rows, join((ps.name, ps.label, ps.unit, repr(v), state, repr(lo), repr(hi),
                          tie, ps.group, ps.kind, ch,
                          replace(ps.doc, "\n" => " ", "\t" => " ")), "\t"))
    end
    return join(rows, "\n")
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
    v = tryparse(Float64, String(value))
    v === nothing && return "not a number: $(value)"
    m.params[Symbol(String(name))] = v
    refresh_model_tab!(sh)
    return ""
end

"""
    shell_set_param_state(name, state) -> String

`free`, `fixed`, or `tied` (with the expression set separately by `shell_set_tie`).
"""
function shell_set_param_state(name, state)
    sh = _sh()
    m = current_model(sh.session)
    m === nothing && return "no model"
    n = Symbol(String(name)); s = String(state)
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
    l = tryparse(Float64, String(lo)); h = tryparse(Float64, String(hi))
    (l === nothing || h === nothing) && return "bounds must be numbers"
    l < h || return "lower bound must be below upper"
    m.bounds[Symbol(String(name))] = (l, h)
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
    n = Symbol(String(name)); e = String(expr)
    m.ties[n] = e
    delete!(m.free, n)
    v = eval_tie(e, m.params)
    v === nothing && return isempty(strip(e)) ? "" : "does not evaluate yet"
    m.params[n] = v
    refresh_model_tab!(sh)
    return Printf.@sprintf("= %.6g", v)
end

"""
    shell_validate_model() -> String

Every problem the schema can see with the current model, one per line. Empty means it builds.
"""
function shell_validate_model()
    sh = _sh()
    m = current_model(sh.session)
    m === nothing && return "no model"
    return join(validate_star_params(star_params(m)), "\n")
end

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
regularisers and the shape gradients are all written against the HEALPix neighbour structure,
so offering it would mean a mesh half the panel does not work on.

`precision` is `Float32` by default. It is what `tessellation_healpix` itself defaults to, it
halves the memory the polygon FT needs, and the visibilities are not measured to seven digits
— but the analytic shape fit needs Float64, and switches to it on its own rather than making
that the default for everything.
"""
function shell_set_tessellation(kind, nside_exp, precision)
    sh = _sh()
    k = Symbol(String(kind))
    k === :healpix || return "only :healpix is wired; :longlat has no regularisers or gradients"
    sh.tessel[] = k
    sh.nside_exp[] = clamp(Int(nside_exp), 2, 6)
    p = String(precision)
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
    "nufft\tQuadrature + NUFFT (fastest)\t" *
    "Gauss-Legendre over each tessel folded into one type-3 NUFFT; 2.6 ms at HEALPix 3 and " *
    "4.1 ms at 5, where the exact kernels take 7 ms and 86 ms — the cost barely moves with " *
    "the mesh. Accurate to 6.8e-7 at HEALPix 3 and 2.5e-9 above it",
    "turbo\tExact, vectorised\t" *
    "the closed-form polygon transform with SIMD transcendentals — no approximation " *
    "parameters at all, 17x the reference at HEALPix 3",
    "scalar\tExact, reference\t" *
    "the same arithmetic in plain Julia; the definition the other two are tested against, " *
    "and what to fall back on if a fit looks wrong"), "\n")

"Which forward kernel is selected."
shell_polyft_backend() = String(ROTIR.POLYFT_BACKEND[])

"""
    shell_set_polyft_backend(kind) -> String

Choose the forward kernel. Invalidates the χ² cache, since the number it holds came from the
other one.
"""
function shell_set_polyft_backend(kind)
    sh = _sh()
    k = Symbol(String(kind))
    k in (:nufft, :turbo, :scalar) || return "backend must be nufft, turbo or scalar"
    ROTIR.POLYFT_BACKEND[] = k
    sh.chi2key[] = nothing
    console!(sh, "polyft backend: $(k)")
    refresh_both!(sh)
    return "polyft backend: $(k)"
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
    la = something(tryparse(Float64, String(lat_deg)), sh.gratlat[])
    lo = something(tryparse(Float64, String(lon_deg)), sh.gratlon[])
    sh.gratlat[] = clamp(la, 5.0, 90.0)
    sh.gratlon[] = clamp(lo, 5.0, 90.0)
    c = String(colour)
    c in GRATICULE_COLORS && (sh.gratcolor[] = c)
    refresh_both!(sh)
    return ""
end

"""
    shell_set_decoration(name, on) -> String

Turn one decoration on or off, everywhere.
"""
function shell_set_decoration(name, on)
    sh = _sh()
    k = Symbol(String(name))
    haskey(sh.decor, k) || return "unknown decoration $(name)"
    sh.decor[k] = String(on) == "1"
    refresh_both!(sh)
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
    sh.intensity[] = String(intensity) == "1"
    m = Symbol(String(model))
    m in (:linear, :planck) || return "unknown intensity model $(model)"
    sh.intensity_model[] = m
    b = something(tryparse(Float64, String(band_um)), 0.0)
    sh.band[] = b > 0 ? b * 1e-6 : _data_band(sh)
    # Fall back to the data's own mean wavelength whenever the field is empty, zero or
    # unparseable — the band the observations were taken in is the only defensible default for
    # a conversion whose whole point is that it depends on the band.
    if m === :planck && sh.band[] <= 0
        return "Planck needs a wavelength: load a dataset or type one in µm"
    end
    refresh_data_tab!(sh)
    refresh_model_tab!(sh)
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
    idx = visible_only ? star.index_quads_visible : eachindex(Imap)
    return Float64.(Imap[idx]) .* Float64.(star.ldmap[idx])
end

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
        star = create_star(model_tessellation(sh), p, t; secondary = m.secondary)
        tmap = parametric_temperature_map(p, star; secondary = m.secondary)
        return (star, tmap)
    catch err
        console!(sh, "preview failed: $(sprint(showerror, err))")
        return nothing
    end
end

"""
    refresh_model_tab!(sh)

Redraw the 3-D preview and the Mollweide from the current model. GUI thread only.
"""
# The decoration flags as keywords for `show_map!`. `:spin` expands into the two the drawing
# layer actually takes.
function _decor(sh::ShellState)
    d = sh.decor
    # Degrees to LINE COUNTS, which is what `graticule_segments` takes: parallels span 180
    # degrees of latitude and meridians 360 of longitude.
    nlat = max(2, round(Int, 180 / sh.gratlat[]))
    nlon = max(2, round(Int, 360 / sh.gratlon[]))
    return (limb = d[:limb], compass = d[:compass], graticules = d[:graticules],
            plotmesh = d[:plotmesh],
            rotation_axis = d[:spin], rotation_arrow = d[:spin],
            graticule_nlat = nlat, graticule_nlon = nlon,
            graticule_color = Symbol(sh.gratcolor[]))
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
    star, tmap = got
    # All three, every time. They are three views of one map, and refreshing only the visible
    # one would mean a stale picture appearing the moment the view is switched — which is
    # worse than the cost, because switching is meant to be instant.
    lbl = sh.intensity[] ? "I (arb.)" : "T (K)"
    sh.msky.cbarlabel[] = lbl
    # ONE evaluation, three views. This called `surface_values` three times — twice with
    # exactly the same arguments — and it runs on every keystroke in the parameter form. With
    # the intensity tick on and Planck selected that is a Planck evaluation over every tessel,
    # done twice for nothing. The visible values are a SLICE of the full ones by construction
    # (`surface_values` differs only in the index it applies), so one call covers all three.
    allv = surface_values(sh, tmap, star; visible_only = false)
    visv = allv[star.index_quads_visible]
    show_map!(sh.msky, star, visv; star_params = star_params(m), _decor(sh)...)
    show_star3d!(sh.star, star, allv)
    show_mollweide!(sh.moll, allv, star)
    return sh
end

"""
    refresh_data_tab!(sh)

Redraw the sky view for the current epoch and the per-epoch χ² panel.
"""
function refresh_data_tab!(sh::ShellState; got = build_epoch_star(sh))
    m = current_model(sh.session)
    if got === nothing
        idle!(sh.sky, current_dataset(sh.session) === nothing ?
                      "no dataset — use Open OIFITS…" :
                      "no model — add one on the Model tab")
    else
        star, tmap = got
        vals = surface_values(sh, tmap, star; visible_only = true)
        sh.sky.cbarlabel[] = sh.intensity[] ? "I (arb.)" : "T (K)"
        d = current_dataset(sh.session)
        ttl = d === nothing ? "" :
              Printf.@sprintf("epoch %d — MJD %.4f  (t = %.3f d)", sh.session.current_epoch,
                              d.mjd[clamp(sh.session.current_epoch, 1, length(d.mjd))],
                              d.tepochs[clamp(sh.session.current_epoch, 1, length(d.tepochs))])
        show_map!(sh.sky, star, vals; title = ttl, star_params = star_params(m),
                  _decor(sh)...)
    end
    c = epoch_chi2(sh)
    c === nothing ? idle!(sh.chi2, "no χ² yet — load a dataset and set a model") :
                    show_chi2!(sh.chi2, c)
    refresh_obs!(sh)
    return sh
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
    p = star_params(m)
    isempty(validate_star_params(p)) || return nothing
    key = (objectid(d), length(d.data), m.secondary, p, sh.nside_exp[], sh.precision[],
           sh.tessel[])
    sh.chi2key[] === key && return sh.chi2cache[]
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
    build_obs_canvas(fig) -> (canvas, model_points)

An OITOOLS `LiveCanvas` plus ONE extra scatter for the model prediction.

The extra series is inserted here, at build time, which is the only moment a plot may be added
— afterwards the window owns the GL context. It is drawn in black rings over OITOOLS' filled
markers so the two read as data and model rather than as two datasets.
"""
function build_obs_canvas(fig)
    O = _oitools_gui()
    O === nothing && return (nothing, Makie.Observable(Makie.Point2f[]))
    ax = Makie.Axis(fig[1, 1])
    c  = O.build_canvas(fig, ax)
    pts = Makie.Observable(Makie.Point2f[])
    Makie.scatter!(ax, pts; color = (:black, 0.0), strokecolor = :black, strokewidth = 1.1,
                   markersize = 9, marker = :circle)
    return (c, pts)
end

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
        push!(rows, "$(strip(kv[1]))\t$(strip(kv[2]))")
    end
    # A file whose tables confuse the check should still be explorable.
    length(rows) == length(always) &&
        append!(rows, ("$(k)\t1" for k in ("v2", "t3amp", "t3phi")))
    return join(rows, "\n")
end

"""
    shell_set_obs_view(kind, color, overlay, logy) -> String

Redraw the observable plot: which quantity, what to colour by, whether to overlay the model's
prediction, and whether the y axis is logarithmic.
"""
function shell_set_obs_view(kind, color, overlay, logy)
    sh = _sh()
    sh.obskind[]    = Symbol(String(kind))
    sh.obscolor[]   = Symbol(String(color))
    sh.obsoverlay[] = String(overlay) == "1"
    sh.obslog[]     = String(logy) == "1"
    return refresh_obs!(sh)
end

"""
    _obs_logscale(sh) -> Bool

Whether the observable plot's y axis is logarithmic for the view it is currently showing.

The tick is remembered across views, but it applies only to V² and T3amp: a closure phase takes
both signs, and uv coverage is a geometry. Decided here rather than only in QML, because the
callback is also reachable from a script and from the tests, and a log axis on a phase is not a
rescaling — it is half the points thrown away.
"""
_obs_logscale(sh::ShellState) = sh.obslog[] && sh.obskind[] in (:v2, :t3amp)

"""
    refresh_obs!(sh) -> String

Push the current epoch's data into the observable plot, and the model's prediction over it.

ONE epoch, not the concatenation: the epochs are separate observations of a rotating star, and
a V² curve stacked across six nights is six curves drawn on top of each other. The epoch is the
one the Data tab's slider selects.
"""
function refresh_obs!(sh::ShellState)
    sh.obs === nothing && return sh.status
    O = _oitools_gui()
    d = current_dataset(sh.session)
    if d === nothing
        sh.obsmodel[] = Makie.Point2f[]
        return sh.status
    end
    i = clamp(sh.session.current_epoch, 1, length(d.data))
    dat = d.data[i]
    # The overlay is emptied BEFORE the data goes in. `update_canvas!` ends by flipping the
    # axis to log10, which transforms every point loaded AT THAT MOMENT — and what is still
    # loaded here is the PREVIOUS view's prediction, a closure phase among them, on which
    # `log10` is a DomainError. Same reason OITOOLS sets `identity` before pushing its series.
    sh.obsmodel[] = Makie.Point2f[]
    try
        # `logscale` is OITOOLS' own keyword: it drops the non-positive points and sets the
        # limits from what is left, so the failure mode it has to survive — a noise-dominated
        # V² at or below zero — is handled once, in the canvas both GUIs share. It throws only
        # when NOTHING is positive, which is a message on the console, not a broken tab.
        O.update_canvas!(sh.obs, dat, sh.obskind[]; color = sh.obscolor[],
                         logscale = _obs_logscale(sh))
    catch err
        console!(sh, "could not draw $(sh.obskind[]): $(sprint(showerror, err))")
        return sh.status
    end
    sh.obsmodel[] = sh.obsoverlay[] ? _model_points(sh, dat, i) : Makie.Point2f[]
    return sh.status
end

"""
    _model_points(sh, dat, epoch) -> Vector{Point2f}

The current model's prediction for the current view, positioned on the same x axis as the data.

Empty when there is no model, when the parameters do not build, or when the view is `uv` —
uv coverage is geometry, not a prediction, so there is nothing for a model to say about it.
"""
function _model_points(sh::ShellState, dat, epoch::Int)
    m = current_model(sh.session)
    m === nothing && return Makie.Point2f[]
    kind = sh.obskind[]
    kind in (:v2, :t3amp, :t3phi) || return Makie.Point2f[]
    # THE CACHED geometry, not a fresh one. This is called on every view switch and every
    # epoch step, and rebuilding the mesh at all epochs and re-running `setup_oi!` for it is
    # the same work `epoch_chi2` has already done with the same parameters.
    st = model_state(sh)
    st === nothing && return Makie.Point2f[]
    try
        v2m, t3am, t3pm = observables(st.x, st.stars[epoch], dat)
        y, xs = kind === :v2    ? (v2m,  dat.v2_baseline) :
                kind === :t3amp ? (t3am, dat.t3_baseline) :
                                  (t3pm, dat.t3_baseline)
        # The same 1e-6 Mλ scaling OITOOLS' ObsSpec applies, so the two series land on one axis.
        pts = [Makie.Point2f(Float32(xs[k]) * 1f-6, Float32(y[k])) for k in eachindex(y)]
        # A log axis transforms the overlay too, and a model V² of exactly zero — a visibility
        # null — transforms to -Inf, which is not a limit Makie's log axis accepts. The data
        # series drops its non-positive points inside `canvas_data`; the overlay drops its own,
        # or the whole panel fails with nothing but "exception in render" to say why.
        return _obs_logscale(sh) ? filter(p -> p[2] > 0, pts) : pts
    catch err
        console!(sh, "model overlay failed: $(sprint(showerror, err))")
        return Makie.Point2f[]
    end
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
    for c in (sh.sky, sh.msky, sh.star, sh.moll, sh.imsky, sh.immoll)
        ok |= set_colormap!(c, String(name))
    end
    # The Mollweide computes its colours from the map, so it has to be redrawn rather than
    # re-tinted; the sky canvases keep their values and `set_colormap!` remaps them in place.
    ok && refresh_model_tab!(sh)
    ok && refresh_data_tab!(sh)
    return ok ? "colormap: $(name)" : "unknown colormap $(name)"
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
    for c in (sh.sky, sh.msky, sh.imsky); reset_zoom!(c); end
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
    for c in (sh.sky, sh.msky, sh.imsky)
        scene = c.axis.scene
        ev = Makie.events(scene)
        downat = Ref((0.0, 0.0))
        Makie.on(ev.mousebutton, priority = 100) do e
            pos = Tuple(Float64.(ev.mouseposition[]))
            if e.button == Makie.Mouse.right
                if e.action == Makie.Mouse.press
                    downat[] = pos
                elseif e.action == Makie.Mouse.release &&
                       hypot(pos[1] - downat[][1], pos[2] - downat[][2]) < 5
                    reset_zoom!(c)
                end
            end
            return Makie.Consume(false)
        end
        # The clamp runs on the LIMITS, not on the scroll event: a drag-zoom never produces a
        # scroll, and clamping only the wheel would leave the other gesture unbounded.
        Makie.on(c.axis.finallimits; priority = 1) do _
            clamp_zoom!(c)
            return Makie.Consume(false)
        end
    end
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

function shell_job_stop()
    sh = _sh()
    j = sh.job
    j === nothing && return "nothing running"
    j.stop[] = true
    return "stopping…"
end

"""
    shell_job_poll() -> String

Called by a QML timer every 200 ms. `running\telapsed\toutput`.

When the worker has finished, this is where the result is taken up and the canvases redrawn —
on the GUI thread, which is the only place a GL call may happen.
"""
function shell_job_poll()
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
    ("gradient",   "VMLMB + analytic gradient (fast)"),
    ("neldermead", "Nelder–Mead (NLopt, local)"),
    ("bobyqa",     "BOBYQA (NLopt, local, quadratic model)"),
    ("hmc",        "NUTS (Hamiltonian, uses the analytic gradient)"),
    ("nautilus",   "Nautilus (importance nested sampling, pure Julia)"),
    ("pigeons",    "Pigeons (parallel tempering — for a multimodal posterior)"),
]
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
    for (k, v) in FIT_METHODS
        k == "nautilus" && !nautilus_available() && continue
        # Tempering runs on the same parametric log-posterior NUTS does, so it is offered in
        # the same place: the rapid rotator, with the gradient stack loaded.
        k == "pigeons" && !(pigeons_available() && m !== nothing && m.surface_type == 2) &&
            continue
        # NUTS needs the analytic gradient, so it is offered exactly where the gradient path
        # is: the rapid rotator, with Zygote loaded.
        k == "hmc" && !(hmc_available() && m !== nothing && m.surface_type == 2) && continue
        # Offered only where the gradient is CONSISTENT with the objective being minimised.
        # See `gradient_fit_kind` — this is not the same set as "has analytic shape
        # derivatives", and the difference is measurable.
        k == "gradient" && (m === nothing || gradient_fit_kind(m) === :none) && continue
        push!(rows, "$(k)\t$(v)")
    end
    return join(rows, "\n")
end

"Whether ROTIRZygoteExt has loaded, i.e. whether `fit_parametric` exists."
fit_parametric_available() = !isempty(methods(ROTIR.fit_parametric))

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
const PARAMETRIC_THETA = Dict(
    :rpole => "rpole", :frac_escapevel => "omega", :inclination => "inc",
    :position_angle => "PA", :beta => "beta", :ld1 => "ld1", :ld2 => "ld2",
    :tpole => "tpole")

"""
    shell_fit(method, maxeval) -> String

Fit the FREE parameters of the current model against every epoch of the current dataset.

Generic over `parametric_chi2` rather than routed to `fit_sphere_ld` / `fit_ellipsoid_ld`:
those take their own fixed parameter orders per shape, while the GUI's model is "whichever
fields the user marked free", which is exactly what this walks. Tied parameters are re-derived
inside the objective, so a tie holds at every trial point rather than only at the start.
"""
function shell_fit(method, maxeval)
    sh = _sh()
    d = current_dataset(sh.session)
    d === nothing && return "no dataset"
    m = current_model(sh.session)
    m === nothing && return "no model"
    names = sort(collect(m.free))
    isempty(names) && return "nothing is free — mark at least one parameter free"
    meth = Symbol(String(method))
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
        m.surface_type == 2 ||
            return "Pigeons runs on the parametric rapid-rotator model (surface_type 2)"
        bad = filter(n -> !haskey(PARAMETRIC_THETA, n), names)
        isempty(bad) || return "not parameters of the tempering path: " * join(bad, ", ")
    end
    if meth === :hmc
        hmc_available() ||
            return "NUTS needs AdvancedHMC, LogDensityProblems and Zygote in this session"
        m.surface_type == 2 ||
            return "NUTS runs on the parametric rapid-rotator model (surface_type 2)"
        bad = filter(n -> !haskey(PARAMETRIC_THETA, n), names)
        isempty(bad) || return "not parameters of the NUTS path: " * join(bad, ", ")
    end
    if meth === :gradient
        kind = gradient_fit_kind(m)
        kind === :none && return _no_gradient_reason(m)
        known = kind === :parametric ? keys(PARAMETRIC_THETA) : keys(SHAPE_THETA[m.surface_type])
        bad = filter(n -> !(n in known), names)
        isempty(bad) || return "not parameters of the gradient path for this surface type: " *
                              join(bad, ", ") * " (free: " * join(sort(collect(known)), ", ") * ")"
    end

    θ0 = [m.params[n] for n in names]
    lb = [m.bounds[n][1] for n in names]
    ub = [m.bounds[n][2] for n in names]
    any(θ0 .< lb) || any(θ0 .> ub) &&
        return "a starting value is outside its bounds"

    # Snapshot the model: the worker must not read a struct the GUI thread can edit under it.
    snap = ModelEntry(m.name, m.surface_type, copy(m.params), copy(m.free), copy(m.bounds),
                      copy(m.ties), m.secondary)
    data = d.data
    tepochs = copy(d.tepochs)
    nd = sum(dd.nv2 + dd.nt3amp + dd.nt3phi for dd in data)
    it = max(10, Int(maxeval))

    # THE PANEL'S MESH, not a hardcoded one. This built `tessellation_healpix(4)` — nside 16,
    # 3072 tessels — while the mesh box above it read "n = 3, 768 tessels", so every fit did
    # four times the work per likelihood that the panel said it would, and a fit that takes
    # seconds from a script took minutes here. Captured before the job starts: the worker must
    # not read a Ref the GUI thread can change under it.
    nexp = sh.nside_exp[]
    prec = sh.precision[]
    console!(sh, "fit $(join(names, ", ")) by $(meth), $(it) evaluations, " *
                 "HEALPix level $(nexp) ($(12 * (2^nexp)^2) tessels, $(prec))"; kind = :cmd)
    return start_job!(sh, :fit, function (stop)
        tess = tessellation_healpix(nexp; T = prec)
        obj = function (θ)
            stop[] && return 1e30
            for (k, n) in enumerate(names); snap.params[n] = θ[k]; end
            p = star_params(snap)
            isempty(validate_star_params(p)) || return 1e30
            c = try
                parametric_chi2(p, tess, data, tepochs)
            catch
                return 1e30
            end
            return isfinite(c) ? c : 1e30
        end
        best, chi2, extra, post = if meth === :hmc
            _run_hmc_fit(snap, data, tepochs, names, θ0, lb, ub, nexp, prec, it)
        elseif meth === :pigeons
            _run_pigeons_fit(snap, data, tepochs, names, θ0, lb, ub, nexp, prec, it)
        elseif meth !== :gradient
            _run_fit(meth, obj, θ0, lb, ub, names, it)
        elseif gradient_fit_kind(snap) === :parametric
            _run_gradient_fit(snap, data, tepochs, names, θ0, lb, ub, nexp, prec, it)
        else
            _run_shape_fit(snap, data, tepochs, names, θ0, lb, ub, nexp, prec, it)
        end
        for (k, n) in enumerate(names); snap.params[n] = best[k]; end
        # The WHOLE result travels back, not a summary of it. `job_succeeded!` turns this into
        # a `FitEntry` on the GUI thread; nothing is reduced here, on the worker, where the
        # draws would otherwise go out of scope and be collected.
        return (; params = Dict(n => best[k] for (k, n) in enumerate(names)),
                table = _fit_table(names, best, extra),
                names = names, best = best, errs = extra, post = post,
                method = meth, model = snap.name, surface_type = snap.surface_type,
                chi2 = chi2, ndata = nd,
                status = Printf.@sprintf("fit done: χ²ᵣ = %.4f over %d points", chi2 / nd, nd))
    end)
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
function _run_gradient_fit(snap, data, tepochs, names, θ0, lb, ub, nexp, prec, maxeval)
    θnames = ROTIR.parametric_param_names(; tpole_free = true)
    # Full theta from the model, INCLUDING the frozen entries: `fit_parametric` wants the whole
    # vector and moves only the indices named in `free`.
    getp(k) = Float64(get(snap.params, k, 0.0))
    θfull = Float64[getp(:rpole), getp(:frac_escapevel), getp(:inclination),
                    getp(:position_angle), getp(:beta), getp(:ld1), getp(:ld2), getp(:tpole)]
    lo = copy(ROTIR.default_parametric_bounds(; tpole_free = true)[1])
    hi = copy(ROTIR.default_parametric_bounds(; tpole_free = true)[2])
    free = String[]
    for (k, n) in enumerate(names)
        haskey(PARAMETRIC_THETA, n) || error("$(n) is not a parameter of the gradient path")
        j = findfirst(==(PARAMETRIC_THETA[n]), θnames)
        push!(free, θnames[j])
        # The panel's bounds win over the defaults; `default_parametric_bounds` carries an
        # infinite upper bound on rpole and tpole, which a box-constrained search accepts but
        # which throws away whatever the user typed.
        lo[j] = lb[k]; hi[j] = ub[k]
    end
    tess = tessellation_healpix(nexp; T = prec)
    base = star_params(snap)
    θ̂, chi2r, info = ROTIR.fit_parametric(data, tess, tepochs, base;
                                          θ0 = θfull, free = free, lb = lo, ub = hi,
                                          tpole_free = true, maxiter = max(20, maxeval ÷ 50),
                                          verb = true)
    best = [θ̂[findfirst(==(PARAMETRIC_THETA[n]), θnames)] for n in names]
    nd = sum(d -> d.nv2 + d.nt3amp + d.nt3phi, data)
    return best, chi2r * nd, Dict{Symbol,Float64}(), nothing
end

"""
    _run_hmc_fit(snap, data, tepochs, names, θ0, lb, ub, nexp, prec, maxeval) -> (best, chi2, extra)

NUTS over the free parametric parameters, returning the posterior MEDIAN as the point estimate
and half the 16–84 interval as the error.

The median, not the maximum: a Hamiltonian sampler is run for the posterior, and the mode of a
skewed one is not the number to quote beside a symmetric error bar. `maxeval` is split into
adaptation and draws in the usual 3:4 ratio — adaptation is where the step size and the mass
matrix are set, and starving it is the commonest way to get divergences.
"""
function _run_hmc_fit(snap, data, tepochs, names, θ0, lb, ub, nexp, prec, maxeval)
    θnames = ROTIR.parametric_param_names(; tpole_free = true)
    getp(k) = Float64(get(snap.params, k, 0.0))
    θfull = Float64[getp(:rpole), getp(:frac_escapevel), getp(:inclination),
                    getp(:position_angle), getp(:beta), getp(:ld1), getp(:ld2), getp(:tpole)]
    lo, hi = ROTIR.default_parametric_bounds(; tpole_free = true)
    lo = copy(lo); hi = copy(hi)
    free = String[]
    for (k, n) in enumerate(names)
        j = findfirst(==(PARAMETRIC_THETA[n]), θnames)
        push!(free, θnames[j])
        lo[j] = lb[k]; hi[j] = ub[k]
    end
    tess = tessellation_healpix(nexp; T = prec)
    base = star_params(snap)
    ndraw = clamp(maxeval ÷ 7, 100, 5000)
    r = ROTIR._fit_hmc(data, tess, tepochs, base; θ0 = θfull, free = free, lb = lo, ub = hi,
                       tpole_free = true, n_samples = 4ndraw ÷ 4, n_adapt = 3ndraw ÷ 4,
                       verb = true)
    # `_fit_hmc` returns the free parameters in the order `parametric_free_indices` sorts them,
    # which is not the order the panel listed them in.
    order = sortperm([findfirst(==(PARAMETRIC_THETA[n]), θnames) for n in names])
    best = zeros(length(names)); errs = Dict{Symbol,Float64}()
    for (slot, k) in enumerate(order)
        best[k] = r.median[slot]
        errs[names[k]] = (r.q84[slot] - r.q16[slot]) / 2
    end
    # `invperm`, not `order`: `order` maps a panel slot to a sampler column, and the columns
    # have to be permuted the other way to land in panel order.
    return best, NaN, errs,
           _posterior(r; order = invperm(order),
                      diagnostics = Printf.@sprintf("%d draws, %d divergences",
                                                    size(r.samples, 1), r.divergences))
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
function _run_pigeons_fit(snap, data, tepochs, names, θ0, lb, ub, nexp, prec, maxeval)
    θnames = ROTIR.parametric_param_names(; tpole_free = true)
    getp(k) = Float64(get(snap.params, k, 0.0))
    θfull = Float64[getp(:rpole), getp(:frac_escapevel), getp(:inclination),
                    getp(:position_angle), getp(:beta), getp(:ld1), getp(:ld2), getp(:tpole)]
    lo, hi = ROTIR.default_parametric_bounds(; tpole_free = true)
    lo = copy(lo); hi = copy(hi)
    free = String[]
    for (k, n) in enumerate(names)
        j = findfirst(==(PARAMETRIC_THETA[n]), θnames)
        push!(free, θnames[j])
        lo[j] = lb[k]; hi[j] = ub[k]
    end
    tess = tessellation_healpix(nexp; T = prec)
    base = star_params(snap)
    # 2^n_rounds scans. Clamped to 6..12: below six the ladder has not tuned itself and the
    # answer is noise, above twelve one run is longer than a working session.
    nr = clamp(round(Int, log2(max(maxeval, 64))) - 5, 6, 12)
    r = ROTIR._fit_pigeons(data, tess, tepochs, base; θ0 = θfull, free = free, lb = lo, ub = hi,
                           tpole_free = true, n_rounds = nr, n_chains = 8,
                           explorer = :slice, multithreaded = Threads.nthreads() > 1,
                           verb = true)
    Printf.@printf("Pigeons: log(Z) = %.4f, %d round trips over 2^%d scans\n",
                   r.logz, r.round_trips, nr)
    order = sortperm([findfirst(==(PARAMETRIC_THETA[n]), θnames) for n in names])
    best = zeros(length(names)); errs = Dict{Symbol,Float64}()
    for (slot, k) in enumerate(order)
        best[k] = r.median[slot]
        errs[names[k]] = (r.q84[slot] - r.q16[slot]) / 2
    end
    return best, NaN, errs,
           _posterior(r; order = invperm(order),
                      diagnostics = Printf.@sprintf("%d draws, %d chains, %d round trips",
                                                    size(r.samples, 1), r.n_chains,
                                                    r.round_trips))
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
function _run_shape_fit(snap, data, tepochs, names, θ0, lb, ub, nexp, prec, maxiter)
    layout = SHAPE_THETA[snap.surface_type]
    n = maximum(values(layout))
    inv = Dict(v => k for (k, v) in layout)
    θ  = Float64[Float64(get(snap.params, inv[j], 0.0)) for j in 1:n]
    lo = fill(-Inf, n); hi = fill(Inf, n)
    # Only the named parameters move. VMLMB has no frozen-variable facility, so the others are
    # pinned by collapsing their bounds onto their current value — which is exact, and cheaper
    # than the scatter-matrix reduction the Zygote path uses.
    for j in 1:n; lo[j] = θ[j]; hi[j] = θ[j]; end
    for (k, nm) in enumerate(names)
        j = layout[nm]
        lo[j] = lb[k]; hi[j] = ub[k]
        θ[j]  = clamp(θ[j], lb[k], ub[k])
    end
    # Float64 throughout: `shape_chi2_fg!` is written with ONE element type `T` shared by the
    # gradients, the map, θ and the tessellation, so the default Float32 mesh and a Float64 θ
    # do not meet a method. The LEVEL still follows the panel — only the element type is
    # forced here, and `prec` is deliberately ignored for that reason.
    tess = tessellation_healpix(nexp; T = Float64)
    base = star_params(snap)
    star = create_star(tess, base, tepochs[1]; secondary = snap.secondary)
    xmap = Float64.(parametric_temperature_map(base, star; secondary = snap.secondary))
    gθ = zeros(Float64, n); gx = zeros(Float64, length(xmap))
    npts = sum(d -> d.nv2 + d.nt3amp + d.nt3phi, data)
    calls = Ref(0)
    function fg!(x, g)
        calls[] += 1
        c = shape_chi2_fg!(gθ, gx, xmap, collect(Float64, x), data, tess, base, tepochs)
        g .= gθ ./ npts
        calls[] % 10 == 1 && Printf.@printf("  it %3d   χ²ᵣ = %.6f\n", calls[], c / npts)
        return c / npts
    end
    θ̂ = OptimPackNextGen.vmlmb(fg!, θ; lower = lo, upper = hi, maxiter = max(20, maxiter ÷ 50),
                               blmvm = false, verb = false)
    best = [θ̂[layout[nm]] for nm in names]
    c = shape_chi2_fg!(gθ, gx, xmap, θ̂, data, tess, base, tepochs)
    return best, c, Dict{Symbol,Float64}(), nothing
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
        un = _fit_nested(meth, obj, String.(names), lb, ub; verb = true,
                         n_eff = clamp(maxeval, 200, 50_000),
                         n_live = clamp(maxeval ÷ 4, 100, 2000),
                         min_num_live_points = clamp(maxeval ÷ 4, 100, 2000))
        return un.median, obj(un.median),
               Dict(n => ((un.q84[k] - un.q16[k]) / 2) for (k, n) in enumerate(names)),
               _posterior(un; diagnostics = Printf.@sprintf("%d samples", size(un.samples, 1)))
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
    return best, chi2, Dict{Symbol,Float64}(), nothing
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

Every regulariser `spheroid_regularization` dispatches on, with what the panel needs to offer
it: a default weight, and the ONE extra scalar (if any) the caller has to choose.

Taken from the dispatch chain itself (src/oichi2_spheroid.jl:995-1019) rather than from the
docs, and it is longer than it looks: `mean`, `bias` and `radialvar` are implemented and were
missing from the first version of this table, which would have made them unreachable from the
GUI while `spheroid_regularization` happily accepted them from a script.

A regulariser entry is FOUR elements — `[name, weight, aux, subset]` — and only the first two
are numbers the panel can supply:

  * `aux` is a precomputed STRUCTURE for most of them: the Sobel operator for `sobel`/`sobel2`,
    the neighbour Laplacian for `tv`/`tv2`, the radial binning for `radflat`/`radialvar`, the
    degenerate direction for `orthold`. Only `bias` takes a plain number there.
  * `subset` is a PIXEL INDEX SET, not a number, and for the radial regularisers it must be
    the exact index set their binning was built from or the call errors.

So the panel chooses names, weights and the one scalar knob, and [`build_regularizers`](@ref)
constructs the rest at run time. Passing four numbers — which is what a naive `name:a:b:c`
string invites — cannot work at all.

Fields: name, default weight, extra-knob label ("" when there is none), extra default, doc.
"""
const REGULARIZER_KINDS = [
    ("sobel",     1e1, "",      0.0, "∫|∇x| dΩ — isotropic L1 on the sphere, edge-preserving"),
    ("sobel2",    1e1, "",      0.0, "∫|∇x|² dΩ — smooth; the usual first choice"),
    ("tv",        1e1, "",      0.0, "‖Lx‖ — curvature, L1"),
    ("tv2",       1e1, "",      0.0, "‖Lx‖² — curvature, L2"),
    ("mem",       1e1, "",      0.0, "maximum entropy: per-pixel contrast, no spatial coupling"),
    ("mean",      1e1, "",      0.0, "departure from the map mean"),
    ("bias",      1e1, "B",     2.0, "harmonic bias for asymmetric brightening; B is the asymmetry"),
    ("radflat",   1e2, "nbins", 6.0, "flatness BETWEEN annuli — single-epoch, non-rotating only"),
    ("radialvar", 1e2, "nbins", 6.0, "variance WITHIN annuli — the complement of radflat"),
    ("orthold",   1e2, "",      0.0, "penalises the limb-darkening-degenerate direction"),
]

"""
    shell_regularizer_kinds() -> String

`name\tweight\textra_label\textra_default\tdoc` per line. The panel shows the extra field only
where `extra_label` is non-empty, which is what keeps eight of the ten rows to two controls.
"""
shell_regularizer_kinds() =
    join(("$(n)\t$(w)\t$(el)\t$(ed)\t$(doc)" for (n, w, el, ed, doc) in REGULARIZER_KINDS), "\n")

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
    # Anything a regulariser has to CONSTRUCT first goes on its own line before the array.
    # Putting it in a trailing comment inside the literal (which the first version did) makes
    # the following comma part of the comment, and the exported script does not parse.
    pre   = String[]
    parts = String[]
    for s in specs
        if s.name in ("sobel", "sobel2")
            push!(parts, "[\"$(s.name)\", $(s.weight), sobel, 1:length(x0)]")
            "sobel = sobel_gradient_healpix($n)" in pre ||
                push!(pre, "sobel = sobel_gradient_healpix($n)")
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

Two subsets are in play and they are not the same. The mesh regularisers act on every pixel;
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
            # several calls down, which names neither the regulariser nor the missing step.
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
    console!(sh, "image_reconstruct_oi(...; nside=2^$n, maxiter=$it)"; kind = :cmd)
    log!(sh.session,
         "tess  = tessellation_healpix($n)\n" *
         "stars = create_star_multiepochs(tess, $(m.name), tepochs)\n" *
         # `setup_oi!` BEFORE the regularisers, and it was missing entirely: the worker does
         # it, the exported script did not, and `radflat_bins` then indexed an empty
         # `star.polyflux` — a script that parsed, read correctly, and died on the fourth line.
         # Nothing caught it until the round trip was actually RUN in the test suite.
         "setup_oi!(data, stars)\n" *
         "x0    = parametric_temperature_map($(m.name), stars[1])\n" *
         _regularizer_source(specs, n) * "\n" *
         "x = image_reconstruct_oi(x0, data, stars; regularizers = regs, maxiter = $it)";
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
                                  callback = cb, callback_every = 25)
        return (; x, n, regs, stars, data)
    end; progress = prog)
end

"""
    job_succeeded!(sh, kind, result)

Take up a finished worker's result. **GUI thread only** — this is where the GL work happens.
"""
function job_succeeded!(sh::ShellState, kind::Symbol, result)
    if kind === :image
        b = chi2_breakdown(result.x, result.stars, result.data)
        e = add_image!(sh.session, "map$(length(sh.session.images) + 1)", result.x,
                       result.n, b.total, b.ndata, result.regs,
                       something(current_dataset(sh.session), (; name = "")).name)
        sh.status = Printf.@sprintf("%s: χ²ᵣ = %.4f over %d points", e.name, b.totalr, b.ndata)
        console!(sh, sh.status)
        show_reconstruction!(sh, result.stars[1], result.x; title = e.name)
        show_chi2!(sh.chi2, b.epochs)
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
    elseif kind === :fit
        sh.lastfit = result.table
        sh.status = result.status
        console!(sh, sh.status)
        _record_fit!(sh, result)
        m = current_model(sh.session)
        if m !== nothing
            for (k, v) in result.params; m.params[k] = v; end
            refresh_model_tab!(sh)
            refresh_data_tab!(sh)
        end
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
    show_map!(sh.imsky, star, Float64.(x[star.index_quads_visible]);
              title = title, star_params = sp, _decor(sh)...)
    show_mollweide!(sh.immoll, Float64.(x), star)
    sh.imstar === nothing || show_star3d!(sh.imstar, star, Float64.(x))
    # And the Model tab's own 3-D scene, which is where the recovered map is compared with the
    # parametric one it started from.
    show_star3d!(sh.star, star, Float64.(x))
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
    refresh_data_tab!(sh; got = got)
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
                 post === nothing ? "" : post.diagnostics)
    add_fit!(sh.session, e)
    isempty(S) ||
        console!(sh, "posterior kept: $(size(S, 1)) draws x $(size(S, 2)) parameters" *
                     (isempty(e.diagnostics) ? "" : "  ($(e.diagnostics))"))
    isfinite(e.logz) &&
        console!(sh, Printf.@sprintf("log(Z) = %.4f ± %.4f", e.logz, e.logzerr))
    return e
end

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
        push!(rows, join((f.name, String(f.method), f.model, string(f.surface_type),
                          f.ndata > 0 ? Printf.@sprintf("%.4f", f.chi2 / f.ndata) : "—",
                          isfinite(f.logz) ? Printf.@sprintf("%.3f", f.logz) : "—",
                          isfinite(f.logzerr) ? Printf.@sprintf("%.3f", f.logzerr) : "—",
                          string(size(f.samples, 1)), f.diagnostics), "\t"))
    end
    return join(rows, "\n")
end

"Which fit the posterior panel is showing, 1-based; `0` when there is none."
shell_current_fit() = string(_sh().session.current_fit)

function shell_select_fit(i)
    sh = _sh()
    n = Int(i)
    n in eachindex(sh.session.fits) || return sh.status
    sh.session.current_fit = n
    refresh_posterior!(sh)
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

shell_images() = join((Printf.@sprintf("%s\t%d\t%.4f\t%d\t%s", e.name, e.nside_exp,
                                       e.ndata > 0 ? e.chi2 / e.ndata : NaN, e.ndata, e.source)
                       for e in _sh().session.images), "\n")

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
    v = set_plot_scale!(something(tryparse(Float64, String(x)), 0.0))
    sh = SHELL[]
    if sh !== nothing
        refresh_both!(sh)
    end
    return v == 0 ? "plot scale: from the screen" : Printf.@sprintf("plot scale: %.2f", v)
end

"""
    shell_plot_scale() -> String

`user_override\tin_force` — the panel shows the first and displays the second beside it, so
"auto" and "1.19" are distinguishable.
"""
shell_plot_scale() = Printf.@sprintf("%.3f\t%.3f", PLOT_SCALE_USER[], live_plot_scale())

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
second save, which silently destroyed the first: reconstruct, look, change a regulariser,
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
