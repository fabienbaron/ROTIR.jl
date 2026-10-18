# The Explore tab's Filter panel: OITOOLS' data filtering, applied epoch by epoch.
#
# The panel and its settings format are OITOOLS' own (`FilterPanel.qml`, `filter_spec_kwargs`),
# so a cut means the same thing in both GUIs. What is ROTIR's is the dataset: a vector of
# epochs rather than OITOOLS' wavelength × time grid, so every epoch is filtered by its own
# discard lists and the result keeps its epoch structure — and its time origin, since the
# rotation phase is measured from it.

# The unfiltered dataset a filter applies to: the current one, or the one it was cut from.
# A filter's settings describe the whole cut, so applying them again replaces the cut rather
# than cutting what an earlier one left.
function _filter_root(s::Session)
    e = current_dataset(s)
    for _ in 1:length(s.datasets)
        (e === nothing || isempty(e.source)) && break
        i = findfirst(x -> x.name == e.source, s.datasets)
        i === nothing && break
        e = s.datasets[i]
    end
    return e
end

_npoints(d) = d.nv2 + d.nt3amp + d.nt3phi + d.nvisamp + d.nvisphi + d.nflux

"""
    filter_dataset!(session, root; read, prep, drop, spec, kwargs...) -> DatasetEntry

Cut dataset `root` into a new dataset, leaving `root` untouched, and make the new one current.

- `kwargs` go to `set_data_filter`, epoch by epoch.
- `read`: `readoifits` options (instrument, target, tables, flags, channel averaging). They are
  decided when a file is read, so with any given the files are read again and split into
  epochs as `root` was.
- `prep`: error-bar floors and multipliers for `oifits_prep`, applied after the cut.
- `drop`: observables `filter_data` removes whole (`:t3amp`, `:visphi`, …).
- `hand`: points picked by hand, per epoch of `root`: `Dict(epoch => Dict(:remove_v2 => [i…]))`
  with the `set_data_filter` keywords as keys. When the files are read again they are found in
  the new epochs by their uv, wavelength and stations; `hand_lost` counts those that were not
  (channel averaging changes a point's wavelength, so it cannot be combined with them), and
  `hand_used` receives the lists as applied, indexed into the epochs that were filtered.

An epoch left with no points is removed. The times stay relative to `root`'s first epoch, so a
rotation phase means the same thing on the filtered dataset as on the one it came from.
"""
function filter_dataset!(s::Session, root::DatasetEntry;
                         read::AbstractDict = Dict{Symbol,Any}(),
                         prep::AbstractDict = Dict{Symbol,Any}(),
                         drop::AbstractVector = Symbol[],
                         hand::AbstractDict = Dict{Int,Dict{Symbol,Vector{Int}}}(),
                         hand_lost::Base.RefValue{Int} = Ref(0),
                         hand_used::Base.RefValue{Any} = Ref{Any}(nothing),
                         spec::AbstractString = "", kwargs...)
    ropts = merge(root.readopts, Dict{Symbol,Any}(read))
    base, mjd = isempty(read) ? (root.data, root.mjd) :
                read_epochs(root.paths; split = root.split, ropts...)
    hb = isempty(read) ? hand : _remap_hand(root, base, mjd, hand, hand_lost)
    hand_used[] = hb
    out = Any[OITOOLS.filter_data(d, OITOOLS.set_data_filter(d; kwargs...,
                                                             get(hb, i, Dict{Symbol,Any}())...);
                                  drop = Symbol.(drop))
              for (i, d) in enumerate(base)]
    isempty(prep) || foreach(d -> OITOOLS.oifits_prep(d; prep...), out)
    keep = findall(d -> _npoints(d) > 0, out)
    isempty(keep) && throw(ArgumentError("the filter removes every point of $(root.name)"))
    out = out[keep]; mjd = Float64.(mjd[keep])
    stem = root.name * " › filtered"
    name, k = stem, 1
    while any(x -> x.name == name, s.datasets)
        k += 1; name = "$stem $k"
    end
    e = DatasetEntry(root.paths, name, out, mjd, mjd .- minimum(root.mjd),
                     [d.nv2 for d in out], [d.nt3phi for d in out], root.split,
                     ropts, root.name, String(spec))
    push!(s.datasets, e)
    s.current_dataset = length(s.datasets)
    s.current_epoch = 1
    return e
end

# The `set_data_filter` keyword each hand-picked kind becomes, and the observable slot of the
# selection tools (`selection_indices`: 1 uv, 2 VIS, 3 V², 4 T3, 5 flux) it removes from.
const HAND_KEYS = ("remove_uv" => :uv_bad, "remove_vis" => :remove_vis, "remove_v2" => :remove_v2,
                   "remove_t3" => :remove_t3, "remove_flux" => :remove_flux)

"""
    _split_hand(spec) -> (spec, hand) or an error String

Take the hand-picked lines out of the Filter panel's settings. ROTIR writes them with the epoch
each point belongs to, `remove_v2\t2:17,2:40,5:3` (epochs of the unfiltered dataset), which
OITOOLS' `filter_spec_kwargs` — one dataset, plain indices — would refuse. What is left is the
spec OITOOLS reads; `hand` is `Dict(epoch => Dict(keyword => indices))`.
"""
function _split_hand(spec::AbstractString)
    keep = String[]
    hand = Dict{Int,Dict{Symbol,Vector{Int}}}()
    for line in split(String(spec), '\n')
        f = split(line, '\t'; limit = 2)
        j = findfirst(p -> p.first == strip(f[1]), HAND_KEYS)
        if j === nothing || length(f) < 2
            push!(keep, line)
            continue
        end
        kw = HAND_KEYS[j].second
        for x in split(f[2], ',')
            isempty(strip(x)) && continue
            ab = tryparse.(Int, split(strip(x), ':'))
            (length(ab) == 2 && !any(isnothing, ab)) ||
                return "! $(strip(f[1])): $(repr(strip(x))) is not an epoch:point pair"
            push!(get!(get!(hand, ab[1], Dict{Symbol,Vector{Int}}()), kw, Int[]), ab[2])
        end
    end
    return join(keep, "\n"), hand
end

# Hand-picked points of `root`'s epochs, found in `base` — the same files read again with other
# options — by their identity: uv, wavelength and stations (OITOOLS' `_point_keys`). An epoch is
# matched by its MJD, since a read option can empty one and shift the rest.
function _remap_hand(root, base, mjd, hand, lost::Base.RefValue{Int})
    O = _oitools_gui()
    out = Dict{Int,Dict{Symbol,Vector{Int}}}()
    for (e, kws) in hand, (kw, idx) in kws
        slot = findfirst(p -> p.second === kw, HAND_KEYS)
        j = findfirst(m -> isapprox(m, root.mjd[e]; atol = 1e-6), mjd)
        if j === nothing
            lost[] += length(idx)
            continue
        end
        rk = O._point_keys(root.data[e], slot)
        where = Dict(k => n for (n, k) in enumerate(O._point_keys(base[j], slot)))
        for i in idx
            n = get(where, rk[i], nothing)
            n === nothing ? (lost[] += 1) :
                push!(get!(get!(out, j, Dict{Symbol,Vector{Int}}()), kw, Int[]), n)
        end
    end
    return out
end

# The filter as script lines, applied to `data` as the log last bound it.
function _log_filter!(s::Session, root::DatasetEntry, e::DatasetEntry, k,
                      hand::AbstractDict = Dict{Int,Dict{Symbol,Vector{Int}}}())
    lines = String[]
    if !isempty(k.read)
        ro = join(("$(a) = $(_literal(v))" for (a, v) in sort(collect(e.readopts), by = first)), ", ")
        push!(lines, length(root.paths) == 1 ?
            "data = split_epochs(readoifits($(_literal(root.paths[1])); $ro, warn = false, verbose = false)[1,1])[1]" :
            "data = [readoifits(f; $ro, warn = false, verbose = false)[1, 1] for f in $(_literal(root.paths))]")
    end
    kw = join(("$(a) = $(_literal(v))" for (a, v) in sort(collect(k.filter), by = first)), ", ")
    dr = isempty(k.drop) ? "" : "; drop = $(_literal(Symbol.(k.drop)))"
    if isempty(hand)
        push!(lines, "data = map(d -> filter_data(d, set_data_filter(d" *
                     (isempty(kw) ? "" : "; $kw") * ")$dr), data)")
    else
        # The points picked by hand, epoch by epoch.
        hs = join(("$(i) => (; " * join(("$(a) = $(_literal(sort(unique(v))))"
                                         for (a, v) in sort(collect(h), by = first)), ", ") * ")"
                   for (i, h) in sort(collect(hand), by = first)), ", ")
        push!(lines, "hand = Dict($hs)")
        push!(lines, "data = [filter_data(d, set_data_filter(d" * (isempty(kw) ? "; " : "; $kw, ") *
                     "get(hand, i, (;))...)$dr) for (i, d) in enumerate(data)]")
    end
    if !isempty(k.prep)
        pr = join(("$(a) = $(_literal(v))" for (a, v) in sort(collect(k.prep), by = first)), ", ")
        push!(lines, "foreach(d -> oifits_prep(d; $pr), data)")
    end
    push!(lines, "data = filter(d -> d.nv2 + d.nt3amp + d.nt3phi + d.nvisamp + d.nvisphi + d.nflux > 0, data)")
    push!(lines, "tepochs = $(_literal(e.tepochs))")
    log!(s, join(lines, "\n"); note = "dataset $(e.name)", binding = "data")
end

"""
    shell_filter_choices() -> String

What the Filter panel offers, from the unfiltered dataset, as the `kind\\t…` lines OITOOLS'
panel reads: `extent_wav` (µm), `extent_bl` (Mλ), `extent_blm` (m), `night\\tlo:hi\\tlabel`,
`insname`, `target`, `station`, `baseline`, `triangle`, and `has\\tobservable`.
"""
function shell_filter_choices()
    sh = _sh()
    O = _oitools_gui()
    e = _filter_root(sh.session)
    (e === nothing || O === nothing) && return ""
    data = e.data
    lines = String[]
    allof(f) = reduce(vcat, (Float64.(vec(getfield(d, f))) for d in data))
    lam, bl, mjd = allof(:uv_lam), allof(:uv_baseline), allof(:uv_mjd)
    g(v) = Printf.@sprintf("%.4g", v)
    if !isempty(lam)
        push!(lines, join(("extent_wav", g(minimum(lam) * 1e6), g(maximum(lam) * 1e6)), "\t"))
        push!(lines, join(("extent_bl", g(minimum(bl) / 1e6), g(maximum(bl) / 1e6)), "\t"))
        blm = bl .* lam
        push!(lines, join(("extent_blm", g(minimum(blm)), g(maximum(blm))), "\t"))
    end
    ns = O._nights(mjd)
    for i in eachindex(ns)
        a, b = O._night_window(ns, i)
        push!(lines, string("night\t", Printf.@sprintf("%.5f:%.5f", a, b), "\t",
                            Dates.Date(1858, 11, 17) + Dates.Day(floor(Int, ns[i].lo)),
                            " UT  ·  ", ns[i].n, " uv points"))
    end
    hdrs = [h for h in (try OITOOLS.load_oifits(f) catch; nothing end for f in e.paths)
            if h !== nothing]
    for n in unique(String[i.insname for h in hdrs for i in h.instr]); push!(lines, "insname\t" * n); end
    for n in unique(String[t.target for h in hdrs for t in h.target]); push!(lines, "target\t" * n);  end
    for (o, n) in (("v2", :nv2), ("t3amp", :nt3amp), ("t3phi", :nt3phi), ("visamp", :nvisamp),
                   ("visphi", :nvisphi), ("flux", :nflux))
        any(d -> getfield(d, n) > 0, data) && push!(lines, "has\t" * o)
    end
    # Station names per epoch: the station tables of different files need not agree.
    names = Dict{Int,String}()
    stations = Set{String}(); pairs = Set{Vector{String}}(); tris = Set{Vector{String}}()
    for d in data
        nm(k) = (j = findfirst(==(k), d.sta_index); j === nothing ? string(k) : String(d.sta_name[j]))
        for c in eachcol(d.v2_sta_index);  p = sort(nm.(c)); push!(pairs, p); union!(stations, p); end
        for c in eachcol(d.vis_sta_index); p = sort(nm.(c)); push!(pairs, p); union!(stations, p); end
        for c in eachcol(d.t3_sta_index);  p = sort(nm.(c)); push!(tris, p);  union!(stations, p); end
    end
    for n in sort(collect(stations)); push!(lines, "station\t" * n); end
    for p in sort(collect(pairs)); push!(lines, "baseline\t" * join(p, "-")); end
    for t in sort(collect(tris));  push!(lines, "triangle\t" * join(t, "-")); end
    return join(lines, "\n")
end

# The parsed settings of the last preview, kept so the points it would remove are drawn again
# when the epoch or the plotted observable changes. `nothing` with no preview.
const FILTER_PREVIEW = Ref{Any}(nothing)

"""
    shell_filter_preview(spec) -> String

What the Filter panel's settings would remove, over every epoch, without applying them:
"Filtered 37/494 V² · 12/260 T3". The points the current epoch would lose are drawn hollow
and grey on the observable plot.

Instrument, target, table, flag and averaging choices are read options, applied when the files
are read again on Apply, so they are not counted here.
"""
function shell_filter_preview(spec)
    sh = _sh()
    O = _oitools_gui()
    O === nothing && return "! filtering needs OITOOLS' GUI extension"
    root = _filter_root(sh.session)
    root === nothing && return ""
    sp = _split_hand(_qmlstr(spec))
    k = sp isa String ? sp : O.filter_spec_kwargs(sp[1])
    if k isa String
        FILTER_PREVIEW[] = nothing
        _draw_filter_preview!(sh)
        return k
    end
    hand = sp[2]
    idxs = try
        [OITOOLS.set_data_filter(d; k.filter..., get(hand, i, Dict{Symbol,Any}())...)
         for (i, d) in enumerate(root.data)]
    catch err
        FILTER_PREVIEW[] = nothing
        _draw_filter_preview!(sh)
        return "! " * sprint(showerror, err)
    end
    parts = String[]
    for (label, field, slot) in (("V²", :v2, 3), ("T3", :t3phi, 4), ("VIS", :visamp, 2), ("flux", :flux, 5))
        n = sum(length(getfield(d, field)) for d in root.data)
        n > 0 || continue
        push!(parts, string(sum(length(ix[slot]) for ix in idxs), "/", n, " ", label))
    end
    # An epoch whose every point would go (slots 2–5: VIS, V², T3, flux) is removed on Apply.
    emptied(d, ix) = all(length(unique(ix[s])) >= length(getfield(d, f))
                         for (s, f) in ((2, :visamp), (3, :v2), (4, :t3phi), (5, :flux)))
    gone = count(((d, ix),) -> emptied(d, ix), zip(root.data, idxs))
    msg = "Filtered " * join(parts, " · ")
    gone > 0 && (msg *= " · $gone of $(length(root.data)) epoch(s) emptied and removed")
    haskey(k.read, :binning) && (msg *= " · read again averaging every $(k.read[:binning]) channels")
    isempty(k.drop) || (msg *= " · without " * join(string.(k.drop), ", "))
    FILTER_PREVIEW[] = (root = root.name, filter = k.filter, idxs = idxs, mjd = root.mjd)
    _draw_filter_preview!(sh)
    return msg
end

"""Remove the filter preview from the plot (the panel was reset)."""
function shell_filter_preview_clear()
    FILTER_PREVIEW[] = nothing
    _draw_filter_preview!(_sh())
    return ""
end

# The points the previewed filter would remove from the epoch on screen, on the observable plot,
# when that plot shows the dataset being filtered.
function _draw_filter_preview!(sh::ShellState)
    c = sh.obs
    (c === nothing || !hasproperty(c, :prevpoints)) && return nothing
    p = FILTER_PREVIEW[]
    e = current_dataset(sh.session)
    O = _oitools_gui()
    if p === nothing || e === nothing || e.name != p.root || O === nothing
        c.prevpoints[] = Makie.Point2f[]
        return nothing
    end
    i = clamp(sh.session.current_epoch, 1, length(e.data))
    d, idx = e.data[i], p.idxs[i]
    kind = sh.obskind[]
    pts = Makie.Point2f[]
    if kind === :uv
        keep(n, bad) = setdiff(1:n, bad)
        kept = Set{Int}()
        union!(kept, d.indx_v2[keep(length(d.indx_v2), idx[3])])
        t3k = keep(length(d.indx_t3_1), idx[4])
        union!(kept, d.indx_t3_1[t3k], d.indx_t3_2[t3k], d.indx_t3_3[t3k])
        union!(kept, d.indx_vis[keep(length(d.indx_vis), idx[2])])
        gone = union(setdiff(1:d.nuv, kept), idx[1])
        u = Float64.(d.uv[1, gone]) ./ 1e6; v = Float64.(d.uv[2, gone]) ./ 1e6
        pts = vcat(Makie.Point2f.(u, v), Makie.Point2f.(-u, -v))
    elseif haskey(O.OBS_SPECS, kind)
        spec = O.OBS_SPECS[kind]
        slot = spec.y === :v2 ? 3 : spec.y in (:t3phi, :t3amp) ? 4 : spec.y === :flux ? 5 : 2
        x = Float64.(getfield(d, spec.x)) .* spec.xscale
        y = Float64.(getfield(d, spec.y))
        r = sort(idx[slot])
        pts = Makie.Point2f.(x[r], y[r])
        _obs_logscale(sh) && filter!(q -> isfinite(q[2]) && q[2] > 0, pts)
    end
    c.prevpoints[] = pts
    return nothing
end

"""
    shell_apply_filter(spec) -> String

Cut the dataset on screen — or, when that is itself a filtered dataset, the one it came from —
with the Filter panel's settings, into a new dataset named "… › filtered", and select it. The
source is kept, so selecting it again is the undo. Logged, so the exported script makes the
same cut.
"""
function shell_apply_filter(spec)
    sh = _sh()
    O = _oitools_gui()
    O === nothing && return "! filtering needs OITOOLS' GUI extension"
    root = _filter_root(sh.session)
    root === nothing && return "! no dataset loaded"
    sp = _qmlstr(spec)
    hs = _split_hand(sp)
    hs isa String && return hs
    k = O.filter_spec_kwargs(hs[1])
    k isa String && return k
    hand = hs[2]
    (isempty(k.read) && isempty(k.filter) && isempty(k.prep) && isempty(k.drop) &&
     isempty(hand)) && return "! nothing to filter: every setting is at its default"
    lost = Ref(0); used = Ref{Any}(nothing)
    e = try
        filter_dataset!(sh.session, root; read = k.read, prep = k.prep, drop = k.drop,
                        hand = hand, hand_lost = lost, hand_used = used, spec = sp, k.filter...)
    catch err
        sh.status = "! could not filter $(root.name): " * sprint(showerror, err)
        console!(sh, sh.status)
        return sh.status
    end
    _log_filter!(sh.session, root, e, k, used[])
    FILTER_PREVIEW[] = nothing
    shell_clear_selection()
    sh.status = "$(e.name): $(length(e.data)) epoch(s), $(sum(e.nv2)) V², $(sum(e.nt3)) T3"
    lost[] > 0 && (sh.status *= " — $(lost[]) hand-picked point(s) not found after reading " *
                                 "the files again, and kept")
    console!(sh, "filter: " * sh.status)
    refresh_data_tab!(sh)
    return sh.status
end

"""The Filter panel's settings that made the current dataset; empty for a dataset as read."""
function shell_filter_spec()
    e = current_dataset(_sh().session)
    return e === nothing ? "" : e.filterspec
end

# ── points by hand: identify, select, add to the filter ──────────────────────────────────
#
# OITOOLS' Explore tools, on the epoch on screen. A click names the nearest point; a box or a
# lasso drawn with a selection tool chosen selects points, and "Add to filter" (in the Filter
# panel) turns them into hand-picked removals — with the epoch they belong to, which is what a
# ROTIR dataset adds to OITOOLS' single one.

# `:none`, `:box` or `:lasso`; and the selection: the dataset, its epoch, the plot kind it was
# made on, the observable slot (`selection_indices`: 1 uv, 2 VIS, 3 V², 4 T3, 5 flux) and the
# indices into it.
const OBS_SEL_TOOL  = Ref(:none)
const OBS_SELECTION = Ref{Any}(nothing)

const _SLOT_NAMES = ("uv point", "VIS point", "V² point", "T3 point", "flux point")

"""
    shell_set_select_tool(mode) -> String

`"none"`, `"box"` or `"lasso"`. With a tool chosen, a left drag on the observable plot selects
instead of zooming; with none, its zoom, pan, identification and reset are untouched.
"""
function shell_set_select_tool(mode)
    sh = _sh()
    t = Symbol(_qmlstr(mode))
    t in (:none, :box, :lasso) || return "! the selection tools are none, box and lasso"
    OBS_SEL_TOOL[] = t
    if sh.obs !== nothing
        t === :none ? Makie.activate_interaction!(sh.obs.axis, :rectanglezoom) :
                      Makie.deactivate_interaction!(sh.obs.axis, :rectanglezoom)
    end
    return ""
end

"""
    select_obs_polygon!(sh, poly) -> Int

Select the points of the current epoch's plot inside `poly` (data coordinates; a box is its four
corners), ring them, and keep them as the selection. What the mouse handler calls on release,
and what a test calls without a mouse. Nothing is selected on the per-group panels or on a
residual plot, whose y is not the data's. Returns how many were selected.
"""
function select_obs_polygon!(sh::ShellState, poly::AbstractVector)
    O = _oitools_gui()
    e = current_dataset(sh.session)
    (O === nothing || e === nothing || sh.obspanels[] || length(poly) < 3 ||
     sh.obscompare[] in (:residmodel, :residimage)) && return 0
    i = clamp(sh.session.current_epoch, 1, length(e.data))
    slot, idx, pts = O.selection_indices(e.data[i], sh.obskind[], poly)
    OBS_SELECTION[] = isempty(idx) ? nothing :
        (; dataset = e.name, epoch = i, kind = sh.obskind[], slot, idx)
    sh.obs === nothing || (sh.obs.selpoints[] = pts)
    return length(idx)
end

"""What is selected: "12 V² points of epoch 2 on lam And", or empty."""
function shell_selection_text()
    s = OBS_SELECTION[]
    s === nothing && return ""
    return string(length(s.idx), " ", _SLOT_NAMES[s.slot], length(s.idx) == 1 ? "" : "s",
                  " of epoch ", s.epoch, " on ", s.dataset)
end

"""Forget the selection and its rings."""
function shell_clear_selection()
    sh = _sh()
    OBS_SELECTION[] = nothing
    if sh.obs !== nothing
        sh.obs.selpoints[] = Makie.Point2f[]
        sh.obs.seloutline[] = Makie.Point2f[]
    end
    return ""
end

"""
    shell_selection_spec() -> String

The selection as a Filter-panel line — `remove_v2\\t2:3,2:17` — in epochs and indices of the
unfiltered dataset the panel filters, so "Add to filter" adds them to its hand-picked points. A
selection made on a filtered dataset is matched to the source by MJD (the epoch) and by uv,
wavelength and stations (the point). `""` with nothing selected.
"""
function shell_selection_spec()
    sh = _sh()
    O = _oitools_gui()
    s = OBS_SELECTION[]
    (s === nothing || O === nothing) && return ""
    k = findfirst(x -> x.name == s.dataset, sh.session.datasets)
    k === nothing && return "! the dataset the selection was made on is gone"
    e = sh.session.datasets[k]
    root = e
    for _ in 1:length(sh.session.datasets)
        isempty(root.source) && break
        j = findfirst(x -> x.name == root.source, sh.session.datasets)
        j === nothing && break
        root = sh.session.datasets[j]
    end
    r = findfirst(m -> isapprox(m, e.mjd[s.epoch]; atol = 1e-6), root.mjd)
    r === nothing && return "! epoch $(s.epoch) is not in $(root.name)"
    idx = if root === e
        s.idx
    else
        where = Dict(q => n for (n, q) in enumerate(O._point_keys(root.data[r], s.slot)))
        sk = O._point_keys(e.data[s.epoch], s.slot)
        sort!([where[sk[i]] for i in s.idx if haskey(where, sk[i])])
    end
    return string(HAND_KEYS[s.slot].first, "\t", join(("$(r):$(i)" for i in idx), ","))
end

"""
    shell_obs_pick() -> String

The plotted point nearest the cursor, described as OITOOLS describes it, or `""` when the cursor
is not near one. Distances are in fractions of the axes, so a baseline in Mλ and a V² weigh the
same; 2 % of the diagonal is roughly a finger's width. Searched over the points held, because
GPU picking needs GLMakie's own render loop, which QMLMakie does not run.
"""
function shell_obs_pick(; tol::Real = 0.02)
    sh = _sh()
    c = sh.obs
    (c === nothing || isempty(sh.obsinfo) || sh.obspanels[]) && return ""
    pts = c.points[]
    isempty(pts) && return ""
    try
        w = Makie.mouseposition(c.axis.scene)
        lims = c.axis.finallimits[]
        x0, y0 = Float64(lims.origin[1]), Float64(lims.origin[2])
        wx, wy = Float64(lims.widths[1]), Float64(lims.widths[2])
        (wx > 0 && wy > 0) || return ""
        # The limits are in data units; on a log axis the distance that matches what is on
        # screen is in log10, across the decades the axis spans.
        logy = c.axis.yscale[] === log10
        ty(y) = logy ? log10(max(Float64(y), floatmin())) : Float64(y)
        hy = logy ? ty(y0 + wy) - ty(y0) : wy
        hy > 0 || return ""
        cx, cy = Float64(w[1]), ty(w[2])
        best, bestd = 0, Inf
        for (i, q) in enumerate(pts)
            dx = (Float64(q[1]) - cx) / wx
            dy = (ty(q[2]) - cy) / hy
            d = dx * dx + dy * dy
            d < bestd && ((bestd, best) = (d, i))
        end
        (best < 1 || best > length(sh.obsinfo)) && return ""
        sqrt(bestd) <= tol || return ""
        d = current_dataset(sh.session)
        ep = d === nothing ? "" : "epoch $(clamp(sh.session.current_epoch, 1, length(d.data))) · "
        return ep * sh.obsinfo[best]
    catch
        return ""
    end
end

"""Last identified point. QML polls this rather than asking for a pick itself."""
shell_pick_text() = _sh().picktext

"""
    install_obs_interactions!(sh)

Left-click identification and the box and lasso, on the observable plot, read from Makie's own
event stream — `MakieArea`'s MouseArea takes an exclusive grab on press, so a QML handler would
stop firing after the first click (OITOOLS' `install_interactions!` says more). Right-click
reset is the axis's own `:limitreset` (`build_obs_canvas`). Nothing here creates a plot.
"""
function install_obs_interactions!(sh::ShellState)
    c = sh.obs
    c === nothing && return nothing
    scene  = c.axis.scene
    events = Makie.events(scene)
    drawing = Ref(false)
    Makie.on(events.mouseposition, priority = 100) do _
        drawing[] || return Makie.Consume(false)
        p = Makie.Point2f(Makie.mouseposition(scene))
        o = c.seloutline[]
        if OBS_SEL_TOOL[] === :box && !isempty(o)
            a = o[1]
            c.seloutline[] = [a, Makie.Point2f(p[1], a[2]), p, Makie.Point2f(a[1], p[2]), a]
        else
            c.seloutline[] = vcat(o, [p])
        end
        return Makie.Consume(true)
    end
    Makie.on(events.mousebutton, priority = 100) do ev
        if OBS_SEL_TOOL[] !== :none && ev.button == Makie.Mouse.left &&
           (ev.action == Makie.Mouse.press ? Makie.is_mouseinside(scene) : drawing[])
            if ev.action == Makie.Mouse.press
                drawing[] = true
                c.seloutline[] = [Makie.Point2f(Makie.mouseposition(scene))]
            else
                drawing[] = false
                poly = c.seloutline[]
                c.seloutline[] = Makie.Point2f[]
                select_obs_polygon!(sh, poly)
            end
            return Makie.Consume(true)
        end
        if ev.button == Makie.Mouse.left && ev.action == Makie.Mouse.press &&
           Makie.is_mouseinside(scene)
            sh.picktext = shell_obs_pick()
        end
        return Makie.Consume(false)      # zoom and pan still need these
    end
    return nothing
end

# ── export ───────────────────────────────────────────────────────────────────────────────

"""
    shell_export_oifits(path) -> String

Write the current dataset — every epoch, filtered as it is — to one OIFITS file through OITOOLS'
`write_oifits`. A file the dataset was read from is refused: writing over it would replace the
data with a copy of itself less whatever was filtered out. Logged.
"""
function shell_export_oifits(path)
    sh = _sh()
    e = current_dataset(sh.session)
    e === nothing && return "! no dataset loaded"
    p = _qmlstr(path)
    isempty(p) && return "! no file name"
    any(f -> isfile(f) && isfile(p) && samefile(p, f), e.paths) &&
        return "! that is a file the dataset was read from; choose another name"
    try
        OITOOLS.write_oifits(p, e.data; overwrite = true)
    catch err
        msg = "! could not write OIFITS: " * sprint(showerror, err)
        console!(sh, msg)
        return msg
    end
    log!(sh.session, "write_oifits($(_literal(p)), data; overwrite = true)";
         note = "export $(e.name)")
    sh.status = "saved $(basename(p)): $(length(e.data)) epoch(s) of $(e.name)"
    console!(sh, sh.status)
    return sh.status
end
