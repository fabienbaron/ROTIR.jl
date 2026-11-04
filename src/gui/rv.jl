# The RV perspective: a table of radial velocities, drawn against the orbit of the Orbit tab.
#
# The table is DATA and the velocity parameters belong to the ORBIT. K1, K2, gamma, the
# per-instrument offsets and the two jitters are stored in `sh.orbit.params` like the
# elements, so the RV tab and the Orbit tab edit one set of numbers, an orbit TOML saves them,
# and a fit writes them back in one place. This tab shows them and edits them; the fit runs in
# the Orbit tab.
#
# The model curve comes from the same `apply_orbit_ties(sh.orbit)` every orbit picture uses, so
# editing an element in the Orbit tab moves the curve here (`refresh_orbit!` calls
# `refresh_rv!`).

"Colours of the two stars, as in OITOOLS' orbit view: a (primary) blue, b (secondary) orange."
const RV_COLOR_A = Makie.RGBAf(0.0, 0.447, 0.698, 1.0)
const RV_COLOR_B = Makie.RGBAf(0.902, 0.624, 0.0, 1.0)
"Marker per instrument, in table order."
const RV_MARKERS = (:circle, :rect, :utriangle, :diamond, :dtriangle, :star5, :hexagon, :cross)

# PLACEHOLDERS, never empty vectors, for the plots with per-point colours or markers: GLMakie
# cannot build a render object for a scatter or a linesegments whose colour (or marker) vector
# is empty ("Failed to resolve gl_renderobject"), QMLMakie reports that only as "exception in
# render", and the area then stays blank even after real data arrive. One NaN point (two for a
# segment) draws nothing and is ignored by the axis limits.
_rv_nopt()  = [Makie.Point2f(NaN, NaN)]
_rv_noseg() = [Makie.Point2f(NaN, NaN), Makie.Point2f(NaN, NaN)]
const RV_NOCOLOR = Makie.RGBAf(0, 0, 0, 0)

const RV_PARAM_DEFAULT_BOUNDS = Dict(:plx => (1e-4, 1e4), :q => (0.01, 100.0),
                                     :gamma => (-1000.0, 1000.0),
                                     :jit_a => (0.0, 50.0), :jit_b => (0.0, 50.0))

"""
    RVCanvas

Velocities against phase (or time) above, residuals O − C below, sharing the x axis.

Every plot is created by [`build_rv_canvas`](@ref) and only reassigned afterwards — the
plot-once rule of livecanvas.jl. Points of both stars and every instrument go through ONE
scatter with per-point colours and markers, and error bars through one `linesegments!`, so a
table of any size and any number of instruments needs no new plots.
"""
struct RVCanvas
    figure::Any
    axis::Any
    resaxis::Any
    curve_a::Makie.Observable{Vector{Makie.Point2f}}
    curve_b::Makie.Observable{Vector{Makie.Point2f}}
    gammaline::Makie.Observable{Vector{Makie.Point2f}}
    pts::Makie.Observable{Vector{Makie.Point2f}}
    ptcolors::Makie.Observable{Vector{Makie.RGBAf}}
    ptmarkers::Makie.Observable{Vector{Symbol}}
    bars::Makie.Observable{Vector{Makie.Point2f}}
    barcolors::Makie.Observable{Vector{Makie.RGBAf}}
    res::Makie.Observable{Vector{Makie.Point2f}}
    rescolors::Makie.Observable{Vector{Makie.RGBAf}}
    resmarkers::Makie.Observable{Vector{Symbol}}
    resbars::Makie.Observable{Vector{Makie.Point2f}}
    resbarcolors::Makie.Observable{Vector{Makie.RGBAf}}
    zeroline::Makie.Observable{Vector{Makie.Point2f}}
    cursor::Makie.Observable{Vector{Float64}}
    message::Makie.Observable{String}
    xlabel::Makie.Observable{String}
    reslabel::Makie.Observable{String}
end

# Every string this canvas can show, so its glyphs can be rasterised before the window exists.
const RV_CANVAS_TEXT = ("orbital phase", "orbital phase, (t − T0)/P", "MJD", "RV (km/s)",
                        "O − C (km/s)", "O − C (σ)", "no RV table — Open RV table…",
                        "the period must be positive", "the table has no velocities",
                        "could not draw the velocities — check the elements",
                        "−0123456789.e+")

"""
Put every glyph the RV canvas can show into Makie's texture atlas, in the axis font, NOW.

A glyph first met after the window is up is rasterised and uploaded to the GPU at render
time, and once a popup or a picker has taken the GL context that upload fails ("GLX: Failed
to make context current") and takes the process down — measured in the click test, on the
x label changing to "orbital phase, (t − T0)/P" after a run of popups. The labels change
with the data and the view, so every string they can take is warmed here, at build time.
"""
function _prewarm_rv_glyphs!()
    atlas = Makie.get_texture_atlas()
    for name in ("TeX Gyre Heros Makie", "TeX Gyre Heros Makie Bold")
        font = try
            Makie.to_font(name)
        catch
            continue
        end
        for str in RV_CANVAS_TEXT, c in str
            try
                Makie.insert_glyph!(atlas, c, font)
            catch
            end
        end
    end
    # The instrument markers too: a marker shape other than :circle is rasterised into the same
    # atlas, and a table with several instruments shows them only after the window is up.
    mm = Makie.default_marker_map()
    for m in RV_MARKERS
        haskey(mm, m) || continue
        try
            Makie.insert_glyph!(atlas, mm[m])
        catch
        end
    end
    return nothing
end

function build_rv_canvas(fig)
    _prewarm_rv_glyphs!()
    xlab = Makie.Observable("orbital phase")
    rlab = Makie.Observable("O − C (km/s)")
    ax  = Makie.Axis(fig[1, 1]; ylabel = "RV (km/s)", xticklabelsvisible = false)
    rax = Makie.Axis(fig[2, 1]; xlabel = xlab, ylabel = rlab)
    Makie.linkxaxes!(ax, rax)
    Makie.rowsize!(fig.layout, 2, Makie.Relative(0.28))
    Makie.rowgap!(fig.layout, 6)
    for a in (ax, rax); a.xgridvisible = false; a.ygridvisible = false; end

    gam = Makie.Observable(Makie.Point2f[])
    Makie.lines!(ax, gam; color = (:grey50, 0.8), linestyle = :dash, linewidth = 1)
    ca  = Makie.Observable(Makie.Point2f[])
    cb  = Makie.Observable(Makie.Point2f[])
    Makie.lines!(ax, ca; color = RV_COLOR_A, linewidth = 1.4)
    Makie.lines!(ax, cb; color = RV_COLOR_B, linewidth = 1.4)
    bars = Makie.Observable(_rv_noseg()); barc = Makie.Observable([RV_NOCOLOR, RV_NOCOLOR])
    Makie.linesegments!(ax, bars; color = barc, linewidth = 1)
    pts  = Makie.Observable(_rv_nopt()); ptc = Makie.Observable([RV_NOCOLOR])
    ptm  = Makie.Observable([:circle])
    Makie.scatter!(ax, pts; color = ptc, marker = ptm, markersize = 8,
                   strokecolor = :black, strokewidth = 0.5)

    zl   = Makie.Observable(Makie.Point2f[])
    Makie.lines!(rax, zl; color = (:grey50, 0.8), linestyle = :dash, linewidth = 1)
    rbar = Makie.Observable(_rv_noseg()); rbarc = Makie.Observable([RV_NOCOLOR, RV_NOCOLOR])
    Makie.linesegments!(rax, rbar; color = rbarc, linewidth = 1)
    res  = Makie.Observable(_rv_nopt()); resc = Makie.Observable([RV_NOCOLOR])
    resm = Makie.Observable([:circle])
    Makie.scatter!(rax, res; color = resc, marker = resm, markersize = 8,
                   strokecolor = :black, strokewidth = 0.5)

    # The Orbit tab's time cursor, in both panels.
    # Excluded from the autolimits: it marks a time, it should not set the x range — and an
    # EMPTY `vlines!` taking part in the limits fails outright (`extrema` of nothing).
    cur = Makie.Observable(Float64[])
    for a in (ax, rax)
        Makie.vlines!(a, cur; color = (:dodgerblue, 0.8), linewidth = 1.2,
                      xautolimits = false, yautolimits = false)
    end

    msg = Makie.Observable("")
    Makie.text!(ax, [Makie.Point2f(0.5, 0)]; text = msg, fontsize = 15, space = :relative,
                color = Makie.RGBAf(0.45, 0.5, 0.55, 1), align = (:center, :bottom))

    c = RVCanvas(fig, ax, rax, ca, cb, gam, pts, ptc, ptm, bars, barc, res, resc, resm,
                 rbar, rbarc, zl, cur, msg, xlab, rlab)
    _rv_clear!(c, "no RV table — Open RV table…")
    return c
end

function _rv_clear!(c::RVCanvas, message::AbstractString)
    for o in (c.curve_a, c.curve_b, c.gammaline, c.zeroline)
        o[] = Makie.Point2f[]
    end
    c.ptcolors[] = [RV_NOCOLOR]; c.ptmarkers[] = [:circle]; c.pts[] = _rv_nopt()
    c.rescolors[] = [RV_NOCOLOR]; c.resmarkers[] = [:circle]; c.res[] = _rv_nopt()
    c.barcolors[] = [RV_NOCOLOR, RV_NOCOLOR]; c.bars[] = _rv_noseg()
    c.resbarcolors[] = [RV_NOCOLOR, RV_NOCOLOR]; c.resbars[] = _rv_noseg()
    c.cursor[] = Float64[]
    c.message[] = String(message)
    return c
end

# ── the velocity parameters, stored in the orbit ─────────────────────────────────────────

"""
    rv_param_names(o, rv) -> Vector{Symbol}

The velocity parameters an RV table brings into the orbit: `plx, q, gamma`, one
`dgamma_<instrument>` per instrument after the first, `jit_a, jit_b`. Empty without a table.
"""
function rv_param_names(rv)
    rv === nothing && return Symbol[]
    return Symbol[:plx, :q, :gamma, ROTIR.rv_offset_name.(rv_instruments(rv)[2:end])...,
                  :jit_a, :jit_b]
end

"""
    seed_rv_params!(o, rv) -> o

Give the orbit the velocity parameters `rv` needs, keeping any it already has: plx, q and
gamma are free, plx and q starting from the amplitudes [`rv_start`](@ref) reads off the data,
at the orbit's current a, i, P and e; offsets start at 0 and are free; jitters
start at 0 and are never free (set in the RV tab). Offsets of instruments the table does not
have are removed.
"""
function seed_rv_params!(o, rv)
    st = rv_start(rv)
    v = apply_orbit_ties(o)
    p0, q0 = rv_physical(v[:a], v[:i], v[:P], v[:e], st.K1, st.K2)
    (isfinite(p0) && p0 > 0) || (p0 = 1.0)
    (isfinite(q0) && q0 > 0) || (q0 = 1.0)
    for (k, v) in ((:plx, p0), (:q, clamp(q0, 0.01, 100.0)), (:gamma, st.gamma))
        haskey(o.params, k) || (o.params[k] = round(v; digits = 3); push!(o.free, k))
    end
    for k in (:jit_a, :jit_b)
        haskey(o.params, k) || (o.params[k] = 0.0)
        delete!(o.free, k); delete!(o.ties, k)
    end
    for (k, b) in RV_PARAM_DEFAULT_BOUNDS; haskey(o.bounds, k) || (o.bounds[k] = b); end
    want = Set(ROTIR.rv_offset_name.(rv_instruments(rv)[2:end]))
    _drop_offsets!(o; keep = want)
    for k in want
        haskey(o.params, k) || (o.params[k] = 0.0; push!(o.free, k))
        haskey(o.bounds, k) || (o.bounds[k] = (-50.0, 50.0))
    end
    return o
end

function _drop_offsets!(o; keep = Set{Symbol}())
    for k in collect(keys(o.params))
        startswith(String(k), "dgamma_") && !(k in keep) || continue
        delete!(o.params, k); delete!(o.free, k); delete!(o.bounds, k); delete!(o.ties, k)
    end
    return o
end

"""
(K1, K2, gamma, offsets, jitter) from the orbit's resolved values, for `rv_model` & co. The
amplitudes follow from plx and q at the orbit's a, i, P and e (`rv_amplitudes`).
"""
function _rv_values(v::AbstractDict, rv)
    ins = rv_instruments(rv)
    off = Dict{String,Float64}(i => Float64(get(v, ROTIR.rv_offset_name(i), 0.0))
                               for i in ins[2:end])
    K1, K2 = rv_amplitudes(v[:a], v[:i], v[:P], v[:e], get(v, :plx, 1.0), get(v, :q, 1.0))
    return (K1 = K1, K2 = K2, gamma = get(v, :gamma, 0.0),
            offsets = off, jitter = (a = get(v, :jit_a, 0.0), b = get(v, :jit_b, 0.0)))
end

"The elements as `binary_RV` reads them, from the orbit's resolved values."
_rv_bp(v::AbstractDict) = (ω = v[:omega], e = v[:e], P = v[:P], T0 = v[:T0], dP = v[:dP],
                           dω = v[:domega])

# ── drawing ──────────────────────────────────────────────────────────────────────────────

_phase(t, T0, P) = mod((t - T0) / P, 1.0)

"""
    show_rv!(c, rv, v; xaxis = :phase, normalized = false, cursor = nothing)

Draw `rv` against the orbit whose resolved values are `v` (`apply_orbit_ties(orbit)`).

`xaxis = :phase` folds every point on (t − T0)/P mod 1 and draws the model over one cycle at
the cursor epoch — or, without a cursor, at the cycle nearest the mean time of the table.
With apsidal motion or a period derivative a folded curve is exact only at that epoch;
`:time` draws the model through every observation. Points are shifted by their
instrument's offset, so all instruments sit on the same curve. Error bars are the tabulated
σ inflated by the star's jitter, √(σ² + s²); `normalized` plots the residuals in those units.
"""
function show_rv!(c::RVCanvas, rv, v::AbstractDict; xaxis::Symbol = :phase,
                  normalized::Bool = false, cursor = nothing)
    p = _rv_values(v, rv)
    bp = _rv_bp(v)
    P, T0 = Float64(v[:P]), Float64(v[:T0])
    P > 0 || return _rv_clear!(c, "the period must be positive")
    length(rv) > 0 || return _rv_clear!(c, "the table has no velocities")
    model = rv_model(bp, rv; K1 = p.K1, K2 = p.K2, gamma = p.gamma, offsets = p.offsets)
    n = length(rv)
    σe = [sqrt(rv.sigma[k]^2 + (rv.star[k] === :a ? p.jitter.a : p.jitter.b)^2) for k in 1:n]
    off = [get(p.offsets, rv.instrument[k], 0.0) for k in 1:n]
    ins = rv_instruments(rv)
    mk  = [RV_MARKERS[mod1(findfirst(==(rv.instrument[k]), ins), length(RV_MARKERS))] for k in 1:n]
    col = [rv.star[k] === :a ? RV_COLOR_A : RV_COLOR_B for k in 1:n]

    mjd = rv.jd .- ROTIR.MJD_TO_JD
    x = xaxis === :phase ? _phase.(rv.jd, T0, P) : mjd
    y = rv.rv .- off                                  # each instrument on the common curve
    r = normalized ? (rv.rv .- model) ./ σe : rv.rv .- model
    re = normalized ? ones(n) : σe

    # The model curves.
    if xaxis === :phase
        tref = cursor === nothing ? sum(rv.jd) / n : Float64(cursor)
        cyc = floor((tref - T0) / P)
        φ = collect(range(0.0, 1.0; length = 400))
        tc = T0 .+ (cyc .+ φ) .* P
        xc = φ
    else
        lo, hi = extrema(rv.jd); pad = 0.03 * max(hi - lo, P)
        npt = clamp(ceil(Int, 40 * (hi - lo + 2pad) / P), 400, 20_000)
        tc = collect(range(lo - pad, hi + pad; length = npt))
        xc = tc .- ROTIR.MJD_TO_JD
    end
    va, vb = binary_RV(bp, tc; K1 = p.K1, K2 = p.K2, γ = p.gamma)
    c.curve_a[] = Makie.Point2f.(xc, va)
    c.curve_b[] = Makie.Point2f.(xc, vb)
    c.gammaline[] = [Makie.Point2f(xc[1], p.gamma), Makie.Point2f(xc[end], p.gamma)]
    c.zeroline[] = [Makie.Point2f(xc[1], 0), Makie.Point2f(xc[end], 0)]

    c.ptcolors[] = col; c.ptmarkers[] = mk
    c.pts[] = Makie.Point2f.(x, y)
    c.barcolors[] = repeat(col; inner = 2)
    c.bars[] = vec(permutedims(hcat(Makie.Point2f.(x, y .- σe), Makie.Point2f.(x, y .+ σe))))
    c.rescolors[] = col; c.resmarkers[] = mk
    c.res[] = Makie.Point2f.(x, r)
    c.resbarcolors[] = repeat(col; inner = 2)
    c.resbars[] = vec(permutedims(hcat(Makie.Point2f.(x, r .- re), Makie.Point2f.(x, r .+ re))))

    c.cursor[] = cursor === nothing ? Float64[] :
                 [xaxis === :phase ? _phase(Float64(cursor), T0, P) :
                  Float64(cursor) - ROTIR.MJD_TO_JD]
    c.xlabel[] = xaxis === :phase ? "orbital phase, (t − T0)/P" : "MJD"
    c.reslabel[] = normalized ? "O − C (σ)" : "O − C (km/s)"
    c.message[] = ""
    Makie.reset_limits!(c.axis); Makie.reset_limits!(c.resaxis)
    return c
end

"""
    refresh_rv!(sh) -> status

Redraw the RV tab from the current table and orbit. Called by `refresh_orbit!`, so every
edit of the orbit reaches this picture too.
"""
function refresh_rv!(sh::ShellState)
    sh.rvcanvas === nothing && return sh.status
    sh.rv === nothing && (_rv_clear!(sh.rvcanvas, "no RV table — Open RV table…"); return sh.status)
    try
        show_rv!(sh.rvcanvas, sh.rv, apply_orbit_ties(sh.orbit); xaxis = sh.rvxaxis,
                 normalized = sh.rvnormalized, cursor = cursor_time(sh))
    catch err
        _rv_clear!(sh.rvcanvas, "could not draw the velocities — check the elements")
        console!(sh, "RV draw failed: $(sprint(showerror, err))")
    end
    return sh.status
end

# ── the callbacks ────────────────────────────────────────────────────────────────────────

"""
    shell_rv_open(path) -> String

Read an RV table (CSV `MJD,RV,sigma,star[,instrument]`, the format OITOOLS reads) and make it
the current one. The orbit gains the velocity parameters it needs (see
[`seed_rv_params!`](@ref)).
"""
function shell_rv_open(path)
    sh = _sh()
    p = strip(replace(_qmlstr(path), r"^file://" => ""))
    isempty(p) && return "no file given"
    rv = try
        read_rv(p)
    catch err
        msg = "could not read $(basename(p)): $(sprint(showerror, err))"
        console!(sh, msg)
        return msg
    end
    sh.rv = rv
    seed_rv_params!(sh.orbit, rv)
    isempty(sh.session.log) || (sh.session.log[end].binding == "rv" && pop!(sh.session.log))
    log!(sh.session, "rv = read_rv($(_literal(String(p))))"; note = "radial velocities",
         binding = "rv")
    console!(sh, "rv = read_rv($(_literal(String(p))))"; kind = :cmd)
    sh.status = shell_rv_summary()
    console!(sh, sh.status)
    refresh_rv!(sh)
    return sh.status
end

"Forget the RV table. plx, q, gamma and the jitters stay in the orbit; offsets go."
function shell_rv_close()
    sh = _sh()
    sh.rv === nothing && return "no RV table open"
    sh.rv = nothing
    _drop_offsets!(sh.orbit)
    refresh_rv!(sh)
    return "RV table closed"
end

"""
    shell_rv_summary() -> String

`file: N velocities (na a, nb b), MJD lo–hi`, plus the instruments; empty without a table.
"""
function shell_rv_summary()
    rv = _sh().rv
    rv === nothing && return ""
    na = count(==(:a), rv.star)
    lo, hi = extrema(rv.jd) .- ROTIR.MJD_TO_JD
    s = Printf.@sprintf("%s: %d velocities (%d a, %d b), MJD %.1f–%.1f", basename(rv.path),
                        length(rv), na, length(rv) - na, lo, hi)
    ins = filter(!isempty, rv_instruments(rv))
    isempty(ins) || (s *= ", instruments " * join(ins, ", "))
    return s
end

"Rows for plx, q, gamma and the instrument offsets, in the orbit form's columns."
function shell_rv_params()
    sh = _sh()
    sh.rv === nothing && return ""
    ns = filter(n -> !(n in (:jit_a, :jit_b)), rv_param_names(sh.rv))
    rows = split(_orbit_param_rows(sh.orbit, ns), '\n')
    # An offset's description names the instrument as the TABLE spells it — the parameter
    # name has it reduced to an identifier — and the reference it is measured from.
    ins = rv_instruments(sh.rv)
    for i in ins[2:end]
        n = String(ROTIR.rv_offset_name(i))
        k = findfirst(r -> startswith(r, n * "\t"), rows)
        k === nothing && continue
        f = split(rows[k], '\t')
        f[8] = "velocity offset of $(i) relative to $(ins[1])"
        rows[k] = join(f, '\t')
    end
    return join(rows, '\n')
end

"""
Rows for the two jitters, in the orbit form's columns. They are SETTINGS: entered here, used
at their values by the fit, never free — the RV tab draws them with a value field only.
"""
function shell_rv_jitter()
    sh = _sh()
    sh.rv === nothing && return ""
    return _orbit_param_rows(sh.orbit, [:jit_a, :jit_b])
end

"""
    shell_rv_rows() -> String

The table, one row per velocity: `MJD RV σ star instrument model O−C (O−C)/σ'`, where σ' is σ
inflated by the star's jitter.
"""
function shell_rv_rows()
    sh = _sh(); rv = sh.rv
    rv === nothing && return ""
    v = apply_orbit_ties(sh.orbit)
    p = _rv_values(v, rv)
    m = try
        rv_model(_rv_bp(v), rv; K1 = p.K1, K2 = p.K2, gamma = p.gamma, offsets = p.offsets)
    catch
        fill(NaN, length(rv))
    end
    rows = String[]
    for k in 1:length(rv)
        s = sqrt(rv.sigma[k]^2 + (rv.star[k] === :a ? p.jitter.a : p.jitter.b)^2)
        push!(rows, join((Printf.@sprintf("%.5f", rv.jd[k] - ROTIR.MJD_TO_JD),
                          Printf.@sprintf("%.3f", rv.rv[k]), Printf.@sprintf("%.3f", rv.sigma[k]),
                          String(rv.star[k]), rv.instrument[k],
                          Printf.@sprintf("%.3f", m[k]), Printf.@sprintf("%.3f", rv.rv[k] - m[k]),
                          Printf.@sprintf("%.2f", (rv.rv[k] - m[k]) / s)), "\t"))
    end
    return join(rows, "\n")
end

"""
    shell_rv_stats() -> String

Per star and in total: the number of velocities, rms(O − C) in km/s, and χ²/N with the
jitter-inflated errors.
"""
function shell_rv_stats()
    sh = _sh(); rv = sh.rv
    rv === nothing && return ""
    v = apply_orbit_ties(sh.orbit)
    p = _rv_values(v, rv)
    r = try
        rv_residuals(_rv_bp(v), rv; K1 = p.K1, K2 = p.K2, gamma = p.gamma, offsets = p.offsets,
                     jitter = p.jitter)
    catch err
        return "model failed: $(sprint(showerror, err))"
    end
    m = rv_model(_rv_bp(v), rv; K1 = p.K1, K2 = p.K2, gamma = p.gamma, offsets = p.offsets)
    lines = String[]
    # What plx and q mean at this orbit: the amplitudes, the distance, the dynamical masses.
    dm = dynamical_masses(v[:a], v[:P], get(v, :plx, 1.0), get(v, :q, 1.0))
    push!(lines, Printf.@sprintf("K1 %.3f  K2 %.3f km/s   d %.2f pc", p.K1, p.K2, dm.d))
    push!(lines, Printf.@sprintf("M %.4g M☉  (M1 %.4g, M2 %.4g)", dm.M, dm.M1, dm.M2))
    for (s, lab) in ((:a, "a (primary)  "), (:b, "b (secondary)"))
        i = findall(==(s), rv.star)
        isempty(i) && continue
        push!(lines, Printf.@sprintf("%s  N %3d   rms O−C %8.3f km/s   χ²/N %8.3f", lab,
                                     length(i), sqrt(sum(abs2, rv.rv[i] .- m[i]) / length(i)),
                                     sum(abs2, r[i]) / length(i)))
    end
    push!(lines, Printf.@sprintf("total          N %3d                              χ²/N %8.3f",
                                 length(rv), sum(abs2, r) / length(rv)))
    return join(lines, "\n")
end

"`xaxis\tnormalized` — `phase` or `time`, and `1` when residuals are in σ."
shell_rv_view() = (sh = _sh(); string(sh.rvxaxis, "\t", sh.rvnormalized ? "1" : "0"))

"Switch the x axis (`phase`/`time`) and the residual units; redraw."
function shell_rv_set_view(xaxis, normalized)
    sh = _sh()
    x = _qmlstr(xaxis)
    x in ("phase", "time") || return "x axis is phase or time, not $(x)"
    sh.rvxaxis = Symbol(x)
    sh.rvnormalized = _qmlstr(normalized) in ("1", "true")
    refresh_rv!(sh)
    return ""
end
