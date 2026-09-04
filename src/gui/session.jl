# What the window is looking at: datasets, models, reconstructed maps.
#
# One mutable `Session` holds everything the three tabs share, and every tab reads it rather
# than keeping a copy. QML holds NO state of its own: an edit goes to Julia, the table is
# re-read, and what is on screen is therefore always what a script would see. That is what
# makes the command log (src/gui/commandlog.jl) able to claim it reproduces the session.
#
# Nothing here touches Makie or Qt. A test can build a Session, load a dataset, edit a model
# and check the numbers with no display at all.

"""
    DatasetEntry

One OIFITS file, already split into epochs.

ROTIR is multi-epoch by construction — a surface map is constrained by how the star looks at
several rotation phases — so the epoch split happens at load time and everything downstream
sees a vector. `tepochs` is in the same units the model's `rotation_period` is (days), because
that is what `create_star_multiepochs` divides.
"""
mutable struct DatasetEntry
    paths::Vector{String}
    name::String
    data::Vector{Any}            # per-epoch OIdata
    mjd::Vector{Float64}         # absolute mean MJD per epoch, for display
    tepochs::Vector{Float64}     # DAYS SINCE THE FIRST EPOCH — see below
    nv2::Vector{Int}
    nt3::Vector{Int}
    # Whether a single file was cut on timestamp gaps. Recorded because the exported script
    # has to reproduce the same epoch structure, and the two produce different ones.
    split::Bool
end

# `tepochs` is relative, and that is not cosmetic. `rotate_star` turns it into a rotation phase
# as `2*pi*t/P`; handed a raw MJD of 55806 against a 54.8-day period that is a thousand
# revolutions, and the phase then depends on the epoch of the calendar rather than on the
# observation. Every demo subtracts the first epoch (`tepochs .- tepochs[1]`), so this does too,
# and the absolute MJD is kept separately for the table.

"""
    ModelEntry

A parametric surface: which `surface_type`, the parameter values, and which of them a fit is
allowed to move.

`free`, `bounds` and `ties` are three views of the same decision and are kept apart on
purpose. A parameter is FIXED (absent from `free`), FREE within `bounds`, or TIED — an
expression in the other parameters, evaluated before the model is built. A tick box cannot
express "expression", which is why this is a three-way choice and not a boolean.

Values are `Float64` even for `surface_type` and `ldtype`; `star_params` converts them back to
`Int` on the way out (see [`star_params`](@ref)). Keeping one numeric type here is what lets
the form, the bounds table and the tie expressions all be plain text round-trips.
"""
mutable struct ModelEntry
    name::String
    surface_type::Int
    params::Dict{Symbol,Float64}
    free::Set{Symbol}
    bounds::Dict{Symbol,Tuple{Float64,Float64}}
    ties::Dict{Symbol,String}
    secondary::Bool              # Roche: which component's potential convention to use
    # The COMPANION, for a binary: another `ModelEntry`, or `nothing` for a single star.
    #
    # A field on the primary rather than a second entry in `session.models`, because the two
    # are not independent models — they share one orbit, they are drawn in one frame, and
    # their χ² is a single number over the pair. `session.models` holds exactly one entry, and
    # this keeps that true for a binary as well.
    #
    # `Any` rather than `Union{Nothing,ModelEntry}`: the struct cannot name itself in its own
    # definition. Everything that reads it checks for `nothing` first.
    companion::Any
end

"""
    ImageEntry

A reconstructed surface map, with what produced it.

`x` is per-tessel, in the HEALPix NESTED order the tessellation uses, so it can be handed
straight to `plot2d_makie`/`plot_mollweide_makie`. `chi2` is the total, `ndata` the point
count, so a reduced value can be shown without re-running the criterion.
"""
mutable struct ImageEntry
    name::String
    x::Vector{Float64}
    nside_exp::Int               # log2(nside), i.e. what `tessellation_healpix` takes
    chi2::Float64
    ndata::Int
    regularizers::Vector{Any}
    source::String               # the dataset it came from
end

"""
    FitEntry

One completed fit, KEPT WHOLE.

Until this existed a sampler's result was reduced to a median and half an interval the moment
it arrived, and everything else — the draws, the evidence, the diagnostics — was discarded on
the worker thread. That threw away most of what a sampler is for: a run of Pigeons or NUTS is
minutes to hours of work whose ANSWER is the shape of the posterior, not two numbers per
parameter. It also made a wrong answer undetectable: a NUTS chain that had not converged
reported a tight interval, and only a separate χ² scan showed it was in the wrong place. One
look at the marginal would have said so.

`samples` is empty for an optimiser, which has no posterior — that is the honest
representation, and the panel says "no posterior" rather than drawing a delta function.

`logz` is `NaN` unless the sampler computes an evidence. Nested sampling and tempering do;
NUTS and the optimisers do not, and a number invented there would be worse than a blank,
because log(Z) is exactly what a model comparison would be read off.
"""
mutable struct FitEntry
    name::String
    method::Symbol
    model::String                  # which model entry it was fitted from
    surface_type::Int              # kept so a comparison across models means something
    names::Vector{Symbol}          # free parameters, in the order the table lists them
    best::Vector{Float64}          # point estimate: the posterior median for a sampler
    err::Vector{Float64}           # half the 16-84 interval, NaN for an optimiser
    q16::Vector{Float64}
    q84::Vector{Float64}
    samples::Matrix{Float64}       # nsamples x nfree; empty for an optimiser
    logz::Float64
    logzerr::Float64
    chi2::Float64
    ndata::Int
    diagnostics::String            # divergences, round trips, chains — sampler-specific
end

"""
    Session

Everything the window is looking at.

`current_*` are 1-based indices into the four vectors, or 0 for "none selected" — 0 rather
than `nothing` because the value crosses into QML, where every index is an integer and a
missing selection is conventionally out of range.
"""
mutable struct Session
    datasets::Vector{DatasetEntry}
    models::Vector{ModelEntry}
    images::Vector{ImageEntry}
    fits::Vector{FitEntry}
    current_dataset::Int
    current_model::Int
    current_image::Int
    current_fit::Int
    current_epoch::Int
    log::Vector{LogEntry}
end

Session() = Session(DatasetEntry[], ModelEntry[], ImageEntry[], FitEntry[],
                    0, 0, 0, 0, 1, LogEntry[])

# ── datasets ────────────────────────────────────────────────────────────────────────────

"""
    split_epochs(data; gap = 0.5) -> (per_epoch_data, mean_mjds)

Split one OIFITS into observing epochs on gaps in the V² timestamps.

`gap` is in days and defaults to half a day: within one night the V² points are minutes apart,
between nights at least several hours, so any threshold between those separates them. This is
the same rule the demo scripts use by hand, moved here so the GUI and the scripts cannot
disagree about what an epoch is.
"""
function split_epochs(data; gap::Real = 0.5)
    mjds = sort(data.v2_mjd)
    isempty(mjds) && return (Any[data], Float64[0.0])
    jumps = findall(diff(mjds) .> gap)
    starts = mjds[[1; jumps .+ 1]]
    ends   = mjds[[jumps; length(mjds)]]
    out  = Any[]
    tmid = Float64[]
    for i in eachindex(starts)
        idx = OITOOLS.set_data_filter(data; mjd_range = [starts[i] - 0.01, ends[i] + 0.01])
        d   = OITOOLS.filter_data(data, idx)
        push!(out, d)
        push!(tmid, Statistics.mean(d.v2_mjd))
    end
    return out, tmid
end

"""
    load_dataset!(session, path; ...) -> DatasetEntry
    load_dataset!(session, paths::AbstractVector; ...) -> DatasetEntry

Read one or more OIFITS files as ONE dataset, and make it current.

Several files is the normal case, not the exception: λ And is six nights in six files, and a
surface map is constrained by how the star looks across them. So a list of files becomes a list
of EPOCHS of one dataset, through `readoifits_multiepochs` — the same call every demo makes —
rather than several unrelated datasets that could not be reconstructed together.

A SINGLE file is split on gaps in its V² timestamps instead (`split_epochs`), because a file
like the Spica set holds three campaigns in one table. Both paths end at the same thing: a
vector of per-epoch `OIdata` and a vector of times relative to the first epoch.

`add = true` appends the files to the CURRENT dataset instead of starting a new one. The epoch
times are recomputed against the earliest epoch of the combined set, since that is the origin
the rotation phase is measured from.

`split = false` takes a single file as ONE epoch rather than cutting it on gaps in its V²
timestamps. It has no effect on a multi-file open, where each file is already one epoch.
"""
function load_dataset!(session::Session, path::AbstractString; kwargs...)
    return load_dataset!(session, [String(path)]; kwargs...)
end

function load_dataset!(session::Session, paths::AbstractVector;
                       gap::Real = 0.5, name::AbstractString = "", add::Bool = false,
                       split::Bool = true)
    ps = String.(paths)
    isempty(ps) && throw(ArgumentError("no files given"))
    if length(ps) == 1
        raw = OITOOLS.readoifits(ps[1]; warn = false, verbose = false)[1, 1]
        # `split = false` takes the file as ONE epoch. The gap rule is right for a file holding
        # several campaigns, and wrong for one night that happens to have a long gap in it —
        # a target lost to cloud and re-acquired two hours later is still one epoch of a star
        # whose rotation period is days.
        per, mjd = split ? split_epochs(raw; gap = gap) :
                   (Any[raw], Float64[Statistics.mean(raw.v2_mjd)])
    else
        # nwav x nepochs; ROTIR works one spectral bin at a time, so row 1 and every column.
        all = OITOOLS.readoifits_multiepochs(ps; warn = false, verbose = false)
        per = collect(Any, all[1, :])
        mjd = Float64[d.mean_mjd for d in per]
    end
    if add && (cur = current_dataset(session)) !== nothing
        append!(cur.data, per); append!(cur.mjd, mjd)
        append!(cur.nv2, [d.nv2 for d in per]); append!(cur.nt3, [d.nt3phi for d in per])
        append!(cur.paths, ps)
        # Sorted by time, so "epoch 3" means the third night rather than the third click.
        o = sortperm(cur.mjd)
        cur.data = cur.data[o]; cur.mjd = cur.mjd[o]
        cur.nv2 = cur.nv2[o];   cur.nt3 = cur.nt3[o]
        cur.tepochs = cur.mjd .- minimum(cur.mjd)
        cur.name = length(cur.paths) == 1 ? basename(cur.paths[1]) :
                   "$(basename(cur.paths[1])) +$(length(cur.paths) - 1)"
        session.current_epoch = 1
        return cur
    end
    nm = !isempty(name) ? String(name) :
         length(ps) == 1 ? basename(ps[1]) : "$(basename(ps[1])) +$(length(ps) - 1)"
    e = DatasetEntry(ps, nm, per, mjd, mjd .- minimum(mjd),
                     [d.nv2 for d in per], [d.nt3phi for d in per], split)
    push!(session.datasets, e)
    session.current_dataset = length(session.datasets)
    session.current_epoch = 1
    return e
end

current_dataset(s::Session) =
    s.current_dataset in eachindex(s.datasets) ? s.datasets[s.current_dataset] : nothing
current_model(s::Session) =
    s.current_model in eachindex(s.models) ? s.models[s.current_model] : nothing
current_image(s::Session) =
    s.current_image in eachindex(s.images) ? s.images[s.current_image] : nothing

current_fit(s::Session) =
    s.current_fit in eachindex(s.fits) ? s.fits[s.current_fit] : nothing

"""
    add_fit!(session, e) -> FitEntry

Keep a completed fit and select it.

Fits ACCUMULATE, unlike models — the GUI holds one model at a time on purpose, but the whole
point of keeping these is comparison: fit a sphere, switch to an ellipsoid, fit again, and the
evidence of the two is what says which the data prefer. Each entry records the surface type it
was fitted with, so the comparison is not between two rows that only look alike.
"""
function add_fit!(session::Session, e::FitEntry)
    push!(session.fits, e)
    session.current_fit = length(session.fits)
    return e
end

# ── models ──────────────────────────────────────────────────────────────────────────────

"""
    add_model!(session, surface_type; name = ..., secondary = false, kwargs...) -> ModelEntry

A new model at the schema defaults for that surface type, with nothing free yet.

Built through `default_star_params`, so the form and the validator agree with the geometry
code by construction, and a field this surface type does not have is rejected here rather than
several calls deep (see src/surface_schema.jl).
"""
function add_model!(session::Session, surface_type;
                    name::AbstractString = "", secondary::Bool = false, kwargs...)
    spec = surface_spec(surface_type)
    p = default_star_params(spec.code; kwargs...)
    params = Dict{Symbol,Float64}(k => Float64(v) for (k, v) in pairs(p) if k !== :surface_type)
    bounds = Dict{Symbol,Tuple{Float64,Float64}}(
        ps.name => (ps.lo, ps.hi) for ps in surface_params(spec.code))
    nm = isempty(name) ? "$(spec.name)_$(length(session.models) + 1)" : String(name)
    m = ModelEntry(nm, spec.code, params, Set{Symbol}(), bounds, Dict{Symbol,String}(),
                   secondary, nothing)
    push!(session.models, m)
    session.current_model = length(session.models)
    return m
end

"""
    star_params(m::ModelEntry) -> NamedTuple

The `star_params` the geometry code wants, with ties applied and the integer fields cast back.

`surface_type` and `ldtype` must be `Int`: `compute_radii` and `compute_ldmap` branch on them
with `==`, and a `Float64` 3.0 works by luck there while reading as a continuous parameter to
anything that iterates the schema.
"""
function star_params(m::ModelEntry)
    vals = apply_model_ties(m)
    pairs_ = Pair{Symbol,Any}[:surface_type => m.surface_type]
    for ps in surface_params(m.surface_type)
        haskey(vals, ps.name) || continue
        v = vals[ps.name]
        push!(pairs_, ps.name => ps.kind === :float ? v : round(Int, v))
    end
    return NamedTuple(pairs_)
end

"""
    apply_model_ties(m) -> Dict{Symbol,Float64}

Parameter values with every tie expression evaluated.

Ties are evaluated in dependency order by iterating to a fixed point rather than by sorting a
graph: the expressions are a handful and the values are cheap, so a few extra passes cost
nothing, and a cyclic tie stops after `length(ties)` rounds with the cycle reported instead of
looping forever.
"""
function apply_model_ties(m::ModelEntry)
    vals = copy(m.params)
    isempty(m.ties) && return vals
    for _ in 1:(length(m.ties) + 1)
        changed = false
        for (name, expr) in m.ties
            v = eval_tie(expr, vals)
            v === nothing && continue
            if !haskey(vals, name) || vals[name] != v
                vals[name] = v
                changed = true
            end
        end
        changed || return vals
    end
    @warn "tie expressions did not settle; a cycle is likely" ties = collect(keys(m.ties))
    return vals
end

"""
    eval_tie(expr, vals) -> Float64 or nothing

Evaluate one tie expression against the current parameter values.

The expression is ordinary Julia (`"P"`, `"2*rpole"`, `"180 - i"`), parsed and evaluated with
the parameter names bound. `nothing` on any failure — an unparseable expression, a name that
is not a parameter — because a half-typed expression is the normal state of a text field
somebody is editing, and it must not raise on every keystroke.
"""
function eval_tie(expr::AbstractString, vals::AbstractDict{Symbol,<:Real})
    s = strip(expr)
    isempty(s) && return nothing
    ex = try
        Meta.parse(s)
    catch
        return nothing
    end
    # Substitute PARAMETER names only, and leave every other symbol for `eval` to resolve.
    # Replacing all of them put the operator itself through the lookup — `Expr(:call, :-, …)`
    # carries `:-` in args[1] — so `"i - 60"` failed as "unknown -" and the tie silently never
    # applied. Leaving them alone is also what makes `sin`, `sqrt` and `pi` usable in a tie.
    sub(e) = e isa Symbol ? get(vals, e, e) :
             e isa Expr   ? Expr(e.head, map(sub, e.args)...) : e
    try
        v = Core.eval(Main, sub(ex))
        return v isa Real && isfinite(v) ? Float64(v) : nothing
    catch
        return nothing
    end
end

# ── images ──────────────────────────────────────────────────────────────────────────────

function add_image!(session::Session, name, x, nside_exp, chi2, ndata, regs, source)
    e = ImageEntry(String(name), Vector{Float64}(x), Int(nside_exp), Float64(chi2),
                   Int(ndata), collect(Any, regs), String(source))
    push!(session.images, e)
    session.current_image = length(session.images)
    return e
end
