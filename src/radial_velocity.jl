# =======================================================================================
# Radial velocities: tables, model, residuals
# =======================================================================================
# The CSV format is OITOOLS' (`OITOOLS/src/radial_velocity.jl`), so one file serves both
# packages: a header naming the columns `MJD,RV,sigma,star`, time in MJD, velocity and error
# in km/s (positive receding), star `a` (primary) or `b` (secondary). Columns are found by
# name, extra columns are ignored, lines starting with `#` are comments. ROTIR also reads a
# `JD` column in place of `MJD`, and an optional `instrument` column.
#
# Times are stored as JD, the unit of every orbit function in src/orbits.jl.
#
# Conventions match OITOOLS' `orbit_rv`: ω is the argument of periapsis of the relative
# orbit (the secondary's), V_a = γ − K1[cos(ν+ω) + e cos ω], V_b = γ + K2[cos(ν+ω) + e cos ω].
# OITOOLS' MJD0 is T0 − 2400000.5.

const MJD_TO_JD = 2400000.5

"""
    RVData

Radial velocities of a two-star system: `jd` (JD), `rv` and `sigma` (km/s), `star`
(`:a` primary, `:b` secondary), `instrument` (`""` when the table has none) and the `path`
read from. Build one with [`read_rv`](@ref).
"""
struct RVData
    jd::Vector{Float64}
    rv::Vector{Float64}
    sigma::Vector{Float64}
    star::Vector{Symbol}
    instrument::Vector{String}
    path::String
    function RVData(jd, rv, sigma, star, instrument, path)
        n = length(jd)
        length(rv) == length(sigma) == length(star) == length(instrument) == n ||
            throw(ArgumentError("RVData: columns differ in length"))
        all(s -> s === :a || s === :b, star) ||
            throw(ArgumentError("RVData: star is :a (primary) or :b (secondary)"))
        all(>(0), sigma) || throw(ArgumentError("RVData: every sigma must be positive"))
        new(Float64.(jd), Float64.(rv), Float64.(sigma), Symbol.(star), String.(instrument),
            String(path))
    end
end

Base.length(rv::RVData) = length(rv.jd)

function Base.show(io::IO, rv::RVData)
    na = count(==(:a), rv.star)
    print(io, "RVData: ", length(rv), " velocities (", na, " a, ", length(rv) - na, " b)")
    ins = filter(!isempty, rv_instruments(rv))
    isempty(ins) || print(io, ", instruments ", join(ins, ", "))
    isempty(rv.path) || print(io, ", ", basename(rv.path))
end

"""
    rv_instruments(rv) -> Vector{String}

The distinct instruments in file order. The first is the reference for the per-instrument
γ offsets. A table without an `instrument` column gives `[""]`.
"""
rv_instruments(rv::RVData) = unique(rv.instrument)

"""
    read_rv(path) -> RVData
    read_rv(path1, path2) -> RVData

Read radial velocities.

With one path: a CSV table with a header, the format OITOOLS reads.

| column | content |
|---|---|
| `MJD` or `JD` | time |
| `RV` | velocity, km/s, positive receding |
| `sigma` | its error, km/s, > 0 |
| `star` | `a` (primary) or `b` (secondary) |
| `instrument` | optional; one γ offset per instrument in a fit |

Columns are found by name (case-insensitive), extra columns are ignored, blank lines and
lines starting with `#` are skipped.

With two paths: whitespace-separated files without header, `JD  RV  [sigma]`, one per star
(primary first). A missing sigma column reads as 1 km/s.
"""
function read_rv(path::AbstractString)
    name = basename(path)
    lines = [strip(l) for l in eachline(path) if !isempty(strip(l)) && !startswith(strip(l), "#")]
    isempty(lines) && throw(ArgumentError("$name: no radial velocities"))
    occursin(',', lines[1]) ||
        throw(ArgumentError("$name: expected a CSV header such as MJD,RV,sigma,star; for " *
                            "two whitespace files (JD RV sigma, one per star) call " *
                            "read_rv(path1, path2)"))
    head = lowercase.(strip.(split(lines[1], ',')))
    function col(n; required = true)
        i = findfirst(==(n), head)
        i === nothing && required &&
            throw(ArgumentError("$name: no column \"$n\" in the header (MJD,RV,sigma,star)"))
        return i
    end
    it_mjd, it_jd = col("mjd"; required = false), col("jd"; required = false)
    (it_mjd === nothing && it_jd === nothing) &&
        throw(ArgumentError("$name: no time column: the header needs MJD or JD"))
    (it_mjd !== nothing && it_jd !== nothing) &&
        throw(ArgumentError("$name: both MJD and JD columns; keep one"))
    itime, offset = it_mjd !== nothing ? (it_mjd, MJD_TO_JD) : (it_jd, 0.0)
    ir, is, ist = col("rv"), col("sigma"), col("star")
    iin = col("instrument"; required = false)
    jd, v, σ, star, inst = Float64[], Float64[], Float64[], Symbol[], String[]
    for (k, l) in enumerate(lines[2:end])
        f = strip.(split(l, ','))
        length(f) >= length(head) ||
            throw(ArgumentError("$name, row $k: $(length(f)) fields, the header has $(length(head))"))
        s = Symbol(lowercase(f[ist]))
        s in (:a, :b) || throw(ArgumentError("$name, row $k: star is a or b, not $(f[ist])"))
        push!(jd, parse(Float64, f[itime]) + offset)
        push!(v, parse(Float64, f[ir])); push!(σ, parse(Float64, f[is])); push!(star, s)
        push!(inst, iin === nothing ? "" : String(f[iin]))
    end
    isempty(jd) && throw(ArgumentError("$name: a header but no velocities"))
    all(>(0), σ) || throw(ArgumentError("$name: every sigma must be positive"))
    return RVData(jd, v, σ, star, inst, String(path))
end

function read_rv(path1::AbstractString, path2::AbstractString)
    function rd(p)
        rows = [split(strip(l)) for l in eachline(p) if !isempty(strip(l)) && !startswith(strip(l), "#")]
        t = [parse(Float64, r[1]) for r in rows]
        v = [parse(Float64, r[2]) for r in rows]
        s = [length(r) >= 3 ? parse(Float64, r[3]) : 1.0 for r in rows]
        return t, v, s
    end
    t1, v1, s1 = rd(path1)
    t2, v2, s2 = rd(path2)
    n1, n2 = length(t1), length(t2)
    return RVData([t1; t2], [v1; v2], [s1; s2], [fill(:a, n1); fill(:b, n2)],
                  fill("", n1 + n2), string(path1, " + ", basename(path2)))
end

"""
    write_rv(path, rv; comment = "")

Write `rv` as a CSV table that [`read_rv`](@ref) and OITOOLS' `read_rv` both read:
`MJD,RV,sigma,star`, plus `instrument` when the table has instruments. `comment` lines are
written first, each prefixed with `#`.
"""
function write_rv(path::AbstractString, rv::RVData; comment::AbstractString = "")
    withinst = any(!isempty, rv.instrument)
    open(path, "w") do io
        for l in split(comment, '\n')
            isempty(strip(l)) || println(io, "# ", l)
        end
        println(io, withinst ? "MJD,RV,sigma,star,instrument" : "MJD,RV,sigma,star")
        for k in eachindex(rv.jd)
            # Rounded to 1e-8 d (~1 ms): JD → MJD subtraction leaves float noise in the last
            # digits that would otherwise be written out verbatim.
            print(io, round(rv.jd[k] - MJD_TO_JD; digits = 8), ",", rv.rv[k], ",",
                  rv.sigma[k], ",", rv.star[k])
            withinst && print(io, ",", rv.instrument[k])
            println(io)
        end
    end
    return path
end

"""
    rv_model(bp, rv; K1, K2, gamma, offsets = Dict()) -> Vector

Model velocity (km/s) at every row of `rv`, for that row's star: `binary_RV` at the
row's time, plus the row's instrument offset `offsets[instrument]` (0 when absent).
`bp` holds the orbital elements as `binary_RV` reads them (`ω`, `e`, `P`, `T0`, `dP`,
optionally `dω`).
"""
function rv_model(bp, rv::RVData; K1, K2, gamma, offsets = Dict{String,Float64}())
    v1, v2 = binary_RV(bp, rv.jd; K1 = K1, K2 = K2, γ = gamma)
    T = promote_type(eltype(v1), eltype(v2))
    m = Vector{T}(undef, length(rv))
    @inbounds for k in eachindex(m)
        m[k] = (rv.star[k] === :a ? v1[k] : v2[k]) + get(offsets, rv.instrument[k], 0.0)
    end
    return m
end

# σ inflated by the star's jitter, √(σ² + s²).
_rv_sigma(rv::RVData, k, jitter) =
    sqrt(rv.sigma[k]^2 + (rv.star[k] === :a ? jitter.a : jitter.b)^2)

"""
    rv_residuals(bp, rv; K1, K2, gamma, offsets = Dict(), jitter = (a = 0, b = 0)) -> Vector

Normalised residuals `(rv − model)/√(σ² + s²)` per row, with `s` the jitter of the row's
star (km/s, added in quadrature to the tabulated error).
"""
function rv_residuals(bp, rv::RVData; K1, K2, gamma, offsets = Dict{String,Float64}(),
                      jitter = (a = 0.0, b = 0.0))
    m = rv_model(bp, rv; K1 = K1, K2 = K2, gamma = gamma, offsets = offsets)
    return [(rv.rv[k] - m[k]) / _rv_sigma(rv, k, jitter) for k in eachindex(m)]
end

"`sum(abs2, rv_residuals(...))`: χ² of the velocities, with jitter-inflated errors."
rv_chi2(bp, rv::RVData; kw...) = sum(abs2, rv_residuals(bp, rv; kw...))

"""
    rv_nll2(bp, rv; K1, K2, gamma, offsets = Dict(), jitter = (a = 0, b = 0))

The velocities' contribution to −2 ln L, offset so that it equals [`rv_chi2`](@ref) when
the jitter is zero:

    Σ (rv − model)²/(σ² + s²) + Σ ln((σ² + s²)/σ²)

The log term is what a fitted jitter needs: without it, χ² falls without limit as `s` grows.
`fit_orbit` minimises this, not the bare χ².
"""
function rv_nll2(bp, rv::RVData; K1, K2, gamma, offsets = Dict{String,Float64}(),
                 jitter = (a = 0.0, b = 0.0))
    m = rv_model(bp, rv; K1 = K1, K2 = K2, gamma = gamma, offsets = offsets)
    s = zero(eltype(m))
    @inbounds for k in eachindex(m)
        σ2 = _rv_sigma(rv, k, jitter)^2
        s += (rv.rv[k] - m[k])^2 / σ2 + log(σ2 / rv.sigma[k]^2)
    end
    return s
end

"""
    rv_start(rv) -> (K1, K2, gamma)

Starting values from the data alone: γ the mean velocity, each K half its star's
peak-to-peak range (50 km/s for a star with fewer than two points).
"""
function rv_start(rv::RVData)
    a = rv.rv[rv.star .=== :a]; b = rv.rv[rv.star .=== :b]
    half(x) = length(x) >= 2 ? (maximum(x) - minimum(x)) / 2 : 50.0
    return (K1 = half(a), K2 = half(b), gamma = sum(rv.rv) / length(rv))
end

# Parameter name of an instrument's γ offset: `dgamma_<name>`, the name reduced to letters,
# digits and `_` so it can appear in a tie expression.
function rv_offset_name(inst::AbstractString)
    s = replace(String(inst), r"[^A-Za-z0-9_]" => "_")
    (isempty(s) || isdigit(first(s))) && (s = "i" * s)
    return Symbol("dgamma_", s)
end

# ── the amplitudes from the orbit: parallax and mass ratio ──────────────────────────────────

"km/s per AU/day."
const AU_PER_DAY_KMS = 1.495978707e8 / 86400.0

"""
    rv_amplitudes(a, i, P, e, plx, q) -> (K1, K2)

Semi-amplitudes (km/s) of the two stars of an orbit of angular semi-major axis `a` (mas),
inclination `i` (deg), period `P` (days) and eccentricity `e`, at parallax `plx` (mas), for a
mass ratio `q = M2/M1`:

    K = 2π (a/plx) sin i / (P √(1 − e²))      (a/plx in AU, P in days)
    K1 = K q/(1 + q),    K2 = K/(1 + q)

With interferometry fitting `a`, `i`, `P` and `e`, the velocities need only `plx` and `q`.
"""
function rv_amplitudes(a, i, P, e, plx, q)
    K = 2π * (a / plx) * abs(sind(i)) / (P * sqrt(1 - e^2)) * AU_PER_DAY_KMS
    return K * q / (1 + q), K / (1 + q)
end

"""
    rv_physical(a, i, P, e, K1, K2) -> (plx, q)

The inverse of [`rv_amplitudes`](@ref): parallax (mas) and mass ratio from the amplitudes.
"""
function rv_physical(a, i, P, e, K1, K2)
    K = K1 + K2
    plx = 2π * a * abs(sind(i)) / (P * sqrt(1 - e^2)) * AU_PER_DAY_KMS / K
    return plx, K1 / K2
end

"""
    dynamical_masses(a, P, plx, q) -> (; M, M1, M2, d)

Total and component masses (M☉) and distance (pc) of an orbit of angular semi-major axis `a`
(mas) and period `P` (days) at parallax `plx` (mas), mass ratio `q = M2/M1`, from Kepler's third
law, M = (a/plx)³ / (P/365.25)².
"""
function dynamical_masses(a, P, plx, q)
    M = (a / plx)^3 / (P / 365.25)^2
    return (; M, M1 = M / (1 + q), M2 = M * q / (1 + q), d = 1000 / plx)
end

"""
    rv_curves(bp, rv; K1, K2, gamma, offsets = Dict(), jitter = (a = 0, b = 0),
              xaxis = :phase, tref = nothing, n = 400) -> NamedTuple

Everything a radial-velocity plot draws, computed once so every backend draws the same thing:

- `x`, `y`, `sigma`, `resid`, `star`, `instrument`: one entry per row of `rv`. `y` is the
  velocity shifted by its instrument's offset, so every instrument sits on the same curve;
  `sigma` is √(σ² + s²) with `s` the star's jitter; `resid` is `rv − model` (km/s).
- `xc`, `va`, `vb`: the model curves of the primary and the secondary.
- `xlabel`.

`xaxis = :phase` folds on (t − T0)/P mod 1 and draws one cycle at `tref` (JD; default the mean
time of the table): with apsidal motion or a period derivative a folded curve is exact only
at that epoch. `xaxis = :time` uses MJD and draws the model through every observation.
"""
function rv_curves(bp, rv::RVData; K1, K2, gamma, offsets = Dict{String,Float64}(),
                   jitter = (a = 0.0, b = 0.0), xaxis::Symbol = :phase, tref = nothing,
                   n::Int = 400)
    xaxis in (:phase, :time) || throw(ArgumentError("xaxis is :phase or :time"))
    P, T0 = Float64(bp.P), Float64(bp.T0)
    m = rv_model(bp, rv; K1 = K1, K2 = K2, gamma = gamma, offsets = offsets)
    k = eachindex(rv.jd)
    sig = [_rv_sigma(rv, i, jitter) for i in k]
    off = [get(offsets, rv.instrument[i], 0.0) for i in k]
    x = xaxis === :phase ? mod.((rv.jd .- T0) ./ P, 1.0) : rv.jd .- MJD_TO_JD
    if xaxis === :phase
        t0 = tref === nothing ? sum(rv.jd) / length(rv) : Float64(tref)
        φ = collect(range(0.0, 1.0; length = n))
        tc = T0 .+ (floor((t0 - T0) / P) .+ φ) .* P
        xc = φ
    else
        lo, hi = extrema(rv.jd); pad = 0.03 * max(hi - lo, P)
        npt = clamp(ceil(Int, n / 10 * (hi - lo + 2pad) / P), n, 50n)
        tc = collect(range(lo - pad, hi + pad; length = npt))
        xc = tc .- MJD_TO_JD
    end
    va, vb = binary_RV(bp, tc; K1 = K1, K2 = K2, γ = gamma)
    return (; x, y = rv.rv .- off, sigma = sig, resid = rv.rv .- m, star = rv.star,
            instrument = rv.instrument, xc, va, vb,
            xlabel = xaxis === :phase ? "Orbital phase, (t − T0)/P" : "MJD")
end
