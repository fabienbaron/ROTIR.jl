# kurucz_intensity.jl
# ---------------------------------------------------------------------------
# Kurucz ATLAS9 specific-intensity packs -> a ROTIR [`RectGrid4`](@ref).
#
# WHY THIS EXISTS. Korg.jl's `interpolate_marcs` stops at Teff = 8000 K, which covers beta Cas
# (7208 K) and nothing hotter. This repo's other targets are Spica (25300 / 20585 K), Vega
# (~9500 K), Regulus (~15000 K) and beta Lyr — all outside it. Kurucz's intensity packs run
# **3500 to 50000 K** at 17 emergent angles, so one reader covers every hot rotator here.
#
# WHICH FILE. `http://kurucz.harvard.edu/grids/gridp00/` holds two solar-metallicity packs and
# the choice matters more than the names suggest:
#
#   ip00k2.pck19   67.5 MB   Teff 3500-50000 K   <- the one you want for hot stars
#   ip00k2new.pck  31.5 MB   Teff 3500- 8750 K   <- newer opacities, cool stars only
#
# `ip00k2new.pck` is the newer ODFNEW grid but tops out at 8750 K, i.e. no better than Korg.
# Other metallicities live in the sibling `gridm05`, `gridm10`, … directories under the same
# naming (`im05k2.pck19` and so on). `grids.html` describes the `I` family as "intensities for
# each model; limbdarkening for 17 angles".
#
# THE FORMAT, from Kurucz's own `intensitypack.for` in the same directory, confirmed against
# the bytes. Per model:
#
#   line 1   "EFF   3500.  GRAVITY 0.00000  SDSC GRID  [+0.0]   VTURB 2.0 KM/S    L/H 1.25"
#   line 2   "  wl(nm)    Inu(ergs/cm**2/s/hz/ster) for 17 mu in 1221 frequency intervals"
#   line 3   "            1.000   .900  .800 … .010"
#   then 1221 lines, Fortran `format(F9.2,1pe10.3,16I6)`:
#              wl(nm)   I(mu=1)   16 integers
#
# Two traps in that, both of which produce a plausible-looking wrong grid:
#
#  1. **The two packs differ by ONE COLUMN, and nothing announces it.** Fortran's old
#     carriage-control convention makes column 1 a printer control character. `ip00k2.pck19`
#     KEEPS it (its header reads `TEFF   3500.`), `ip00k2new.pck` has had it STRIPPED (its
#     header reads `EFF   3500.` and every data line is shifted one column left). Since the
#     record is fixed-width `format(F9.2,1pe10.3,16I6)` and its fields abut — a 160000.00 nm
#     wavelength fills F9.2 exactly, and a ratio of 100000 fills I6 — the columns cannot be
#     recovered by splitting on whitespace. The parser therefore DETECTS the shift from the
#     header and offsets every field by it. Hard-coding either layout silently misreads the
#     other: one column left turns the continuum field into the wavelength's last digit plus
#     nine characters of the mantissa, which still parses as a float.
#  2. **The 16 integers are RATIOS, not intensities.** `intensitypack.for` stores
#     `INU(mu) = I(mu)/I(mu=1) * 100000 + 0.5`, so the absolute intensity is
#     `I(mu) = centi * INU(mu) / 100000` with `centi` the `1pe10.3` column. Reading them as
#     intensities gives numbers ~1e5 times too large, uniformly — which cancels out of a
#     flux-normalised visibility and so survives every test that does not check absolute flux.
#
# And one convention: the 17 mu values DESCEND (1.0 … 0.01), while `RectGrid4` requires
# ascending axes. [`read_kurucz_intensity`](@ref) reverses them.
# ---------------------------------------------------------------------------

"The 17 emergent angles of a Kurucz intensity pack, as written (descending)."
const KURUCZ_MU = (1.000, 0.900, 0.800, 0.700, 0.600, 0.500, 0.400, 0.300, 0.250,
                   0.200, 0.150, 0.125, 0.100, 0.075, 0.050, 0.025, 0.010)

"""
    KuruczModel

One model atmosphere out of an intensity pack: its `Teff`, `logg`, the wavelength grid in
**metres** (converted from the file's nm), and `I` as `(17, nλ)` in erg/cm²/s/Hz/sr with the
μ axis still in the file's descending order.
"""
struct KuruczModel{T<:AbstractFloat}
    Teff::T
    logg::T
    λ::Vector{T}
    I::Matrix{T}
end

"""
    read_kurucz_models(path; T = Float64, λrange = nothing, verbose = true) -> Vector{KuruczModel}

Parse every model in a Kurucz intensity pack.

`λrange` is `(λmin, λmax)` in **metres** and is strongly recommended: a full pack is 1221
wavelengths × 17 angles × ~500 models, and an H-band grid needs perhaps 30 of those
wavelengths. Restricting here keeps the result small and the parse fast.

Reads the `1pe10.3` continuum column and the 16 packed ratio integers and multiplies them out,
so `I` is absolute — see the note on trap 2 in this file's header.
"""
function read_kurucz_models(path::AbstractString; T::Type = Float64,
                            λrange = nothing, verbose::Bool = true)
    models = KuruczModel{T}[]
    # Matches both layouts: `TEFF` contains `EFF`, so the regex is layout-agnostic and the
    # shift is deduced separately from whether the `T` survived.
    hdr = r"EFF\s+([0-9.]+)\s+GRAVITY\s+([0-9.]+)"
    nμ = length(KURUCZ_MU)
    shift = -1          # 0 = carriage-control column present (.pck19), 1 = stripped (new.pck)
    open(path, "r") do io
        Teff = logg = zero(T)
        λbuf = T[]; Ibuf = Vector{T}[]
        inmodel = false
        flush!() = begin
            if inmodel && !isempty(λbuf)
                # (nλ, nμ) rows -> (nμ, nλ), and reverse μ so the axis ascends.
                M = Matrix{T}(undef, nμ, length(λbuf))
                @inbounds for (j, row) in enumerate(Ibuf), i in 1:nμ
                    M[i, j] = row[i]
                end
                push!(models, KuruczModel{T}(Teff, logg, copy(λbuf), M))
            end
            empty!(λbuf); empty!(Ibuf)
        end
        for line in eachline(io)
            m = match(hdr, line)
            if m !== nothing
                flush!()
                Teff = parse(T, m.captures[1]); logg = parse(T, m.captures[2])
                # Deduce the column shift ONCE, from the first header seen.
                shift < 0 && (shift = occursin("TEFF", line) ? 0 : 1)
                inmodel = true
                continue
            end
            inmodel || continue
            w1 = 9 - shift                       # width of the F9.2 wavelength field
            length(line) < w1 + 10 && continue
            # `format(F9.2,1pe10.3,16I6)`: fixed columns, so slice rather than split — the
            # fields abut at the extremes of their ranges (see trap 1 in the header).
            wl = tryparse(T, line[1:w1]);              wl === nothing && continue
            ct = tryparse(T, line[w1+1:w1+10]);        ct === nothing && continue
            λm = wl * T(1e-9)                              # nm -> m
            λrange !== nothing && !(λrange[1] <= λm <= λrange[2]) && continue
            row = Vector{T}(undef, nμ)
            row[1] = ct                                    # mu = 1 is absolute
            ok = true
            @inbounds for k in 2:nμ
                a = w1 + 11 + 6*(k-2); b = a + 5
                b > length(line) && (ok = false; break)
                r = tryparse(Int, strip(line[a:b]))
                r === nothing && (ok = false; break)
                row[k] = ct * T(r) / T(100_000)            # RATIO, not intensity
            end
            ok || continue
            push!(λbuf, λm); push!(Ibuf, row)
        end
        flush!()
    end
    verbose && @info "read_kurucz_models" file=basename(path) models=length(models) λ_per_model=(isempty(models) ? 0 : length(first(models).λ))
    isempty(models) && error("read_kurucz_models: no models parsed from $path. If `λrange` " *
                             "is set, check it overlaps the pack's coverage (9.09 nm to " *
                             "160 um); otherwise the file may not be a Kurucz intensity " *
                             "pack — the header must contain `EFF <Teff> GRAVITY <logg>`.")
    return models
end

"""
    read_kurucz_intensity(path; λrange, Teff = nothing, logg = nothing, T = Float64)
        -> RectGrid4

A Kurucz intensity pack as a ROTIR grid, ready for [`TabulatedProvider`](@ref).

`λrange` is `(λmin, λmax)` in metres and is required — see [`read_kurucz_models`](@ref).
`Teff`/`logg` optionally restrict the axes to the models you want; the default takes every
(Teff, logg) the pack contains that shares the full λ sampling.

**The pack is not a complete rectangle.** Kurucz computes each Teff only over the gravities
where a model converges, so the hot end has no low-gravity entries and the cool end no high
ones. A `RectGrid4` must be rectangular, so this fills missing (Teff, logg) nodes by copying
the nearest available `logg` at the same `Teff`. That is a real approximation and it is
deliberate rather than silent: restricting `logg` to the span your surface actually needs
avoids relying on it. A warning names how many nodes were filled.

μ comes out ASCENDING (0.01 … 1.0), reversed from the file, because `RectGrid4` requires it.
"""
function read_kurucz_intensity(path::AbstractString; λrange,
                               Teff = nothing, logg = nothing,
                               T::Type = Float64, verbose::Bool = true)
    ms = read_kurucz_models(path; T = T, λrange = λrange, verbose = verbose)
    nλ = length(first(ms).λ)
    full = filter(m -> length(m.λ) == nλ, ms)
    length(full) < length(ms) && @warn "read_kurucz_intensity: dropped $(length(ms)-length(full)) model(s) with a different λ sampling"
    Ta = sort!(unique(m.Teff for m in full)); ga = sort!(unique(m.logg for m in full))
    Teff !== nothing && (Ta = filter(t -> Teff[1] <= t <= Teff[2], Ta))
    logg !== nothing && (ga = filter(g -> logg[1] <= g <= logg[2], ga))
    (length(Ta) >= 2 && length(ga) >= 2) ||
        error("read_kurucz_intensity: need >=2 Teff and >=2 logg nodes after filtering " *
              "(got $(length(Ta)) x $(length(ga))); widen `Teff`/`logg`")
    bykey = Dict((m.Teff, m.logg) => m for m in full)
    nμ = length(KURUCZ_MU)
    vals = Array{T,4}(undef, length(Ta), length(ga), nμ, nλ)
    nfill = 0
    for (i, t) in enumerate(Ta), (j, g) in enumerate(ga)
        m = get(bykey, (t, g), nothing)
        if m === nothing
            # Nearest available gravity AT THE SAME Teff — never a different Teff, which would
            # mix in a different temperature's limb darkening.
            cand = [gg for (tt, gg) in keys(bykey) if tt == t]
            isempty(cand) && error("read_kurucz_intensity: Teff = $t has no models at all")
            m = bykey[(t, cand[argmin(abs.(cand .- g))])]
            nfill += 1
        end
        @inbounds for k in 1:nμ, l in 1:nλ
            vals[i, j, k, l] = m.I[nμ + 1 - k, l]    # reverse μ to ascending
        end
    end
    nfill > 0 && @warn "read_kurucz_intensity: filled $nfill of $(length(Ta)*length(ga)) " *
                       "(Teff, logg) nodes from the nearest gravity at the same Teff — the " *
                       "pack is not a complete rectangle. Narrow `logg` to the span your " *
                       "surface needs if this matters."
    return RectGrid4(T.(Ta), T.(ga), T[KURUCZ_MU[nμ + 1 - k] for k in 1:nμ],
                     first(full).λ, vals)
end
