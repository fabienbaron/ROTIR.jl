# Reading and writing a surface map as FITS.
#
# A ROTIR map is a bare `Vector` of per-tessel values — nothing in it says which tessellation
# it belongs to, whether it is a temperature or an intensity, or what geometry it was fitted
# against. Saved on its own it is unusable a week later: HEALPix order is not recoverable from
# the length alone (nside 8 and a 768-element long-lat grid are the same number), and without
# the star parameters the map cannot be turned back into a χ².
#
# So the file carries the map AND everything needed to rebuild the model that produced it:
#
#   Primary HDU   the map itself, as an image, plus the tessellation in the header
#   PARAMS        one row per `star_params` field: name, value, and the Julia type it had
#   EPOCHS        the observation times, when the caller has them
#
# The PARAMS type column is what makes the round trip exact. `surface_type` and `ldtype` are
# `Int` and the code branches on them with `==`; a Float64 3.0 compares equal but propagates
# as a float through `create_star_multiepochs` and changes what `compute_radii` builds. FITS
# has no integer column type that survives a mixed table, so the type travels beside the value
# and `load_surface_map` restores it.

"""
    save_surface_map(path, x, params; nside_exp, kwargs...) -> path

Write a surface map to FITS, together with everything needed to reproduce its χ².

`x` is the per-tessel map, `params` the `star_params` NamedTuple it was computed with, and
`nside_exp` the HEALPix level as `tessellation_healpix` takes it — log₂(nside), so 3 means
nside 8 and 768 tessels.

Optional:

  * `tepochs`, `mjd` — the observation times. `tepochs` is what `create_star_multiepochs`
    takes (days since the first epoch), `mjd` the absolute dates; both are written when given.
  * `field` — what the values ARE (`:temperature` by default, `:intensity` for a map with limb
    darkening already multiplied in). A map of one kind fed to a χ² expecting the other is
    silently wrong, which is why it is recorded.
  * `secondary`, `chi2`, `ndata`, `comment` — provenance, written to the primary header.

To get the χ² back:

```julia
m = load_surface_map("map.fits")
tess  = tessellation_healpix(m.nside_exp)
stars = create_star_multiepochs(tess, m.params, m.tepochs; secondary = m.secondary)
setup_oi!(data, stars)
chi2_breakdown(m.x, data, stars)
```
"""
function save_surface_map(path::AbstractString, x::AbstractVector, params::NamedTuple;
                          nside_exp::Integer,
                          tessellation::Symbol = :healpix,
                          tepochs::Union{Nothing,AbstractVector} = nothing,
                          mjd::Union{Nothing,AbstractVector} = nothing,
                          secondary::Bool = false,
                          field::Symbol = :temperature,
                          chi2::Real = NaN, ndata::Integer = 0,
                          comment::AbstractString = "")
    tessellation === :healpix || tessellation === :longlat ||
        throw(ArgumentError("tessellation must be :healpix or :longlat (got $tessellation)"))
    npix = length(x)
    if tessellation === :healpix
        exp_npix = 12 * (2^Int(nside_exp))^2
        npix == exp_npix ||
            throw(ArgumentError("a level-$(nside_exp) HEALPix map has $(exp_npix) tessels, " *
                                "not $(npix) — nside_exp is log2(nside)"))
    end

    hdr = FITSIO.FITSHeader(
        String[], Any[], String[])
    set!(k, v, c) = (hdr[k] = v; FITSIO.set_comment!(hdr, k, c))
    set!("ORIGIN",   "ROTIR",                   "written by ROTIR.save_surface_map")
    set!("DATE",     Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS"), "file creation date")
    set!("TESSEL",   String(tessellation),      "tessellation the map is on")
    set!("NSIDEEXP", Int(nside_exp),            "log2(nside), what tessellation_healpix takes")
    set!("NSIDE",    2^Int(nside_exp),          "HEALPix nside")
    set!("NPIX",     npix,                      "number of tessels")
    set!("ORDERING", "NESTED",                  "HEALPix ordering ROTIR uses")
    set!("FIELD",    String(field),             "temperature or intensity")
    set!("SECONDRY", secondary,                 "map is the secondary of a binary")
    isfinite(chi2) && set!("CHI2", Float64(chi2), "chi2 of this map against its data")
    ndata > 0      && set!("NDATA", Int(ndata),   "number of data points in CHI2")
    isempty(comment) || set!("NOTE", String(comment), "")

    names  = String[]; values = Float64[]; types = String[]
    for k in keys(params)
        v = getfield(params, k)
        v isa Real || continue          # a NamedTuple field that is not a number is not a
                                        # fit parameter; nothing downstream reads one.
        push!(names, String(k)); push!(values, Float64(v))
        push!(types, v isa Integer ? "I" : v isa Bool ? "B" : "D")
    end

    FITSIO.FITS(String(path), "w") do f
        write(f, collect(Float64, x); header = hdr)
        write(f, Dict("NAME" => names, "VALUE" => values, "TYPE" => types); name = "PARAMS")
        if tepochs !== nothing || mjd !== nothing
            n = tepochs === nothing ? length(mjd) : length(tepochs)
            cols = Dict{String,Any}()
            cols["TDAY"] = tepochs === nothing ? fill(NaN, n) : collect(Float64, tepochs)
            cols["MJD"]  = mjd     === nothing ? fill(NaN, n) : collect(Float64, mjd)
            write(f, cols; name = "EPOCHS")
        end
    end
    return String(path)
end

"""
    load_surface_map(path) -> NamedTuple

Read back what [`save_surface_map`](@ref) wrote.

Returns `x`, `params`, `nside_exp`, `tessellation`, `field`, `secondary`, `tepochs`, `mjd`,
`chi2`, `ndata` and `header`. `params` comes back as a NamedTuple with the field types it was
saved with, so it can go straight into `create_star_multiepochs`.

A file without a PARAMS table still loads — `params` is then empty and only the map and its
tessellation are available, which is enough to draw it but not to refit it.
"""
function load_surface_map(path::AbstractString)
    FITSIO.FITS(String(path), "r") do f
        hdr = read_header(f[1])
        x   = Float64.(vec(read(f[1])))
        g(k, d) = haskey(hdr, k) ? hdr[k] : d
        nside_exp = Int(g("NSIDEEXP", round(Int, log2(sqrt(length(x) / 12)))))
        tessel    = Symbol(strip(String(g("TESSEL", "healpix"))))
        field     = Symbol(strip(String(g("FIELD", "temperature"))))

        nms = Symbol[]; vls = Any[]
        if _has_hdu(f, "PARAMS")
            t  = f["PARAMS"]
            pn = read(t, "NAME"); pv = read(t, "VALUE")
            pt = _has_col(t, "TYPE") ? read(t, "TYPE") : fill("D", length(pn))
            for j in eachindex(pn)
                push!(nms, Symbol(strip(String(pn[j]))))
                s = strip(String(pt[j]))
                push!(vls, s == "I" ? round(Int, pv[j]) :
                           s == "B" ? (pv[j] != 0)      : Float64(pv[j]))
            end
        end
        tep = nothing; mjd = nothing
        if _has_hdu(f, "EPOCHS")
            t = f["EPOCHS"]
            _has_col(t, "TDAY") && (tep = Float64.(read(t, "TDAY")))
            _has_col(t, "MJD")  && (mjd = Float64.(read(t, "MJD")))
            tep !== nothing && all(isnan, tep) && (tep = nothing)
            mjd !== nothing && all(isnan, mjd) && (mjd = nothing)
        end
        return (x = x, params = NamedTuple{Tuple(nms)}(Tuple(vls)),
                nside_exp = nside_exp, tessellation = tessel, field = field,
                secondary = Bool(g("SECONDRY", false)),
                tepochs = tep, mjd = mjd,
                chi2 = Float64(g("CHI2", NaN)), ndata = Int(g("NDATA", 0)),
                header = hdr)
    end
end

# `FITS` indexing by name throws rather than returning nothing, and an older file simply has
# no such table — a missing PARAMS is a thinner file, not a corrupt one.
function _has_hdu(f, name::AbstractString)
    try
        f[name]; return true
    catch
        return false
    end
end
_has_col(t, name::AbstractString) = name in FITSIO.colnames(t)
