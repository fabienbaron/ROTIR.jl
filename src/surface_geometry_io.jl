# Reading and writing a model's GEOMETRY as FITS.
#
# The sibling of `surface_map_io.jl`, and deliberately the other half of it: that file saves
# the MAP — what the surface is emitting — and this one saves the SURFACE, with no map values
# in it at all. Two files rather than one because they have different lifetimes: a map belongs
# to a reconstruction, a geometry belongs to a model, and the same geometry is reused across
# every map fitted on it.
#
# WHY THE VERTICES AND NOT JUST THE PARAMETERS. `create_star(tess, params, t)` is
# deterministic, so in principle the parameters plus the tessellation are the geometry. In
# practice they are the geometry *of the version of the code that reads them*: the Roche
# equipotential is a root solve, the rapid rotator's radii come from a cubic, the visibility
# clip is a sigmoid with a κ, and the element type can be Float32 or Float64. Writing the mesh
# means a geometry can be compared against what it was, rather than rebuilt and assumed equal.
# `load_star_geometry` returns both, and the GUI reports the largest disagreement when it
# rebuilds — see `shell_load_geometry`.
#
# WHAT IS NOT HERE, on purpose:
#
#   * `polyft`. That matrix is the polygon Fourier transform against ONE dataset's uv points —
#     a property of the data, not of the star — and it is the largest array in the struct by
#     far (nuv × npix complex). `setup_oi!(data, stars)` rebuilds it in the time it would take
#     to read.
#   * `polyflux`, which `setup_oi!` fills at the same time, is small and IS written: the radial
#     regularizers need it and their failure without it is an unhelpful BoundsError.
#   * the map. `save_surface_map` has it.
#   * free / fixed / tied and bounds. Those describe a FIT, not a surface.
#
# Layout:
#
#   Primary HDU   the primary component's `vertices_xyz`, with everything scalar in the header
#   VSPHERE …     the rest of that component's mesh, one image HDU each
#   PARAMS        one row per `star_params` field: name, value, and the Julia type it had
#   *2            the same names with a `2` suffix, for the secondary of a binary

# The mesh fields written as image HDUs, in the order they are written. `polyft` is absent for
# the reason given above; `t`, `npix` and the two type codes are scalars and live in the header.
const _GEOM_ARRAYS = (
    (:vertices_spherical, "VSPHERE"),
    (:normals,            "NORMALS"),
    (:proj_west,          "PROJWEST"),
    (:proj_north,         "PROJNRTH"),
    (:ldmap,              "LDMAP"),
    (:vis_weights,        "VISW"),
    (:sig_args,           "SIGARGS"),
    (:center_offsets,     "CENTOFF"),
    (:polyflux,           "POLYFLUX"),
)

"""
    save_star_geometry(path, star, params; nside_exp, kwargs...) -> path

Write a model's geometry to FITS: the mesh, the parameters that built it, and nothing about
what it is emitting.

`star` is a `stellar_geometry` and `params` the `star_params` NamedTuple it was built from.
`nside_exp` is the HEALPix level as `tessellation_healpix` takes it — log₂(nside).

For a binary, pass the companion as `star2`/`params2` together with how it is placed
(`place`, one of `:offset` or `:orbit`, and `offset` in `(West, North, toward-observer)` mas).
Both components' meshes are then written, the secondary's under the same HDU names with a `2`
suffix.

`tepochs`/`mjd` record the observation times the model was set up against, and `comment` is
free text for the primary header.

To get the geometry back:

```julia
g = load_star_geometry("model.fits")
g.star.vertices_xyz          # exactly the mesh that was written
tess  = tessellation_healpix(g.nside_exp)
star  = create_star(tess, g.params, g.t)      # and the rebuild, to compare against it
```
"""
function save_star_geometry(path::AbstractString, star, params::NamedTuple;
                            nside_exp::Integer,
                            tessellation::Symbol = :healpix,
                            star2 = nothing,
                            params2::Union{Nothing,NamedTuple} = nothing,
                            place::Symbol = :offset,
                            offset = (0.0, 0.0, 0.0),
                            tepochs::Union{Nothing,AbstractVector} = nothing,
                            mjd::Union{Nothing,AbstractVector} = nothing,
                            secondary::Bool = false,
                            comment::AbstractString = "")
    tessellation === :healpix || tessellation === :longlat ||
        throw(ArgumentError("tessellation must be :healpix or :longlat (got $tessellation)"))
    place === :offset || place === :orbit ||
        throw(ArgumentError("place must be :offset or :orbit (got $place)"))
    (star2 === nothing) == (params2 === nothing) ||
        throw(ArgumentError("a companion needs both star2 and params2"))

    hdr = FITSIO.FITSHeader(String[], Any[], String[])
    set!(k, v, c) = (hdr[k] = v; FITSIO.set_comment!(hdr, k, c))
    set!("ORIGIN",   "ROTIR",                 "written by ROTIR.save_star_geometry")
    set!("DATE",     Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS"), "file creation date")
    set!("CONTENT",  "geometry",              "the surface, NOT the map on it")
    set!("TESSEL",   String(tessellation),    "tessellation the mesh is on")
    set!("NSIDEEXP", Int(nside_exp),          "log2(nside), what tessellation_healpix takes")
    set!("TESSTYPE", Int(star.tessellation_type), "0 healpix, 1 longitude/latitude")
    set!("SURFTYPE", Int(star.surface_type),  "surface_type compute_radii branches on")
    set!("NPIX",     Int(star.npix),          "number of tessels")
    set!("NVISIBLE", Int(star.nquads_visible), "tessels facing the observer at this epoch")
    set!("EPOCHT",   Float64(star.t),         "epoch time the mesh was rotated to, days")
    set!("PRECISON", string(eltype(star.vertices_xyz)), "element type of the mesh")
    set!("SECONDRY", secondary,               "this component uses the secondary convention")
    set!("BINARY",   star2 !== nothing,       "a companion mesh is present")
    if star2 !== nothing
        set!("PLACE",  String(place),         "orbit or offset places the secondary")
        set!("OFFX",   Float64(offset[1]),    "secondary offset, West, mas")
        set!("OFFY",   Float64(offset[2]),    "secondary offset, North, mas")
        set!("OFFZ",   Float64(offset[3]),    "secondary offset, toward observer, mas")
        set!("NPIX2",  Int(star2.npix),       "tessels in the secondary mesh")
    end
    isempty(comment) || set!("NOTE", String(comment), "")

    FITSIO.FITS(String(path), "w") do f
        write(f, collect(star.vertices_xyz); header = hdr)
        _write_geom_arrays(f, star, "")
        write(f, _params_columns(params); name = "PARAMS")
        if star2 !== nothing
            write(f, collect(star2.vertices_xyz); name = "VXYZ2")
            _write_geom_arrays(f, star2, "2")
            write(f, _params_columns(params2); name = "PARAMS2")
        end
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

function _write_geom_arrays(f, star, suffix::AbstractString)
    for (fld, nm) in _GEOM_ARRAYS
        a = getfield(star, fld)
        # FITS has no zero-length image, and `polyflux` IS zero-length until `setup_oi!` has
        # run — a model built for a preview has never seen a uv point. An absent HDU is how
        # the reader learns that; it fills the field with an empty vector, which is what
        # `create_star` leaves there too.
        isempty(a) && continue
        write(f, collect(a); name = nm * suffix)
    end
    # The visible index set is Int64 and stays Int64: it indexes into the mesh, and a float
    # round trip through a 32-bit image would silently renumber a large tessellation.
    isempty(star.index_quads_visible) ||
        write(f, collect(Int64, star.index_quads_visible); name = "IDXVIS" * suffix)
    return f
end

# `star_params` as the three-column table `load_star_geometry` reads back. The TYPE column is
# what makes the round trip exact: `surface_type` and `ldtype` are `Int` and the geometry
# branches on them with `==`, and a Float64 3.0 compares equal while propagating as a float.
function _params_columns(params::NamedTuple)
    names = String[]; values = Float64[]; types = String[]
    for k in keys(params)
        v = getfield(params, k)
        v isa Real || continue
        # ASCII aliases for the three Unicode orbital elements; see `_param_ascii` in
        # src/surface_map_io.jl for why a FITS table cannot carry them verbatim.
        push!(names, _param_ascii(k)); push!(values, Float64(v))
        push!(types, v isa Integer ? "I" : v isa Bool ? "B" : "D")
    end
    return Dict("NAME" => names, "VALUE" => values, "TYPE" => types)
end

"""
    load_star_geometry(path) -> NamedTuple

Read back what [`save_star_geometry`](@ref) wrote.

Returns `star` (a `stellar_geometry` rebuilt from the stored arrays, with an empty `polyft` —
run `setup_oi!` if a χ² is wanted), `params`, `nside_exp`, `tessellation`, `t`, `secondary`,
and for a binary `star2`, `params2`, `place` and `offset`. `tepochs`/`mjd` come back when the
file carried them, `nothing` otherwise.

The returned `star` is the mesh AS SAVED, not a rebuild — that is the point of the file. To
check a rebuild against it, build one from `params` and compare `vertices_xyz`.
"""
function load_star_geometry(path::AbstractString)
    FITSIO.FITS(String(path), "r") do f
        hdr = read_header(f[1])
        g(k, d) = haskey(hdr, k) ? hdr[k] : d
        T = _precision_type(String(strip(String(g("PRECISON", "Float64")))))
        vxyz = Array{T}(read(f[1]))
        st   = Int(g("SURFTYPE", 0))
        tt   = Int(g("TESSTYPE", 0))
        # The header value stays Float64 in the returned `t` and only the struct's own field
        # narrows to the mesh type, which is what `create_star` does with the `t` it is handed.
        # The epoch a Float32 mesh CARRIES is already Float32 — 0.7 is stored as 0.69999999 —
        # so this preserves the star's own value rather than inventing precision it never had.
        tday = Float64(g("EPOCHT", 0.0))
        star = _read_star(f, "", vxyz, st, tt, T(tday), T)

        tep = nothing; mjd = nothing
        if _has_hdu(f, "EPOCHS")
            t = f["EPOCHS"]
            _has_col(t, "TDAY") && (tep = Float64.(read(t, "TDAY")))
            _has_col(t, "MJD")  && (mjd = Float64.(read(t, "MJD")))
            tep !== nothing && all(isnan, tep) && (tep = nothing)
            mjd !== nothing && all(isnan, mjd) && (mjd = nothing)
        end

        star2 = nothing; params2 = nothing
        if Bool(g("BINARY", false)) && _has_hdu(f, "VXYZ2")
            v2 = Array{T}(read(f["VXYZ2"]))
            star2 = _read_star(f, "2", v2, st, tt, T(tday), T)
            params2 = _read_params(f, "PARAMS2")
        end

        return (star = star, params = _read_params(f, "PARAMS"),
                nside_exp = Int(g("NSIDEEXP", round(Int, log2(sqrt(star.npix / 12))))),
                tessellation = Symbol(strip(String(g("TESSEL", "healpix")))),
                t = tday, secondary = Bool(g("SECONDRY", false)),
                star2 = star2, params2 = params2,
                place = Symbol(strip(String(g("PLACE", "offset")))),
                offset = (Float64(g("OFFX", 0.0)), Float64(g("OFFY", 0.0)),
                          Float64(g("OFFZ", 0.0))),
                tepochs = tep, mjd = mjd, header = hdr)
    end
end

# Only the two element types the mesh is ever written in. Anything else is a file from a
# future ROTIR and reading it as Float64 would be a quiet lie about its precision.
function _precision_type(s::AbstractString)
    s == "Float32" && return Float32
    s == "Float64" && return Float64
    throw(ArgumentError("unsupported mesh element type $(s) — expected Float32 or Float64"))
end

function _read_star(f, suffix::AbstractString, vxyz, st::Int, tt::Int, tval, ::Type{T}) where {T}
    npix = size(vxyz, 1)
    a(nm, dims...) = _has_hdu(f, nm * suffix) ? Array{T}(read(f[nm * suffix])) :
                     zeros(T, dims...)
    vsph = a("VSPHERE", npix, 5, 3)
    nrm  = a("NORMALS", npix, 3)
    pw   = a("PROJWEST", npix, 5)
    pn   = a("PROJNRTH", npix, 5)
    ld   = vec(a("LDMAP", npix))
    vw   = vec(a("VISW", npix))
    sa   = vec(a("SIGARGS", npix))
    co   = vec(a("CENTOFF", 3))
    pf   = vec(a("POLYFLUX", 0))
    idx  = _has_hdu(f, "IDXVIS" * suffix) ?
           vec(Int64.(read(f["IDXVIS" * suffix]))) : Int64[]
    # `polyft` is EMPTY by construction, not missing data: it belongs to a dataset's uv points
    # and `observables` takes the matrix-free route when it finds none.
    return stellar_geometry{T}(st, tt, npix, vxyz, vsph, nrm, idx, length(idx),
                               pw, pn, ld, vw, sa, co, pf,
                               Matrix{Complex{T}}(undef, 0, 0), T(tval))
end

function _read_params(f, name::AbstractString)
    _has_hdu(f, name) || return NamedTuple()
    t  = f[name]
    pn = read(t, "NAME"); pv = read(t, "VALUE")
    pt = _has_col(t, "TYPE") ? read(t, "TYPE") : fill("D", length(pn))
    nms = Symbol[]; vls = Any[]
    for j in eachindex(pn)
        push!(nms, _param_unicode(strip(String(pn[j]))))
        s = strip(String(pt[j]))
        push!(vls, s == "I" ? round(Int, pv[j]) :
                   s == "B" ? (pv[j] != 0)      : Float64(pv[j]))
    end
    return NamedTuple{Tuple(nms)}(Tuple(vls))
end
