# sky_export.jl — the model as SKY IMAGES, for use outside ROTIR (OITOOLS above all)
#
# ROTIR's own forward model never makes an image: it Fourier-transforms the projected facets
# exactly. An image cube is for a DIFFERENT program — OITOOLS computing visibilities at the
# epochs and channels of some dataset — so the one requirement that matters is that the cube,
# fed to OITOOLS, reproduces ROTIR's visibilities. Everything below is organised around that.
#
# THE LAYOUT IS OITOOLS'S, measured rather than assumed (`check_sky_cube` re-measures it):
#
#   cube[i, j, ichan, iepoch], with i increasing EAST and j increasing NORTH, the phase centre
#   at pixel (nx÷2+1, nx÷2+1) — where NFFT.jl puts the zero offset and so where OITOOLS's
#   `setup_nfft` does. Axes 3 and 4 follow `readoifits_multiepochs(...; polychromatic = true)`,
#   which returns `data[nwav, nepoch]`, so plane `[:, :, iλ, it]` goes with `data[iλ, it]`.
#
# A layout that is one transpose or one reflection off still gives the right V² — and closure
# phases of the wrong sign; one pixel off gives the right closure phases and a phase gradient
# on every complex visibility. Both were measured while writing this (the wrong candidates
# disagree by 1e-1, the right one by 4.5e-6), which is why the layout is ONE function.
#
# RENDERING. OITOOLS treats a pixel as a point sample. Two ways to make one:
#
#   :bandlimited (default) — the exact polygon Fourier transform on the FFT grid
#       (`polyft_nfft_image`), inverse-transformed. Its point-sample transform IS the exact
#       one on that grid, so OITOOLS reproduces ROTIR's visibilities to ~1e-5. The price is
#       ringing at the sharp limb: small negative pixels, which are part of what makes it exact.
#   :area — each pixel the exact facet-overlap area times intensity (`rasterize_polygon_image`):
#       positive and flux-conserving, the picture to LOOK at, but an area average is a box
#       filter, so OITOOLS sees visibilities attenuated by the pixel's sinc (measured 2.7e-3).
#
# The intensity is the forward model's own: `intensity(x, model, λ)` at each channel's λ, times
# limb darkening and the soft visibility weight, at the epoch's rotation phase — the `xw` of
# `fused_cvis_parts`, so the cube cannot disagree with ROTIR's χ² about what the star looks like.

"""
    sky_layout(img_north_west) -> Matrix

Convert the `(north, west)` image of `polyft_nfft_image` to OITOOLS' `(east, north)` layout,
centre at pixel `nx÷2+1`: a transpose, a reflection of the first axis, and a one-pixel shift
that puts the centre back after the reflection.
"""
sky_layout(img) = circshift(reverse(permutedims(img); dims = 1), (1, 0))

# ── components ────────────────────────────────────────────────────────────────────────────
#
# A sky cube draws one star or several. Each COMPONENT is a NamedTuple
#
#     (; params, x, tess, secondary = false, offsets = nothing)
#
# with `offsets` one `(west, north)` displacement in mas per epoch (`nothing` = at the origin).
# That is exactly ROTIR's binary model: `binary_cvis` multiplies the companion's polygon
# transform by `binary_phase_shift(uv, west, north)`, which IS the transform of the same
# polygons translated by (west, north) — so drawing the translated facets into one band-limited
# image reproduces the binary's visibilities exactly, with no special case. Mutual occultation
# is not drawn (the GUI's binary χ² does not apply it either); see `check_binary_overlap`.

_component(c) = (; params = c.params, x = Float64.(c.x), tess = c.tess,
                 secondary = hasproperty(c, :secondary) ? c.secondary : false,
                 offsets = hasproperty(c, :offsets) ? c.offsets : nothing)

_stars_of(c, ts) = c.secondary === false || c.secondary === nothing ?
    create_star_multiepochs(c.tess, c.params, ts) :
    create_star_multiepochs(c.tess, c.params, ts; secondary = c.secondary)

_offset(c, it) = c.offsets === nothing ? (0.0, 0.0) : (Float64(c.offsets[it][1]), Float64(c.offsets[it][2]))

"""
    sky_grid(params, tess, tepochs; data, oversample = 4, fov_factor = 6.0, min_nx = 256, pixsize, nx)
    sky_grid(components, tepochs; ...)
        -> (; pixsize, nx, diameter, extent, mb_per_plane, reason)

A pixel size and image size for a sky cube, with the reason for each.

**Pixel size** (mas): `1/(2ρmax)/oversample`, with `ρmax` the longest baseline in `data`
(cycles per radian). Without data, 1/64 of the largest projected stellar diameter.

**Image size**: the whole system plus `(fov_factor − 1)` diameters of the largest star (for one
star, `fov_factor` diameters), at least `min_nx`, rounded up to an even 2·3·5-smooth size.

On the ζ And coverage, OITOOLS reproduces ROTIR's visibilities to 1.3e-4 with a field of 2
diameters, 3e-5 with 4 (closure phases 0.28° → 0.05° rms) and about 5e-6 with 12; 256 pixels
is about 11 diameters there.

A binary's field grows with its separation while the pixel size is set by the baselines;
`mb_per_plane` gives the memory per Float32 plane. `pixsize` or `nx` given explicitly are
used as is.
"""
sky_grid(params, tess, tepochs; kw...) = sky_grid([(; params, x = Float64[], tess)], tepochs; kw...)

function sky_grid(components::AbstractVector, tepochs; data = nothing, oversample::Real = 4,
                  fov_factor::Real = 6.0, min_nx::Integer = 256, pixsize = nothing, nx = nothing)
    ts = collect(Float64, tepochs)
    comps = [_component(merge((; x = Float64[]), c)) for c in components]
    diameter = 0.0; extent = 0.0
    for c in comps
        for (it, st) in enumerate(_stars_of(c, ts))
            ix = st.index_quads_visible
            pw = st.proj_west[ix, :]; pn = st.proj_north[ix, :]
            diameter = max(diameter, 2maximum(abs.(vcat(vec(pw), vec(pn)))))
            ox, oy = _offset(c, it)
            extent = max(extent, maximum(abs.(vec(pw) .+ ox)), maximum(abs.(vec(pn) .+ oy)))
        end
    end
    reason = String[]
    if pixsize === nothing
        if data === nothing
            pixsize = diameter / 64
            push!(reason, "pixel = diameter/64 (no data to set it from)")
        else
            ds = data isa AbstractArray ? vec(data) : [data]
            ρ = maximum(d -> maximum(sqrt.(sum(abs2, d.uv; dims = 1))), ds)   # cycles / rad
            pixsize = 206264806.2 / (2 * ρ) / oversample                       # mas
            push!(reason, @sprintf("pixel = Nyquist/%g at the longest baseline", oversample))
        end
    end
    if nx === nothing
        fov = 2extent + (fov_factor - 1) * diameter
        n = ceil(Int, fov / pixsize)
        nx = nextprod([2, 3, 5], max(n, min_nx))
        isodd(nx) && (nx = nextprod([2, 3, 5], nx + 1))
        push!(reason, length(comps) == 1 ?
              @sprintf("field = %.3g × the projected diameter (%.4g mas)", nx * pixsize / diameter, diameter) :
              @sprintf("field = the system (%.4g mas across) + %.3g diameters", 2extent,
                       (nx * pixsize - 2extent) / diameter))
    end
    mb = nx^2 * 4 / 2^20
    return (; pixsize = Float64(pixsize), nx = Int(nx), diameter, extent,
            mb_per_plane = mb, reason = join(reason, "; "))
end

# The physical λ dependence the dimensionless conversions leave out, so planes at different
# wavelengths can share one flux scale: B_λ ∝ λ⁻⁵/(e^{c₂/λT} − 1) for Planck (`intensity`
# returns only the second factor) and λ⁻⁴·T in the Rayleigh–Jeans limit (`:linear` is T).
_lambda_factor(model::Symbol, λ) = model === :planck ? λ^-5 : λ^-4

"""
    render_sky_cube(params, xmap, tess; tepochs, wavelengths, pixsize, nx,
                    rendering = :bandlimited, intensity_model = :linear,
                    normalize = :plane, secondary = nothing, mjd = nothing,
                    bandwidths = nothing, T = Float32) -> (cube, meta)
    render_sky_cube(components; tepochs, wavelengths, ...) -> (cube, meta)

Render the surface `xmap` (on `tess`, with geometry `params`), or several components, as sky
images: an `(nx, nx, nλ, nt)` array in OITOOLS' layout (see [`sky_layout`](@ref)), one plane
per wavelength and epoch. `nx` must be even.

Each component is a NamedTuple `(; params, x, tess, secondary = false, offsets = nothing)`;
`offsets` holds one `(west, north)` displacement in mas per epoch (`nothing`: at the origin).
A binary is two components, the companion offset by its orbit. Mutual occultation is not drawn.

  * `tepochs` are ROTIR epoch times (the same numbers `create_star_multiepochs` takes); `mjd`,
    when given, is recorded alongside them.
  * `wavelengths` in metres. With `intensity_model = :linear` (Rayleigh–Jeans) every channel
    has the same SHAPE; with `:planck` the contrast — and a binary's flux ratio — change with λ.
  * `rendering`: `:bandlimited` (default) is the exact polygon Fourier transform on the FFT
    grid, inverse-transformed; OITOOLS reproduces ROTIR's visibilities from it to about 1e-5,
    and it has small negative pixels near the limb. `:area` gives each pixel the exact facet
    overlap times intensity: positive and flux-conserving, for display, with visibilities
    attenuated by the pixel's sinc (2.7e-3 measured).
  * `normalize`: `:plane`, every plane sums to 1 (as OITOOLS imaging expects); or `:relative`,
    relative flux across channels and epochs (the λ⁻⁵ or λ⁻⁴ factor of the intensity is
    restored), scaled so the brightest plane sums to 1.
  * `bandwidths` (m), optional: recorded as `EFF_BAND` in the `CHANNELS` table, like
    `OI_WAVELENGTH`. The planes themselves are monochromatic at `wavelengths`.

`meta` carries everything [`write_sky_cube`](@ref) records and [`check_sky_cube`](@ref) needs:
pixel size, wavelengths, epochs, per-plane flux before normalisation, rendering, intensity
model, the components and their parameters (`meta.params` is the first component's).
"""
render_sky_cube(params, xmap, tess; secondary = nothing, kw...) =
    render_sky_cube([(; params, x = xmap, tess, secondary = something(secondary, false))]; kw...)

function render_sky_cube(components::AbstractVector; tepochs, wavelengths, pixsize::Real,
                         nx::Integer, rendering::Symbol = :bandlimited,
                         intensity_model::Symbol = :linear, normalize::Symbol = :plane,
                         mjd = nothing, bandwidths = nothing,
                         T::Type{<:AbstractFloat} = Float32)
    rendering in (:bandlimited, :area) ||
        throw(ArgumentError("rendering must be :bandlimited or :area, not :$(rendering)"))
    normalize in (:plane, :relative) ||
        throw(ArgumentError("normalize must be :plane or :relative, not :$(normalize)"))
    intensity_model in (:linear, :planck) ||
        throw(ArgumentError("intensity_model must be :linear (Rayleigh–Jeans) or :planck"))
    iseven(nx) || throw(ArgumentError("nx must be even (the centre is pixel nx÷2+1)"))
    λs = collect(Float64, wavelengths); ts = collect(Float64, tepochs)
    isempty(λs) && throw(ArgumentError("no wavelengths"))
    isempty(ts) && throw(ArgumentError("no epochs"))
    comps = [_component(c) for c in components]
    isempty(comps) && throw(ArgumentError("no components"))
    for c in comps
        c.offsets === nothing || length(c.offsets) == length(ts) ||
            throw(ArgumentError("a component's offsets must have one (west, north) per epoch"))
    end
    stars = [_stars_of(c, ts) for c in comps]
    cube = zeros(T, nx, nx, length(λs), length(ts))
    flux = zeros(Float64, length(λs), length(ts))
    # The quadrature is chosen for the HIGHEST frequency on the grid, 1/(2·pixsize) on each
    # axis, so it resolves every frequency the image is made of (`quadrature_for`).
    kg = [1 / (2 * pixsize)]
    for it in eachindex(ts)
        # Every component's visible facets, translated to where it is at this epoch, as ONE
        # polygon set: the transform is linear, so this is the sum of the components.
        pws = Matrix{Float64}[]; pns = Matrix{Float64}[]; ws = Vector{Float64}[]; ids = Vector{Int}[]
        for (k, c) in enumerate(comps)
            st = stars[k][it]; ix = st.index_quads_visible
            ox, oy = _offset(c, it)
            push!(pws, Matrix{Float64}(st.proj_west[ix, :]) .+ ox)
            push!(pns, Matrix{Float64}(st.proj_north[ix, :]) .+ oy)
            push!(ws, Float64.(st.vis_weights[ix] .* st.ldmap[ix])); push!(ids, ix)
        end
        pw = reduce(vcat, pws); pn = reduce(vcat, pns)
        # Quadrature for the widest QUAD, not the widest system: the rule resolves the phase
        # across one facet, and translating a facet does not change that.
        ng, ns = quadrature_for(pw .- pw[:, 1:1], pn .- pn[:, 1:1], kg, kg)
        for (iλ, λ) in enumerate(λs)
            xw = reduce(vcat, [(intensity_model === :linear ? c.x : intensity(c.x, intensity_model, λ))[ids[k]] .* ws[k]
                               for (k, c) in enumerate(comps)])
            plane = rendering === :bandlimited ?
                sky_layout(polyft_nfft_image(pw, pn, xw, pixsize, nx; ngauss = ng, nsub = ns,
                                             T = Float64, fftflags = FFTW.ESTIMATE)) :
                permutedims(rasterize_polygon_image(-pw, pn, xw, pixsize, nx; T = Float64))
            f = sum(plane)
            flux[iλ, it] = f * _lambda_factor(intensity_model, λ)
            cube[:, :, iλ, it] .= T.(normalize === :plane ? plane ./ f : plane .* _lambda_factor(intensity_model, λ))
        end
    end
    if normalize === :relative
        cube ./= T(maximum(flux))
    end
    meta = (; pixsize = Float64(pixsize), nx = Int(nx), wavelengths = λs, tepochs = ts,
            mjd = mjd === nothing ? fill(NaN, length(ts)) : collect(Float64, mjd),
            bandwidths = bandwidths === nothing ? nothing : collect(Float64, bandwidths),
            flux, rendering, intensity_model, normalize, params = comps[1].params,
            components = [(; params = c.params, secondary = c.secondary, offsets = c.offsets) for c in comps])
    return cube, meta
end

"""
    write_sky_cube(path, cube, meta; comment = "") -> String

Write a cube from [`render_sky_cube`](@ref) as FITS: the `(nx, nx, nλ, nt)` primary array with
a celestial WCS, plus `CHANNELS` and `EPOCHS` tables and the model's parameters (`PARAMS`).

The first axis increases East (`CDELT1 > 0`) and the reference pixel is `nx÷2+1`, the layout
OITOOLS uses (OITOOLS' own `writefits` writes `CDELT1 < 0` and `CRPIX = nx/2`). Pixel scales
are in degrees with `CUNIT` set; OITOOLS' `fits_pixsize` reads them in mas.

Axes 3 and 4 get a linear WCS (`WAVE` in metres, `TIME` as MJD) when the channels or epochs are
evenly spaced to 1 %; otherwise `WAVEIRR`/`TIMEIRR` are set and the `CHANNELS`/`EPOCHS` tables
give the values.
"""
function write_sky_cube(path::AbstractString, cube::AbstractArray{<:Real,4}, meta;
                        comment::AbstractString = "")
    nx = size(cube, 1); nλ = size(cube, 3); nt = size(cube, 4)
    hdr = FITSIO.FITSHeader(String[], Any[], String[])
    set!(k, v, c) = (hdr[k] = v; FITSIO.set_comment!(hdr, k, c))
    pixdeg = meta.pixsize / 3.6e6
    set!("ORIGIN",   "ROTIR", "written by ROTIR.write_sky_cube")
    set!("DATE",     Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS"), "file creation date")
    set!("CTYPE1",   "RA---SIN", "axis 1 increases EAST (OITOOLS layout)")
    set!("CUNIT1",   "deg",      "")
    set!("CDELT1",   pixdeg,     "positive: RA increases with pixel index")
    set!("CRPIX1",   Float64(nx ÷ 2 + 1), "phase centre: NFFT zero offset, pixel nx/2+1")
    set!("CRVAL1",   0.0,        "offset from the star's centre")
    set!("CTYPE2",   "DEC--SIN", "axis 2 increases NORTH")
    set!("CUNIT2",   "deg",      "")
    set!("CDELT2",   pixdeg,     "")
    set!("CRPIX2",   Float64(nx ÷ 2 + 1), "phase centre")
    set!("CRVAL2",   0.0,        "")
    set!("PIXSIZE",  meta.pixsize, "pixel size in mas")
    regular(v) = length(v) < 2 || (s = (v[end] - v[1]) / (length(v) - 1);
                                   s != 0 && maximum(abs.(diff(v) .- s)) <= 0.01 * abs(s))
    λ = meta.wavelengths
    set!("CTYPE3", "WAVE", "wavelength; the CHANNELS table is the authority")
    set!("CUNIT3", "m", "")
    if regular(λ)
        set!("CRPIX3", 1.0, ""); set!("CRVAL3", λ[1], "centre wavelength of channel 1")
        set!("CDELT3", nλ > 1 ? (λ[end] - λ[1]) / (nλ - 1) : 0.0, "metres per channel")
    else
        set!("WAVEIRR", true, "channels unevenly spaced: see the CHANNELS table")
    end
    mjd = meta.mjd
    set!("CTYPE4", "TIME", "epoch; the EPOCHS table is the authority")
    set!("CUNIT4", "d", "")
    set!("TIMESYS", "UTC", "MJD as in the OIFITS data")
    if all(isfinite, mjd) && regular(mjd)
        set!("MJDREF", 0.0, "TIME axis values are MJD")
        set!("CRPIX4", 1.0, ""); set!("CRVAL4", mjd[1], "MJD of epoch 1")
        set!("CDELT4", nt > 1 ? (mjd[end] - mjd[1]) / (nt - 1) : 0.0, "days per epoch")
    else
        set!("TIMEIRR", true, "epochs unevenly spaced or without MJD: see the EPOCHS table")
    end
    set!("BUNIT",    meta.normalize === :plane ? "plane-normalized" : "relative flux",
         meta.normalize === :plane ? "each plane sums to 1" :
                                     "physical relative flux; brightest plane sums to 1")
    set!("NORMALIZ", meta.normalize === :plane ? "PLANE" : "RELATIVE",
         "PLANE: each plane sums to 1; RELATIVE: SED and light curve kept")
    set!("INTMODEL", meta.intensity_model === :linear ? "Rayleigh-Jeans" : "Planck",
         "temperature -> intensity conversion")
    set!("RENDER",   String(meta.rendering),
         meta.rendering === :bandlimited ? "exact polygon FT on the grid; rings at the limb" :
                                           "exact facet-overlap area per pixel")
    isempty(comment) || set!("NOTE", String(comment), "")
    comps = hasproperty(meta, :components) ? meta.components : [(; params = meta.params, offsets = nothing)]
    ncomp = length(comps)
    set!("NCOMP", ncomp, ncomp == 1 ? "one star" : "components drawn into one image; see EPOCHS")
    # PARAMS: every numeric parameter, prefixed `c2.`, `c3.` … for the components after the
    # first, so a binary's two inclinations stay distinguishable.
    names = String[]; values = Float64[]
    for (k, c) in enumerate(comps), key in keys(c.params)
        v = getfield(c.params, key)
        v isa Real || continue
        push!(names, (k == 1 ? "" : "c$(k).") * _param_ascii(key)); push!(values, Float64(v))
    end
    # EPOCHS: MJD and ROTIR's epoch time, plus where each component after the first sits —
    # `OFF<k>_W` / `OFF<k>_N` in mas, west and north of the reference pixel (the first
    # component's centre). That is the offset `binary_phase_shift` takes.
    ecols = Dict{String,Any}("MJD" => collect(Float64, mjd), "TDAY" => collect(Float64, meta.tepochs))
    for (k, c) in enumerate(comps)
        k == 1 && continue
        offs = c.offsets === nothing ? fill((0.0, 0.0), nt) : c.offsets
        ecols["OFF$(k)_W"] = Float64[o[1] for o in offs]
        ecols["OFF$(k)_N"] = Float64[o[2] for o in offs]
    end
    ccols = Dict{String,Any}("EFF_WAVE" => collect(Float64, λ),
                             "FLUX" => vec(sum(meta.flux; dims = 2)) ./ max(nt, 1))
    hasproperty(meta, :bandwidths) && meta.bandwidths !== nothing &&
        (ccols["EFF_BAND"] = collect(Float64, meta.bandwidths))
    FITSIO.FITS(String(path), "w") do f
        write(f, collect(cube); header = hdr)
        write(f, ccols; name = "CHANNELS")
        write(f, ecols; name = "EPOCHS")
        write(f, Dict("NAME" => names, "VALUE" => values); name = "PARAMS")
    end
    return String(path)
end

"""
    read_sky_cube(path) -> (; cube, pixsize, wavelengths, mjd, tepochs)

The cube and its axes back from a file [`write_sky_cube`](@ref) wrote.
"""
function read_sky_cube(path::AbstractString)
    FITSIO.FITS(String(path), "r") do f
        cube = read(f[1]); hdr = read_header(f[1])
        ch = f["CHANNELS"]; ep = f["EPOCHS"]
        return (; cube, pixsize = Float64(hdr["PIXSIZE"]),
                wavelengths = read(ch, "EFF_WAVE"), mjd = read(ep, "MJD"),
                tepochs = read(ep, "TDAY"))
    end
end

"""
    check_sky_cube(cube, meta, xmap, tess, data) -> NamedTuple
    check_sky_cube(cube, meta, components, data) -> NamedTuple

Compare the visibilities OITOOLS computes from the cube with ROTIR's own, on the uv coverage
of `data` and in units of its error bars.

`data` is one `OIdata` per epoch (a vector, matched to the cube's epochs in order) or a
`(nλ, nt)` matrix as `readoifits_multiepochs(...; polychromatic = true)` returns. For each
plane, OITOOLS's `setup_nfft`/`image_to_vis` and ROTIR's exact model are both evaluated at the
same uv with the same intensity at the plane's λ. For several components the exact model is ROTIR's binary one,
`(Σ F_k · binary_phase_shift(uv, offset_k)) / Σ flux_k`. With a vector of single-band epochs
the plane whose λ is closest to the data's mean wavelength is used.

Returns `chi2dist`, Σ(Δ/σ)²/n over V², T3amp and T3φ (the χ² the image adds), and `rel_dV`,
`max_dV2` and `rms_dT3phi_deg`, per plane and the worst over all planes.
"""
check_sky_cube(cube::AbstractArray{<:Real,4}, meta, xmap, tess, data) =
    check_sky_cube(cube, meta, [(; params = meta.params, x = xmap, tess,
                                  secondary = hasproperty(meta, :components) ? meta.components[1].secondary : false)], data)

function check_sky_cube(cube::AbstractArray{<:Real,4}, meta, components::AbstractVector, data)
    nλ, nt = size(cube, 3), size(cube, 4)
    comps = [_component(c) for c in components]
    stars = [_stars_of(c, meta.tepochs) for c in comps]
    rows = NamedTuple[]
    for it in 1:nt
        planes = data isa AbstractMatrix ? [(iλ, data[iλ, it]) for iλ in 1:min(nλ, size(data, 1))] :
                 [(argmin(abs.(meta.wavelengths .- Statistics.mean(data[it].uv_lam))), data[it])]
        for (iλ, d) in planes
            λ = meta.wavelengths[iλ]
            Tp = eltype(d.uv)
            ft = OITOOLS.setup_nfft(d, meta.nx, meta.pixsize)
            cO = ComplexF64.(OITOOLS.image_to_vis(Tp.(cube[:, :, iλ, it]), ft[1]))
            Fsum = zeros(ComplexF64, size(d.uv, 2)); flsum = 0.0
            for (k, c) in enumerate(comps)
                F, fl = fused_cvis_parts(c.x, stars[k][it], d; intensity_model = meta.intensity_model,
                                         band = meta.intensity_model === :linear ? nothing : λ)
                ox, oy = _offset(c, it)
                Fsum .+= ComplexF64.(F) .* (ox == 0 && oy == 0 ? 1.0 :
                                            ComplexF64.(binary_phase_shift(Float64.(d.uv), ox, oy)))
                flsum += Float64(fl)
            end
            cR = Fsum ./ flsum
            v2O = abs2.(cO[d.indx_v2]); v2R = abs2.(cR[d.indx_v2])
            t3O = cO[d.indx_t3_1] .* cO[d.indx_t3_2] .* cO[d.indx_t3_3]
            t3R = cR[d.indx_t3_1] .* cR[d.indx_t3_2] .* cR[d.indx_t3_3]
            dφ = rad2deg.(mod.(angle.(t3O) .- angle.(t3R) .+ π, 2π) .- π)
            r = vcat((v2O .- v2R) ./ d.v2_err, (abs.(t3O) .- abs.(t3R)) ./ d.t3amp_err,
                     dφ ./ d.t3phi_err)
            push!(rows, (; iλ, it, chi2dist = sum(abs2, r) / max(length(r), 1),
                         rel_dV = LinearAlgebra.norm(cO .- cR) / LinearAlgebra.norm(cR),
                         max_dV2 = isempty(v2O) ? 0.0 : maximum(abs.(v2O .- v2R)),
                         rms_dT3phi_deg = isempty(dφ) ? 0.0 : sqrt(Statistics.mean(abs2, dφ))))
        end
    end
    worst = (; chi2dist = maximum(r.chi2dist for r in rows), rel_dV = maximum(r.rel_dV for r in rows),
             max_dV2 = maximum(r.max_dV2 for r in rows),
             rms_dT3phi_deg = maximum(r.rms_dT3phi_deg for r in rows))
    return (; planes = rows, worst)
end
