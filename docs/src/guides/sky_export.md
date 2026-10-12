# Exporting sky images for OITOOLS

ROTIR's forward model never makes an image: it Fourier-transforms the projected facets
exactly. A sky image is for a **different program** — OITOOLS computing visibilities at the
epochs and channels of some dataset — so the one requirement that matters is that the image,
fed to OITOOLS, reproduces ROTIR's own visibilities. Everything below is organised around
that, and [`check_sky_cube`](@ref) measures it.

## The cube

[`render_sky_cube`](@ref) returns an `(nx, nx, nλ, nt)` array — one plane per wavelength
channel and epoch, no axes, no graticule — in OITOOLS's own layout:

* `cube[i, j, iλ, it]`, `i` increasing **east**, `j` increasing **north**;
* the phase centre at pixel `(nx÷2+1, nx÷2+1)`, where NFFT.jl (and so OITOOLS's `setup_nfft`)
  puts the zero offset;
* axes 3 and 4 in the order `readoifits_multiepochs(...; polychromatic = true)` returns
  `data[nwav, nepoch]`, so plane `[:, :, iλ, it]` goes with `data[iλ, it]`.

The layout is one function, [`sky_layout`](@ref), because the wrong candidates fail quietly:
a transpose or a reflection still gives the right V² with closure phases of the wrong sign,
and a one-pixel shift gives the right closure phases with a phase gradient on every complex
visibility.

```julia
tep = [minimum(d.v2_mjd) for d in data]                 # epoch times, as for create_star_multiepochs
g   = sky_grid(params, tess, tep; data = data)           # pixel size and field, with reasons
cube, meta = render_sky_cube(params, x, tess; tepochs = tep, wavelengths = [1.55e-6, 1.65e-6],
                             pixsize = g.pixsize, nx = g.nx, mjd = tep)
write_sky_cube("star_sky.fits", cube, meta)
check_sky_cube(cube, meta, x, tess, data).worst           # rel_dV, rms_dT3phi_deg, chi2dist
```

## Pixel and field

[`sky_grid`](@ref) picks both and says why (`g.reason`):

* **pixel**: `1/(2ρmax)/4` — a quarter of Nyquist at the longest baseline in the data;
* **field**: the whole system plus a margin, by default six stellar diameters for one star,
  never below 256 pixels. Six, not two, because a band-limited image rings beyond the limb and
  a tight field folds that ringing back in (measured on ζ And: 1.3e-4 → 3e-5).

Large images are fine — 256×256 is ordinary in interferometry — and `g.mb_per_plane` reports
the cost per plane.

## Rendering

| `rendering` | what a pixel is | agrees with ROTIR in OITOOLS | use |
|---|---|---|---|
| `:bandlimited` (default) | exact polygon transform on the FFT grid, inverted | ~1e-5 | computing visibilities |
| `:area` | exact facet-overlap area × intensity | ~3e-3 (box-filtered) | looking at |

The band-limited image has small negative pixels at the limb; they are part of what makes it
exact. The area image is positive and flux-conserving, but an area average is a box filter,
so its visibilities carry the pixel's sinc.

## Intensity and normalisation

Each channel's intensity is the model's own conversion at that channel's wavelength
(`intensity_model = :linear` for Rayleigh–Jeans, `:planck` for Planck — see
[From temperature to intensity](@ref)). Under Rayleigh–Jeans every channel has the same
shape; under Planck the contrast, and a binary's flux ratio, change with λ.

`normalize = :plane` (default) makes every plane sum to 1, which is what OITOOLS imaging
expects. `normalize = :relative` keeps the physical flux across channels and epochs — the
``\lambda^{-5}`` or ``\lambda^{-4}`` the dimensionless conversions drop is put back — scaled
so the brightest plane sums to 1, so the cube carries the spectrum and the light curve.

Doppler-shifted lines are not drawn: each plane is the continuum surface at its wavelength.

## Binaries

A cube can hold several components. Each is a NamedTuple
`(; params, x, tess, secondary = false, offsets = nothing)`, with `offsets` one
`(west, north)` displacement in mas per epoch:

```julia
comps = [(; params = p1, x = x1, tess),
         (; params = p2, x = x2, tess, secondary = true, offsets = offs)]
g = sky_grid(comps, tep; data = data)                   # the field grows with the separation
cube, meta = render_sky_cube(comps; tepochs = tep, wavelengths = [1.6e-6],
                             pixsize = g.pixsize, nx = g.nx)
```

This is exactly ROTIR's binary model: drawing the companion's facets translated by its offset
is the same as multiplying its transform by `binary_phase_shift`, so the binary's
visibilities are reproduced with no special case. Two limbs ring instead of one, so the
agreement is ~1e-4 rather than ~1e-5. Mutual occultation is not drawn (the binary χ² does not
apply it either). The pixel is still fixed by the baselines while the field grows with the
separation, so a wide pair costs `nx²` per plane — `sky_grid` reports it and the caller
decides.

## The FITS file

[`write_sky_cube`](@ref) writes the cube as the primary array with a celestial WCS
(`RA---SIN`/`DEC--SIN`, degrees, `CDELT1 > 0` because the first axis increases east,
`CRPIX = nx÷2+1`), `WAVE` on axis 3 and `TIME` (MJD) on axis 4 when those are evenly spaced —
otherwise the header says so and the tables are the authority. Three binary tables follow:
`CHANNELS` (`EFF_WAVE`, `FLUX`, `EFF_BAND`), `EPOCHS` (`MJD`, `TDAY`, and the companion's
`OFF2_W`/`OFF2_N`) and `PARAMS` (the model, companion fields prefixed `c2.`). OITOOLS's
`fits_pixsize` reads the pixel scale back in mas; [`read_sky_cube`](@ref) reads everything.

In the GUI, **Export sky FITS…** opens a dialog with the same choices: model or a
reconstruction, data epochs or MJDs or rotation phases, data channels or custom wavelengths,
pixel and field (with an *auto* button for `sky_grid`'s choice), rendering, normalisation,
and an optional check against ROTIR on the loaded data.
