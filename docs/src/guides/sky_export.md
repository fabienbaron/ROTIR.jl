# Exporting sky images for OITOOLS

ROTIR computes visibilities from the projected facets directly, without an image. A sky cube
lets OITOOLS compute the same visibilities at the epochs and channels of a dataset;
[`check_sky_cube`](@ref) measures how closely they agree.

## The cube

[`render_sky_cube`](@ref) returns an `(nx, nx, nλ, nt)` array, one plane per wavelength
channel and epoch, in OITOOLS's layout:

* `cube[i, j, iλ, it]`, `i` increasing east, `j` increasing north;
* the phase centre at pixel `(nx÷2+1, nx÷2+1)`, the zero offset of NFFT.jl and OITOOLS's
  `setup_nfft`;
* axes 3 and 4 in the order of `readoifits_multiepochs(...; polychromatic = true)`, which
  returns `data[nwav, nepoch]`: plane `[:, :, iλ, it]` goes with `data[iλ, it]`.

[`sky_layout`](@ref) applies this layout. A transpose or reflection would give the right V²
and closure phases of the wrong sign; a one-pixel shift would add a phase gradient to every
complex visibility.

```julia
tep = [minimum(d.v2_mjd) for d in data]                 # epoch times, as for create_star_multiepochs
g   = sky_grid(params, tess, tep; data = data)           # pixel size and field
cube, meta = render_sky_cube(params, x, tess; tepochs = tep, wavelengths = [1.55e-6, 1.65e-6],
                             pixsize = g.pixsize, nx = g.nx, mjd = tep)
write_sky_cube("star_sky.fits", cube, meta)
check_sky_cube(cube, meta, x, tess, data).worst           # rel_dV, rms_dT3phi_deg, chi2dist
```

## Pixel and field

[`sky_grid`](@ref) chooses both and gives the reason in `g.reason`:

* pixel: `1/(2ρmax)/oversample`, with `oversample = 4` (a quarter of Nyquist at the longest
  baseline); without data, the diameter / 64;
* field: the whole system plus a margin, `fov_factor = 6` stellar diameters for one star, at
  least `min_nx = 256` pixels. The band-limited image rings beyond the limb; a smaller field
  folds the ringing back in (ζ And: error 1.3e-4 at two diameters, 3e-5 at six).

`g.mb_per_plane` is the memory per plane (Float32). `pixsize` and `nx` override the choice.

## Rendering

| `rendering` | what a pixel is | agrees with ROTIR in OITOOLS | use |
|---|---|---|---|
| `:bandlimited` (default) | exact polygon transform on the FFT grid, inverted | ~1e-5 | computing visibilities |
| `:area` | exact facet-overlap area × intensity | ~3e-3 (box-filtered) | display |

The band-limited image has small negative pixels at the limb. The area image is positive and
conserves flux; as a box filter, it multiplies the visibilities by the pixel's sinc.

## Intensity and normalisation

Each channel's intensity uses the model's conversion at that channel's wavelength:
`intensity_model = :linear` (Rayleigh–Jeans, default) or `:planck` (see
[From temperature to intensity](@ref)). Under Rayleigh–Jeans every channel has the same
shape; under Planck the contrast and a binary's flux ratio change with λ.

`normalize = :plane` (default) makes every plane sum to 1, as OITOOLS imaging expects.
`normalize = :relative` keeps the relative flux across channels and epochs (restoring the
``\lambda^{-5}`` or ``\lambda^{-4}`` factor), scaled so the brightest plane sums to 1; the
cube then carries the spectrum and the light curve.

Lines are not drawn: each plane is the continuum surface at its wavelength.

## Binaries

A cube can hold several components. Each is a NamedTuple
`(; params, x, tess, secondary = false, offsets = nothing)`, with `offsets` one
`(west, north)` displacement in mas per epoch:

```julia
comps = [(; params = p1, x = x1, tess),
         (; params = p2, x = x2, tess, secondary = true, offsets = offs)]
g = sky_grid(comps, tep; data = data)                   # field includes the separation
cube, meta = render_sky_cube(comps; tepochs = tep, wavelengths = [1.6e-6],
                             pixsize = g.pixsize, nx = g.nx)
```

Translating the companion's facets by its offset is equivalent to ROTIR's
`binary_phase_shift`, so the cube reproduces the binary model. With two limbs ringing, the
agreement is ~1e-4. Mutual occultation is not drawn (the binary χ² does not include it
either). The pixel is set by the baselines and the field by the separation, so a wide pair
needs a large `nx`; `g.mb_per_plane` gives the memory.

## The FITS file

[`write_sky_cube`](@ref) writes the cube as the primary array with a celestial WCS
(`RA---SIN`/`DEC--SIN`, degrees, `CDELT1 > 0` because the first axis increases east,
`CRPIX = nx÷2+1`), `WAVE` on axis 3 and `TIME` (MJD) on axis 4. When channels or epochs are
not evenly spaced (to 1 %), `WAVEIRR` or `TIMEIRR` is set and the tables hold the values.
Three binary tables follow: `CHANNELS` (`EFF_WAVE`, `FLUX`, and `EFF_BAND` if `bandwidths`
was given), `EPOCHS` (`MJD`, `TDAY`, and `OFF2_W`/`OFF2_N` for a companion) and `PARAMS` (the
numeric model parameters, companion fields prefixed `c2.`). OITOOLS's `fits_pixsize` reads the
pixel scale in mas. [`read_sky_cube`](@ref) returns the cube, `pixsize`, `wavelengths`, `mjd`
and `tepochs`.

In the GUI, **Export sky FITS…** opens a dialog with the same options: model or
reconstruction; data epochs, MJDs or rotation phases; data channels or custom wavelengths;
pixel and field (*auto* uses `sky_grid`); rendering; normalisation; and an optional check
against the loaded data.
