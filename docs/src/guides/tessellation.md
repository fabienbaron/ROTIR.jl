# Tessellation

ROTIR discretizes the stellar surface into quadrilateral pixels using one of
two tessellation schemes. The choice affects pixel uniformity, neighbor
structure, and multi-resolution capabilities.

## HEALPix scheme

The HEALPix (Hierarchical Equal Area isoLatitude Pixelization) scheme divides
the sphere into equal-area pixels arranged in a nested hierarchy. This is the
recommended scheme for most applications.

```julia
n = 3                                   # resolution parameter
tessels = tessellation_healpix(n)       # nside = 2^n = 8, npix = 768
```

The resolution parameter `n` determines the number of pixels:

| `n` | nside | npix | Typical use |
|-----|-------|------|-------------|
| 1 | 2 | 48 | Quick tests |
| 2 | 4 | 192 | Low resolution |
| 3 | 8 | 768 | Standard reconstruction |
| 4 | 16 | 3072 | High resolution |
| 5 | 32 | 12288 | Very high resolution |

Conversion functions:
- `nside2npix(nside)` -- returns `12 * nside^2`
- `npix2n(npix)` -- returns the resolution parameter `n`

### Pixel structure

Each pixel has 5 associated points: 4 corner vertices and 1 center. These are
stored in the tessellation as:

- `tessels.unit_xyz` -- Cartesian coordinates `(npix, 5, 3)` on the unit sphere
- `tessels.unit_spherical` -- spherical coordinates `(r, theta, phi)` in `(npix, 5, 3)`

where `theta` is the colatitude (0 at north pole, pi at south pole) and `phi` is
the longitude.

### Neighbors and total variation

HEALPix pixels have 8 neighbors (6 or 7 at special locations). The neighbor
structure is used for total variation regularization:

```julia
tv_info = tv_neighbours_healpix(n)
```

This returns a tuple containing the sparse difference matrix and Hessian used
by the TV regularization terms. For reconstructions where some pixels are never
visible, use:

```julia
tv_info = tv_neighbours_healpix_visible(n, stars)
```

This excludes invisible pixels from the regularization, preventing artifacts at
the limb.

## Longitude/latitude scheme

The longitude/latitude grid is a regular grid in colatitude and longitude:

```julia
ntheta = 50                                     # latitude rings
nphi = 50                                       # pixels per ring
tessels = tessellation_latlong(ntheta, nphi)     # npix = 2500
```

All rings have the same number of pixels (`nphi`), giving a total of
`ntheta * nphi` pixels. The grid ranges over colatitude `[0, pi)` and longitude
`[0, 2*pi)`.

Neighbor information for total variation:

```julia
tv_info = tv_neighbours_longlat(ntheta, nphi)
```

### Spots

The longitude/latitude scheme provides convenience functions for creating and
moving temperature spots:

```julia
tmap = make_circ_spot(tmap, ntheta, nphi, spot_radius, latitude, longitude;
                       bright_frac=0.8)
tmap = make_spot_move(tmap, ntheta, nphi, period, B_rot, tepochs)
```

## Choosing a scheme

| Feature | HEALPix | Lon/Lat |
|---------|---------|---------|
| Equal area pixels | Yes | No (poles oversampled) |
| Multi-resolution | Yes (factor-of-4 up/downsampling) | No |
| Hierarchical nesting | Yes | No |
| Uniform neighbor count | ~8 | 4 |
| Spot creation | Manual | Built-in `make_circ_spot` |

For image reconstruction from interferometry, **HEALPix is recommended** due to
equal-area pixels and multi-resolution support.
