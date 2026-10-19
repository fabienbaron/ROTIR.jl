# Multi-resolution imaging

Coarse-to-fine reconstruction on HEALPix: reconstruct at a low level, split each pixel into its
four children, and reconstruct again at the next level. Fewer pixels at the coarse levels mean
faster iterations, and the large-scale structure is fitted before the small-scale one.

## Using `multires_reconstruct_oi`

```julia
using ROTIR

# Load data
oifitsfiles = ["epoch1.oifits", "epoch2.oifits"]
data_all = readoifits_multiepochs(oifitsfiles)
data = data_all[1, :]
tepochs = [d.mean_mjd for d in data] .- data[1].mean_mjd

star_params = (
    surface_type=2, rpole=1.37, tpole=4800.0, ldtype=3, ld1=0.23, ld2=0.0,
    inclination=78.0, position_angle=24.0, rotation_period=54.8,
    beta=0.08, frac_escapevel=0.9, B_rot=0.0,
)

tmap, stars = multires_reconstruct_oi(data, star_params, tepochs;
                                       n_start=2, n_end=4,
                                       maxiter=500, reg_weight=10.0,
                                       reg_type="sobel2", verbose=true)
```

This reconstructs at HEALPix levels n=2, 3, 4:

| Level | nside | npix |
|-------|-------|------|
| n=2 | 4 | 192 |
| n=3 | 8 | 768 |
| n=4 | 16 | 3072 |

Each level starts from the previous map, upsampled (each pixel copied into its 4 children).
The first level starts from `parametric_temperature_map`.

### Keywords

| Keyword | Default | Description |
|---------|---------|-------------|
| `n_start` | `2` | Starting HEALPix level |
| `n_end` | `4` | Final HEALPix level |
| `maxiter` | `500` | Max iterations per level |
| `reg_weight` | `10.0` | regularization weight, the same at every level |
| `reg_type` | `"sobel2"` | `"sobel"`, `"sobel2"`, `"tv"` or `"tv2"` |
| `verbose` | `true` | Print diagnostics |

Other keywords go to `image_reconstruct_oi` (e.g. `lower`, `upper`).

## Manual multi-resolution

The same ladder, step by step:

```julia
# Level n=2
n = 2
tessels = tessellation_healpix(n)
stars = create_star_multiepochs(tessels, star_params, tepochs)
tmap = parametric_temperature_map(star_params, stars[1])
setup_oi!(data, stars)
regularizers = [["sobel2", 10.0, sobel_gradient_healpix(n), 1:length(tmap)]]
tmap = image_reconstruct_oi(tmap, data, stars; maxiter=500, regularizers=regularizers)

# Upsample to level n=3
tmap, stars = upsample_map_stars(tmap, stars, star_params, tepochs)
setup_oi!(data, stars)
n = 3
regularizers = [["sobel2", 10.0, sobel_gradient_healpix(n), 1:length(tmap)]]
tmap = image_reconstruct_oi(tmap, data, stars; maxiter=500, regularizers=regularizers)

# Upsample to level n=4
tmap, stars = upsample_map_stars(tmap, stars, star_params, tepochs)
setup_oi!(data, stars)
n = 4
regularizers = [["sobel2", 10.0, sobel_gradient_healpix(n), 1:length(tmap)]]
tmap = image_reconstruct_oi(tmap, data, stars; maxiter=500, regularizers=regularizers)
```

`"sobel"` and `"sobel2"` include the `4π/npix` solid angle, so one weight has the same effect at
every level. A `"tv2"` term on a smooth map falls by about ×0.4 per doubling of `nside`, so a
fixed `"tv2"` weight regularizes the finest level least. Sobel weights are typically O(1–10²);
`"tv2"` weights O(10⁻⁵–10⁻⁴).

## Downsampling

To go down one level:

```julia
tmap_coarse, stars_coarse = downsample_map_stars(tmap, stars, star_params, tepochs)
```

Each group of 4 child pixels is averaged into its parent. For maps alone,
`upsample_healpix_xmap(x, n_coarse, n_fine)` and `downsample_healpix_xmap(x, n_fine, n_coarse)`
change resolution by any number of levels.
