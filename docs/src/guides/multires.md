# Multi-resolution imaging

ROTIR supports a coarse-to-fine reconstruction strategy using the hierarchical
structure of HEALPix. Starting at a low resolution, the reconstruction is
progressively refined by upsampling and re-optimizing at each level.

## Why multi-resolution?

At high HEALPix resolution (n=4 or 5), the optimizer faces a high-dimensional
problem with many local minima. Starting from a coarse solution and refining
helps:

- **Faster convergence**: fewer pixels at low resolution means fewer iterations
- **Better global minimum**: coarse features are established first
- **Robust regularization**: the effective regularization adapts to each scale

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
                                       maxiter=500, reg_weight=1e-5,
                                       reg_type="tv2", verbose=true)
```

This runs reconstruction at HEALPix levels n=2, 3, 4:

| Level | nside | npix | Description |
|-------|-------|------|-------------|
| n=2 | 4 | 192 | Coarse: captures large-scale features |
| n=3 | 8 | 768 | Medium: refines spot structure |
| n=4 | 16 | 3072 | Fine: resolves sharp features |

At each level, the map from the previous level is upsampled (each HEALPix pixel
splits into 4 children) and used as the starting point.

### Keywords

| Keyword | Default | Description |
|---------|---------|-------------|
| `n_start` | `2` | Starting HEALPix level |
| `n_end` | `4` | Final HEALPix level |
| `maxiter` | `500` | Max iterations per level |
| `reg_weight` | `1e-5` | TV regularization weight |
| `reg_type` | `"tv2"` | Regularization type (`"tv"` or `"tv2"`) |
| `verbose` | `true` | Print diagnostics |

Additional keywords are passed through to `image_reconstruct_oi` (e.g., `lower`,
`upper`).

## Manual multi-resolution

For more control, you can manually upsample and reconstruct at each level:

```julia
# Level n=2
n = 2
tessels = tessellation_healpix(n)
stars = create_star_multiepochs(tessels, star_params, tepochs)
tmap = parametric_temperature_map(star_params, stars[1])
setup_oi!(data, stars)
regularizers = [["tv2", 1e-5, tv_neighbors_healpix(n), 1:length(tmap)]]
tmap = image_reconstruct_oi(tmap, data, stars; maxiter=500, regularizers=regularizers)

# Upsample to level n=3
tmap, stars = upsample_map_stars(tmap, stars, star_params, tepochs)
setup_oi!(data, stars)
n = 3
regularizers = [["tv2", 1e-5, tv_neighbors_healpix(n), 1:length(tmap)]]
tmap = image_reconstruct_oi(tmap, data, stars; maxiter=500, regularizers=regularizers)

# Upsample to level n=4
tmap, stars = upsample_map_stars(tmap, stars, star_params, tepochs)
setup_oi!(data, stars)
n = 4
regularizers = [["tv2", 1e-5, tv_neighbors_healpix(n), 1:length(tmap)]]
tmap = image_reconstruct_oi(tmap, data, stars; maxiter=500, regularizers=regularizers)
```

## Downsampling

To reduce resolution (e.g., for quick visualization or comparison):

```julia
tmap_coarse, stars_coarse = downsample_map_stars(tmap, stars, star_params, tepochs)
```

Each group of 4 child pixels is averaged into one parent pixel.
