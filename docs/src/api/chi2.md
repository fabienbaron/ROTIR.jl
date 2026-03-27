# Chi-squared & imaging

## Setup

| Function | Description |
|----------|-------------|
| `setup_oi!(data, stars)` | Precompute polyflux and polyft matrices for all epochs (threaded) |
| `setup_polygon_ft(data, star)` | Return `(polyflux, polyft)` arrays for all epochs |
| `setup_polyflux_single(proj_west, proj_north)` | Shoelace polygon areas for one epoch |
| `setup_polyft_single(uv, proj_west, proj_north)` | Complex visibility matrix (nuv x npix) for one epoch |

## Forward model

| Function | Description |
|----------|-------------|
| `poly_to_cvis(x, star)` | Temperature map to flux-normalized complex visibilities |
| `poly_to_flux(x, star)` | Temperature map to total flux |
| `observables(x, star, data)` | Returns `(v2_model, t3amp_model, t3phi_model)` |
| `cvis_to_v2(cvis, indx)` | Complex visibilities to squared visibilities |
| `cvis_to_t3(cvis, i1, i2, i3)` | Complex visibilities to triple product, T3amp, T3phi |
| `mod360(x)` | Wrap angle to [-180, 180] |

## Chi-squared

| Function | Description |
|----------|-------------|
| `chi2s(x, star, data; verbose)` | Per-observable chi2: returns `(chi2_v2, chi2_t3amp, chi2_t3phi)` |
| `spheroid_chi2_f(x, star, data)` | Single-epoch chi2 (value only) |
| `spheroid_chi2_fg(x, g, star, data)` | Single-epoch chi2 + gradient (matrix-based) |
| `spheroid_chi2_allepochs_f(x, stars, data)` | Multi-epoch chi2 (value only) |
| `spheroid_crit_multiepochs_fg(x, g, stars, data; regularizers)` | Multi-epoch chi2 + regularization + gradient |

## Reconstruction

| Function | Description |
|----------|-------------|
| `image_reconstruct_oi(x, data, stars; kwargs...)` | Main reconstruction: VMLMB optimizer with bounds and regularization |
| `image_reconstruct_oi_crit(x, data, stars; regularizers)` | Evaluate criterion at fixed point |
| `image_reconstruct_oi_chi2(x, data, stars)` | Evaluate chi2 at fixed point |
| `image_reconstruct_oi_chi2_fg(x, data, stars)` | Evaluate chi2 + gradient at fixed point |
| `multires_reconstruct_oi(data, star_params, tepochs; n_start, n_end, kwargs...)` | Multi-resolution HEALPix pyramid reconstruction |

### `image_reconstruct_oi` keywords

| Keyword | Default | Description |
|---------|---------|-------------|
| `maxiter` | `200` | Maximum VMLMB iterations |
| `lower` | `0` | Lower bound on pixel values |
| `upper` | `Inf` | Upper bound on pixel values |
| `regularizers` | `[]` | List of regularization terms |
| `epochs_weights` | `[]` | Per-epoch weights (empty = uniform) |
| `verbose` | `true` | Print per-iteration diagnostics |

### `multires_reconstruct_oi` keywords

| Keyword | Default | Description |
|---------|---------|-------------|
| `n_start` | `2` | Starting HEALPix level |
| `n_end` | `4` | Final HEALPix level |
| `maxiter` | `500` | Max iterations per level |
| `reg_weight` | `1e-5` | TV regularization weight |
| `reg_type` | `"tv2"` | Regularization type |
| `verbose` | `true` | Print diagnostics |

## Regularization

| Function | Description |
|----------|-------------|
| `spheroid_regularization(x, g; regularizers, verbose)` | Dispatch to regularization functions, accumulates into gradient `g` |
| `spheroid_total_variation(x, g, tvinfo)` | L1 total variation |
| `spheroid_l2_fg(x, g, tvinfo)` | Quadratic total variation (TV2) |
| `spheroid_harmon_bias_fg(x, g, B)` | Harmonic bias regularization |

## Parametric fitting

| Function | Description |
|----------|-------------|
| `parametric_temperature_map(params, star)` | Generate von Zeipel temperature map for any surface type |
| `spheroid_parametric_f(params, tessels, data, tepochs)` | End-to-end: parameters to chi2 |
