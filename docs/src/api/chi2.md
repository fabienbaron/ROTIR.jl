# Chi-squared & imaging

## Two forward-model paths

ROTIR provides two ways to compute complex visibilities from a
temperature map.  Both produce identical results; the choice depends on
what is being optimized.

### Matrix path (precomputed polyft)

`setup_oi!()` precomputes a dense complex matrix `polyft` (Nuv x Npix)
and a flux vector `polyflux` (Npix) for each epoch, stored inside the
`stellar_geometry` struct.  Complex visibilities are then a single
matrix-vector multiply: `cvis = polyft * xw / flux`.  The gradient is
the transpose multiply.

This is the default path used by `spheroid_chi2_fg`,
`image_reconstruct_oi`, and the standard reconstruction pipeline.

### Fused matrix-free path

`fused_spheroid_chi2_fg` (in `fused_polyft.jl`) computes visibilities
on-the-fly in a loop over pixels and UV points, without ever forming the
dense polyft matrix.  A second adjoint pass computes the gradient.  This
path also supports `compute_adjoint_vertices!`, which backpropagates
through the vertex positions — needed for joint shape + map optimization
(`shape_chi2_fg!`).

### When to use which

|                     | Matrix path          | Fused path                  |
|---------------------|----------------------|-----------------------------|
| Memory              | O(Nuv x Npix) dense  | O(Nuv + Npix)               |
| Setup cost          | One-time `setup_oi!` | None                        |
| Per-iteration cost  | Mat-vec multiply     | Loop over pixels x UV pts   |
| Gradient            | Transpose multiply   | Adjoint loop                |
| Vertex gradients    | No                   | Yes                         |
| Used by             | `image_reconstruct_oi`, `spheroid_chi2_fg` | `shape_chi2_fg!`, `joint_reconstruct_oi` |
| Best when           | Map-only optimization (fixed geometry) | Shape optimization, or large Nuv x Npix |

For most reconstructions where only the temperature map is optimized,
the matrix path is simpler and fast (a single BLAS call per epoch).
Switch to the fused path when optimizing shape parameters (inclination,
radii, position angle) jointly with the map, or when the polyft matrix
is too large to fit in memory.

## Setup (matrix path)

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

## Binary forward model

| Function | Description |
|----------|-------------|
| `orbit_to_rotir_offset(bparams, tepoch_jd)` | Convert orbital position to ROTIR's (West, North) projected frame; returns `(offset_x, offset_y)` in mas |
| `binary_phase_shift(uv, offset_x, offset_y)` | Per-baseline phase shift from binary separation |
| `binary_cvis(x1, star1, x2, star2, phase_shift)` | Combined complex visibilities for both stars, flux-normalized |
| `binary_observables(x1, star1, x2, star2, data, phase_shift)` | Returns `(v2, t3amp, t3phi)` for a binary model |
| `binary_chi2_f(x1, star1, x2, star2, data, phase_shift; verbose)` | Binary chi-squared (value only) |

## Parametric fitting

| Function | Description |
|----------|-------------|
| `parametric_temperature_map(params, star)` | Generate von Zeipel temperature map for any surface type |
| `spheroid_parametric_f(params, tessels, data, tepochs)` | End-to-end: parameters to chi2 |
