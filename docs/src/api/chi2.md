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

!!! note
    For producing real-space images (visualization, image-plane fitting),
    ROTIR also provides **rasterization** and **NFFT** methods that avoid
    the dense polyft matrix entirely. See [Rasterization & NFFT](@ref)
    for the API and [Direct imaging methods](@ref) for a usage guide.

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
| `reg_type` | `"sobel2"` | Regularization type |
| `verbose` | `true` | Print diagnostics |

## Regularization

| Function | Description |
|----------|-------------|
| `spheroid_regularization(x, g; regularizers, verbose)` | Dispatch to regularization functions, accumulates into gradient `g` |
| `spheroid_total_variation(x, g, tvinfo)` | L1 total variation |
| `spheroid_l2_fg(x, g, tvinfo)` | Quadratic total variation (TV2) |
| `spheroid_harmon_bias_fg(x, g, B)` | Harmonic bias regularization |
| `spheroid_radflat_fg(x, g, bins)` | RADFLAT — flatten the azimuthally averaged radial profile |
| `spheroid_radialvar_fg(x, g, bins)` | RADIALVAR — remove azimuthal scatter within each annulus |
| `radflat_bins(star; nbins)` | Projected-radius binning both of the above need |

Regularizers are passed as `["name", weight, aux, pixel_subset]`. Recognised names:
`mem`, `tv`, `tv2`, `mean`, `bias`, `radflat`, `radialvar` — an unknown name raises rather
than silently contributing zero.

The two radial regularizers have their own page, including when they help and when they
destroy real structure: [Radial Regularizers](radial_regularizers.md).

## Binary forward model

| Function | Description |
|----------|-------------|
| `orbit_to_rotir_offset(bparams, tepoch_jd)` | Convert orbital position to ROTIR's (West, North) projected frame; returns `(offset_x, offset_y)` in mas |
| `binary_phase_shift(uv, offset_x, offset_y)` | Per-baseline phase shift from binary separation |
| `binary_cvis(x1, star1, x2, star2, phase_shift)` | Combined complex visibilities for both stars, flux-normalized |
| `binary_observables(x1, star1, x2, star2, data, phase_shift)` | Returns `(v2, t3amp, t3phi)` for a binary model |
| `binary_chi2_f(x1, star1, x2, star2, data, phase_shift; verbose)` | Binary chi-squared (value only) |

## Binary imaging

Both surfaces of a binary, reconstructed at a separation the orbit (or a fitted offset)
already fixes. The unknown is the two maps **concatenated**, `[x1; x2]`, split at
`stars1[1].npix`; `split_binary_map` takes them apart again.

| Function | Description |
|----------|-------------|
| `binary_chi2_fg(x1, g1, star1, x2, g2, star2, data, phase_shift)` | One epoch's chi-squared and the gradient w.r.t. **both** maps |
| `binary_crit_allepochs_fg(x, g, stars1, stars2, data, phase_shifts; regularizers1, regularizers2, epochs_weights)` | The criterion a binary reconstruction minimises |
| `binary_reconstruct_oi(x_start, data, stars1, stars2, phase_shifts; maxiter, regularizers1, regularizers2, callback)` | VMLMB over both maps |
| `split_binary_map(x, stars1)` | Views of the two halves |

```julia
stars1 = create_star_multiepochs(tess, p1, tepochs; secondary = false)
stars2 = create_star_multiepochs(tess, p2, tepochs; secondary = true)
setup_oi!(data, stars1); setup_oi!(data, stars2)      # both: the dense route is what runs
# The separation is fixed through the run — it comes from the orbit, it is not fitted here.
shifts = [binary_phase_shift(data[i].uv, offs[i]...) for i in eachindex(data)]
x0 = vcat(parametric_temperature_map(p1, stars1[1]),
          parametric_temperature_map(p2, stars2[1]; secondary = true))
x  = binary_reconstruct_oi(x0, data, stars1, stars2, shifts;
                           regularizers1 = regs1, regularizers2 = regs2, maxiter = 200)
map1, map2 = split_binary_map(x, stars1)
```

Two regularizer lists, and neither defaults to the other: an entry carries a structure built
from *its* star (`radflat_bins` bins that component's radii, `orthold_direction` is that
component's degenerate direction), so the primary's list applied to the secondary regularizes
it against the wrong shape. A subset in element 4 indexes its own component's map, not the
concatenated vector.

`intensity_model = :planck` is not available on this path: it makes the map a temperature and
the brightness a nonlinear function of it, while this derivative is the linear one — the same
restriction `spheroid_chi2_fg` carries.

The GUI's Imaging tab does not drive this yet; it reconstructs a single component and refuses
a model with a companion rather than silently imaging the primary alone.

## Parametric fitting

| Function | Description |
|----------|-------------|
| `parametric_temperature_map(params, star)` | Generate von Zeipel temperature map for any surface type |
| `spheroid_parametric_f(params, tessels, data, tepochs)` | End-to-end: parameters to chi2 |

## Saving a surface map

A ROTIR map is a bare `Vector` of per-tessel values: on its own it does not record which
tessellation it belongs to, whether it is a temperature or an intensity, or what geometry it
was fitted against — so it cannot be turned back into a chi-squared a week later. These write
it to FITS together with everything needed to reproduce that, and read it back with the
parameter types it was saved with (`surface_type` and `ldtype` stay `Int`, which matters
because the code branches on them with `==`).

See [Saving a reconstruction](../guides/reconstruction.md) for the round trip.

```@docs
save_surface_map
load_surface_map
```
