# Chi-squared & imaging

## Two forward-model paths

ROTIR computes complex visibilities from a surface map in two ways. They agree to rounding.

### Matrix path

`setup_oi!` stores in each `stellar_geometry` a dense complex matrix `polyft`
(Nuv × visible tessels) and the projected tessel areas `polyflux`. The visibilities are one
matrix-vector product, `cvis = polyft * xw / dot(polyflux, xw)`, with `xw` the visible map
multiplied by `vis_weights` and `ldmap`; the gradient is the transposed product.

`spheroid_chi2_fg`, `image_reconstruct_oi` and `binary_reconstruct_oi` use this path and
require `setup_oi!`. `observables`, `chi2s` and `chi2_breakdown` use it when `polyft` is set
and the fused path otherwise.

### Fused path

`fused_spheroid_chi2_fg` computes the visibilities without forming `polyft`, with the kernel
selected by `POLYFT_BACKEND` (see [Fused polygon Fourier transform](fused_polyft.md)); a
second, adjoint pass gives the gradient. `compute_adjoint_vertices!` also gives the gradient
with respect to the vertex positions, which the joint shape and map fit (`shape_chi2_fg!`,
`joint_reconstruct_oi`) uses.

|                     | Matrix path          | Fused path                  |
|---------------------|----------------------|-----------------------------|
| Memory              | O(Nuv × Npix)        | O(Nuv + Npix)               |
| Setup               | `setup_oi!` per geometry | none                    |
| Per evaluation      | matrix-vector product | transform kernel           |
| Vertex gradients    | no                   | yes                         |
| Used by             | `image_reconstruct_oi`, `spheroid_chi2_fg` | `shape_chi2_fg!`, `joint_reconstruct_oi` |
| Suited to           | fixed geometry, map only | geometry that changes, or a matrix too large for memory |

!!! note
    For real-space images, see [Rasterization & NFFT](@ref) and
    [Direct imaging methods](@ref).

## Setup (matrix path)

| Function | Description |
|----------|-------------|
| `setup_oi!(data, stars)` | `polyflux` and `polyft` for every epoch, stored in `stars` (threaded over epochs) |
| `setup_polygon_ft(data, stars)` | the same arrays, returned as `(polyflux, polyft)` vectors |
| `setup_polyflux_single(proj_west, proj_north)` | projected tessel areas (shoelace formula), one epoch |
| `setup_polyft_single(uv, proj_west, proj_north)` | `polyft` (Nuv × tessels), one epoch |

## Forward model

| Function | Description |
|----------|-------------|
| `poly_to_cvis(x, star; intensity_model, band)` | map → flux-normalised complex visibilities |
| `poly_to_flux(x, star; intensity_model, band)` | map → total flux |
| `observables(x, star, data)` | `(v2, t3amp, t3phi)` |
| `cvis_to_v2(cvis, indx)` | complex visibilities → V² |
| `cvis_to_t3(cvis, i1, i2, i3)` | complex visibilities → `(t3, t3amp, t3phi)`, phase in degrees |
| `mod360(x)` | wrap an angle (degrees) to [−180, 180) |

`intensity_model = :linear` (default) uses the map as surface brightness; `:planck` treats it
as a temperature and converts it with `intensity` at wavelength `band` (metres).

## Chi-squared

| Function | Description |
|----------|-------------|
| `chi2s(x, star, data; verbose)` | `(chi2_v2, chi2_t3amp, chi2_t3phi)`; prints the reduced values |
| `chi2_breakdown(x, stars, data)` | absolute, counts and reduced χ², total and per epoch, as a NamedTuple |
| `spheroid_chi2_f(x, star, data; weights)` | one epoch, value |
| `spheroid_chi2_fg(x, g, star, data; weights)` | one epoch, value and gradient (matrix path) |
| `spheroid_chi2_allepochs_f(x, stars, data; epochs_weights, weights)` | all epochs, value |
| `ROTIR.spheroid_crit_allepochs_fg(x, g, stars, data; regularizers, epochs_weights, weights)` | all epochs, χ² + regularisation, value and gradient |

`weights = (w_V², w_T3amp, w_T3φ)`, default `(1, 1, 1)`, scales each observable's χ² and its
gradient. `epochs_weights` scales each epoch's χ²; it must have one finite, non-negative entry
per epoch.

## Reconstruction

| Function | Description |
|----------|-------------|
| `image_reconstruct_oi(x, data, stars; kwargs...)` | VMLMB with bounds and regularisation |
| `image_reconstruct_oi_crit(x, data, stars; regularizers)` | criterion at `x` |
| `image_reconstruct_oi_chi2(x, data, stars)` | χ² at `x` |
| `image_reconstruct_oi_chi2_fg(x, data, stars)` | `(χ², gradient)` at `x` |
| `multires_reconstruct_oi(data, star_params, tepochs; kwargs...)` | coarse-to-fine HEALPix reconstruction; returns `(tmap, stars)` |

### `image_reconstruct_oi` keywords

| Keyword | Default | Description |
|---------|---------|-------------|
| `maxiter` | `100` | maximum VMLMB iterations |
| `lower` | `0` | lower bound on pixel values |
| `upper` | `Inf` | upper bound on pixel values |
| `regularizers` | `[]` | regularisation terms |
| `epochs_weights` | `[]` | per-epoch weights (empty: uniform) |
| `weights` | `(1, 1, 1)` | `(w_V², w_T3amp, w_T3φ)` |
| `intensity_model`, `band` | `:linear`, `nothing` | see [Forward model](@ref) |
| `callback` | `nothing` | `callback(x, nevals, f)`, called on the first evaluation and every `callback_every` |
| `callback_every` | `25` | evaluations between callbacks |
| `verbose` | `true` | print per-iteration diagnostics |

`callback` receives VMLMB's working vector (copy it to keep it) and counts evaluations, not
iterations. A callback that throws is dropped with a warning; the reconstruction continues.

### `multires_reconstruct_oi` keywords

| Keyword | Default | Description |
|---------|---------|-------------|
| `n_start` | `2` | first HEALPix level |
| `n_end` | `4` | last HEALPix level |
| `maxiter` | `500` | iterations per level |
| `reg_weight` | `10.0` | regularisation weight, the same at every level |
| `reg_type` | `"sobel2"` | `"sobel"` or `"sobel2"` (`sobel_gradient_healpix`); any other name uses `tv_neighbors_healpix` |
| `verbose` | `true` | print diagnostics |

Other keywords go to `image_reconstruct_oi`.

## Regularisation

Each regulariser is an entry `[name, weight, aux, pixel_subset]`:

| Name | Penalty | `aux` |
|------|---------|-------|
| `"sobel"` | `∫|∇x|dΩ / x̄` (edge-preserving) | `sobel_gradient_healpix(n)` |
| `"sobel2"` | `∫|∇x|²dΩ / x̄²` | `sobel_gradient_healpix(n)` |
| `"laplace_beltrami"` | `∫|∇x|²dΩ / x̄²` (cotangent Laplacian) | `laplace_beltrami_healpix(n)` |
| `"tv"` | `‖Lx‖`, `L` the graph Laplacian | `tv_neighbors_healpix(n)` |
| `"tv2"` | `‖Lx‖²` | `tv_neighbors_healpix(n)` |
| `"mem"` | `Σ xm log xm`, `xm = x/x̄` | `nothing` |
| `"mean"` | `Σ |x − x̄|` | `nothing` |
| `"bias"` | weighted variance about the mean, weight `B` above the mean | `B` |
| `"radflat"` | flat radial profile | `radflat_bins(star)` |
| `"radialvar"` | no azimuthal scatter | `radflat_bins(star)` |
| `"orthold"` | no component along the limb-darkening-degenerate direction | `orthold_direction(star, x0, params)` |

`pixel_subset` (e.g. `1:npix`) selects the pixels the term acts on; `aux` must be built for
the same subset. An unknown name raises an error.

| Function | Description |
|----------|-------------|
| `spheroid_regularization(x, g; regularizers, verbose)` | sum of the terms; adds their gradient into `g` |
| `ROTIR.spheroid_total_variation_fg(x, g, tvinfo)` | `"tv"` |
| `ROTIR.spheroid_total_variation2_fg(x, g, tvinfo)` | `"tv2"` |
| `max_entropy_fg(x, g)` | `"mem"` |
| `spheroid_harmon_bias_fg(x, g, B)` | `"bias"` |
| `spheroid_radflat_fg(x, g, bins)` | `"radflat"` |
| `spheroid_radialvar_fg(x, g, bins)` | `"radialvar"` |
| `spheroid_orthold_fg(x, g, od)` | `"orthold"` |

```@docs
sobel_gradient_healpix
laplace_beltrami_healpix
spheroid_sobel_fg
spheroid_sobel2_fg
spheroid_laplace_beltrami_fg
```

`"radflat"`, `"radialvar"` and `"orthold"` are described in
[Radial Regularizers](radial_regularizers.md).

## Binary forward model

| Function | Description |
|----------|-------------|
| `orbit_to_rotir_offset(bparams, tepoch_jd)` | orbital position as ROTIR's `(West, North)` offset, mas |
| `binary_phase_shift(uv, offset_x, offset_y)` | per-baseline phase factor for that offset |
| `binary_cvis(x1, star1, x2, star2, phase_shift)` | combined, flux-normalised complex visibilities |
| `binary_observables(x1, star1, x2, star2, data, phase_shift)` | `(v2, t3amp, t3phi)` |
| `binary_chi2_f(x1, star1, x2, star2, data, phase_shift)` | χ², value |

## Binary imaging

Both surfaces of a binary, reconstructed at a fixed separation (from the orbit or a fitted
offset). The unknown is the two maps concatenated, `[x1; x2]`, split at `stars1[1].npix`.

| Function | Description |
|----------|-------------|
| `binary_chi2_fg(x1, g1, star1, x2, g2, star2, data, phase_shift)` | one epoch's χ² and the gradient with respect to both maps |
| `binary_crit_allepochs_fg(x, g, stars1, stars2, data, phase_shifts; regularizers1, regularizers2, epochs_weights, weights)` | the criterion `binary_reconstruct_oi` minimises |
| `binary_reconstruct_oi(x_start, data, stars1, stars2, phase_shifts; maxiter, regularizers1, regularizers2, weights, callback)` | VMLMB over both maps |
| `split_binary_map(x, stars1)` | views of the two halves |

```julia
stars1 = create_star_multiepochs(tess, p1, tepochs; secondary = false)
stars2 = create_star_multiepochs(tess, p2, tepochs; secondary = true)
setup_oi!(data, stars1); setup_oi!(data, stars2)      # matrix path
shifts = [binary_phase_shift(data[i].uv, offs[i]...) for i in eachindex(data)]  # fixed
x0 = vcat(parametric_temperature_map(p1, stars1[1]),
          parametric_temperature_map(p2, stars2[1]; secondary = true))
x  = binary_reconstruct_oi(x0, data, stars1, stars2, shifts;
                           regularizers1 = regs1, regularizers2 = regs2, maxiter = 200)
map1, map2 = split_binary_map(x, stars1)
```

Each component has its own regulariser list. Structures in an entry (`radflat_bins`,
`orthold_direction`) are built from that component's star, and a pixel subset in element 4
indexes that component's map, not the concatenated vector.

`binary_cvis`, `binary_observables` and `binary_chi2_f` take `intensity_model`, `band` and
`occultation` (default `false`). `binary_chi2_fg` takes `occultation` only;
`binary_crit_allepochs_fg` and `binary_reconstruct_oi` take neither, so a binary
reconstruction uses the map as surface brightness and has no mutual occultation.

The GUI's Imaging tab reconstructs one component and refuses a model with a companion.

## Parametric model

| Function | Description |
|----------|-------------|
| `parametric_temperature_map(params, star; secondary)` | gravity-darkened temperature map for any `surface_type` |
| `spheroid_parametric_f(params, tessels, data, tepochs)` | parameters → χ² over all epochs |

## Saving a surface map

`save_surface_map` writes a map to FITS with its tessellation, the `star_params` it was
computed with (integer fields stay `Int`), and optionally the epochs, so its χ² can be
reproduced. See [Saving a reconstruction](../guides/reconstruction.md).

```@docs
save_surface_map
load_surface_map
```

## Saving a model geometry

`save_star_geometry` writes the surface instead of the map: the tessellated mesh, the
parameters it was built from, and for a binary the second component and its placement. The
mesh is stored, so it is reproduced exactly, including a Roche component whose shape comes
from a root solve at the instantaneous separation.

```julia
g = load_star_geometry("spica_primary.fits")
g.star                       # the mesh as saved
g.params                     # the parameters that built it
setup_oi!(data, [g.star])    # the Fourier setup is not stored
```

```@docs
save_star_geometry
load_star_geometry
```
