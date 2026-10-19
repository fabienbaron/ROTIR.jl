# Image reconstruction

ROTIR reconstructs a surface temperature map by minimising χ² plus regularisation terms, with
the VMLMB bounded quasi-Newton optimiser from OptimPackNextGen.jl.

## Basic reconstruction

```julia
using ROTIR

oifitsfiles = ["epoch1.oifits", "epoch2.oifits"]
data_all = readoifits_multiepochs(oifitsfiles)   # data_all[wavelength_bin, epoch]
data = data_all[1, :]
tepochs = [d.mean_mjd for d in data] .- data[1].mean_mjd   # days

n = 3                                            # HEALPix level: nside = 8, 768 pixels
tessels = tessellation_healpix(n)

star_params = (
    surface_type=2,           # 0=sphere, 1=ellipsoid, 2=rapid rotator, 3=Roche lobe
    rpole=1.37,               # polar radius in milliarcseconds
    tpole=4800.0,             # polar temperature in Kelvin
    ldtype=3, ld1=0.23, ld2=0.0,  # Hestroffer power-law limb darkening
    inclination=78.0,         # inclination in degrees (90°=edge-on)
    position_angle=24.0,      # position angle of rotation axis on sky (degrees)
    rotation_period=54.8,     # rotation period in days
    beta=0.08,                # von Zeipel exponent: T ∝ g^β (e.g. 0.25 radiative, 0.08 convective)
    frac_escapevel=0.9,       # Ω / Ω_crit
    B_rot=0.0,                # differential rotation (0 = solid body)
)

stars = create_star_multiepochs(tessels, star_params, tepochs)   # geometry at each epoch
tmap_start = parametric_temperature_map(star_params, stars[1])   # gravity-darkened start

# Polygon fluxes and Fourier matrices, stored in each star (polyflux, polyft)
setup_oi!(data, stars)

# ["type", weight, operator, pixel_range]
regularizers = [["sobel2", 10.0, sobel_gradient_healpix(n), 1:length(tmap_start)]]

tmap = image_reconstruct_oi(tmap_start, data, stars;
                             maxiter=500, regularizers=regularizers, verbose=true)
```

## Regularizations

`regularizers` is a list; each entry is `[name, weight, operator, pixel_range]`. The
regulariser applies to `x[pixel_range]`.

| Name | Syntax | Penalises | Response |
|------|--------|-----------|----------|
| `"sobel"` | `["sobel", mu, sobel_info, pixel_range]` | `∫\|∇x\|dΩ` (isotropic L1, edge-preserving) | `k²` |
| `"sobel2"` | `["sobel2", mu, sobel_info, pixel_range]` | `∫\|∇x\|²dΩ` | `k²` |
| `"laplace_beltrami"` | `["laplace_beltrami", mu, lb_info, pixel_range]` | `∫\|∇x\|²dΩ` (converges to the continuum) | `k²` |
| `"tv2"` | `["tv2", mu, tv_info, pixel_range]` | `‖Lx‖²` (curvature) | `k⁴` |
| `"tv"` | `["tv", mu, tv_info, pixel_range]` | `‖Lx‖` (curvature) | `k⁴` |
| `"mem"` | `["mem", mu, nothing, pixel_range]` | maximum entropy of `x/mean(x)` (no spatial coupling) | — |
| `"mean"` | `["mean", mu, nothing, pixel_range]` | `Σ\|x − mean(x)\|` | — |
| `"bias"` | `["bias", mu, B, pixel_range]` | `(x − mean(x))²`, weighted by `B` above the mean | — |

Radial regularizers for the single-epoch, non-rotating case (`"radflat"`, `"radialvar"`,
`"orthold"`) are described in [Radial regularizers](../api/radial_regularizers.md).

### Which to use

`"sobel"` and `"sobel2"` use `sobel_gradient_healpix`, a tangent-plane least-squares gradient
(the spherical analogue of the 2-D Sobel stencil). They penalise the first derivative: each
neighbour difference is squared before summing, so `‖Dx‖²` on a tessel approximates `|∇x|²`
there and only a constant map has zero cost.

`"laplace_beltrami"` uses `laplace_beltrami_healpix`, the cotangent (P1 finite-element)
Laplace–Beltrami operator on the triangulation of the tessel centres. Its `‖Dx‖²` is the
Dirichlet energy `∫|∇x|²dΩ` of the piecewise-linear interpolant. It converges at O(h²) and
isotropically (on `∫|∇z|²`: 1.6 %, 0.4 %, 0.1 % low at n = 2, 3, 4), so a weight has the same
meaning at every resolution. `"sobel2"` is within a few per cent of the same limit. At one
resolution the two give the same map (ζ And, 768 tessels, VarInf MAP: r = 0.999, each at its
own roughness-matched weight). Use `"laplace_beltrami"` when a weight is carried between
resolutions (coarse-to-fine, multires).

`"tv"` and `"tv2"` use `tv_neighbors_healpix`, the graph Laplacian (diagonal = neighbour count,
`−1` off-diagonal), and penalise the second derivative. On HEALPix, `"tv"` is `‖Lx‖`, a single
global norm and a monotone transform of `"tv2"`; it is not an edge-preserving total variation
(that is `"sobel"`).

The `k⁴` response of `"tv"`/`"tv2"` leaves large scales almost untouched and suppresses small
ones strongly, so between a spiky map and a nearly constant one there is little usable range.
The `k²` response suppresses tessel-scale noise and keeps larger structure.

When transferring a weight between targets or resolutions:

* `"sobel"`, `"sobel2"`, `"laplace_beltrami"` and `"mem"` are scale-invariant (normalised by
  `mean(x)`), as is χ² (the visibilities are flux-normalised). `"tv"` and `"tv2"` are not: the
  same weight is about 42× stronger on a 25000 K star than on a 3900 K one.
* `"sobel"` and `"sobel2"` include the `4π/npix` solid angle, so on a smooth map their value
  converges as the mesh refines. `"tv2"` falls by about 0.4× per `nside` doubling.

Build the operator argument with `sobel_gradient_healpix(n)` for the Sobel forms,
`laplace_beltrami_healpix(n)` for `"laplace_beltrami"`, or `tv_neighbors_healpix(n)` /
`tv_neighbors_longlat(ntheta, nphi)` for the Laplacian ones. `"laplace_beltrami"` is
normalised by `mean(x)²` like `"sobel2"`, and its energy includes the solid angle.

### Measured comparison

Each prior at its own best weight, from the ladder in `demos/rwcep_radflat.jl`
(`ALPHAS=0 TVWS=0,1e1,1e2,1e3,1e4`, RW Cep at nside 3, 2819 data against 428 patches, no
radial regularizer):

| prior | best weight | χ²ᵥ₂/n | χ²ₜ₃ₚ/n |
|---|---|---|---|
| `"tv2"` | 1e−4 | 2.29 | 1.7 |
| *(none)* | — | 2.00 | 1.5 |
| `"sobel2"` | 1e1 | 1.93 | 1.5 |
| `"sobel"` | 1e1 | 1.86 | 1.4 |

`"tv2"` fits worse than no prior: strong enough to suppress tessel-scale noise, it also
flattens structure the closure phases constrain. `"sobel"` fits best on both observables.

The weights are orders of magnitude. Sobel weights are `O(1–10²)` and independent of `tpole`
and `nside`, but the balance against χ² scales with the number of data points: re-run the
ladder for a new dataset.

!!! warning "Weights differ between the two families"
    On the same data a `"tv2"` weight is `O(10⁻⁵–10⁻⁴)` and a `"sobel2"` weight `O(1–10²)`.
    Swapping one name for the other without re-tuning gives either no regularisation or a
    prior that dominates the likelihood.

## Evaluating the fit

```julia
# χ² + regularisation
crit = image_reconstruct_oi_crit(tmap, data, stars; regularizers=regularizers, verbose=true)

# χ² only
chi2 = image_reconstruct_oi_chi2(tmap, data, stars; verbose=true)

# χ² per observable, one epoch (prints the reduced values)
chi2_v2, chi2_t3amp, chi2_t3phi = chi2s(tmap, stars[1], data[1]; verbose=true)

# The same as a NamedTuple, with counts and reduced values
b = chi2_breakdown(tmap, stars[1], data[1])

# Model observables
v2_model, t3amp_model, t3phi_model = observables(tmap, stars[1], data[1])
```

## Keyword reference

`image_reconstruct_oi` keywords:

| Keyword | Default | Description |
|---------|---------|-------------|
| `maxiter` | `100` | maximum optimiser iterations |
| `lower` | `0` | lower bound on pixel values |
| `upper` | `Inf` | upper bound on pixel values |
| `regularizers` | `[]` | regularisation terms |
| `epochs_weights` | `[]` | per-epoch weights (empty = equal) |
| `weights` | `(1, 1, 1)` | weights on the V², T3amp and T3φ χ² |
| `intensity_model` | `:linear` | temperature-to-intensity law (see `poly_to_cvis`) |
| `band` | `nothing` | photometric band for the intensity law |
| `callback` | `nothing` | `callback(x, neval, f)`, called every `callback_every` (25) evaluations |
| `verbose` | `true` | print per-iteration diagnostics |

## Chi-squared functions

| Function | Method |
|----------|--------|
| `spheroid_chi2_fg(x, g, star, data)` | uses the `polyft` matrix built by `setup_oi!` (fast for small problems) |
| `fused_spheroid_chi2_fg(x, g, star, data)` | matrix-free (less memory for large problems) |

Both take the same arguments and return the same χ² and gradient to rounding.

## Soft visibility

Shape fits weight pixels near the limb with a sigmoid instead of the hard mask
`normals_z > 0`, which makes χ² differentiable with respect to the shape parameters
(`soft_visibility`):

```
w(p) = sigmoid(kappa * nz(p))
```

where `nz` is the z-component of the surface normal and `kappa` (default 50) sets the
sharpness.

## Joint shape + map optimization

`joint_reconstruct_oi` alternates between optimising the map at fixed shape and the shape
(radii, inclination, position angle) at fixed map, with analytic gradients.

!!! note
    `tessellation_healpix` defaults to `Float32`. The map and `θ` must have the same element
    type; for Float64 use `tessellation_healpix(n, T=Float64)`.

```julia
# [rpole (mas), frac_escapevel, inclination (°), PA (°)]
θ_start = [1.37, 0.9, 78.0, 24.0]

tmap_final, θ_final = joint_reconstruct_oi(
    tmap_start, θ_start, data, tessels, star_params, tepochs;
    maxiter_xmap=200,     # iterations per map step
    maxiter_θ=50,         # iterations per shape step
    nouter=5,             # alternating cycles
    regularizer="sobel",  # map regulariser (default)
    reg_weight=10.0,      # its weight (default)
    κ=50.0,               # soft-visibility sharpness
    θ_lower=[0.5, 0.0, 0.0, -180.0],    # [rpole, ω, inc, PA] lower bounds
    θ_upper=[3.0, 1.0, 180.0, 180.0],   # [rpole, ω, inc, PA] upper bounds
)
```

`θ` depends on `surface_type`:

| Surface type | Parameters |
|-------------|------------|
| 0 (Sphere) | `[radius, inclination, PA]` |
| 1 (Ellipsoid) | `[rx, ry, rz, inclination, PA]` |
| 2 (Rapid Rotator) | `[rpole, frac_escapevel, inclination, PA]` |

Surface type 3 (Roche) has no analytic shape gradient and is refused.

## Single-epoch imaging

Single-epoch data (slow rotators such as Cepheids and supergiants, one-night snapshots, or
diameter and limb-darkening fits) use the same workflow with `tepochs = [0.0f0]`:

```julia
using ROTIR

data_all = readoifits_multiepochs(["polaris.oifits"], T=Float32)
data = data_all[1, :]          # length 1
tepochs = [0.0f0]

n = 4                          # nside = 16, 3072 pixels
tessels = tessellation_healpix(n)

star_params = (
    surface_type    = 0,       # sphere
    radius          = 1.6,     # angular radius in milliarcseconds
    tpole           = 6000.0,  # effective temperature in Kelvin
    ldtype          = 3,       # Hestroffer power-law limb darkening
    ld1             = 0.24,    # limb-darkening coefficient
    ld2             = 0.0,     # unused for Hestroffer
    inclination     = 90.0,    # equator-on; arbitrary for a sphere
    position_angle  = 0.0,     # arbitrary for a sphere
    rotation_period = 1.0,     # irrelevant at t = 0
)

stars = create_star_multiepochs(tessels, star_params, tepochs)
setup_oi!(data, stars)
tmap_start = parametric_temperature_map(star_params, stars[1])   # uniform

# A single snapshot has sparser uv coverage: higher weight than multi-epoch
regularizers = [["sobel2", 100.0, sobel_gradient_healpix(n), 1:length(tmap_start)]]

tmap = image_reconstruct_oi(tmap_start, data, stars;
    maxiter=500, regularizers=regularizers, verbose=true)

chi2 = image_reconstruct_oi_chi2(tmap, data, stars; verbose=true)
plot2d(tmap, stars[1], intensity=true, compass=true)
```

### Key differences from multi-epoch

| Aspect | Single epoch | Multi-epoch |
|--------|-------------|-------------|
| `tepochs` | `[0.0f0]` | `[0.0f0, 3.5f0, ...]` |
| Rotation | not used | sets the phase coverage |
| Inclination/PA | set the projected geometry | same |
| `stars` array | Length 1 | Length N |
| `data` array | Length 1 | Length N |
| uv coverage | one night | grows with epochs |

### Choosing parameters for non-rotating targets

For a single-epoch sphere:

- `rotation_period`: any nonzero value. The phase `2π · t / P` is zero at `t = 0`.
- `position_angle`: rotates only the projected coordinate frame; `0.0` unless known.
- `inclination`: does not change the projected shape; `90.0` (equator-on), or a known value.

For an ellipsoid or rapid rotator, inclination and PA set the projected shape and the limb
darkening even at one epoch.

### Regularization tips for single-epoch data

- Use `"sobel"` or `"sobel2"` with a higher weight than for multi-epoch data (the example
  above uses 100 against 10).
- `lower` / `upper` bound the pixel values, e.g. `lower=3000` for a 6000 K star.

See [`demos/polaris_imaging.jl`](https://github.com/fabienbaron/ROTIR.jl/blob/main/demos/polaris_imaging.jl)
reconstructs Polaris from one epoch.

## Saving a reconstruction

A map is a `Vector` of per-tessel values without its tessellation or geometry.
[`save_surface_map`](@ref) writes it to FITS with the HEALPix level, every `star_params` field
and the epoch times:

```julia
save_surface_map("polaris.fits", x, params;
                 nside_exp = 4, tepochs = tepochs, mjd = mjd,
                 chi2 = chi2, ndata = ndata)
```

[`load_surface_map`](@ref) reads it back with the parameter types it was saved with
(`surface_type` and `ldtype` stay `Int`), so the χ² can be reproduced from the file:

```julia
m     = load_surface_map("polaris.fits")
stars = create_star_multiepochs(tessellation_healpix(m.nside_exp), m.params, m.tepochs)
setup_oi!(data, stars)
chi2_breakdown(m.x, stars[1], data[1])
```

The GUI's **Save model map** and **Load model map** buttons call these. Save writes the
reconstruction if there is one, otherwise the model's parametric map; load restores the
parameters as a model and the values as an image.
