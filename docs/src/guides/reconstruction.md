# Image reconstruction

ROTIR reconstructs stellar surface temperature maps by minimizing the
chi-squared between model and observed interferometric data, subject to
regularization. The optimizer is VMLMB (variable metric limited-memory
quasi-Newton with bounds) from OptimPackNextGen.jl.

## Basic reconstruction

```julia
using ROTIR

# Read OIFITS files — returns a 2D array: data_all[wavelength_bin, epoch]
oifitsfiles = ["epoch1.oifits", "epoch2.oifits"]
data_all = readoifits_multiepochs(oifitsfiles)

# Select the first wavelength bin across all epochs
data = data_all[1, :]

# Compute relative epoch times (days since first observation) for tracking rotation
tepochs = [d.mean_mjd for d in data] .- data[1].mean_mjd

# Tessellate the stellar surface using nested HEALPix with resolution level 3
# (nside=2^3=8, giving 768 equal-area pixels)
n = 3
tessels = tessellation_healpix(n)

# Define the stellar model: a rapid rotator (surface_type=2)
star_params = (
    surface_type=2,           # 0=sphere, 1=ellipsoid, 2=rapid rotator, 3=Roche lobe
    rpole=1.37,               # polar radius in milliarcseconds
    tpole=4800.0,             # polar temperature in Kelvin
    ldtype=3, ld1=0.23, ld2=0.0,  # Hestroffer power-law limb darkening
    inclination=78.0,         # inclination in degrees (90°=edge-on)
    position_angle=24.0,      # position angle of rotation axis on sky (degrees)
    rotation_period=54.8,     # rotation period in days
    beta=0.08,                # von Zeipel exponent: T ∝ g^β (e.g. 0.25 radiative, 0.08 convective)
    frac_escapevel=0.9,       # rotational velocity as fraction of escape velocity
    B_rot=0.0,                # differential rotation coefficient (0=solid body)
)

# Build projected, rotated stellar geometry for each epoch
# (computes surface shape, normals, visible pixels, and limb darkening)
stars = create_star_multiepochs(tessels, star_params, tepochs)

# Generate initial temperature map from the von Zeipel gravity-darkening law
# (serves as the starting point for the optimizer)
tmap_start = parametric_temperature_map(star_params, stars[1])

# Pre-compute polygon flux and Fourier transform matrices for each epoch.
# This stores polyflux and polyft in-place in each stellar_geometry struct.
# Must be called before image_reconstruct_oi or any chi2 evaluation.
setup_oi!(data, stars)

# Gradient-based smoothing (see "Regularizations" below for why this is preferred to "tv2")
# Format: ["type", weight, operator_info, pixel_range]
regularizers = [["sobel2", 10.0, sobel_gradient_healpix(n), 1:length(tmap_start)]]

# Run the reconstruction: iteratively adjusts pixel temperatures to fit
# the interferometric observables (V², closure phases, triple amplitudes)
tmap = image_reconstruct_oi(tmap_start, data, stars;
                             maxiter=500, regularizers=regularizers, verbose=true)
```

## Regularizations

Multiple regularizations can be combined in the `regularizers` list. Each entry
is a vector specifying the type, weight, and any additional arguments:

| Name | Syntax | Penalises | Response |
|------|--------|-----------|----------|
| `"sobel"` | `["sobel", mu, sobel_info, pixel_range]` | `∫\|∇x\|dΩ` — isotropic L1, edge-preserving | `k²` |
| `"sobel2"` | `["sobel2", mu, sobel_info, pixel_range]` | `∫\|∇x\|²dΩ` — smooth | `k²` |
| `"tv2"` | `["tv2", mu, tv_info, pixel_range]` | `‖Lx‖²` — curvature | `k⁴` |
| `"tv"` | `["tv", mu, tv_info, pixel_range]` | `‖Lx‖` — curvature | `k⁴` |
| `"mem"` | `["mem", mu]` | per-pixel contrast (no spatial coupling) | — |
| `"mean"` | `["mean", mu]` | departure from the mean | — |
| `"bias"` | `["bias", mu, B]` | harmonic bias for asymmetric brightening | — |

Radial regularizers for the single-epoch, non-rotating case (`"radflat"`, `"radialvar"`,
`"orthold"`) are covered separately in [Radial regularizers](../api/radial_regularizers.md).

### Which to use

`"sobel"` / `"sobel2"` are built on `sobel_gradient_healpix`, a tangent-plane
least-squares gradient — the sphere's analogue of the 2-D Sobel stencil, which is itself the
inverse-distance²-weighted least-squares gradient on a 3×3 grid. They penalise the **first**
derivative.

`"tv"` / `"tv2"` are built on `tv_neighbors_healpix`, whose operator is the graph
Laplacian (diagonal = neighbour count, `−1` off-diagonal), so both penalise the **second**
derivative. Note that despite the name, `"tv"` on HEALPix is `‖Lx‖`, a single global norm —
a monotone transform of `"tv2"`, not an edge-preserving L1 total variation. For that, use
`"sobel"`.

The `k⁴` response of `"tv"`/`"tv2"` discriminates very sharply by scale: it barely touches
large structure and crushes fine structure, so there is little usable middle ground between
a spiky map and a nearly constant one. `k²` rolls off gently enough to suppress
tessel-scale noise while leaving genuine structure.

Two further differences matter when transferring a weight between targets or resolutions:

* `"sobel"`, `"sobel2"` and `"mem"` are **scale-invariant** (normalised by `mean(x)`), which
  matches χ² — exactly invariant under `x → c·x`, since the visibilities are flux-normalised.
  `"tv"` and `"tv2"` are not, so their effective strength rides on the map's arbitrary
  overall level: the same weight is ~42× stronger on a 25000 K star than on a 3900 K one.
* `"sobel"` and `"sobel2"` carry the `4π/npix` solid angle, so on a smooth map their value
  converges as the mesh refines. `"tv2"` falls by ~×0.4 per `nside` doubling, so a weight
  tuned at one resolution does not transfer to another.

Build the operator argument with `sobel_gradient_healpix(n)` for the Sobel forms, or
`tv_neighbors_healpix(n)` / `tv_neighbors_longlat(ntheta, nphi)` for the Laplacian ones.

### Measured comparison

Each prior at its own best weight, from the ladder in `demos/rwcep_radflat.jl`
(`ALPHAS=0 TVWS=0,1e1,1e2,1e3,1e4`, RW Cep at nside 3, 2819 data against 428 patches, no
radial regularizer so the smoothing prior is the only thing acting):

| prior | best weight | χ²ᵥ₂/n | χ²ₜ₃ₚ/n |
|---|---|---|---|
| `"tv2"` | 1e−4 | 2.29 | 1.7 |
| *(none)* | — | 2.00 | 1.5 |
| `"sobel2"` | 1e1 | 1.93 | 1.5 |
| `"sobel"` | 1e1 | **1.86** | **1.4** |

`"tv2"` fits **worse than using no prior at all**, which is the practical consequence of the
`k⁴` response: by the time it is strong enough to suppress tessel-scale noise it is already
flattening structure the closure phases can see. `"sobel"` is best on both observables — the
L1 form suppresses noise without penalising a genuine edge in proportion to its height.

Read the weights as orders of magnitude, not as constants. They are `O(1–10²)` for the Sobel
forms because those are normalised by `mean(x)²` and by solid angle — so `tpole` and `nside`
drop out — but the balance against χ² still scales with the **number of data points**, so
re-run the ladder for a new dataset rather than carrying a weight across.

!!! warning "Weights are not comparable between the two families"
    A `"tv2"` weight is `O(10⁻⁵–10⁻⁴)` and a `"sobel2"` weight `O(1–10²)` on the same data —
    four to six orders of magnitude apart, because `"tv2"` rides on the map's absolute
    temperature. Substituting one name for the other without re-tuning silently produces
    either no regularization or a prior that dominates the likelihood.

## Evaluating the fit

After reconstruction, evaluate the fit quality:

```julia
# Total criterion (chi2 + regularization)
crit = image_reconstruct_oi_crit(tmap, data, stars; regularizers=regularizers, verbose=true)

# Chi2 only (no regularization)
chi2 = image_reconstruct_oi_chi2(tmap, data, stars; verbose=true)

# Per-observable chi2 for a single epoch
chi2_v2, chi2_t3amp, chi2_t3phi = chi2s(tmap, stars[1], data[1]; verbose=true)

# Model observables
v2_model, t3amp_model, t3phi_model = observables(tmap, stars[1], data[1])
```

## Keyword reference

`image_reconstruct_oi` keywords:

| Keyword | Default | Description |
|---------|---------|-------------|
| `maxiter` | `200` | Maximum optimizer iterations |
| `lower` | `0` | Lower bound on pixel values |
| `upper` | `Inf` | Upper bound on pixel values |
| `regularizers` | `[]` | List of regularization terms |
| `epochs_weights` | `[]` | Per-epoch weights (empty = equal) |
| `verbose` | `true` | Print per-iteration diagnostics |

## Chi-squared functions

Two implementations of the forward model and gradient are available:

| Function | Description |
|----------|-------------|
| `spheroid_chi2_fg(x, g, star, data)` | Matrix-based: precomputes polyft matrix (fast for small problems) |
| `fused_spheroid_chi2_fg(x, g, star, data)` | Matrix-free: computes FT on-the-fly (memory-efficient for large problems) |

Both are drop-in replacements for each other and produce identical results.

## Soft visibility

ROTIR uses a sigmoid-based soft visibility to smoothly weight pixels near the
limb, replacing the traditional hard mask (`normals_z > 0`). This makes the
chi-squared differentiable with respect to shape parameters.

```
w(p) = sigmoid(kappa * nz(p))
```

where `nz` is the z-component of the surface normal and `kappa` (default 50)
controls the sharpness of the transition.

## Joint shape + map optimization

ROTIR can simultaneously optimize the surface map and shape parameters
(radii, inclination, position angle) using analytical gradients.

!!! note
    `tessellation_healpix` defaults to `Float32`. For joint reconstruction,
    pass `T=Float64` to match the optimizer's precision:
    `tessellation_healpix(n, T=Float64)`.

```julia
# Starting shape parameters: [rpole (mas), omega, inclination (°), PA (°)]
θ_start = [1.37, 0.9, 78.0, 24.0]

tmap_final, θ_final = joint_reconstruct_oi(
    tmap_start, θ_start, data, tessels, star_params, tepochs;
    maxiter_xmap=200,     # iterations per temperature-map step
    maxiter_θ=50,         # iterations per shape-parameter step
    nouter=5,             # number of alternating cycles
    reg_weight=1e-5,      # TV regularization weight
    κ=50.0,               # soft-visibility sigmoid sharpness
    θ_lower=[0.5, 0.0, 0.0, -180.0],    # [rpole, ω, inc, PA] lower bounds
    θ_upper=[3.0, 1.0, 180.0, 180.0],   # [rpole, ω, inc, PA] upper bounds
)
```

This alternates between:
1. Optimizing the temperature map with fixed shape
2. Optimizing shape parameters with fixed map

The shape parameter vector depends on `surface_type`:

| Surface type | Parameters |
|-------------|------------|
| 0 (Sphere) | `[radius, inclination, PA]` |
| 1 (Ellipsoid) | `[rx, ry, rz, inclination, PA]` |
| 2 (Rapid Rotator) | `[rpole, omega, inclination, PA]` |

## Single-epoch imaging

Most ROTIR examples use multiple epochs to exploit stellar rotation, but many
science cases only have — or need — a single interferometric snapshot:

- **Slow rotators** (Cepheids, supergiants) where the rotation period is much
  longer than the observing baseline.
- **Snapshot surveys** where only one night of data is available.
- **Symmetric stars** where the goal is limb-darkening or diameter fitting
  rather than surface mapping.

The workflow is identical to multi-epoch imaging, with `tepochs = [0.0f0]`:

```julia
using ROTIR

# Read a single OIFITS file — returns a 2D array: data_all[wavelength_bin, epoch]
data_all = readoifits_multiepochs(["polaris.oifits"], T=Float32)

# Select the first wavelength bin (single epoch → data and stars are length-1 arrays)
data = data_all[1, :]

# Single epoch: no rotation phase, so tepochs is just [0.0]
tepochs = [0.0f0]

# Tessellate the stellar surface using nested HEALPix with resolution level 4
# (nside=2^4=16, giving 3072 equal-area pixels)
n = 4
tessels = tessellation_healpix(n)

# Define a spherical stellar model (surface_type=0)
star_params = (
    surface_type    = 0,       # 0=sphere (no oblateness or tidal distortion)
    radius          = 1.6,     # angular radius in milliarcseconds
    tpole           = 6000.0,  # effective temperature in Kelvin
    ldtype          = 3,       # Hestroffer power-law limb darkening
    ld1             = 0.24,    # limb-darkening coefficient
    ld2             = 0.0,     # second coefficient (unused for Hestroffer)
    inclination     = 90.0,    # pole-on; arbitrary for a sphere
    position_angle  = 0.0,     # rotation axis PA; arbitrary for a sphere
    rotation_period = 1.0,     # irrelevant for single epoch (phase = 2π·t/P = 0)
)

# Build projected, rotated stellar geometry
# (for a sphere at t=0, this just sets up the surface normals and limb darkening)
stars = create_star_multiepochs(tessels, star_params, tepochs)

# Pre-compute polygon Fourier transform matrices for fast χ² evaluation
setup_oi!(data, stars)

# Generate initial temperature map (uniform for a sphere with no gravity darkening)
tmap_start = parametric_temperature_map(star_params, stars[1])

# Gradient smoothing — use a higher weight than multi-epoch, since a single snapshot
# has sparser UV coverage and needs more regularization
regularizers = [["sobel2", 100.0, sobel_gradient_healpix(n), 1:length(tmap_start)]]

# Run the reconstruction
tmap = image_reconstruct_oi(tmap_start, data, stars;
    maxiter=500, regularizers=regularizers, verbose=true)

# Evaluate fit quality and plot the result
chi2 = image_reconstruct_oi_chi2(tmap, data, stars; verbose=true)
plot2d(tmap, stars[1], intensity=true, compass=true)
```

### Key differences from multi-epoch

| Aspect | Single epoch | Multi-epoch |
|--------|-------------|-------------|
| `tepochs` | `[0.0f0]` | `[0.0f0, 3.5f0, ...]` |
| Rotation | Not used — `rotation_period` is irrelevant | Drives phase coverage |
| Inclination/PA | Still define the projected geometry | Same |
| `stars` array | Length 1 | Length N |
| `data` array | Length 1 | Length N |
| UV coverage | Limited to single night | Improves with epochs |

### Choosing parameters for non-rotating targets

For a single-epoch sphere, several parameters are arbitrary:

- **`rotation_period`**: Set to any nonzero value (e.g. `1.0`). It only affects
  the rotation phase `2π · t / P`, and with `t = 0` this is always zero.
- **`position_angle`**: For a sphere, PA only rotates the projected
  coordinate frame. Set it to `0.0` unless you have prior knowledge.
- **`inclination`**: For a sphere, inclination does not change the projected
  shape. Set it to `90.0` (equator-on) for interpretability, or to a known
  value if the star has a measured spin axis.

For non-spherical surfaces (ellipsoid, rapid rotator), inclination and PA
are physically meaningful even in single-epoch mode — they determine the
projected shape and limb-darkening pattern.

### Regularization tips for single-epoch data

With only one epoch, UV coverage is sparser than multi-epoch datasets.
Regularization plays a larger role:

- **`tv2`** (quadratic TV) is a good default — it smooths the map while
  preserving large-scale structure like hot/cool spots.
- Increase the regularization weight (e.g. `1e-3` to `1e-2`) compared to
  multi-epoch reconstructions, since there is less data to constrain the map.
- **`lower` / `upper` bounds** on pixel values can prevent unphysical
  temperatures. For example, `lower=3000` for a 6000 K star.

See [`demos/polaris_imaging.jl`](https://github.com/fabienbaron/ROTIR.jl/blob/main/demos/polaris_imaging.jl)
for a complete single-epoch reconstruction of Polaris.
