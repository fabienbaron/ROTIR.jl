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
    beta=0.08,                # gravity-darkening exponent (T ∝ g^β)
    frac_escapevel=0.9,       # rotational velocity as fraction of escape velocity
    B_rot=0.0,                # differential rotation coefficient (0=solid body)
)

# Build projected, rotated stellar geometry for each epoch
# (computes surface shape, normals, visible pixels, and limb darkening)
stars = create_star_multiepochs(tessels, star_params, tepochs)

# Generate initial temperature map from the von Zeipel gravity-darkening law
# (serves as the starting point for the optimizer)
tmap_start = parametric_temperature_map(star_params, stars[1])

# Pre-compute polygon flux and Fourier transform matrices for each epoch
# (stored in-place in the stars objects for fast χ² evaluation)
setup_oi!(data, stars)

# Set up quadratic total-variation regularization to enforce smooth maps
# Format: ["type", weight, neighbor_info, pixel_range]
regularizers = [["tv2", 1e-5, tv_neighbors_healpix(n), 1:length(tmap_start)]]

# Run the reconstruction: iteratively adjusts pixel temperatures to fit
# the interferometric observables (V², closure phases, triple amplitudes)
tmap = image_reconstruct_oi(tmap_start, data, stars;
                             maxiter=500, regularizers=regularizers, verbose=true)
```

## Regularizations

Multiple regularizations can be combined in the `regularizers` list. Each entry
is a vector specifying the type, weight, and any additional arguments:

| Name | Syntax | Description |
|------|--------|-------------|
| `"tv"` | `["tv", mu, tv_info, pixel_range]` | Total variation (L1) |
| `"tv2"` | `["tv2", mu, tv_info, pixel_range]` | Total variation squared (quadratic) |
| `"mem"` | `["mem", mu]` | Maximum entropy |
| `"mean"` | `["mean", mu]` | Mean constraint |
| `"bias"` | `["bias", mu, B]` | Harmonic bias for asymmetric brightening |

The `tv_info` argument comes from `tv_neighbors_healpix(n)` or
`tv_neighbors_longlat(ntheta, nphi)`.

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
(radii, inclination, position angle) using analytical gradients:

```julia
# Starting shape parameters: [rpole (mas), omega, inclination (°), PA (°)]
theta_start = [1.37, 0.9, 78.0, 24.0]

tmap_final, theta_final = joint_reconstruct_oi(
    tmap_start, theta_start, data, tessels, star_params, tepochs;
    maxiter_xmap=200,     # iterations per temperature-map step
    maxiter_theta=50,     # iterations per shape-parameter step
    nouter=5,             # number of alternating cycles
    reg_weight=1e-5,      # TV regularization weight
    kappa=50.0,           # soft-visibility sigmoid sharpness
    theta_lower=[0.5, 0.0, 0.0, -180.0],    # [rpole, omega, inc, PA] lower bounds
    theta_upper=[3.0, 1.0, 180.0, 180.0],   # [rpole, omega, inc, PA] upper bounds
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

# Quadratic total-variation regularization — use a higher weight than multi-epoch
# since a single snapshot has sparser UV coverage and needs more regularization
regularizers = [["tv2", 1e-4, tv_neighbors_healpix(n), 1:length(tmap_start)]]

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
