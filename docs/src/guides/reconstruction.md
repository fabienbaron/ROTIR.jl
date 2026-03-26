# Image reconstruction

ROTIR reconstructs stellar surface temperature maps by minimizing the
chi-squared between model and observed interferometric data, subject to
regularization. The optimizer is VMLMB (variable metric limited-memory
quasi-Newton with bounds) from OptimPackNextGen.jl.

## Basic reconstruction

```julia
using ROTIR

# Load data
oifitsfiles = ["epoch1.oifits", "epoch2.oifits"]
data_all = readoifits_multiepochs(oifitsfiles)
data = data_all[1, :]
tepochs = [d.mean_mjd for d in data] .- data[1].mean_mjd

# Create geometry
n = 3
tessels = tessellation_healpix(n)
star_params = (
    surface_type=2, rpole=1.37, tpole=4800.0, ldtype=3, ld1=0.23, ld2=0.0,
    inclination=78.0, position_angle=24.0, rotation_period=54.8,
    beta=0.08, frac_escapevel=0.9, B_rot=0.0,
)
stars = create_star_multiepochs(tessels, star_params, tepochs)

# Starting map and setup
tmap_start = parametric_temperature_map(star_params, stars[1])
setup_oi!(data, stars)

# Regularization
regularizers = [["tv2", 1e-5, tv_neighbours_healpix(n), 1:length(tmap_start)]]

# Reconstruct
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

The `tv_info` argument comes from `tv_neighbours_healpix(n)` or
`tv_neighbours_longlat(ntheta, nphi)`.

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
theta_start = [1.37, 0.9, 78.0, 24.0]   # [rpole, omega, inc, PA] for rapid rotator

tmap_final, theta_final = joint_reconstruct_oi(
    tmap_start, theta_start, data, tessels, star_params, tepochs;
    maxiter_xmap=200, maxiter_theta=50, nouter=5,
    reg_weight=1e-5, kappa=50.0,
    theta_lower=[0.5, 0.0, 0.0, -180.0],
    theta_upper=[3.0, 1.0, 180.0, 180.0],
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
