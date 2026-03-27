# Single-epoch imaging

Most ROTIR examples use multiple epochs to exploit stellar rotation, but many
science cases only have — or need — a single interferometric snapshot:

- **Slow rotators** (Cepheids, supergiants) where the rotation period is much
  longer than the observing baseline.
- **Snapshot surveys** where only one night of data is available.
- **Symmetric stars** where the goal is limb-darkening or diameter fitting
  rather than surface mapping.

The workflow is the same as multi-epoch imaging, but with `tepochs = [0.0f0]`.

## Minimal example

```julia
using ROTIR

# 1. Load a single OIFITS file
data_all = readoifits_multiepochs(["polaris.oifits"], T=Float32)
data = data_all[1, :]         # first wavelength bin
tepochs = [0.0f0]             # single epoch

# 2. Create geometry (sphere)
n = 4
tessels = tessellation_healpix(n)
star_params = (
    surface_type    = 0,       # sphere
    radius          = 1.6,     # mas
    tpole           = 6000.0,  # K
    ldtype          = 3,       # Hestroffer power-law
    ld1             = 0.24,
    ld2             = 0.0,
    inclination     = 90.0,    # pole-on (arbitrary for sphere)
    position_angle  = 0.0,
    rotation_period = 1.0      # irrelevant for single epoch
)

stars = create_star_multiepochs(tessels, star_params, tepochs)
setup_oi!(data, stars)

# 3. Starting map and regularization
tmap_start = parametric_temperature_map(star_params, stars[1])
regularizers = [["tv2", 1e-4, tv_neighbours_healpix(n), 1:length(tmap_start)]]

# 4. Reconstruct
tmap = image_reconstruct_oi(tmap_start, data, stars,
    maxiter=500, regularizers=regularizers, verbose=true)

# 5. Evaluate and plot
chi2 = image_reconstruct_oi_chi2(tmap, data, stars, verbose=true)
plot2d(tmap, stars[1], intensity=true, compass=true)
```

## Key differences from multi-epoch

| Aspect | Single epoch | Multi-epoch |
|--------|-------------|-------------|
| `tepochs` | `[0.0f0]` | `[0.0f0, 3.5f0, ...]` |
| Rotation | Not used — `rotation_period` is irrelevant | Drives phase coverage |
| Inclination/PA | Still define the projected geometry | Same |
| `stars` array | Length 1 | Length N |
| `data` array | Length 1 | Length N |
| UV coverage | Limited to single night | Improves with epochs |

## Choosing parameters for non-rotating targets

For a single-epoch sphere, several parameters are arbitrary:

- **`rotation_period`**: Set to any nonzero value (e.g. `1.0`). It only affects
  the rotation phase `2pi * t / P`, and with `t = 0` this is always zero.
- **`position_angle`**: For a sphere, PA only rotates the projected
  coordinate frame. Set it to `0.0` unless you have prior knowledge.
- **`inclination`**: For a sphere, inclination does not change the projected
  shape. Set it to `90.0` (equator-on) for interpretability, or to a known
  value if the star has a measured spin axis.

For non-spherical surfaces (ellipsoid, rapid rotator), inclination and PA
are physically meaningful even in single-epoch mode — they determine the
projected shape and limb-darkening pattern.

## Regularization tips

With only one epoch, UV coverage is sparser than multi-epoch datasets.
Regularization plays a larger role:

- **`tv2`** (quadratic TV) is a good default — it smooths the map while
  preserving large-scale structure like hot/cool spots.
- Increase the regularization weight (e.g. `1e-3` to `1e-2`) compared to
  multi-epoch reconstructions, since there is less data to constrain the map.
- **`lower` / `upper` bounds** on pixel values can prevent unphysical
  temperatures. For example, `lower=3000` for a 6000 K star.

## Working with the results

All analysis functions work identically for single and multi-epoch data:

```julia
# Per-observable reduced chi2
chi2_v2, chi2_t3amp, chi2_t3phi = chi2s(tmap, stars[1], data[1]; verbose=true)

# Model observables for comparison with data
v2_model, t3amp_model, t3phi_model = observables(tmap, stars[1], data[1])

# Residual plots (from OITOOLS)
plot_v2_residuals(v2_model, data[1])
plot_t3phi_residuals(t3phi_model, data[1])

# Surface plots
plot2d(tmap, stars[1], intensity=true, graticules=true, compass=true)
plot_mollweide(tmap, stars[1])
```

## Full demo

See [`demos/polaris_imaging.jl`](https://github.com/fabienbaron/ROTIR.jl/blob/main/demos/polaris_imaging.jl)
for a complete single-epoch reconstruction of Polaris.
