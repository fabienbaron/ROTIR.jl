# Overview

ROTIR reconstructs stellar surface temperature maps from optical interferometry data.

## Concepts

**Tessellation.** The surface is divided into quadrilateral pixels on a HEALPix or
longitude/latitude grid. Each pixel has 4 corners and a centre, stored in Cartesian and
spherical coordinates on the unit sphere.

**Surface geometry.** The unit-sphere tessellation is scaled to the stellar shape: sphere,
triaxial ellipsoid, rapid rotator, or Roche-lobe-filling binary component. The geometry sets
the gravity darkening, by the von Zeipel or Espinosa Lara–Rieutord law (`gravity_law`).

**Rotation and projection.** At each epoch the star is rotated by its spin phase, inclination
and position angle, and projected onto the sky. Pixels facing away from the observer are
masked with a sigmoid weight on the normal, so the mask is differentiable.

**Polygon Fourier transform.** Each projected pixel is a quadrilateral, whose Fourier
transform at each uv point is analytic (edge-based sinc integrals). The complex visibility is
the flux-weighted sum over visible pixels.

**Direct imaging.** Two methods produce sky images without the dense polyft matrix:
rasterisation (exact polygon–pixel clipping, Sutherland–Hodgman) and NFFT (Gauss–Legendre
quadrature on a non-uniform FFT). See [Direct imaging methods](@ref).

**Reconstruction.** The temperature map minimises the χ² of V², closure phases and triple
amplitudes plus regularisation, with the VMLMB quasi-Newton optimiser and analytic gradients.

## Minimal workflow

```julia
using ROTIR

# 1. Load multi-epoch OIFITS data
oifitsfiles = ["epoch1.oifits", "epoch2.oifits", "epoch3.oifits"]
data_all = readoifits_multiepochs(oifitsfiles)
data = data_all[1, :]                            # first wavelength bin, all epochs
tepochs = [d.mean_mjd for d in data]
tepochs = tepochs .- tepochs[1]                  # relative MJDs

# 2. Create tessellation
n = 3                                            # HEALPix level (npix = 768)
tessels = tessellation_healpix(n)

# 3. Define stellar parameters
star_params = (
    surface_type    = 2,        # 0=sphere, 1=ellipsoid, 2=rapid rotator, 3=Roche
    rpole           = 1.37,     # polar radius (mas)
    tpole           = 4800.0,   # polar temperature (K)
    ldtype          = 3,        # limb darkening: 0 none, 1 linear, 2 quadratic, 3 power law, 4 Claret
    ld1             = 0.23,     # first LD coefficient
    ld2             = 0.0,      # second LD coefficient
    inclination     = 78.0,     # degrees
    position_angle  = 24.0,     # degrees
    rotation_period = 54.8,     # days
    beta            = 0.08,     # gravity darkening, T ∝ g^β (0.25 radiative, 0.08 convective)
    frac_escapevel  = 0.9,      # Ω / Ω_crit (rapid rotator)
    B_rot           = 0.0,      # differential rotation (0 = solid body)
)

# 4. Build geometry for all epochs
stars = create_star_multiepochs(tessels, star_params, tepochs)

# 5. Starting map: the gravity-darkened temperature
tmap_start = parametric_temperature_map(star_params, stars[1])

# 6. Polygon flux and Fourier transform, stored in each star
setup_oi!(data, stars)

# 7. Set up regularization
regularizers = [["sobel2", 10.0, sobel_gradient_healpix(n), 1:length(tmap_start)]]

# 8. Reconstruct
tmap = image_reconstruct_oi(tmap_start, data, stars;
                             maxiter=500, regularizers=regularizers, verbose=true)

# 9. Plot
plot2d_allepochs(tmap, stars)
plot_mollweide(tmap, stars[1])
```

## Data flow

```
OIFITS files
    |
    v
readoifits_multiepochs()  -->  Matrix{OIdata}  [wavelength bin, epoch]
    |
    v
tessellation_healpix(n)   -->  tessellation     (unit sphere grid)
    |
    v
create_star_multiepochs() -->  Vector{stellar_geometry}
    |                           (rotated, projected, visibility-masked)
    v
setup_oi!()               -->  polyflux, polyft matrices stored in stars
    |
    v
image_reconstruct_oi()    -->  temperature map (Vector, length npix)
    |
    v
plot2d / plot_mollweide   -->  visualization
```

## Example output

Gravity-darkened temperature map of a rapid rotator, Mollweide projection:

![Mollweide output](../assets/plot_mollweide.png)
