# ROTIR: Regularized Imaging of Stellar Surfaces

|     **Status**                  | **Documentation**               | **License**                     |**Build**                      |
|:--------------------------------|:--------------------------------|:--------------------------------|:------------------------------|
| [![][proj-img]][proj-url] | [![][doc-dev-img]][doc-dev-url] | [![][license-img]][license-url] | [![][build-img]][build-url] |

[proj-img]: http://www.repostatus.org/badges/latest/active.svg
[proj-url]: http://www.repostatus.org/#active

[doc-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[doc-dev-url]: https://fabienbaron.github.io/ROTIR.jl/dev

[license-url]: ./LICENSE
[license-img]: http://img.shields.io/badge/license-GPL3-brightgreen.svg?style=flat

[build-img]: https://github.com/fabienbaron/ROTIR.jl/workflows/Documentation/badge.svg
[build-url]: https://github.com/fabienbaron/ROTIR.jl/actions

ROTIR is a Julia package for regularized imaging of stellar surfaces from
optical interferometry data, developed by Prof. Fabien Baron (Georgia State
University) and collaborators. It reconstructs temperature maps on tessellated
stellar surfaces by fitting interferometric observables (V², closure phases,
triple amplitudes).

### **[:book: Full Documentation](https://fabienbaron.github.io/ROTIR.jl/dev)**

## Installation

```julia
using Pkg
pkg"registry add General"
pkg"registry add https://github.com/emmt/EmmtRegistry"
Pkg.add(url="https://github.com/fabienbaron/OITOOLS.jl.git")
Pkg.add(url="https://github.com/fabienbaron/ROTIR.jl.git")
using ROTIR
```

See the [installation guide](https://fabienbaron.github.io/ROTIR.jl/dev/install/) for details.

## Quick start

```julia
using ROTIR

# Load multi-epoch OIFITS data
data_all = readoifits_multiepochs(["epoch1.oifits", "epoch2.oifits"])
data = data_all[1, :]
tepochs = [d.mean_mjd for d in data] .- data[1].mean_mjd

# Create HEALPix tessellation and rapid rotator geometry
tessels = tessellation_healpix(3)
star_params = (
    surface_type=2, rpole=1.37, tpole=4800.0, ldtype=3, ld1=0.23, ld2=0.0,
    inclination=78.0, position_angle=24.0, rotation_period=54.8,
    beta=0.08, frac_escapevel=0.9, B_rot=0.0,
)
stars = create_star_multiepochs(tessels, star_params, tepochs)
tmap_start = parametric_temperature_map(star_params, stars[1])
setup_oi!(data, stars)

# Reconstruct temperature map
regularizers = [["tv2", 1e-5, tv_neighbours_healpix(3), 1:length(tmap_start)]]
tmap = image_reconstruct_oi(tmap_start, data, stars;
                             maxiter=500, regularizers=regularizers)

# Visualize
plot2d_allepochs(tmap, stars)
plot_mollweide(tmap, stars[1])
```

## Features

- **Multiple surface geometries**: spheres, triaxial ellipsoids, rapid rotators
  (centrifugally distorted), and Roche-lobe-filling stars in binaries
- **Two tessellation schemes**: HEALPix (equal-area, hierarchical) and
  longitude/latitude grids
- **Multi-epoch reconstruction**: simultaneously fit data from multiple rotation
  phases to recover the full surface map
- **Multi-resolution imaging**: coarse-to-fine HEALPix pyramid for robust
  convergence
- **Joint shape + map optimization**: analytical gradients for shape parameters
  (radii, inclination, position angle) alongside the surface map
- **Matrix-free polygon Fourier transform**: fused forward/adjoint passes with
  O(Nuv + Npix) memory instead of O(Nuv * Npix)
- **Regularization**: total variation (L1, quadratic), maximum entropy, mean
  constraint, and harmonic bias
- **Gradient-based optimization**: VMLMB quasi-Newton with bounds
  (OptimPackNextGen)

ROTIR uses [OITOOLS.jl](https://github.com/fabienbaron/OITOOLS.jl) for OIFITS
I/O and data handling.

## Documentation

Full documentation is available at
[fabienbaron.github.io/ROTIR.jl](https://fabienbaron.github.io/ROTIR.jl/dev),
including:

- [Overview and workflow](https://fabienbaron.github.io/ROTIR.jl/dev/guides/overview/)
- [Tessellation schemes](https://fabienbaron.github.io/ROTIR.jl/dev/guides/tessellation/)
- [Surface types](https://fabienbaron.github.io/ROTIR.jl/dev/guides/surfaces/)
- [Image reconstruction](https://fabienbaron.github.io/ROTIR.jl/dev/guides/reconstruction/)
- [Multi-resolution imaging](https://fabienbaron.github.io/ROTIR.jl/dev/guides/multires/)
- [API reference](https://fabienbaron.github.io/ROTIR.jl/dev/api/chi2/)

## Development install

```julia
using Pkg
pkg"registry add General"
pkg"registry add https://github.com/emmt/EmmtRegistry"
Pkg.develop(url="https://github.com/fabienbaron/ROTIR.jl.git")
```
