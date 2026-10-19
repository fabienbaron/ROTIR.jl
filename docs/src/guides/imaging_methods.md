# Direct imaging methods

Three methods produce a sky image from the projected surface. All take the projected polygon
vertices and their intensity weights.

## Overview

| Method | Function | Approach | Cost |
|--------|----------|----------|------|
| **Polyft + irfft** | `setup_polyft_single` + `irfft` | Dense matrix multiply in Fourier space | O(Nuv x Npix) memory |
| **Rasterization** | `rasterize_polygon_image` | Exact polygon-pixel clipping in image space | O(Npix x patch_area) |
| **NFFT** | `polyft_nfft_image` | Gauss-Legendre quadrature + adjoint NFFT | O(Ns log N) |

The **polyft matrix** holds the exact polygon Fourier transform at every uv point as a dense
matrix: one mat-vec per epoch for χ², O(Nuv × Npix) memory.

**Rasterization** clips each projected quadrilateral against the pixels it overlaps
(Sutherland–Hodgman) and adds the intersection areas. It is geometrically exact (no Fourier
ringing or truncation) and multithreaded over polygons.

**NFFT** approximates each polygon's Fourier integral by Gauss–Legendre quadrature and sums
all quadrature points in one adjoint NFFT. It gives the visibilities on a regular grid, which
an inverse FFT turns into the image. It is the fastest for large grids.

## When to use which

| Scenario | Recommended |
|----------|-------------|
| Reconstruction from uv data | polyft matrix or fused path |
| Quick image preview | `rasterize_polygon_image` |
| Large image grids (nx > 256) | `polyft_nfft_image` |
| Need visibility amplitudes on a grid | `polyft_nfft_forward` |
| Image-plane fitting (e.g. direct imaging) | `rasterize_polygon_image` + `rasterize_adjoint!` |

## Rasterization

```julia
using ROTIR

# star and tmap as in the Overview
indx = star.index_quads_visible
pw = star.proj_west[indx, :]     # (nvis, 4), mas
pn = star.proj_north[indx, :]    # (nvis, 4), mas
x_weighted = tmap[indx] .* star.vis_weights[indx] .* star.ldmap[indx]   # with limb darkening

# 128 × 128 pixels of 0.015 mas
nx = 128
pixsize = 0.015f0
img = rasterize_polygon_image(pw, pn, x_weighted, pixsize, nx)
```

`rasterize_polygon_image!` writes into a pre-allocated image. Both take the centre pixel as
`cx`, `cy` (default `nx ÷ 2 + 1`).

`rasterize_adjoint!` is the exact transpose: it maps a residual image back to polygon
intensities with the same overlap areas.

## NFFT-based imaging

### One-step image

```julia
img = polyft_nfft_image(pw, pn, x_weighted, pixsize, nx; ngauss=6)
```

### Step-by-step Fourier pipeline

To get the visibilities as well:

```julia
using FFTW

# Frequency grid (cycles/mas)
u_freq = rfftfreq(nx, 1 / pixsize)   # non-negative, length nx÷2+1
v_freq = fftfreq(nx, 1 / pixsize)    # length nx

# Complex visibilities, (nx÷2+1, nx), rfft layout: rows follow rfftfreq, columns fftfreq
F = polyft_nfft_forward(pw, pn, x_weighted, pixsize, nx; ngauss=6)

using PythonPlot
imshow(log10.(abs.(fftshift(F, 2))), origin="lower", cmap="inferno")

# Image
img = fftshift(irfft(F, nx))
```

The grid spacing is `1 / (nx * pixsize)` cycles/mas, up to the Nyquist frequency
`1 / (2 * pixsize)`.

### Quadrature accuracy

`ngauss` is the number of Gauss–Legendre points per axis in each polygon (default 4):

| `ngauss` | Points per polygon | When to use |
|----------|-------------------|-------------|
| 2 | 4 | Sub-pixel polygons (HEALPix n >= 6) |
| 4 | 16 | default |
| 6 | 36 | Larger polygons or high accuracy |
| 8-10 | 64-100 | Very large polygons or high Nyquist |

The accuracy depends on the phase span across a polygon (highest frequency × polygon
diagonal). `ngauss=4` is accurate for HEALPix level n ≥ 4 at typical pixel scales.

`nsub` (default 1) splits each polygon into `nsub × nsub` sub-quads before quadrature. For
polygons large relative to the pixel scale this converges faster than raising `ngauss`.

## Complete example

`demos/demo_rasterize_nfft.jl` compares rasterization and NFFT on a gravity-darkened rapid rotator,
including the Fourier steps above.
