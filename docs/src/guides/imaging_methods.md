# Direct imaging methods

ROTIR provides three ways to produce a real-space image from the projected
stellar surface model. All three take the same inputs — polygon vertex
coordinates and intensity weights — but differ in how they compute the image.

## Overview

| Method | Function | Approach | Cost |
|--------|----------|----------|------|
| **Polyft + irfft** | `setup_polyft_single` + `irfft` | Dense matrix multiply in Fourier space | O(Nuv x Npix) memory |
| **Rasterization** | `rasterize_polygon_image` | Exact polygon-pixel clipping in image space | O(Npix x patch_area) |
| **NFFT** | `polyft_nfft_image` | Gauss-Legendre quadrature + adjoint NFFT | O(Ns log N) |

The **polyft matrix** is ROTIR's original method. It precomputes the exact
polygon Fourier transform at every UV point and stores it as a dense matrix.
This is fast for chi-squared evaluation (one BLAS mat-vec per epoch) but
requires O(Nuv x Npix) memory and is not needed when only a real-space image
is desired.

**Rasterization** works entirely in image space. Each projected quadrilateral
is clipped against every pixel it overlaps using Sutherland-Hodgman clipping,
and the exact intersection area determines the pixel's contribution. This is
geometrically exact — no Fourier ringing or truncation — and parallelized
across polygons.

The **NFFT** method approximates each polygon's Fourier integral using
Gauss-Legendre quadrature and folds all the non-uniform point evaluations
into a single adjoint NFFT. The result is the visibility function on a
regular Fourier grid, which can be inverse-FFT'd to produce the image. This
approach is fastest for large grids and also provides access to the
intermediate Fourier-domain representation.

## When to use which

| Scenario | Recommended |
|----------|-------------|
| Interferometric reconstruction (fitting UV data) | polyft matrix or fused path |
| Quick image preview | `rasterize_polygon_image` |
| Large image grids (nx > 256) | `polyft_nfft_image` |
| Need visibility amplitudes on a grid | `polyft_nfft_forward` |
| Image-plane fitting (e.g. direct imaging) | `rasterize_polygon_image` + `rasterize_adjoint!` |

## Rasterization

Rasterization produces a geometrically exact image by computing the
intersection area between each projected polygon and each image pixel.

```julia
using ROTIR

# After creating the star and temperature map (see Overview)
indx = star.index_quads_visible
pw = star.proj_west[indx, :]     # (nvis, 4) in mas
pn = star.proj_north[indx, :]    # (nvis, 4) in mas
x_weighted = tmap[indx] .* star.vis_weights[indx]

# Rasterize onto a 128x128 grid with 0.015 mas/pixel
nx = 128
pixsize = 0.015f0
img = rasterize_polygon_image(pw, pn, x_weighted, pixsize, nx)
```

The in-place variant `rasterize_polygon_image!` lets you pre-allocate the
output and control the center pixel via `cx` and `cy` keyword arguments.

For gradient-based optimization, `rasterize_adjoint!` provides the exact
transpose operation: given a residual image, it maps the gradient back to
polygon intensities weighted by the same overlap areas as the forward pass.

## NFFT-based imaging

The NFFT method approximates the polygon Fourier transform using
Gauss-Legendre quadrature inside each quadrilateral, then evaluates all
the resulting non-uniform exponential sums in a single adjoint NFFT.

### One-step image

```julia
img = polyft_nfft_image(pw, pn, x_weighted, pixsize, nx; ngauss=6)
```

### Step-by-step Fourier pipeline

For more control, or to inspect the visibilities, use the two-step approach
with `rfftfreq` / `fftfreq`:

```julia
using FFTW

# Step 1: Build the frequency grid (cycles/mas)
u_freq = rfftfreq(nx, 1 / pixsize)   # non-negative freqs, length nx÷2+1
v_freq = fftfreq(nx, 1 / pixsize)    # full freqs,         length nx

# Step 2: Compute complex visibilities on the rfft grid
F = polyft_nfft_forward(pw, pn, x_weighted, pixsize, nx; ngauss=6)
# F is (nx÷2+1, nx) complex — standard rfft layout:
#   rows match rfftfreq (non-negative u)
#   cols match fftfreq  (standard FFT order)

# Step 3: Inspect visibilities (fftshift for centered display)
using PyPlot
imshow(log10.(abs.(fftshift(F, 2))), origin="lower", cmap="inferno")

# Step 4: Inverse real-FFT to recover the image
img = fftshift(irfft(F, nx))
```

The frequency grid has spacing `1 / (nx * pixsize)` cycles/mas and extends
to the Nyquist frequency `1 / (2 * pixsize)` cycles/mas.

### Quadrature accuracy

The `ngauss` parameter controls the number of Gauss-Legendre quadrature
points per axis inside each polygon. Higher orders are more accurate but
slower:

| `ngauss` | Points per polygon | When to use |
|----------|-------------------|-------------|
| 2 | 4 | Sub-pixel polygons (HEALPix n >= 6) |
| 4 | 16 | Default — good for most cases |
| 6 | 36 | Larger polygons or high accuracy |
| 8-10 | 64-100 | Very large polygons or high Nyquist |

The accuracy depends on the phase span across each polygon: if the highest
spatial frequency times the polygon diagonal is small, low-order quadrature
suffices. As a rule of thumb, `ngauss=4` is accurate for HEALPix level n >= 4
with typical pixel scales.

For cases where polygons are large relative to the pixel scale, the `nsub`
parameter subdivides each polygon into `nsub x nsub` sub-squares before
applying quadrature. Subdivision converges faster than increasing `ngauss`
alone for oscillatory integrands.

## Complete example

See `demos/demo_rasterize_nfft.jl` for a complete script that compares all
three methods on a gravity-darkened rapid rotator, including the step-by-step
Fourier pipeline with `rfftfreq`/`fftfreq`.
