# Rasterization & NFFT

Two alternative forward-model implementations that avoid building the dense
polyft matrix. Both produce real-space images from polygon vertex coordinates
and intensity weights.

See [Direct imaging methods](@ref) for a guide with usage examples and
guidance on when to use each approach.

## Rasterization

Exact polygon-to-pixel rendering via Sutherland-Hodgman clipping. Works
entirely in image space — no Fourier transforms involved.

| Function | Description |
|----------|-------------|
| `rasterize_polygon_image(proj_west, proj_north, x_weighted, pixsize, nx; cx, cy, T)` | Allocate and return an `(nx, nx)` rasterized image |
| `rasterize_polygon_image!(img, proj_west, proj_north, x_weighted, pixsize; cx, cy)` | In-place rasterization into preallocated `img` |
| `rasterize_adjoint!(grad, res_img, proj_west, proj_north, pixel_weight, pixsize; cx, cy)` | Adjoint (transpose) of rasterization for gradient computation |

### `rasterize_polygon_image` arguments

| Argument | Type | Description |
|----------|------|-------------|
| `proj_west` | `(Npix, 4)` matrix | Polygon vertex West coordinates (mas) |
| `proj_north` | `(Npix, 4)` matrix | Polygon vertex North coordinates (mas) |
| `x_weighted` | `(Npix,)` vector | Intensity weights (e.g. `tmap .* vis_weights`) |
| `pixsize` | scalar | Pixel scale (mas/pixel) |
| `nx` | integer | Image size (pixels per side) |

| Keyword | Default | Description |
|---------|---------|-------------|
| `cx` | `nx÷2+1` | Column index of the projection origin |
| `cy` | `nx÷2+1` | Row index of the projection origin |
| `T` | `Float32` | Element type of output image |

### `rasterize_adjoint!` details

The adjoint operation is the exact transpose of the forward rasterization:
the forward spreads polygon intensities into pixels weighted by overlap area;
the adjoint gathers pixel values back to polygons weighted by the same
overlap area.

```
grad[p] = pixel_weight[p] * pixsize^2 * sum_{iy,ix} res_img[iy,ix] * Area(tile_p ∩ pixel)
```

Both forward and adjoint are threaded over polygons using per-thread local
buffers, reduced at the end.

### Internal helpers

| Function | Description |
|----------|-------------|
| `quad_box_area(q1x,q1y,...,q4x,q4y, xmin,xmax,ymin,ymax, Ax,Ay,Bx,By)` | Exact area of intersection between a quadrilateral and an axis-aligned box via Sutherland-Hodgman clipping + shoelace formula |

## NFFT

Fast polygon Fourier transform via Gauss-Legendre quadrature and the
Non-uniform Fast Fourier Transform. Approximates the continuous polygon FT
integral by placing quadrature points inside each quadrilateral and
evaluating the resulting non-uniform exponential sum with a single adjoint
NFFT.

| Function | Description |
|----------|-------------|
| `polyft_nfft_forward(proj_west, proj_north, x_weighted, pixsize, nx; ngauss, nsub, T)` | Complex visibilities on an rfft grid — returns `(nx÷2+1, nx)` array |
| `polyft_nfft_image(proj_west, proj_north, x_weighted, pixsize, nx; ngauss, nsub, T)` | Convenience: `fftshift(irfft(polyft_nfft_forward(...), nx))` |
| `build_gauss_samples(proj_west, proj_north, x_weighted; ngauss, nsub, T)` | Build Gauss-Legendre quadrature samples for all polygons |

### `polyft_nfft_forward` arguments

| Argument | Type | Description |
|----------|------|-------------|
| `proj_west` | `(Npix, 4)` matrix | Polygon vertex West coordinates (mas) |
| `proj_north` | `(Npix, 4)` matrix | Polygon vertex North coordinates (mas) |
| `x_weighted` | `(Npix,)` vector | Intensity weights |
| `pixsize` | scalar | Pixel scale (mas/pixel) |
| `nx` | integer | Grid size (pixels per side) |

| Keyword | Default | Description |
|---------|---------|-------------|
| `ngauss` | `4` | Gauss-Legendre order per axis (2, 3, 4, 5, 6, 8, or 10) |
| `nsub` | `1` | Subdivision level per axis (total quadrature points = `(ngauss*nsub)^2` per polygon) |
| `T` | `Float32` | Floating-point type |

### Return value layout

`polyft_nfft_forward` returns a `Complex{T}` array of shape `(nx÷2+1, nx)` in
the standard rfft layout:

- **Rows** correspond to `rfftfreq(nx, 1/pixsize)` — non-negative spatial
  frequencies from 0 to the Nyquist frequency `1/(2*pixsize)` cycles/mas.
- **Columns** correspond to `fftfreq(nx, 1/pixsize)` — the full frequency
  axis in standard FFT order (zero-frequency at index 1).

This layout is compatible with Julia's `irfft(F, nx)` to recover the
real-space image. Apply `fftshift` to center the image:

```julia
img = fftshift(irfft(F, nx))
```

### Convention notes

The NFFT.jl adjoint computes `sum_j f_j exp(+2*pi*i * k . x_j)`. By
negating the sample positions (`pos = -r/L` where `L = nx*pixsize`), the
adjoint directly produces the standard DFT convention
`F[k] = sum f_j exp(-2*pi*i * k . r_j)`. The `ifftshift` then converts
from NFFT natural order (`k = -N/2, ..., N/2-1`) to FFT order
(`k = 0, 1, ..., N/2-1, -N/2, ..., -1`), and extracting the first
`nx÷2+1` rows gives the rfft layout.

### `build_gauss_samples` details

For each visible polygon, `build_gauss_samples` parameterizes the
quadrilateral using bilinear shape functions on `[-1,1]^2`, optionally
subdivides into `nsub x nsub` sub-squares, and places `ngauss x ngauss`
tensor-product Gauss-Legendre nodes in each sub-square.

Returns `(xs, ys, fs)` where:
- `xs, ys`: spatial positions (mas) of the quadrature points
- `fs`: complex weights including the intensity, Jacobian, and Gauss weight

Total samples: `(ngauss * nsub)^2 * Npix`.
