# Rasterization & NFFT

Two ways to make a real-space image from the projected tessels and their intensity weights,
without the dense polyft matrix. See [Direct imaging methods](@ref) for examples.

## Rasterization

Each pixel receives the exact overlap area of each tessel (Sutherland–Hodgman clipping)
times its intensity. No Fourier transform is involved.

| Function | Description |
|----------|-------------|
| `rasterize_polygon_image(proj_west, proj_north, x_weighted, pixsize, nx; cx, cy, T)` | Allocate and return an `(nx, nx)` rasterized image |
| `rasterize_polygon_image!(img, proj_west, proj_north, x_weighted, pixsize; cx, cy)` | In-place rasterization into preallocated `img` |
| `rasterize_adjoint!(grad, res_img, proj_west, proj_north, pixel_weight, pixsize; cx, cy)` | Adjoint (transpose) of the rasterization, for gradients |

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
| `T` | `float(real(eltype(x_weighted)))` | Element type of the image |

The image is indexed `img[iy, ix]`; `proj_west = proj_north = 0` falls on pixel `(cx, cy)`.

### `rasterize_adjoint!`

The exact transpose of the rasterization: it gathers pixel values back to the tessels,
weighted by the same overlap areas.

```
grad[p] = pixel_weight[p] * pixsize^2 * sum_{iy,ix} res_img[iy,ix] * Area(tile_p ∩ pixel)
```

Both are threaded over tessels.

### Internal helpers

| Function | Description |
|----------|-------------|
| `quad_box_area(q1x,q1y,...,q4x,q4y, xmin,xmax,ymin,ymax, Ax,Ay,Bx,By)` | Area of the intersection of a quadrilateral and an axis-aligned box (Sutherland–Hodgman clipping, shoelace formula); `Ax,Ay,Bx,By` are scratch buffers of at least 12 elements |

## NFFT

The polygon Fourier transform on a regular frequency grid: Gauss–Legendre quadrature points
inside each tessel, summed with one adjoint NFFT.

| Function | Description |
|----------|-------------|
| `polyft_nfft_forward(proj_west, proj_north, x_weighted, pixsize, nx; ngauss, nsub, T, fftflags)` | Complex visibilities on an rfft grid, a `(nx÷2+1, nx)` array |
| `polyft_nfft_image(proj_west, proj_north, x_weighted, pixsize, nx; ngauss, nsub, T, fftflags)` | `fftshift(irfft(polyft_nfft_forward(...), nx))` |
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
| `ngauss` | `4` | Gauss–Legendre order per axis: 2, 3, 4, 5, 6, 8 or 10 |
| `nsub` | `1` | Sub-squares per axis; `(ngauss*nsub)^2` points per tessel |
| `T` | `float(real(eltype(proj_west)))` | Floating-point type |
| `fftflags` | `FFTW.MEASURE` | FFTW planner flags; `FFTW.ESTIMATE` for a single call |

### Return value layout

`polyft_nfft_forward` returns a `Complex{T}` array of shape `(nx÷2+1, nx)` in
the standard rfft layout:

- rows: `rfftfreq(nx, 1/pixsize)`, from 0 to the Nyquist frequency `1/(2*pixsize)`
  cycles/mas;
- columns: `fftfreq(nx, 1/pixsize)`, in FFT order (zero frequency at index 1).

`irfft(F, nx)` gives the real-space image; `fftshift` centres it:

```julia
img = fftshift(irfft(F, nx))
```

### Sign convention

`F[k] = Σ f_j exp(-2πi k·r_j)`. NFFT.jl's adjoint has the opposite sign, so the sample
positions are negated (`pos = -r/L`, `L = nx*pixsize`); `ifftshift` then takes NFFT order
(`k = -N/2, …, N/2-1`) to FFT order, and the first `nx÷2+1` rows form the rfft layout.

### `build_gauss_samples`

Maps `[-1,1]^2` onto each tessel bilinearly, splits it into `nsub × nsub` sub-squares, and puts
an `ngauss × ngauss` Gauss–Legendre rule in each. Returns `(xs, ys, fs)`: the positions (mas)
of the points and their weights (intensity × Jacobian × Gauss weight).

Total samples: `(ngauss * nsub)^2 * Npix`.

## Sky images for OITOOLS

FITS image cubes from which OITOOLS computes visibilities; see
[Exporting sky images for OITOOLS](@ref).

```@docs
sky_layout
sky_grid
render_sky_cube
write_sky_cube
read_sky_cube
check_sky_cube
```
