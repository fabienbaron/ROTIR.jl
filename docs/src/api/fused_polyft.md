# Fused polygon Fourier transform

Matrix-free polygon Fourier transform: the visibilities are summed directly instead of
through the dense (Nuv × Npix) polyft matrix. Memory is O(Nuv + Npix) instead of
O(Nuv × Npix).

See [Chi-squared & imaging](@ref) for when to use the matrix path and when the fused path.

## Forward pass

| Function | Description |
|----------|-------------|
| `compute_polyflux_and_cvis!(F, polyflux, kx, ky, k2_inv_im, proj_west, proj_north, xw)` | complex visibilities `F[k]` and projected pixel areas `polyflux[p]` |
| `precompute_k2_inv_im(kx, ky)` | `-im / (2π(kx² + ky²))` for each uv point (0 at the origin) |
| `fused_cvis(x, star, data; intensity_model=:linear, band=nothing)` | normalised complex visibilities, the matrix-free counterpart of `poly_to_cvis` |
| `fused_cvis_parts(x, star, data; intensity_model, band, extra_weights)` | `(F, flux)`, the unnormalised transform and the total flux, used for binaries |

For each visible pixel p and uv point k, each quad edge contributes

```
F[k] += k2_inv_im[k] * sinc(k.d) * exp(-i*pi*k.c) * k_perp * xw[p]
```

with
- `d = (x2-x1, y2-y1)` the edge vector,
- `c = (x2+x1, y2+y1)` twice the edge midpoint,
- `k_perp = ky*dx - kx*dy`.

The pixel area is the shoelace sum `polyflux[p] = 0.5 * sum(x[j]*y[j+1] - x[j+1]*y[j])`.

## Adjoint passes

| Function | Description |
|----------|-------------|
| `compute_adjoint_cvis!(grad_xw, adj, kx, ky, k2_inv_im, proj_west, proj_north, polyflux)` | gradient with respect to the weighted pixel values |
| `compute_adjoint_vertices!(grad_proj_west, grad_proj_north, adj, kx, ky, k2_inv_im, proj_west, proj_north, xw)` | gradient with respect to the vertex positions |

The vertex adjoint includes the derivatives of the sinc, phase and `k_perp` terms, and of
the shoelace area through the flux normalisation.

## χ² and gradient

`fused_spheroid_chi2_fg(x, g, star, data; verbose=true)` returns χ² and fills `g`, without a
polyft matrix. It takes the same arguments as `spheroid_chi2_fg`. Steps:

1. forward pass: `F[k]` and `polyflux`;
2. V², T3amp, T3φ and χ²;
3. adjoint signal from the residuals;
4. adjoint pass to the weighted pixel values;
5. flux normalisation and the soft-visibility × limb-darkening weights.

## Kernels

`POLYFT_BACKEND[]` selects the forward kernel for the process; `with_polyft_backend(f, kind)`
sets it for one task and the tasks it starts; `polyft_backend()` returns the kernel in force.

| backend | method | HEALPix 3 | 4 | 5 |
|---------|--------|-----------|---|---|
| `:scalar` | exact closed form, plain Julia | 119 ms | 425 ms | 1472 ms |
| `:turbo` | exact closed form, SIMD (needs LoopVectorization) | 6.9 ms | 21.8 ms | 86 ms |
| `:nufft` | Gauss quadrature + FINUFFT type-3 NUFFT (needs FINUFFT) | 2.6 ms | 3.8 ms | 4.1 ms |
| `:t3` (default) | Gauss quadrature + ROTIR's type-3 NUFFT | 0.20 ms | 0.28 ms | 0.62 ms |

`:t3` and `:nufft` apply to `fused_cvis` and `fused_cvis_parts`. The gradient paths
(`compute_polyflux_and_cvis!` and the two adjoints, used by `shape_chi2_fg!`,
`interferometric_chi2` and `fused_spheroid_chi2_fg`) run the exact kernel: `:turbo` when it
is selected and loaded, `:scalar` otherwise. `turbo_available()` reports whether
LoopVectorization is loaded.
