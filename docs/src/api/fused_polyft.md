# Fused polygon Fourier transform

Matrix-free two-pass polygon Fourier transform that eliminates the dense
(Nuv x Npix) polyft matrix. Memory scales as O(Nuv + Npix) instead of
O(Nuv * Npix).

## Forward pass

| Function | Description |
|----------|-------------|
| `compute_polyflux_and_cvis!(F, polyflux, kx, ky, k2_inv_im, projx, projy, xw)` | Compute complex visibilities `F[k]` and polygon areas `polyflux[p]` simultaneously |
| `precompute_k2_inv_im(kx, ky)` | Precompute `-im / (2*pi*(kx^2 + ky^2))` for all UV points |

### Forward pass details

For each visible pixel p and each UV point k, the contribution from each quad
edge is:

```
F[k] += k2_inv_im[k] * sinc(k.d) * exp(-i*pi*k.c) * k_perp * xw[p]
```

where:
- `d = (x2-x1, y2-y1)` is the edge vector
- `c = (x2+x1, y2+y1)` is twice the edge midpoint
- `k_perp = ky*dx - kx*dy`

The polygon area uses the shoelace formula:
`polyflux[p] = 0.5 * sum(x[j]*y[j+1] - x[j+1]*y[j])`.

## Adjoint passes

| Function | Description |
|----------|-------------|
| `compute_adjoint_cvis!(grad_xw, adj, kx, ky, k2_inv_im, projx, projy, polyflux)` | Backpropagate adjoint signal to pixel value gradients |
| `compute_adjoint_vertices!(grad_projx, grad_projy, adj, kx, ky, k2_inv_im, projx, projy, xw, polyflux)` | Backpropagate adjoint signal to vertex position gradients |

The adjoint pass recomputes the same edge contributions as the forward pass but
accumulates the gradient signal instead of visibilities. The vertex adjoint
includes derivatives of the sinc, phase, and perpendicular terms, plus the
shoelace area gradient for flux normalization.

## Drop-in chi2 function

| Function | Description |
|----------|-------------|
| `fused_spheroid_chi2_fg(x, g, star, data; verbose)` | Matrix-free chi2 + gradient, drop-in replacement for `spheroid_chi2_fg` |

This function:
1. Forward pass: computes F[k] and polyflux
2. Computes observables (V2, T3amp, T3phi) and chi2
3. Builds adjoint signal from observable residuals
4. Adjoint pass: computes gradient w.r.t. pixel temperatures
5. Applies flux normalization and soft visibility chain rule
