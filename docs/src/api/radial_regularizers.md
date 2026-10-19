# Radial Regularizers — RADFLAT and RADIALVAR

Two regularisers that act on the projected geometry of the visible disk. Both use the same
radial binning and act on the surface map without limb darkening.

!!! note "Attribution"
    RADFLAT is an idea of Prof. John Monnier (University of Michigan), communicated privately
    (2026-07-01) and implemented in his SURFING reconstruction code. ROTIR's implementation
    works on a HEALPix tessellation rather than an image grid.

    RADIALVAR is the tessellated counterpart of the `"radialvar"` regulariser in
    [OITOOLS](https://github.com/fabienbaron/OITOOLS.jl) (`radial_variance`).

## Purpose

On a non-rotating star observed at one epoch, a reconstruction can place dark or bright
patches on the limb that trade off against the limb-darkening coefficient: χ² barely changes
while the coefficient wanders. RADFLAT removes that freedom. Several epochs of a rotating star
resample the limb and break the degeneracy without it; there RADFLAT suppresses real
structure.

## Definitions

[`radflat_bins`](@ref) sorts the visible tessels into `nbins` bins of projected radius ρ
(default 6, as in SURFING).

RADFLAT flattens the azimuthally averaged radial profile:

```
f = scale · Σ_b w_b · (I_b/I_mean − 1)²,      w_b = ρ_b² · (projected area in bin b)
```

with `I_b` the mean brightness in bin `b` and `I_mean` the mean over the visible disk.

RADIALVAR removes the azimuthal scatter within each bin:

```
f = scale · Σ_b var_b / I_mean²,      var_b = (1/(n_b−1)) Σ_{j∈b} (x_j − x̄_b)²
```

The two are the between-bin and within-bin terms of the variance of the visible disk (one-way
ANOVA):

```
Σ_j (x_j − x̄)²  =  Σ_b n_b (x̄_b − x̄)²   +   Σ_b Σ_{j∈b} (x_j − x̄_b)²
                    └──── RADFLAT ────┘       └──── RADIALVAR ────┘
                     between annuli               within annuli
```

`test/test_radial_regularizers.jl` checks this identity to machine precision. Each one alone
leaves the other term free: on α Cen A, RADFLAT drives the radial profile rms to 1e-4 and
raises the azimuthal scatter 17×; RADIALVAR does the reverse. Both together drive the map to
a constant.

## Usage

```julia
bins = radflat_bins(star; nbins = 6)          # from the projected geometry
regs = Any[["sobel2", 10.0, sobel_gradient_healpix(n), 1:length(x0)],
           ["radflat",   1.0, bins, bins.idx]]   # 4th element: bins.idx
x = image_reconstruct_oi(x0, data, stars; maxiter = 20000, regularizers = regs)
```

The pixel subset (4th element) must be the index set `radflat_bins` was built from; a length
mismatch raises an error.

```@docs
radflat_bins
spheroid_radflat_fg
spheroid_radialvar_fg
```

## Choosing the weight

The weight does not transfer between stars. On α Cen A at `tessellation_healpix(3)` (parametric ld1 = 0.1399):

| weight | ld1* | a/σ_a | profrms |
|---|---|---|---|
| 1 | 0.1852 ± 0.0077 | 18 | 1.4e-2 |
| 10 | 0.1488 ± 0.0030 | 58 | 4e-4 |
| **100** | **0.1478 ± 0.0021** | **81** | **1e-4** |
| 1000 | 0.1621 ± 0.0067 | 24 | 2e-4 |
| 10⁴ | 0.1645 ± 0.0077 | 21 | 3e-5 |

Too weak, and residual map structure biases ld1. Too strong, and VMLMB stops converging: past
10³ the significance falls while the profile keeps flattening.

α Cen has no surface structure. On RW Cep at `tessellation_healpix(3)`, with `"tv2"` off, the sweep has no
plateau:

| RADFLAT weight | a/σ_a | ld1* | profrms |
|---|---|---|---|
| 0 | 1.4 | unconstrained | 0.485 |
| 1e-3 | −2.1 | unconstrained | 0.267 |
| 1e-2 | 2.2 | 1.67 ± 0.17 | 0.140 |
| 1e-1 | −0.5 | unconstrained | 0.447 |
| 1 | 0.6 | unconstrained | 0.056 |
| 100 | 3.9 | 1.73 ± 0.09 | 7e-4 |

The significance changes sign with no trend in the weight, and changes again with the
`"tv2"` weight (weight 100: 3.9 at `"tv2"` = 0, 2.1 at 1e-7 with a negative vertex, 1.4 at
1e-4; weight 1: 0.6 at 0, 6.2 at 1e-7). The scatter along each row, 0.03–0.08 in χ²ᵣ, is
comparable to the curvature being fitted, so a/σ_a of 2–4 here is noise. Do not tune the
weight against a/σ_a.

## Measured behaviour

Benchmarks: `demos/alphacen_ld.jl` (α Cen A/B, no resolved surface structure, published limb
darkening) and `demos/rwcep_radflat.jl` (RW Cep, asymmetric surface).

**RADFLAT leaves the fit unchanged.** On RW Cep with `"tv2"` off it takes the radial profile
rms from 0.49 to 7e-4 while χ²ᵥ₂ goes 2.00 → 2.07, χ²ₜ₃ₐ 1.25 → 1.28, χ²ₜ₃ₚ 1.5 → 1.5, at
every weight and `"tv2"` setting tried.

**The recovered ld1 depends on the smoothing prior.** RW Cep, parametric Hestroffer fit
(V² + T3amp) ld1 = 1.58:

| RADFLAT weight | with `"tv2"` @1e−4 | with `"sobel"` @1e1 |
|---|---|---|
| 1e−1 | — | 1.21 ± 0.23 (2.6σ), χ²ᵥ₂ 1.90 |
| 1e0 | 0.76 ± 0.38 (2.7σ), χ²ᵥ₂ 1.99 | 1.32 ± 0.11 (5.0σ), χ²ᵥ₂ 1.95 |
| 1e1 | — | 1.55 ± 0.11 (2.9σ), χ²ᵥ₂ 1.95 |

With `"sobel"` the coefficient approaches the parametric value as the weight rises; with
`"tv2"` (`k⁴` response, not scale- or resolution-normalised; see
[Image reconstruction](../guides/reconstruction.md)) it does not. A flat profile leaves the
model no radial freedom, so χ²(ld1) tends to the parametric curve and 1.55 against 1.58 is
expected. What the reconstruction adds is the error bar (±0.11), which includes the
azimuthal structure the map still carries (azimuthal rms 3.5).

**orthoLD gave no reliable constraint.** With `"tv2"` it gave ld1 = 1.95 ± 0.20 at 1e−3
(χ²ᵥ₂ 2.04), tightening to 2.05 ± 0.01 (36.6σ) at 1e−1 with χ²ᵥ₂ 5.93. With `"sobel"` it
gave no constraint at any weight from 1e−3 to 1e1, with χ² unchanged (1.84–1.85). On α Cen its
43.5σ minimum was 41 % above the published coefficient.

**RADIALVAR removes real structure.** On RW Cep it raises χ²ₜ₃ₚ from 1.4 to 8.0 and χ²ᵥ₂ from
1.8 to 8.5, at three `"tv2"` weights. Surface features at assorted longitudes are azimuthal
structure, which is what RADIALVAR penalises. Do not use it on a star with surface structure.

**A reconstruction does not measure limb darkening.** On α Cen, neither regulariser alone
constrains ld1 (|a/σ_a| ≲ 1). Both together recover the published coefficient to 3.3 % (A) and
0.3 % (B), because the map is then constant and χ²(ld1) is the parametric curve; χ²ᵣ at the
minimum equals the parametric value. Fit limb darkening parametrically
([`fit_sphere_ld`](@ref)).

!!! warning "Convergence"
    Every number above is a difference of χ² between reconstructions, each of which must be
    converged. On RW Cep at `tessellation_healpix(3)`, χ²/n changes by 0.14–0.66 between 300 and 3000
    iterations, more than the spread across a limb-darkening scan. `image_reconstruct_oi`
    uses a tight `gtol`, so runs stop at `maxiter`; choose it from a convergence ladder.

!!! warning "Under-populated bins"
    Equal-ρ bins hold a number of tessels ∝ ρ dρ, so the innermost bin empties first.
    `radflat_bins` warns when a bin is empty or holds fewer than 3 tessels. nside ≥ 8
    (`tessellation_healpix(3)`) is stable to about 1 %.

## Effective map dimension

Each tessel is weighted by `vis_weights × ldmap`, and limb darkening suppresses the outer
ones, so the data constrain fewer pixels than are visible, and the number depends on ld1:

| ld1 | effective pixels (Σw)²/Σw² | of 428 visible |
|---|---|---|
| 0.4 | 351 | 82 % |
| 1.6 | 239 | 56 % |
| 3.0 | 168 | 39 % |

At `tessellation_healpix(3)` (nside 8, 768 tessels) about half the visible map is
unconstrained: ~185 of 428 directions have no measurable curvature, in Float32 and Float64
alike. An uncertainty on a coefficient that changes the effective dimension is therefore
unreliable. A map coarse enough for the data to constrain all of it makes a joint fit of map,
radius and limb darkening tractable, as in SURFING.
