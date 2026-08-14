# Radial Regularizers — RADFLAT and RADIALVAR

Two regularizers that act on the *projected* geometry of the visible disk rather than on
neighbouring tessels. Both are built on the same radial binning and both operate on the
surface map **without** limb darkening.

!!! note "Attribution"
    **RADFLAT is an idea of Prof. John Monnier (University of Michigan)**, communicated
    privately (2026-07-01) and implemented in his SURFING reconstruction code. The
    implementation here is ROTIR's own, on a HEALPix tessellation rather than an image
    grid, but the regularizer and the reasoning behind it are his.

    RADIALVAR is the tessellated counterpart of the `"radialvar"` regularizer in
    [OITOOLS](https://github.com/fabienbaron/OITOOLS.jl) (`radial_variance`).

## The problem RADFLAT solves

On a **non-rotating star observed at a single epoch**, the surface regularizer is weak
enough that the reconstruction can park dark or bright patches on the limb. Those trade off
almost exactly against the limb-darkening coefficient: χ² barely moves while the LD
coefficient wanders, because a limb-darkened disk and a uniform disk with a dark rim are
nearly the same object to an interferometer.

Several epochs of a *rotating* star break this by themselves — spots move, the limb is
resampled — so RADFLAT is deliberately wrong there and would suppress real structure. It is
for the single-epoch, non-rotating case it was designed for.

## Definitions

Both use [`radflat_bins`](@ref), which sorts the visible tessels into `nbins` bins of
projected radius ρ (6 by default, matching SURFING).

**RADFLAT** — force the azimuthally averaged radial profile flat:

```
f = scale · Σ_b w_b · (I_b/I_mean − 1)²,      w_b = ρ_b² · (flux weight in bin b)
```

with `I_b` the mean brightness in bin `b` and `I_mean` the mean over the whole visible disk.

**RADIALVAR** — remove azimuthal scatter *within* each annulus, i.e. push the disk toward
circular symmetry:

```
f = scale · Σ_b var_b / I_mean²,      var_b = (1/(n_b−1)) Σ_{j∈b} (x_j − x̄_b)²
```

### They are two halves of one decomposition

Split the variance of the visible disk over the radial bins and the two regularizers are
exactly the between-group and within-group terms of a one-way ANOVA:

```
Σ_j (x_j − x̄)²  =  Σ_b n_b (x̄_b − x̄)²   +   Σ_b Σ_{j∈b} (x_j − x̄_b)²
                    └──── RADFLAT ────┘       └──── RADIALVAR ────┘
                     between annuli               within annuli
```

This is an identity, verified to machine precision in `test/test_radial_regularizers.jl`.
Its practical consequence: **each regularizer alone constrains only half the map, and using
one alone pushes structure into the other half.** On α Cen A, RADFLAT drives the radial
profile rms to 1e-4 while raising the azimuthal scatter 17×; RADIALVAR does the mirror
image. Driving *both* to zero drives the map to a constant.

## Usage

```julia
bins = radflat_bins(star; nbins = 6)          # bins follow the PROJECTED geometry
regs = Any[["sobel2", 10.0, sobel_gradient_healpix(n), 1:length(x0)],
           ["radflat",   1.0, bins, bins.idx]]   # 4th element MUST be bins.idx
x = image_reconstruct_oi(x0, data, stars; maxiter = 20000, regularizers = regs)
```

The regularizer's pixel subset (4th element) must be the same index set `radflat_bins` was
built from, or the call errors rather than silently misaligning.

```@docs
radflat_bins
spheroid_radflat_fg
spheroid_radialvar_fg
```

## Choosing the weight

The weight matters more than the mesh resolution, and it does **not** transfer between
stars. A ladder on α Cen A at nside 3 (target: the parametric ld1 = 0.1399) gives

| weight | ld1* | a/σ_a | profrms |
|---|---|---|---|
| 1 | 0.1852 ± 0.0077 | 18 | 1.4e-2 |
| 10 | 0.1488 ± 0.0030 | 58 | 4e-4 |
| **100** | **0.1478 ± 0.0021** | **81** | **1e-4** |
| 1000 | 0.1621 ± 0.0067 | 24 | 2e-4 |
| 10⁴ | 0.1645 ± 0.0077 | 21 | 3e-5 |

Too weak and residual map structure biases the answer; too strong and the quadratic penalty
stiffens the problem until VMLMB stops converging — note the significance *collapsing* past
10³ while the structure metrics keep improving, which is an optimizer failing rather than a
prior succeeding.

That ladder is well behaved because α Cen has no surface structure. **On a star that does,
the same sweep has no plateau at all.** RW Cep at nside 3 with TV2 off:

| RADFLAT weight | a/σ_a | ld1* | profrms |
|---|---|---|---|
| 0 | 1.4 | unconstrained | 0.485 |
| 1e-3 | −2.1 | unconstrained | 0.267 |
| 1e-2 | 2.2 | 1.67 ± 0.17 | 0.140 |
| 1e-1 | −0.5 | unconstrained | 0.447 |
| 1 | 0.6 | unconstrained | 0.056 |
| 100 | 3.9 | 1.73 ± 0.09 | 7e-4 |

The significance oscillates in sign with no systematic weight dependence, and the same
configuration at a different TV2 weight moves it again (weight 100 gives 3.9 at TV2 = 0,
2.1 at 1e-7 with a *negative* vertex, 1.4 at 1e-4; weight 1 gives 0.6 at TV2 = 0 but 6.2 at
TV2 = 1e-7). The point-to-point scatter along each row, 0.03–0.08 in χ²ᵣ, is comparable to
the curvature being fitted, so values of a/σ_a in the 2–4 range here are within the noise of
the estimator rather than detections. Do not tune the weight against a/σ_a.

## What these regularizers will and will not do

Measured on two benchmarks, `demos/alphacen_ld.jl` (α Cen A/B, a star with no resolved
surface structure and a published limb darkening) and `demos/rwcep_radflat.jl` (RW Cep, a
star with real asymmetry).

**RADFLAT costs nothing in fit quality.** On RW Cep with TV2 off it takes the radial profile
rms from 0.49 to 7e-4 while leaving every observable alone: χ²ᵥ₂ 2.00 → 2.07, χ²ₜ₃ₐ
1.25 → 1.28, χ²ₜ₃ₚ 1.5 → 1.5. This is reproduced at every weight and every TV2 setting
tried, and it is the claim RADFLAT was designed to make good on.

**What RADFLAT recovers depends on the smoothing prior it is paired with — and `"tv2"` is
the wrong partner.** Repeating the RW Cep weight sweep against both smoothing priors, with
the parametric Hestroffer fit (V² + T3amp, no map freedom at all) giving ld1 = 1.58:

| RADFLAT weight | with `"tv2"` @1e−4 | with `"sobel"` @1e1 |
|---|---|---|
| 1e−1 | — | 1.21 ± 0.23 (2.6σ), χ²ᵥ₂ 1.90 |
| 1e0 | 0.76 ± 0.38 (2.7σ), χ²ᵥ₂ 1.99 | 1.32 ± 0.11 (5.0σ), χ²ᵥ₂ 1.95 |
| 1e1 | — | 1.55 ± 0.11 (2.9σ), χ²ᵥ₂ 1.95 |

Against `"sobel"` the recovered coefficient climbs toward the parametric value as the weight
rises (1.21 → 1.32 → 1.55 vs 1.58), with the radial profile rms driven to 7e-4 — RADFLAT
behaving exactly as designed. Against `"tv2"` it lands at 0.76, inconsistent with the
parametric answer. `"tv2"` penalises curvature with a `k⁴` response and is neither scale- nor
resolution-normalised (see [Image reconstruction](../guides/reconstruction.md)); it is
evidently distorting the radial profile that RADFLAT then measures.

Read the `"sobel"` column carefully, though: driving the profile flat is *how* RADFLAT works,
and in the flat limit the model no longer has any radial freedom, so χ²(ld1) necessarily
approaches the parametric curve. Recovering 1.55 against a parametric 1.58 is therefore close
to tautological. What RADFLAT adds over a plain parametric fit is the **error bar** — ±0.11
here — which accounts for the azimuthal structure the map still carries (azimuthal rms 3.5)
and which a parametric fit assumes away.

**orthoLD did not survive the same test.** With `"tv2"` it produced an apparently sharp
minimum (ld1 = 1.95 ± 0.20 at 1e−3 with χ²ᵥ₂ 2.04, tightening to 2.05 ± 0.01 at 36.6σ by
1e−1, at the cost of χ²ᵥ₂ 5.93). Against `"sobel"` it produces **no constraint at any weight
from 1e−3 to 1e1**, while leaving χ² untouched (1.84–1.85). A minimum that vanishes when the
smoothing prior changes was a property of the prior pair, not of the data. On α Cen, where
the answer is known, orthoLD's 43.5σ minimum was biased +41 % against the published
coefficient — consistent with the same warning.

**RADIALVAR destroys real structure — do not use it on a star that has any.** On RW Cep it
sends χ²ₜ₃ₚ from 1.4 to 8.0, and χ²ᵥ₂ from 1.8 to 8.5, reproduced at three TV2 weights.
Closure phases are the signature of genuine asymmetry, and stellar surface features are
*azimuthal* — spots at assorted longitudes — which is precisely what RADIALVAR penalises.
This is why Monnier's regularizer is the radial one: the radial profile is dominated by limb
darkening, a model degeneracy, rather than by astrophysics.

**Do not read a limb-darkening coefficient off a reconstruction.** On α Cen, where the
answer is known, neither regularizer alone constrains ld1 (|a/σ_a| ≲ 1); both together
recover the published coefficient to 3.3 % (A) and 0.3 % (B) — but only because zeroing
both halves of the variance drives the map to a constant, at which point the model *is* the
parametric limb-darkened disk and χ²(ld1) is the parametric curve. The tell is that χ²ᵣ at
the minimum sits at the parametric floor. It is a consistency check on the forward model,
not a measurement. Fit limb darkening parametrically ([`fit_sphere_ld`](@ref)); it is far
better determined that way.

!!! warning "Convergence is part of the measurement"
    Every number above is a *difference* of χ² between reconstructions, so each one must be
    at its own minimum. On RW Cep at nside 3 the drift in χ²/n from 300 to 3000 iterations
    is 0.14–0.66 — larger than the entire spread being measured across a limb-darkening
    scan. `image_reconstruct_oi` uses a tight `gtol`, so runs stop on `maxiter`: set it from
    a convergence ladder, not from what is quick.

!!! warning "Under-populated bins"
    Equal-ρ bins on a projected disk hold ∝ ρ dρ tessels, so the *innermost* bin starves
    first. `radflat_bins` warns when a bin is empty or holds fewer than 3 tessels. nside ≥ 3
    is stable to about 1 %.

## Effective map dimension

A caveat that limits what any of this can measure. The weight a tessel carries is
`vis_weights × ldmap`, and limb darkening suppresses the outer ones, so the number of map
pixels the data actually constrain is well below the number of visible pixels — and it
*depends on the limb-darkening coefficient being fitted*:

| ld1 | effective pixels (Σw)²/Σw² | of 428 visible |
|---|---|---|
| 0.4 | 351 | 82 % |
| 1.6 | 239 | 56 % |
| 3.0 | 168 | 39 % |

At `tessellation_healpix(3)` (nside 8, 768 tessels) roughly half the visible map is
therefore unconstrained by the data, and the map curvature confirms it: ~185 of 428
directions have no measurable curvature at all — a real rank deficiency, not a numerical
one, unchanged when the whole calculation is redone in Float64.

Two practical consequences. Any uncertainty quoted for a coefficient whose value changes
the map's effective dimension should be treated with suspicion. And the natural remedy is a
map coarse enough that the data constrain all of it — which is also what makes a joint
fit of map, radius and limb darkening tractable, as SURFING does.
