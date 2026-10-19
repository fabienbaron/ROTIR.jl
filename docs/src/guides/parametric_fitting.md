# Parametric Fitting

A surface reconstruction needs the star's angular size first: the tessellation is built at a
fixed `radius`, and the map cannot correct a wrong one. On one dataset a 55 % radius error gave
`χ²/n ≈ 4000` after reconstruction, against `≈ 4` at the correct radius.

This page covers stars whose surface has no structure beyond limb darkening. For a
gravity-darkened rapid rotator, use `fit_parametric` (see
[Parametric Model & AD](../api/parametric_gradient.md)). Its `gravity_law` parameter selects
von Zeipel (1924, `gravity_law = 1`) or Espinosa Lara & Rieutord (2011, A&A 533, A43;
`gravity_law = 2` or `:elr`). Gravity darkening and limb darkening both take flux out of the
limb, so a law with too much equatorial cooling pushes the limb-darkening coefficient down,
possibly below zero. See [Choosing a gravity-darkening law](surfaces.md#Choosing-a-gravity-darkening-law)
and `demos/gravity_law_comparison.jl`.

## The three functions

| function | model | free parameters |
|---|---|---|
| [`parametric_chi2`](@ref) | any `surface_type` | none; evaluates χ² |
| [`fit_sphere_ld`](@ref) | round star, uniform surface | `radius`, `ld1` (+`ld2`) |
| [`fit_ellipsoid_ld`](@ref) | triaxial ellipsoid, uniform surface | `radius_x`, `radius_y`, `radius_z`, `inc`, `PA`, `ld1` |

Both fitters default to NLopt's Nelder–Mead (`algorithm = :LN_NELDERMEAD`; `:LN_BOBYQA` and
`:LN_SBPLX` also work). `method = :nautilus` or `:ultranest` returns a posterior and `log(Z)`;
see [Choosing the optimiser](@ref).

## Fitting a round star

```julia
using ROTIR
data    = readoifits("RW_Cep.oifits")[1,1]     # readoifits returns Array{OIdata}(nwavbin, ntimebin)
tessels = tessellation_healpix(4)

fit = fit_sphere_ld([data], tessels; ldtype = 3, tepochs = [0.0])

fit.radius          # polar radius in mas; the angular diameter is 2*radius
fit.ld1             # Hestroffer exponent α for ldtype = 3
fit.chi2_per_datum
fit.params          # NamedTuple for create_star / create_star_multiepochs
```

`fit.params` is the input to the reconstruction. `tpole` (default 5000 K) only sets the flux
scale, which cancels in normalised visibilities, so it is not fitted. The starting radius
defaults to a quarter of `λ/B_max`.

### Limb-darkening laws

| `ldtype` | law | `ld1` | default bounds |
|---|---|---|---|
| 1 | linear, `1 − u(1−μ)` | `u` | `0 ≤ u ≤ 1` |
| 2 | quadratic, `1 − u(1−μ) − v(1−μ)²` | `u`; `ld2 = v` is fitted with `fit_ld2 = true` | `0 ≤ u ≤ 1`, `−1 ≤ v ≤ 1` |
| 3 | Hestroffer power law, `μ^α` | `α` | `0 ≤ α ≤ 4` |

`ld1` values are not comparable across laws. On RW Cep the linear fit stops at `u = 1.0000`,
the physical limit (`u > 1` gives negative intensity at the limb), with `χ²/n = 6.53`; the
power law reaches `α = 1.575` at `χ²/n = 4.82`. A coefficient at its bound means the law cannot
express the limb darkening the data require.

## Closure phases are excluded by default

Both fitters default to `weights = [1, 1, 0]`: V² and T3amp, not T3φ.

A centrosymmetric model (uniform disc, limb-darkened disc, uniform ellipsoid) has closure
phases of 0° or 180° only. RW Cep's T3φ has an rms of 118°. Including T3φ in a symmetric fit
biases the size and limb darkening and worsens the amplitude fits:

| fitted on | D (mas) | α | χ²v2/n | χ²t3a/n |
|---|---|---|---|---|
| V² only | 3.00 | 1.80 | 8.6 | 6.1 |
| V² + T3amp (default) | 2.93 | 1.58 | 8.9 | 5.4 |
| V² + T3amp + T3φ | 2.55 | 0.89 | 24.8 | 9.5 |

With T3φ, the diameter is 13 % smaller and the LD coefficient 44 % smaller, and the amplitude
χ² is about 3× larger. The T3φ residual is the asymmetry the imaging step reconstructs. Use
`weights = [1,1,1]` only for a star known to be symmetric.

## Accuracy versus tessellation resolution

OITOOLS has closed-form visibilities for the same limb-darkened discs. Fitting RW Cep with
both:

| | D (mas) | α | ΔD | Δα |
|---|---|---|---|---|
| OITOOLS analytic (power-law LD disc) | 2.9309 | 1.5799 | — | — |
| ROTIR nside 3 | 2.9271 | 1.5604 | −0.13 % | −1.2 % |
| ROTIR nside 4 | 2.9299 | 1.5749 | −0.03 % | −0.3 % |
| ROTIR nside 5 | 2.9308 | 1.5792 | −0.003 % | −0.04 % |

The discretisation error in V² is an rms of `3.3e-4` at nside 3, against typical error bars of
`2e-3`. nside 3 is enough for a size; use 4–5 when the limb-darkening coefficient is the
measurement.

## Fitting a non-spherical star

```julia
fit = fit_ellipsoid_ld(data, tessels; tepochs = ts, rotation_period = 17.769426, ldtype = 3)
fit.radius_x, fit.radius_y, fit.radius_z, fit.inclination, fit.position_angle, fit.ld1
fit.equatorial_elongation, fit.polar_flattening, fit.mean_radius
fit.chi2, fit.chi2_sphere                            # compare: see below
```

The free parameters are the three semi-axes of `surface_type = 1`, so the fit can represent a
tidally elongated star (`radius_x > radius_y ≈ radius_z`) as well as an oblate one.
`equatorial_elongation = radius_x/radius_y` measures tidal distortion;
`polar_flattening = 1 − 2·radius_z/(radius_x+radius_y)` measures rotational flattening.
`mean_radius` is `(radius_x·radius_y·radius_z)^(1/3)`.

**ζ And** (six MIRC-X epochs over 1.13 rotations, `P = 17.769426 d`):

| | nside 3 | nside 2 |
|---|---|---|
| `radius_x`, `radius_y`, `radius_z` [mas] | 1.3275, 1.2663, 1.2120 | 1.3364, 1.2728, 1.2237 |
| equatorial elongation | 1.048 | 1.050 |
| polar flattening | 0.066 | 0.062 |
| equatorial/polar | 1.070 | 1.066 |
| `inc`, `PA` [deg] | 64.7, 122.8 | 68.0, 123.7 |
| χ²/n (sphere) | 3.63 (9.13) | 3.60 (9.26) |

Roettenbacher et al. find an equatorial/polar ratio of 1.06 with a different model. The 0.4 %
difference between the two meshes is the discretisation error.

!!! warning "Several minima"
    Twelve starts on ζ And at nside 3 found three minima: the one above at `χ²/n = 3.626`, a
    pole-on one at 4.325 (`radius_z` at its upper bound, `inc = 9°`) and an edge-on one at
    5.561. Four of the twelve, including the default start, ended pole-on with an
    equatorial/polar ratio of 0.41. `restarts = 6` (default) runs Nelder–Mead from a grid of
    starting inclinations and position angles and keeps the best; `restarts = 1` runs once from
    `inc0` (default 60°) and `pa0` (default 45°).

!!! warning "χ² improvement and shape"
    The ellipsoid lowers ζ And's χ²/n from 9.13 to 3.63, but at `χ²/n = 3.6` the model is still
    incomplete, and V² and T3amp alone cannot separate surface brightness asymmetry from
    elongation. Synthetic data from a figure with equatorial/polar ratio 1.060 and three spots
    at ζ And's 20 % contrast return 1.072 (`demos/zetand_triaxial.jl`): spots bias the ratio by
    about 1 %. Shapes at the few-percent level are measurable here; tens of percent would be
    spot signal.

    Uncertainties from Δχ² are not valid when `χ²/n ≫ 1`.

!!! note "Rotation period and phase"
    A triaxial figure projects to a different ellipse at each rotation phase, which is what
    makes the three radii separately measurable. `rotation_period` (same units as `tepochs`)
    has no default and is required when the epochs differ; `t0` sets phase zero. For a
    synchronised binary, set `t0` to conjunction so that phase 0 faces the companion. With all
    epochs at one phase, the outline is one ellipse (three numbers for five shape and
    orientation parameters): the fit runs and warns that the radii are degenerate.

!!! note "Uniform surface"
    `fit_ellipsoid_ld` holds the surface uniform. `temperature_map_vonZeipel_ellipsoid` uses
    `g ∝ 1/r²` with `rpole = radius_x` and is not validated. For a gravity-darkened oblate star
    use `surface_type = 2` with `fit_parametric`; for a tidally distorted star in a close binary,
    `surface_type = 3` (Roche lobe).

!!! warning "Degeneracies and checks"
    Pole-on (`inc → 0`), the projection is circular for any shape, and the radii, `inc` and
    `PA` are not separately identifiable. Bounds are `inc ∈ [0°, 90°]` and `PA ∈ [0°, 180°]`:
    `180° − inc` and `PA + 180°` are almost indistinguishable (for a triaxial body, flipping the
    inclination reverses the direction the long axis sweeps with phase).

    - The radii and `ld1` start from a [`fit_sphere_ld`](@ref) fit (`seed_from_sphere = true`).
    - The sphere is the equal-radii case, so the ellipsoid χ² cannot be higher; if it is, the
      function warns and `chi2_sphere` shows the comparison.
    - It warns when a parameter is within 0.1 % of the range from a bound, and when `inc < 15°`.

## Choosing the optimiser

```julia
fit_sphere_ld(data, tessels; method = :neldermead)   # default, point estimate
fit_sphere_ld(data, tessels; method = :nautilus)     # posterior, uncertainties, log(Z); `using Nautilus`
fit_sphere_ld(data, tessels; method = :ultranest)    # the same through Python; `using PythonCall`
```

`:nautilus` is pure-Julia importance nested sampling. `:ultranest` runs UltraNest through
PythonCall. Loading PythonCall invalidates OITOOLS' precompiled plotting code: one canvas build
goes from 338 ms to 2477 ms.

Parallel tempering (Pigeons) is offered in the GUI's Model tab for sphere, ellipsoid and
rapid-rotator fits, with `using Pigeons, Distributions, LogDensityProblems, ADTypes, Zygote`.
Use it for multimodal posteriors: `demos/rho_cas_basins.jl` finds several minima between 2.2
and 3.7 mas diameter for one star. Check the round-trip count; with fewer than three, the
chains did not cross the temperature ladder and the result covers one mode only.

The samplers use the same box prior as the optimiser bounds. The result has the same fields
as the point estimate, with the posterior median as the value, plus:

| field | meaning |
|---|---|
| `radius_err`, `ld1_err` | half the 68 % interval (`fit_sphere_ld` only) |
| `q16`, `q84` | the 68 % interval, per parameter |
| `logz`, `logzerr` | log evidence and its uncertainty |
| `samples` | posterior samples, `nsamples × nfree` |

### When to sample

- **Uncertainties.** On a spotted star the surface can absorb part of the limb signal, so the
  limb-darkening coefficient may be much less constrained than the χ² minimum suggests.
- **Model comparison.** `fit_ellipsoid_ld` contains `fit_sphere_ld` as a special case, but Δχ²
  is meaningful only when `χ²/n ≈ 1`. On RW Cep `χ²/n = 4.4`, and a Δχ² of about 1200 does not
  show oblateness. Compare `log(Z)` instead.

### Cost

Nested sampling uses no gradient and its cost grows quickly with dimension: about 10⁴
likelihood evaluations for two parameters, many more for six. Each evaluation builds a
geometry and a polygon Fourier transform; use nside 3 for sampling runs.

For UltraNest, `use_stepsampler` defaults to `false` for the sphere and `true` for the ellipsoid,
and `min_num_live_points = 400`. Nautilus evaluates likelihoods in parallel (`threaded = true`).

!!! tip "Cross-check with `:neldermead`"
    The optimiser returns the χ² minimum and the samplers the posterior median; for a skewed
    posterior these differ. They should agree within a fraction of the posterior width. A
    poorly chosen prior gives a narrow posterior around a wrong value.

## Evaluating χ² without fitting

```julia
parametric_chi2(params, tessels, data, tepochs; weights = [1,1,0], uniform = true)
```

`weights` defaults to `(1, 1, 1)`. `uniform = true` uses a uniform surface instead of
`parametric_temperature_map` (identical for `surface_type = 0`).

## Next

[Image Reconstruction](reconstruction.md) recovers the surface structure left in the residuals
(on RW Cep, a residual `χ²t3φ/n` of several hundred).

```@docs
parametric_chi2
fit_sphere_ld
fit_ellipsoid_ld
```
