# Parametric Fitting

Before you can reconstruct a surface map you need the star's angular size, because the
tessellation is built at a fixed `radius`. This is not a detail to postpone: on a real
dataset a 55 % radius error gave a reconstruction with `χ²/n ≈ 4000` where the correct radius
gives `≈ 4`. The map cannot absorb a wrong size — it will try, and produce nonsense.

This page covers fitting a star whose surface carries no structure beyond limb darkening.
For a gravity-darkened rapid rotator see `fit_parametric` and
[Parametric Model & AD](../api/parametric_gradient.md), which solve the Roche problem
properly; the wrappers here are deliberately simpler models for stars that do not rotate fast
enough to be oblate.

A rapid-rotator fit also has to choose a gravity-darkening law, and the choice is not free of
consequences for the limb darkening it is fitted alongside: both take flux out of the limb, so
a law that asks for too much equatorial cooling drives the limb-darkening coefficient down and
can push it negative. Set `gravity_law` on the model to pick between von Zeipel (1924) and
Espinosa Lara & Rieutord (2011, A&A 533, A43), and rank them by log-evidence when the data can
tell. See [Choosing a gravity-darkening law](surfaces.md#Choosing-a-gravity-darkening-law) for
the contrast tables and the references, and `demos/gravity_law_comparison.jl` for a worked
comparison.

## The three functions

| function | model | free parameters |
|---|---|---|
| [`parametric_chi2`](@ref) | any `surface_type` | none — evaluates χ² at fixed parameters |
| [`fit_sphere_ld`](@ref) | round star, uniform surface | `radius`, `ld1` (+`ld2`) |
| [`fit_ellipsoid_ld`](@ref) | triaxial ellipsoid, uniform surface | `radius_x`, `radius_y`, `radius_z`, `inc`, `PA`, `ld1` |

Both fitters default to NLopt's Nelder–Mead (`:LN_NELDERMEAD`) — derivative-free and
bound-constrained, which is the right trade at two to five parameters where the χ² is smooth
but each evaluation costs a geometry build and a polygon Fourier transform. Pass
`method = :ultranest` for a posterior with uncertainties and `log(Z)` instead of a point
estimate; see [Choosing the optimiser](#Choosing-the-optimiser:-`method`) below.

## Fitting a round star

```julia
using ROTIR
data    = readoifits("RW_Cep.oifits")[1,1]     # note the [1,1]: readoifits returns
tessels = tessellation_healpix(4)              # Array{OIdata}(nwavbin, ntimebin)

fit = fit_sphere_ld([data], tessels; ldtype = 3, tepochs = [0.0])

fit.radius          # POLAR radius in mas — angular DIAMETER is 2*radius
fit.ld1             # Hestroffer exponent α for ldtype = 3
fit.chi2_per_datum
fit.params          # NamedTuple ready for create_star / create_star_multiepochs
```

`fit.params` is the point of the whole exercise — hand it straight to the reconstruction.

### Limb-darkening laws

`ldtype` selects the law, and `ld1` means something different in each:

| `ldtype` | law | `ld1` | range |
|---|---|---|---|
| 1 | linear, `1 − u(1−μ)` | `u` | `0 ≤ u ≤ 1` |
| 2 | quadratic, `1 − u(1−μ) − v(1−μ)²` | `u`, with `ld2 = v` | — |
| 3 | Hestroffer power law, `μ^α` | `α` | `α ≥ 0`, unbounded above |

Do not compare `ld1` across laws. On RW Cep the linear fit **pins at `u = 1.0000`**, its
physical ceiling — `u > 1` implies negative intensity at the limb — with `χ²/n = 6.53`,
while the power law reaches `α = 1.575` at `χ²/n = 4.82`. The data want more limb darkening
than a linear law can express, and the pinned coefficient is the symptom.

## Why closure phases are excluded by default

Both fitters default to `weights = [1, 1, 0]`: V² and T3amp, **not** T3φ.

Any centrosymmetric model — a uniform disc, a limb-darkened disc, a uniform ellipsoid —
produces closure phases of exactly 0° or 180°. A spotted star does not. RW Cep's T3φ has an
rms of **118°**. Fitting it with a symmetric model does not improve the fit; it drags the
model toward a compromise that is worse in the amplitudes:

| fitted on | D (mas) | α | χ²v2/n | χ²t3a/n |
|---|---|---|---|---|
| V² only | 3.00 | 1.80 | 8.6 | 6.1 |
| **V² + T3amp** (default) | **2.93** | **1.58** | 8.9 | 5.4 |
| V² + T3amp + T3φ | 2.55 | 0.89 | 24.8 | 9.5 |

Including T3φ shrinks the diameter by 13 % and the LD coefficient by 44 % while making the
amplitude fits ~3× worse. The large residual T3φ is not a failure — it **is** the asymmetry,
and it is exactly what the imaging step reconstructs. Pass `weights = [1,1,1]` only if you
have independent reason to believe the star is symmetric.

## Accuracy versus tessellation resolution

ROTIR evaluates the model on a tessellated sphere; OITOOLS has closed-form visibilities for
the same limb-darkened discs. They agree, and the agreement improves with resolution.
Fitting the same RW Cep data with both:

| | D (mas) | α | ΔD | Δα |
|---|---|---|---|---|
| OITOOLS analytic (`visibility_ldpow`) | 2.9309 | 1.5799 | — | — |
| ROTIR nside 3 | 2.9271 | 1.5604 | −0.13 % | −1.2 % |
| ROTIR nside 4 | 2.9299 | 1.5749 | −0.03 % | −0.3 % |
| ROTIR nside 5 | 2.9308 | 1.5792 | −0.003 % | −0.04 % |

The discretisation-induced χ² is negligible at every resolution (rms `ΔV² = 3.3e-4` at
nside 3 against typical error bars of `2e-3`). nside 3 is adequate for a size estimate;
use 4–5 if the limb-darkening coefficient itself is the measurement.

## Fitting a non-spherical star

```julia
fit = fit_ellipsoid_ld(data, tessels; tepochs = ts, rotation_period = 17.769426, ldtype = 3)
fit.radius_x, fit.radius_y, fit.radius_z, fit.inclination, fit.position_angle, fit.ld1
fit.equatorial_elongation, fit.polar_flattening      # the two shape numbers worth reading
fit.chi2, fit.chi2_sphere                            # ALWAYS compare these two — see below
```

The free parameters are the model's own three semi-axes. This matters more than it sounds:
until 2026-10-03 this function imposed `radius_x = radius_y` and fitted an **oblate** figure,
which cannot represent a tidally distorted star — elongated *towards its companion*, so
`radius_x > radius_y ≈ radius_z` — at any parameter value. Asked for one it drives the shape
to its bound.

The two reported shape numbers separate the two physical causes: `equatorial_elongation`
= `radius_x/radius_y` is tidal, `polar_flattening` = `1 − 2·radius_z/(radius_x+radius_y)` is
rotational.

**Worked example, ζ And** (six MIRC-X epochs over 1.13 rotations, `P = 17.769426 d`):

| | nside 3 | nside 2 |
|---|---|---|
| `radius_x`, `radius_y`, `radius_z` [mas] | 1.3275, 1.2663, 1.2120 | 1.3364, 1.2728, 1.2237 |
| equatorial elongation | 1.048 | 1.050 |
| polar flattening | 0.066 | 0.062 |
| equatorial/polar | **1.070** | **1.066** |
| `inc`, `PA` [deg] | 64.7, 122.8 | 68.0, 123.7 |
| χ²/n (sphere) | 3.63 (9.13) | 3.60 (9.26) |

against the published equatorial-to-polar ratio of **1.06**. The 0.4 % difference between the
two meshes is the discretisation error on this quantity.

!!! warning "Restart, or you will get a different answer"
    This objective has several basins and the wrong ones look converged. Twelve starts on ζ And
    at nside 3 gave three: the physical one above at `χ²/n = 3.626`, a **pole-on** corner at
    4.325 with `radius_z` pinned to its ceiling and `inc = 9°`, and an edge-on one at 5.561.
    Four of the twelve — including the default start — landed in the pole-on corner and
    reported an equatorial/polar ratio of 0.41 where the answer is 1.07. Hence `restarts = 6`
    by default, and hence the warnings the function emits when a parameter sits at a bound or
    the inclination falls below 15°.

!!! warning "A better χ² is not evidence of a shape"
    The ellipsoid improves ζ And's χ²/n from 9.13 to 3.63, but at `χ²/n = 3.6` the model is
    still wrong and three extra parameters are free to absorb model error. With T3φ excluded
    the fit sees only V² and T3amp, where an asymmetric brightness distribution readily mimics
    an elongated disc. The size of that effect is measurable: synthetic data from a *known*
    figure carrying three spots at the 20 % contrast of ζ And's published map returns the
    equatorial/polar ratio **1 % too large** (1.070 against a truth of 1.060). So a shape at
    the few-percent level is believable here and one at the tens-of-percent level is the spots
    talking — and either way the imaging step, not a bound, is what should absorb them.

    Parameter significance computed from Δχ² is meaningless whenever `χ²/n ≫ 1`, because the
    error bars no longer describe the residuals.

!!! note "The rotation phase is now load-bearing"
    An oblate figure is axisymmetric, so the epoch times did not affect it. A triaxial one
    presents a different projected ellipse at every rotation phase — which is exactly what
    makes the three radii separately identifiable — so `rotation_period` is **required** for
    multi-epoch data and has no default, and `t0` sets the phase zero. From a single phase the
    projected outline is one ellipse, three numbers against five shape and orientation
    parameters; the fit runs and warns, and the radii are one point on a degenerate surface.

!!! note "Uniform surface on purpose"
    `fit_ellipsoid_ld` holds the surface uniform rather than calling
    `temperature_map_vonZeipel_ellipsoid`. That function is a placeholder — `g ∝ 1/r²` with
    an unverified `rpole = radius_x` and a `# to check` comment in the source — so folding it
    into a shape fit would attribute its errors to the shape. For a physically
    gravity-darkened oblate star use `surface_type = 2` and `fit_parametric`.

!!! warning "Degeneracies, and the convergence guard"
    Pole-on (`inc → 0`) the projection is circular whatever the shape, so the radii, `inc` and
    `PA` become jointly unidentifiable. Started cold, Nelder-Mead walks straight into that
    corner: an early version of this fitter returned a flattening of 0.354 with `inc` and `PA`
    pinned at their bounds and a χ² *worse* than the sphere.

    Four things prevent it now. `r0`/`ld0` are seeded from [`fit_sphere_ld`](@ref) by default
    (`seed_from_sphere = true`); `restarts` sweeps the starting orientation and keeps the best
    χ²; the bounds are `inc ∈ [0°, 90°]` and `PA ∈ [0°, 180°)`, removing the two-fold
    degeneracies `i ↔ 180° − i` and `PA ↔ PA + 180°` (all but exact for a triaxial figure,
    since flipping the inclination reverses the direction the long axis sweeps with phase); and
    because the sphere is the equal-radii special case, a converged ellipsoid can never fit
    worse — if it does, the function **warns** and returns `chi2_sphere` so you can check
    rather than receive a plausible-looking parameter set. It also warns when any parameter
    lands within 0.1 % of a bound, or when `inc` falls below 15°.

## Choosing the optimiser: `method`

Both fitters take `method`, selecting how the χ² surface is explored:

```julia
fit_sphere_ld(data, tessels; method = :neldermead)   # default — fast point estimate
fit_sphere_ld(data, tessels; method = :nautilus)     # posterior, uncertainties, log(Z)
fit_sphere_ld(data, tessels; method = :ultranest)    # the same, through Python
fit_sphere_ld(data, tessels; method = :pigeons)      # tempering — for a MULTIMODAL posterior
```

!!! tip "Which sampler, and when"
    `:nautilus` is the default choice for a posterior with one mode: pure Julia, importance
    nested sampling, and it reuses every likelihood evaluation. `:ultranest` answers the same
    question through Python and is kept for continuity with published runs.

    `:pigeons` is for the case the other two cannot see. A ROTIR χ² is routinely
    **multimodal** — `demos/rho_cas_basins.jl` finds several minima between 2.2 and 3.7 mas
    diameter on one star — and a single NUTS chain samples whichever basin it started in and
    reports a tight, confident interval that says nothing about the others. Non-reversible
    parallel tempering moves between basins along the annealed chain ladder and returns
    `log(Z)` by stepping stone, so nothing is given up for the coverage.

    Read the **round-trip count**, not the sample count: a run whose chains never traversed
    the ladder has not visited the other modes, and its posterior is the single-basin answer a
    cheaper sampler would have given. `_fit_pigeons` warns when it is below three.

!!! note "`:pigeons` needs `using Pigeons, Distributions, LogDensityProblems, ADTypes, Zygote`"
    All four besides Pigeons are Pigeons' own dependencies, so nothing extra is installed.
    Like PythonCall, Pigeons is a **weak** dependency and the GUI launcher does not load it:
    measured, `using Pigeons` invalidates the plot-construction code OITOOLS precompiles for
    its live canvas and takes one canvas build from 370 ms to 3166 ms, i.e. 2.8 s onto every
    GUI start. Load it in a script, or before `gui()` when you want that fit method.

!!! note "`:ultranest` needs `using PythonCall`"
    UltraNest is a Python package, and ROTIR reaches it through
    `ext/ROTIRUltraNestExt.jl`, which PythonCall triggers. PythonCall is a **weak**
    dependency: loading it invalidates the plot-construction code OITOOLS precompiles for
    its live canvas — measured at 338 ms → 2477 ms for one canvas build — which put 1.2 s on
    every ROTIR GUI start whether or not anything sampled. Add `using PythonCall` to a script
    that wants `:ultranest`, or use `:nautilus`, which is the pure-Julia nested sampler and
    needs nothing extra.

`:neldermead` uses NLopt's simplex and returns a point. `:ultranest` runs nested sampling
over the *same* uniform box prior — deliberately the same bounds, so the two methods are
answering the same question and `:ultranest` cannot quietly explore somewhere the optimiser
could not.

Both return the same fields, so code reading `fit.params` or `fit.ld1` works either way.
`:ultranest` adds:

| field | meaning |
|---|---|
| `radius_err`, `ld1_err` | half the 68 % interval |
| `q16`, `q84` | the 68 % interval itself, per parameter |
| `logz`, `logzerr` | log evidence and its uncertainty |
| `samples` | the full posterior, `nsamples × nfree` |

### When the point estimate is not enough

Two situations make `:ultranest` worth its cost.

**You need an uncertainty, not a value.** A limb-darkening coefficient quoted without a width
says nothing about whether the data constrain it — and on a spotted star the surface map can
absorb much of the limb signal, leaving the coefficient far less determined than a single
χ² minimum suggests. If the width is the result you care about, the optimiser cannot give it
to you.

**You want to compare models.** `fit_ellipsoid_ld` nests `fit_sphere_ld` at equal radii, so Δχ²
looks like a clean test — but Δχ² only means what you think when `χ²/n ≈ 1`. On RW Cep it is
4.4, the error bars do not describe the residuals, and an apparently decisive Δχ² ≈ 1200 is
not evidence of oblateness. `log(Z)` from both fits is the honest comparison, because it
penalises the extra parameters through the prior volume rather than assuming Gaussian errors.

### Cost, and the thread requirement

Nested sampling ignores gradients and its cost grows steeply with dimension. Two parameters
need order 10⁴ likelihood evaluations; five with a step sampler can need far more. Each
evaluation here builds a star geometry and a polygon Fourier transform, so budget accordingly
— and prefer a coarse tessellation (nside 3) for the sampling run, since the size estimate
converges long before the LD coefficient does.

`use_stepsampler` defaults differ by model: `false` for the two-parameter sphere, `true` for
the five-parameter ellipsoid. UltraNest's region sampler degrades above roughly four
dimensions, and the ellipsoid additionally has two near-degenerate directions.

!!! warning "UltraNest requires a single Julia thread"
    `:ultranest` errors immediately unless Julia is running single-threaded — the Python
    interface it uses is not thread-safe, and a multi-threaded run would crash part-way
    through, losing the sampling done so far.

    ```
    JULIA_NUM_THREADS=1 julia --project=demos yourscript.jl
    ```

    `JULIA_NUM_THREADS` in your environment counts, so passing `-t 1` on the command line is
    not enough if that variable is set to something else. `:neldermead` has no such
    restriction.

!!! tip "Cross-check against `:neldermead`"
    Run both and confirm they agree to within the posterior width. Nested sampling fails
    quietly — a badly chosen prior gives a confident, tightly-constrained posterior around a
    wrong answer, which is indistinguishable from a good fit unless you have something to
    compare against.

    The two do not answer identical questions: `:neldermead` returns the χ² minimum,
    `:ultranest` the posterior median, and these differ for a skewed posterior. Agreement to
    within a fraction of 1σ is the criterion, not exact equality.

## Evaluating χ² without fitting

```julia
parametric_chi2(params, tessels, data, tepochs; weights = [1,1,0], uniform = true)
```

Useful for scanning a parameter by hand, or for checking a published model against your data.
`uniform = true` forces uniform surface brightness instead of calling
`parametric_temperature_map`; for `surface_type = 0` that is already what the map is,
so it changes nothing.

## Where to go next

With `fit.params` in hand, see [Image Reconstruction](reconstruction.md) for recovering the
surface structure that the parametric model leaves in the residuals — for RW Cep, a residual
`χ²t3φ/n` of several hundred, which is the spot signal.

```@docs
parametric_chi2
fit_sphere_ld
fit_ellipsoid_ld
```
