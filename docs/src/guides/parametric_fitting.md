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

## The three functions

| function | model | free parameters |
|---|---|---|
| [`parametric_chi2`](@ref) | any `surface_type` | none — evaluates χ² at fixed parameters |
| [`fit_sphere_ld`](@ref) | round star, uniform surface | `radius`, `ld1` (+`ld2`) |
| [`fit_ellipsoid_ld`](@ref) | oblate spheroid, uniform surface | `req`, flattening, `inc`, `PA`, `ld1` |

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

## Fitting an oblate star

```julia
fit = fit_ellipsoid_ld([data], tessels; ldtype = 3, tepochs = [0.0])
fit.req, fit.rpol, fit.flattening, fit.inclination, fit.position_angle, fit.ld1
fit.chi2, fit.chi2_sphere     # ALWAYS compare these two — see below
```

On RW Cep this returns `req = 1.504 mas`, `f = 0.047`, `inc = 87.9°`, `PA = 0°`,
`ld1 = 1.615`, at `χ²/n = 4.397` against the sphere's `4.820`. Both parameters are genuinely
constrained — scanning them at the optimum gives

| PA | 0° | 40° | 80° | 120° | 160° | 180° |
|---|---|---|---|---|---|---|
| χ²/n | **4.397** | 4.954 | 5.501 | 5.879 | 4.431 | **4.397** |

| f | 0.000 | 0.020 | 0.047 | 0.080 | 0.120 |
|---|---|---|---|---|---|
| χ²/n | 5.128 | 4.695 | **4.397** | 5.014 | 8.701 |

(PA = 0° and 180° agree exactly, as they must.)

!!! warning "A better χ² is not evidence of oblateness"
    RW Cep is a slow-rotating supergiant with no reason to be measurably oblate, yet the
    ellipsoid improves χ² by Δχ² ≈ 1200. That improvement is real but its *interpretation*
    is not: at `χ²/n = 4.4` the model is still badly wrong, and the three extra parameters
    are free to absorb model error. With T3φ excluded the fit sees only V² and T3amp, and an
    asymmetric brightness distribution readily mimics an elongated disc in those. Treat a
    nonzero flattening here as "the spots are not being modelled" rather than "the star is
    oblate" — and let the imaging step reconstruct them instead.

    Parameter significance computed from Δχ² is meaningless whenever `χ²/n ≫ 1`, because the
    error bars no longer describe the residuals.

The geometry sets `radius_x = radius_y = req` and `radius_z = req·(1−f)`, i.e. rotational
symmetry about the polar axis. `surface_type = 1` permits three independent semi-axes, but
fitting a triaxial ellipsoid to interferometry is rarely justified.

!!! note "Uniform surface on purpose"
    `fit_ellipsoid_ld` holds the surface uniform rather than calling
    `temperature_map_vonZeipel_ellipsoid`. That function is a placeholder — `g ∝ 1/r²` with
    an unverified `rpole = radius_x` and a `# to check` comment in the source — so folding it
    into a shape fit would attribute its errors to the shape. For a physically
    gravity-darkened oblate star use `surface_type = 2` and `fit_parametric`.

!!! warning "Degeneracies, and the convergence guard"
    Pole-on (`inc → 0`) the projection is circular whatever the flattening, so `f`, `inc` and
    `PA` become jointly unidentifiable. Started cold, Nelder-Mead walks straight into that
    corner: an early version of this fitter returned `f = 0.354` with `inc` and `PA` pinned
    at their bounds and a χ² *worse* than the sphere.

    Three things prevent it now. `req0`/`ld0` are seeded from [`fit_sphere_ld`](@ref) by
    default (`seed_from_sphere = true`); the bounds are `inc ∈ [0°, 90°]` and
    `PA ∈ [0°, 180°)`, removing the exact two-fold degeneracies `i ↔ 180° − i` and
    `PA ↔ PA + 180°`; and because the sphere is the `f = 0` special case, a converged
    ellipsoid can never fit worse — if it does, the function **warns** and returns
    `chi2_sphere` so you can check rather than receive a plausible-looking parameter set.

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

**You want to compare models.** `fit_ellipsoid_ld` nests `fit_sphere_ld` at `f = 0`, so Δχ²
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
