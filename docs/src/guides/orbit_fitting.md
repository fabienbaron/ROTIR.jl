# Orbit Fitting

Fit the relative orbit of a two-component system directly to interferometric observables,
starting from nothing but `readoifits` output and a guess.

```julia
using ROTIR

data = readoifits("betlyr_2007Jul03.oifits")[1, 1]

res = fit_orbit(data,
                UniformDisk(diameter = 0.57),                        # component 1
                EllipticalGaussian(fwhm = 0.90, ratio = 0.61, pa = 73.7);  # component 2
                elements = (P = 12.9414, e = 0.0, a = 0.865,
                            i = 93.5, Omega = 73.7, T0 = 2454283.043),
                flux_ratio = 0.81,
                method = :neldermead)

res.elements.a      # fitted semi-major axis, mas
res.chi2_red        # reduced χ²
```

Component 1 sits at the origin and component 2 is displaced by the orbit, so the model is

```
cvis = (V₁ + f·V₂·exp(-2πi(u·Δα + v·Δδ))) / (1 + f)
```

with `f` the flux of component 2 relative to component 1.

## Data

`data` is a single `OIdata` or a vector of them, one per epoch:

```julia
nights = [readoifits(f)[1, 1] for f in filter(endswith(".oifits"), readdir(dir, join=true))]
res = fit_orbit(nights, comp1, comp2; elements = (...))
```

Epochs are fitted jointly against one set of orbital elements — that is the point, since a
single night rarely constrains an orbit.

## Components

| Type | Parameters | Notes |
|------|------------|-------|
| `PointSource()` | — | unresolved |
| `UniformDisk(diameter = d)` | `diameter` | mas |
| `LimbDarkenedDisk(diameter = d, law = :linear, ld1 = 0.3)` | `diameter, ld1[, ld2]` | `law` ∈ `:linear`, `:quadratic`, `:power` |
| `GaussianDisk(fwhm = w)` | `fwhm` | circular |
| `EllipticalGaussian(fwhm = w, ratio = r, pa = θ)` | `fwhm, ratio, pa` | `pa` in degrees |

The disk and Gaussian visibilities come from OITOOLS (`visibility_ud`, `visibility_ldlin`,
…), so a ROTIR orbit fit and an OITOOLS model fit use the same analytic profiles.

## Choosing what to fit

Parameters are named: the eight orbital elements `a, i, Omega, omega, e, P, T0, dP`, each
component's own parameters prefixed `c1_` / `c2_`, and the flux ratio `f`.

```julia
res = fit_orbit(data, comp1, comp2;
                elements = (...),
                free   = [:a, :i, :Omega, :T0, :f, :c1_diameter, :c2_fwhm],
                bounds = Dict(:a => (0.5, 1.5), :f => (0.1, 3.0)))
```

The default frees `a, i, Omega, T0, f` and every component parameter. Deliberately **not**
free by default:

* **`e` and `omega`** — for a circularised system `e = 0` and `omega` is undefined, so
  fitting it samples a direction the likelihood is exactly flat in.
* **`P`** — an eclipse or spectroscopic ephemeris usually pins the period far better than a
  few nights of visibilities can. See *Fitting the period* below.
* **`dP`** — needs a baseline long enough for the phase drift to accumulate.

## Degeneracies

Three are handled by the default bounds, and all three are exact rather than approximate:

**`(Ω, T0) → (Ω + 180°, T0 + P/2)` at `e = 0`.** Both give identical observables. `Omega`
is restricted to `[0, 180)` while `T0` spans one full period — together exactly one
fundamental domain, each physical solution appearing once. Halving *both* would exclude
real solutions.

**Elliptical Gaussian position angle.** The profile is symmetric under a half turn, so
`pa ∈ [0, 360)` would hold two copies of every solution. The default bound is 180°.

**T0 and disc PA are circular.** With a nested sampler, mark them so a mode sitting on a
prior boundary is not split in two:

```julia
# names of the free parameters, in order
wrapped = [n in (:T0, :c2_pa) for n in res.spec.names[res.spec.free]]
```

## Fitters

```julia
fit_orbit(...; method = :neldermead)   # NLopt, point estimate, fast
fit_orbit(...; method = :ultranest)    # nested sampling, posterior + log Z
```

`:ultranest` additionally returns `posterior` (nsamples × nfree), `logz`, `logzerr`, `q16`,
`q84`. For a peaked posterior — a poorly fitting model, or many parameters — keep
`use_stepsampler = true` (the default): plain region rejection stalls once the likelihood
occupies a vanishing fraction of the prior box.

!!! note "A step sampler makes `vectorized` moot"
    UltraNest's step samplers propose one point at a time by construction, so with a step
    sampler every likelihood call receives a batch of exactly one row. Any parallelism
    across the batch is then inert; parallelise *inside* the likelihood instead, or use
    `use_stepsampler = false` where rejection sampling is viable.

## Fitting the period

`P` can be freed as a consistency check:

```julia
res = fit_orbit(data, comp1, comp2;
                elements = (...),
                free   = [:a, :i, :Omega, :T0, :f, :P],
                bounds = Dict(:P => (12.9404, 12.9424)))
```

Scale the prior to the physics, not to what the data might detect. The observable is
accumulated phase drift, `n·δP/P` over `n` orbits — so a short baseline simply cannot see a
period error, and a flat posterior on `P` is the correct answer rather than a failure. A
prior much wider than the amount `P` itself varies over the baseline (via `dP`) asserts an
uncertainty the ephemeris contradicts.

What the check *does* establish: whether the other elements shift when `P` is freed. If they
do not, the orbital solution does not secretly depend on the assumed period.

## Performance

One Kepler solve per **distinct time**, not per uv point. An OIFITS night has thousands of
uv samples across tens of exposures, so `orbit_fit_data` reduces the times to their unique
values and gathers back — 32× fewer solves on β Lyr, and ~0.05 ms per likelihood over seven
nights. That is what makes nested sampling over ten parameters practical.

## Evaluating without fitting

Pass `free = Symbol[]` to compute χ² at the given elements:

```julia
res = fit_orbit(data, comp1, comp2; elements = (...), free = Symbol[], verbose = false)
res.chi2, res.chi2_split      # total, and (V², T3amp, T3φ)
```

## Tessellated components

`model = :tessellated` is reserved for fitting orbits with resolved ROTIR surfaces (Roche
shapes, gravity darkening, irradiation) and is **not yet implemented** — it needs its own
parameterisation rather than the analytic component library. For a tessellated binary today,
drive `create_binary_geometry` and `binary_chi2_f` directly; see
`demos/spica_binary_roche.jl`.

## API

```@docs
fit_orbit
OrbitComponent
PointSource
UniformDisk
LimbDarkenedDisk
GaussianDisk
EllipticalGaussian
orbit_fit_data
orbit_fit_spec
orbit_model_cvis
orbit_chi2
```
