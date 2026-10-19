# Orbit Fitting

`fit_orbit` fits the relative orbit of a two-component system to interferometric observables.

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

Component 1 is at the origin and component 2 is displaced by the orbit:

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

All epochs are fitted with one set of orbital elements.

## Components

| Type | Parameters | Notes |
|------|------------|-------|
| `PointSource()` | — | unresolved |
| `UniformDisk(diameter = d)` | `diameter` | mas |
| `LimbDarkenedDisk(diameter = d, law = :linear, ld1 = 0.3)` | `diameter, ld1[, ld2]` | `law` ∈ `:linear`, `:quadratic`, `:power` |
| `GaussianDisk(fwhm = w)` | `fwhm` | circular |
| `EllipticalGaussian(fwhm = w, ratio = r, pa = θ)` | `fwhm, ratio, pa` | `pa` in degrees |

Default bounds: diameters and FWHM 0.01–20 mas, `ratio` 0.05–1, `pa` 0–180°, `ld1` 0–1,
`ld2` −0.5–1. The disk visibilities come from the OITOOLS model API (`ud`, `ldlin`, `ldquad`,
`ldpow`), so ROTIR and OITOOLS use the same profiles.

## Choosing what to fit

The parameters are the nine orbital elements `a, i, Omega, omega, e, P, T0, dP, domega`, the
component parameters prefixed `c1_` / `c2_`, and the flux ratio `f`. Each is free, fixed or
tied:

| | how | use |
|---|---|---|
| free | listed in `free`, fitted within `bounds` | constrained by the data |
| fixed | not listed; keeps its `elements` value | known from elsewhere (an ephemeris) |
| tied | listed in `ties` as an expression | determined by other parameters |

`P` is required in `elements`; the others default to `a = 1, i = 90, Omega = 0, omega = 90,
e = 0, T0 = 0, dP = 0, domega = 0`. Units: `a` in mas, angles in degrees, `P` and `T0` in days,
`dP` in days per day, `domega` in degrees per day.

```julia
res = fit_orbit(data, comp1, comp2;
                elements = (a = 0.9, i = 91.2, Omega = 73.5, P = 12.96, T0 = 2454283.0),
                free   = [:a, :i, :Omega, :T0, :f, :c1_diameter, :c2_fwhm],
                bounds = Dict(:a => (0.8, 1.0), :f => (0.1, 3.0)),
                ties   = Dict(:c2_pa => "-Omega"),
                method = :ultranest)
```

The default frees `a, i, Omega, T0, f` and every component parameter. Not free by default:

- `e` and `omega`: at `e = 0`, `omega` is undefined.
- `P`: usually known better from an eclipse or spectroscopic ephemeris. See
  [Fitting the period](@ref).
- `dP`: needs a long baseline for the phase drift to accumulate.
- `domega` (apsidal motion): needs a baseline of years. Spica's periapsis advances about
  2.6°/yr, so `omega` drifts about 21° across the 2007–2015 CHARA data.

## Tied parameters

A tie computes one parameter from the others. Example: a disc in the orbital plane has its
projected major axis along the line of nodes, so its position angle follows from `Omega`.

```julia
ties = Dict(:c2_pa => "-Omega")
```

Expressions are Julia code over the parameter names. A tie may reference another tie;
ties are evaluated in dependency order, and a cycle is an error.

```julia
ties = Dict(:c2_pa    => "-Omega",
            :c2_ratio => "abs(cosd(i))")     # a flat disc's projected axis ratio
```

Ties are compiled once with `RuntimeGeneratedFunctions`; resolving them costs one allocation
per likelihood evaluation.

A parameter listed as both free and tied is an error (`orbit_fit_spec` raises).

!!! warning "Position-angle conventions"
    `EllipticalGaussian`'s `pa` is the angle of the **minor** axis (the one scaled by
    `ratio < 1`), measured from +RA toward +Dec. `Omega` is measured North through East. The
    major axis lies on the line of nodes when `pa = -Omega`. A wrong tie does not show in the
    χ²: the fit adjusts the other parameters around it.

To test a tie, fit the parameter freely once and compare. On β Lyrae the free position angle
is 1.9° from the tied value, and the tie raises χ²/n from 3.154 to 3.189.

## Degeneracies

The default bounds remove two exact degeneracies:

- `(Omega, T0) → (Omega + 180°, T0 + P/2)` at `e = 0`: `Omega` is bounded to `[0, 180)` and
  `T0` to one period centred on its starting value, so each solution appears once.
- `EllipticalGaussian` `pa`: the profile is symmetric under a half turn, so `pa` is bounded to
  `[0, 180)`.

!!! warning "The `Omega` bound assumes `e = 0`"
    Eccentricity breaks the `(Omega, T0)` degeneracy, and so do closure phases. For an
    eccentric orbit, widen the bound:

    ```julia
    bounds = Dict(:Omega => (0.0, 360.0))
    ```

    With the default, a published `Omega = 309.9°` folds to `129.9°` and the fit finds the
    mirrored orbit. `demos/spica_orbit_fit_tied.jl` bounds `Omega` to 290–330°.

`T0` and `pa` are periodic, but `fit_orbit` does not tell the sampler so; a mode at a bound is
split in two. `demos/betlyr/betlyr_orbit_fit_ultranest.jl` calls UltraNest directly with
`wrapped_params` for these two.

## Fitters

```julia
fit_orbit(...; method = :neldermead)   # NLopt, point estimate (default)
fit_orbit(...; method = :ultranest)    # nested sampling, needs `using PythonCall`
fit_orbit(...; method = :nautilus)     # nested sampling, needs `using Nautilus`
```

The nested samplers also return `posterior` (nsamples × nfree), `logz`, `logzerr`, `q16` and
`q84`; the point estimate is the posterior median. `maxeval = 20_000` bounds Nelder–Mead;
`min_num_live_points = 400` and `use_stepsampler = true` go to UltraNest. Plain region
rejection (`use_stepsampler = false`) stalls on a sharply peaked posterior.

With a step sampler, UltraNest passes the likelihood one point at a time, so parallelism goes
inside the likelihood (`threaded = true`).

## Fitting the period

Free `P` as a consistency check:

```julia
res = fit_orbit(data, comp1, comp2;
                elements = (...),
                free   = [:a, :i, :Omega, :T0, :f, :P],
                bounds = Dict(:P => (12.9404, 12.9424)))
```

The data measure accumulated phase drift, `n·δP/P` over `n` orbits, so a short baseline gives
a flat posterior on `P`. Keep the bound on `P` near the range the ephemeris allows over the
baseline (including `dP`). The check shows whether the other elements move when `P` is free.

## Performance

`orbit_fit_data` solves Kepler's equation once per distinct time, not per uv point. On β Lyr
that is 50 times for 1613 uv points, and one likelihood over seven nights takes about 0.05 ms.

## Evaluating without fitting

`free = Symbol[]` computes χ² at the given elements:

```julia
res = fit_orbit(data, comp1, comp2; elements = (...), free = Symbol[], verbose = false)
res.chi2, res.chi2_split      # total, and (V², T3amp, T3φ)
```

## Tessellated components

`model = :tessellated` (ROTIR surfaces: Roche shapes, gravity darkening, irradiation) is not
implemented and raises an error. For a tessellated binary, use `create_binary_geometry` and
`binary_chi2_f`; see `demos/spica_binary_roche.jl`.

## Worked examples

| Script | System | Shows |
|---|---|---|
| `demos/betlyr/betlyr_orbit_fit_tied.jl` | β Lyrae — donor + accretion disc | a tie (`c2_pa = -Omega`) on a circular orbit |
| `demos/spica_orbit_fit_tied.jl` | Spica — two limb-darkened discs | eccentric orbit, `e`/`omega` free, apsidal motion with `APSIDAL=1` |

Both default to UltraNest; `METHOD=neldermead` runs a quick Nelder–Mead fit. Neither models
occultation or mutual irradiation.

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
OrbitFitSpec
orbit_fit_spec
orbit_model_cvis
orbit_chi2
resolve_params
compile_ties
apply_ties
OrbitTies
```
