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

- `(Omega, omega) → (Omega + 180°, omega + 180°)`: the two orbits give the same sky positions
  at every epoch and differ only in which star is in front. Interferometry, closure phases
  included, cannot separate them at any eccentricity. `Omega` is bounded to `[0, 180)` and
  `T0` to one period centred on its starting value, so each solution appears once.
- `EllipticalGaussian` `pa`: the profile is symmetric under a half turn, so `pa` is bounded to
  `[0, 180)`.

A published `Omega = 309.9°` therefore corresponds to `Omega = 129.9°` with `omega` shifted by
180° in an interferometry-only fit. The Orbit tab's fit makes this change itself when `Omega`
is free, its box is `[0, 180)` and its value lies outside it; a box set to `[0, 360)` is left as
it is. To keep a published `Omega` in a script, widen its box:

```julia
bounds = Dict(:Omega => (0.0, 360.0))
```

With radial velocities (`rv = ...`) the mirror is broken and `Omega` spans `[0, 360)`; see
[Radial velocities](@ref "Radial velocities in the fit").

`T0`, `omega` and `pa` are periodic, but `fit_orbit` does not tell the sampler so; a mode at a
bound is split in two. When the solution lies near `omega = 0°`, move its box, e.g.
`bounds = Dict(:omega => (180.0, 540.0))`. `demos/betlyr/betlyr_orbit_fit_ultranest.jl` calls UltraNest directly with
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

## Radial velocities in the fit

Radial velocities are read from a CSV table, the format OITOOLS reads:

```
# HD 210763, Gallenne et al. 2023
MJD,RV,sigma,star
59510.16368,6.805,0.072,a
59510.16368,24.530,0.172,b
```

| column | content |
|---|---|
| `MJD` (or `JD`) | time |
| `RV`, `sigma` | velocity and its error, km/s, positive receding |
| `star` | `a` primary, `b` secondary |
| `instrument` | optional |

Columns are found by name in any order; extra columns are ignored; `#` lines are comments.
`read_rv(path1, path2)` reads the older format of two whitespace files `JD RV sigma`, one
per star, and `write_rv` writes the CSV.

```julia
rv  = read_rv("demos/data/hd210763_rv_gallenne2023.csv")
res = fit_orbit(data, comp1, comp2; elements = (...), rv = rv,
                rv_params = (plx = 4.5, q = 0.86, gamma = 15.0))
res.rvpar          # (plx, q, gamma, jit_a, jit_b, K1, K2, M, M1, M2, d)
res.chi2_split     # (V², T3amp, T3φ, RV)
```

The velocity model is `binary_RV`, with ω the argument of periapsis of the relative
orbit:

    V_a = γ − K1 [cos(ν + ω) + e cos ω],    V_b = γ + K2 [cos(ν + ω) + e cos ω]

This is the convention of OITOOLS' `orbit_rv`; its `MJD0` is `T0 − 2400000.5`.

`rv` adds these parameters after `f`:

| name | meaning | starting value | box |
|---|---|---|---|
| `plx` | parallax, mas | from the data's amplitudes at the starting elements | [1e-4, 1e4] |
| `q` | mass ratio M2/M1 | idem | [0.01, 100] |
| `gamma` | systemic velocity, km/s | mean velocity | [−1000, 1000] |
| `dgamma_<instrument>` | offset of each instrument after the first, km/s | 0 | [−50, 50] |
| `jit_a`, `jit_b` | jitter of each star, km/s | 0 | [0, 50] |

The semi-amplitudes are not parameters. The interferometry fits `a`, `i`, `P` and `e`, and the
amplitudes follow from them, the parallax and the mass ratio ([`rv_amplitudes`](@ref)):

    K = 2π (a/plx) sin i / (P √(1 − e²)),    K1 = K q/(1 + q),    K2 = K/(1 + q)

with `a/plx` in AU and `P` in days. `res.rvpar` also reports `K1`, `K2`, the distance `d = 1000/plx`
(pc) and the dynamical masses `M = (a/plx)³/(P/365.25)²`, `M1`, `M2` (M☉,
[`dynamical_masses`](@ref)). A known parallax (Gaia, say) is fixed by leaving `plx` out of
`free`. `rv_parametrization = :amplitudes` fits `K1`, `K2` (km/s) instead, which is what a fit
with the interferometry switched off needs.

`plx`, `q`, `gamma` and the offsets are free by default; the jitters are fitted only when
listed in `free`. An instrument offset exists only when the table has an `instrument`
column with more than one value; the first instrument in the file is the reference.

The jitter `s` is added in quadrature to the tabulated error, `σ² + s²`. The fit minimises
the interferometric χ² plus [`rv_nll2`](@ref),

    Σ (rv − model)²/(σ² + s²) + Σ ln((σ² + s²)/σ²),

which equals χ²_RV when the jitters are 0. `chi2_split[4]` and `chi2_red` report χ² with the
inflated errors; `chi2_red` divides by the interferometric plus RV point count.

`weights = [0, 0, 0]` switches the interferometry off: the visibilities are then not
computed, and the fit uses the velocities alone.

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
| `demos/spica_global_orbit_fit.jl` | Spica — two limb-darkened discs, with radial velocities | Nautilus over the whole prior box, apsidal motion, periodic `Omega` and `omega` |

The β Lyrae script defaults to UltraNest; `METHOD=neldermead` runs a quick Nelder–Mead fit.
Neither script models occultation or mutual irradiation. The Spica result is the starting
point of `demos/spica_refine_orbit_fit.jl`, which refines it with Roche surfaces outside
`fit_orbit`.

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
orbit_periodic
orbit_residuals
resolve_params
compile_ties
apply_ties
OrbitTies
RVData
read_rv
write_rv
rv_instruments
rv_model
rv_residuals
rv_chi2
rv_nll2
rv_start
rv_amplitudes
rv_physical
dynamical_masses
```
