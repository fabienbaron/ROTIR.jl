# Bootstrap uncertainties

Parameter uncertainties for **parametric** fits by nonparametric block bootstrap: refit
many replicates in which the observations are resampled in blocks, and report the scatter
of the best-fit parameters.

This is deliberately not offered for image reconstruction. A bootstrap around a regularized
inversion measures how much the data move the image while saying nothing about the
regularization bias, and per-pixel percentiles from such a run read as error bars without
being any.

The resampling units and the replicate statistics come from OITOOLS
(`data_blocks`, `resample_blocks`, `bootstrap_driver`), so ROTIR and OITOOLS quote
uncertainties computed the same way. The calibration study behind the schemes lives in
`OITOOLS/demos/bootstrap_validation`.

## What a block is

A replicate resamples *which observations you have*, not their values. The unit is a
block, and all wavelength channels of a block stay together.

| `granularity` | one block per | typical count (ρ Cas example) |
|---------------|---------------|-------------------------------|
| `:config` (default) | (MJD, telescope configuration) — baseline for V², triangle for T3 | 746 |
| `:epoch` | MJD, i.e. one exposure with all its baselines | 113 |
| `:point` | single data point (i.i.d. bootstrap) | 6763 |

`:point` destroys the correlation structure and **underestimates** uncertainties on real
data; it is there for comparison studies. MJDs are stored in `Float64` by `readoifits`
whatever `T` is, so blocks are identical for `Float32` and `Float64` data.

## Resampling schemes

| `mode` | what it does | notes |
|--------|--------------|-------|
| `:replacement` (default) | draw `nblocks` blocks with replacement | textbook nonparametric bootstrap |
| `:halfsample` | keep a random half of the blocks | conservative when blocks are few |
| `:weights` | multiplier ("Bayesian") bootstrap: divide each block's error bars by `√w` | never drops a block; the most optimistic of the calibrated schemes |

`:replacement` is the default. `:weights` would keep every array length fixed across
replicates, but that buys nothing in ROTIR — `resample_blocks` leaves `data.uv` untouched
and the polygon FT runs over every uv row regardless — while being ~10% optimistic at ~100
blocks.

## Multi-epoch data

ROTIR data is a `Vector{OIdata}`, one entry per epoch, each with its own geometry.
`resample_epochs` resamples all of them and preserves their order.

- `stratify=true` (default) draws blocks independently within each epoch. The epochs *are*
  the rotational phase coverage that constrains inclination, position angle and ω, so an
  unstratified draw that empties one inflates the uncertainties for a reason unrelated to
  data quality.
- `stratify=false` draws once over the pooled block set; an epoch can then lose most of its
  data, and a warning is issued when it does.

## Fitting

```julia
using ROTIR, Zygote          # the fitter lives in an extension: no AD, no fit

data    = readoifits_multiepochs(["rho_Cas_example.oifits"]; T=Float64)[1, :]
tessels = tessellation_healpix(3, T=Float64)
base    = (surface_type=2, rpole=1.25, tpole=4000.0, ldtype=3, ld1=1.75, ld2=0.0,
           inclination=0.0, position_angle=0.0, rotation_period=1e6,
           beta=0.08, frac_escapevel=0.0, B_rot=0.0)

θ0 = [1.25, 0.0, 0.0, 0.0, 0.08, 1.75, 0.0]     # rpole ω inc PA β ld1 ld2
θ̂, chi2r, info = fit_parametric(data, tessels, [0.0], base;
                                θ0 = θ0, free = ["rpole", "ld1"])

b = bootstrap_parametric(data, tessels, [0.0], base;
                         θ0 = θ0, free = ["rpole", "ld1"],
                         nboot = 200, seed = 42, sigma_clipping = 4.5)
b.median, b.sigma, b.correlation
```

**Fit a subset.** `free` takes names or indices; everything else is held at its `θ0` value.
This is the normal case: a resolved single star constrains its angular size and its limb
darkening and little else. With `ω` frozen at 0 the rapid-rotator model *is* the sphere +
limb-darkening model (the radius reduces to `rpole` everywhere and a uniform von Zeipel map
makes inclination, position angle, β and tpole inert). Leaving inert parameters free lets
the optimizer wander a flat subspace and pile replicates onto the bounds.

**Warm start.** Every replicate starts at the full-data θ̂, so the spread measures the data
rather than the optimizer's path. This matters more than it sounds: on ρ Cas a fit started
far away lands in a different basin entirely (a diameter of 5.2 mas instead of 2.5 mas).

## Reading the result

`ParametricBootstrap` wraps OITOOLS' `BootstrapResult` — `median`, `sigma`,
`sigma_plus`/`sigma_minus`, `covar`, `correlation`, `samples`, `mask`, `nfailed` are all
reachable directly — and adds:

| field | meaning |
|-------|---------|
| `natbound` | kept replicates that converged onto a box bound, per parameter |
| `ndegenerate` | replicates whose resampled data lost all closure quantities |
| `nmaxiter` | replicates whose fit stopped on the iteration limit |
| `nblocks_per_epoch`, `stratify` | the resampling setup |

A large `natbound` means the percentiles are shaped by the bounds rather than by the data;
a nonzero `nmaxiter` means some replicates are not minimisers and are widening the
distribution with optimiser noise.

`sigma` is a **percentile half-width**, not a standard deviation: with a multimodal χ² a
handful of replicates can land in a distant basin, which inflates `std(samples)` while the
16/84 range stays put. Use `sigma_clipping` (e.g. 4.5) to reject them explicitly and check
how many were dropped via `count(b.mask)`.

## Bootstrap versus posterior

`demos/rapid_rotator_betCas_nuts.jl` samples the posterior of the same model with NUTS.
The two answer different questions and it is worth running both:

- the posterior takes the quoted error bars at face value and assumes they are independent;
- the block bootstrap does not — correlated calibration errors and mis-stated errors show
  up as extra scatter between replicates.

Disagreement beyond about √2 is a diagnostic that the quoted errors are wrong, not a bug.
Neither is a defence against a misspecified model: when χ²ᵣ is far from 1 (ρ Cas fits at
2.8 even with oblateness and gravity darkening free), both quantify parameters of a model
the star does not obey.

## Without an AD backend

`bootstrap_parametric` also takes any closure `replicate_data -> θ̂`, so you can bootstrap a
fitter of your own without loading Zygote:

```julia
b = bootstrap_parametric(rdata -> my_fit(rdata), data; nboot=200, seed=1)
```
It may also return `(θ̂, chi2r)` or `(θ̂, chi2r, extra)`; `chi2r` enables the `chi2r_max`
rejection, and `extra` is carried through for diagnostics.
