# Bootstrap uncertainties

Parameter uncertainties for parametric fits by block bootstrap: refit replicates in which
the observations are resampled in blocks, and report the scatter of the best-fit
parameters. It is not available for image reconstruction, where the scatter would not
include the regularisation bias.

The blocks and the replicate statistics come from OITOOLS (`data_blocks`,
`resample_blocks`, `bootstrap_driver`). The calibration of the resampling schemes is in
`OITOOLS/demos/bootstrap_validation`.

## What a block is

A replicate resamples which observations are present, not their values. All wavelength
channels of a block stay together.

| `granularity` | one block per | typical count (ρ Cas example) |
|---------------|---------------|-------------------------------|
| `:config` (default) | (MJD, telescope configuration) — baseline for V², triangle for T3 | 746 |
| `:epoch` | MJD, i.e. one exposure with all its baselines | 113 |
| `:point` | single data point (i.i.d. bootstrap) | 6763 |

`:point` ignores correlations between points and underestimates uncertainties on real data.
`readoifits` stores MJDs in `Float64` for any `T`, so blocks are the same for `Float32` and
`Float64` data.

## Resampling schemes

| `mode` | what it does | notes |
|--------|--------------|-------|
| `:replacement` (default) | draw `nblocks` blocks with replacement | textbook nonparametric bootstrap |
| `:halfsample` | keep a random half of the blocks | conservative when blocks are few |
| `:weights` | multiplier ("Bayesian") bootstrap: divide each block's error bars by `√w` | never drops a block; the most optimistic of the calibrated schemes |

`:weights` is about 10 % optimistic at ~100 blocks.

## Multi-epoch data

ROTIR data is a `Vector{OIdata}`, one entry per epoch. `resample_epochs` resamples every
epoch and keeps their order.

- `stratify=true` (default) draws blocks independently within each epoch, so each rotation
  phase keeps about its data volume.
- `stratify=false` draws once over all blocks; an epoch can lose most of its data, and a
  warning is issued when one drops below a tenth of its points.

## Fitting

```julia
using ROTIR, Zygote          # fit_parametric needs Zygote

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

`free` takes names or indices; other entries stay at their `θ0` value. With `ω = 0` the
rapid-rotator model is a limb-darkened sphere, and inclination, position angle, β and
`tpole` have no effect; leave them fixed, or replicates collect on the bounds.

Every replicate starts at the full-data θ̂. On ρ Cas, a fit started far from it finds a
different minimum (diameter 5.2 mas instead of 2.5 mas).

`bootstrap_parametric` with data and a model also takes the `fit_parametric` keywords (`lb`,
`ub`, `intensity_model`, `band`, `κ`, `GM`, `tpole_free`, `gravity_law`, `logprior`,
`provider`, `layout`, `maxiter`, `gtol`, `mem`) and `refit_full` (default `true`; `false`
when `θ0` is already the full-data fit).

## Reading the result

`ParametricBootstrap` wraps OITOOLS' `BootstrapResult` (`median`, `sigma`,
`sigma_plus`/`sigma_minus`, `covar`, `correlation`, `samples`, `mask`, `nfailed`, read
directly from the result) and adds:

| field | meaning |
|-------|---------|
| `natbound` | kept replicates that converged onto a box bound, per parameter |
| `ndegenerate` | replicates whose resampled data lost all closure quantities |
| `nmaxiter` | replicates whose fit stopped on the iteration limit |
| `nblocks_per_epoch`, `stratify` | the resampling setup |

A large `natbound` means the bounds shape the percentiles. A non-zero `nmaxiter` means some
replicates did not reach a minimum; raise `maxiter`.

`sigma` is the half-width of the 16–84 % range, not a standard deviation, so a few
replicates in a distant minimum do not inflate it. `sigma_clipping` (e.g. 4.5) rejects such
replicates; `count(b.mask)` gives the number kept.

## Bootstrap versus posterior

`demos/rapid_rotator_betCas_nuts.jl` samples the posterior of the same model with NUTS.

- The posterior assumes the quoted error bars are correct and independent.
- The block bootstrap does not: correlated or mis-stated errors appear as extra scatter
  between replicates.

A ratio beyond about √2 between the two indicates wrong error bars. Neither accounts for a
wrong model: ρ Cas fits at χ²ᵣ = 2.8 even with oblateness and gravity darkening free.

## Without an AD backend

`bootstrap_parametric` also takes a closure `replicate_data -> θ̂`, which needs no Zygote:

```julia
b = bootstrap_parametric(rdata -> my_fit(rdata), data; nboot=200, seed=1)
```
It may return `(θ̂, chi2r)` or `(θ̂, chi2r, extra)`; `chi2r` enables the `chi2r_max`
rejection, and `extra` is kept for diagnostics.

`parallel = :outer` (default) runs replicates on all threads and requires a closure that is
safe to call concurrently; a failing replicate is dropped without an error. `parallel =
:inner` runs replicates one at a time. Other keywords: `nboot` (200), `mode`, `granularity`,
`mjd_digits` (5), `stratify`, `lb`/`ub` (to count replicates at a bound), `list_free_params`,
`x_opt`, `sigma_clipping`, `chi2r_max`, `seed`, `verb`.
