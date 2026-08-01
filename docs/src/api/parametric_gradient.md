# Parametric model and AD

A fully parametric forward model of a rapid rotator, built from primitives that each
carry a hand-coded `ChainRulesCore.rrule`. Unlike [Shape gradients](shape_gradient.md),
where the surface map is a free per-tessel vector, here the map is *determined* by the
physical parameters through von Zeipel gravity darkening — so the whole likelihood is a
function of a handful of numbers and is suitable for parameter estimation and MCMC.

ROTIR does **not** depend on an AD backend: the primitives only define ChainRules
rrules, so any ChainRules-compatible AD composes them. Load Zygote (or another backend)
in your own environment:

```julia
using ROTIR, Zygote
logπ = build_parametric_logπ(data, tessels, tepochs, star_params)
g    = Zygote.gradient(logπ, θ)[1]
```

The sampling stack used by `demos/rapid_rotator_betCas_nuts.jl` (AdvancedHMC,
Pathfinder, LogDensityProblems) lives in `demos/Project.toml`.

## Parameter vector

| Index | Parameter | Units | Notes |
|-------|-----------|-------|-------|
| 1 | `rpole` | mas | Polar radius |
| 2 | `ω` | – | Fractional escape velocity, 0 ≤ ω < 1 (break-up at 1) |
| 3 | `inc` | degrees | Spin-axis inclination |
| 4 | `PA` | degrees | Spin-axis position angle |
| 5 | `β` | – | von Zeipel gravity-darkening exponent (0.25 radiative, 0.08 convective) |
| 6 | `ld1` | – | First limb-darkening coefficient |
| 7 | `ld2` | – | Second limb-darkening coefficient (used by `ldtype = 2`) |
| 8 | `tpole` | K | Only present when `tpole_free = true` |

`ldtype` is taken from `base_params` and is not fitted (it selects the law, it is not a
continuous parameter).

## Composed log-density

| Function | Description |
|----------|-------------|
| `build_parametric_logπ(data_epochs, tessels, tepochs, base_params; kwargs...)` | Returns a differentiable closure `θ -> -0.5·χ²(θ) + logprior(θ)` |

### Keywords

| Keyword | Default | Description |
|---------|---------|-------------|
| `intensity_model` | `:linear` | `:linear` (I = T, historical proxy) or `:planck` |
| `band` | `nothing` | Wavelength in metres, required for `:planck` |
| `κ` | `50` | Sigmoid sharpness of the soft visibility |
| `GM` | `1` | Gravitational parameter used by the von Zeipel map |
| `tpole_free` | `false` | Append `tpole` to θ as an 8th parameter |
| `logprior` | `nothing` | Optional `θ -> log p(θ)` added to the log-likelihood |

All per-epoch constants (UV frequencies, `k2_inv_im`) are captured once when the closure
is built, so repeated evaluations do no setup work.

!!! note "tpole is degenerate under the linear intensity model"
    With `intensity_model = :linear` the visibilities are normalized by the total flux,
    so a global temperature scale cancels and `∂logπ/∂tpole ≡ 0`. Use
    `intensity_model = :planck` with a `band` to make `tpole` identifiable. This
    null-space is asserted in `test/test_parametric_gradient.jl`.

## Primitives

Each forward function is paired with a `*_and_derivs` routine that returns the analytic
derivatives; the rrule wraps the latter.

| Function | Description |
|----------|-------------|
| `vonzeipel_map(rpole, fev, β, tpole, sinθ, cosθ; GM)` | Per-tessel temperature map from von Zeipel gravity darkening |
| `vonzeipel_map_and_derivs(...)` | Map plus derivatives w.r.t. `rpole`, `fev`, `β`, `tpole` |
| `intensity(x, model, band)` | Temperature map → intensity (`:linear` or `:planck`) |
| `planck_and_dT(T, λ)` | Planck function and its temperature derivative |
| `ld_weight(nz, ldtype, ld1, ld2)` | Limb-darkening map from the normal z-component |
| `ld_and_derivs(nz, ldtype, ld1, ld2)` | LD map plus `∂/∂nz`, `∂/∂ld1`, `∂/∂ld2` |
| `visibility_weight(nz, κ)` | Soft visibility `σ(κ·nz)` |
| `project_geometry(rpole, fev, inc, PA, tessels, t, base_params)` | Projected quad vertices and normal z-components |
| `interferometric_chi2(xw, proj_west, proj_north, kx, ky, k2_inv_im, data)` | Single-epoch χ² (V² + T3amp + T3phi) from weighted tessels |
| `limb_mu(nz)`, `mu_and_dmu(nz)` | Limb cosine `μ = max(nz, 0)` and its derivative |

`limb_mu` is the single source of truth for the limb cosine, shared with
`compute_ldmap` in the forward geometry — so the parametric model and
`create_star` apply exactly the same limb-darkening law.

## How the pieces compose

```julia
x    = vonzeipel_map(rpole, ω, β, tpole, sinθ, cosθ)   # temperature map
Imap = intensity(x, intensity_model, band)             # → intensity
pw, pn, nz = project_geometry(rpole, ω, inc, PA, tessels, t, base)
xw   = Imap .* visibility_weight(nz, κ) .* ld_weight(nz, ldtype, ld1, ld2)
chi2 = interferometric_chi2(xw, pw, pn, kx, ky, k2_inv_im, data)
```

The `interferometric_chi2` pullback reuses the same adjoint kernels as
`shape_chi2_fg!` (`compute_adjoint_cvis!`, `compute_adjoint_vertices!`), including the
flux-normalization and shoelace corrections, so both paths share one validated
implementation.

## Validation

`test/test_parametric_gradient.jl` checks, against central finite differences:

1. the leaf derivatives (von Zeipel map, all three LD laws w.r.t. `nz`, `ld1`, `ld2`),
2. forward consistency — `logπ ≈ -0.5·spheroid_parametric_f`, i.e. the parametric path
   agrees with the standard χ² path on the same model,
3. the end-to-end AD gradient for every LD type,
4. the `tpole` null-space gate above, and that Planck breaks it.

Run the suite with `julia --project=. -e 'using Pkg; Pkg.test()'`.
