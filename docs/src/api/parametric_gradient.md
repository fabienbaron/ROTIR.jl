# Parametric model and AD

A parametric forward model of a rapid rotator, built from primitives that each carry a
hand-coded `ChainRulesCore.rrule`. The temperature map is set by the physical parameters
through gravity darkening, so the likelihood is a function of a few numbers and suits
parameter estimation and sampling. In [Shape gradients](shape_gradient.md) the map is a free
per-tessel vector instead.

ROTIR defines only the rrules and has no AD dependency; any ChainRules-compatible AD
backend composes them. Load one yourself:

```julia
using ROTIR, Zygote
logπ = build_parametric_logπ(data, tessels, tepochs, star_params)
g    = Zygote.gradient(logπ, θ)[1]
```

`demos/rapid_rotator_betCas_nuts.jl` samples this model with Pathfinder and AdvancedHMC;
its packages are in `demos/Project.toml`.

## Parameter vector

| Index | Parameter | Units | Notes |
|-------|-----------|-------|-------|
| 1 | `rpole` | mas | polar radius |
| 2 | `ω` | – | `frac_escapevel`, Ω/Ω_crit, 0 ≤ ω < 1 (break-up at 1) |
| 3 | `inc` | degrees | spin-axis inclination |
| 4 | `PA` | degrees | spin-axis position angle |
| 5 | `β` | – | gravity-darkening exponent (0.25 radiative, 0.08 convective) |
| 6 | `ld1` | – | first limb-darkening coefficient |
| 7 | `ld2` | – | second limb-darkening coefficient (read by `ldtype = 2` and `4`) |
| 8 | `tpole` | K | only with `tpole_free = true` |

`ldtype` comes from `base_params` and is not fitted. A `ParametricLayout` (from
`parametric_layout`) changes the vector: `ldtype = 0` drops `ld1` and `ld2`, and
`distance_free = true` appends the distance `d` (pc).

## Composed log-density

`build_parametric_logπ(data_epochs, tessels, tepochs, base_params; kwargs...)` returns a
differentiable closure `θ -> -0.5·χ²(θ) + logprior(θ)`. Per-epoch constants (uv frequencies,
`k2_inv_im`) are computed once, when the closure is built.

| Keyword | Default | Description |
|---------|---------|-------------|
| `intensity_model` | `:linear` | `:linear` (I = T) or `:planck` |
| `band` | `nothing` | wavelength in metres, required for `:planck` and with a `provider` |
| `κ` | `50` | sharpness of the soft visibility σ(κ·nz) |
| `GM` | `1` | gravitational parameter of the gravity-darkening map |
| `tpole_free` | `false` | append `tpole` to θ |
| `gravity_law` | `nothing` | `:vonzeipel` or `:elr`; `nothing` reads `base_params.gravity_law` |
| `logprior` | `nothing` | `θ -> log p(θ)`, added to the log-likelihood |
| `provider` | `nothing` | an [`IntensityProvider`](@ref): intensity from a model atmosphere, `I(Teff, logg, μ, band)` |
| `layout` | `nothing` | a `ParametricLayout`; the default is the table above |
| `weights` | `(1, 1, 1)` | weights of the V², T3amp and T3φ χ² terms |

With a `provider`, `ldtype` must be 0: the provider supplies the μ dependence and no
limb-darkening law is applied. `logg` is computed per tessel by [`logg_map`](@ref) from
`rpole`, `d`, `ω` and `rotation_period`, so fitting `d` needs `distance_free = true` in the
layout.

!!! note "tpole and the linear intensity model"
    With `intensity_model = :linear` the visibilities are normalised by the total flux, so a
    global temperature scale cancels and `∂logπ/∂tpole = 0`. Use `intensity_model = :planck`
    with a `band` to constrain `tpole`. `test/test_parametric_gradient.jl` checks this.

## Primitives

Each forward function has a `*_and_derivs` routine that returns the analytic derivatives;
the rrule calls it.

| Function | Description |
|----------|-------------|
| `vonzeipel_map(rpole, fev, β, tpole, sinθ, cosθ; GM)` | per-tessel von Zeipel temperature map |
| `vonzeipel_map_and_derivs(...)` | map and derivatives with respect to `rpole`, `fev`, `β`, `tpole` |
| `intensity(x, model, band)` | temperature to intensity (`:linear` or `:planck`) |
| `planck_and_dT(T, λ)` | Planck function and its temperature derivative |
| `ld_weight(nz, ldtype, ld1, ld2, ld3=0, ld4=0)` | limb-darkening map from the normal z-component |
| `ld_and_derivs(nz, ldtype, ld1, ld2, ld3=0, ld4=0)` | LD map and `∂/∂nz`, `∂/∂ld1`, `∂/∂ld2` |
| `visibility_weight(nz, κ)` | soft visibility `σ(κ·nz)` |
| `project_geometry(rpole, fev, inc, PA, tessels, t, base_params)` | projected quad vertices and normal z-components |
| `interferometric_chi2(xw, proj_west, proj_north, kx, ky, k2_inv_im, data, weights=(1,1,1), V0=1)` | single-epoch χ² (V², T3amp, T3φ) from weighted tessels |
| `limb_mu(nz)`, `mu_and_dmu(nz)` | limb cosine `μ = max(nz, 0)` and its derivative |

`compute_ldmap` in the forward geometry uses the same `limb_mu`, so the parametric model and
`create_star` apply the same limb-darkening law. `ld3` and `ld4` (Claret, `ldtype = 4`) are
read from `base_params` and held fixed.

`V0` in `interferometric_chi2` multiplies the complex visibility: V² scales by `V0²`, T3amp
by `V0³`, T3φ is unchanged. It absorbs a uniform scale on the measured visibilities
(calibrator diameter, coherence loss). It trades against the angular diameter unless the
data reach the first visibility null; check `corr(V0, rpole)`. On ζ And, Roettenbacher et al.
report 1.02 ± 0.02 and a free `V0` here gives 1.020–1.022.

## How the pieces compose

```julia
x    = vonzeipel_map(rpole, ω, β, tpole, sinθ, cosθ)   # temperature map
Imap = intensity(x, intensity_model, band)             # intensity
pw, pn, nz = project_geometry(rpole, ω, inc, PA, tessels, t, base)
xw   = Imap .* visibility_weight(nz, κ) .* ld_weight(nz, ldtype, ld1, ld2)
chi2 = interferometric_chi2(xw, pw, pn, kx, ky, k2_inv_im, data)
```

The `interferometric_chi2` pullback uses the same adjoint kernels as `shape_chi2_fg!`
(`compute_adjoint_cvis!`, `compute_adjoint_vertices!`), including the flux-normalisation
and shoelace terms.

## Validation

`test/test_parametric_gradient.jl` checks against central finite differences:

1. the leaf derivatives (von Zeipel map; the LD laws with respect to `nz`, `ld1`, `ld2`);
2. forward consistency: `logπ ≈ -0.5·spheroid_parametric_f` on the same model;
3. the end-to-end AD gradient for every LD type;
4. `∂logπ/∂tpole = 0` under the linear model, and non-zero under Planck.

## The ellipsoid's θ layout

[`fit_ellipsoid_ld`](@ref) fits the first five entries of this layout (the three semi-axes,
inclination, position angle) and `ld1`.

```@docs
ellipsoid_param_names
```
