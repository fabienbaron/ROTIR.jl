# Model atmospheres and limb darkening

A limb-darkening *law* — linear, quadratic, Hestroffer, Claret — is a fitted
approximation. A model *atmosphere* gives the specific intensity
``I(T_\mathrm{eff}, \log g, \mu, \lambda)`` directly, so the limb darkening is
**predicted rather than fitted** and its coefficients leave the parameter vector
entirely.

That matters for a gravity-darkened star. Gravity darkening and limb darkening
both remove flux from the limb, an interferometer measures only their sum, and a
fit that is free in both buys the difference out of the limb-darkening
coefficient. On β Cas the correlation between the gravity-darkening exponent β
and `ld1` runs from +0.67 to +0.99 depending on rotation rate; with the intensity
predicted it disappears, because `ld1` is no longer there to absorb anything.

## Where `log g` comes from: the mass is derived, not fitted

## Where `log g` comes from: the mass is derived, not fitted

An atmosphere is indexed by ``\log g``, and a surface gravity needs a mass. ROTIR
does not fit one. For a Roche rotator the rotation rate, the polar radius and the
mass are not independent — ``\Omega = f_\mathrm{ev}\sqrt{8GM/27R_p^3}`` — so with
``\Omega = 2\pi/P`` known from the rotation period the mass follows:

```math
M = \frac{27\,R_p^3\,\Omega^2}{8\,G\,f_\mathrm{ev}^2}
```

`rpole` is an *angular* radius, so the only new parameter this needs is a
**distance** `d` in parsecs, which converts it to a length. That is a quantity
interferometry cannot measure but Gaia already has: the data constrain `d` to
only ~12 pc through its effect on ``\log g``, against ~0.1 pc from a parallax,
so it is pinned externally rather than fitted against the star.

```julia
p = (surface_type = 2, rpole = 0.849, d = 16.8, frac_escapevel = 0.92,
     rotation_period = 0.893, tpole = 7208.0, beta = 0.25, inclination = 20.0,
     position_angle = 115.0, gravity_law = 2, ldtype = 0,
     ld1 = 0.0, ld2 = 0.0, ld3 = 0.0, ld4 = 0.0)

derived_quantities(p)      # mass, R_p, R_eq, logg_pole, v_eq, vsini
println(derived_summary_text(p))
```

For β Cas those four inputs give **M = 1.936 M⊙, R_p = 3.067 R⊙,
R_eq = 3.796 R⊙, ``\log g_\mathrm{pole}`` = 3.752 and vsini = 73.2 km/s** — every
one matching the literature for an F2III-IV star, and **vsini was never an
input**. It falls out of the geometry, which is the cross-check that the
parameterization is self-consistent rather than merely sufficient.

Propagating a Fisher covariance through `derived_quantities` at that solution
gives **M = 1.936 ± 0.042 M⊙ (2.1 %)**, with the polar radius to 0.63 % and the
equatorial to 0.60 % — a per-cent-level mass for a *single* star, from
interferometry plus a parallax and a known rotation period.

Two consequences worth stating plainly:

* **``\log g`` varies over the surface, by 0.466 dex pole to equator on β Cas.**
  That spread is the whole reason a scalar `logg` cannot serve this model, and it
  is why the grid must cover the range the *surface* spans rather than the polar
  values. [`logg_map`](@ref) returns the per-tessel values and carries its own
  adjoint, so `d` is a genuinely differentiable parameter.
* **`rotation_period` and `frac_escapevel` can disagree**, because before this
  they were read by different functions and never compared. An implausible
  derived mass is the visible symptom, and `advise_star_params` reports it — but
  no optimiser consults that, so a fit is free to walk to 7 M⊙. For a production
  fit, either pass a `logprior` penalising an implausible mass or free
  `rotation_period`.

Without a distance the derived quantities are simply unavailable:
`has_physical_scale(p)` is false, `derived_summary_text` returns an empty string,
and an atmosphere provider cannot be used at all. A `PlanckProvider` still can —
it ignores ``\log g``, and ``\partial\log\pi/\partial d`` is then exactly zero,
which the test suite asserts as a null-space gate.

## The provider interface

Everything below produces an [`IntensityProvider`](@ref), which the fitting and
sampling paths consume identically:

```julia
I = provider_intensity(provider, Teff, logg, μ, λ)
```

`owns_mu(provider)` says whether the provider already supplies the μ dependence.
When it does, set **`ldtype = 0`** on the star so no limb-darkening factor is
applied on top — applying both counts the limb twice, and nothing downstream
raises if you do. `check_provider_consistency` reports the mismatch.

A provider needs a per-tessel `logg`, which needs a mass, which ROTIR *derives*
from `(rpole, d, frac_escapevel, rotation_period)`. So an atmosphere backend
requires a distance `d` on the star. See [`derive_mass`](@ref).

## Choosing a backend

| | Teff | per-μ? | geometry |
|---|---|---|---|
| Korg / MARCS | 2800–8000 K | synthesized at any μ | plane-parallel or spherical |
| Kurucz ATLAS9 | 3500–50000 K | 17 fixed angles | plane-parallel |

Build the grid over the range the **whole surface** spans, not the polar values:
on β Cas the surface runs Teff 5888–7206 K and ``\log g`` 3.29–3.75, a 0.466 dex
spread in gravity, which is why a scalar `logg` cannot serve this model.

The two backends have a page each: [Korg and MARCS](@ref) and [Kurucz ATLAS9](@ref).
