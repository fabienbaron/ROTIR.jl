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

## Korg — MARCS atmospheres, 2800–8000 K

[Korg.jl](https://github.com/ajwheeler/Korg.jl) synthesizes spectra from **MARCS**
model atmospheres, which it downloads as an artifact on first use. ROTIR uses it
as a grid *builder*: it emits a `RectGrid4` that `TabulatedProvider` evaluates,
so there is one evaluation path and one adjoint regardless of where the numbers
came from.

```julia
using Korg                      # loads ext/ROTIRKorgExt.jl
grid = build_korg_grid(Teff = range(5800, 7400, length = 9),
                       logg = range(3.2,  3.9,  length = 5),
                       λ    = (1.5e-6, 1.8e-6),
                       μ    = limb_mu_vec(12),
                       spherical = false)
save_intensity_grid("betcas_H_korg.fits", grid)
prov = TabulatedProvider(grid; name = "Korg/MARCS H band")
```

Four things measured against Korg 1.3.1 that its documentation does not state,
each of which silently produces a plausible wrong grid:

* `I_scheme = "linear"` is mandatory; the default ignores `mu_values` entirely.
* `result.intensity` is `(nμ, nλ, nlayers)` and **layer 1 is the surface**.
* Its first axis is rays ordered **inward then outward**, so the emergent
  intensity is the *last* `nμ` rows. Taking the first `nμ` gives mostly zeros
  that are still monotone in μ — and for a plane-parallel atmosphere the two
  coincide, so the mistake survives testing on dwarfs and corrupts every giant.
* Wavelengths are **vacuum**. Hα is 6562.79 Å in air and 6564.60 in vacuum —
  83 km/s, which is not a rounding detail for velocity-resolved work.

`spherical` is an explicit keyword and Korg's own `logg`-dependent default is
never used: a grid must use one geometry throughout, since a spherical model's
intensity genuinely vanishes below a μ_min that moves with (Teff, logg), and
interpolating across that discontinuity is meaningless. Plane-parallel is the
default because it is what is self-consistent with a tessellated Roche surface —
the mesh already carries the global geometry and each tessel is a local patch.
Use `spherical = true` for a single extended star (ρ Cas, RW Cep near log g 0.5).

## Kurucz ATLAS9 — 3500–50000 K

Korg's MARCS grid stops at 8000 K, which covers β Cas (7208 K) and nothing
hotter. Kurucz's ATLAS9 **intensity packs** run to 50000 K at 17 emergent
angles, which covers Spica, Vega, Regulus and β Lyr.

Download one pack from `http://kurucz.harvard.edu/grids/gridp00/` — nothing is
bundled:

```bash
curl -O http://kurucz.harvard.edu/grids/gridp00/ip00k2.pck19     # 67 MB
```

```julia
grid = read_kurucz_intensity("ip00k2.pck19";
                             λrange = (1.5e-6, 1.8e-6),
                             Teff = (18000.0, 28000.0), logg = (3.0, 4.5))
prov = TabulatedProvider(grid; name = "Kurucz ATLAS9 H band")
```

**Pick the right pack.** `ip00k2.pck19` (67 MB) reaches 50000 K.
`ip00k2new.pck` (31 MB) has newer ODFNEW opacities but tops out at **8750 K**,
so it is no better than MARCS for this purpose.

Two format traps, both of which yield a plausible wrong grid rather than an
error — the reader handles both, and they are documented here because the files
themselves announce neither:

1. **The two packs differ by one column.** Fortran's carriage-control convention
   makes column 1 a printer control character; `.pck19` keeps it (header
   `TEFF   3500.`) and `new.pck` has it stripped (`EFF   3500.`, every line
   shifted one left). The record is fixed-width `format(F9.2,1pe10.3,16I6)` and
   its fields abut at the extremes — a 160000.00 nm wavelength fills F9.2
   exactly — so whitespace splitting cannot recover them. The reader detects the
   shift from the header.
2. **The 16 integers are RATIOS**, `I(μ)/I(μ=1) × 10⁵`, not intensities. Reading
   them as intensities is uniformly 10⁵ too large, which cancels out of a
   flux-normalised visibility and so survives every relative-only test.

`λrange` is required and worth restricting hard: a full pack is 1221 wavelengths
× 17 angles × ~500 models, while an H-band grid needs perhaps 30 wavelengths.

A result worth knowing before choosing a target: **a hot star is *less*
limb-darkened in H than a warm one** — ``I(1)/I(0.01)`` is 1.44 at 25000 K
against 1.95 at 7000 K, because at 1.65 μm a 25000 K photosphere sits far down
the Planck tail where the source function varies slowly with depth. So *fitting*
limb darkening on a hot star is even less constrained than on β Cas, and
predicting it matters more there, not less.

## SATLAS — spherical ATLAS, available but not yet read

SATLAS (Neilson & Lester 2013) provides **spherically extended** ATLAS models
with per-μ intensities, for red giants (A&A 554, A98) and FGK dwarfs
(A&A 556, A86), Teff 3000–8000 K. They are distributed through CDS:

* `https://cdsarc.cds.unistra.fr/ftp/J/A+A/554/A98/` — giants, ``\log g < 4``
* `https://cdsarc.cds.unistra.fr/ftp/J/A+A/556/A86/` — dwarfs, ``\log g ≥ 4``

ROTIR has **no SATLAS reader yet**. Adding one belongs beside
`read_kurucz_intensity` and would produce the same `RectGrid4`. The one point
needing care is that a spherical model's tabulated μ runs over the *extended*
limb, so it has to be rescaled onto the Rosseland radius —
``μ' = (μ - μ_0)/(1 - μ_0)`` with ``μ_0 = \sqrt{1 - (R_\mathrm{Ross}/R_\mathrm{LD})^2}``
— before it can be pasted onto a tessellated figure, or the extension is counted
twice.

## TLUSTY — not currently obtainable

TLUSTY's NLTE grids (OSTAR2002, 27.5–55 kK; BSTAR2006, 15–30 kK) would be the
right choice for O and early-B stars. **The canonical host is gone**:
`nova.astro.umd.edu` redirects to `tlusty.oca.eu`, which serves an unrelated
page, and the Arizona mirror 404s. SVO and CDS mirror TLUSTY *fluxes*, not
per-μ intensities, which is what a limb-darkened surface needs.

It matters less than it sounds. Kurucz is LTE and LTE is adequate through
Spica's 25300 K; NLTE becomes important for O stars, which this package has no
target for. If the intensity files turn up locally, a reader belongs beside
`read_kurucz_intensity`.

## Choosing a backend

| | Teff | per-μ? | geometry | status |
|---|---|---|---|---|
| Korg / MARCS | 2800–8000 K | synthesized at any μ | either | built on demand |
| Kurucz ATLAS9 | 3500–50000 K | 17 fixed angles | plane-parallel | one download |
| SATLAS | 3000–8000 K | tabulated | spherical | no reader yet |
| TLUSTY | 15–55 kK | 20 angles | plane-parallel | host gone |

Build the grid over the range the **whole surface** spans, not the polar values:
on β Cas the surface runs Teff 5888–7206 K and ``\log g`` 3.29–3.75, a 0.466 dex
spread in gravity, which is why a scalar `logg` cannot serve this model.
