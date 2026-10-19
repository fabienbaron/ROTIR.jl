# Model atmospheres and limb darkening

A limb-darkening law (linear, quadratic, Hestroffer, Claret) is a fitted approximation. A
model atmosphere gives the specific intensity ``I(T_\mathrm{eff}, \log g, \mu, \lambda)``
directly: the limb darkening is predicted, and its coefficients are not fitted.

Gravity darkening and limb darkening both remove flux from the limb, and an interferometer
measures their sum. On β Cas, with both free, the correlation between the gravity-darkening
exponent β and `ld1` is +0.67 to +0.99 depending on rotation rate. With the intensity
predicted, `ld1` is not a parameter and the correlation is gone.

## From temperature to intensity

ROTIR's maps are temperatures. The conversion to intensity is set once per model and applies
to any temperature map (von Zeipel, ELR, spotted, reconstructed) and any limb-darkening law.
In the GUI it is the **Intensity from** row of the Model form; in a script, the
`intensity_source` field of [`INTENSITY_PARAMS`](@ref):

| code | source | intensity | μ dependence | reads |
|---|---|---|---|---|
| 5 | Rayleigh–Jeans (default) | ``I = \lambda T / c_2``, linear in ``T`` | from the LD law | nothing |
| 1 | Planck | ``1/(e^{c_2/\lambda T} - 1)`` | from the LD law | nothing |
| 2 | ROTIR grid | a grid shipped with ROTIR, or a custom FITS grid | its own | `grid_choice` (+ `grid_path` if custom) |
| 3 | Kurucz pack | the ATLAS9 intensity pack, downloaded once | its own | λ, Teff, log g windows |
| 4 | Korg (MARCS) | synthesized, Teff ≤ 8000 K, needs `using Korg` | its own | windows + `grid_nodes` |

Rayleigh–Jeans is `intensity_model = :linear`: intensity proportional to the map. This is the
conversion made by imaging codes that fit a brightness per patch, such as SURFING, so
reproducing their maps needs it. Planck and Rayleigh–Jeans share a normalisation and agree
when ``c_2/\lambda T \ll 1``. Both drop the ``\lambda^{-5}`` (``\lambda^{-4}``) prefactor,
which is exact when every component in a χ² is at the same wavelength.

The three atmosphere sources carry their own μ dependence and require `ldtype = 0`. They also
need a per-tessel ``\log g``, hence a distance `d` (next section).

The shipped grids are listed in [`SHIPPED_INTENSITY_GRIDS`](@ref) and selected by
`grid_choice` (0 is a custom file):

| `grid_choice` | file | covers |
|---|---|---|
| 1 | `spica_H_kurucz.fits` | Kurucz ATLAS9, H band, Teff 17–30 kK, log g 3.0–4.5 |
| 2 | `betcas_H_korg.fits` | Korg/MARCS, H band, Teff 5.4–7.6 kK, log g 2.9–4.1 |
| 3 | `betcas_Br11_korg.fits` | Korg/MARCS, Br-11 window, Teff 5.4–7.6 kK, log g 2.9–4.1 |

The Kurucz pack (67.5 MB) is downloaded once and kept by ROTIR; see
[`kurucz_pack_path`](@ref). The form shows only the rows the chosen source reads
([`intensity_params_used`](@ref)); an analytic source shows none.

From a script:

```julia
cfg = Dict(:intensity_source => 2.0, :grid_choice => 2.0)       # shipped Korg H-band grid
isempty(validate_intensity_source(0, cfg, "")) || error("bad intensity configuration")
prov = ROTIR.build_intensity_provider(cfg)
```

The Imaging tab uses the model's **Intensity from**: Rayleigh–Jeans reconstructs with
`intensity_model = :linear`, Planck with `:planck` at the data's mean wavelength. The
atmosphere sources are refused there, because `image_reconstruct_oi` takes an
`intensity_model` symbol, not an `IntensityProvider`. A sky-image export uses the same
conversion at each channel (see [Exporting sky images for OITOOLS](@ref)).

## Where `log g` comes from: the mass is derived

An atmosphere is indexed by ``\log g``, which needs a mass. ROTIR does not fit one. For a
Roche rotator, ``\Omega = f_\mathrm{ev}\sqrt{8GM/27R_p^3}``, so with ``\Omega = 2\pi/P``
from the rotation period:

```math
M = \frac{27\,R_p^3\,\Omega^2}{8\,G\,f_\mathrm{ev}^2}
```

`rpole` is an angular radius, so the one added parameter is the distance `d` in parsecs. The
interferometric data constrain `d` to about 12 pc through ``\log g``; a parallax gives about
0.1 pc, so `d` is fixed from the parallax rather than fitted.

```julia
p = (surface_type = 2, rpole = 0.849, d = 16.8, frac_escapevel = 0.92,
     rotation_period = 0.893, tpole = 7208.0, beta = 0.25, inclination = 20.0,
     position_angle = 115.0, gravity_law = 2, ldtype = 0,
     ld1 = 0.0, ld2 = 0.0, ld3 = 0.0, ld4 = 0.0)

derived_quantities(p)      # mass, R_p, R_eq, logg_pole, v_eq, vsini
println(derived_summary_text(p))
```

For β Cas these give M = 1.936 M⊙, R_p = 3.067 R⊙, R_eq = 3.796 R⊙,
``\log g_\mathrm{pole}`` = 3.752 and vsini = 73.2 km/s, consistent with the literature for an
F2III-IV star. vsini is not an input; it follows from the geometry.

A Fisher covariance propagated through `derived_quantities` at that solution gives
M = 1.936 ± 0.042 M⊙ (2.1 %), the polar radius to 0.63 % and the equatorial radius to 0.60 %.

* ``\log g`` varies over the surface: 0.466 dex from pole to equator on β Cas. The grid must
  cover the range the whole surface spans. [`logg_map`](@ref) returns the per-tessel values
  and has its own adjoint, so `d` is differentiable.
* `rotation_period`, `frac_escapevel`, `rpole` and `d` over-determine the mass, and can
  imply an implausible one. `advise_star_params` reports it; optimisers do not check it, and a
  fit can reach 7 M⊙. For a production fit, pass a `logprior` that penalises an implausible
  mass, or free `rotation_period`.

Without a distance, `has_physical_scale(p)` is false, `derived_quantities` returns `nothing`,
`derived_summary_text` returns an empty string, and atmosphere providers cannot be used. A
`PlanckProvider` can: it ignores ``\log g``, and ``\partial\log\pi/\partial d`` is then zero
(checked by the test suite).

## The provider interface

Each source is an [`IntensityProvider`](@ref), used the same way by the fitting and sampling
code:

```julia
I = provider_intensity(provider, Teff, logg, μ, λ)
```

`owns_mu(provider)` is true when the provider supplies the μ dependence; the star then needs
`ldtype = 0`, otherwise the limb darkening is applied twice. `check_provider_consistency`
reports that combination, and `build_parametric_logπ` refuses it.

A provider needs a per-tessel `logg`, hence a mass, which ROTIR derives from
`(rpole, d, frac_escapevel, rotation_period)`; see [`derive_mass`](@ref).

## Choosing a backend

| | Teff | μ dependence | geometry |
|---|---|---|---|
| Korg / MARCS | 2800–8000 K | synthesized at any μ | plane-parallel or spherical |
| Kurucz ATLAS9 | 3500–50000 K | 17 fixed angles | plane-parallel |
| Rayleigh–Jeans / Planck | any | from the LD law | — |

Build the grid over the range of the whole surface: on β Cas, Teff 5888–7206 K and
``\log g`` 3.29–3.75.

Each backend has its own page: [Korg and MARCS](@ref) and [Kurucz ATLAS9](@ref).
