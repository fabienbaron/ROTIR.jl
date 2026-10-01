# Model Atmospheres, Stellar Physics & Spectra

API for the intensity providers, the derived stellar quantities that index them,
the surface velocity field and the velocity-resolved observables built on it.
See the [Model Atmospheres](@ref "Model atmospheres and limb darkening") guide
for how these fit together.

## Derived stellar quantities

The mass is **derived** from `(rpole, d, frac_escapevel, rotation_period)` rather
than fitted, which is what makes a per-tessel surface gravity available at all.

```@docs
derive_mass
derived_quantities
derived_summary_text
has_physical_scale
logg_map
logg_pole
equatorial_velocity
projected_veq
```

## Intensity providers

```@docs
IntensityProvider
PlanckProvider
TabulatedProvider
RectGrid4
owns_mu
needs_logg
provider_support
check_provider_consistency
provider_map
analytic_test_grid
limb_mu_vec
```

## Grid persistence

```@docs
save_intensity_grid
load_intensity_grid
```

## Backends

`build_korg_grid` lives in a package extension and needs `using Korg`. The
Kurucz readers need a downloaded intensity pack; nothing is bundled.

```@docs
build_korg_grid
read_kurucz_models
read_kurucz_intensity
KuruczModel
```

## Velocity field

```@docs
los_velocity
los_velocity_from_proj
doppler_lambda
velocity_field_summary
```

## Velocity-resolved observables

```@docs
surface_state
channel_intensity
spectral_cvis
line_profile
differential_observables
continuum_mask
normalize_profile
line_equivalent_width
```
