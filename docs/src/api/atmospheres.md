# Model Atmospheres, Stellar Physics & Spectra

Intensity providers, the derived stellar quantities they use, the surface velocity field and
the velocity-resolved observables. See the guide
[Model atmospheres and limb darkening](@ref).

## Derived stellar quantities

The mass is derived from `(rpole, d, frac_escapevel, rotation_period)`; with it, each
tessel has an absolute surface gravity.

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

## Temperature to intensity: choosing a source

The GUI's "Intensity from" setting; see [From temperature to intensity](@ref).

```@docs
INTENSITY_PARAMS
INTENSITY_SOURCES
SHIPPED_INTENSITY_GRIDS
ROTIR.build_intensity_provider
validate_intensity_source
intensity_params_used
intensity_spec
intensity_grid_path
intensity_lambda_range
kurucz_pack_path
intensity_cache_dir
```

## Intensity providers

```@docs
IntensityProvider
RayleighJeansProvider
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

`build_korg_grid` and `korg_provider` need `using Korg`. The Kurucz readers need the
intensity pack ([`kurucz_pack_path`](@ref) downloads it).

```@docs
build_korg_grid
korg_provider
read_kurucz_models
read_kurucz_intensity
KuruczModel
```

## Velocity field

```@docs
los_velocity
los_velocity_from_proj
doppler_lambda
rest_lambda
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
