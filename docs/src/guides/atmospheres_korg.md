# Korg and MARCS

[Korg.jl](https://github.com/ajwheeler/Korg.jl) synthesizes spectra from MARCS model
atmospheres, which it downloads on first use. ROTIR uses it to build a `RectGrid4`, which
`TabulatedProvider` evaluates like any other grid. See
[Model atmospheres and limb darkening](@ref) for the provider interface.

## Building a grid

```julia
using Korg                      # loads ext/ROTIRKorgExt.jl
grid = build_korg_grid(Teff = range(5800, 7400, length = 9),
                       logg = range(3.2,  3.9,  length = 5),
                       λ    = (15000.0, 18000.0),      # Å
                       μ    = (1e-3, 1.0, 12),         # (lo, hi, n) or a vector
                       spherical = false)
save_intensity_grid("betcas_H_korg.fits", grid)
prov = TabulatedProvider(grid; name = "Korg/MARCS H band")
```

`Teff`, `logg` and `μ` are each `(lo, hi, n)` or a vector; μ must be increasing and start
above 0. `λ` is `(start, stop)` in Å, with `λ_step` (default 0.05 Å). The grid is returned
with λ in metres and intensity in erg/s/cm²/Å/sr. MARCS covers Teff 2800–8000 K and
log g −0.5 to 5.5.

Korg works in vacuum wavelengths: Hα is 6562.79 Å in air and 6564.60 Å in vacuum, a
difference of 83 km/s.

`build_korg_grid` handles three behaviours of Korg 1.3.1 that otherwise give a wrong grid
without an error: it passes `I_scheme = "linear"` (the default ignores `mu_values`), takes
layer 1 of `result.intensity` as the surface, and takes the emergent intensity from the last
`nμ` rays.

## `spherical`

`spherical` is always explicit; Korg's own `logg`-dependent switch is not used, so a grid has
one geometry throughout. A spherical model's intensity is zero below a μ_min that depends on
(Teff, logg), which cannot be interpolated across a mixed grid.

- `spherical = false` (default): plane-parallel, for a tessellated surface, where each tessel
  is a local patch and the mesh carries the global shape.
- `spherical = true`: for a single extended star (ρ Cas, RW Cep near log g 0.5). Keep the
  `logg` axis fine enough to resolve the moving μ_min.
