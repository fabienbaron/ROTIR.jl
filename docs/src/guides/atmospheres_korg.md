# Korg and MARCS

[Korg.jl](https://github.com/ajwheeler/Korg.jl) synthesizes spectra from **MARCS** model
atmospheres. See [Model atmospheres and limb darkening](@ref) for the interface these feed.

## Building a grid

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
