# Kurucz ATLAS9

Kurucz's ATLAS9 intensity packs reach 50000 K at 17 emergent angles, which covers the hot
rotators (Spica, Vega, Regulus, β Lyr); Korg/MARCS stops at 8000 K. See
[Model atmospheres and limb darkening](@ref) for the provider interface.

## Getting a pack

[`kurucz_pack_path`](@ref) with `download = true` downloads `ip00k2.pck19` (67.5 MB) once from
`http://kurucz.harvard.edu/grids/gridp00/` into `intensity_cache_dir()` (under the first Julia
depot) and returns the cached path afterwards. Choosing "Kurucz pack" in the GUI does the same.
Set `ROTIR_KURUCZ_PACK` to use an existing copy.

```julia
grid = read_kurucz_intensity(kurucz_pack_path(; download = true);
                             λrange = (1.5e-6, 1.8e-6),          # metres
                             Teff = (18000.0, 28000.0), logg = (3.0, 4.5))
prov = TabulatedProvider(grid; name = "Kurucz ATLAS9 H band")
```

`λrange` is required. A full pack is 1221 wavelengths × 17 angles × ~500 models; an H-band
grid needs about 30 wavelengths. `Teff` and `logg` restrict the axes.

The pack is not a complete (Teff, logg) rectangle: missing nodes are filled from the nearest
`logg` at the same Teff, with a warning giving the count. Restricting `logg` to what the
surface needs avoids the filled nodes.

`ip00k2new.pck` (31.5 MB) has the newer ODFNEW opacities but stops at 8750 K. The reader reads
both formats (they differ by one column) and converts the stored values, which are
`I(μ)/I(1) × 10⁵`, to intensities.

## Limb darkening of hot stars

A hot star is less limb-darkened in H than a warm one: ``I(1)/I(0.01)`` is 1.44 at 25000 K
and 1.95 at 7000 K. At 1.65 μm a 25000 K photosphere is far down the Planck tail, where the
source function varies slowly with depth. Fitted limb darkening is therefore even less
constrained on a hot star than on β Cas.
