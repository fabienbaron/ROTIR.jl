# Kurucz ATLAS9

Korg's MARCS grid stops at 8000 K. Kurucz's ATLAS9 **intensity packs** reach 50000 K at 17
emergent angles, which is what covers the hot rotators. See
[Model atmospheres and limb darkening](@ref) for the interface these feed.

## Getting a pack

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
