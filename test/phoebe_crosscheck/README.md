# PHOEBE cross-check for the radiosity solver

Validates `src/reflection.jl` against PHOEBE 2.5's own C++ implementation
(`libphoebe.mesh_radiosity_problem_nbody_convex`, from `phoebe/lib/reflection.h`) by
running **the same triangulation** through both and comparing element by element.

This is not part of `Pkg.test()` — it needs a Python environment with PHOEBE's compiled
extension. The Julia-side unit tests in `test/test_reflection.jl` (analytic limits,
energy conservation, Horvat ≡ Wilson) are the gate; this is the external confirmation.

## Result

Agreement to machine precision on 1536 + 384 triangles, albedo 0.6:

```
  uniform  Horvat  body1: max|Δ(reflected)|/max = 1.6e-14 | max Fout/F0  julia 1.010247172  phoebe 1.010247172
  uniform  Horvat  body2: max|Δ(reflected)|/max = 5.5e-15 | max Fout/F0  julia 1.135802378  phoebe 1.135802378
  uniform  Wilson  body1: max|Δ(reflected)|/max = 2.4e-13 | max Fout/F0  julia 1.010247172  phoebe 1.010247172
  uniform  Wilson  body2: max|Δ(reflected)|/max = 2.3e-14 | max Fout/F0  julia 1.135802378  phoebe 1.135802378
  linear   Horvat  body1: max|Δ(reflected)|/max = 1.6e-14 | max Fout/F0  julia 1.010230881  phoebe 1.010230881
  linear   Horvat  body2: max|Δ(reflected)|/max = 2.8e-15 | max Fout/F0  julia 1.135782746  phoebe 1.135782746
  linear   Wilson  body1: max|Δ(reflected)|/max = 1.8e-13 | max Fout/F0  julia 1.010284219  phoebe 1.010284219
  linear   Wilson  body2: max|Δ(reflected)|/max = 3.0e-14 | max Fout/F0  julia 1.135848157  phoebe 1.135848157
```

The **linear** rows are the important ones: a non-uniform bolometric law is what exercises
the `D0 = 2π∫₀¹D(μ)μdμ` normalisation in `ld_bol_D0`. With uniform LD the limb-darkened
and Lambert kernels coincide and the test is much weaker (it is also the case where Horvat
and Wilson must agree exactly, which they do).

## Running it

```bash
julia --project=. test/phoebe_crosscheck/export_geometry.jl     # writes data/ + Julia results
<phoebe-python> test/phoebe_crosscheck/compare_with_phoebe.py   # runs PHOEBE, prints the table
```

## Getting a working PHOEBE

PHOEBE 2.5 is not on PyPI in a form that builds cleanly everywhere. From a source clone:

```bash
python -m venv ~/SOFTWARE/phoebe_env
~/SOFTWARE/phoebe_env/bin/pip install numpy
cd ~/SOFTWARE/phoebe2
~/SOFTWARE/phoebe_env/bin/pip install --no-build-isolation --no-deps .
```

`--no-deps` is deliberate. On recent toolchains the `ndpolator` dependency fails to build
(`-Werror` plus a `_POSIX_C_SOURCE` redefinition against Python 3.14 headers), which aborts
the whole install even though PHOEBE itself compiles. `ndpolator` is only used by
`phoebe/atmospheres/passbands.py`, so `import phoebe` will fail — but `import libphoebe`
works, and `libphoebe` is where the radiosity solver lives. This script only uses
`libphoebe`.

## API notes (easy to get wrong)

* `model` and `support` must be **bytes**, not `str` — `libphoebe.cpp:107` defines
  `PyString_Type` as `PyBytes_Type`.
* Limb-darkening coefficients must be **numpy float64 arrays**, not Python lists:
  `LDmodelFromTuple` does an unchecked `PyArray_DATA((PyArrayObject*)...)` cast, so a list
  segfaults rather than raising.
* Connectivity must be C-contiguous `int32`; vertices, normals, areas, albedos and `F0`
  C-contiguous `float64`. `PyArray_To3DPointVector` reads raw memory with no dtype check.
* `support="triangles"` expects `N`, `A`, `R`, `F0` per **triangle** and `V` per vertex.
* PHOEBE's `gravb_bol` is 4× ROTIR's `beta` (`T = Tpole·((g/gpole)^gravb_bol)^0.25` versus
  `T = tpole·(g/gpole)^beta`). Not used by this script, but it bites when transcribing
  parameters.
