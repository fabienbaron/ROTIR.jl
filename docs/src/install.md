# Installation

## Dependencies

| Package | Purpose |
|---------|---------|
| [OITOOLS.jl](https://github.com/fabienbaron/OITOOLS.jl) | OIFITS I/O, data handling, UV frequency setup |
| OptimPackNextGen | VMLMB gradient optimizer for image reconstruction |
| LinearAlgebra | Matrix operations, dot products |
| SparseArrays | Sparse total variation matrices |
| Statistics | Mean, standard deviation |
| FITSIO | Reading/writing FITS images |
| PythonCall / PythonPlot | Matplotlib-based plotting |

## Step 1: Julia registry

ROTIR depends on OptimPackNextGen from the EmmtRegistry:

```julia
using Pkg
pkg"registry add General"
pkg"registry add https://github.com/emmt/EmmtRegistry"
```

## Step 2: Install OITOOLS

```julia
Pkg.add(url="https://github.com/fabienbaron/OITOOLS.jl.git")
```

## Step 3: Install ROTIR

```julia
Pkg.add(url="https://github.com/fabienbaron/ROTIR.jl.git")
```

## Python / matplotlib

ROTIR uses PythonPlot.jl for plotting, which needs a Python environment with
matplotlib. You normally do not have to do anything: PythonCall manages a private
environment through CondaPkg, and installs matplotlib (plus `ultranest` and
`astroquery`, declared by OITOOLS) on first use.

Two notes:

* To point PythonCall at an interpreter you already have, set `JULIA_PYTHONCALL_EXE` to its
  absolute path and `JULIA_CONDAPKG_BACKEND=Null`; that environment must then supply
  matplotlib itself.
* OITOOLS ≥ 0.11.1 requires **numpy ≥ 2**. PythonCall's `__array__` passes
  `copy=None`, the numpy-2 spelling of "copy only if needed", which numpy 1.x rejects
  with `ValueError: NoneType copy mode not allowed`. The failure is data-dependent —
  dense arrays convert fine while nested vertex lists do not — so an old numpy tends to
  surface as plots failing rather than as an import error.

## Verify

```julia
using ROTIR
```

This should load without errors and re-export key OITOOLS functions
(`readoifits`, `readoifits_multiepochs`, `readfits`, `writefits`).
