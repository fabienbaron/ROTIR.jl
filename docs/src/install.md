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
| PyCall / PyPlot | Matplotlib-based plotting |

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

## Verify

```julia
using ROTIR
```

This should load without errors and re-export key OITOOLS functions
(`readoifits`, `readoifits_multiepochs`, `readfits`, `writefits`).
