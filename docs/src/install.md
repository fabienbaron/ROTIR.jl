# Installation

## Dependencies

These are always installed with ROTIR:

| Package | Purpose |
|---------|---------|
| [OITOOLS.jl](https://github.com/fabienbaron/OITOOLS.jl) | OIFITS I/O, data handling, UV frequency setup |
| OptimPackNextGen | VMLMB gradient optimizer for image reconstruction |
| FINUFFT | Type-3 non-uniform FFT — the default polygon-transform backend |
| NFFT, FFTW | Gridded transforms for the rasterised route |
| NLopt | Derivative-free optimisers for parametric fits |
| FITSIO | Reading and writing FITS images and surface maps |
| LinearAlgebra, SparseArrays, Statistics | Matrix operations, sparse regularisers, summaries |

Plotting, the graphical interface, the samplers and the vectorised kernel are **not**
here — they are weak dependencies, described next.

## Weak dependencies and extensions

`using ROTIR` loads no plotting toolkit, no Qt, no sampler and no Python. Everything
optional lives behind a [package
extension](https://pkgdocs.julialang.org/v1/creating-packages/#Conditional-loading-of-code-in-packages):
the code ships inside ROTIR but only becomes active once you load its trigger package
yourself.

This is not tidiness. It is measured: `using PythonCall` alone takes one OITOOLS canvas
build from 338 ms to 2477 ms by invalidating precompiled plotting code, and `using
Pigeons` and `using LoopVectorization` cost 2.8 s and 2.3 s of GUI startup the same way.
Anything a script does not need should not be loaded by it.

| Load this | and you get |
|-----------|-------------|
| `Makie` | `plot2d_makie`, `plot3d_makie`, `plot_mollweide_makie` and the decorations |
| `GLMakie, QMLMakie, QML` | the graphical interface — [`rotirgui`](@ref) and `gui` |
| `PythonPlot` | the matplotlib plotting layer (`plot2d`, `plot3d`, `plot_mollweide`, …) |
| `Zygote` | `fit_parametric`, `bootstrap_parametric`, and the gradient fit in the GUI |
| `AdvancedHMC, LogDensityProblems` (with `Zygote`) | `method = :nuts` |
| `Nautilus` | `method = :nautilus` — nested sampling, with an evidence |
| `Pigeons, Distributions, ADTypes` (with `Zygote`) | `method = :pigeons` — parallel tempering |
| `PythonCall` | `method = :ultranest` |
| `LoopVectorization` | `POLYFT_BACKEND[] = :turbo`, the vectorised exact kernel |

Each of these must be installed in **your own** project — a weak dependency of ROTIR is
not automatically available to you:

```julia
using Pkg; Pkg.add(["GLMakie", "QMLMakie", "QML"])
```

Calling something whose extension is not loaded raises a `MethodError` naming the
function rather than an `UndefVarError`, because the names are declared in ROTIR itself.
Several also have a predicate you can ask first: `nautilus_available()`,
`hmc_available()`, `pigeons_available()`, `ultranest_available()`, `turbo_available()`.

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

## The graphical interface

Three tabs over one session — data exploring, model fitting and imaging. From a plain
Julia session:

```julia
using ROTIR
rotirgui()                             # empty session
rotirgui("mydata.oifits")              # with a dataset loaded
rotirgui("ep1.oifits", "ep2.oifits")   # several epochs
```

[`rotirgui`](@ref) does the whole startup itself: it sets the Mesa, GLFW and Qt hints
(which must be applied *before* the first OpenGL context exists, so no extension can do
it), loads GLMakie, then QMLMakie and QML, then the optional sampler stack, and opens the
window.

If GLMakie, QMLMakie or QML are missing it offers to install them, and lists the optional
samplers at the same time — pressing **Enter** installs only what the window needs and
goes straight to the GUI:

```
rotirgui needs GLMakie.
  These will be added to …/v1.12/Project.toml, and precompiled — a few minutes, once.

  Optional. Each adds one fit method to the Model panel, and the panel reads
  that list when it opens — so a method not installed now is not offered for
  the rest of the session. All pure Julia; nothing here loads Python.

   1  Zygote                             gradient fit — analytic gradients
   2  AdvancedHMC + LogDensityProblems   :nuts — a posterior (wants Zygote too)
   3  Nautilus                           :nautilus — a posterior and an evidence
   4  Pigeons + Distributions + ADTypes  :pigeons — multimodal posteriors; ~2.8 s per start

  ENTER  install GLMakie only, then open the window
  a      install the optional ones as well
  1 3    install only those, by number
  n      cancel
```

`install = true` skips the question and takes everything, `install = :required` takes only
the toolkit, and `install = false` explains and stops. `samplers = false` starts faster at
the cost of those fit methods not being offered.

To open the window with a stack you have loaded yourself — which is what
[`rotirgui`](@ref) ultimately calls — use `gui`:

```julia
using ROTIR, GLMakie, QMLMakie, QML
gui()
```

From a clone there is also a pinned launcher environment, which is the reproducible way to
run it:

```
julia --project=bin bin/rotirgui.jl [file.oifits …]
```

## Python / matplotlib

Python is entirely optional, and the GUI never loads it: `:ultranest` is not offered
there, and matplotlib's backend probe would map a second Qt into a process where QML.jl
needs its own to be the only one. `method = :nautilus` covers the same ground in pure
Julia.

For the matplotlib plotting layer (`using PythonPlot`) or `method = :ultranest` (`using
PythonCall`) you need a Python environment with matplotlib. You normally do not have to do
anything: PythonCall manages a private environment through CondaPkg, and installs
matplotlib (plus `ultranest` and `astroquery`, declared by OITOOLS) on first use.

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

That alone should be quick and should pull in neither Makie nor Qt nor Python. To check the
window as well:

```julia
using ROTIR; rotirgui()
```

This should load without errors and re-export key OITOOLS functions
(`readoifits`, `readoifits_multiepochs`, `readfits`, `writefits`).

## Reference

```@docs
rotirgui
```
