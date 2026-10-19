# Installation

Julia 1.13 or later.

```julia
using Pkg
pkg"registry add General"
pkg"registry add https://github.com/emmt/EmmtRegistry"   # OptimPackNextGen
Pkg.add(url="https://github.com/fabienbaron/OITOOLS.jl.git")
Pkg.add(url="https://github.com/fabienbaron/ROTIR.jl.git")
```

`using ROTIR` should load quickly, without Makie, Qt or Python. It re-exports `readoifits`,
`readoifits_multiepochs`, `readfits` and `writefits` from OITOOLS.

## Dependencies

| Package | Use |
|---------|-----|
| [OITOOLS.jl](https://github.com/fabienbaron/OITOOLS.jl) | OIFITS reading, data handling, uv setup |
| OptimPackNextGen | VMLMB optimiser for reconstructions |
| NFFT, FFTW | transforms for the rasterised route |
| NLopt | derivative-free optimisers for parametric fits |
| FITSIO | FITS images and surface maps |

## Optional packages

Optional features are [package
extensions](https://pkgdocs.julialang.org/v1/creating-packages/#Conditional-loading-of-code-in-packages),
active once their trigger packages are loaded. Each must be installed in your own project.

| Load | Enables |
|------|---------|
| `Makie` | `plot2d_makie`, `plot3d_makie`, `plot_mollweide_makie`, the decorations |
| `GLMakie, QMLMakie, QML` | the GUI: [`rotirgui`](@ref), `gui` |
| `PythonPlot` | the matplotlib plots (`plot2d`, `plot3d`, `plot_mollweide`, …) |
| `Zygote` | `fit_parametric`, `bootstrap_parametric`, the GUI's gradient fit |
| `AdvancedHMC, LogDensityProblems` (with `Zygote`) | `method = :nuts` |
| `Nautilus` | `method = :nautilus` (nested sampling, with evidence) |
| `Pigeons, Distributions, ADTypes` (with `Zygote`) | `method = :pigeons` (parallel tempering) |
| `PythonCall` | `method = :ultranest` |
| `FINUFFT` | `POLYFT_BACKEND[] = :nufft` |
| `LoopVectorization` | `POLYFT_BACKEND[] = :turbo` (vectorised exact kernel) |
| `Korg` | MARCS model-atmosphere intensities ([Korg and MARCS](@ref)) |

The default polygon transform, `:t3`, needs neither FINUFFT nor LoopVectorization.

Calling a function whose extension is not loaded raises a `MethodError`.
`nautilus_available()`, `hmc_available()`, `pigeons_available()`, `ultranest_available()` and
`turbo_available()` report whether an extension is loaded.

Loading PythonCall, Pigeons or LoopVectorization invalidates precompiled plotting code and
adds about 2 s each to GUI startup.

## The graphical interface

```julia
using ROTIR
rotirgui()                             # empty session
rotirgui("mydata.oifits")              # one dataset
rotirgui("ep1.oifits", "ep2.oifits")   # several epochs
```

If GLMakie, QMLMakie or QML are missing, `rotirgui` lists them with the optional samplers and
asks what to install; Enter installs only the required ones. `install = true` installs
everything, `install = :required` only the required packages, `install = false` stops with
instructions. `samplers = false` does not load the optional samplers: startup is faster and
their fit methods are not offered.

With the packages already loaded:

```julia
using ROTIR, GLMakie, QMLMakie, QML
gui()
```

From a clone, with the pinned launcher environment:

```
julia --project=bin bin/rotirgui.jl [file.oifits …]
```

## Python

Python is needed only for PythonPlot and `method = :ultranest`; the GUI does not load it
(`method = :nautilus` is the pure-Julia alternative). PythonCall installs a private
environment through CondaPkg on first use, with matplotlib, `ultranest` and `astroquery`.

- To use an existing interpreter, set `JULIA_PYTHONCALL_EXE` to its path and
  `JULIA_CONDAPKG_BACKEND=Null`; that environment must provide matplotlib.
- numpy 2 or later is required. With numpy 1.x, PythonCall's array conversion fails with
  `ValueError: NoneType copy mode not allowed` on some arrays, which shows up as plots failing.

## Reference

```@docs
rotirgui
```
