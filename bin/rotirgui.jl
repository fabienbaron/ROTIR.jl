#!/usr/bin/env julia
#
# Launcher.
#
#     julia --project=bin bin/rotirgui.jl [file.oifits]
#
# Use --project=bin, not --project=. — bin/Project.toml is where GLMakie, QMLMakie and QML are
# real dependencies, so bin/Manifest.toml pins them. Under --project=. those three are WEAK
# dependencies and are not loadable at all; if they resolve anyway it is from your default
# environment, at whatever versions happen to be there.
#
# The ordering below is the whole point of this file, and one fact forces it: Mesa and GLFW both
# read their configuration when the first OpenGL context is created, so both hints have to be
# set above `using GLMakie`. That is why neither call lives in the GUI extension — loading
# GLMakie is what creates the extension, so anything inside it is already too late.
#
# `configure_graphics!` and `configure_qt_platform!` come from OITOOLS, which ROTIR already
# depends on. They are not reimplemented here: the questions they answer (which Mesa driver,
# which Qt platform plugin, do GLFW and Qt agree about X11 vs Wayland) have nothing to do with
# which of the two packages is asking.
#
# MPLBACKEND does not need forcing: ROTIR keeps matplotlib in an extension
# (ROTIRPythonPlotExt), so `using ROTIR` starts no backend probe and maps no second Qt.
# `check_qt_conflict` runs inside `gui()` in case something else in the session loaded one.

# A CondaPkg resolve that gets interrupted leaves its lock file behind, and every later start
# then blocks on it FOREVER — the process stays alive, no window appears, and the only clue is
# one CondaPkg info line, which is indistinguishable from a graphics hang. So: notice it, say
# what it is, and say how to clear it.
let lock = joinpath(@__DIR__, "..", ".CondaPkg", "lock")
    if isfile(lock)
        age = round(Int, time() - mtime(lock))
        @warn """
            A CondaPkg lock file exists and startup will block on it until it is released.

                $lock   (age $(age)s)

            If no other Julia process is resolving a Conda environment right now, it is stale —
            delete it and start again:

                rm -f $lock
            """
    end
end

using ROTIR
using OITOOLS: configure_graphics!, configure_qt_platform!, prefer_native_wayland!

configure_graphics!()          # before the first GL context

using GLFW_jll                 # activates OITOOLSGLFWExt; dlopens libglfw, needs no display

# GLMakie's platform is decided FIRST, then Qt is made to match it. The two must agree — one on
# Wayland and one on X11 means two EGL display connections in one process — and only GLMakie's
# choice can fail: GLFW.jl HARD-CODES X11 on Linux, so under a Wayland compositor `using
# GLMakie` silently lands on XWayland. `prefer_native_wayland!` initialises GLFW with the
# Wayland hint first; `glfwInit` is idempotent, so GLFW.jl's later `Init()` becomes a no-op and
# its X11 hint is irrelevant. **This must run BEFORE `using GLMakie`** — which is why it lives
# in OITOOLS' GLFW_jll extension rather than in a GUI extension that GLMakie itself triggers.
#
# Without it the launcher reports "GLMakie is on X11, so Qt is pinned to xcb", i.e. the whole
# application on XWayland — which works, but is the arrangement whose popup surface lifetimes
# the GL-context bug lives in.
wl = prefer_native_wayland!()
configure_qt_platform!(; match_x11 = !wl.applied)

using GLMakie

using QMLMakie, QML
# Nautilus is optional; without it the GUI simply does not offer that fit method.
# Loading it here rather than lazily keeps the method list stable for the whole
# session, since the list is read once when the panel first refreshes.
# Zygote gives the Model panel its fast gradient fit (`fit_parametric`); without it that method
# is simply not offered. Loaded here rather than lazily so the method list is stable for the
# whole session — the panel reads it once, when it first refreshes.
try
    @eval using Zygote
catch err
    @info "Zygote not available; the gradient fit method will not be offered" exception = err
end
# AdvancedHMC + LogDensityProblems + Zygote give the Model panel its NUTS method. Without them
# it is simply not offered.
try
    @eval using AdvancedHMC, LogDensityProblems
catch err
    @info "AdvancedHMC not available; the NUTS fit method will not be offered" exception = err
end
try
    @eval using Nautilus
catch err
    @info "Nautilus not available; the :nautilus fit method will not be offered" exception = err
end
# Pigeons gives the Model panel `:pigeons`, the sampler for a MULTIMODAL posterior — the shape
# a real ROTIR χ² has, and the one no other method here can see: NUTS samples the basin it
# started in and reports a confident interval that says nothing about the others.
#
# IT IS NOT FREE, and the cost is at startup rather than at use. MEASURED: `using Pigeons`
# invalidates the plot-construction code OITOOLS precompiles for its live canvas, and one
# `build_canvas` goes from ~370 ms to 3166 ms — about 2.8 s onto every session, sampler or no
# sampler. Nothing else loaded here does that; Zygote, AdvancedHMC, LogDensityProblems and
# Nautilus were each measured at 370-380 ms. It is loaded anyway because a fit method that is
# listed but unusable is worse than a slower start, and because the panel reads the method list
# once when it first refreshes — a lazy import would not appear in it.
try
    @eval using Pigeons, Distributions, LogDensityProblems, ADTypes
catch err
    @info "Pigeons not available; the :pigeons fit method will not be offered" exception = err
end

# PYTHON IS NOT LOADED, AND THE GUI IS PYTHON-FREE BY CONSTRUCTION. `:ultranest` is not offered
# at all — not gated, not listed — so nothing here can pull PythonCall in. `:nautilus` covers
# the same ground (importance nested sampling, an evidence, uncertainties) in pure Julia and
# agrees with UltraNest to 1.4e-6 mas on a λ And sphere fit, so nothing is given up. UltraNest
# is still there for a script: `using PythonCall` loads ROTIRUltraNestExt and
# `fit_sphere_ld(...; method = :ultranest)` works as it always did.
#
# Beyond the sampler, matplotlib's backend probe maps a second Qt into this process, and QML.jl
# needs its own Qt to be the only one — so Python in a session with a window is not merely
# slow. PythonCall on its own also costs 1.2 s here for the same invalidation reason as
# Pigeons: 338 ms to 2477 ms for one canvas build, measured with nothing else changed.

# The GUI's own types live in ROTIRGUIExt, so they are reached after the fact.
const GUI = Base.get_extension(ROTIR, :ROTIRGUIExt)
GUI === nothing && error("ROTIRGUIExt did not load — check that GLMakie, QML and QMLMakie " *
                         "are all present in this environment")

session = GUI.Session()
for arg in ARGS
    isfile(arg) || (@warn "not a file, skipping" path = arg; continue)
    GUI.load_dataset!(session, arg)
end

gui(session)
