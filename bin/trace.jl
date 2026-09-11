# What `app/build.jl` traces to build the application sysimage.
#
#     xvfb-run -a julia --project=bin bin/trace.jl        # standalone, to check it runs
#
# PackageCompiler records every method this script compiles and bakes those specialisations
# into the image, so what belongs here is exactly the work a user waits for on a cold start.
#
# THIS TRACE DOES NOT OPEN THE WINDOW, and that is deliberate — the same choice OITOOLS makes.
#
# It used to. ROTIR's `gui` takes `autoquit_ms`, so unlike OITOOLS' event loop it does return,
# and tracing the real startup path looked like a free win. It is not: the window runs Makie on
# Qt's own render thread, and when that thread dies the whole build dies with it. MEASURED —
# `QSGRenderThread::syncAndRender -> JuliaRenderer::render -> QMLMakie renderfunction` took
# SIGSEGV under `--compile=all`, after twenty minutes, with the bundle already half written. A
# segfault in another thread cannot be caught, so there is no way to make it survivable here.
#
# What the window would have contributed — `loadqml`, the QML engine, the Qt/GLMakie bridge and
# the first frame — comes from `app/precompile_statements.jl` instead, which is generated from
# a REAL session (see its header) and therefore covers those paths better than one scripted
# launch ever did.
#
# NEEDS A DISPLAY. Xvfb is enough; `app/build.jl` says so and does not supply one itself.
#
# NO PYTHON. PythonCall and PythonPlot are weak dependencies of ROTIR and are not in `app/`:
# an application that dragged in a conda environment would map a second Qt into the process
# and defeat the point. Nested sampling therefore comes from Nautilus, not UltraNest.

using ROTIR
import OITOOLS
using OITOOLS: configure_graphics!, configure_qt_platform!, prefer_native_wayland!

# The same order as bin/rotirgui.jl, and for the same reason: Mesa and GLFW read their
# configuration when the first GL context is created, so both hints must be set above
# `using GLMakie`. See that file's header for why neither can live in the GUI extension.
configure_graphics!()
using GLFW_jll
wl = prefer_native_wayland!()
configure_qt_platform!(; match_x11 = !wl.applied)

using GLMakie, QMLMakie, QML

# Optional engines, loaded rather than exercised: `create_app` bundles them as ordinary
# dependencies of `app/`, and the fit-method list is read once when a panel first refreshes,
# so loading them here keeps that list the same as the application's.
for pkg in (:Nautilus, :Zygote)
    try
        @eval using $pkg
    catch err
        @debug "not in this environment; not traced" pkg err
    end
end

const GUI = Base.get_extension(ROTIR, :ROTIRGUIExt)
const MK  = Base.get_extension(ROTIR, :ROTIRMakieExt)
(GUI === nothing || MK === nothing) && error("extensions did not load; nothing to trace")

const FILES = filter(!isnothing,
                     [ROTIR.resource("demos", "data", "2011Sep02.lam_And_prepped.oifits"),
                      ROTIR.resource("demos", "data", "2011Sep06.lam_And_prepped.oifits")])
isempty(FILES) && error("no demo data to trace against")
const OUT = joinpath(mktempdir(), "trace.png")

"""
    traced(f, what)

Run one section, reporting rather than throwing if it fails.

A bundle build is tens of minutes and everything before a failure is discarded, so one unlucky
step must not destroy it. The cost of catching is a thinner image, so a failure is reported
loudly: a warning here means whatever `what` names gets compiled in the user's first session
instead.
"""
function traced(f, what)
    try
        f()
    catch err
        @warn "trace section failed; its methods will NOT be in the image" section = what exception = err
    end
    return nothing
end

# ── the numerics, before any window ──────────────────────────────────────────
#
# These are what a fit and a reconstruction spend their first call compiling. Deterministic,
# so they are not wrapped: a failure here is a real fault and should stop the build.

# `readoifits_multiepochs` returns nwav x nepochs; ROTIR works one spectral bin at a time,
# so row 1 and every column — the same slice `load_dataset!` takes.
data = collect(OITOOLS.readoifits_multiepochs(FILES; warn = false, verbose = false)[1, :])
tess = tessellation_healpix(3)
p    = default_star_params(2; rpole = 1.37, tpole = 4800.0, inclination = 78.0,
                           rotation_period = 54.0)
tepochs = [d.mean_mjd for d in data]; tepochs .-= minimum(tepochs)
stars = create_star_multiepochs(tess, p, tepochs)
setup_oi!(data, stars)
x0 = Float64.(parametric_temperature_map(p, stars[1]))
chi2_breakdown(x0, stars, data)

# One reconstruction, short: `spheroid_crit_allepochs_fg` plus the regularizers is the whole
# imaging inner loop, and VMLMB's own specialisations come with it.
regs = Any[Any["sobel2", 10.0, sobel_gradient_healpix(3), 1:stars[1].npix]]
image_reconstruct_oi(x0, data, stars; regularizers = regs, maxiter = 3, verbose = false)

# The binary path, which the Model tab reaches whenever a companion is ticked, and which
# imaging now has a gradient for.
traced("binary") do
    p2 = default_star_params(3; rpole = 0.5, tpole = 4500.0)
    s2 = create_star_multiepochs(tess, p2, tepochs; secondary = true)
    setup_oi!(data, s2)
    x2 = Float64.(parametric_temperature_map(p2, s2[1]; secondary = true))
    ph = [binary_phase_shift(d.uv, 3.0, 1.5) for d in data]
    binary_chi2_f(x0, stars[1], x2, s2[1], data[1], ph[1])
    g = zeros(length(x0) + length(x2))
    binary_crit_allepochs_fg(vcat(x0, x2), g, stars, s2, data, ph)
end

# ── the plots ────────────────────────────────────────────────────────────────
#
# `Makie.save` is what forces the render, and rendering is where most of the compilation is.
# Building a Figure and stopping leaves the expensive half untraced.
traced("plots") do
    # Each returns `(fig, ax)`; `save` wants the figure.
    for (fig, _) in (plot2d_makie(x0, stars[1]),
                     plot2d_makie(x0, stars[1]; intensity = true, plotmesh = true),
                     plot3d_makie(x0, stars[1]),
                     plot_mollweide_makie(x0, stars[1]))
        Makie.save(OUT, fig)
    end
end

# ── everything below the window ───────────────────────────────────────────────
#
# The canvases and both tab refreshes: the work `gui()` does BEFORE `loadqml`, which is the
# part a user waits through on a cold start and the part that can be traced without handing a
# thread to Qt. The window itself is not opened here — see the header.
traced("the canvases") do
    session = GUI.Session()
    GUI.load_dataset!(session, FILES)
end

@info "trace complete"
