# The window: Qt callback registration and `gui()`.
#
# Split from the rest of src/gui/ because everything else describes what to draw, while this
# builds the thing that draws it and hands control to Qt's event loop.

function __init__()
    # Callbacks must be registered before any QML file that calls them is loaded.
    QML.@qmlfunction(shell_ready, shell_console, shell_status, shell_ui_scale,
                     shell_refresh,
                     shell_open, shell_open_many, shell_datasets, shell_epochs,
                     shell_remove_epoch, shell_close_dataset,
                     shell_select_dataset, shell_select_epoch, shell_current_dataset,
                     shell_current_epoch,
                     shell_obs_kinds, shell_set_obs_view,
                     shell_surface_field, shell_set_surface_field,
                     shell_decorations, shell_set_decoration, shell_clear_model,
                     shell_graticule, shell_set_graticule, shell_graticule_colors,
                     shell_tessellation, shell_set_tessellation,
                     shell_polyft_backends, shell_polyft_backend,
                     shell_set_polyft_backend,
                     shell_orbit_params, shell_set_orbit_param, shell_set_orbit_state,
                     shell_set_orbit_bound, shell_set_orbit_tie,
                     shell_orbit_component_kinds, shell_orbit_components,
                     shell_set_orbit_component, shell_orbit_options,
                     shell_orbit_star_models, shell_orbit_star_model,
                     shell_set_orbit_star_model,
                     shell_save_map, shell_load_map, shell_save_figure,
                     shell_fits, shell_current_fit, shell_select_fit, shell_fit_params,
                     shell_posterior_pair, shell_set_posterior_pair,
                     shell_set_orbit_option, shell_orbit_render_params,
                     shell_set_orbit_render_param, shell_fit_orbit,
                     shell_orbit_result, shell_save_orbit, shell_load_orbit,
                     shell_surface_types, shell_add_model, shell_models, shell_select_model,
                     shell_params, shell_set_param, shell_set_param_state,
                     shell_set_bound, shell_set_tie, shell_validate_model,
                     shell_fit_methods, shell_fit, shell_last_fit,
                     shell_colormaps, shell_set_colormap, shell_reset_zoom,
                     shell_regularizer_kinds, shell_reconstruct, shell_images,
                     shell_imaging_context,
                     shell_job_poll, shell_job_stop, shell_job_running,
                     shell_export, shell_script,
                     shell_save_settings, shell_load_settings,
                     shell_set_plot_scale, shell_plot_scale,
                     # The in-window file picker: QtQuick.Dialogs can leave its window mapped,
                     # and a native dialog gets its own GL surface, which is one of the two
                     # known ways to break the shared context permanently.
                     picker_list, picker_places, picker_parent, picker_join, picker_start)
    return nothing
end

"""
Folder the file dialog opens in.

Next to the last dataset loaded if there is one, then `\$ROTIRGUI_DATA_DIR`, then the package's
`demos/data`, which is the user-facing data directory — λ And, Spica and β Lyr all live there.
Opening in the filesystem root is a small tax paid on every single open, and it makes automated
clicking guesswork.
"""
function _initial_folder(session::Session)
    if !isempty(session.datasets)
        d = dirname(abspath(session.datasets[end].paths[end]))
        isdir(d) && return "file://" * d
    end
    forced = get(ENV, "ROTIRGUI_DATA_DIR", "")
    isempty(forced) || (isdir(forced) && return "file://" * abspath(forced))
    p = joinpath(pkgdir(ROTIR), "demos", "data")
    isdir(p) && return "file://" * p
    return "file://" * pwd()
end

# What the title bar says about a session that was populated before the window opened — by the
# launcher's command-line arguments, say. A fixed "no dataset loaded" was wrong exactly when a
# dataset HAD been loaded, which is the case a user notices.
function _session_status(s::Session)
    isempty(s.datasets) && return "no dataset loaded"
    d = s.datasets[end]
    return "$(d.name): $(length(d.data)) epoch(s)"
end

"Tab names, in the order Main.qml lists them."
const TAB_NAMES = ("data", "model", "imaging", "orbit")

"""
Which perspective the window opens on, from `\$ROTIRGUI_TAB` (default Data).

Its reason for existing is testing: a tab that is never current is constructed but never laid
out, so any layout warning it would emit never appears. Naming a tab is what makes an automated
run render it. An unknown name is reported rather than silently ignored, since a typo would
otherwise look like the setting having no effect.
"""
function _initial_tab()
    want = lowercase(strip(get(ENV, "ROTIRGUI_TAB", "")))
    isempty(want) && return 0
    i = findfirst(==(want), TAB_NAMES)
    if i === nothing
        @warn "ROTIRGUI_TAB is not a tab name; opening on Data" got = want options = TAB_NAMES
        return 0
    end
    return i - 1        # QML indexes from zero
end

"""
    gui(session = Session(); qmlfile, autoquit_ms = 0) -> Session

Open the ROTIR window.

    using ROTIR, GLMakie, QMLMakie, QML
    gui()

Everything is built BEFORE `loadqml`: four figures, every plot inside them, and the glyph
atlas. Once Qt owns the GL context, inserting a plot allocates buffers with none bound — see
the note at the top of src/gui/livecanvas.jl — so `loadqml` is the point after which only
Observable assignment is allowed.

`autoquit_ms` closes the window after that many milliseconds, which is how an automated run
gets to exercise the layout without a human closing it.
"""
function ROTIR.gui(session::Session = Session();
                   qmlfile::AbstractString = joinpath(pkgdir(ROTIR), "src", "gui", "qml",
                                                      "Main.qml"),
                   autoquit_ms::Integer = 0)
    check_qt_conflict()
    # A backstop for sessions started by hand rather than through bin/rotirgui.jl. GLMakie is
    # imported by now but has not built a context yet, and Mesa reads these at context
    # creation, so this is still early enough. A no-op when the launcher already ran it.
    gfx = OITOOLS.configure_graphics!()
    on_x11 = GLMakie.GLFW.GetPlatform() == GLMakie.GLFW.PLATFORM_X11
    qt = OITOOLS.configure_qt_platform!(; match_x11 = on_x11, verbose = false)

    uiscale = ui_scale_override()

    # Fill the glyph atlas before ANY figure exists. Text first rasterised after the Qt window
    # is up renders corrupted; the same text set beforehand is clean. A plain statement,
    # deliberately: inside a `@debug` string this would never run, because Julia's logging
    # macros skip their interpolations when the level is off.
    nglyphs = prewarm_glyphs!()

    # SIX figures, and every one of them is separate on purpose.
    #
    # Four kinds of axis — a 2-D sky projection, a 3-D scene, a bar chart, a whole-sphere
    # projection — and no axis can be more than one of them. On top of that, a Makie Figure can
    # be attached to ONE MakieArea at a time: pointing the Imaging tab's areas at the Data and
    # Model tabs' figures leaves one area blank and draws the other detached over the panel
    # beside it. So Imaging gets its own sky and Mollweide, which also means a reconstruction
    # cannot wipe the plot being read on another tab.
    skyfig  = Makie.Figure(size = (760, 700), backgroundcolor = :white)
    starfig = Makie.Figure(size = (760, 700), backgroundcolor = :white)
    mollfig = Makie.Figure(size = (900, 560), backgroundcolor = :white)
    chi2fig = Makie.Figure(size = (760, 360), backgroundcolor = :white)
    imskyfig  = Makie.Figure(size = (760, 620), backgroundcolor = :white)
    immollfig = Makie.Figure(size = (900, 560), backgroundcolor = :white)
    # The Model tab's own sky view, and the reason it exists rather than sharing the Data
    # tab's: a Figure belongs to one MakieArea at a time. This is the DEFAULT view there —
    # the sky projection is what the data actually constrains, and the 3-D scene is a
    # deliberate step away from it.
    mskyfig = Makie.Figure(size = (760, 700), backgroundcolor = :white)
    # The Data tab's observable plot — V² against baseline, closure phase, uv coverage. It uses
    # OITOOLS' font set, because it IS an OITOOLS canvas.
    obsfig  = Makie.Figure(size = (860, 640), backgroundcolor = :white)
    orbfig  = Makie.Figure(size = (820, 760), backgroundcolor = :white)

    sky    = build_sky_canvas(skyfig)
    star   = build_star_canvas(starfig)
    moll   = build_moll_canvas(mollfig)
    chi2   = build_chi2_canvas(chi2fig)
    imsky  = build_sky_canvas(imskyfig)
    immoll = build_moll_canvas(immollfig)
    msky   = build_sky_canvas(mskyfig)
    obs, obsmodel = build_obs_canvas(obsfig)
    orbitcanvas   = build_orbit_canvas(orbfig)
    postfig = Makie.Figure(size = (900, 520), backgroundcolor = :white)
    post    = build_post_canvas(postfig)
    # The Imaging tab's own 3-D scene: a Figure belongs to ONE MakieArea, so the Model tab's
    # cannot be shared with it.
    imstarfig = Makie.Figure(size = (760, 700), backgroundcolor = :white)
    imstar    = build_star_canvas(imstarfig)
    obs === nothing &&
        @warn "OITOOLSGUIExt is not loaded, so the Data tab's observable plot is unavailable"

    sh = ShellState(session, sky, star, moll, chi2, imsky, immoll, msky,
                    obs, obsmodel, Ref(:v2), Ref(:baseline), Ref(true),
                    # Intensity, not temperature, by DEFAULT: the map is a temperature but what
                    # the interferometer measures is the emergent intensity, so the intensity
                    # is the picture a χ² can be reasoned about from.
                    Ref(true), Ref(:linear), Ref(0.0),
                    Dict(:limb => true, :compass => true, :graticules => false,
                         :spin => false, :plotmesh => false),
                    Ref(30.0), Ref(30.0), Ref("black"),
                    Ref(:healpix), Ref(3), Ref{DataType}(Float32),
                    default_orbit(), orbitcanvas, "",
                    String[],
                    _session_status(session), nothing, :none, "", Dict{Symbol,Any}(),
                    Ref{Any}(nothing), Ref{Any}(nothing),
                    post, Ref(1), Ref(2), Ref(false), imstar,
                    Ref{Any}(nothing))
    SHELL[] = sh
    # Right-click to reset the zoom, and the zoom bound itself. Both read Makie's event
    # stream, so they attach to the canvases rather than to anything in QML.
    install_interactions!(sh)

    # Seed the console before the window exists, so it opens with a transcript rather than an
    # empty box — including the graphics decisions, which are exactly what you want to see when
    # the window looks wrong.
    # The orbit picture is drawn once here, before the window: it needs no data, so there is
    # no reason for the tab to open empty and fill in on the first click.
    refresh_orbit!(sh)
    console!(sh, "ROTIRGUI ready")
    console!(sh, "graphics: $(gfx.reason)")
    console!(sh, "qt platform: $(qt.reason)")
    console!(sh, "glyph atlas: $(nglyphs) glyphs pre-warmed")
    console!(sh, "ui scale: $(uiscale.reason)")
    isempty(session.datasets) && console!(sh, "no dataset loaded — use Open OIFITS…")

    # `skyfig`, not `skyfig.scene`: QMLMakie's MakieArea resolves the scene itself. Handing it
    # the root Scene directly skips whatever the figure-level path sets up.
    QML.loadqml(qmlfile;
                skyPlot         = skyfig,
                starPlot        = starfig,
                mollPlot        = mollfig,
                chi2Plot        = chi2fig,
                imSkyPlot       = imskyfig,
                imMollPlot      = immollfig,
                modelSkyPlot    = mskyfig,
                obsPlot         = obsfig,
                orbitPlot       = orbfig,
                postPlot        = postfig,
                imStarPlot      = imstarfig,
                autoQuitMs      = Int(autoquit_ms),
                initialTab      = _initial_tab(),
                uiScaleOverride = uiscale.scale,
                initialFolder   = _initial_folder(session),
                initialStatus   = sh.status)
    QML.exec()

    # An automated run has no way to read the console pane once the window is gone, so hand it
    # over on the way out. Only active when asked for; production runs never write a file.
    dump = get(ENV, "ROTIRGUI_CONSOLE_DUMP", "")
    if !isempty(dump)
        try
            open(dump, "w") do io
                for l in sh.console; println(io, l); end
            end
        catch err
            @warn "could not write console dump" path = dump exception = err
        end
    end
    return session
end
