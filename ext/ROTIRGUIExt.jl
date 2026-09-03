# The ROTIR GUI: three perspectives over one session, in a QML window drawing with Makie.
#
#     using ROTIR, GLMakie, QMLMakie, QML
#     gui()
#
# The whole GUI lives here rather than in the core package because Makie and Qt together cost
# several seconds of load time, and a reduction script, a CI run, or any session that never
# opens a window should not pay it.
#
# The trigger lists four packages but a caller only ever names three: GLMakie pulls in Makie,
# so `using GLMakie, QMLMakie, QML` loads all four and fires this extension. Makie is named
# because the sources below use it directly, and an extension may only `using` the parent's
# dependencies and its own triggers.
#
# `configure_graphics!` and `configure_qt_platform!` are NOT here, and not reimplemented
# either — they come from OITOOLS, which ROTIR already depends on. Both must run before the
# first GL context exists, i.e. before `using GLMakie`, i.e. before this extension can exist at
# all, which is why bin/rotirgui.jl calls them itself and `gui()` only repeats them as a
# backstop.
#
# Sources live in `src/gui/` and are included from here, so they stay editable in the normal
# place. `Session`, `ShellState` and the canvases are therefore extension types: functions can
# be forward-declared in the parent (`function gui end`), types cannot. `gui()` builds its own
# `Session`, so callers never need to name one; tests reach the rest through
# `Base.get_extension(ROTIR, :ROTIRGUIExt)`.
module ROTIRGUIExt

using ROTIR
import ROTIR: gui
using Printf, Statistics, TOML
using NLopt
using OptimPackNextGen
using LinearAlgebra
import OITOOLS

# The Makie drawing layer belongs to ROTIRMakieExt. `using GLMakie` activates it, so it is
# always present by the time this extension exists — but one extension cannot import from
# another, so these arrive through the functions the parent package declares.
using ROTIR: tessel_polygons, map_colors, sky_axis_max, star_mesh, relative_orbit_track,
             mollweide_xy, mollweide_grid
# Not exported by the parent, and needed by the canvases and the shell.
using ROTIR: _padded_cmap, _map_range, _body_radius, graticule_segments, convex_hull_2d,
             _spin_axis,
             # The nested-sampler dispatcher. Not imported until a test actually ran
             # `:nautilus` through the panel, at which point both it and `:ultranest` failed
             # with `UndefVarError: _fit_nested` — the fit methods were listed, gated on
             # availability, and unreachable.
             _fit_nested,
             # Orbit-tab internals: the element list and its default boxes are
             # not exported, and the tab's form is generated from them.
             ORBIT_ELEMENT_BOUNDS

"""
    prewarm_glyphs!() -> Int

Fill Makie's texture atlas before any Figure exists.

Text first rasterised after the Qt window is up renders corrupted; the same text set beforehand
is clean. OITOOLS measured that and carries the glyph list, so this delegates rather than
keeping a second copy that would drift. It is reached at RUNTIME through the OITOOLS module
rather than imported, because the implementation lives in OITOOLSMakieExt — a sibling of
OITOOLS' own extensions, which this one may not import from. By the time this is called
GLMakie is loaded, so that extension exists; the guard covers the case where it does not.
"""
prewarm_glyphs!() = isempty(methods(OITOOLS.prewarm_glyphs!)) ? 0 :
                    OITOOLS.prewarm_glyphs!()

using Makie
using QML, QMLMakie, GLMakie
using PrecompileTools

# Exported so tests and scripts can `using .ROTIRGUIExt` after fetching the module with
# `Base.get_extension`. `gui` itself is exported by ROTIR, since it is declared there.
export Session, DatasetEntry, ModelEntry, ImageEntry, LogEntry, OrbitEntry,
       load_dataset!, add_model!, add_image!, star_params, export_script

include(joinpath(pkgdir(ROTIR), "src", "gui", "commandlog.jl"))
include(joinpath(pkgdir(ROTIR), "src", "gui", "scaling.jl"))
include(joinpath(pkgdir(ROTIR), "src", "gui", "session.jl"))
include(joinpath(pkgdir(ROTIR), "src", "gui", "livecanvas.jl"))
include(joinpath(pkgdir(ROTIR), "src", "gui", "filepicker.jl"))
include(joinpath(pkgdir(ROTIR), "src", "gui", "shell.jl"))
include(joinpath(pkgdir(ROTIR), "src", "gui", "orbit.jl"))
include(joinpath(pkgdir(ROTIR), "src", "gui", "snapshot.jl"))
include(joinpath(pkgdir(ROTIR), "src", "gui", "window.jl"))

# ── precompilation ───────────────────────────────────────────────────────────────────────
#
# The core package's workload cannot reach any of this: these methods exist only once GLMakie,
# QML and QMLMakie are loaded. Everything below runs headless — building a Figure and drawing
# into it needs no backend, and no window is created — so it is safe in any environment.
#
# What is worth precompiling here is what happens BEFORE the first window appears, because that
# is time the user spends looking at nothing: building the four canvases (one of every plot the
# shell can ever show), reading a dataset, evaluating a model on it, and the first refresh of
# each tab. The GUI's own callbacks are cheap by comparison but they are also where the first
# CLICK goes, so the table-rendering ones are included too.
@setup_workload begin
    _pcfile = joinpath(pkgdir(ROTIR), "demos", "data", "2011Sep02.lam_And_prepped.oifits")
    @compile_workload begin
        # The pure-Julia half: no plotting, so this is compiled whatever the environment.
        _s = Session()
        _m = add_model!(_s, 2; rpole = 0.6, tpole = 4800.0, inclination = 78.0,
                        rotation_period = 54.0, frac_escapevel = 0.2)
        star_params(_m)
        apply_model_ties(_m)
        eval_tie("2*rpole", _m.params)
        export_script(_s)

        # The canvases. `build_*_canvas` creates one of every plot the shell can show, and that
        # construction is a large part of the delay before the first window appears — it is
        # also the part that cannot be deferred, because inserting a plot after the window
        # exists is the GL-context failure this whole file is arranged around.
        #
        # `build_obs_canvas` is NOT reachable from here and the call below always takes the
        # early return. MEASURED: `_oitools_gui()` is `nothing` during this workload, every
        # time — Julia does not activate a dependency's extensions while precompiling one of
        # this package's, and `Base.retry_load_extensions()` does not change that. So the
        # OITOOLS half of the observable canvas can only be precompiled by OITOOLS, which
        # does precompile it; ROTIR's job is to not RUIN that, which is why PythonCall is a
        # weak dependency here (see src/fit_ultranest.jl for the 338 ms → 2477 ms
        # measurement). The line stays because it costs nothing and becomes useful the day
        # Julia loads sibling extensions.
        _oitools_gui() === nothing || build_obs_canvas(Makie.Figure())
        _skyfig  = Makie.Figure(); _sky  = build_sky_canvas(_skyfig)
        _starfig = Makie.Figure(); _star = build_star_canvas(_starfig)
        _mollfig = Makie.Figure(); _moll = build_moll_canvas(_mollfig)
        _chifig  = Makie.Figure(); _chi  = build_chi2_canvas(_chifig)

        if isfile(_pcfile)
            load_dataset!(_s, _pcfile)
            _sh = ShellState(_s, _sky, _star, _moll, _chi, _sky, _moll, _sky,
                             nothing, Makie.Observable(Makie.Point2f[]),
                             Ref(:v2), Ref(:baseline), Ref(true),
                             Ref(true), Ref(:linear), Ref(0.0),
                             Dict(:limb => true, :compass => true, :graticules => false,
                                  :spin => false, :plotmesh => false),
                             Ref(30.0), Ref(30.0), Ref("black"),
                             Ref(:healpix), Ref(3), Ref{DataType}(Float32),
                             default_orbit(), nothing, "", String[], "",
                             nothing, :none, "", Dict{Symbol,Any}(),
                             Ref{Any}(nothing), Ref{Any}(nothing),
                             build_post_canvas(Makie.Figure()), Ref(1), Ref(2), Ref(false),
                             build_star_canvas(Makie.Figure()))
            SHELL[] = _sh
            # The two refreshes are the whole redraw path: geometry, temperature map, colours,
            # polygons, the 3-D mesh, the Mollweide resampling and the per-epoch χ². MEASURED
            # at 1.21 s on first call and 0.00 s afterwards — all of it compilation, and all of
            # it spent while the user is looking at a window that has not drawn yet.
            refresh_data_tab!(_sh)
            refresh_model_tab!(_sh)
            # Again with every DECORATION on: the graticule and spin-axis polylines go through
            # `graticule_segments` and `_spin_axis`, which the plain refresh never touches, so
            # the first tick would otherwise pay for them.
            for k in keys(_sh.decor); _sh.decor[k] = true; end
            refresh_data_tab!(_sh)
            refresh_model_tab!(_sh)
            # And with the colour coming from the intensity rather than the temperature, which
            # is a different code path through `surface_values` and the Planck conversion.
            shell_set_surface_field("1", "planck", "1.65")
            shell_set_surface_field("0", "linear", "0")
            shell_set_colormap("viridis")
            shell_set_colormap("gist_heat")
            _sh.chi2key[] = nothing
            model_state(_sh); epoch_chi2(_sh)

            # THROUGH THE SHELL CALLBACKS, not the functions under them. `shell_add_model`
            # calls `add_model!` with NO keywords while the line above calls it with five, and
            # those are different specialisations of the keyword body: measured at 43 ms for
            # `add_model!` and 23 ms for `default_star_params` still compiling on the first
            # "+ model" click with the low-level call alone in this workload.
            for _code in SURFACE_TYPE_ORDER
                shell_add_model(_code)
                refresh_data_tab!(_sh)
            end
            shell_add_model(2)
            shell_clear_model()
            shell_add_model(2)

            # The model prediction drawn over the data. `refresh_obs!` returns early here
            # because there is no observable canvas (see above), so the ROTIR half — geometry
            # at every epoch, `setup_oi!`, `observables` — is reached directly. Each view is
            # its own specialisation and each was ~100 ms on first switch.
            let _d1 = _s.datasets[1].data[1]
                for _k in (:v2, :t3amp, :t3phi, :uv)
                    _sh.obskind[] = _k
                    _model_points(_sh, _d1, 1)
                end
                _sh.obskind[] = :v2
            end
            # The callbacks a first click lands in.
            shell_epochs(); shell_datasets(); shell_params(); shell_surface_types()
            shell_regularizer_kinds(); shell_fit_methods(); shell_validate_model()
            shell_console(); shell_images(); shell_job_poll()
            shell_obs_kinds(); shell_current_epoch(); shell_select_epoch(1)
            shell_decorations(); shell_graticule(); shell_tessellation()
            picker_places()
            for _purpose in ("data", "orbit", "map")
                picker_list(dirname(_pcfile), "0", _purpose)
            end
            # The Orbit tab's own tables and its two renderings, which are a different path
            # from the model tabs entirely.
            shell_orbit_params(); shell_orbit_options(); shell_orbit_component_kinds()
            shell_orbit_star_models(); shell_orbit_render_params()
            # The posterior panel and the snapshot path, both of which build a Figure of
            # their own — the dearest thing either does and all of it compilation.
            shell_fits(); shell_fit_params(); shell_posterior_pair(); refresh_posterior!(_sh)
            _snapshot_figure(_sh, "sky", (400, 300))
            _snapshot_figure(_sh, "chi2", (400, 300))
            # The regulariser construction, which the Reconstruct button runs before the
            # engine sees anything — and which builds a sparse operator over the whole mesh.
            _specs = parse_regularizers("sobel2:10.0:0;radflat:100.0:6;mem:1.0:0")
            let _p = star_params(_m), _st = create_star(tessellation_healpix(3), _p, 0.0)
                # `setup_oi!` first: `radflat_bins` weights its annuli by `star.polyflux`,
                # which stays empty until the Fourier setup has run.
                setup_oi!([_s.datasets[1].data[1]], [_st])
                _x0 = Float64.(parametric_temperature_map(_p, _st))
                # The Reconstruct button runs this before the engine sees anything, and it
                # assembles a sparse operator over the whole mesh.
                _regs = build_regularizers(_specs, 3, _st, _x0, _p)
                _g = zeros(length(_x0))
                spheroid_regularization(_x0, _g; regularizers = _regs, verbose = false)
            end
            # The analytic shape fit, in the Float64 specialisation the GUI uses. Measured at
            # 9.6 s on first call against ~1 s of actual work, all of it compiling
            # `shape_chi2_fg!` for this element type — which happens while the user is
            # watching a Fit button that has not done anything yet.
            let _t64 = tessellation_healpix(2; T = Float64),
                _p64 = star_params(_m),
                _s64 = create_star(_t64, _p64, 0.0),
                _d1  = _s.datasets[1].data[1]
                _x64 = Float64.(parametric_temperature_map(_p64, _s64))
                _gθ  = zeros(Float64, 4); _gx = zeros(Float64, length(_x64))
                _θ   = Float64[_p64.rpole, _p64.frac_escapevel, _p64.inclination,
                               _p64.position_angle]
                shape_chi2_fg!(_gθ, _gx, _x64, _θ, [_d1], _t64, _p64, [0.0])
            end
            SHELL[] = nothing
        end
    end
end

end # module
