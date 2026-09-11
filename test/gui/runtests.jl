# Headless GUI tests: every shell callback, driven from Julia, with no window.
#
#     julia --project=bin test/gui/runtests.jl
#
# `--project=bin`, not `--project=.` — GLMakie, QML and QMLMakie are WEAK dependencies of
# ROTIR, so they are not loadable from the package environment at all, and neither is the
# extension that defines everything below.
#
# WHAT THIS COVERS, AND WHAT IT CANNOT.
#
# Every callback QML can reach, the session model underneath them, and the canvases — which are
# built for real, because building them is where the plot-once constraint lives. What it does
# NOT cover is anything that needs a WINDOW: opening a popup, the GL context that a popup
# takes over, and whether a click lands on the control it appears to. `gui_click.sh` does that,
# and this file passes on a machine where the GUI is visibly broken. Both are needed; neither
# replaces the other.

using Test
using ROTIR
using GLMakie, QMLMakie, QML          # activates ROTIRGUIExt
using Makie

const G = Base.get_extension(ROTIR, :ROTIRGUIExt)
G === nothing && error("ROTIRGUIExt did not load; GLMakie, QMLMakie and QML are all needed")

const DATA = joinpath(pkgdir(ROTIR), "demos", "data")
const LAM  = sort([joinpath(DATA, f) for f in readdir(DATA) if occursin("lam_And", f)])

"A ShellState with every canvas built, as `gui()` builds them."
function fresh_shell()
    s = G.Session()
    figs = [Makie.Figure() for _ in 1:10]
    # BY KEYWORD: only the canvases this harness builds. Everything else — the view state,
    # the decorations, the caches — takes the default that lives with the field in
    # `ShellState`, so a field added there does not have to be added here too. It used to be
    # positional, and adding one meant editing this, `gui()`, and the `@compile_workload` in
    # ext/ROTIRGUIExt.jl, with a 39-argument `MethodError` as the only warning.
    #
    # No observable canvas: it needs OITOOLSGUIExt's `build_canvas`, and `refresh_obs!` is a
    # no-op without it. The observable plot is covered by gui_click.sh, which has a window.
    sh = G.ShellState(; session = s,
        sky    = G.build_sky_canvas(figs[1]),  star  = G.build_star_canvas(figs[2]),
        moll   = G.build_moll_canvas(figs[3]), chi2  = G.build_chi2_canvas(figs[4]),
        imsky  = G.build_sky_canvas(figs[5]),  immoll = G.build_moll_canvas(figs[6]),
        msky   = G.build_sky_canvas(figs[7]),
        obsmodel = Makie.Observable(Makie.Point2f[]),
        orbitcanvas = G.build_orbit_canvas(figs[8]),
        post   = G.build_post_canvas(figs[9]),
        imstar = G.build_star_canvas(figs[10]))
    G.SHELL[] = sh
    return sh
end

"Run a job callback to completion, the way the QML poll timer does."
function drain!(; limit = 600)
    for _ in 1:limit
        G.shell_job_running() == "1" || break
        sleep(0.2); G.shell_job_poll()
    end
    G.shell_job_poll()
    return nothing
end

# The nine orbital elements, by the name the Orbit tab spells them with. They are always
# listed, whatever the star model is.
const ORBIT_ELEMENTS_NAMES = ["a", "i", "Omega", "omega", "e", "P", "T0", "dP", "domega"]

rows(s) = isempty(s) ? String[] : split(s, '\n')
cols(r) = split(r, '\t')

@testset "ROTIR GUI (headless)" begin

@testset "the GUI session loads no Python" begin
    # MEASURED, and the reason PythonCall is a weak dependency: with it loaded, OITOOLS'
    # precompiled `build_canvas` is invalidated and one canvas build goes from 338 ms to
    # 2477 ms — 1.2 s on every GUI start, sampler or no sampler. Nothing in the GUI path may
    # pull it back in, and this is what says so before a user notices the second.
    @test !any(k -> k.name == "PythonCall", keys(Base.loaded_modules))
    @test !ultranest_available()
    # Pigeons is the SAME rule and a bigger number: 370 ms → 3166 ms for one canvas build,
    # so the launcher does not load it either and `:pigeons` is offered only to a session
    # that asked for it. Nothing the launcher does load costs anything — Zygote, AdvancedHMC,
    # LogDensityProblems and Nautilus were each measured at 370-380 ms.
    @test !any(k -> k.name == "Pigeons", keys(Base.loaded_modules))
    @test !pigeons_available()
    # ...and the sampler that needs it is not offered, rather than offered and failing.
    sh = fresh_shell()
    G.shell_open(LAM[1], "0")
    G.shell_add_model(2)
    G.shell_set_param_state("rpole", "free")
    # `:ultranest` is ABSENT from the list, not merely gated — the GUI must have no method
    # whose availability depends on loading Python — and asking for it by name explains that
    # rather than reporting an unknown method.
    @test !occursin("ultranest", G.shell_fit_methods())
    @test occursin("does not offer :ultranest", G.shell_fit("ultranest", 100))
    # Pigeons IS gated rather than absent: it costs nothing to have present, and the launcher
    # loads it. This test session does not, so it must not be offered here.
    @test !occursin("pigeons", G.shell_fit_methods())
    @test occursin("Pigeons", G.shell_fit("pigeons", 100))
end

@testset "canvases build" begin
    # Every plot must exist before a window does — that is the whole arrangement of
    # src/gui/livecanvas.jl — so a build failure here is a build failure in `gui()`.
    sh = fresh_shell()
    @test sh.sky !== nothing && sh.star !== nothing
    @test sh.moll !== nothing && sh.chi2 !== nothing
    # Idle from the start, rather than an empty axis beside a 0-1 colorbar.
    @test !isempty(sh.sky.message[])
    @test !isempty(sh.star.message[])
end

@testset "loading" begin
    sh = fresh_shell()
    @test G.shell_datasets() == ""
    @test G.shell_epochs() == ""

    G.shell_open(LAM[1], "0")
    d = G.current_dataset(sh.session)
    @test d !== nothing
    @test length(d.data) == 1
    @test d.tepochs == [0.0]                       # relative to the first epoch, always

    # Several files at once become several EPOCHS of one dataset, not several datasets.
    sh2 = fresh_shell()
    G.shell_open_many(join(LAM, "\n"), "0")
    d2 = G.current_dataset(sh2.session)
    @test length(sh2.session.datasets) == 1
    @test length(d2.data) == length(LAM)
    @test d2.tepochs[1] == 0.0
    @test issorted(d2.mjd)                          # sorted by time, not by click order
    @test all(d2.tepochs .≈ d2.mjd .- minimum(d2.mjd))

    # One row per epoch, seven fields, and the observables as separate columns.
    r = rows(G.shell_epochs())
    @test length(r) == length(LAM)
    @test all(length(cols(x)) == 7 for x in r)
    @test cols(r[1])[7] == "0"                      # no model yet -> counts, not χ²
end

@testset "adding epochs and removing them" begin
    sh = fresh_shell()
    G.shell_open(LAM[1], "0")
    G.shell_open(LAM[2], "1")                       # add mode
    d = G.current_dataset(sh.session)
    @test length(sh.session.datasets) == 1
    @test length(d.data) == 2

    # Removing the FIRST epoch re-bases the origin: the phase is measured from the earliest
    # epoch that is still loaded, and a stale origin would silently rotate every map.
    G.shell_open(LAM[3], "1")
    before = copy(G.current_dataset(sh.session).mjd)
    G.shell_remove_epoch(1)
    d = G.current_dataset(sh.session)
    @test length(d.data) == 2
    @test d.mjd == before[2:end]
    @test d.tepochs[1] == 0.0

    # Removing the last epoch closes the dataset rather than leaving an empty one.
    G.shell_remove_epoch(1); G.shell_remove_epoch(1)
    @test isempty(sh.session.datasets)
    @test G.shell_epochs() == ""

    # Closing selects the previous dataset.
    G.shell_open(LAM[1], "0"); G.shell_open(LAM[2], "0")
    @test length(sh.session.datasets) == 2
    G.shell_close_dataset()
    @test length(sh.session.datasets) == 1
    @test G.shell_current_dataset() == "1"
end

@testset "the surface schema drives the form" begin
    sh = fresh_shell()
    # Every implemented surface type is offered, plus BINARY — which is not a surface type
    # and has no schema entry, because it is two components each with a type of its own. Its
    # code is negative so it cannot collide with a real `surface_type`.
    st = rows(G.shell_surface_types())
    @test length(st) == length(SURFACE_TYPE_ORDER) + 1
    @test Set(parse(Int, cols(r)[1]) for r in st) ==
          union(Set(SURFACE_TYPE_ORDER), Set([G.BINARY_CODE]))
    @test G.BINARY_CODE < 0
    @test !(G.BINARY_CODE in SURFACE_TYPE_ORDER)
    bin_row = only(filter(r -> parse(Int, cols(r)[1]) == G.BINARY_CODE, st))
    @test cols(bin_row)[2] == "binary"
    @test cols(bin_row)[3] == "Roche Binary"           # what it builds, said in the name
    # Choosing it builds BOTH components in one step, so a binary is one action rather than
    # "add a model, then remember to tick something".
    G.shell_add_model(G.BINARY_CODE)
    let m = G.current_model(sh.session)
        @test m.companion !== nothing
        @test m.surface_type == 3 && m.companion.surface_type == 3   # Roche both
        @test G.shell_binary() == "1"
        # PLACED BY THE ORBIT, and that is the difference from ticking `secondary` on some
        # other surface type. A Roche shape is computed FROM the instantaneous separation, so
        # this entry already has an orbit in it; placing the pair by a hand-typed offset while
        # shaping it from an orbit would say two things at once.
        @test m.companion.place === :orbit
        @test m.companion.offset == G.DEFAULT_COMPANION_OFFSET   # kept, not destroyed
        # The positions are a set of PARAMETERS, in the same 13 columns as every other form —
        # inert here, because the orbit is what decides where the secondary is.
        pr = rows(G.shell_position_params())
        @test length(pr) == 3
        @test all(length(cols(x)) == 13 for x in pr)
        @test [cols(x)[1] for x in pr] == ["pos_x", "pos_y", "pos_z"]
        @test all(cols(x)[13] == "1" for x in pr)
        @test occursin("orbit places", G.shell_set_position_state("pos_y", "free"))
        # Switch to a fixed offset and they come alive, with the separation intact.
        G.shell_set_binary_placement("offset")
        @test all(cols(x)[13] == "0" for x in rows(G.shell_position_params()))
        @test G.shell_set_position_state("pos_x", "free") == ""
        @test :pos_x in m.companion.free
    end
    G.shell_clear_model()

    # And the OTHER path — a model of some surface type, then the `secondary` tick — still
    # defaults to a fixed offset: one displacement is what a snapshot constrains, and an orbit
    # needs elements set up on another tab before it places anything at all.
    G.shell_add_model(0)
    G.shell_set_binary("1", 0)
    @test G.current_model(sh.session).companion.place === :offset
    G.shell_clear_model()

    for code in SURFACE_TYPE_ORDER
        G.shell_add_model(code)
        m = G.current_model(sh.session)
        @test m.surface_type == code
        @test isempty(G.shell_validate_model())      # defaults must build
        pr = rows(G.shell_params())
        @test length(pr) == length(surface_params(code))
        # THIRTEEN columns: the twelve the form draws, plus `inert` — whether the row is a
        # limb-darkening coefficient the current law does not read. The panel greys those and
        # `shell_set_param_state` refuses to free them.
        @test all(length(cols(x)) == 13 for x in pr)
        @test all(cols(x)[13] in ("0", "1") for x in pr)
        # Only limb-darkening coefficients are ever inert.
        @test all(cols(x)[1] in ("ld1", "ld2", "ld3", "ld4")
                  for x in pr if cols(x)[13] == "1")
        # The `ldtype` row is a choice field carrying its options.
        ld = findfirst(x -> cols(x)[1] == "ldtype", pr)
        @test ld !== nothing
        @test cols(pr[ld])[10] == "choice"
        @test occursin("=", cols(pr[ld])[11])
        # And it has no free/fixed/tied to make: a law index is not a coordinate an optimiser
        # can walk, so both the panel and the shell refuse it. Refused, not silently ignored —
        # a state that was accepted and then dropped would show the wrong thing in the form.
        @test occursin("discrete choice", G.shell_set_param_state("ldtype", "free"))
        @test occursin("discrete choice", G.shell_set_param_state("ldtype", "tied"))
        @test cols(rows(G.shell_params())[ld])[5] == "fixed"
        @test G.shell_set_param_state("ldtype", "fixed") == ""
        # Switching the LD law zeroes the coefficients it does not read, releases them, and
        # marks them inert. A value left behind is invisible in the fit, reappears when the
        # law is switched back, and if it was FREE it stays in the parameter vector as a
        # direction the χ² is exactly flat along.
        G.shell_set_param("ldtype", "4")
        for (n, v) in (("ld1", "0.3"), ("ld2", "0.4"), ("ld3", "0.5"), ("ld4", "0.6"))
            G.shell_set_param(n, v)
        end
        G.shell_set_param_state("ld2", "free")
        ldrow(nm) = only(filter(x -> cols(x)[1] == nm, rows(G.shell_params())))
        @test all(cols(ldrow(n))[13] == "0" for n in ("ld1", "ld2", "ld3", "ld4"))
        G.shell_set_param("ldtype", "1")                 # linear reads ld1 only
        @test cols(ldrow("ld1"))[13] == "0"
        @test all(cols(ldrow(n))[13] == "1" for n in ("ld2", "ld3", "ld4"))
        @test all(parse(Float64, cols(ldrow(n))[4]) == 0.0 for n in ("ld2", "ld3", "ld4"))
        @test cols(ldrow("ld2"))[5] == "fixed"           # released, not left free
        @test occursin("not used", G.shell_set_param_state("ld2", "free"))
        G.shell_set_param("ldtype", "2")                 # quadratic brings ld2 back
        @test cols(ldrow("ld2"))[13] == "0"
        @test G.shell_set_param_state("ld2", "free") == ""
        # Every label has to FIT the form's column, which elides silently rather than
        # wrapping: a truncated "Limb-darkening l…" is how this was noticed.
        for x in pr
            c = cols(x)
            @test length(c[2]) + (isempty(c[3]) ? 0 : length(c[3]) + 3) <= 20
        end
    end
    # ONE model, however many were added: "+ model" REPLACES rather than appends, so that
    # nothing but the visible model can decide what the χ² column is about.
    @test length(rows(G.shell_models())) == 1
    @test parse(Int, cols(rows(G.shell_models())[1])[2]) == last(SURFACE_TYPE_ORDER)
end

@testset "editing parameters, states, bounds and ties" begin
    sh = fresh_shell()
    G.shell_add_model(2)
    m = G.current_model(sh.session)

    @test G.shell_set_param("rpole", "1.5") == ""
    @test m.params[:rpole] == 1.5
    # A half-typed number is refused rather than silently zeroing the field.
    @test occursin("not a number", G.shell_set_param("rpole", "1.5e"))
    @test m.params[:rpole] == 1.5

    G.shell_set_param_state("rpole", "free")
    @test :rpole in m.free
    @test G.shell_set_bound("rpole", "0.5", "3.0") == ""
    @test m.bounds[:rpole] == (0.5, 3.0)
    @test occursin("below", G.shell_set_bound("rpole", "3.0", "0.5"))

    # A tie is an EXPRESSION, which is why the state is three-way and not a tick box.
    G.shell_set_param("inclination", "78.0")
    @test occursin("=", G.shell_set_tie("position_angle", "inclination - 60"))
    @test G.current_model(sh.session).params[:position_angle] ≈ 18.0
    @test !(:position_angle in m.free)
    # Operators and Base functions resolve; only the parameter names are substituted.
    G.shell_set_tie("beta", "0.25*sqrt(4)")
    @test G.apply_model_ties(m)[:beta] ≈ 0.5
    # A half-typed expression reports rather than raising.
    @test occursin("does not evaluate", G.shell_set_tie("beta", "0.25*sqrt("))

    G.shell_set_param_state("beta", "fixed")
    @test !haskey(m.ties, :beta)
end

@testset "per-epoch χ² and its cache" begin
    sh = fresh_shell()
    G.shell_open_many(join(LAM[1:3], "\n"), "0")
    G.shell_add_model(2)
    G.shell_set_param("rpole", "1.37131")
    G.shell_set_param("tpole", "4800.0")
    G.shell_set_param("rotation_period", "54.8")

    b = G.epoch_chi2(sh)
    @test b !== nothing && length(b) == 3
    @test all(isfinite(e.v2r) for e in b)
    # The columns now carry χ², and the flag says so.
    r = rows(G.shell_epochs())
    @test cols(r[1])[7] == "1"

    # Cached: the same key must not recompute. This is what made a tab switch cost 184 ms.
    t1 = @elapsed G.epoch_chi2(sh)
    sh.chi2key[] = nothing
    t2 = @elapsed G.epoch_chi2(sh)
    @test t1 < t2 / 10

    # Editing a parameter must INVALIDATE it, or the panel shows a fit that is no longer the
    # model on screen.
    v2_before = G.epoch_chi2(sh)[1].v2r
    G.shell_set_param("rpole", "2.0")
    @test G.epoch_chi2(sh)[1].v2r != v2_before
end

@testset "regularisers" begin
    sh = fresh_shell()
    # Every regulariser `spheroid_regularization` dispatches on is offered, and nothing else.
    kinds = rows(G.shell_regularizer_kinds())
    @test length(kinds) == 10
    names = Set(cols(k)[1] for k in kinds)
    @test names == Set(["sobel","sobel2","tv","tv2","mem","mean","bias",
                        "radflat","radialvar","orthold"])
    # SIX columns: `name weight extra_label extra_default short doc`. The row draws `short`
    # (the formula) and the tooltip carries `doc` (the advice) — one string could not do both
    # without eliding mid-sentence in the narrow rows, which are the three that also carry an
    # extra knob.
    @test all(length(cols(k)) == 6 for k in kinds)
    # The short form has to FIT, which is the whole point of splitting them.
    @test all(length(cols(k)[5]) <= 18 for k in kinds)
    @test all(!isempty(cols(k)[6]) for k in kinds)

    spec = join(("$(cols(k)[1]):$(cols(k)[2]):$(cols(k)[4])" for k in kinds), ";")
    specs = G.parse_regularizers(spec)
    @test length(specs) == 10
    @test G.parse_regularizers("bogus:1:0") == []          # unknown names dropped
    @test G.parse_regularizers("sobel:notanumber:0") == []

    # Building them needs a star with its Fourier setup — every entry is
    # [name, weight, aux, subset], and aux is a STRUCTURE for most of them.
    G.shell_open(LAM[1], "0")
    p = default_star_params(:rapid_rotator; rpole = 0.6, tpole = 4800.0)
    star = create_star(tessellation_healpix(3), p, 0.0)
    setup_oi!([sh.session.datasets[1].data[1]], [star])
    x0 = Float64.(parametric_temperature_map(p, star))
    regs = G.build_regularizers(specs, 3, star, x0, p)
    @test length(regs) == 10
    @test all(length(r) == 4 for r in regs)
    for r in regs
        g = zeros(length(x0))
        f = spheroid_regularization(x0, g; regularizers = Any[r], verbose = false)
        @test isfinite(f)
    end

    # radflat/radialvar weight their annuli by polyflux, so they need setup_oi! and say so.
    bare = create_star(tessellation_healpix(3), p, 0.0)
    @test_throws ErrorException G.build_regularizers(
        G.parse_regularizers("radflat:100:6"), 3, bare, x0, p)
end

# Shipped resources through `ROTIR.resource`, not `pkgdir`.
#
# This is what makes an application bundle possible: `create_app` copies no package source, so
# `pkgdir(ROTIR)` inside a bundle names a directory on the machine that BUILT it. Every lookup
# that reaches a shipped file has to go through `resource`, and the failure mode of a
# reintroduced `pkgdir` is silent on this machine — the checkout is still there — and only
# shows up in a bundle. `$ROTIR_RESOURCE_DIR` is what lets the whole mechanism be driven
# without building one.
@testset "shipped resources relocate" begin
    @test ROTIR.resource("no", "such", "thing") === nothing
    @test ROTIR.resource_dir() !== nothing              # a checkout is a valid root

    stage = mktempdir()
    mkpath(joinpath(stage, "src", "gui"))
    mkpath(joinpath(stage, "demos"))
    cp(joinpath(pkgdir(ROTIR), "src", "gui", "qml"), joinpath(stage, "src", "gui", "qml"))
    mkpath(joinpath(stage, "demos", "data"))
    mkpath(joinpath(stage, "demos", "orbits"))

    withenv("ROTIR_RESOURCE_DIR" => stage) do
        @test ROTIR.resource_dir() == stage
        # The forced root WINS over the checkout, or a bundle could be shadowed by whatever
        # happened to be sitting at pkgdir.
        for parts in (("src", "gui", "qml", "Main.qml"), ("demos", "data"), ("demos", "orbits"))
            p = ROTIR.resource(parts...)
            @test p !== nothing
            @test startswith(something(p, ""), stage)
        end
        # ...and the call sites follow it. These are the four that a bundle breaks.
        @test occursin(stage, G._initial_folder(G.Session()))
        @test any(r -> occursin(stage, r), rows(G.picker_places()))
    end

    # Outside the block the checkout answers again, so nothing leaked into global state.
    @test !startswith(something(ROTIR.resource("demos", "data"), ""), stage)
end

@testset "reconstruction" begin
    sh = fresh_shell()
    G.shell_open(LAM[1], "0")
    @test occursin("no model", G.shell_reconstruct(3, "sobel2:10:0", 5))
    # A BINARY is refused, not reconstructed as its primary. This panel builds one `stars`
    # vector from the primary's parameters, so with a companion in the model it fitted one
    # surface against data containing two and reported a χ² for it — a meaningless map that
    # nothing on screen distinguished from a good one.
    G.shell_add_model(0)
    G.shell_set_binary("1", 0)
    let msg = G.shell_reconstruct(3, "sobel2:10:0", 5)
        @test occursin("one component", msg)
        @test occursin("binary_reconstruct_oi", msg)
    end
    @test sh.job === nothing          # refused BEFORE a worker was started
    G.shell_clear_model()
    G.shell_add_model(2)
    G.shell_set_param("rpole", "1.37131")
    G.shell_set_param("tpole", "4800.0")

    @test G.shell_reconstruct(3, "radflat:100.0:6;sobel2:10.0:0", 10) == ""
    drain!()
    @test sh.job === nothing
    im = rows(G.shell_images())
    @test length(im) == 1
    @test length(sh.session.images) == 1
    e = sh.session.images[1]
    @test length(e.x) == 12 * (2^3)^2
    @test isfinite(e.chi2) && e.ndata > 0
    @test occursin("χ²", sh.status)
end

@testset "fitting" begin
    sh = fresh_shell()
    G.shell_open(LAM[1], "0")
    G.shell_add_model(0)
    @test occursin("nothing is free", G.shell_fit("neldermead", 100))

    G.shell_set_param("radius", "1.2")
    G.shell_set_param_state("radius", "free")
    G.shell_set_bound("radius", "0.5", "3.0")
    @test G.shell_fit("neldermead", 200) == ""
    drain!()
    @test occursin("χ²", sh.status)
    @test occursin("radius", sh.lastfit)
    @test G.current_model(sh.session).params[:radius] != 1.2

    @test occursin("unknown method", G.shell_fit("nosuchmethod", 10))

    # The gradient path is offered only where the gradient is CONSISTENT with the objective.
    ms = Set(cols(r)[1] for r in rows(G.shell_fit_methods()))
    @test "gradient" in ms                       # surface_type 0: uniform map, exact
    G.shell_add_model(3)                         # Roche: no analytic gradient at all
    @test !("gradient" in Set(cols(r)[1] for r in rows(G.shell_fit_methods())))
    @test G.gradient_fit_kind(G.current_model(sh.session)) === :none
end

@testset "temperature or intensity" begin
    sh = fresh_shell()
    G.shell_open(LAM[1], "0")
    G.shell_add_model(2)
    G.shell_set_param("rpole", "1.2")
    G.shell_set_param("frac_escapevel", "0.7")     # a real gravity-darkening gradient
    got = G.build_epoch_star(sh)
    @test got !== nothing
    star, tmap = got

    # INTENSITY is the default: the map is a temperature, but what an interferometer measures
    # is the emergent intensity, so that is the picture a χ² can be reasoned about from.
    @test sh.intensity[]

    # Temperature: the map itself, one value per visible tessel.
    G.shell_set_surface_field("0", "linear", "0")
    @test !sh.intensity[]
    tv = G.surface_values(sh, tmap, star; visible_only = true)
    @test length(tv) == length(star.index_quads_visible)
    @test tv ≈ Float64.(tmap[star.index_quads_visible])

    # Intensity: limb darkening multiplied in, so the limb is DARKER than the temperature says
    # and the values are no longer the map.
    G.shell_set_surface_field("1", "linear", "0")
    @test sh.intensity[]
    iv = G.surface_values(sh, tmap, star; visible_only = true)
    @test length(iv) == length(tv)
    @test !(iv ≈ tv)
    @test all(iv .<= tv .* (1 + 1e-9))             # ld <= 1 everywhere
    @test sh.sky.cbarlabel[] == "I (arb.)"

    # Planck: a real surface brightness, strongly non-linear in T, so the CONTRAST across the
    # star is larger than the linear proxy gives.
    @test occursin("planck", G.shell_set_surface_field("1", "planck", "1.65"))
    @test sh.band[] ≈ 1.65e-6
    pv = G.surface_values(sh, tmap, star; visible_only = true)
    @test all(isfinite, pv)
    @test !(pv ≈ iv)
    # The physical claim, stated without limb darkening in the way: at 1.65 µm the Planck
    # function compresses the COOL end far harder than a linear proxy, so the pole-to-equator
    # brightness ratio is larger than the temperature ratio. Asserting it on `pv` directly
    # cannot work — `ldmap` is exactly 0 at the limb, so both contrasts are Inf.
    Tlo, Thi = extrema(Float64.(tmap))
    Ilo, Ihi = ROTIR.intensity([Tlo, Thi], :planck, 1.65e-6)
    @test Ihi / Ilo > Thi / Tlo

    # Band 0 takes the wavelength from the data rather than needing one typed in.
    G.shell_set_surface_field("1", "planck", "0")
    @test sh.band[] > 0

    @test occursin("unknown", G.shell_set_surface_field("1", "nosuchmodel", "0"))
    G.shell_set_surface_field("0", "linear", "0")
    @test !sh.intensity[]
    @test sh.sky.cbarlabel[] == "T (K)"
    # The whole-surface label rule, stated once where the field is set. `linear` intensity of
    # a temperature map is that map, so on a view with no limb it is still "T (K)"; `planck`
    # is a different physical quantity and is labelled as one.
    G.shell_set_surface_field("0", "linear", "0")
    @test G.whole_surface_label(sh) == "T (K)"
    G.shell_set_surface_field("1", "linear", "0")
    @test G.whole_surface_label(sh) == "T (K)"
    @test sh.moll.cbarlabel[] == "T (K)" && sh.star.cbarlabel[] == "T (K)"
    @test sh.msky.cbarlabel[] == "I (arb.)"
    G.shell_set_surface_field("1", "planck", "1.65")
    @test G.whole_surface_label(sh) == "I (arb.)"
    @test sh.moll.cbarlabel[] == "I (arb.)" && sh.star.cbarlabel[] == "I (arb.)"
    G.shell_set_surface_field("0", "linear", "0")

    # Per-pixel indexing for the 3-D and Mollweide views, not per visible tessel.
    @test length(G.surface_values(sh, tmap, star; visible_only = false)) == star.npix

    # A WHOLE-SURFACE view has no limb, so limb darkening must not be applied to it. `ldmap`
    # is a viewing quantity and is exactly 0 on the half of the surface facing away, so
    # multiplying it in blanked half of every 3-D and Mollweide view.
    @test count(iszero, star.ldmap) >= div(star.npix, 2)     # the half that would be lost
    @test !any(iszero, G.surface_values(sh, tmap, star; visible_only = false))
    G.shell_set_surface_field("1", "linear", "0")
    @test !any(iszero, G.surface_values(sh, tmap, star; visible_only = false))
    G.shell_set_surface_field("0", "linear", "0")
end

@testset "views, colormaps and zoom" begin
    sh = fresh_shell()
    G.shell_open(LAM[1], "0")
    G.shell_add_model(2)
    G.refresh_data_tab!(sh)

    @test isempty(sh.sky.message[])                      # busy once something is drawn
    @test !isempty(sh.sky.lastvalues[])
    @test sh.sky.homespan[] > 0

    @test occursin("colormap", G.shell_set_colormap("viridis"))
    @test occursin("unknown", G.shell_set_colormap("not-a-colormap"))
    @test length(rows(G.shell_colormaps())) == length(G.SURFACE_COLORMAPS)

    # Zoom is bounded: overzooming is refused, not clamped into a loop.
    a = sh.sky.homespan[] / 2
    Makie.xlims!(sh.sky.axis, -a/100, a/100); Makie.ylims!(sh.sky.axis, -a/100, a/100)
    @test G.clamp_zoom!(sh.sky)
    G.reset_zoom!(sh.sky)
    @test !G.clamp_zoom!(sh.sky)

    # THE ORBIT VIEW IS BOUNDED TOO, and it is the one that was not.
    #
    # Without a limiter it kept Makie's own `ScrollZoom`, which receives QMLMakie's 120-unit
    # wheel event undivided — one notch drove the limits somewhere the renderer could not
    # survive, and it took the machine with it. The bound is what makes a wheel gesture here
    # safe, so it is pinned rather than left to the next person to rediscover.
    let c = sh.orbitcanvas
        @test c !== nothing
        G.show_orbit!(c, G.default_orbit(), [0.0, 1.0])
        @test c.homespan[] > 0                       # the frame was recorded
        home = c.homespan[]
        G.zoom_step!(c, 400)                         # far more than any real gesture
        span = abs(c.axis.finallimits[].widths[1])
        @test span >= G.ZOOM_MIN_SPAN * home * 0.99  # refused, not driven to zero
        G.zoom_step!(c, -400)
        span = abs(c.axis.finallimits[].widths[1])
        @test span <= G.ZOOM_MAX_SPAN * home * 1.01  # and bounded on the way out
        G.reset_zoom!(c)
        @test isapprox(abs(c.axis.finallimits[].widths[1]), home; rtol = 1e-3)
    end

    # THE REQUESTED GRATICULE SPACING IS WHAT GETS DRAWN.
    #
    # The panel used to convert degrees to a line COUNT, and the conversion could not express
    # the request: `nlat = round(180/75) = 2` parallels were then placed `180/(nlat+1) = 60`
    # apart, so asking for 75 drew 60. Degrees now go through untouched and parallels sit at
    # multiples of the spacing from the equator — 75 gives -75, 0, +75.
    #
    # Asserted through the run COUNT, which is what the spacing controls: a coarser spacing
    # must draw strictly fewer curves than a finer one, and both must draw something.
    G.shell_add_model(0)
    sh.decor[:graticules] = true
    G.shell_set_graticule(75.0, 45.0, "black")
    coarse = count(p -> isnan(p[1]), sh.msky.grat[])
    G.shell_set_graticule(30.0, 30.0, "black")
    fine = count(p -> isnan(p[1]), sh.msky.grat[])
    @test coarse > 0
    @test fine > coarse
    # Int arguments, as a QML SpinBox actually sends them — `String(Int32)` has no method and
    # the exception used to escape through `julia_call` and freeze the window.
    @test G.shell_set_graticule(Int32(45), Int32(60), "black") == ""
    @test sh.gratlat[] == 45.0 && sh.gratlon[] == 60.0
    G.shell_clear_model()

    # CLEARING A MODEL CLEARS THE PICTURE, not just the caption. `idle!` used to set the
    # message alone, so "− model" left the old star on screen with "no model" written over it
    # — measured at 422 polygons and 3072 mesh vertices still drawn. The caption was already
    # right, which is why the click test's own check did not catch this.
    G.shell_add_model(0)
    G.refresh_both!(sh)
    @test !isempty(sh.msky.polys[])
    G.shell_clear_model()
    @test isempty(sh.msky.polys[])
    @test isempty(sh.msky.limb[]) && isempty(sh.msky.grat[])
    @test length(Makie.GeometryBasics.coordinates(sh.star.mesh[])) == 3   # the placeholder
    @test !isempty(sh.msky.message[])

    # Back to idle when there is nothing to draw.
    G.shell_close_dataset()
    sh.session.current_model = 0
    G.refresh_data_tab!(sh)
    @test !isempty(sh.sky.message[])
end

@testset "the observable plot's log y axis" begin
    sh = fresh_shell()
    G.shell_open(LAM[1], "0")
    G.shell_add_model(2)

    # "1"/"0", like every other boolean crossing this boundary.
    @test G.shell_set_obs_view("v2", "baseline", "1", "1") isa String
    @test sh.obslog[] && G._obs_logscale(sh)
    G.shell_set_obs_view("v2", "baseline", "1", "0")
    @test !sh.obslog[] && !G._obs_logscale(sh)

    # The tick is remembered across views but applies to neither a phase nor a geometry, so
    # switching to one cannot ask for log10 of a negative.
    G.shell_set_obs_view("t3phi", "baseline", "1", "1")
    @test sh.obslog[] && !G._obs_logscale(sh)
    G.shell_set_obs_view("uv", "baseline", "1", "1")
    @test !G._obs_logscale(sh)
    G.shell_set_obs_view("t3amp", "baseline", "1", "1")
    @test G._obs_logscale(sh)

    # The rest needs the real canvas, which is OITOOLS' — `fresh_shell` has none.
    c, pts = G.build_obs_canvas(Makie.Figure())
    if c === nothing
        @info "OITOOLSGUIExt is not loaded; the log axis itself is untested"
    else
        sh.obs = c
        sh.obsmodel = pts
        for k in ("v2", "t3amp", "t3phi", "uv"), lg in ("0", "1")
            G.shell_set_obs_view(k, "baseline", "1", lg)
        end
        @test !any(l -> occursin("could not draw", l), sh.console)

        G.shell_set_obs_view("v2", "baseline", "1", "1")
        @test sh.obs.axis.yscale[] === log10
        # The overlay is transformed by the same axis, so a model null at exactly zero is
        # dropped with the data's non-positive points rather than reaching log10.
        @test all(p -> p[2] > 0, sh.obsmodel[])
        G.shell_set_obs_view("v2", "baseline", "1", "0")
        @test sh.obs.axis.yscale[] === identity
        G.shell_set_obs_view("t3phi", "baseline", "1", "1")
        @test sh.obs.axis.yscale[] === identity     # the tick is on, the view is signed

        # A NOISE-DOMINATED V²: zero and negative points are what real data has, and a log
        # axis on them threw before — invisibly, because QMLMakie swallows it into
        # "exception in render".
        dat = G.current_dataset(sh.session).data[1]
        v2 = copy(dat.v2)
        dat.v2[1] = -0.02; dat.v2[2] = 0.0
        G.shell_set_obs_view("v2", "baseline", "1", "1")
        @test sh.obs.axis.yscale[] === log10
        @test !any(l -> occursin("could not draw", l), sh.console)

        # Nothing positive at all is reported rather than thrown, and the axis is left in a
        # state the next view can use.
        dat.v2 .= -1.0
        @test G.shell_set_obs_view("v2", "baseline", "1", "1") isa String
        @test any(l -> occursin("could not draw", l), sh.console)
        dat.v2 .= v2
        @test G.shell_set_obs_view("v2", "baseline", "1", "1") isa String
    end
end

@testset "the command log reproduces the session" begin
    sh = fresh_shell()
    G.shell_open_many(join(LAM[1:2], "\n"), "0")
    G.shell_add_model(2)
    G.shell_reconstruct(3, "radflat:100.0:6;sobel2:10.0:0;mem:1.0:0", 5)
    drain!()

    src = G.export_script(sh.session)
    @test occursin("using ROTIR", src)
    @test occursin("readoifits_multiepochs", src)
    @test occursin("image_reconstruct_oi(x0, data, stars", src)   # data BEFORE stars
    @test occursin("radflat_bins", src)
    @test occursin("sobel_gradient_healpix", src)
    # It has to PARSE. The first version put the `bins = …` construction in a trailing comment
    # inside the array literal, which swallowed the following comma.
    @test Meta.parseall(src) isa Expr

    # AND IT HAS TO RUN. Parsing only proves the text is Julia; the claim the command log
    # makes is that the script REPRODUCES the session, and nothing checked that a name it
    # binds is a name it later uses, that the argument order is right, or that a regulariser
    # it constructs is one the reconstruction accepts. Every one of those has been wrong here
    # at least once — `image_reconstruct_oi(x0, stars, data)` shipped with the arguments
    # swapped, and no test above would have caught it.
    #
    # Run in a MODULE of its own so the script's bindings (`data`, `stars`, `x0`, `x`) cannot
    # collide with the test file's, and with `maxiter` cut to keep it seconds rather than
    # minutes — the point is that the path executes, not that it converges.
    sandbox = Module(:CommandLogSandbox)
    Core.eval(sandbox, :(using ROTIR))
    runnable = replace(src, "maxiter = 5" => "maxiter = 2")
    @test (Core.eval(sandbox, Meta.parseall(runnable)); true)
    # The script's own final map, not the session's: this is what a reader would get.
    xs = Core.eval(sandbox, :x)
    @test length(xs) == 12 * (2^3)^2          # a level-3 HEALPix map
    @test all(isfinite, xs)

    # One dataset line, however many files were opened — and it survives a close.
    @test count(e -> e.binding == "data", sh.session.log) == 1
    G.shell_close_dataset()
    @test count(e -> e.binding == "data", sh.session.log) == 0
end

@testset "the current epoch is Julia's to decide" begin
    # The Data tab's table READS this; it used to remember its own row. Adding an epoch moves
    # the session onto it, and the two disagreeing is what made the plot and the table show
    # different nights.
    sh = fresh_shell()
    G.shell_open(LAM[1], "0")
    @test G.shell_current_epoch() == "1"
    G.shell_open_many(LAM[2], "1", "1")             # add, split
    @test length(G.current_dataset(sh.session).data) > 1
    # ONE epoch, shared by every tab. The Model tab's arrows and the Data tab's table are two
    # views of it, and the table drifting off it is what drew one night beside another's
    # numbers — so what the table must do is read this back, not remember its own row.
    G.shell_select_epoch(2)
    @test G.shell_current_epoch() == "2"
    G.shell_select_epoch(1)
    @test G.shell_current_epoch() == "1"
    # Out of range clamps rather than throwing: the table can ask for a row that has just gone.
    G.shell_select_epoch(99)
    @test parse(Int, G.shell_current_epoch()) == length(G.current_dataset(sh.session).data)
    G.shell_close_dataset()
    @test G.shell_current_epoch() == "1"
end

@testset "a running job reports as it goes" begin
    # Before this, a twenty-minute reconstruction showed a spinner and a scrolling trace: no
    # count, no criterion, and the map only at the end. The engine now reports through the
    # job's slot and `shell_job_poll` — which already runs on the GUI thread every 200 ms —
    # draws it. Nothing is drawn on the worker; that is the whole arrangement.
    sh = fresh_shell()
    G.shell_open(LAM[1], "0")
    G.shell_add_model(2)
    # Idle: four fields, and the last one empty.
    idle = cols(G.shell_job_poll())
    @test length(idle) == 4 && idle[1] == "0" && isempty(idle[4])

    G.shell_reconstruct(3, "sobel2:10.0:0", 200)
    reports = String[]
    polys = 0
    for _ in 1:600
        G.shell_job_running() == "1" || break
        sleep(0.05)
        f = cols(G.shell_job_poll())
        if length(f) >= 4 && !isempty(f[4])
            f[4] in reports || push!(reports, f[4])
            polys = max(polys, length(sh.imsky.polys[]))
        end
    end
    G.shell_job_poll()
    @test !isempty(reports)
    @test all(r -> occursin("evaluations", r) && occursin("f = ", r), reports)
    # THE LIVE MAP: polygons on the imaging canvas while the engine was still running, which
    # is the difference between a progress number and watching it converge.
    @test polys > 0
    # The counts advance rather than repeating one report.
    ns = [parse(Int, match(r"^(\d+) evaluations", r).captures[1]) for r in reports]
    @test issorted(ns)
    @test occursin("χ²ᵣ", sh.status)        # and it still finishes normally
end

@testset "a binary keeps its separation across placements" begin
    # Switching to the orbit and back must not move the secondary. It did: the placement
    # setter also took x, y and z, and once those became parameters in their own table the
    # panel stopped passing them — so every call wrote the argument DEFAULTS of zero, put the
    # secondary exactly on the primary, and a binary drew as one star.
    #
    # It produced a plausible picture rather than an error, and only on the way BACK, since
    # switching to the orbit ignores the offset. That is why the round trip is the test rather
    # than a single direction.
    sh = fresh_shell()
    G.shell_add_model(0)
    G.shell_set_binary("1", 0)
    m = G.current_model(sh.session)
    npoly() = (G.refresh_both!(sh); length(sh.msky.polys[]))

    solo = let                       # one star, for the count a binary must exceed
        G.shell_set_binary("0")
        n = npoly()
        G.shell_set_binary("1", 0)
        n
    end
    @test npoly() > solo             # both components drawn from the start

    G.shell_set_position_param("pos_x", "3.0")
    G.shell_set_position_param("pos_y", "1.5")
    @test G.current_model(sh.session).companion.offset == (3.0, 1.5, 0.0)
    both = npoly()
    @test both > solo

    G.shell_set_binary_placement("orbit")
    # The offset SURVIVES being switched away from — a placement is a choice between two
    # sources, and picking one must not destroy the other's numbers.
    @test G.current_model(sh.session).companion.offset == (3.0, 1.5, 0.0)
    @test npoly() == both

    G.shell_set_binary_placement("offset")
    @test G.current_model(sh.session).companion.offset == (3.0, 1.5, 0.0)
    @test npoly() == both
    # And the offset actually places it: the drawn separation is the one that was set.
    got = G.build_epoch_star(sh)
    @test length(got) >= 5
    @test got[5][1] ≈ 3.0 && got[5][2] ≈ 1.5
end

@testset "a binary has ONE orbit, shown once" begin
    # The orbital elements are part of a Roche component's SURFACE definition — its shape
    # follows the instantaneous separation — so the schema lists them on surface type 3 and a
    # binary of two Roche components had them twice: an editable copy in the Primary frame and
    # another in the Secondary frame, free to disagree about the orbit of one pair.
    sh = fresh_shell()
    G.shell_open(LAM[1], "0")
    G.shell_add_model(G.BINARY_CODE)
    m = G.current_model(sh.session)

    groups(t) = unique([cols(r)[9] for r in rows(t)])
    @test !("orbit" in groups(G.shell_params()))        # gone from the primary frame
    @test !("orbit" in groups(G.shell_params2()))       # and from the secondary's
    @test "orbit" in groups(G.shell_binary_orbit_params())

    # A single Roche star is NOT a binary and keeps its elements: there is no Positions frame
    # for them to move to, and they are still what shapes it.
    G.shell_add_model(3)
    @test "orbit" in groups(G.shell_params())
    G.shell_add_model(G.BINARY_CODE)
    m = G.current_model(sh.session)

    # The Roche Binary entry is placed BY THE ORBIT: a Roche shape is computed from the
    # instantaneous separation, so the orbit is already part of the model. (Ticking `secondary`
    # on some other surface type still defaults to a fixed offset — see the schema testset.)
    @test cols(G.shell_binary_placement())[1] == "orbit"
    # The rest of this testset is about the SHARED orbit rows being editable, which they are
    # only when nothing else owns them.
    G.shell_set_binary_placement("offset")
    @test cols(G.shell_binary_placement())[1] == "offset"

    # `q` is the one element whose NUMBER differs between the components — the Roche potential
    # reads M_companion/M_self — so the companion's copy is inverted. Both were 0.5, which is
    # a different system on each side of the same model.
    @test m.params[:q] ≈ 0.5
    @test m.companion.params[:q] ≈ 2.0

    # One edit, both components, each in its own sense — AND the Orbit tab, so the two tabs
    # cannot end up describing different systems. The sync already ran the other way; under a
    # fixed offset nothing owned these, so an element edited here left the orbit plot on the
    # old one.
    @test G.shell_set_binary_orbit_param("a", "2.5") == ""
    @test m.params[:a] ≈ 2.5 && m.companion.params[:a] ≈ 2.5
    @test sh.orbit.params[:a] ≈ 2.5
    @test G.shell_set_binary_orbit_param("Ω", "37") == ""
    @test sh.orbit.params[:Omega] ≈ 37          # ASCII there, Unicode in the schema
    @test m.params[:Ω] ≈ 37
    @test G.shell_set_binary_orbit_param("q", "0.25") == ""
    @test m.params[:q] ≈ 0.25 && m.companion.params[:q] ≈ 4.0
    @test occursin("zero", G.shell_set_binary_orbit_param("q", "0"))

    # Freeing a shared element adds ONE coordinate, on the owner.
    @test G.shell_set_binary_orbit_state("a", "free") == ""
    @test G.shell_free_count() == "1"

    # Under ORBIT placement the Orbit tab owns the elements: the rows show its values, are
    # marked inert, refuse edits, and drop out of the fit — while `q` and `d`, which are not
    # elements, stay editable here.
    G.shell_set_binary_placement("orbit")
    tab = Dict(cols(r)[1] => cols(r) for r in rows(G.shell_binary_orbit_params()))
    @test tab["a"][13] == "1" && tab["P"][13] == "1"
    @test tab["q"][13] == "0" && tab["d"][13] == "0"
    @test parse(Float64, tab["a"][4]) ≈ sh.orbit.params[:a]
    @test parse(Float64, tab["ω"][4]) ≈ sh.orbit.params[:omega]
    @test occursin("Orbit tab", G.shell_set_binary_orbit_param("a", "9"))
    @test occursin("Orbit tab", G.shell_set_binary_orbit_state("a", "free"))
    @test G.shell_free_count() == "0"                   # `a` is the Orbit tab's to fit
    @test G.shell_set_binary_orbit_param("q", "0.3") == ""

    # And the orbit that SHAPES the components is the one that MOVES them: the elements the
    # placement uses are copied into both, so the star cannot be shaped on a 10-day period
    # while being moved on the tab's.
    G.shell_set_orbit_param("P", "4.5")
    @test G.epoch_chi2(sh) !== nothing
    @test m.params[:P] ≈ 4.5 && m.companion.params[:P] ≈ 4.5
end

@testset "the Imaging views honour the view options" begin
    # The three Imaging views were drawn ONCE, when a run finished, straight from the map —
    # so the intensity tick did nothing there while it worked on the Model tab, and neither
    # did a decoration. The visible difference is limb darkening: a temperature map is flat
    # across the disk, the emergent intensity falls towards the limb, and that is what the
    # interferometer actually measures.
    sh = fresh_shell()
    G.shell_open(LAM[1], "0")
    G.shell_add_model(0)

    # Nothing yet: the tab says so rather than showing an empty frame.
    @test sh.lastmap[] === nothing
    G.refresh_image_tab!(sh)
    @test occursin("no reconstruction", sh.imsky.message[])

    got = G.build_epoch_star(sh)
    @test got !== nothing
    star, tmap = got
    G.show_reconstruction!(sh, star, Float64.(tmap))

    # The map is REMEMBERED, which is what makes a later redraw possible at all.
    @test sh.lastmap[] !== nothing
    @test length(sh.lastmap[].x) == length(tmap)

    G.shell_set_surface_field("1", "linear", "0")     # intensity
    ci = copy(sh.imsky.colors[])
    mi = copy(sh.immoll.colors[])
    @test sh.imsky.cbarlabel[] == "I (arb.)"
    G.shell_set_surface_field("0", "linear", "0")     # temperature
    ct = copy(sh.imsky.colors[])
    mt = copy(sh.immoll.colors[])
    @test sh.imsky.cbarlabel[] == "T (K)"

    # The ORTHOGRAPHIC view must actually CHANGE. The map here is a uniform-tpole sphere, so
    # the temperature view is one flat colour and the intensity view is a limb-darkened ramp —
    # if the tick were still ignored these would be identical arrays.
    @test length(ci) == length(ct) > 0
    @test ci != ct
    # And the temperature view of a uniform map IS flat, which is the other half of the
    # check: a difference could otherwise come from anything.
    @test length(unique(ct)) == 1
    @test length(unique(ci)) > 1

    # The MOLLWEIDE has no limb, so the tick reaches it through the intensity LAW alone — and
    # `:linear` intensity of a uniform temperature map is that same uniform map, so identical
    # colours here are the correct answer, not a tick being ignored. This DID differ once, by
    # applying `ldmap` to a whole-surface view, which blanked the far side of every Mollweide
    # and 3-D picture (see `surface_values`).
    @test mi == mt
    @test length(unique(mt)) == 1

    # AND THE LABEL SAYS SO. These two colour bars read "T (K)" as a CONSTANT, so a Planck map
    # — a real surface brightness, not a temperature — was drawn under a temperature label.
    # They follow the field now, by their own rule: `linear` intensity IS the temperature on a
    # view with no limb, so it is still labelled as one.
    G.shell_set_surface_field("1", "linear", "0")
    @test sh.immoll.cbarlabel[] == "T (K)"
    @test sh.imstar === nothing || sh.imstar.cbarlabel[] == "T (K)"
    @test sh.imsky.cbarlabel[] == "I (arb.)"          # the orthographic view HAS a limb
    G.shell_set_surface_field("1", "planck", "1.65")
    @test sh.immoll.cbarlabel[] == "I (arb.)"
    G.shell_set_surface_field("0", "linear", "0")
    @test sh.immoll.cbarlabel[] == "T (K)"

    # Where the whole-surface view does respond is a non-uniform map under Planck, which
    # compresses the cool end and so gives a contrast the temperature does not have.
    G.shell_add_model(2)
    G.shell_set_param("frac_escapevel", "0.9")
    g2 = G.build_epoch_star(sh)
    @test g2 !== nothing
    G.show_reconstruction!(sh, g2[1], Float64.(g2[2]))
    G.shell_set_surface_field("0", "linear", "0")
    mt2 = copy(sh.immoll.colors[])
    G.shell_set_surface_field("1", "planck", "1.65")
    mi2 = copy(sh.immoll.colors[])
    @test length(unique(mt2)) > 1
    @test mi2 != mt2
    # And the whole surface is drawn either way: the far side is part of the map.
    @test !any(c -> Makie.alpha(c) == 0, mt2)
    @test !any(c -> Makie.alpha(c) == 0, mi2)
    G.shell_set_surface_field("0", "linear", "0")

    # A decoration reaches the tab too, by the same route — the graticule polyline is its own
    # Observable, empty until the tick is on.
    @test isempty(sh.imsky.grat[])
    G.shell_set_decoration("graticules", "1")
    @test !isempty(sh.imsky.grat[])
    G.shell_set_decoration("graticules", "0")
    @test isempty(sh.imsky.grat[])
end

@testset "a fit is kept whole" begin
    sh = fresh_shell()
    G.shell_open(LAM[1], "0")
    G.shell_add_model(0)
    G.shell_set_param("radius", "1.35"); G.shell_set_param("tpole", "4800")
    G.shell_set_param_state("radius", "free"); G.shell_set_param_state("ld1", "free")
    G.shell_set_bound("radius", "1.0", "1.8"); G.shell_set_bound("ld1", "0.0", "0.9")
    @test G.shell_fits() == ""
    @test G.shell_current_fit() == "0"

    # An OPTIMISER has no posterior, and the entry has to say so rather than carry an empty
    # matrix that the panel would draw as a spike at the point estimate.
    G.shell_fit("neldermead", 300); drain!()
    frows = rows(G.shell_fits())
    @test length(frows) == 1
    f1 = cols(frows[1])
    @test length(f1) == 9
    @test f1[2] == "neldermead" && f1[6] == "—"                   # no evidence
    # The draws/evals column holds whichever number the METHOD produced: a sampler reports
    # draws, a local optimiser reports how many times it evaluated the criterion. It used to
    # read "0" here, which said nothing — an optimiser has no draws, and the count it does
    # have is the thing that says whether it converged or ran out of budget.
    @test f1[8] != "0" && f1[8] != "—"
    nev = parse(Int, f1[8])
    @test 0 < nev <= 300                                          # NLopt's own count, ≤ budget
    @test G.current_fit(sh.session).nevals == nev
    @test isempty(G.current_fit(sh.session).samples)              # and still no posterior
    G.refresh_posterior!(sh)
    @test occursin("no posterior", sh.post.message[])

    # THE COLUMN ORDER, which is the part that fails silently. `_fit_hmc` and `_fit_pigeons`
    # return the free parameters sorted by `parametric_free_indices`, which is not the order
    # the form lists them in — so a posterior can land under the wrong parameter's name and
    # look entirely reasonable. Distinct per-column values are what make a swap visible.
    n = 40
    S = hcat(fill(1.0, n) .+ 0.01 .* (1:n), fill(9.0, n) .+ 0.01 .* (1:n))
    raw = (samples = S, q16 = [1.0, 9.0], q84 = [1.4, 9.4], logz = -12.5, logzerr = 0.25)
    straight = G._posterior(raw; diagnostics = "40 draws")
    @test straight.samples[1, 1] ≈ 1.01 && straight.samples[1, 2] ≈ 9.01
    swapped = G._posterior(raw; order = invperm([2, 1]))
    @test swapped.samples[1, 1] ≈ 9.01 && swapped.samples[1, 2] ≈ 1.01
    @test swapped.q16 == [9.0, 1.0] && swapped.q84 == [9.4, 1.4]
    @test straight.logz ≈ -12.5 && straight.logzerr ≈ 0.25

    # A sampler's entry, through the same recorder the job path uses.
    res = (names = [:radius, :ld1], best = [1.2, 9.2], errs = Dict(:radius => 0.2, :ld1 => 0.2),
           post = straight, method = :nautilus, model = "sphere_1", surface_type = 0,
           chi2 = 1500.0, ndata = 1000, table = "", params = Dict{Symbol,Float64}(),
           status = "")
    e = G._record_fit!(sh, res)
    @test size(e.samples) == (n, 2)
    @test e.logz ≈ -12.5
    @test length(sh.session.fits) == 2 && G.shell_current_fit() == "2"
    # log(Z) reaches the comparison table, with its error beside it — a difference smaller
    # than the error is not a preference for either model.
    r2 = cols(rows(G.shell_fits())[2])
    @test r2[6] == "-12.500" && r2[7] == "0.250" && r2[8] == "40"
    @test length(rows(G.shell_fit_params())) == 2

    # And it draws: the marginal, the band and the pair.
    G.refresh_posterior!(sh)
    @test isempty(sh.post.message[])
    @test length(sh.post.dens[]) > 4
    @test length(sh.post.pair[]) == n
    @test sh.post.axis1.xlabel[] == "radius"
    G.shell_set_posterior_pair(2, 1)
    @test sh.post.axis1.xlabel[] == "ld1" && sh.post.axis2.ylabel[] == "radius"
    @test G.shell_posterior_pair() == "2\t1"

    # Selecting the optimiser's fit puts the panel back to saying there is nothing to show.
    G.shell_select_fit(1)
    @test occursin("no posterior", sh.post.message[])
end

@testset "saving a view" begin
    sh = fresh_shell()
    dir = mktempdir()
    # Nothing loaded: a message, not an exception, and no file.
    @test startswith(G.shell_save_figure("sky", joinpath(dir, "a.png"), 400, 300), "!")
    @test startswith(G.shell_save_figure("nosuch", joinpath(dir, "b.png"), 400, 300), "!")

    G.shell_open(LAM[1], "0")
    G.shell_add_model(2)
    for w in ("sky", "mollweide", "star3d", "chi2", "orbit")
        p = joinpath(dir, "$(w).png")
        # Rebuilt offscreen through the window's own builders — the live framebuffer cannot be
        # read under QMLMakie at all. See src/gui/snapshot.jl.
        @test G.shell_save_figure(w, p, 500, 400) == ""
        @test isfile(p) && filesize(p) > 5_000
    end
    # An extension is appended when the caller leaves it off, rather than writing a file
    # nothing will open.
    G.shell_save_figure("sky", joinpath(dir, "noext"), 400, 300)
    @test isfile(joinpath(dir, "noext.png"))

    # AND IT LEAVES NO SCREEN BEHIND. `Makie.save` attaches a hidden GLFW window to the figure
    # and registers it; left open it keeps the process alive after the Qt window closes, which
    # showed up as a click-test run that did everything right and then never exited.
    before = length(GLMakie.ALL_SCREENS)
    G.shell_save_figure("chi2", joinpath(dir, "screens.png"), 400, 300)
    @test length(GLMakie.ALL_SCREENS) <= before
end

@testset "the polyft kernel is selectable" begin
    sh = fresh_shell()
    rows_ = rows(G.shell_polyft_backends())
    @test length(rows_) == 3
    # Fastest first: the panel offers them in the order the measurements put them.
    @test [cols(r)[1] for r in rows_] == ["nufft", "turbo", "scalar"]
    @test G.shell_polyft_backend() == "nufft"
    @test occursin("scalar", G.shell_set_polyft_backend("scalar"))
    @test G.shell_polyft_backend() == "scalar"
    @test G.shell_set_polyft_backend("rasterize") ==
          "backend must be nufft, turbo or scalar"
    # ALL THREE agree on a real χ². That is the only thing that makes offering a choice safe:
    # a backend that is fast and slightly wrong would bias every fit run through it.
    G.shell_open(LAM[1], "0"); G.shell_add_model(0)
    c_scalar = G.epoch_chi2(sh)[1].total
    for b in ("turbo", "nufft")
        # The RETURN and the selection, not only the number. This compared χ² alone, and
        # since all three agree, a backend that silently failed to be selected passed the
        # test — which is exactly what happened twice: once when `:turbo` moved into an
        # extension and the launcher manifest had not recorded it, and once because loading
        # LoopVectorization mid-session leaves its methods invisible to the frame that
        # loaded it (world age), so the χ² came back `nothing` and only the selection said
        # anything was wrong.
        #
        # The ORDER here is the regression test for the second: `:turbo` is selected and its
        # χ² taken inside ONE top-level block, so the world age never advances between the
        # load and the call — which is the GUI's situation, where every callback runs in the
        # world age fixed when `QML.exec()` was entered.
        msg = G.shell_set_polyft_backend(b)
        @test occursin(b, msg)
        @test G.shell_polyft_backend() == b
        @test abs(G.epoch_chi2(sh)[1].total - c_scalar) / c_scalar < 1e-4
    end
    # `:turbo` needs LoopVectorization, which the GUI does NOT load at startup — selecting it
    # loads it on demand. By here that has happened, so the extension must be live: if the
    # lazy load had failed, `shell_set_polyft_backend` would have said so above.
    @test ROTIR.turbo_available()
    G.shell_set_polyft_backend("nufft")
end

@testset "the orbit tab" begin
    sh = fresh_shell()

    # The two frames the tab keeps apart. The elements are the orbit; the star model is what
    # sits at the two positions it puts them in.
    els = rows(G.shell_orbit_params())
    @test !isempty(els)
    @test all(length(cols(r)) == 8 for r in els)
    names = [cols(r)[1] for r in els]
    @test "a" in names && "i" in names && "Omega" in names && "P" in names
    # `q`, `rpole` and `tpole` are star-model quantities. A relative astrometric orbit says
    # nothing about a mass ratio, so one appearing among the elements would be a claim the
    # data cannot support.
    @test !("q" in names) && !("rpole1" in names) && !("tpole1" in names)

    mdls = rows(G.shell_orbit_star_models())
    @test length(mdls) == 2
    @test [cols(r)[1] for r in mdls] == ["analytic", "tessellated"]
    @test G.shell_orbit_star_model() == "analytic"

    # Under the analytic model only "show stars" is offered: the other three are surface
    # physics that an analytic profile has no surface for.
    @test [cols(r)[1] for r in rows(G.shell_orbit_options())] == ["render"]
    @test cols(rows(G.shell_orbit_options())[1])[2] == "show stars"

    G.shell_set_orbit_param("a", "3.0"); G.shell_set_orbit_param("P", "10.0")
    G.shell_set_orbit_param("e", "0.3")
    # The ELEMENT values, keyed by name: the row SET changes with the star model — `f` and the
    # per-component profile parameters are analytic-only — so a positional comparison would be
    # comparing two different tables.
    elvals() = Dict(cols(r)[1] => cols(r)[3] for r in rows(G.shell_orbit_params()))
    before = elvals()
    # THE 3-D COMPONENTS ARE THE MODEL TAB'S BINARY, so without one there is nothing to put at
    # the two orbital positions and this refuses rather than falling back to a hardcoded pair.
    @test occursin("define a binary", G.shell_set_orbit_star_model("tessellated"))
    @test G.shell_orbit_star_model() == "analytic"      # and it does not half-switch
    G.shell_add_model(3)
    G.shell_set_binary("1", 3)
    @test G.shell_binary() == "1"
    @test occursin("tessellated", G.shell_set_orbit_star_model("tessellated"))
    @test G.shell_orbit_star_model() == "tessellated"
    # SWITCHING THE STAR MODEL MUST NOT MOVE THE SECONDARY. The elements describe the orbit,
    # not the stars, and this is the assertion that keeps the two frames separate.
    after = elvals()
    @test all(after[k] == v for (k, v) in before if haskey(after, k))
    @test all(k in keys(after) for k in ORBIT_ELEMENTS_NAMES)
    # What DOES go away is the component half: analytic-only, and the tessellated path reads
    # none of it.
    @test !haskey(after, "c1_diameter") && !haskey(after, "f")
    @test [cols(r)[1] for r in rows(G.shell_orbit_options())] ==
          ["render", "roche", "irradiation", "occultation"]
    @test !occursin("dearest", G.shell_orbit_options())
    @test G.shell_set_orbit_star_model("bogus") == "star model must be analytic or tessellated"

    # Both star models draw. The analytic one draws the profiles the fit uses — rings, so a
    # limb-darkened disk is distinguishable from a uniform one — and the tessellated one real
    # surfaces, which is why it produces far more polygons.
    G.shell_set_orbit_star_model("analytic")
    G.shell_set_orbit_option("render", "1")
    G.shell_set_orbit_param("c1_diameter", "0.8"); G.shell_set_orbit_param("c2_diameter", "0.5")
    G.show_orbit!(sh.orbitcanvas, sh.orbit, [0.0, 2.5, 5.0])
    nanalytic = length(sh.orbitcanvas.polys[])
    @test nanalytic == 3 * 2 * 10                      # epochs x components x rings
    @test length(sh.orbitcanvas.polycolors[]) == nanalytic
    @test sh.orbitcanvas.cbarlabel[] == "relative brightness"
    G.shell_set_orbit_star_model("tessellated")
    G.show_orbit!(sh.orbitcanvas, sh.orbit, [0.0, 2.5, 5.0];
                  binary = G._orbit_binary(sh))
    @test length(sh.orbitcanvas.polys[]) > nanalytic
    @test sh.orbitcanvas.cbarlabel[] == "T (K)"
    # And the components carry the MODEL TAB's parameters, not a hardcoded pair: change the
    # companion's polar radius and the bparams the renderer builds must follow.
    m = G.current_model(sh.session)
    m.companion.params[:rpole] = 0.83
    @test G.orbit_bparams(sh.orbit; binary = m).star2.rpole ≈ 0.83
    m.params[:ldtype] = 2.0; m.params[:ld1] = 0.37
    bp = G.orbit_bparams(sh.orbit; binary = m)
    @test bp.star1.ldtype == 2 && bp.star1.ld1 ≈ 0.37
    # Without the binary it falls back to the tab's own thin description, which is different.
    @test G.orbit_bparams(sh.orbit).star1.ldtype == 3
    # Off is off, in both.
    G.shell_set_orbit_option("render", "0")
    G.show_orbit!(sh.orbitcanvas, sh.orbit, [0.0, 2.5])
    @test isempty(sh.orbitcanvas.polys[])

    # A round trip through TOML, star model included.
    G.shell_set_orbit_state("e", "free"); G.shell_set_orbit_tie("omega", "-Omega")
    f = joinpath(mktempdir(), "orbit.toml")
    G.shell_save_orbit(f)
    o2 = G.load_orbit(f)
    @test o2.model === sh.orbit.model
    @test o2.params[:a] ≈ 3.0 && o2.params[:e] ≈ 0.3
    @test :e in o2.free && o2.ties[:omega] == "-Omega"
end

@testset "the orbit orients both components" begin
    # A binary placed BY ITS ORBIT has no orientation of its own to set. `orbit_bparams` has
    # taken the spin axis from the orbital inclination and node since the Orbit tab was
    # written — `create_binary_geometry` does the same — while the Model tab passed each
    # component's OWN `inclination`/`position_angle` to `create_star_multiepochs`, so the two
    # halves of one picture disagreed. `graticule_segments` records the same trap from the
    # other side, measuring the answers 102.6 degrees apart on a Spica-like binary.
    sh = fresh_shell()
    G.shell_open(LAM[1], "0")
    G.shell_add_model(G.BINARY_CODE)
    m = G.current_model(sh.session)
    @test m.companion.place === :orbit

    G.shell_set_orbit_param("i", "116")
    G.shell_set_orbit_param("Omega", "310")
    tbl(f) = Dict(cols(r)[1] => cols(r) for r in rows(f()))

    p1 = tbl(G.shell_params)
    @test parse(Float64, p1["inclination"][4]) ≈ 180 - 116
    @test parse(Float64, p1["position_angle"][4]) ≈ 310 - 180
    @test p1["inclination"][13] == "1" && p1["position_angle"][13] == "1"   # inert
    @test occursin("from the orbit", p1["inclination"][12])                 # and it says why

    # BOTH components, not only the primary.
    p2 = tbl(G.shell_params2)
    @test parse(Float64, p2["inclination"][4]) ≈ 180 - 116
    @test parse(Float64, p2["position_angle"][4]) ≈ 310 - 180
    @test p2["inclination"][13] == "1"
    @test m.params[:inclination] ≈ 64.0 && m.companion.params[:inclination] ≈ 64.0
    @test m.params[:position_angle] ≈ 130.0 && m.companion.params[:position_angle] ≈ 130.0

    # A derived value is not editable and is not a fit coordinate: the next build overwrites
    # it, so a fit that moved it would be moving a direction the chi-squared only appears to
    # respond to.
    @test occursin("orbit orients", G.shell_set_param("inclination", "20"))
    @test occursin("orbit orients", G.shell_set_param_state("position_angle", "free"))
    @test occursin("orbit orients", G.shell_set_param2("inclination", "20"))
    @test m.params[:inclination] ≈ 64.0
    @test G.shell_free_count() == "0"

    # Under a FIXED OFFSET the orientation is the component's own again — there is no orbit
    # deciding where it points, so nothing else owns it.
    G.shell_set_binary_placement("offset")
    q = tbl(G.shell_params)
    @test q["inclination"][13] == "0"
    @test G.shell_set_param("inclination", "20") == ""
    @test m.params[:inclination] ≈ 20.0
    @test G.shell_set_param_state("inclination", "free") == ""
    @test G.shell_free_count() == "1"

    # And going back to the orbit takes it away again, free set included.
    G.shell_set_binary_placement("orbit")
    G.shell_params()                       # the sync runs where the rows are built
    @test m.params[:inclination] ≈ 64.0
    @test !(:inclination in m.free)
    @test G.shell_free_count() == "0"
end

@testset "the Orbit tab's component parameters are 2-D only" begin
    # `c1_diameter` and its siblings come from `orbit_param_names`, which appends them for
    # whichever component KINDS are selected — so they landed in the same flat table as the
    # elements and were drawn whatever the star model was. The tessellated path never reads
    # them: it sizes its components from `rpole1`/`rpole2` and takes their relative brightness
    # from the temperatures.
    sh = fresh_shell()
    names() = [cols(r)[1] for r in rows(G.shell_orbit_params())]
    @test G.shell_orbit_star_model() == "analytic"
    @test "c1_diameter" in names() && "c2_diameter" in names() && "f" in names()

    # The 3-D components ARE the Model tab's binary, so one has to exist first.
    G.shell_add_model(G.BINARY_CODE)
    @test occursin("tessellated", G.shell_set_orbit_star_model("tessellated"))
    n2 = names()
    @test !("c1_diameter" in n2) && !("c2_diameter" in n2) && !("f" in n2)
    @test "a" in n2 && "i" in n2 && "Omega" in n2 && "P" in n2     # the elements stay

    # UNLISTED, not deleted: switching back brings them as they were, which is what makes the
    # two star models comparable at all.
    @test haskey(sh.orbit.params, :c1_diameter)
    G.shell_set_orbit_star_model("analytic")
    @test "c1_diameter" in names()
end

@testset "the epoch marks fall on the orbit" begin
    # They did not, on Spica, and the cause was a UNIT: `refresh_orbit!` handed the dataset's
    # MJDs straight to `orbit_to_rotir_offset`, which reads an epoch in the same JD the
    # elements are quoted in — `companion_offsets` and `fit_orbit`'s `_uv_times` both add
    # 2400000.5 for exactly that reason.
    #
    # A wrong time alone would only put a mark at the wrong PHASE of the right ellipse, still
    # on the curve. What took it OFF the curve is apsidal motion: `omega_at` advances ω by
    # `dω·(t − T0)`, and Spica ships ω̇ = 0.0071 deg/day, so being 2.4 million days early
    # rotated the apsidal line by ~96 degrees against the track. Spica is the test for that
    # reason — a `dω` of zero cannot show the bug at all.
    sh = fresh_shell()
    G.shell_open(LAM[1], "0")
    @test occursin("loaded spica",
                   G.shell_load_orbit(joinpath(pkgdir(ROTIR), "demos", "orbits",
                                               "spica.toml")))
    @test sh.orbit.params[:domega] > 0                    # the orbit that can show it
    G.refresh_orbit!(sh)

    trk = sh.orbitcanvas.track[]
    marks = sh.orbitcanvas.marks[]
    @test length(trk) > 100
    @test length(marks) == length(rows(G.shell_epochs()))

    # Every mark within a track-sampling step of the curve. 360 points over one period puts
    # neighbouring samples about 1.7 % of the radius apart, so 3 % is the sampling itself and
    # nothing more; the bug put them of order the semi-major axis away.
    scale = maximum(hypot(p[1], p[2]) for p in trk)
    @test scale > 0
    worst = maximum(minimum(hypot(mk[1] - p[1], mk[2] - p[2]) for p in trk) for mk in marks)
    @test worst < 0.03 * scale

    # And the track is the one the marks live on: anchored to the period containing the middle
    # of the data, not to a `T0` that can be years away from it.
    bp = G.orbit_bparams(sh.orbit)
    @test G._track_epoch(bp, Float64[]) ≈ bp.T0
    tep = G.current_dataset(sh.session).mjd .+ 2_400_000.5
    @test abs(G._track_epoch(bp, tep) - sum(tep) / length(tep)) <= bp.P
end

@testset "the analytic profiles are whole, and the size of the 3-D ones" begin
    # A SLICE was missing from every analytic disc. An annulus is drawn as one polygon with a
    # slit — outer arc, then inner arc reversed — and both arcs were sampled over the OPEN
    # range 0 to 2π-Δ, so the polygon jumped inwards at 2π-Δ and left a Δ-wide wedge out. At
    # `nang = 72` that is 5 degrees missing from nine of the ten rings; only the central disc,
    # which is a filled polygon and closes itself, came out whole.
    nang = 72
    polys, vals = G._profile_rings(:uniform, Dict(:c1_diameter => 1.0), :c1, 0.0, 0.0;
                                  nring = 10, nang = nang)
    @test length(polys) == 10 == length(vals)
    R = 0.5
    edges = range(0, R, length = 11)
    for k in 2:10                                # every annulus, not just the outermost
        rin, rout = edges[k], edges[k + 1]
        ring = polys[k]
        # The OUTER arc of this ring: the half of its points beyond its own mid-radius.
        angs = sort([atan(p[2], p[1]) for p in ring
                     if hypot(p[1], p[2]) > (rin + rout) / 2])
        @test length(angs) >= nang               # the 2π endpoint is there
        gaps = diff(vcat(angs, angs[1] + 2π))
        @test maximum(gaps) <= 2π / nang + 1e-4  # one sampling step, not two
    end

    # And the two star models draw the same STAR. The tessellated polar radii are rendering
    # knobs, independent of the analytic component sizes that the fit moves, and both shipped
    # orbits carried `default_orbit`'s old 0.2 and 0.12 — so Spica's 3-D surfaces came out at
    # less than half the size of its 0.894 mas analytic primary.
    o = G.load_orbit(joinpath(pkgdir(ROTIR), "demos", "orbits", "spica.toml"))
    @test o.rpole1 ≈ o.params[:c1_diameter] / 2
    @test o.rpole2 ≈ o.params[:c2_diameter] / 2
    d = G.default_orbit()
    @test d.rpole1 ≈ d.params[:c1_diameter] / 2
    @test d.rpole2 ≈ d.params[:c2_diameter] / 2
    # A file that gives a radius still wins: a component can render at a size its analytic
    # profile does not have.
    f = joinpath(mktempdir(), "explicit.toml")
    open(f, "w") do io
        println(io, "[parameters]\nc1_diameter = 1.0\nc2_diameter = 1.0")
        println(io, "[rendering]\nrpole1 = 0.07")
    end
    o2 = G.load_orbit(f)
    @test o2.rpole1 ≈ 0.07
    @test o2.rpole2 ≈ 0.5
end

@testset "the spin axis emerges from the pole that is visible" begin
    # It was ONE line, tip to tip, drawn after the polygons and so straight across the disk.
    # That reads as an axis passing THROUGH the star, and it draws the hidden pole's half as
    # confidently as the visible one — which is also why the axis looked as though its pole
    # were in the wrong place: what the eye takes for the pole is where the line crosses the
    # limb, not where the pole is.
    #
    # MEASURED, on the two surface types whose mesh actually carries an orientation: the
    # axis's pole and the graticule's are the SAME to 0.000 degrees. `_spin_axis`'s mesh
    # branch, `_spin_axis`'s analytic branch and `_mesh_rotation`'s `R[3,:]` all agree at
    # inclinations 60, 30 and 120 on a rapid rotator and on a Roche component. There was
    # never a geometry error to fix here; see the note in plans/gui_todo.md.
    tess = tessellation_healpix(3)
    function stub(inc)
        p = merge(default_star_params(2),
                  (inclination = inc, position_angle = 0.0))
        star = create_star(tess, p, 0.0)
        n, s = ROTIR._spin_axis(star, p, NaN, NaN)
        return (G._axis_polyline(star, p, 0.0, 0.0), n, s)
    end

    # Below 90 the NORTH pole faces us: one stub, starting AT the north pole — not at the
    # south, and not at the centre — and running outward from it.
    pts, north, south = stub(60.0)
    @test north[3] > 0 && south[3] < 0
    @test length(pts) == 2
    @test !any(q -> isnan(q[1]), pts)
    @test pts[1][1] ≈ -north[1] atol = 1e-4          # the plot negates west
    @test pts[1][2] ≈ north[2] atol = 1e-4
    @test hypot(pts[2]...) > hypot(pts[1]...)

    # Past 90 the south pole is the visible one and the stub swaps ends, with no rule to keep
    # in step with the geometry: the test is the pole's own sky z.
    pts2, n2, s2 = stub(120.0)
    @test n2[3] < 0 && s2[3] > 0
    @test length(pts2) == 2
    @test pts2[1][1] ≈ -s2[1] atol = 1e-4
    @test pts2[1][2] ≈ s2[2] atol = 1e-4
    @test hypot(pts2[2]...) > hypot(pts2[1]...)

    # Never both, away from edge-on. (AT 90 degrees both poles sit on the limb and which one
    # the mesh calls nearer is Float32 noise, so that case is not asserted either way — the
    # NaN break exists to draw two pieces from one polyline if it comes up.)
    for inc in (5.0, 45.0, 89.0, 91.0, 135.0, 175.0)
        q, _, _ = stub(inc)
        @test length(q) == 2
    end

    # And the tick reaches the canvas the Model tab draws on.
    sh = fresh_shell()
    G.shell_open(LAM[1], "0")
    G.shell_add_model(2)
    G.shell_set_decoration("spin", "1")
    @test length(sh.msky.axis3d[]) == 2
    G.shell_set_decoration("spin", "0")
    @test isempty(sh.msky.axis3d[])

    # HEAVIER THAN THE GRATICULE, AND OVER IT. Two thin black lines are one line to the eye —
    # a meridian passing near the pole read as the axis. And "drawn after" only means anything
    # if the axis overdraws too: the graticule must overdraw (or z-fighting with the tessels
    # renders it dashed), and an overdrawn line beats an ordinary one whatever the insertion
    # order. With both overdrawn, order decides and the axis is created later.
    @test sh.msky.axisplot.linewidth[] > sh.msky.gratplot.linewidth[]
    @test sh.msky.spinplot.linewidth[] > sh.msky.gratplot.linewidth[]
    @test sh.msky.axisplot.overdraw[] && sh.msky.spinplot.overdraw[]
    @test sh.msky.gratplot.overdraw[]
end

@testset "decorations cover both components of a binary" begin
    # They followed the PRIMARY alone. A limb around one star and not the other reads as the
    # second having no edge, and a spin axis on one says the companion does not rotate — and a
    # binary is exactly where the comparison between the two is the point.
    sh = fresh_shell()
    G.shell_open(LAM[1], "0")
    G.shell_add_model(3)                       # ONE Roche star, for the baseline
    for d in ("limb", "graticules", "spin"); G.shell_set_decoration(d, "1"); end
    G.refresh_both!(sh)
    solo = (limb = length(sh.msky.limb[]), grat = length(sh.msky.grat[]),
            axis = length(sh.msky.axis3d[]), spin = length(sh.msky.spin[]))
    @test all(v -> v > 0, values(solo))

    G.shell_set_binary("1", 3)
    G.shell_set_binary_placement("offset")     # each component keeps its own orientation
    G.shell_set_position_param("pos_x", "4.0")
    G.refresh_both!(sh)
    both = (limb = length(sh.msky.limb[]), grat = length(sh.msky.grat[]),
            axis = length(sh.msky.axis3d[]), spin = length(sh.msky.spin[]))
    for k in keys(solo)
        @test both[k] > solo[k]                # the companion's half is actually there
    end

    # STILL ONE POLYLINE each, with the pen lifted between the components: this canvas cannot
    # insert a plot after the window exists, so a NaN break is how two pieces are drawn.
    breaks(v) = count(p -> isnan(p[1]), v)
    @test breaks(sh.msky.limb[]) >= 1
    @test breaks(sh.msky.axis3d[]) >= 1

    # And `_break` itself: either side empty joins nothing and leaves no stray NaN at an end.
    a = [Makie.Point2f(0, 0), Makie.Point2f(1, 1)]
    b = [Makie.Point2f(2, 2)]
    @test length(G._break(a, b)) == 4 && breaks(G._break(a, b)) == 1
    @test G._break(a, Makie.Point2f[]) == a
    @test G._break(Makie.Point2f[], b) == b
    @test isempty(G._break(Makie.Point2f[], Makie.Point2f[]))
end

@testset "the Orbit tab's time cursor" begin
    # VISUALIZATION ONLY: the numbered marks stay the observations', and what the cursor
    # changes is which single time the surfaces are rendered at. Nothing here reaches the
    # dataset, the model or a chi-squared.
    sh = fresh_shell()
    t = cols(G.shell_times())
    @test length(t) == 7
    @test t[1] == "0"                          # off by default
    @test G.cursor_times(sh) === nothing
    @test G.cursor_time(sh) === nothing

    # The resolved defaults show even while it is off — one period from periastron, in 48
    # steps, which is the span `_orbit_track_2d` already samples.
    v = G.apply_orbit_ties(sh.orbit)
    @test parse(Float64, t[2]) ≈ v[:T0]
    @test parse(Float64, t[3]) ≈ v[:T0] + v[:P]
    @test parse(Int, t[6]) == 49

    @test occursin("frames", G.shell_set_times("1", "", "", ""))
    ts = G.cursor_times(sh)
    @test ts !== nothing && length(ts) == 49
    @test G.cursor_time(sh) ≈ first(ts)
    @test length(sh.orbitcanvas.cursor[]) == 1
    @test length(sh.orbitcanvas.cursorlabel[]) == 1
    @test occursin("JD", sh.orbitcanvas.cursorlabel[][1])
    @test isempty(sh.orbitcanvas.marks[])      # no dataset, so no observed epochs

    # Stepping CLAMPS rather than wraps: an orbit is periodic but a time range is not, and
    # jumping from the last frame back to the first would hide that the end was reached.
    G.shell_step_time(1)
    @test sh.times.index == 2
    G.shell_step_time(-5)
    @test sh.times.index == 1
    G.shell_set_time_index(999)
    @test sh.times.index == 49
    @test G.cursor_time(sh) ≈ last(ts)

    # An explicit range, and it tracks the elements rather than being frozen at the defaults.
    @test occursin("frames", G.shell_set_times("1", "2450000", "2450010", "2.5"))
    @test length(G.cursor_times(sh)) == 5
    G.shell_set_times("1", "", "", "")
    G.shell_set_orbit_param("P", "20")
    @test length(G.cursor_times(sh)) == 49
    @test last(G.cursor_times(sh)) ≈ v[:T0] + 20

    # Unticking takes the cursor off the plot entirely.
    @test occursin("times off", G.shell_set_times("0", "", "", ""))
    @test G.cursor_time(sh) === nothing
    @test isempty(sh.orbitcanvas.cursor[])
    @test isempty(sh.orbitcanvas.cursorlabel[])
end

@testset "no callback throws on the shapes QML sends" begin
    # THE FREEZE CLASS. An exception thrown inside a Julia callback escapes through
    # `QML.julia_call` and stops the whole window responding — no error on screen, nothing in
    # the console, just a dead GUI. So for these entry points "does not throw" is the property
    # that matters, ahead of what they return.
    #
    # The shapes are the ones QML really sends: a TextField gives a String, a SpinBox an
    # Int32, a Slider a Float64, and a RECYCLED DELEGATE gives `undefined`, which arrives as
    # `nothing` — `onEditingFinished` fires on focus loss and can run after its row's model
    # context is gone. `String(::Int32)`, `String(::Float64)` and `String(::Nothing)` all have
    # no method, and every one of those was a reachable freeze.
    shapes = Any["1", "0", "", "abc", Int32(1), Int32(-1), 0.0, 1.0, -5.0, nothing]

    sh = fresh_shell()
    G.shell_open(LAM[1], "0")
    G.shell_add_model(G.BINARY_CODE)

    # `f(args...)` must return, whatever it returns. A refusal is fine; a throw is not.
    function survives(f, args...)
        try
            f(args...)
            return true
        catch err
            @error "callback threw — this freezes the window" f args err
            return false
        end
    end

    for a in shapes
        # name-taking setters: the recycled-delegate case is `nothing` here
        @test survives(G.shell_set_param, a, "1.0")
        @test survives(G.shell_set_param_state, a, "free")
        @test survives(G.shell_set_bound, a, "0", "1")
        @test survives(G.shell_set_tie, a, "x")
        @test survives(G.shell_set_param2, a, "1.0")
        @test survives(G.shell_set_param_state2, a, "free")
        @test survives(G.shell_set_bound2, a, "0", "1")
        @test survives(G.shell_set_tie2, a, "x")
        @test survives(G.shell_set_position_param, a, "1.0")
        @test survives(G.shell_set_position_state, a, "free")
        @test survives(G.shell_set_binary_orbit_param, a, "1.0")
        @test survives(G.shell_set_binary_orbit_state, a, "free")
        @test survives(G.shell_set_binary_orbit_bound, a, "0", "1")
        @test survives(G.shell_set_binary_orbit_tie, a, "x")
        @test survives(G.shell_set_orbit_param, a, "1.0")
        @test survives(G.shell_set_orbit_state, a, "free")
        @test survives(G.shell_set_orbit_bound, a, "0", "1")
        @test survives(G.shell_set_orbit_tie, a, "x")
        @test survives(G.shell_set_orbit_render_param, a, "1.0")
        @test survives(G.shell_set_decoration, a, "1")
        # and the VALUE side of the same calls
        @test survives(G.shell_set_param, "rpole", a)
        @test survives(G.shell_set_binary_orbit_param, "a", a)
        @test survives(G.shell_set_orbit_param, "a", a)
        @test survives(G.shell_set_binary_orbit_state, "a", a)
        @test survives(G.shell_set_binary_orbit_bound, "a", a, a)
        @test survives(G.shell_set_binary_orbit_tie, "a", a)
        # flags and choices, which arrive from ticks and combos
        @test survives(G.shell_set_decoration, "limb", a)
        @test survives(G.shell_set_binary, a, 3)
        @test survives(G.shell_set_binary_placement, a)
        @test survives(G.shell_set_surface_field, a, "linear", "0")
        @test survives(G.shell_set_surface_field, "1", a, "0")
        @test survives(G.shell_set_orbit_option, a, "1")
        @test survives(G.shell_set_orbit_option, "render", a)
        @test survives(G.shell_set_orbit_star_model, a)
        @test survives(G.shell_set_orbit_component, a, "uniform")
        @test survives(G.shell_set_times, a, a, a, a)
        @test survives(G.shell_set_time_index, a)
        @test survives(G.shell_step_time, a)
        @test survives(G.shell_set_graticule, a, a, "black")
        @test survives(G.shell_set_colormap, a)
    end

    # And the READERS, which QML calls on every refresh: one of these throwing freezes the
    # window on a tab switch rather than on an edit, which is harder to place.
    for f in (G.shell_params, G.shell_params2, G.shell_position_params,
              G.shell_binary_orbit_params, G.shell_orbit_params, G.shell_times,
              G.shell_binary_placement, G.shell_free_count, G.shell_surface_field,
              G.shell_epochs, G.shell_models, G.shell_validate_model)
        @test survives(f)
    end
end

@testset "the surface map as a file" begin
    sh = fresh_shell()
    # Nothing to save is a message, not an exception.
    @test occursin("nothing to save", G.shell_save_map(joinpath(mktempdir(), "x.fits")))

    G.shell_open(LAM[1], "0")
    G.shell_add_model(2)
    G.shell_set_param("rpole", "1.37"); G.shell_set_param("tpole", "4800")
    G.shell_set_param("inclination", "78"); G.shell_set_param("rotation_period", "54.8")
    f = joinpath(mktempdir(), "map.fits")
    @test occursin("768 tessels", G.shell_save_map(f))
    @test isfile(f)

    # The point of the file: the χ² is reproducible from it alone, with no session.
    m = load_surface_map(f)
    @test m.nside_exp == 3 && m.tessellation === :healpix && m.field === :temperature
    @test m.params.surface_type === 2 && m.params.surface_type isa Int
    @test m.params.rpole ≈ 1.37 && m.params.tpole ≈ 4800
    @test m.tepochs !== nothing && length(m.tepochs) == 1
    stars = create_star_multiepochs(tessellation_healpix(m.nside_exp), m.params, m.tepochs)
    data  = readoifits(LAM[1])[1, 1]
    setup_oi!([data], stars)
    b = chi2_breakdown(m.x, stars[1], data)
    @test isfinite(b.total) && b.total > 0 && b.ndata > 0

    @test occursin("768 tessels", G.shell_load_map(f))
    @test length(sh.session.models) == 1            # the saved parameters came back as a model
    @test length(sh.session.images) == 1            # and the values as an image
    @test sh.session.models[end].params[:rpole] ≈ 1.37
    @test occursin("could not read", G.shell_load_map(joinpath(DATA, "no_such.fits")))

    # An extension appended when the user leaves it off — a map written as "map" is a file
    # nothing will open.
    g = joinpath(mktempdir(), "noext")
    G.shell_save_map(g)
    @test isfile(g * ".fits")
end

@testset "the model geometry as a file" begin
    sh = fresh_shell()
    @test occursin("no model", G.shell_save_geometry(joinpath(mktempdir(), "x.fits")))

    G.shell_open(LAM[1], "0")
    G.shell_add_model(G.BINARY_CODE)
    # A FIXED OFFSET for this one: the Roche Binary entry is placed by the orbit, and under the
    # orbit the elements are the Orbit tab's and the offset is not what decides the position —
    # neither of which this test is about. It is about the file carrying both back.
    G.shell_set_binary_placement("offset")
    G.shell_set_param("rpole", "0.7")
    G.shell_set_binary_orbit_param("a", "4.2")
    G.shell_set_position_param("pos_x", "3.0")
    f = joinpath(mktempdir(), "geom.fits")
    st = G.shell_save_geometry(f)
    @test occursin("wrote", st)
    @test occursin("secondary", st)                 # both components are in the file
    @test isfile(f)

    # A fresh session, loaded from the file alone.
    sh2 = fresh_shell()
    st2 = G.shell_load_geometry(f)
    @test occursin("loaded", st2) && occursin("binary", st2)
    m = G.current_model(sh2.session)
    @test m !== nothing
    @test m.surface_type == 3
    @test m.params[:rpole] ≈ 0.7
    @test m.params[:a] ≈ 4.2
    @test m.companion !== nothing
    @test m.companion.surface_type == 3
    @test m.companion.place === :offset
    @test m.companion.offset[1] ≈ 3.0
    # `q` still inverted on the secondary, which is the convention the file carried.
    @test m.companion.params[:q] ≈ 1 / m.params[:q]

    # The point of storing the MESH: the rebuild is measured, not assumed. A Roche surface is
    # a root solve, so the honest answer is Float32 round-off rather than zero.
    msg = filter(l -> occursin("geometry:", l), sh2.console)
    @test length(msg) == 1
    @test occursin("identical", msg[1]) || occursin("differs", msg[1])

    # A single star writes no companion, and loading it leaves none behind.
    G.shell_add_model(0)
    f2 = joinpath(mktempdir(), "solo.fits")
    @test occursin("wrote", G.shell_save_geometry(f2))
    fresh_shell()
    G.shell_load_geometry(f2)
    @test G.current_model(G.SHELL[].session).companion === nothing

    # Not a geometry file at all: a message, not an exception.
    bad = joinpath(mktempdir(), "nope.fits")
    write(bad, "not fits")
    @test occursin("could not read", G.shell_load_geometry(bad))
end

@testset "the file picker" begin
    @test !isempty(G.picker_places())
    @test all(length(cols(r)) == 2 for r in rows(G.picker_places()))
    listing = rows(G.picker_list(DATA))
    @test !isempty(listing)
    @test all(length(cols(r)) == 3 for r in listing)
    @test any(occursin(".oifits", cols(r)[2]) for r in listing)
    # Only the interesting extensions unless asked.
    @test length(rows(G.picker_list(DATA, "1"))) >= length(listing)
    @test G.picker_parent("/") == "/"                    # never an empty string
    @test G.picker_start("file://" * DATA) == DATA
    @test G.picker_start("/no/such/place") == pwd()
    # One picker opened for three purposes, each with its own extensions: hunting an orbit
    # TOML among forty OIFITS files is the listing this avoids.
    @test !occursin(".toml", G.picker_list(DATA, "0", "data"))
    @test all(occursin(".toml", cols(r)[2]) || cols(r)[1] == "dir"
              for r in rows(G.picker_list(DATA, "0", "orbit")))
    # `.oifits` does not end in `.fits`, so the map listing excludes the data files and
    # shows a saved map — which has to be written somewhere to be listed.
    mapdir = mktempdir(); touch(joinpath(mapdir, "a_map.fits"))
    touch(joinpath(mapdir, "obs.oifits"))
    maplist = rows(G.picker_list(mapdir, "0", "map"))
    @test [cols(r)[2] for r in maplist] == ["a_map.fits"]
end

end
