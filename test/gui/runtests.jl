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
    sh = G.ShellState(s,
        G.build_sky_canvas(figs[1]), G.build_star_canvas(figs[2]),
        G.build_moll_canvas(figs[3]), G.build_chi2_canvas(figs[4]),
        G.build_sky_canvas(figs[5]), G.build_moll_canvas(figs[6]),
        G.build_sky_canvas(figs[7]),
        # No observable canvas: it needs OITOOLSGUIExt's `build_canvas`, and `refresh_obs!`
        # is a no-op without it. The observable plot is covered by gui_click.sh, which has a
        # window.
        nothing, Makie.Observable(Makie.Point2f[]), Ref(:v2), Ref(:baseline), Ref(true),
        Ref(true), Ref(:linear), Ref(0.0),
        Dict(:limb => true, :compass => true, :graticules => false,
             :spin => false, :plotmesh => false),
        Ref(30.0), Ref(30.0), Ref("black"),
        Ref(:healpix), Ref(3), Ref{DataType}(Float32),
        G.default_orbit(), G.build_orbit_canvas(figs[8]), "",
        String[], "", nothing, :none, "", Dict{Symbol,Any}(),
        Ref{Any}(nothing), Ref{Any}(nothing),
        G.build_post_canvas(figs[9]), Ref(1), Ref(2), Ref(false),
        G.build_star_canvas(figs[10]), Ref{Any}(nothing))
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
    # Every implemented surface type is offered, and nothing else.
    st = rows(G.shell_surface_types())
    @test length(st) == length(SURFACE_TYPE_ORDER)
    @test Set(parse(Int, cols(r)[1]) for r in st) == Set(SURFACE_TYPE_ORDER)

    for code in SURFACE_TYPE_ORDER
        G.shell_add_model(code)
        m = G.current_model(sh.session)
        @test m.surface_type == code
        @test isempty(G.shell_validate_model())      # defaults must build
        pr = rows(G.shell_params())
        @test length(pr) == length(surface_params(code))
        @test all(length(cols(x)) == 12 for x in pr)
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
    @test all(length(cols(k)) == 5 for k in kinds)

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

@testset "reconstruction" begin
    sh = fresh_shell()
    G.shell_open(LAM[1], "0")
    @test occursin("no model", G.shell_reconstruct(3, "sobel2:10:0", 5))
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
    # Per-pixel indexing for the 3-D and Mollweide views, not per visible tessel.
    @test length(G.surface_values(sh, tmap, star; visible_only = false)) == star.npix
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

    # Both views must actually CHANGE. The map here is a uniform-tpole sphere, so the
    # temperature view is one flat colour and the intensity view is a limb-darkened ramp —
    # if the tick were still ignored these would be identical arrays.
    @test length(ci) == length(ct) > 0
    @test ci != ct
    @test mi != mt
    # And the temperature view of a uniform map IS flat, which is the other half of the
    # check: a difference could otherwise come from anything.
    @test length(unique(ct)) == 1
    @test length(unique(ci)) > 1

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
    @test f1[2] == "neldermead" && f1[6] == "—" && f1[8] == "0"   # no evidence, no draws
    @test isempty(G.current_fit(sh.session).samples)
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
        G.shell_set_polyft_backend(b)
        @test abs(G.epoch_chi2(sh)[1].total - c_scalar) / c_scalar < 1e-4
    end
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
    before = [cols(r)[3] for r in rows(G.shell_orbit_params())]
    @test occursin("tessellated", G.shell_set_orbit_star_model("tessellated"))
    @test G.shell_orbit_star_model() == "tessellated"
    # SWITCHING THE STAR MODEL MUST NOT MOVE THE SECONDARY. The elements describe the orbit,
    # not the stars, and this is the assertion that keeps the two frames separate.
    @test [cols(r)[3] for r in rows(G.shell_orbit_params())] == before
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
    G.show_orbit!(sh.orbitcanvas, sh.orbit, [0.0, 2.5, 5.0])
    @test length(sh.orbitcanvas.polys[]) > nanalytic
    @test sh.orbitcanvas.cbarlabel[] == "T (K)"
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
