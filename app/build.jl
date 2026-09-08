# Build the application bundle.
#
#     xvfb-run -a julia --project=bin app/build.jl          # headless
#     julia --project=bin app/build.jl                      # with a screen
#
# RUN IT WITH --project=bin, not --project=app. PackageCompiler is needed by this script and
# must NOT be a dependency of the application, or it would be bundled into it; `bin/` already
# carries it. What gets built is `app/`, which this script names explicitly.
#
# A DISPLAY IS NEEDED, because the trace opens a GL context and, in ROTIR's case, an actual Qt
# window (see bin/trace.jl). Xvfb is enough.
#
# WHAT create_app DOES NOT DO, and this script therefore does:
#
#   * ship the package's data files. `bundle_project` writes a stub Project.toml and nothing
#     else — no source, no QML, no demo data. They are staged into share/rotir/, mirroring the
#     repository's own relative paths, which is what `ROTIR.resource` expects to find under
#     `<Sys.BINDIR>/../share/rotir`.
#   * seed Makie's font atlas. Without it the first launch tries to DOWNLOAD the atlas and,
#     failing that, re-renders every glyph.
#   * install a launcher. On Linux the windowing system has to be chosen before the process
#     starts; see app/launcher.sh.

using Pkg
using PackageCompiler

const ROOT   = normpath(joinpath(@__DIR__, ".."))
const APPSRC = joinpath(ROOT, "app")
const OUT    = get(ENV, "ROTIR_APP_DIR", joinpath(ROOT, "build", "ROTIR"))
const TRACE  = joinpath(ROOT, "bin", "trace.jl")
const STMTS  = joinpath(APPSRC, "precompile_statements.jl")

isfile(TRACE) || error("missing precompile trace: $TRACE")

# ── no console window on Windows ─────────────────────────────────────────────
#
# PackageCompiler links the C driver as a CONSOLE subsystem executable, so Windows opens a
# terminal and the GUI appears out of it. The subsystem is a two-byte field in the PE header,
# so it is patched on the finished executable rather than passed as a link flag — not via
# `JULIA_CC`, which `get_compiler_cmd` parses with `Base.shell_split`, and which treats a
# backslash as an escape so a Windows compiler path does not survive it.
#
# Every field is checked before the write and the current value must be 3 (console), so a
# layout this does not understand is left alone with a warning rather than corrupted.
"""
    set_windows_gui_subsystem!(exe) -> Bool

Flip a PE executable from the console subsystem to the windows one, in place.

`e_lfanew` at 0x3C gives the PE signature; the optional header follows the 4-byte signature and
the 20-byte COFF header, and `Subsystem` sits at offset 68 within it — the same place in PE32
and PE32+, since everything before it is fixed width in both.
"""
function set_windows_gui_subsystem!(exe::AbstractString)
    isfile(exe) || return false
    open(exe, "r+") do io
        read(io, 2) == b"MZ" || (@warn "not a PE file; console subsystem left as it is" exe; return false)
        seek(io, 0x3C); pe = Int(read(io, UInt32))
        seek(io, pe);   read(io, 4) == b"PE\0\0" || (@warn "no PE signature; left alone" exe; return false)
        opt = pe + 4 + 20
        seek(io, opt); magic = read(io, UInt16)
        magic in (0x10b, 0x20b) || (@warn "unknown optional header magic; left alone" exe magic; return false)
        seek(io, opt + 68); sub = read(io, UInt16)
        sub == 2 && return true                       # already a GUI binary
        sub == 3 || (@warn "unexpected subsystem; left alone" exe subsystem = sub; return false)
        seek(io, opt + 68); write(io, UInt16(2))
        return true
    end
end

@info "Building the application" out = OUT source = APPSRC
@info "This takes tens of minutes and several GB of scratch space."

t = @elapsed create_app(APPSRC, OUT;
                        force = true,
                        precompile_execution_file = TRACE,

                        # A raw `--trace-compile` union, when one has been generated. The
                        # execution trace above already opens the window — unlike OITOOLS',
                        # whose event loop never returns — so the QML bridge and the first
                        # frame are covered by it. What a statements file adds here is the
                        # paths ONE run cannot reach: the other three tabs, the fit engines,
                        # the file picker. Missing is not an error, only a pause on first use.
                        #
                        # Regenerate by unioning traces (see app/precompile_statements.jl):
                        #   xvfb-run -a julia --project=bin --trace-compile=a.trace bin/trace.jl
                        #   xvfb-run -a julia --project=bin --trace-compile=b.trace test/gui/runtests.jl
                        #   cat *.trace | sort -u | grep -vE '\\b(Main|Test)\\.' > app/precompile_statements.jl
                        precompile_statements_file = isfile(STMTS) ? [STMTS] : String[],

                        # THESE TWO GO TOGETHER, and getting the pair wrong is fatal rather
                        # than merely wasteful.
                        #
                        # A JLL compiled into the image runs its `__init__` whether or not
                        # anything calls it. With `include_transitive_dependencies = true` a
                        # JLL that nothing loads is in the image anyway, so omitting its lazy
                        # artifact makes that `__init__` call `find_artifact_dir` on a
                        # directory that was never bundled — and the application dies before
                        # `julia_main`. Invisible on the build machine, where the artifact is
                        # still sitting in ~/.julia.
                        #
                        # `false` is what the option is for: the manual says it "only makes a
                        # difference if some packages do not load all their dependencies when
                        # themselves are loaded". The lazy artifacts can then be left out too.
                        include_transitive_dependencies = false,
                        include_lazy_artifacts = false,

                        # Only the stdlibs this project actually names. The risk the manual
                        # warns about is depending on one WITHOUT naming it — `rand()` needs
                        # Random, `A * B` needs LinearAlgebra and Random both, because those
                        # stdlibs practise type piracy and merely loading them changes
                        # behaviour. A sandbox run is what checks it, since the failure would
                        # be a MethodError at run time rather than anything the build notices.
                        filter_stdlibs = true,

                        # -g0 stops DWARF being GENERATED. Julia's default is -g1 and
                        # PackageCompiler does not override it, which is where a few hundred MB
                        # of .debug_* comes from. It has to be suppressed here because it
                        # cannot be removed afterwards: `strip --strip-debug` produces an image
                        # that then SIGSEGVs.
                        #
                        # Two heavier options are deliberately NOT used.
                        #
                        # --strip-metadata takes more off the serialised heap, but it removes
                        # source locations from backtraces, and `julia_main`'s crash log is the
                        # only diagnostic a user can send back. Not worth a report that says
                        # `none:11`.
                        #
                        # --strip-ir would take more still, and it is the one that breaks this
                        # package: tie expressions in the Model tab are compiled at RUN time
                        # through `eval` (see `eval_tie`), and inlining into that new code
                        # needs the IR of what it calls.
                        sysimage_build_args = `-g0`,

                        # The default is a multiversioned target, which is what lets the binary
                        # run on a CPU other than this one. Do not narrow it: the native code
                        # is a small fraction of the image, so there is little to win here and
                        # portability to lose.
                        )

# ── the shipped resources ────────────────────────────────────────────────────
#
# Repository-relative paths, mirrored. `resource` is then a root substitution and nothing else,
# with no layout mapping to keep in step with the file picker's places.

const SHARE = joinpath(OUT, "share", "rotir")

"""
    treesize(dir) -> Int

Bytes a directory really occupies.

`filesize` FOLLOWS a symlink and reports its target, so summing it over `walkdir` counts every
linked file twice. A bundle is full of versioned `.so` links, which overstates it by more than
half.
"""
treesize(dir) = sum(islink(p) ? 0 : filesize(p)
                    for (r, _, fs) in walkdir(dir) for p in (joinpath(r, f) for f in fs);
                    init = 0)

function stage(rel...)
    src = joinpath(ROOT, rel...)
    ispath(src) || (@warn "resource missing, not staged" src; return 0)
    dst = joinpath(SHARE, rel...)
    mkpath(dirname(dst))
    cp(src, dst; force = true)
    return treesize(dst)
end

# QMLMakie's own QML module ("Makie", supplying MakieArea) is registered by its `__init__` with
# `QML.add_import_path(joinpath(@__DIR__, "qml"))` — and `@__DIR__` is expanded at PRECOMPILE
# time, so the path baked into the image is this machine's package directory. In a bundle it
# does not exist, Qt cannot find the module, and every .qml that imports it fails with
# `module "Makie" is not installed` before the window appears. A few kB, and the application
# registers it from here at startup.
function stage_qmlmakie()
    pid = Base.PkgId(Base.UUID("08f9cac3-3b11-4f1c-9d88-d0e81c500f64"), "QMLMakie")
    src = Base.locate_package(pid)
    src === nothing && (@warn "QMLMakie not found; the Makie QML module is not staged"; return 0)
    from = joinpath(dirname(src), "qml")
    isdir(from) || (@warn "QMLMakie has no qml directory" from; return 0)
    to = joinpath(SHARE, "qml-modules")
    mkpath(dirname(to))
    cp(from, to; force = true)
    return treesize(to)
end

staged = 0
staged += stage_qmlmakie()
staged += stage("src", "gui", "qml")       # FATAL if missing: loadqml has no fallback
staged += stage("demos", "data")           # the picker's "ROTIR data" place, and gui()'s default dir
staged += stage("demos", "orbits")         # the beta Lyrae and Spica presets orbit_dir() seeds
staged += stage("app", "assets")           # the icon, for the desktop entry and the installer
@info "Resources staged" dir = SHARE MB = round(staged / 1e6, digits = 1)

# ── Makie's font atlas ───────────────────────────────────────────────────────
#
# The trace has already built it, into this machine's scratch space. Copy it beside the
# resources; the application seeds a per-user cache from there on first run, because the
# bundle's own depot is read-only once installed.

try
    Makie = Base.require(Base.PkgId(Base.UUID("ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"), "Makie"))
    cache = Base.invokelatest(Makie.get_cache_path)
    bins  = filter(f -> endswith(f, ".bin"), readdir(cache; join = true))
    if isempty(bins)
        @warn "no font atlas found; the first launch will render one" cache
    else
        dst = joinpath(SHARE, "makie-cache")
        mkpath(dst)
        for b in bins
            cp(b, joinpath(dst, basename(b)); force = true)
        end
        @info "Font atlas staged" files = length(bins) dir = dst
    end
catch err
    @warn "could not stage the font atlas; the first launch will render one" err
end

# ── the console window ───────────────────────────────────────────────────────

if Sys.iswindows()
    exe = joinpath(OUT, "bin", "ROTIRApp.exe")
    if set_windows_gui_subsystem!(exe)
        @info "Windows: linked for the windows subsystem; no console window" exe
        @info "stdout and stderr therefore go nowhere — crash.log is the diagnostic" *
              " (see crash_log_path() in app/src/ROTIRApp.jl)"
    end
end

# ── the launcher ─────────────────────────────────────────────────────────────

if Sys.islinux()
    cp(joinpath(APPSRC, "launcher.sh"), joinpath(OUT, "ROTIR"); force = true)
    chmod(joinpath(OUT, "ROTIR"), 0o755)
    # Menu entry, icon and the .oifits association. Run once after unpacking; it computes its
    # own paths, so the bundle can live anywhere and be moved by re-running it.
    cp(joinpath(APPSRC, "install-desktop.sh"), joinpath(OUT, "install-desktop.sh"); force = true)
    chmod(joinpath(OUT, "install-desktop.sh"), 0o755)
end

total = treesize(OUT)
@info "Done" seconds = round(t, digits = 1) total_MB = round(total / 1e6, digits = 1)
println("""

Start it with:

    $(joinpath(OUT, Sys.islinux() ? "ROTIR" :
                    joinpath("bin", Sys.iswindows() ? "ROTIRApp.exe" : "ROTIRApp"))) [file.oifits]

$(Sys.islinux() ? """
On Linux go through the launcher rather than bin/ROTIRApp: a bundle cannot choose its own
windowing system, because GLFW and Qt are both initialised from the sysimage before any code
in this package runs, and left alone on a Wayland session they land on different ones.
""" : """
There is no launcher on this platform and none is needed: GLFW defaults to the native
windowing system here, so the two halves agree without being told to.
""")""")
