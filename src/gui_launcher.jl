# `rotirgui()` — the GUI from a plain `using ROTIR` session.
#
# WHY THIS IS A FUNCTION AND NOT A SCRIPT. Opening the window needs four things done in one
# specific order, and only the first is obvious:
#
#   1. the Mesa/GLFW/Qt hints, which are read when the first OpenGL context is created and are
#      therefore unreachable from anywhere inside the GUI extension — loading GLMakie is what
#      creates that extension, so by then it is already too late;
#   2. GLMakie, whose platform choice (Wayland vs X11) has to be made before Qt is told to
#      match it;
#   3. QMLMakie and QML;
#   4. the optional sampler stack, which decides which fit methods the Model panel lists.
#
# `bin/rotirgui.jl` did all of that as a script, which meant it only worked from the pinned
# `bin` environment. Everything here is the same sequence, in a function, so a user with
# ROTIR + GLMakie + QMLMakie + QML in their own project gets the window with one call — and
# `bin/rotirgui.jl` now calls this too, so the two paths cannot drift.

# GLFW_jll is NOT a dependency of ROTIR — it arrives underneath GLMakie — so it cannot be
# named in a `using`. Loading it by UUID works regardless of whose dependency it is, which is
# what lets `prefer_native_wayland!` (in OITOOLSGLFWExt) exist before GLMakie is loaded.
const _GLFW_JLL_UUID = Base.UUID("0656b61e-2033-5cc2-a64a-77c0f6c09b89")

# Pkg is a STDLIB but not a dependency of ROTIR, and it should not become one: declaring it
# would load it on every `using ROTIR`, and it is only ever wanted on the one path where
# something is missing. By UUID, on demand, like GLFW_jll above.
const _PKG_UUID = Base.UUID("44cfe95a-1eb2-52ea-b672-e2afdf69b78f")

"The toolkit packages the window needs, in the order they must be loaded."
const GUI_PACKAGES = [:GLMakie, :QMLMakie, :QML]

# The optional stack: each entry adds ONE fit method to the Model panel, and the panel reads
# that list once when it first refreshes, so anything not loaded before the window opens is not
# offered for the rest of the session. That is why they are proposed up front rather than
# installed when someone reaches for the method.
#
# EVERY ONE IS PURE JULIA. `:ultranest` is deliberately absent — not listed, not gated — and
# nothing here can pull PythonCall in. Beyond the sampler, matplotlib's backend probe maps a
# second Qt into the process and QML.jl needs its own Qt to be the only one, so Python in a
# session with a window is not merely slow. `:nautilus` covers the same ground (importance
# nested sampling, an evidence, uncertainties) and agrees with UltraNest to 1.4e-6 mas on a
# λ And sphere fit, so nothing is given up. UltraNest remains available to a SCRIPT:
# `using PythonCall` loads ROTIRUltraNestExt and `method = :ultranest` works as it always did.
const OPTIONAL_PACKAGES = [
    (pkgs = [:Zygote],
     gives = "gradient fit — analytic gradients"),
    (pkgs = [:AdvancedHMC, :LogDensityProblems],
     gives = ":nuts — a posterior (wants Zygote too)"),
    (pkgs = [:Nautilus],
     gives = ":nautilus — a posterior and an evidence"),
    (pkgs = [:Pigeons, :Distributions, :ADTypes],
     gives = ":pigeons — multimodal posteriors; ~2.8 s per start"),
]

"""
    _source_url(name) -> String or nothing

The git URL ROTIR's own `[sources]` gives for `name`, if it has one.

NOT every package here is registered: Nautilus lives at a URL, so `Pkg.add("Nautilus")` fails
with "not found in a registry" — and `[sources]` does NOT propagate from a dependency, so a
user's own project cannot learn the URL from ROTIR's Project.toml by itself. Reading it back
out at install time keeps one copy of the address: change `[sources]` and this follows.
"""
function _source_url(name::Symbol)
    proj = joinpath(pkgdir(ROTIR), "Project.toml")
    isfile(proj) || return nothing
    entry = get(get(TOML.parsefile(proj), "sources", Dict{String,Any}()), String(name), nothing)
    entry isa AbstractDict || return nothing
    url = get(entry, "url", nothing)
    return url isa AbstractString ? url : nothing
end

"""
    _pkg_add(Pkg, names) -> nothing

`Pkg.add` the named packages, by URL for the unregistered ones.

Called through `invokelatest` because `Pkg` was loaded moments ago: `Pkg.PackageSpec` is a
constructor in a module newer than the frame that wants it.
"""
function _pkg_add(Pkg, names::Vector{Symbol})
    specs = map(names) do n
        url = _source_url(n)
        url === nothing ? Pkg.PackageSpec(name = String(n)) :
                          Pkg.PackageSpec(name = String(n), url = url)
    end
    Pkg.add(specs)
    return nothing
end

"The optional groups that are not fully installed, as indices into `OPTIONAL_PACKAGES`."
_missing_optional() =
    [i for (i, o) in enumerate(OPTIONAL_PACKAGES) if !isempty(_not_installed(o.pkgs))]

"""
    _not_installed(names) -> Vector{Symbol}

Which of `names` the active project cannot resolve. Asked BEFORE loading anything, so that a
missing set can be installed in one `Pkg.add` rather than discovered one `using` at a time.
"""
_not_installed(names) =
    filter(n -> Base.identify_package(Main, String(n)) === nothing, names)

# Relative to `Main`, because `Main` is where the `using` happens. And note that a package can
# be resolvable WITHOUT being a direct dependency of the active project: ROTIR's weakdeps are
# recorded in the manifest so its extensions can load, so in a project that only has ROTIR,
# `using QML` already works and only GLMakie is genuinely missing. Asking the loader, rather
# than reading the project file, is what gets that right.

"\"X\", \"X and Y\", \"X, Y and Z\" — a list a sentence can contain."
function _and_list(xs)
    n = length(xs)
    n == 0 && return ""
    n == 1 && return string(xs[1])
    return join(xs[1:(end - 1)], ", ") * " and " * string(xs[end])
end

"""
    _ask_what_to_install(missing_required, optional_idx, proj) -> Vector{Int} or nothing

The interactive offer. Returns which optional groups to add, or `nothing` if cancelled.

ENTER is the important key: it installs only what the window cannot open without and goes
straight to the GUI. Anything else has to be typed, so the default answer costs nothing and
the fast path stays fast.
"""
function _ask_what_to_install(missing_required::Vector{Symbol}, optional_idx::Vector{Int},
                              proj::AbstractString)
    println()
    if isempty(missing_required)
        printstyled("rotirgui: optional packages are not installed.\n"; bold = true)
    else
        printstyled("rotirgui needs $(_and_list(missing_required)).\n"; bold = true)
        println("  These will be added to $(proj), and precompiled — a few minutes, once.")
    end
    if !isempty(optional_idx)
        println()
        println("  Optional. Each adds one fit method to the Model panel, and the panel reads")
        println("  that list when it opens — so a method not installed now is not offered for")
        println("  the rest of the session. All pure Julia; nothing here loads Python.")
        println()
        for (k, i) in enumerate(optional_idx)
            o = OPTIONAL_PACKAGES[i]
            @printf("   %d  %-34s %s\n", k, join(o.pkgs, " + "), o.gives)
        end
    end
    println()
    if isempty(missing_required)
        println("  ENTER  open the window now, without them")
    else
        println("  ENTER  install $(_and_list(missing_required)) only, then open the window")
    end
    isempty(optional_idx) || println("  a      install the optional ones as well")
    isempty(optional_idx) || println("  1 3    install only those, by number")
    println("  n      cancel")
    print("> ")
    flush(stdout)
    reply = try
        strip(lowercase(readline()))
    catch
        "n"
    end
    isempty(reply) && return Int[]
    reply == "n" && return nothing
    reply == "a" && return copy(optional_idx)
    picked = Int[]
    for w in split(reply)
        k = tryparse(Int, w)
        if k === nothing || !(1 <= k <= length(optional_idx))
            println("  did not understand \"$(w)\" — installing nothing optional")
            return Int[]
        end
        push!(picked, optional_idx[k])
    end
    return picked
end

"""
    _provision(install, samplers) -> nothing

Make sure the window can open, and offer the optional stack while asking.

`install` is `:ask` (the default), `true` (everything, no prompt), `:required` (the toolkit
only, no prompt) or `false` (explain and stop).

It only PROMPTS when something required is missing, which is the moment the user is already
blocked. When the window can open, a missing optional package is one `@info` line rather than
a question on every launch — nagging someone who has deliberately not installed Pigeons is a
good way to make them stop reading the output.
"""
function _provision(install, samplers::Bool)
    install in (true, false, :ask, :required) ||
        throw(ArgumentError("install must be true, false, :ask or :required; " *
                            "got $(repr(install))"))
    req = _not_installed(GUI_PACKAGES)
    opt = samplers ? _missing_optional() : Int[]
    isempty(req) && isempty(opt) && return nothing

    proj = something(Base.active_project(), "the active project")
    listed = join(("\"$(n)\"" for n in req), ", ")
    explain = """
        the GUI needs $(_and_list(req)), which $(length(req) == 1 ? "is" : "are") not \
        installed in this environment.

        $(length(req) == 1 ? "It is a WEAK dependency" : "They are WEAK dependencies") of \
        ROTIR — deliberately, so that `using ROTIR` in a script costs no Makie and no Qt — \
        so they have to be in your own project:

            using Pkg; Pkg.add([$(listed)])

        or let the launcher do it:

            rotirgui(; install = true)

        From a ROTIR clone there is also a pinned environment that already has them:

            julia --project=bin bin/rotirgui.jl
        """

    chosen = Int[]
    if install === false
        isempty(req) || error(explain)
    elseif install === :required
        # `chosen` stays empty; the required set is installed below.
    elseif install === true
        chosen = opt
    else                                # :ask
        if isempty(req)
            # Nothing is blocked, so say what is not offered and get out of the way.
            missing_names = unique(reduce(vcat, (OPTIONAL_PACKAGES[i].pkgs for i in opt)))
            @info """
                  not installed, so those fit methods will not be offered:
                  $(join(("  " * join(OPTIONAL_PACKAGES[i].pkgs, " + ") * " → " *
                          first(split(OPTIONAL_PACKAGES[i].gives, " — ")) for i in opt), "\n"))
                  `rotirgui(install = true)` adds them.
                  """
            return nothing
        end
        isinteractive() || error(explain)
        picked = _ask_what_to_install(req, opt, proj)
        picked === nothing && error("cancelled. " * explain)
        chosen = picked
    end

    add = copy(req)
    for i in chosen; append!(add, _not_installed(OPTIONAL_PACKAGES[i].pkgs)); end
    unique!(add)
    isempty(add) && return nothing

    Base.require(Base.PkgId(_PKG_UUID, "Pkg"))
    Pkg = Base.loaded_modules[Base.PkgId(_PKG_UUID, "Pkg")]
    @info "installing $(_and_list(add)) into $(proj)"
    # `invokelatest`: Pkg was loaded on the line above, so its methods are newer than this
    # frame — the same rule that governs every other load in this file.
    Base.invokelatest(_pkg_add, Pkg, add)
    return nothing
end

"""
    _load_or_explain(names) -> nothing

`using` the named packages in `Main`, or raise an error naming them together.

Reached only after `_provision`, so a failure here is not a missing package but
a broken one — a precompilation error, or a version that cannot resolve.
"""
function _load_or_explain(names::Vector{Symbol})
    try
        Core.eval(Main, Expr(:using, (Expr(:., n) for n in names)...))
    catch err
        error("""
              the GUI needs $(_and_list(names)), which \
              $(length(names) == 1 ? "is" : "are") installed but failed to load:

                  $(sprint(showerror, err))
              """)
    end
    return nothing
end

"Load a package by UUID, returning whether it worked. Used for a package nobody declares."
function _try_require(uuid::Base.UUID, name::AbstractString)
    try
        Base.require(Base.PkgId(uuid, String(name)))
        return true
    catch
        return false
    end
end

"""
    rotirgui(files...; samplers = true, session = nothing)

Open the ROTIR GUI, loading everything it needs.

```julia
using ROTIR
rotirgui()                                  # an empty session
rotirgui("mydata.oifits")                   # with a dataset loaded
rotirgui("ep1.oifits", "ep2.oifits")        # several epochs
```

GLMakie, QMLMakie and QML must be installed in the active project — they are weak dependencies
of ROTIR, so that `using ROTIR` in a script costs no Makie and no Qt — and this says how to add
them if they are missing. Everything else it does itself, in the order the graphics stack
requires.

# Keywords
- `samplers = true` — also load Zygote, AdvancedHMC/LogDensityProblems, Nautilus and Pigeons,
  each optional and skipped with a note if absent. They decide which fit methods the Model
  panel lists, and the panel reads that list ONCE when it first refreshes, so a method not
  loaded here is not offered for the rest of the session. `false` starts faster —
  MEASURED, Pigeons alone is worth about 2.8 s, because loading it invalidates the plot
  code OITOOLS precompiles for its live canvas — at the cost of `:pigeons`, `:nautilus`,
  `:nuts` and the gradient fit not appearing.
- `install = :ask` — what to do about packages that are not installed. `:ask` (interactive
  only) lists what is missing, offers the optional stack alongside it, and treats ENTER as
  "install only what the window needs, then open it". `true` installs everything without
  asking, `:required` installs only the toolkit, `false` explains and stops. When the window
  can already open, `:ask` does not prompt at all — a missing optional package becomes one
  `@info` line, not a question on every launch.
- `session` — an existing `Session` to reopen instead of a fresh one.

Anything else is passed to `gui` — `autoquit_ms` and `qmlfile`, which is how an
automated run drives the window without a person to close it.

Returns when the window closes.
"""
function rotirgui(files::AbstractString...; samplers::Bool = true, session = nothing,
                  install = :ask, kwargs...)
    # FIRST, before the graphics hints: if something has to be installed, the question — and
    # the several minutes of precompilation that follow a yes — should come before any of the
    # environment fiddling below, not be buried underneath it.
    _provision(install, samplers)
    # A CondaPkg resolve that was interrupted leaves its lock behind, and every later start
    # then blocks on it FOREVER: the process stays alive, no window appears, and the only clue
    # is one CondaPkg info line, indistinguishable from a graphics hang.
    let lock = joinpath(pkgdir(ROTIR), ".CondaPkg", "lock")
        if isfile(lock)
            @warn """
                A CondaPkg lock file exists and startup will block on it until it is released.

                    $lock   (age $(round(Int, time() - mtime(lock)))s)

                If no other Julia process is resolving a Conda environment, it is stale:

                    rm -f $lock
                """
        end
    end

    # The hints, before anything can create a GL context. If GLMakie is ALREADY loaded they
    # are too late to have any effect, and saying so is better than setting them silently and
    # leaving the user to wonder why the window still landed on XWayland.
    late = haskey(Base.loaded_modules,
                  Base.PkgId(Base.UUID("e9467ef8-e4e7-5192-8a1a-b1aee30e663a"), "GLMakie"))
    if late
        @info("GLMakie was already loaded, so the graphics hints cannot be applied — the " *
              "window will work, but on whatever platform GLMakie already chose. Call " *
              "`rotirgui()` from a fresh session to let it pick.")
    else
        OITOOLS.configure_graphics!()
        # GLMakie's platform is decided FIRST, then Qt is made to match: one on Wayland and
        # one on X11 means two EGL display connections in one process. Only GLMakie's choice
        # can go wrong on its own — GLFW.jl hard-codes X11 on Linux, so under a Wayland
        # compositor `using GLMakie` silently lands on XWayland.
        wl = if _try_require(_GLFW_JLL_UUID, "GLFW_jll")
            # `invokelatest` for the same reason as the `gui` call at the end, one step
            # earlier: `prefer_native_wayland!` is a stub in OITOOLS whose method arrives with
            # OITOOLSGLFWExt, which the `_try_require` on the line above has just loaded — into
            # a world newer than this frame's. Without it: "no method matching
            # prefer_native_wayland!(); the applicable method may be too new".
            Base.invokelatest(OITOOLS.prefer_native_wayland!)
        else
            (; applied = false)
        end
        OITOOLS.configure_qt_platform!(; match_x11 = !wl.applied)
    end

    # Separately, and in this order: GLMakie's platform choice has to be made — and Qt made
    # to match it — before QML is loaded.
    # The Qt Quick Controls style, before QML is imported. See `gui` in src/gui/window.jl for
    # why macOS needs asking: its native style makes the tab bar a third too tall and runs a
    # ComboBox's dropdown arrow over its own text.
    Sys.isapple() && get!(ENV, "QT_QUICK_CONTROLS_STYLE", "Fusion")

    _load_or_explain([:GLMakie])
    _load_or_explain([:QMLMakie, :QML])

    if samplers
        for o in OPTIONAL_PACKAGES
            # SILENTLY skipped when not installed: `_provision` has already listed exactly
            # these, and saying it twice trains the reader to skip the whole block. A failure
            # to load something that IS installed is different — that is a broken package,
            # and it gets a warning.
            isempty(_not_installed(o.pkgs)) || continue
            try
                Core.eval(Main, Expr(:using, (Expr(:., p) for p in o.pkgs)...))
            catch err
                @warn "$(join(o.pkgs, ", ")) installed but failed to load; " *
                      "$(first(split(o.gives, " — "))) will not be offered" exception = err
            end
        end
    end

    # `invokelatest`, and it is REQUIRED. Every method above — `gui` itself, and everything in
    # ROTIRGUIExt it reaches — was added to the running session by the `using` calls in this
    # very frame, and methods added after a frame's world age are invisible to it. Without
    # this the call lands on the `function gui end` stub and fails with "no method matching",
    # having just successfully loaded the extension that defines it.
    return Base.invokelatest(_open_gui, collect(String, files), session; kwargs...)
end

"The half of `rotirgui` that runs in the new world: everything here needs the extension."
function _open_gui(files::Vector{String}, session; kwargs...)
    ext = Base.get_extension(ROTIR, :ROTIRGUIExt)
    ext === nothing && error("GLMakie, QMLMakie and QML all loaded but ROTIRGUIExt did not; " *
                             "check for a precompilation error above")
    s = session === nothing ? ext.Session() : session
    for f in files
        isfile(f) || (@warn "not a file, skipping" path = f; continue)
        ext.load_dataset!(s, f)
    end
    return gui(s; kwargs...)
end
