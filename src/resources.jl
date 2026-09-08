# Where the shipped resources live.
#
# The QML, the demo OIFITS and the shipped orbit presets are files on disk that the package
# reads at RUN time, not code a sysimage can carry. Everything else in this package is
# compiled; these are not.
#
# This file exists because `pkgdir` stops being an answer the moment the package is inside an
# application bundle. `PackageCompiler.create_app` bundles the sysimage, the artifacts and
# Julia's libraries, and writes only a stub `Project.toml` — no package source is copied. So
# `pkgdir(ROTIR)` then names a directory on the machine that BUILT the bundle, which does not
# exist on the machine running it, and every `joinpath(pkgdir(ROTIR), ...)` silently misses.
#
# What each missing resource costs, read off the call sites rather than guessed:
#
#   src/gui/qml   gone -> FATAL. `gui()` passes it to `loadqml`, which has no fallback.
#   demos/data    gone -> silent. The picker loses its "Demo data" place (filepicker.jl) and
#                         `gui()`'s default open directory falls back (window.jl).
#   demos/orbits  gone -> silent. `orbit_dir()` creates the user folder but seeds no presets,
#                         so beta Lyrae and Spica are simply absent from the Orbit tab.
#
# The fatal one masks the silent two, which is why all four call sites move together rather
# than one at a time in response to whatever the bundle complains about first.
#
# A BUNDLE MIRRORS THE REPOSITORY'S OWN RELATIVE PATHS — `src/gui/qml`, `demos/data`,
# `demos/orbits` under the resource root. That makes `resource` a root substitution and
# nothing else. A mapping from repository layout to bundle layout would be a second list to
# keep in step with the file picker's places, and it would drift the first time a place was
# added.

"Environment variable that names the resource root outright, overriding every other candidate."
const RESOURCE_DIR_VAR = "ROTIR_RESOURCE_DIR"

"""
    resource_dir() -> String or nothing

The directory holding the shipped resources, or `nothing` if there is none.

Three candidates, in this order, first existing one wins:

 1. `\$ROTIR_RESOURCE_DIR`. This is what lets a relocation test drive the whole mechanism
    without building an application.
 2. `<Sys.BINDIR>/../share/rotir`, where an application bundle keeps them. `create_app` copies
    the executable into `<app>/bin`, so this anchor holds wherever the bundle is unpacked,
    needs nothing writable, and costs no precompile-time capture.
 3. `pkgdir(ROTIR)`, which is a development checkout.

The bundle is tried BEFORE the checkout so that a bundle cannot be shadowed. In practice the
two are mutually exclusive: in a bundle the checkout path does not exist, and in a checkout
Julia's own `share/` has no `rotir` directory.
"""
function resource_dir()
    for r in _resource_roots()
        isdir(r) && return r
    end
    return nothing
end

function _resource_roots()
    roots = String[]
    forced = get(ENV, RESOURCE_DIR_VAR, "")
    isempty(forced) || push!(roots, abspath(expanduser(forced)))
    push!(roots, normpath(joinpath(Sys.BINDIR, "..", "share", "rotir")))
    root = pkgdir(@__MODULE__)
    root === nothing || push!(roots, root)
    return roots
end

"""
    resource(parts...) -> String or nothing

One shipped file or directory, named by its path relative to the repository root, or `nothing`
when no root has it.

    resource("src", "gui", "qml", "Main.qml")
    resource("demos", "data")

`nothing` rather than a non-existent path, because every caller but one treats a missing
resource as "this feature has nothing to offer" and carries on — a picker place is simply not
listed, the shipped orbits are not seeded. The exception is `Main.qml`, whose absence has no
sensible fallback and which therefore reports the path it looked for.
"""
function resource(parts::AbstractString...)
    sub = joinpath(parts...)
    for r in _resource_roots()
        p = joinpath(r, sub)
        ispath(p) && return p
    end
    return nothing
end
