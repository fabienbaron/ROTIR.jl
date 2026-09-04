# The in-window file picker.
#
# NOT `QtQuick.Dialogs.FileDialog`. On some systems that leaves its own window mapped after it
# is dismissed, and — worse for this GUI — a native dialog gets its own GL surface, which is
# one of the two known ways to permanently break the shared context (see the note at the top of
# src/gui/livecanvas.jl). Drawing the picker inside our own window avoids both, at the cost of
# Julia having to supply everything the toolkit would: the listing, the parent directory, the
# shortcut places, and path joining that behaves the same on every platform.
#
# Everything crosses as text, like the rest of the shell.

# Extensions per PURPOSE. The picker is one widget opened for three different things, and a
# listing that shows every one of them makes each opening harder: hunting an orbit TOML among
# forty OIFITS files, or the reverse. Case-insensitive.
const PICKER_EXTENSIONS = (".oifits", ".fits", ".fit", ".oifit")
const PICKER_EXTENSIONS_ORBIT = (".toml",)
const PICKER_EXTENSIONS_MAP   = (".fits", ".fit")

"Extensions the picker offers for a purpose."
picker_extensions(purpose::AbstractString) =
    purpose == "orbit" ? PICKER_EXTENSIONS_ORBIT :
    purpose == "map"   ? PICKER_EXTENSIONS_MAP   : PICKER_EXTENSIONS

"""
    picker_places() -> String

Shortcut destinations, `label\tpath` per line.

The package's `demos/data` is listed first and by name, because it is where every dataset that
ships with ROTIR lives and it is otherwise buried inside a depot path nobody would navigate to.
The orbit folder is here for the same reason — it is under the platform's per-user config
directory, which is not somewhere anyone would think to browse to by hand.
"""
function picker_places()
    rows = String[]
    d = joinpath(pkgdir(ROTIR), "demos", "data")
    isdir(d) && push!(rows, "ROTIR data\t$(d)")
    # Whether or not this picker is being opened FOR an orbit: a shortcut that appears and
    # disappears depending on which button opened the dialog is harder to learn than one that
    # is always in the same place. `orbit_dir` creates the folder, which is what makes the
    # shortcut land somewhere rather than on a path that does not exist yet.
    o = try
        orbit_dir()
    catch
        ""
    end
    isempty(o) || !isdir(o) || push!(rows, "Orbits\t$(o)")
    push!(rows, "Home\t$(homedir())")
    push!(rows, "Working dir\t$(pwd())")
    forced = get(ENV, "ROTIRGUI_DATA_DIR", "")
    isempty(forced) || (isdir(forced) && push!(rows, "ROTIRGUI_DATA_DIR\t$(abspath(forced))"))
    return join(unique(rows), "\n")
end

"""
    picker_list(dir, show_all) -> String

One entry per line: `kind\tname\tsize`, directories first and each group sorted by name.

`kind` is `dir` or `file`. `size` is bytes for a file and empty for a directory. Hidden entries
are skipped unless `show_all` is "1"; an unreadable entry is skipped rather than raising, since
one file with no permissions must not empty the whole listing.

`purpose` picks the extension filter: `data` for OIFITS, `orbit` for a saved orbit TOML, `map`
for a surface-map FITS.
"""
function picker_list(dir, show_all = "0", purpose = "data")
    d = String(dir)
    isdir(d) || return ""
    all = String(show_all) == "1"
    exts = picker_extensions(String(purpose))
    dirs = String[]; files = String[]
    for name in try; sort(readdir(d)); catch; String[]; end
        (all || !startswith(name, ".")) || continue
        p = joinpath(d, name)
        try
            if isdir(p)
                push!(dirs, "dir\t$(name)\t")
            elseif all || _picker_wanted(name, exts)
                push!(files, "file\t$(name)\t$(filesize(p))")
            end
        catch
            continue
        end
    end
    return join(vcat(dirs, files), "\n")
end

_picker_wanted(name, exts = PICKER_EXTENSIONS) =
    any(endswith(lowercase(name), e) for e in exts)

"""
    picker_parent(dir) -> String

The parent directory, or `dir` itself at the filesystem root — never an empty string, which
QML would render as a listing of the working directory and make "up" look like a jump.
"""
function picker_parent(dir)
    d = abspath(String(dir))
    p = dirname(d)
    return isempty(p) || p == d ? d : p
end

"Join a directory and a name the way this platform does."
picker_join(dir, name) = joinpath(String(dir), String(name))

"""
    picker_start(folder) -> String

Turn the `initialFolder` QML was handed (a `file://` URL) into a directory that exists,
falling back to the working directory. QML calls this once, on open.
"""
function picker_start(folder)
    f = String(folder)
    startswith(f, "file://") && (f = f[8:end])
    isdir(f) && return abspath(f)
    isfile(f) && return dirname(abspath(f))
    return pwd()
end
