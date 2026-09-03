# The extensions' contract with the core package, checked without loading any of them.
#
# An extension reaches into ROTIR for two kinds of name: the STUBS it adds methods to
# (`import ROTIR: plot2d`), and the internals it calls (`using ROTIR: graticule_segments`).
# Both are written as plain import lists in ext/*.jl and neither is checked by anything until
# the extension is loaded — and even then the failure is a warning, not an error:
#
#     WARNING: Imported binding ROTIR.list_free_params was undeclared at import time
#              during import to ROTIRUltraNestExt
#
# which is what `ROTIRUltraNestExt` did for as long as it existed. `list_free_params` is a
# keyword argument and a field of the returned NamedTuple; there is no such function, and the
# import asked for a binding that does not exist. It happened to be harmless, but the same
# warning is what a RENAMED internal would produce — and that one breaks the extension on a
# machine that has the trigger package installed, which CI here does not for all of them.
#
# So: parse every import list and require that ROTIR actually defines what it promises. This
# needs no trigger package and costs milliseconds.

using Test
using ROTIR

@testset "extension imports resolve" begin
    extdir = joinpath(pkgdir(ROTIR), "ext")
    @test isdir(extdir)
    files = sort(filter(f -> endswith(f, ".jl"), readdir(extdir)))
    @test !isempty(files)

    "Every name in this file's `using ROTIR:` / `import ROTIR:` lists."
    function rotir_imports(path)
        names = Symbol[]
        lines = readlines(path)
        i = 1
        while i <= length(lines)
            m = match(r"^\s*(?:using|import)\s+ROTIR:\s*(.*)$", lines[i])
            if m !== nothing
                # An import list continues while the line ends in a comma — which is how all
                # of these are wrapped, since they are longer than the line limit.
                body = strip(m.captures[1])
                while endswith(body, ",") && i < length(lines)
                    i += 1
                    body *= " " * strip(lines[i])
                end
                for n in split(body, ",")
                    n = strip(first(split(n, "#")))            # drop a trailing comment
                    isempty(n) && continue
                    push!(names, Symbol(n))
                end
            end
            i += 1
        end
        return names
    end

    total = 0
    for f in files
        ns = rotir_imports(joinpath(extdir, f))
        total += length(ns)
        @testset "$(f)" begin
            # Each file must import SOMETHING from ROTIR — an extension that reaches for
            # nothing is either dead or accessing the parent unqualified, which resolves to
            # whatever `using ROTIR` exported and breaks the moment an export is dropped.
            @test !isempty(ns)
            missing_names = filter(n -> !isdefined(ROTIR, n), ns)
            @test missing_names == Symbol[]
        end
    end
    # A wrong regex that matched nothing would pass every test above, so pin the count too.
    @test total > 60
end
