#!/usr/bin/env julia
# The Kurucz ATLAS9 intensity-pack reader (src/kurucz_intensity.jl).
#
# Standalone script in the style of test_parametric_gradient.jl: prints its own table, never
# throws, exposes `nfail[]` for runtests.jl. No AD package and NO DOWNLOAD needed.
#
#     julia --project=. test/test_kurucz_intensity.jl
#
# WHY SYNTHETIC FILES. A real pack is 67 MB and lives on a web server that has already moved
# once in this project's lifetime (TLUSTY's canonical host is gone entirely). So the fixtures
# here are written in-test, byte-for-byte in Kurucz's `format(F9.2,1pe10.3,16I6)`, with KNOWN
# intensities — which also lets the ratio decoding be checked exactly rather than for
# plausibility.
#
# THE TRAP THESE GUARD. The two distributed packs differ by ONE COLUMN and nothing announces
# it: `ip00k2.pck19` keeps Fortran's carriage-control column (header `TEFF   3500.`) while
# `ip00k2new.pck` has had it stripped (header `EFF   3500.`, every data line shifted one left).
# The record is fixed-width and its fields abut at the extremes of their ranges — a 160000.00
# nm wavelength fills F9.2 exactly, a ratio of 100000 fills I6 — so whitespace splitting
# cannot recover them, and hard-coding either layout misreads the other into a float that
# still parses. Both layouts are therefore written and read here.
#
# Verified against the real files when they were present: the two packs agree on I(mu)/I(1) to
# three decimals at the same (Teff, logg), with I(mu=1) differing only in the fourth digit
# (4.025e-5 vs 4.026e-5) from the different opacity tables. That cross-check is what
# established the shift detection is right, and it cannot be reproduced without the downloads.

using ROTIR, Printf, LinearAlgebra

npass = Ref(0); nfail = Ref(0)
function cb(label, ok)
    ok ? (npass[] += 1) : (nfail[] += 1)
    @printf("  %-58s %s\n", label, ok ? "✓" : "✗")
    return ok
end

# A model's intensities: I(mu=1) absolute, and the 16 other angles as ratios of it. Chosen to
# descend monotonically toward the limb like a real atmosphere, and to be exactly
# representable in I6 after the x100000 scaling so the decode can be checked to the digit.
const RATIOS = [1.0, 0.98, 0.95, 0.92, 0.89, 0.86, 0.82, 0.78, 0.75,
                0.71, 0.67, 0.64, 0.60, 0.56, 0.50, 0.43, 0.35]

"""
    write_pack(path, models, λnm; shift = 0)

A synthetic intensity pack. `shift = 0` is `ip00k2.pck19`'s layout; `shift = 1` is
`ip00k2new.pck`'s, which is *literally the same file with the first character of every line
removed* — that is all Fortran's dropped carriage-control column amounts to, and writing it
that way is what makes this fixture mirror the real difference rather than an idea of it.
"""
function write_pack(path, models, λnm; shift::Int = 0)
    drop(l) = shift == 0 ? l : l[2:end]
    open(path, "w") do io
        for (Teff, logg, centi) in models
            println(io, drop(@sprintf("TEFF %7.0f.  GRAVITY %7.5f  SDSC GRID  [+0.0]",
                                      Teff, logg)))
            println(io, drop("  wl(nm)    Inu(ergs/cm**2/s/hz/ster) for 17 mu in 1221 " *
                             "frequency intervals"))
            println(io, drop("            1.000   .900  .800"))
            for λ in λnm
                # `format(F9.2,1pe10.3,16I6)`, exactly.
                s = @sprintf("%9.2f%10.3e", λ, centi)
                for k in 2:17
                    s *= @sprintf("%6d", round(Int, RATIOS[k] * 100_000))
                end
                println(io, drop(s))
            end
        end
    end
    return path
end

const MODELS = [(20000.0, 3.0, 1.0e-4), (20000.0, 3.5, 1.1e-4),
                (25000.0, 3.0, 2.0e-4), (25000.0, 3.5, 2.2e-4)]
const LAMNM  = [1600.0, 1650.0, 1700.0]          # nm -> 1.60, 1.65, 1.70 um

println("\n[1] both column layouts parse, and parse IDENTICALLY")
let d = mktempdir()
    p0 = write_pack(joinpath(d, "keep.pck"), MODELS, LAMNM; shift = 0)
    p1 = write_pack(joinpath(d, "strip.pck"), MODELS, LAMNM; shift = 1)
    m0 = read_kurucz_models(p0; verbose = false)
    m1 = read_kurucz_models(p1; verbose = false)
    cb("carriage-control KEPT (TEFF): 4 models", length(m0) == 4)
    cb("carriage-control STRIPPED (EFF): 4 models", length(m1) == 4)
    # THE decisive check: the one-column shift must not change a single number.
    cb("the two layouts agree exactly on Teff/logg",
       [(m.Teff, m.logg) for m in m0] == [(m.Teff, m.logg) for m in m1])
    cb("the two layouts agree exactly on every intensity",
       all(m0[i].I == m1[i].I for i in eachindex(m0)))
    cb("the two layouts agree exactly on lambda", all(m0[i].λ == m1[i].λ for i in eachindex(m0)))

    println("\n[2] the packed RATIOS are decoded, not read as intensities")
    m = m0[1]
    cb("Teff and logg recovered", m.Teff == 20000.0 && m.logg == 3.0)
    cb("lambda converted nm -> metres", m.λ ≈ [1.60e-6, 1.65e-6, 1.70e-6])
    cb("I(mu=1) is the absolute continuum column", m.I[1, 1] ≈ 1.0e-4)
    # I(mu_k) = centi * round(ratio*1e5)/1e5 — exact, because the fixture's ratios survive I6.
    want = [1.0e-4 * round(Int, RATIOS[k] * 100_000) / 100_000 for k in 1:17]
    cb("all 17 angles decode to centi x ratio", m.I[:, 1] ≈ want)
    # Reading the integers AS intensities would give ~1e5 times too much, uniformly — which
    # cancels out of a flux-normalised visibility and so survives any relative-only test.
    cb("not off by the 1e5 packing factor", maximum(m.I) < 1.0e-3)
    cb("monotone decreasing toward the limb (file order)", issorted(m.I[:, 1], rev = true))

    println("\n[3] read_kurucz_intensity -> a usable RectGrid4")
    g = read_kurucz_intensity(p0; λrange = (1.55e-6, 1.75e-6), verbose = false)
    cb("shape is (nTeff, nlogg, 17, nlambda)", size(g) == (2, 2, 17, 3))
    cb("Teff axis ascending", g.Teff == [20000.0, 25000.0])
    cb("logg axis ascending", g.logg == [3.0, 3.5])
    # RectGrid4 REQUIRES ascending axes; the file's mu descends, so it must be reversed.
    cb("mu axis REVERSED to ascending", issorted(g.μ) && g.μ[1] ≈ 0.010 && g.μ[end] ≈ 1.0)
    cb("mu = 0 is absent (RectGrid4 refuses it)", first(g.μ) > 0)
    # And the reversal must carry the VALUES with it, not just relabel the axis.
    cb("values follow the reversal: disc centre is brightest",
       argmax(g.values[1, 1, :, 1]) == length(g.μ))
    cb("the brightest value is I(mu=1)", g.values[1, 1, end, 1] ≈ 1.0e-4)

    println("\n[4] it works as a provider, and survives a FITS round trip")
    prov = TabulatedProvider(g; name = "synthetic Kurucz")
    I = provider_intensity(prov, [22000.0, 24000.0], [3.2, 3.4], [0.5, 0.95], 1.65e-6)
    cb("provider_intensity returns finite positives", all(isfinite, I) && all(>(0), I))
    cb("owns_mu and needs_logg are both true", owns_mu(prov) && needs_logg(prov))
    f = joinpath(d, "g.fits"); save_intensity_grid(f, g)
    h = load_intensity_grid(f)
    cb("FITS round trip is exact", h.values == g.values && h.μ == g.μ && h.Teff == g.Teff)

    println("\n[5] lambda filtering, and the failure messages")
    gf = read_kurucz_intensity(p0; λrange = (1.62e-6, 1.72e-6), verbose = false)
    cb("lambda filter narrows the grid", size(gf, 4) == 2 && gf.λ ≈ [1.65e-6, 1.70e-6])
    # A single wavelength is refused, by RectGrid4 rather than here — interpolating in λ needs
    # two nodes, and a grid that cannot be interpolated is not a grid.
    cb("filtering down to ONE lambda is refused",
       (try read_kurucz_intensity(p0; λrange = (1.62e-6, 1.68e-6), verbose = false); false
        catch; true; end))
    cb("a lambda range with no overlap is an error",
       (try read_kurucz_models(p0; λrange = (1.0e-9, 2.0e-9), verbose = false); false
        catch; true; end))
    cb("too few nodes after filtering is an error",
       (try read_kurucz_intensity(p0; λrange = (1.55e-6, 1.75e-6),
                                  Teff = (19000.0, 21000.0), verbose = false); false
        catch; true; end))
end

println("\n[6] the rectangle is filled from the nearest gravity at the SAME Teff")
# Kurucz computes each Teff only where a model converges, so the pack is ragged while a
# RectGrid4 must be rectangular. Filling across Teff would mix in another temperature's limb
# darkening, so the fill is restricted to logg — this checks it is.
let d = mktempdir()
    ragged = [(20000.0, 3.0, 1.0e-4), (20000.0, 3.5, 1.1e-4),
              (25000.0, 3.0, 2.0e-4)]            # 25000/3.5 deliberately missing
    p = write_pack(joinpath(d, "ragged.pck"), ragged, LAMNM; shift = 0)
    g = read_kurucz_intensity(p; λrange = (1.55e-6, 1.75e-6), verbose = false)
    cb("a ragged pack still yields a rectangle", size(g) == (2, 2, 17, 3))
    # The missing (25000, 3.5) must be filled from (25000, 3.0) — same Teff — whose continuum
    # is 2.0e-4, NOT from (20000, 3.5) whose continuum is 1.1e-4.
    cb("filled from the same Teff, not the same logg",
       g.values[2, 2, end, 1] ≈ 2.0e-4)
end

@printf("\n=== %d passed, %d failed ===\n", npass[], nfail[])
