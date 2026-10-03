#!/usr/bin/env julia
# Nested samplers must RUN on a multithreaded session.
#
# THE BUG THIS EXISTS FOR. `fit_sphere_ld`, `fit_ellipsoid_ld` and `fit_orbit` each passed
# `min_num_live_points` (and two of them `use_stepsampler`) down to `_fit_nested` while
# SILENTLY DROPPING `threaded`, so every one inherited that function's `threaded = true`.
# Nautilus then parallelises its likelihood calls with `Threads.@threads` while ROTIR's χ²
# threads internally with `Threads.@threads :static` (src/type3_nufft.jl, reached by the
# default `:t3` kernel), and Julia refuses:
#
#     nested task error: `@threads :static` cannot be used concurrently or nested
#
# So `method = :nautilus` could not run AT ALL on any session with more than one thread —
# on all three fitters, including the ellipsoid, which is the one whose own docstring says a
# nested sampler is most needed ("5 parameters with two near-degenerate directions").
#
# Nothing caught it: `test_orbit_fit.jl` exercises `:neldermead` only, and no other test runs a
# nested sampler. The static check below is the cheap guard that would have; the functional
# ones are the proof.

using Test
using ROTIR

@testset "nested samplers under threads" begin

# ── 1. the static guard: every call site forwards `threaded` ────────────────────────────
#
# Costs nothing, needs no sampler, and runs even on one thread. It catches the regression
# class directly — a new fitter that forgets to forward the keyword fails here rather than
# months later on someone's 16-thread machine.
@testset "every _fit_nested call site forwards `threaded`" begin
    srcs = [joinpath(PKGDIR, "src", f) for f in readdir(joinpath(PKGDIR, "src"))
            if endswith(f, ".jl")]
    sites = Tuple{String,Int,String}[]
    for f in srcs
        lines = readlines(f)
        for (i, l) in pairs(lines)
            occursin("_fit_nested(", l) || continue
            startswith(strip(l), "#") && continue
            occursin("function _fit_nested", l) && continue
            occursin("_fit_nested(method, chi2_of, names, lo, hi;", l) && continue  # docstring
            # the call spans a few lines; look at the whole argument list
            blob = join(lines[i:min(i + 4, end)], " ")
            push!(sites, (basename(f), i, blob))
        end
    end
    @test length(sites) >= 3          # sphere, ellipsoid, orbit — the walk found them
    for (f, i, blob) in sites
        ok = occursin("threaded", blob)
        ok || @error "a _fit_nested call site does not forward `threaded`" file = f line = i
        @test ok
    end
end

# ── 2. the functional proof ─────────────────────────────────────────────────────────────
#
# Only meaningful with more than one thread: on a single thread the nesting cannot occur and
# the bug is invisible, which is precisely why it survived. Skipped rather than silently
# passed, so the log says which it was.
LAM1 = joinpath(PKGDIR, "demos", "data", "2011Sep02.lam_And_prepped.oifits")

if Threads.nthreads() < 2
    @info "single-threaded session: the nesting bug cannot manifest, skipping the live runs"
elseif !nautilus_available()
    @info "Nautilus not loaded, skipping the live nested-sampler runs"
else
    d = readoifits(LAM1)[1]
    tess = tessellation_healpix(2)      # 192 tessels: enough to exercise the kernel, cheap

    @testset "fit_sphere_ld with :nautilus" begin
        r = fit_sphere_ld([d], tess; tepochs = [0.0], ldtype = 3,
                          radius0 = 0.6, ld0 = 0.3,
                          radius_bounds = (0.3, 1.2), ld_bounds = (0.0, 2.0),
                          weights = [1.0, 1.0, 0.0], method = :nautilus,
                          min_num_live_points = 100, verbose = false)
        @test isfinite(r.radius) && 0.3 <= r.radius <= 1.2
        @test isfinite(r.ld1)
        @test isfinite(r.logz)
        @test size(r.samples, 1) > 0        # it really sampled, rather than returning a point
    end

    @testset "fit_ellipsoid_ld with :nautilus" begin
        # The case the bug hurt most: its own docstring says the region sampler struggles on
        # these five parameters, so a nested sampler is the point rather than a luxury.
        r = fit_ellipsoid_ld([d], tess; tepochs = [0.0], ldtype = 3,
                             req0 = 0.6, f0 = 0.05, inc0 = 45.0, pa0 = 45.0, ld0 = 0.3,
                             seed_from_sphere = false,
                             req_bounds = (0.3, 1.2), f_bounds = (0.0, 0.3),
                             weights = [1.0, 1.0, 0.0], method = :nautilus,
                             min_num_live_points = 100, verbose = false)
        @test isfinite(r.req) && isfinite(r.flattening) && isfinite(r.ld1)
        @test isfinite(r.logz)
        @test size(r.samples, 1) > 0
    end
end
end
