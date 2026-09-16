#!/usr/bin/env julia
# ROTIR test suite —  run with:  julia --project=. -e 'using Pkg; Pkg.test()'
#
# The individual files are standalone scripts (they print their own tables and can be run
# directly), so each is included into its own throw-away module to keep top-level names
# from colliding. `test_parametric_gradient.jl` tracks its own pass/fail counters, which
# are asserted here — it never throws on a numerical mismatch by itself.
#
# Figure-generating tests (test_spot_euclidean.jl, ~5 PNGs into test/figures/) are opt-in:
#     ROTIR_TEST_FIGURES=1 julia --project=. -e 'using Pkg; Pkg.test()'

# Headless-safe matplotlib: must be set before ROTIR (hence PyPlot) is loaded anywhere.
get!(ENV, "MPLBACKEND", "Agg")

using Test

const TESTDIR = @__DIR__
const PKGROOT = dirname(TESTDIR)

"""
    run_script(file) -> Module

Include `file` in a fresh module, with the working directory at the package root (the
scripts use paths like `./demos/data/...`). Returns the module so counters can be read.
"""
function run_script(file)
    m = Module(Symbol(:Script_, replace(basename(file), r"\.jl$" => "")))
    cd(PKGROOT) do
        Base.include(m, joinpath(TESTDIR, file))
    end
    return m
end

@testset "ROTIR" begin
    @testset "parametric gradient" begin
        m = run_script("test_parametric_gradient.jl")
        # The script skips its Zygote section with a @warn if the AD load fails; require
        # that the section actually ran, otherwise a silent skip would look like success.
        @test m.zygote_ran[]
        @test m.nfail[] == 0
    end

    @testset "bootstrap" begin
        m = run_script("test_bootstrap.jl")
        @test m.zygote_ran[]     # the extension section must run, not silently skip
    end

    # Binary frame + mutual irradiation. Self-contained @test assertions (analytic limits
    # and conservation laws, no recorded reference values), so it is included directly
    # rather than through run_script.
    include(joinpath(TESTDIR, "test_reflection.jl"))

    # Every plotting entry point with every decoration, asserting on artist structure
    # rather than pixels. Headless (Agg) and a few seconds, so it runs by default —
    # structural checks are what catch a corrupted Julia->Python array conversion.
    include(joinpath(TESTDIR, "test_plotting.jl"))

    # Generic orbit fitter: component library, spec assembly, and a synthetic round trip.
    include(joinpath(TESTDIR, "test_orbit_fit.jl"))

    # RADFLAT / RADIALVAR / orthoLD: hand-written gradients against FiniteDifferences, plus
    # the ANOVA identity that makes RADFLAT and RADIALVAR complementary.
    include(joinpath(TESTDIR, "test_radial_regularizers.jl"))

    # The reconstruction's callable surface: the mid-run callback, per-epoch weights, and the
    # joint fitter's refusal to guess a θ layout it does not have.
    include(joinpath(TESTDIR, "test_reconstruct_hooks.jl"))

    # Every name the extensions import from the core package, without loading any of them —
    # a phantom import is only a WARNING at load time, and one had been there since the file
    # was written.
    include(joinpath(TESTDIR, "test_extensions.jl"))

    # The @turbo forward kernel against the scalar reference it rewrites — both precisions,
    # three mesh levels, three surface types, three datasets, and the degenerate cases.
    include(joinpath(TESTDIR, "test_fused_polyft.jl"))

    include(joinpath(TESTDIR, "test_type3_nufft.jl"))

    # The binary forward model: the matrix-free route against the dense one it replaced.
    # AFTER the file above, which loads LoopVectorization — otherwise the `:turbo` case here
    # skips itself and the backend it is meant to pin goes untested.
    include(joinpath(TESTDIR, "test_binary_cvis.jl"))

    # The binary IMAGING gradient: finite differences on both components, and the reduction
    # to `spheroid_chi2_fg` when the secondary is dark.
    include(joinpath(TESTDIR, "test_binary_imaging.jl"))

    # The geometry and map FITS files: the mesh back bit-for-bit, and the ASCII aliasing that
    # lets a Roche model's Unicode orbital elements into a FITS string column at all.
    include(joinpath(TESTDIR, "test_geometry_io.jl"))

    # The SPHERE's log-posterior, which is what lets NUTS run on a single star: its gradient
    # against FiniteDifferences, and — the part that matters — that every parameter in θ
    # actually moves the posterior. Needs Zygote, so it is guarded like the other AD tests.
    include(joinpath(TESTDIR, "test_sphere_logpi.jl"))

    # The ELLIPSOID's von Zeipel map differentiated: the piece a consistent gradient fit for
    # surface type 1 was missing, and the proof that its map does not depend on orientation —
    # which is what makes a shape-only gradient exact when only the angles are free.
    include(joinpath(TESTDIR, "test_ellipsoid_map_derivs.jl"))

    # The ELLIPSOID's log-posterior: the two analytic rrules composed, and the measurement
    # behind refusing to free `tpole` under a linear intensity law.
    include(joinpath(TESTDIR, "test_ellipsoid_logpi.jl"))

    # The two GRAVITY-DARKENING laws and the registry that selects between them: ELR against
    # the closed forms in the paper, the equatorial symmetry a finite-difference check cannot
    # see, and both laws through `build_parametric_logπ`.
    include(joinpath(TESTDIR, "test_gravity_darkening.jl"))

    if get(ENV, "ROTIR_TEST_FIGURES", "0") == "1"
        @testset "spot placement (figures)" begin
            run_script("test_spot_euclidean.jl")   # contains its own @test assertions
        end
    else
        @info "Skipping figure tests; set ROTIR_TEST_FIGURES=1 to enable."
    end
end
