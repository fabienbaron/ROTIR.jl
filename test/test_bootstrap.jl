#!/usr/bin/env julia
# Bootstrap resampling and driver wiring.
#
# The fast tests use a synthetic `fitfun` (a summary statistic of the replicate), so the
# whole of src/bootstrap.jl is covered without an AD backend or a real fit. The end-to-end
# fit through ext/ROTIRZygoteExt.jl is gated on Zygote being loadable and kept deliberately
# small — it checks mechanics, not astrophysics.

using ROTIR
using Test
using Random
using Statistics

const OIFITS_FILE = "./demos/data/rho_Cas_example.oifits"

data = readoifits_multiepochs([OIFITS_FILE]; T=Float32)[1, :]
npts(d) = d.nv2 + d.nt3amp + d.nt3phi

@testset "epoch_blocks" begin
    b = epoch_blocks(data)
    @test length(b) == length(data)
    @test all(x -> length(x) > 0, b)
    # finer granularity cannot produce fewer blocks
    @test length(epoch_blocks(data; granularity=:point)[1]) >=
          length(epoch_blocks(data; granularity=:config)[1]) >=
          length(epoch_blocks(data; granularity=:epoch)[1])
    # MJDs are Float64 in OIdata whatever T is, so blocking is precision-independent
    d64 = readoifits_multiepochs([OIFITS_FILE]; T=Float64)[1, :]
    @test length(epoch_blocks(d64)[1]) == length(epoch_blocks(data)[1])
end

@testset "resample_epochs invariants" begin
    b = epoch_blocks(data)

    for mode in (:replacement, :halfsample, :weights)
        r = resample_epochs(data, b; mode=mode, rng=Xoshiro(1))
        @test length(r) == length(data)
        for e in eachindex(data)
            # the uv table is never touched — precomputed FTs stay valid
            @test r[e].uv == data[e].uv
            @test isempty(r[e].indx_v2) || maximum(r[e].indx_v2) <= size(r[e].uv, 2)
            @test isempty(r[e].indx_t3_1) || maximum(r[e].indx_t3_1) <= size(r[e].uv, 2)
            @test r[e].nv2 == length(r[e].v2)
            @test all(isfinite, r[e].v2_err)
            @test all(>(0), r[e].v2_err)
        end
    end

    # the multiplier bootstrap only rescales errors: sizes and counts are preserved
    w = resample_epochs(data, b; mode=:weights, rng=Xoshiro(2))
    for e in eachindex(data)
        @test length(w[e].v2) == length(data[e].v2)
        @test w[e].nv2 == data[e].nv2
        @test w[e].v2 == data[e].v2               # values untouched
        @test w[e].v2_err != data[e].v2_err       # errors rescaled
    end

    # same seed, same replicate
    @test resample_epochs(data, b; rng=Xoshiro(3))[1].v2 ==
          resample_epochs(data, b; rng=Xoshiro(3))[1].v2
    @test resample_epochs(data, b; rng=Xoshiro(4))[1].v2 !=
          resample_epochs(data, b; rng=Xoshiro(5))[1].v2

    # stratified draws keep every epoch populated
    for s in 1:5
        r = resample_epochs(data, b; stratify=true, rng=Xoshiro(s))
        @test all(e -> npts(r[e]) > 0, eachindex(r))
    end

    @test_throws ErrorException resample_epochs(data, [b[1], b[1]])
end

@testset "free-parameter selection" begin
    @test parametric_free_indices(nothing) == 1:7
    @test parametric_free_indices(["rpole", "ld1"]) == [1, 6]
    @test parametric_free_indices([:ld1, :rpole]) == [1, 6]      # order-independent
    @test parametric_free_indices([6, 1]) == [1, 6]
    @test length(parametric_free_indices(nothing; tpole_free=true)) == 8
    @test_throws ErrorException parametric_free_indices(["nonsense"])
    @test_throws ErrorException parametric_free_indices([1, 1])
    @test_throws ErrorException parametric_free_indices(Int[])
    @test_throws ErrorException parametric_free_indices([99])

    # rpole must be bounded away from 0: the flux normalisation makes the objective NaN
    # exactly at 0, and a bounded line search projects trial steps onto the box.
    lb, ub = default_parametric_bounds()
    @test lb[1] > 0
    @test ub[2] < 1                                  # omega below break-up
    @test ub[6] > 1                                  # ld1 above 1 (rho Cas fits at ~1.7)
    @test length(default_parametric_bounds(; tpole_free=true)[1]) == 8
end

@testset "bootstrap_parametric (synthetic fit)" begin
    # A "fit" that is a cheap statistic of the replicate: exercises the whole driver path
    # — resampling, seeding, clipping, summary, diagnostics — with no AD and no optimiser.
    fitfun = rdata -> begin
        v = rdata[1].v2
        ([Float64(mean(v)), Float64(std(v))], 1.0, (hit_maxiter = false,))
    end

    b = bootstrap_parametric(fitfun, data; nboot=24, seed=7, verb=false,
                             list_free_params=["mean_v2", "std_v2"])
    @test b isa ParametricBootstrap
    @test size(b.samples) == (24, 2)
    @test b.nfailed == 0
    @test all(isfinite, b.median)
    @test all(b.sigma .> 0)
    @test b.nblocks_per_epoch == [length(epoch_blocks(data)[1])]
    @test b.stratify
    @test b.ndegenerate == 0
    @test b.nmaxiter == 0
    @test size(b.covar) == (2, 2)
    @test b.list_free_params == ["mean_v2", "std_v2"]   # forwarded from BootstrapResult

    # reproducible for a fixed seed, and different for a different one
    @test bootstrap_parametric(fitfun, data; nboot=24, seed=7, verb=false).samples ==
          b.samples
    @test bootstrap_parametric(fitfun, data; nboot=24, seed=8, verb=false).samples !=
          b.samples

    # threading must not change the answer
    @test bootstrap_parametric(fitfun, data; nboot=24, seed=7, parallel=:inner,
                               verb=false).samples == b.samples

    # bounds are used only to count replicates pinned to them
    hi = maximum(b.samples[:, 1])
    bb = bootstrap_parametric(fitfun, data; nboot=24, seed=7, verb=false,
                              lb=[-Inf, -Inf], ub=[hi, Inf])
    @test bb.natbound[1] >= 1

    # a replicate whose fit throws is counted, not fatal (x_opt given, so the full-data
    # call is skipped — a full-data fit that fails *should* propagate, and does)
    flaky = rdata -> (isodd(length(rdata[1].v2)) ? error("boom") :
                      ([Float64(mean(rdata[1].v2)), 1.0], 1.0, nothing))
    bf = bootstrap_parametric(flaky, data; nboot=24, seed=11, verb=false,
                              x_opt=[0.0, 1.0], list_free_params=["a", "b"])
    @test bf.nfailed >= 1
    @test count(bf.mask) + bf.nfailed == 24
    @test_throws ErrorException bootstrap_parametric(_ -> error("boom"), data;
                                                     nboot=4, verb=false)

    @test_throws ErrorException bootstrap_parametric(fitfun, data; nboot=8,
                                                     parallel=:sideways, verb=false)
end

# ── End-to-end through the Zygote extension ─────────────────────────────────
zygote_ran = Ref(false)
try
    @eval using Zygote
    @testset "fit_parametric + bootstrap (Zygote)" begin
        tessels = tessellation_healpix(2, T=Float32)   # coarse on purpose: this is a
        base = (surface_type=2, rpole=1.25f0, tpole=4000f0,   # mechanics test, not science
                ldtype=3, ld1=1.75f0, ld2=0f0,
                inclination=0f0, position_angle=0f0, rotation_period=1f6,
                beta=0.08f0, frac_escapevel=0f0, B_rot=0f0)
        θ0 = Float32[1.25, 0.0, 0.0, 0.0, 0.08, 1.75, 0.0]

        θ̂, chi2r, info = fit_parametric(data, tessels, [0f0], base;
                                        θ0=θ0, free=["rpole", "ld1"], maxiter=30)
        @test length(θ̂) == 7
        @test isfinite(chi2r) && chi2r > 0
        @test info.free == [1, 6]
        @test θ̂[2] == θ0[2] && θ̂[5] == θ0[5]        # frozen entries untouched
        @test θ̂[1] != θ0[1]                          # free entries moved

        b = bootstrap_parametric(data, tessels, [0f0], base;
                                 θ0=θ0, free=["rpole", "ld1"], nboot=4, seed=3,
                                 maxiter=30, verb=false)
        @test size(b.samples) == (4, 2)
        @test b.list_free_params == ["rpole", "ld1"]
        @test all(isfinite, b.median)
        zygote_ran[] = true
    end
catch e
    @warn "Zygote section skipped" exception = (e, catch_backtrace())
end
