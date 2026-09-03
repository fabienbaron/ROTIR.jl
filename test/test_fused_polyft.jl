# The `@turbo` forward kernel against the scalar reference.
#
# `_cvis_turbo!` is a REWRITE of `_cvis_scalar!`: real accumulators instead of complex,
# `sin(πa)/(πa)` with an `ifelse` instead of `sinc`, `sin`/`cos` instead of `cis`, and the loop
# order inverted so the UV index is the vectorised one. Every one of those is a place the
# arithmetic could silently change, and one of them — the sign folded into the `−i` constant —
# did, in the first draft, producing a uniform relative error of exactly 2.0.
#
# So the reference is kept in the package and this asserts they agree: across both precisions,
# across mesh levels, across surface types (which is what changes the shape of the projected
# quads), and on the degenerate cases the reference special-cases.

using Test
using ROTIR
# `_cvis_turbo!` is declared in ROTIR but DEFINED in ROTIRLoopVectorizationExt, so the kernel
# under test does not exist until LoopVectorization is loaded — it is a weak dependency
# because loading it costs 1.8 s of GUI startup (see src/turbo_polyft.jl). The `import` below
# resolves either way, since the stub is ROTIR's; only the methods arrive with the extension.
using LoopVectorization
using ROTIR: _cvis_scalar!, _cvis_turbo!, POLYFT_BACKEND

@testset "fused polyft: turbo vs reference" begin
    # Which is to say: the extension is loaded. Without this the whole file would test a
    # stub, and `_cvis_turbo!` would fail with a MethodError rather than a wrong number.
    @test ROTIR.turbo_available()

    D = joinpath(pkgdir(ROTIR), "demos", "data")

    "Both kernels on one geometry, returning `(F_scalar, F_turbo)`."
    function both(T, n, params, file)
        data = readoifits(joinpath(D, file); verbose = false)[1, 1]
        tess = tessellation_healpix(n; T = T)
        star = create_star(tess, params, zero(T))
        idx = star.index_quads_visible
        pjx = Array(star.proj_west[idx, :]); pjy = Array(star.proj_north[idx, :])
        x = T.(parametric_temperature_map(params, star))
        xw = x[idx] .* (star.vis_weights[idx] .* star.ldmap[idx])
        kx = T.(data.uv[1, :]) * T(-π / (180 * 3600000))
        ky = T.(data.uv[2, :]) * T( π / (180 * 3600000))
        k2 = precompute_k2_inv_im(kx, ky)
        Fs = Vector{Complex{T}}(undef, length(kx))
        Ft = similar(Fs)
        _cvis_scalar!(Fs, kx, ky, k2, pjx, pjy, xw)
        _cvis_turbo!(Ft, kx, ky, k2, pjx, pjy, xw)
        return Fs, Ft
    end

    relerr(a, b) = maximum(abs.(a .- b)) / max(maximum(abs.(b)), eps(Float64))

    @testset "$(T), level $(n), surface_type $(st)" for T in (Float32, Float64),
                                                        n in (2, 3, 4),
                                                        st in (0, 1, 2)
        # A different surface type is a different projected quad SHAPE, which is what the
        # kernel actually sees: a sphere gives near-square quads, an ellipsoid stretched ones,
        # a rapid rotator strongly sheared ones near the equator.
        p = st == 0 ? default_star_params(0; T = T, radius = 3.2, tpole = 5000.0,
                                          ldtype = 1, ld1 = 0.3) :
            st == 1 ? default_star_params(1; T = T, radius_x = 3.6, radius_y = 3.0,
                                          radius_z = 2.8, inclination = 62.0) :
                      default_star_params(2; T = T, rpole = 3.0, frac_escapevel = 0.7,
                                          inclination = 71.0, beta = 0.15)
        Fs, Ft = both(T, n, p, "polaris.oifits")
        # The tolerance follows the PRECISION, not the level: this is a rearrangement of the
        # same sum, so the only difference is rounding order and the transcendental library.
        tol = T === Float32 ? 2e-5 : 1e-12
        @test relerr(Ft, Fs) < tol
        @test all(isfinite, Ft)
        @test length(Ft) == length(Fs)
    end

    @testset "other datasets: $(f)" for f in ("2011Sep02.lam_And_prepped.oifits",
                                              "2007_2012_2015.Spica.oifits")
        # A different uv distribution, and in Spica's case many more epochs folded into one
        # table — the kernel sees a different `kx`/`ky` spread, which is what its `1/(kx²+ky²)`
        # factor is sensitive to.
        p = default_star_params(0; T = Float32, radius = 1.0, tpole = 5000.0,
                                ldtype = 1, ld1 = 0.3)
        Fs, Ft = both(Float32, 3, p, f)
        @test relerr(Ft, Fs) < 2e-5
        @test all(isfinite, Ft)
    end

    @testset "degenerate inputs" begin
        T = Float64
        data = readoifits(joinpath(D, "polaris.oifits"); verbose = false)[1, 1]
        tess = tessellation_healpix(2; T = T)
        p = default_star_params(0; T = T, radius = 3.2, tpole = 5000.0, ldtype = 1, ld1 = 0.3)
        star = create_star(tess, p, zero(T))
        idx = star.index_quads_visible
        pjx = Array(star.proj_west[idx, :]); pjy = Array(star.proj_north[idx, :])
        kx = T.(data.uv[1, :]) * T(-π / (180 * 3600000))
        ky = T.(data.uv[2, :]) * T( π / (180 * 3600000))
        nvis = length(idx); nuv = length(kx)
        Fs = Vector{Complex{T}}(undef, nuv); Ft = similar(Fs)

        # A ZERO BASELINE. `precompute_k2_inv_im` returns 0 rather than Inf there, and the
        # rewrite carries that through `c = -imag(k2_inv_im)` rather than recomputing
        # `1/(kx²+ky²)` — which would have been Inf.
        kx0 = copy(kx); ky0 = copy(ky); kx0[1] = 0; ky0[1] = 0
        k20 = precompute_k2_inv_im(kx0, ky0)
        xw = ones(T, nvis)
        _cvis_scalar!(Fs, kx0, ky0, k20, pjx, pjy, xw)
        _cvis_turbo!(Ft, kx0, ky0, k20, pjx, pjy, xw)
        @test Fs[1] == 0 && Ft[1] == 0
        @test all(isfinite, Ft)
        @test relerr(Ft, Fs) < 1e-12

        # ZERO-WEIGHT tessels: the reference `continue`s past them, so the rewrite must skip
        # them the same way rather than adding 0·NaN.
        k2 = precompute_k2_inv_im(kx, ky)
        xz = zeros(T, nvis); xz[1:2:end] .= 1
        _cvis_scalar!(Fs, kx, ky, k2, pjx, pjy, xz)
        _cvis_turbo!(Ft, kx, ky, k2, pjx, pjy, xz)
        @test relerr(Ft, Fs) < 1e-12

        # ALL zero: both must give exactly zero, not NaN.
        _cvis_turbo!(Ft, kx, ky, k2, pjx, pjy, zeros(T, nvis))
        @test all(iszero, Ft)

        # A DEGENERATE quad (all four corners coincident) contributes nothing, and must not
        # produce NaN through the sinc singularity.
        pjd = copy(pjx); pjnd = copy(pjy)
        pjd[1, :] .= pjd[1, 1]; pjnd[1, :] .= pjnd[1, 1]
        _cvis_scalar!(Fs, kx, ky, k2, pjd, pjnd, ones(T, nvis))
        _cvis_turbo!(Ft, kx, ky, k2, pjd, pjnd, ones(T, nvis))
        @test all(isfinite, Ft)
        @test relerr(Ft, Fs) < 1e-12
    end

    @testset "the type-3 NUFFT backend" begin
        # A QUADRATURE, unlike the other two, so it is tested against them rather than the
        # other way round — and its accuracy depends on the phase span across one tessel,
        # which is why `ngauss` matters at coarse meshes and stops mattering at fine ones.
        T = Float64
        data = readoifits(joinpath(D, "polaris.oifits"); verbose = false)[1, 1]
        kx = T.(data.uv[1, :]) * T(-π / (180 * 3600000))
        ky = T.(data.uv[2, :]) * T( π / (180 * 3600000))
        k2 = precompute_k2_inv_im(kx, ky)
        @testset "level $(n)" for n in (2, 3, 4)
            tess = tessellation_healpix(n; T = T)
            p = default_star_params(0; T = T, radius = 3.2, tpole = 5000.0,
                                    ldtype = 1, ld1 = 0.3)
            star = create_star(tess, p, zero(T))
            idx = star.index_quads_visible
            pjx = Matrix(star.proj_west[idx, :]); pjy = Matrix(star.proj_north[idx, :])
            xw = T.(parametric_temperature_map(p, star)[idx] .*
                    star.vis_weights[idx] .* star.ldmap[idx])
            Fex = Vector{Complex{T}}(undef, length(kx))
            _cvis_turbo!(Fex, kx, ky, k2, pjx, pjy, xw)
            nrm = maximum(abs, Fex)
            F4 = polyft_cvis_nufft(pjx, pjy, xw, kx, ky; ngauss = 4)
            @test maximum(abs.(F4 .- Fex)) / nrm < 1e-5
            @test all(isfinite, F4)
            @test length(F4) == length(Fex)
            # Refining the rule must IMPROVE it — that is what distinguishes a quadrature
            # from the rasterised route, whose error floored at 5e-3 whatever the grid.
            e2 = maximum(abs.(polyft_cvis_nufft(pjx, pjy, xw, kx, ky; ngauss = 2) .- Fex)) / nrm
            e4 = maximum(abs.(F4 .- Fex)) / nrm
            @test e4 <= e2
        end
    end

    @testset "observables pick the route from the star" begin
        # `observables` reads `star.polyft` when `setup_oi!` has filled it and computes the
        # visibilities matrix-free when it has not. That is what lets imaging keep the matrix
        # (fixed geometry, many maps) while a χ² table skips it (geometry changed with the
        # parameters, one evaluation) — with no caller having to say which.
        data = [readoifits(joinpath(D, "polaris.oifits"); verbose = false)[1, 1]]
        tess = tessellation_healpix(3)
        p = default_star_params(0; radius = 3.2, tpole = 5000.0, ldtype = 1, ld1 = 0.3)

        free = create_star_multiepochs(tess, p, [0.0])
        @test isempty(free[1].polyft)                 # nothing built it
        x = parametric_temperature_map(p, free[1])
        bf = chi2_breakdown(x, free, data)

        dense = create_star_multiepochs(tess, p, [0.0])
        setup_oi!(data, dense)
        @test !isempty(dense[1].polyft)
        bd = chi2_breakdown(x, dense, data)

        # The two routes must give the same χ², per observable, not just in total.
        @test abs(bf.total - bd.total) / bd.total < 1e-6
        @test abs(bf.v2    - bd.v2)    / bd.v2    < 1e-6
        @test abs(bf.t3amp - bd.t3amp) / bd.t3amp < 1e-6
        @test abs(bf.t3phi - bd.t3phi) / bd.t3phi < 1e-6
        @test bf.ndata == bd.ndata
        # And so must the observables themselves.
        vf, af, pf_ = observables(x, free[1], data[1])
        vd, ad, pd  = observables(x, dense[1], data[1])
        @test maximum(abs.(vf .- vd)) / maximum(abs, vd) < 1e-6
        @test maximum(abs.(af .- ad)) / maximum(abs, ad) < 1e-6
    end

    @testset "the backend switch reaches the χ²" begin
        # The switch has to change which kernel runs and nothing else about the answer.
        data = [readoifits(joinpath(D, "polaris.oifits"); verbose = false)[1, 1]]
        tess = tessellation_healpix(3)
        p = default_star_params(0; radius = 3.2, tpole = 5000.0, ldtype = 1, ld1 = 0.3)
        old = POLYFT_BACKEND[]
        try
            POLYFT_BACKEND[] = :scalar
            c_scalar = parametric_chi2(p, tess, data, [0.0])
            @test isfinite(c_scalar) && c_scalar > 0
            # All three must agree on the χ², which is the only thing that makes offering a
            # choice safe: a backend that is fast and slightly wrong would bias every fit.
            for b in (:turbo, :nufft)
                POLYFT_BACKEND[] = b
                c = parametric_chi2(p, tess, data, [0.0])
                @test abs(c - c_scalar) / c_scalar < 1e-4
            end
        finally
            POLYFT_BACKEND[] = old
        end
    end
end
