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
# FINUFFT is a weak dependency too now — `:nufft` needs ROTIRFINUFFTExt, exactly as `:turbo`
# needs ROTIRLoopVectorizationExt. It is in the test target for this reason.
using FINUFFT
using ROTIR: _cvis_scalar!, _cvis_turbo!, POLYFT_BACKEND,
             compute_adjoint_cvis!, compute_adjoint_vertices!,
             compute_polyflux_and_cvis!

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

    # THE TWO ADJOINTS, which are where a gradient's time actually goes. MEASURED per
    # evaluation on one lam And epoch, scalar against turbo:
    #
    #     nside   forward          adj_cvis         adj_verts        all three
    #       3     1.5 -> 0.2 ms    3.6 -> 0.2 ms    6.3 -> 0.2 ms    11.4 -> 0.7 ms  16.3x
    #       4     5.3 -> 0.7       9.0 -> 0.5      18.1 -> 0.9       32.3 -> 2.2     15.0x
    #       5    18.1 -> 2.4      39.0 -> 3.0      67.8 -> 3.6      124.9 -> 9.0     13.8x
    #
    # They vectorise BETTER than the forward (13-25x against 7x) because of a structural
    # difference: in the forward every `F[k]` is touched by every tessel, which forces chunked
    # per-thread accumulators; in both adjoints each `p` writes only its own outputs, so the
    # thread split over `p` is free and `k` is a clean vectorised reduction inside it.
    #
    # The scalar kernel is the DEFINITION and the turbo one a rewrite in real arithmetic —
    # `@turbo` will not take `Complex` — so this is the only thing standing between that
    # algebra and a silently wrong gradient.
    @testset "adjoints, $(T), level $(n)" for T in (Float32, Float64), n in (2, 3)
        data = readoifits(joinpath(D, "2011Sep02.lam_And_prepped.oifits");
                          verbose = false)[1, 1]
        tess = tessellation_healpix(n; T = T)
        prm  = default_star_params(2; T = T, rpole = 3.0, frac_escapevel = 0.7,
                                   tpole = 5000.0, ldtype = 1, ld1 = 0.3)
        star = create_star(tess, prm, zero(T))
        idx  = star.index_quads_visible
        pjx  = Array(star.proj_west[idx, :]); pjy = Array(star.proj_north[idx, :])
        x    = T.(parametric_temperature_map(prm, star))
        xw   = x[idx] .* (star.vis_weights[idx] .* star.ldmap[idx])
        kx   = T.(data.uv[1, :]) * T(-π / (180 * 3600000))
        ky   = T.(data.uv[2, :]) * T( π / (180 * 3600000))
        k2   = precompute_k2_inv_im(kx, ky)
        nuv  = length(kx); npx = length(xw)
        F    = Vector{Complex{T}}(undef, nuv); pf = zeros(T, npx)
        compute_polyflux_and_cvis!(F, pf, kx, ky, k2, pjx, pjy, xw)
        # A COMPLEX cotangent with both parts populated: a real-only one would leave half the
        # expanded algebra untested, and it is the imaginary half that carries the sign work.
        adj = Complex{T}.(range(T(-1), T(1), length = nuv),
                          range(T(0.7), T(-0.4), length = nuv))
        gs = Vector{T}(undef, npx); gt = Vector{T}(undef, npx)
        ws = zeros(T, npx, 4); ns = zeros(T, npx, 4)
        wt = zeros(T, npx, 4); nt = zeros(T, npx, 4)
        old = POLYFT_BACKEND[]
        try
            POLYFT_BACKEND[] = :scalar
            compute_adjoint_cvis!(gs, adj, kx, ky, k2, pjx, pjy, pf)
            compute_adjoint_vertices!(ws, ns, adj, kx, ky, k2, pjx, pjy, xw, pf)
            POLYFT_BACKEND[] = :turbo
            compute_adjoint_cvis!(gt, adj, kx, ky, k2, pjx, pjy, pf)
            compute_adjoint_vertices!(wt, nt, adj, kx, ky, k2, pjx, pjy, xw, pf)
        finally
            POLYFT_BACKEND[] = old
        end
        # The tolerance is the float type's own floor for a reduction this long, not the
        # kernel's: measured 1.9e-14 at Float64 and 6.5e-6 at Float32 for `adj_cvis`.
        tol = T === Float64 ? 1e-10 : 2e-3
        @test relerr(gt, gs) < tol
        @test relerr(wt, ws) < tol
        @test relerr(nt, ns) < tol
        # And not trivially zero, which would pass every comparison above.
        @test maximum(abs, gs) > 0 && maximum(abs, ws) > 0 && maximum(abs, ns) > 0
    end

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

    @testset "the type-3 kernel follows the mesh precision" begin
        # IT USED TO PROMOTE. `polyft_cvis_nufft` built its samples with a hardcoded
        # `T = Float64`, cast to `ComplexF64`, passed `Float64(tol)` and wrote its targets as
        # `2π .* kx` — and `2π` is a Float64, so even a Float32 `kx` came out Float64. The
        # whole default kernel therefore ran in double and narrowed on return, whatever the
        # panel's precision box said.
        data = [readoifits(joinpath(D, "polaris.oifits"); verbose = false)[1, 1]]
        p = default_star_params(0; radius = 3.2, tpole = 5000.0, ldtype = 1, ld1 = 0.3)

        "The quad geometry, the weights and the uv, all in `T`."
        function inputs(T)
            tess = tessellation_healpix(3; T = T)
            star = create_star(tess, p, 0.0)
            i = star.index_quads_visible
            xw = T.(star.vis_weights[i] .* star.ldmap[i])
            (Matrix(star.proj_west[i, :]), Matrix(star.proj_north[i, :]), xw,
             T.(data[1].uv[1, :]) * T(-π / (180 * 3600000)),
             T.(data[1].uv[2, :]) * T(π / (180 * 3600000)))
        end

        @test nufft_work_type(zeros(Float32, 2, 4), zeros(Float32, 3)) === Float32
        @test nufft_work_type(zeros(Float64, 2, 4), zeros(Float64, 3)) === Float64
        # A MIXED pair is Float64: the transform is only as good as its worse argument, and
        # silently dropping the Float64 half to single would be a downgrade the caller did not
        # ask for.
        @test nufft_work_type(zeros(Float32, 2, 4), zeros(Float64, 3)) === Float64
        # And a type FINUFFT cannot take at all resolves to Float64 rather than failing inside
        # the C call — it supports exactly Float32 and Float64.
        @test nufft_work_type(zeros(Float16, 2, 4), zeros(Float16, 3)) === Float64
        @test nufft_tol(Float32) == 1.0e-6
        @test nufft_tol(Float64) == 1.0e-9

        a32 = inputs(Float32); a64 = inputs(Float64)
        # THE DEFAULT IS DOUBLE, deliberately and whatever the mesh is — see the docstring for
        # the measurement. Single is opt-in, and asking for it is what this checks.
        @test eltype(polyft_cvis_nufft(a32...)) === ComplexF64
        F32 = polyft_cvis_nufft(a32...; T = Float32)
        F64 = polyft_cvis_nufft(a64...)
        @test eltype(F32) === ComplexF32
        @test eltype(F64) === ComplexF64
        # And following the inputs is one keyword away for a caller that wants it.
        @test eltype(polyft_cvis_nufft(a32...; T = nufft_work_type(a32[1], a32[4]))) === ComplexF32
        # The pinned reference stays double whatever it is handed, which is what makes it
        # usable as the thing the single path is measured against.
        @test eltype(polyft_cvis_nufft_f64(a32...)) === ComplexF64

        # ACCURACY AGAINST THE EXACT KERNEL, which is neither of these: `:scalar` evaluates the
        # closed-form polygon transform with no quadrature and no tolerance.
        exact = let (pw, pn, xw, kx, ky) = a64
            F = Vector{ComplexF64}(undef, length(kx)); pf = zeros(length(xw))
            old = POLYFT_BACKEND[]
            try
                POLYFT_BACKEND[] = :scalar
                compute_polyflux_and_cvis!(F, pf, kx, ky, ROTIR.precompute_k2_inv_im(kx, ky),
                                           pw, pn, xw)
            finally
                POLYFT_BACKEND[] = old
            end
            F
        end
        err(F) = maximum(abs, ComplexF64.(F) .- exact) / maximum(abs, exact)
        # MEASURED: 5.8e-10 double against 1.9e-6 single on lam And and 1.0e-5 here on
        # polaris. The single figure is the transform's own floor rather than slack in the
        # bound — `nufft_tol(Float32)` is already FINUFFT's single-precision limit, and double
        # at that same 1e-6 gives 7.2e-8. Note it is LARGER than the quadrature's 6.8e-7, so
        # in single precision the transform becomes the limiting term rather than the rule,
        # which is the measurement behind not making it the default.
        @test err(F64) < 1e-8
        @test err(F32) < 3e-5
        @test err(F32) > 100 * err(F64)
    end

    @testset "with_polyft_backend scopes the kernel" begin
        # WHY THE SCOPE EXISTS. The kernel a fit wants is not the one the rest of the process
        # wants: `fused_cvis` has a `:nufft` branch and `interferometric_chi2` does not, so a
        # gradient fit needs `:turbo` while the interactive χ² beside it keeps `:nufft`. The
        # GUI runs the fit on a worker thread while its event loop keeps answering on another,
        # so this cannot be a global assignment: that would change the objective function
        # under a running line search.
        @test polyft_backend() === POLYFT_BACKEND[]          # no scope: the Ref answers
        old = POLYFT_BACKEND[]
        try
            POLYFT_BACKEND[] = :nufft
            with_polyft_backend(:turbo) do
                @test polyft_backend() === :turbo
                @test POLYFT_BACKEND[] === :nufft            # the Ref is NOT written
                # A task started inside inherits it, which is what makes the `Threads.@threads`
                # loops in the two adjoint kernels see the same choice as their caller.
                @test fetch(Threads.@spawn polyft_backend()) === :turbo
            end
            @test polyft_backend() === :nufft                # and it is restored on exit
            # Restored THROUGH A THROW as well; a fit that errors must not leave the kernel
            # changed for the rest of the session.
            @test_throws ErrorException with_polyft_backend(:scalar) do
                @test polyft_backend() === :scalar
                error("boom")
            end
            @test polyft_backend() === :nufft

            # And the scope reaches the kernels, not just the accessor: the same χ² under a
            # scoped `:turbo` as under a globally assigned one.
            data = [readoifits(joinpath(D, "polaris.oifits"); verbose = false)[1, 1]]
            # FLOAT64, deliberately. At Float32 this χ² is 1.11e7, so the quadrature kernel's
            # 6.8e-7 relative difference from the exact one is below the last bit and the two
            # come back bit-identical — which makes the "they differ" check below a tautology
            # that passes for the wrong reason.
            tess = tessellation_healpix(3; T = Float64)
            p = default_star_params(0; radius = 3.2, tpole = 5000.0, ldtype = 1, ld1 = 0.3)
            POLYFT_BACKEND[] = :turbo
            c_global = parametric_chi2(p, tess, data, [0.0])
            POLYFT_BACKEND[] = :nufft
            c_scoped = with_polyft_backend(:turbo) do
                parametric_chi2(p, tess, data, [0.0])
            end
            @test c_scoped == c_global
            # Not a tautology: the unscoped call gives the OTHER kernel's number, which agrees
            # to the quadrature's accuracy rather than exactly.
            c_plain = parametric_chi2(p, tess, data, [0.0])
            @test c_plain != c_global
            @test 0 < abs(c_plain - c_global) / c_global < 1e-4
        finally
            POLYFT_BACKEND[] = old
        end
    end
end
