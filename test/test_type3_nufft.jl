# ROTIR's own 2-D type-3 NUFFT (src/type3_nufft.jl).
#
# THE TEST THAT MATTERS MOST IS THE DOT PRODUCT. The forward can be checked against the exact
# closed-form kernel, but an adjoint has no independent reference — so it is checked against
# the forward it is supposed to be the transpose of:
#
#     Re( Σ_k adj[k] · F[k] )  ==  Σ_p xw[p] · grad[p]
#
# which holds to machine precision for a true transpose and fails immediately for anything
# else. It is the assertion that caught this adjoint reading uninitialised memory after the
# stencil buffers moved into tuples: the timings looked BETTER, the forward was untouched, and
# only the dot product said anything was wrong.

using Test, ROTIR, LinearAlgebra, Random
# `:nufft` is compared against here, and it lives behind ROTIRFINUFFTExt.
using FINUFFT
using ROTIR: POLYFT_BACKEND, compute_polyflux_and_cvis!, compute_adjoint_cvis!,
             precompute_k2_inv_im, TYPE3_PLANS, t3_extent_bin

@testset "the type-3 NUFFT" begin
    D = joinpath(pkgdir(ROTIR), "demos", "data")
    d = readoifits(joinpath(D, "polaris.oifits"); verbose = false)[1, 1]

    "Geometry, weights and uv coordinates in `T`, plus the visible-quad count."
    function bits(nexp, T)
        tess = tessellation_healpix(nexp; T = T)
        star = create_star(tess, default_star_params(0; radius = 3.2, tpole = 5000.0,
                                                     ldtype = 1, ld1 = 0.3), 0.0)
        i = star.index_quads_visible
        (Matrix(star.proj_west[i, :]), Matrix(star.proj_north[i, :]),
         T.(star.vis_weights[i] .* star.ldmap[i]),
         T.(d.uv[1, :]) * T(-π / (180 * 3600000)),
         T.(d.uv[2, :]) * T(π / (180 * 3600000)), length(i))
    end

    "The exact closed-form forward and adjoint, in Float64 — the reference for both precisions."
    function reference(nexp)
        pw, pn, xw, kx, ky, npx = bits(nexp, Float64)
        k2 = precompute_k2_inv_im(kx, ky); J = length(kx)
        F = Vector{ComplexF64}(undef, J); pf = zeros(npx); g = Vector{Float64}(undef, npx)
        adj = ComplexF64.(range(-1, 1, length = J), range(0.7, -0.4, length = J))
        old = POLYFT_BACKEND[]
        try
            POLYFT_BACKEND[] = :scalar
            compute_polyflux_and_cvis!(F, pf, kx, ky, k2, pw, pn, xw)
            compute_adjoint_cvis!(g, adj, kx, ky, k2, pw, pn, fill(1.0, npx))
        finally
            POLYFT_BACKEND[] = old
        end
        return F, g, adj
    end
    rel(a, b) = maximum(abs, ComplexF64.(a) .- b) / maximum(abs, b)
    relr(a, b) = maximum(abs, Float64.(a) .- b) / maximum(abs, b)

    @testset "the forward against the exact closed form" begin
        for nexp in (3, 4)
            Fe, _, _ = reference(nexp)
            for (T, tol) in ((Float64, 1e-9), (Float32, 3e-5))
                pw, pn, xw, kx, ky, _ = bits(nexp, T)
                F = type3_cvis(pw, pn, xw, kx, ky)
                # THE PLAN IS Float64 WHATEVER THE MESH IS — single buys ~10 % and costs four
                # orders of accuracy — so the transform returns double and the caller narrows.
                @test eltype(F) === ComplexF64
                # MEASURED at 2.5e-11 (HEALPix 3) and 3.3e-10 (4) in double, against FINUFFT's
                # 5.8e-10 at its own 1e-9 tolerance — so this route is the more accurate one.
                @test rel(F, Fe) < tol
                @test all(isfinite, F)
            end
        end
    end

    @testset "the adjoint is the transpose of the forward" begin
        for nexp in (3, 4), T in (Float64, Float32)
            _, ge, adj64 = reference(nexp)
            pw, pn, xw, kx, ky, npx = bits(nexp, T)
            adj = Complex{T}.(adj64)
            F = type3_cvis(pw, pn, xw, kx, ky)
            g = Vector{T}(undef, npx)
            type3_cvis_adj!(g, pw, pn, adj, kx, ky)
            @test eltype(g) === T                 # narrowed into the caller's buffer
            # THE DOT TEST. Exact to machine precision in double; in single it can only hold to
            # Float32 eps, which is itself worth knowing — a gradient check run in single
            # precision cannot resolve better than that.
            lhs = real(sum(adj[k] * F[k] for k in eachindex(F)))
            rhs = dot(xw, g)
            @test abs(lhs - rhs) / abs(rhs) < (T === Float64 ? 1e-12 : 1e-4)
            # A GRID THAT ACTUALLY HOLDS THE STENCIL. `nf` is floored at `w/(1 − 1/γ)` because
            # the spread loop has no bounds check; a geometry with a small space-bandwidth
            # product used to give nf = 16 against w = 11 and index off the end of the grid.
            pl = type3_plan_for(pw, pn, kx, ky)
            @test pl.nf >= 2 * pl.w
            # And against the EXACT operator's adjoint, which differs by the quadrature.
            @test relr(g, ge) < (T === Float64 ? 1e-7 : 1e-3)
            @test !all(iszero, g)
        end
    end

    @testset "type3_points! against a direct sum" begin
        # The generic point-set entry point, checked against the definition rather than against
        # another approximation — this is what pins the sign convention and the scaling.
        Random.seed!(3)
        A = 3.2; N = 3000; J = 200
        xs = A * (2rand(N) .- 1); ys = A * (2rand(N) .- 1); cs = randn(N)
        kx = 2rand(J) .- 1; ky = 2rand(J) .- 1
        direct = [sum(cs[l] * cis(-2π * (kx[j]*xs[l] + ky[j]*ys[l])) for l in 1:N)
                  for j in 1:J]
        p = plan_type3(A, kx, ky; w = 13)
        out = Vector{ComplexF64}(undef, J)
        type3_points!(out, p, xs, ys, cs)
        @test maximum(abs, out .- direct) / maximum(abs, direct) < 1e-9
    end

    @testset "the plan is cached, and re-made when the star outgrows it" begin
        pw, pn, xw, kx, ky, _ = bits(3, Float64)
        empty!(TYPE3_PLANS)
        p1 = type3_plan_for(pw, pn, kx, ky)
        p2 = type3_plan_for(pw, pn, kx, ky)
        @test p1 === p2                      # same geometry, same plan object
        @test length(TYPE3_PLANS) == 1
        # SHRINKING reuses it — the extent bin is an upper bound, so a smaller star still fits.
        p3 = type3_plan_for(pw .* 0.96, pn .* 0.96, kx, ky)
        @test p3 === p1
        # GROWING past the bin gets a new one, which is the point: MEASURED, a plan built at
        # radius 3.2 degrades to 7.0e-8 at 4.0 and 1.5e-3 at 5.6 because the k-side
        # oversampling falls below 2.
        p4 = type3_plan_for(pw .* 2, pn .* 2, kx, ky)
        @test p4 !== p1
        @test length(TYPE3_PLANS) == 2
        # Different uv targets are a different plan even at the same extent.
        @test type3_plan_for(pw, pn, kx .* 1.3, ky) !== p1
        # The bin is an upper bound on the extent it was asked for, and monotone.
        @test all(t3_extent_bin(a) >= a for a in (0.3, 1.0, 3.2, 3.3, 7.9, 12.0))
        @test t3_extent_bin(3.2) <= t3_extent_bin(4.0) <= t3_extent_bin(5.6)
        empty!(TYPE3_PLANS)
    end

    @testset "the quadrature rule prefers order to subdivision" begin
        pw, pn, _, kx, ky, _ = bits(3, Float64)
        ng3, ns3 = quadrature_for_type3(pw, pn, kx, ky)
        ng0, ns0 = quadrature_for(pw, pn, kx, ky)
        # The existing rule holds ngauss at 4 and subdivides; this one raises the order. At
        # HEALPix 3 that is 6/1 against 4/3 — four times fewer nodes for the same accuracy,
        # because Gauss-Legendre is spectrally accurate on an analytic integrand and
        # subdivision is only algebraic.
        @test ns3 == 1
        @test ng3 > ng0
        @test ng3 * ns3 < ng0 * ns0
        # A finer mesh needs less: the phase span across a quad shrinks with the quad.
        ng5, ns5 = quadrature_for_type3(bits(5, Float64)[1], bits(5, Float64)[2], kx, ky)
        @test ng5 <= ng3
        en, ew = t3_gauss_rule(ng3, ns3, Float64)
        @test length(en) == length(ew) == ng3 * ns3
        @test sum(ew) ≈ 2 rtol=1e-12          # the rule integrates 1 over [-1,1]
        @test_throws ErrorException t3_gauss_rule(7, 1, Float64)   # no 7-point table
        # SO THE RULE MUST NEVER ASK FOR ONE. `_GL_NODES` has no 7- or 9-point entry, and the
        # raw width formula can land on either — which errored the moment a geometry outside
        # the benchmark set came along (the precompile workload's synthetic uv). Every order
        # the rule can return must be buildable, over the whole plausible span range.
        @test all(o in keys(ROTIR._GL_NODES) for o in ROTIR.T3_GAUSS_ORDERS)
        let pwq = [0.0 1.0 1.0 0.0], pnq = [0.0 0.0 1.0 1.0]
            for s in exp10.(range(-3, 1.6, length = 60))
                ng, ns = quadrature_for_type3(pwq .* s, pnq .* s, [1.0, 0.5], [0.3, -0.7])
                @test ng in keys(ROTIR._GL_NODES)
                @test ns >= 1
                @test length(first(t3_gauss_rule(ng, ns, Float64))) == ng * ns
            end
        end
    end

    @testset "the 1/ψ̂ fit is not the limiting term" begin
        _, _, _, kx, ky, _ = bits(3, Float64)
        # A DEGENERATE-COVERAGE CASE, which is what exposed the missing grid floor: very short
        # baselines make the space-bandwidth product tiny and `nf` collapses.
        for (kxs, kys) in ((kx, ky), (kx .* 1e-3, ky .* 1e-3))
            pl = plan_type3(3.2, kxs, kys; w = 11)
            @test pl.nf >= 2 * pl.w
            out = Vector{ComplexF64}(undef, length(kxs))
            @test all(isfinite, type3_points!(out, pl, [0.4, -1.2], [0.7, 2.0], [1.0, 0.5]))
        end
        for T in (Float64, Float32)
            p = plan_type3(3.2, kx, ky; T = T, w = 11)
            e, _ = type3_psihat_fit_error(p)
            # MEASURED 1.3e-13 in double and 7.1e-7 in single, both under the transform's own
            # error at that precision. The fit and its evaluation are in Float64 whatever `T`
            # is, which is what buys the single-precision figure — in Float32 Clenshaw floors
            # at 1.4e-5 and would dominate.
            @test e < (T === Float64 ? 1e-11 : 5e-6)
        end
    end

    @testset ":t3 reaches the χ² and the adjoint through the backend switch" begin
        data = [readoifits(joinpath(D, "polaris.oifits"); verbose = false)[1, 1]]
        tess = tessellation_healpix(3; T = Float64)
        sp = default_star_params(0; radius = 3.2, tpole = 5000.0, ldtype = 1, ld1 = 0.3)
        old = POLYFT_BACKEND[]
        try
            POLYFT_BACKEND[] = :scalar
            c_exact = parametric_chi2(sp, tess, data, [0.0])
            POLYFT_BACKEND[] = :t3
            c_t3 = parametric_chi2(sp, tess, data, [0.0])
            @test abs(c_t3 - c_exact) / c_exact < 1e-6
            # AND THE ADJOINT BRANCH, which `:nufft` does not have — the reason this exists.
            pw, pn, xw, kx, ky, npx = bits(3, Float64)
            k2 = precompute_k2_inv_im(kx, ky); J = length(kx)
            adj = ComplexF64.(range(-1, 1, length = J), range(0.7, -0.4, length = J))
            one_ = fill(1.0, npx)
            gs = Vector{Float64}(undef, npx); gt = Vector{Float64}(undef, npx)
            POLYFT_BACKEND[] = :scalar
            compute_adjoint_cvis!(gs, adj, kx, ky, k2, pw, pn, one_)
            POLYFT_BACKEND[] = :t3
            compute_adjoint_cvis!(gt, adj, kx, ky, k2, pw, pn, one_)
            @test relr(gt, gs) < 1e-7
            @test !all(iszero, gt)
        finally
            POLYFT_BACKEND[] = old
        end
    end
end
