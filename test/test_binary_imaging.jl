# The binary IMAGING gradient: `binary_chi2_fg` and the criterion built on it.
#
# What makes this worth pinning rather than trusting: the gradient is `ℜ(Jᵀ·a)` with a complex
# adjoint, and getting the conjugation backwards is not uniformly wrong. It leaves flux-like
# directions — brighten a tessel — nearly right, while inverting phase-like ones, which here
# means everything the SEPARATION carries. A reconstruction with the sign flipped on that half
# still descends, to the wrong map. Finite differences along random directions catch it; a
# convergence check does not.
#
# Two independent references, because each one alone has a blind spot:
#   * FINITE DIFFERENCES on the full concatenated map — the general case, but conditioned by
#     a χ² of ~10⁵ against directional derivatives of ~10⁻², so it is good to ~1e-6.
#   * REDUCTION TO ONE STAR — with the secondary's map at zero the pair IS the primary, and
#     `spheroid_chi2_fg` is the independently written expression to compare against. Exact to
#     the accumulation-order difference between the two χ² summations, so it pins the
#     conjugation and the normalization to 1e-8 where FD cannot.

using Test
using LinearAlgebra
using FiniteDifferences
using Random
using ROTIR
using ROTIR: binary_chi2_f, binary_chi2_fg, binary_crit_allepochs_fg, binary_reconstruct_oi,
             binary_phase_shift, split_binary_map, spheroid_chi2_fg, spheroid_chi2_f,
             spheroid_regularization, sobel_gradient_healpix

@testset "binary imaging gradient" begin
    D = joinpath(pkgdir(ROTIR), "demos", "data")
    # Two epochs from ONE file: the criterion sums over epochs and weights them, and a
    # single-epoch test cannot tell a sum from an assignment.
    raw = readoifits(joinpath(D, "2007_2012_2015.Spica.oifits"); verbose = false)[1, 1]
    data = [raw, raw]
    tess = tessellation_healpix(2)          # 192 tessels per component; FD affordable
    p1 = default_star_params(3; rpole = 0.447, tpole = 25300.0, q = 0.6188)
    p2 = default_star_params(3; rpole = 0.227, tpole = 20585.0, q = 1 / 0.6188)
    te = zeros(length(data))

    stars1 = create_star_multiepochs(tess, p1, te; secondary = false)
    stars2 = create_star_multiepochs(tess, p2, te; secondary = true)
    setup_oi!(data, stars1)
    setup_oi!(data, stars2)
    x1 = Float64.(parametric_temperature_map(p1, stars1[1]; secondary = false))
    x2 = Float64.(parametric_temperature_map(p2, stars2[1]; secondary = true))
    n1 = length(x1); n2 = length(x2)
    # Two different separations, so an epoch's phase shift cannot be silently reused for the
    # other — which is the mistake a single-epoch or equal-offset test would let through.
    shifts = [binary_phase_shift(data[1].uv, 0.6, 0.4),
              binary_phase_shift(data[2].uv, -0.3, 0.9)]

    @testset "value agrees with binary_chi2_f" begin
        g1 = zeros(n1); g2 = zeros(n2)
        f = binary_chi2_fg(x1, g1, stars1[1], x2, g2, stars2[1], data[1], shifts[1])
        ref = binary_chi2_f(x1, stars1[1], x2, stars2[1], data[1], shifts[1])
        # Not equal to the last bit: the two sum the same three terms in a different order
        # over ~10⁴ points at χ² ~ 10⁵.
        @test isapprox(f, ref; rtol = 1e-7)
        @test all(isfinite, g1) && all(isfinite, g2)
        # A tessel the epoch cannot see does not enter the transform and must not pick up a
        # gradient, or the optimiser will push on the far side of the star.
        hidden1 = setdiff(1:n1, stars1[1].index_quads_visible)
        @test isempty(hidden1) || all(iszero, g1[hidden1])
    end

    @testset "finite differences, both components" begin
        # The criterion as a function of the concatenated map, through the SAME forward path
        # the gradient differentiates — value-only, so this is not testing the fg against
        # itself but against the model it claims to differentiate.
        function f_of(v)
            a = v[1:n1]; b = v[n1+1:end]
            ga = zeros(n1); gb = zeros(n2)
            binary_chi2_fg(a, ga, stars1[1], b, gb, stars2[1], data[1], shifts[1])
        end
        g1 = zeros(n1); g2 = zeros(n2)
        binary_chi2_fg(x1, g1, stars1[1], x2, g2, stars2[1], data[1], shifts[1])
        g = vcat(g1, g2)
        x = vcat(x1, x2)

        Random.seed!(20260904)
        fdm = central_fdm(5, 1)
        # Directional derivatives, not the full 384-column Jacobian: the direction is what
        # mixes the two components, and a random direction with a component on each half
        # exercises the shared `S` term that a coordinate-wise sweep also would, at 5
        # evaluations instead of 1920.
        for _ in 1:6
            # Scaled to the map's own magnitude (temperatures ~2e4). A unit-norm direction
            # makes the FD step relatively tiny and the difference of two ~1e5 numbers
            # cancels to noise.
            v = randn(n1 + n2); v .*= norm(x) / norm(v)
            an = dot(g, v)
            fd = fdm(t -> f_of(x .+ t .* v), 0.0)
            @test isapprox(an, fd; rtol = 2e-5)
        end
        # And one direction on each component alone, so a gradient that is right on the pair
        # only because the two halves' errors cancel does not pass.
        for half in (1:n1, n1+1:n1+n2)
            v = zeros(n1 + n2); v[half] .= randn(length(half))
            v .*= norm(x) / norm(v)
            @test isapprox(dot(g, v), fdm(t -> f_of(x .+ t .* v), 0.0); rtol = 2e-5)
        end
    end

    @testset "reduces to the single star when the secondary is dark" begin
        # A zero secondary contributes no flux and no transform, so the pair IS the primary
        # and the primary's gradient must be the one `spheroid_chi2_fg` writes — a completely
        # separate expression, with its adjoint written out by hand.
        z2 = zeros(n2)
        g1 = zeros(n1); g2 = zeros(n2)
        f = binary_chi2_fg(x1, g1, stars1[1], z2, g2, stars2[1], data[1], shifts[1])
        gref = zeros(n1)
        fref = spheroid_chi2_fg(x1, gref, stars1[1], data[1]; verbose = false)
        @test isapprox(f, fref; rtol = 1e-7)
        @test isapprox(g1, gref; rtol = 1e-6, atol = 1e-9 * maximum(abs, gref))
        # The secondary still has a gradient: lighting it up from zero changes the χ².
        @test any(!iszero, g2)
    end

    @testset "criterion over epochs" begin
        x = vcat(x1, x2)
        g = zeros(n1 + n2)
        f = binary_crit_allepochs_fg(x, g, stars1, stars2, data, shifts)
        perepoch = sum(binary_chi2_f(x1, stars1[i], x2, stars2[i], data[i], shifts[i])
                       for i in eachindex(data))
        @test isapprox(f, perepoch; rtol = 1e-7)

        # Per-epoch weights go on the terms AND on the gradient, or the map is optimised
        # against a different criterion from the one reported.
        gw = zeros(n1 + n2)
        fw = binary_crit_allepochs_fg(x, gw, stars1, stars2, data, shifts;
                                      epochs_weights = [2.0, 0.0])
        g1only = zeros(n1 + n2)
        f1only = binary_crit_allepochs_fg(x, g1only, stars1[1:1], stars2[1:1], data[1:1],
                                          shifts[1:1])
        @test isapprox(fw, 2 * f1only; rtol = 1e-10)
        @test isapprox(gw, 2 .* g1only; rtol = 1e-10)

        @test_throws DimensionMismatch binary_crit_allepochs_fg(x[1:end-1], g[1:end-1], stars1,
                                                                stars2, data, shifts)
        @test_throws DimensionMismatch binary_crit_allepochs_fg(x, g, stars1, stars2, data,
                                                                shifts[1:1])
    end

    @testset "a regularizer stays on its own component" begin
        x = vcat(x1, x2)
        S = sobel_gradient_healpix(2)
        # Element 4 is a PIXEL INDEX SET, not a flag: `spheroid_regularization` reads
        # `x[reg[4]]` unconditionally.
        regs = Any[Any["sobel2", 1.0, S, 1:n1]]
        g0 = zeros(n1 + n2)
        f0 = binary_crit_allepochs_fg(x, g0, stars1, stars2, data, shifts)
        g2r = zeros(n1 + n2)
        f2r = binary_crit_allepochs_fg(x, g2r, stars1, stars2, data, shifts;
                                       regularizers2 = regs)
        @test f2r > f0
        # The primary's half is untouched; only the secondary's picks up the penalty.
        @test g2r[1:n1] == g0[1:n1]
        @test g2r[n1+1:end] != g0[n1+1:end]
        # ...and the same list on the primary moves the other half instead.
        g1r = zeros(n1 + n2)
        binary_crit_allepochs_fg(x, g1r, stars1, stars2, data, shifts; regularizers1 = regs)
        @test g1r[n1+1:end] == g0[n1+1:end]
        @test g1r[1:n1] != g0[1:n1]
        # The gradient the regularizer adds is the regularizer's own, evaluated on that
        # component's map alone — the indices in element 4 are into `x2`, not into `[x1; x2]`.
        rg = zeros(n2)
        spheroid_regularization(x2, rg; regularizers = regs)
        # A loose tolerance on purpose: the left side is a difference of two O(1)
        # gradients whose regularizer part is O(1e-9), so eight digits are lost to
        # cancellation before the comparison starts.
        @test isapprox(g2r[n1+1:end] .- g0[n1+1:end], rg; rtol = 1e-6)
    end

    @testset "a short reconstruction descends" begin
        x0 = vcat(x1, x2) .* 0.8      # off the truth, so there is something to recover
        g = zeros(n1 + n2)
        f0 = binary_crit_allepochs_fg(x0, g, stars1, stars2, data, shifts)
        xs = binary_reconstruct_oi(x0, data, stars1, stars2, shifts;
                                   maxiter = 15, verbose = false)
        f1 = binary_crit_allepochs_fg(xs, g, stars1, stars2, data, shifts)
        @test f1 < f0
        @test length(xs) == n1 + n2
        @test all(>=(0), xs)          # `lower = 0`: a negative surface brightness is not a map

        a, b = split_binary_map(xs, stars1)
        @test length(a) == n1 && length(b) == n2
        @test vcat(a, b) == xs
        @test split_binary_map(xs, n1)[1] == a
    end

    @testset "Float32 carries through" begin
        # The maps are the unknown and may be single precision; nothing in the gradient may
        # widen them back to Float64. The χ² itself comes back double because the data's
        # errors are, which is OITOOLS' promotion and not this path's business.
        a = Float32.(x1); b = Float32.(x2)
        ga = zeros(Float32, n1); gb = zeros(Float32, n2)
        f32 = binary_chi2_fg(a, ga, stars1[1], b, gb, stars2[1], data[1], shifts[1])
        g64a = zeros(n1); g64b = zeros(n2)
        binary_chi2_fg(x1, g64a, stars1[1], x2, g64b, stars2[1], data[1], shifts[1])
        @test eltype(ga) === Float32 && eltype(gb) === Float32
        @test isapprox(ga, Float32.(g64a); rtol = 1e-4)
        @test isapprox(gb, Float32.(g64b); rtol = 1e-4)
        g32 = zeros(Float32, n1 + n2)
        @test isfinite(binary_crit_allepochs_fg(Float32.(vcat(x1, x2)), g32, stars1, stars2,
                                                data, shifts))
        @test eltype(g32) === Float32
        @test isfinite(f32)
    end

    @testset "the callback and the verbose trace run" begin
        seen = Ref(0)
        x0 = vcat(x1, x2)
        binary_reconstruct_oi(x0, data, stars1, stars2, shifts; maxiter = 3, verbose = false,
                              callback = (x, n, f) -> (seen[] += 1), callback_every = 1)
        @test seen[] > 0
        # A callback that throws is a progress report, not the run: it must be dropped and
        # the reconstruction must finish.
        xs = binary_reconstruct_oi(x0, data, stars1, stars2, shifts; maxiter = 3,
                                   verbose = false, callback = (x, n, f) -> error("boom"),
                                   callback_every = 1)
        @test length(xs) == n1 + n2
        # The per-observable trace has its own arithmetic and is only ever exercised by
        # printing it.
        g1 = zeros(n1); g2 = zeros(n2)
        @test isfinite(binary_chi2_fg(x1, g1, stars1[1], x2, g2, stars2[1], data[1], shifts[1];
                                      verbose = true))
    end
end
