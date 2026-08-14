# Tests for the two radial-binning regularizers in src/oichi2_spheroid.jl:
#   RADFLAT   (`spheroid_radflat_fg`)   — flattens the azimuthally averaged radial profile
#   RADIALVAR (`spheroid_radialvar_fg`) — removes azimuthal scatter within each annulus
#
# They are the BETWEEN-annuli and WITHIN-annuli terms of a one-way ANOVA split of the
# visible disk's variance, and the identity below pins that down rather than leaving it as
# a claim in a docstring. Everything here is analytic or a finite-difference check — no
# recorded reference numbers — so a failure means the regularizer is wrong, not that a
# reference drifted.

using Test, ROTIR, FiniteDifferences, Statistics, Random, LinearAlgebra

const FDM = central_fdm(5, 1)

# A real star with a populated `polyflux`: `radflat_bins` needs the projected areas, which
# only exist after `setup_oi!`. α Cen A is the benchmark dataset (demos/alphacen_ld.jl).
const ACFILE = joinpath(@__DIR__, "..", "demos", "data", "AlphaCenA.oifits")

function rr_star(nside)
    d = readoifits(ACFILE; use_t3 = false, filter_bad_data = true,
                   verbose = false, warn = false)[1, 1]
    p = (surface_type = 0, radius = 4.25, tpole = 5000.0, ldtype = 3, ld1 = 0.14,
         ld2 = 0.0, inclination = 90.0, position_angle = 0.0, rotation_period = 1.0)
    stars = create_star_multiepochs(tessellation_healpix(nside), p, [0.0])
    setup_oi!([d], stars)
    return stars[1], d
end

@testset "radial regularizers" begin
    isfile(ACFILE) || error("test_radial_regularizers: missing $ACFILE")

    @testset "nside $nside" for nside in (3, 4)
        star, _ = rr_star(nside)
        bins = radflat_bins(star; nbins = 6)
        n    = length(bins.idx)
        Random.seed!(20260814 + nside)
        x = 4000.0 .+ 400.0 .* randn(n)

        # --- gradients ------------------------------------------------------------
        # Both regularizers carry a hand-written gradient with a global coupling term from
        # the 1/I_mean normalisation; that term is easy to drop and impossible to notice
        # from the forward value alone, so check it directly.
        for (nm, f!) in (("radflat",   (v, g) -> spheroid_radflat_fg(v, g, bins)),
                         ("radialvar", (v, g) -> spheroid_radialvar_fg(v, g, bins)),
                         ("radialvar unnormalised",
                          (v, g) -> spheroid_radialvar_fg(v, g, bins; normalize = false)))
            g   = similar(x)
            f!(x, g)
            gfd = grad(FDM, v -> f!(v, similar(v)), x)[1]
            # 1e-6 rather than tighter: the unnormalised variant evaluates to ~1e11 at
            # nside 4, so the central difference itself loses relative precision there.
            # The normalised forms come in around 5e-10.
            @test maximum(abs.(g .- gfd)) / maximum(abs.(gfd)) < 1e-6
        end

        # --- the ANOVA identity ---------------------------------------------------
        # Σ(x−x̄)² = Σ_b n_b(x̄_b−x̄)² + Σ_b Σ_{j∈b}(x_j−x̄_b)². The first term is what
        # RADFLAT penalises, the second what RADIALVAR does. Exact, not approximate.
        let m = mean(x), nb = bins.nbins
            S = zeros(nb); cnt = zeros(Int, nb)
            for i in 1:n; S[bins.bin[i]] += x[i]; cnt[bins.bin[i]] += 1; end
            mb      = S ./ cnt
            within  = sum((x[i] - mb[bins.bin[i]])^2 for i in 1:n)
            between = sum(cnt[k] * (mb[k] - m)^2 for k in 1:nb)
            @test isapprox(within + between, sum(abs2, x .- m); rtol = 1e-10)
            # On a random map the azimuthal term dominates by orders of magnitude — which
            # is exactly the freedom RADFLAT leaves untouched.
            @test within > 10 * between
        end

        # --- what each one actually annihilates -----------------------------------
        # A map that is constant within every annulus but steps between annuli: RADIALVAR
        # must see nothing, RADFLAT must see a lot. And vice versa for a map with the same
        # mean in every annulus but scatter inside them.
        let stepped = [1000.0 * bins.bin[i] for i in 1:n]
            g = similar(stepped)
            @test spheroid_radialvar_fg(stepped, g, bins) ≈ 0 atol = 1e-8
            @test all(abs.(g) .< 1e-8)
            @test spheroid_radflat_fg(stepped, similar(stepped), bins) > 1.0
        end
        let balanced = copy(x)                     # force every annulus to a common mean
            for k in 1:bins.nbins
                m = findall(==(k), bins.bin)
                isempty(m) && continue
                balanced[m] .+= 4000.0 - mean(balanced[m])
            end
            @test spheroid_radflat_fg(balanced, similar(balanced), bins) ≈ 0 atol = 1e-6
            @test spheroid_radialvar_fg(balanced, similar(balanced), bins) > 1.0
        end

        # --- both vanish on a constant map, and only there ------------------------
        let flat = fill(4200.0, n)
            @test spheroid_radflat_fg(flat, similar(flat), bins)   ≈ 0 atol = 1e-8
            @test spheroid_radialvar_fg(flat, similar(flat), bins) ≈ 0 atol = 1e-8
        end

        # --- scale invariance -----------------------------------------------------
        # Both are normalised by I_mean, so doubling the map must not change the penalty.
        # `normalize = false` is the OITOOLS convention and must scale as x² instead.
        let g = similar(x)
            @test spheroid_radflat_fg(2x, g, bins)   ≈ spheroid_radflat_fg(x, g, bins)   rtol = 1e-10
            @test spheroid_radialvar_fg(2x, g, bins) ≈ spheroid_radialvar_fg(x, g, bins) rtol = 1e-10
            @test spheroid_radialvar_fg(2x, g, bins; normalize = false) ≈
                  4 * spheroid_radialvar_fg(x, g, bins; normalize = false) rtol = 1e-10
        end

        # --- wiring through the regularizer dispatch ------------------------------
        # The name must be routed, and a typo must still raise rather than contribute zero.
        let xf = fill(4200.0, star.npix), g = zeros(star.npix)
            xf[bins.idx] .= x
            regs = Any[["radflat", 1.0, bins, bins.idx], ["radialvar", 1.0, bins, bins.idx]]
            f = ROTIR.spheroid_regularization(xf, g; regularizers = regs)
            @test f ≈ spheroid_radflat_fg(x, similar(x), bins) +
                      spheroid_radialvar_fg(x, similar(x), bins) rtol = 1e-10
            @test_throws ErrorException ROTIR.spheroid_regularization(
                xf, g; regularizers = Any[["radialvarr", 1.0, bins, bins.idx]])
        end
    end

    # --- orthoLD ------------------------------------------------------------------
    # Rank-1 penalty along the map direction that is exactly degenerate with ld1.
    @testset "orthoLD" begin
        star, _ = rr_star(3)
        p = (surface_type = 0, radius = 4.25, tpole = 5000.0, ldtype = 3, ld1 = 1.6,
             ld2 = 0.0, inclination = 90.0, position_angle = 0.0, rotation_period = 1.0)
        x0 = fill(4000.0, star.npix)
        od = orthold_direction(star, x0, p)
        n  = length(od.idx)
        Random.seed!(11)
        x = od.x0 .+ 300.0 .* randn(n)

        # gradient
        g = similar(x); spheroid_orthold_fg(x, g, od)
        gfd = grad(FDM, v -> spheroid_orthold_fg(v, similar(v), od), x)[1]
        @test maximum(abs.(g .- gfd)) / maximum(abs.(gfd)) < 1e-8

        # rank one: the penalty depends on x ONLY through ⟨x−x₀, ĉ⟩, so any move
        # orthogonal to ĉ is free. This is the whole point — it removes one direction,
        # not a subspace.
        ĉ = od.c ./ sqrt(sum(abs2, od.c))
        let r = randn(MersenneTwister(5), n)
            r .-= dot(r, ĉ) .* ĉ                       # project out the degenerate mode
            @test spheroid_orthold_fg(od.x0 .+ 100 .* r, similar(x), od) ≈ 0 atol = 1e-6
        end
        # ...and a move ALONG it costs, quadratically
        f1 = spheroid_orthold_fg(od.x0 .+ 1.0 .* ĉ, similar(x), od)
        f2 = spheroid_orthold_fg(od.x0 .+ 2.0 .* ĉ, similar(x), od)
        @test f1 > 0 && isapprox(f2/f1, 4.0; rtol = 1e-8)

        # the reference map itself is free
        @test spheroid_orthold_fg(copy(od.x0), similar(x), od) ≈ 0 atol = 1e-10

        # The direction is finite everywhere despite ∂ln(ldmap)/∂ld1 = ln μ diverging at
        # the limb: the model-space weighting cancels it. A pixel-space version would not.
        @test all(isfinite, od.c)
        # Where it peaks is an analytic prediction, and a sharp check on the whole
        # construction: |c| ∝ vis²·μ^{2α}·|ln μ| for the Hestroffer law, maximised at
        # μ = exp(−1/2α). NOT at the limb — μ^{2α} extinguishes the limb faster than ln μ
        # diverges — even though the pixel-space direction x₀·ln μ does diverge there.
        let μ = Float64.(star.normals[od.idx, 3]), α = 1.6
            μpeak = μ[argmax(abs.(od.c))]
            @test isapprox(μpeak, exp(-1/(2α)); atol = 0.12)
            @test 0.5 < μpeak < 0.95
        end

        # dispatch wiring, and a typo must still raise
        let xf = copy(x0), g2 = zeros(length(x0))
            xf[od.idx] .= x
            f = ROTIR.spheroid_regularization(xf, g2;
                    regularizers = Any[["orthold", 1.0, od, od.idx]])
            @test f ≈ spheroid_orthold_fg(x, similar(x), od) rtol = 1e-12
            @test_throws ErrorException ROTIR.spheroid_regularization(
                xf, g2; regularizers = Any[["ortholdd", 1.0, od, od.idx]])
        end
    end

    # Maximum entropy. Pointwise (no neighbour coupling), so unlike `tv2` it imposes no
    # correlation length — but its gradient carries a GLOBAL term, and dropping that for a
    # per-pixel one leaves a gradient almost orthogonal to the truth while still looking
    # plausible. Only a finite-difference check catches it.
    @testset "max entropy" begin
        rng = MersenneTwister(11)
        for n in (20, 200), spread in (0.3, 0.9)
            x = 1.0 .+ spread .* rand(rng, n)
            g = similar(x)
            f = ROTIR.max_entropy_fg(x, g)
            gfd = grad(FDM, y -> ROTIR.max_entropy_fg(y, similar(y)), x)[1]
            @test isfinite(f)
            @test g ≈ gfd rtol = 1e-7 atol = 1e-9
            @test dot(g, gfd)/(norm(g)*norm(gfd)) ≈ 1.0 atol = 1e-10
        end

        # Scale invariance: `xm = x/mean(x)` is untouched by x → c·x, and χ² is invariant
        # too (flux-normalised visibilities), so the weight means the same thing whatever
        # the map level. `tv2` does NOT have this property.
        let x = 1.0 .+ 0.3 .* rand(rng, 50), g = zeros(50)
            f1 = ROTIR.max_entropy_fg(x, g)
            f2 = ROTIR.max_entropy_fg(1000 .* x, similar(x))
            @test f1 ≈ f2 rtol = 1e-12
        end

        # dispatch wiring
        let x = 1.0 .+ 0.3 .* rand(rng, 40), g = zeros(40)
            f = ROTIR.spheroid_regularization(x, g;
                    regularizers = Any[["mem", 1.0, nothing, 1:40]])
            @test f ≈ ROTIR.max_entropy_fg(x, similar(x)) rtol = 1e-12
        end
    end

    # Spherical Sobel gradient. The 2-D Sobel stencil is the inverse-distance²-weighted
    # least-squares gradient on a 3×3 grid; `sobel_gradient_healpix` is that same
    # construction on an irregular sphere. Validated against analytic integrals rather than
    # recorded numbers, so a failure means the operator is wrong.
    @testset "sobel gradient" begin
        # ∫|∇cosθ|²dΩ = ∫sin²θ dΩ = 8π/3, and ∫|∇(3z²−1)/2|²dΩ = l(l+1)∫Y² = 6·4π/5.
        # Both must converge as the mesh refines — that is what makes it a real gradient.
        errs_z = Float64[]; errs_y = Float64[]
        for n in (2, 3, 4)
            S = sobel_gradient_healpix(n; T = Float64)
            u, _ = ROTIR.pix2vec_nest(2^n, collect(1:S.npix); T = Float64)
            @test S.npix == 12*(2^n)^2
            @test S.area ≈ 4π/S.npix

            z = u[:, 3]
            push!(errs_z, abs(S.area*(sum(abs2, S.Gx*z) + sum(abs2, S.Gy*z)) - 8π/3))
            y = (3 .* u[:,3].^2 .- 1) ./ 2
            push!(errs_y, abs(S.area*(sum(abs2, S.Gx*y) + sum(abs2, S.Gy*y)) - 6*4π/5))

            # constants are in the null space exactly: the diagonal carries −Σⱼ
            c = ones(S.npix)
            @test maximum(abs, S.Gx*c) < 1e-12
            @test maximum(abs, S.Gy*c) < 1e-12
        end
        @test errs_z[end] < errs_z[1] && errs_z[end] < 0.01
        @test errs_y[end] < errs_y[1] && errs_y[end] < 0.10

        S = sobel_gradient_healpix(3; T = Float64)
        rng = MersenneTwister(5)
        for (fn, nm) in ((spheroid_sobel2_fg, "sobel2"), (spheroid_sobel_fg, "sobel"))
            x = 4000.0 .+ 400.0 .* randn(rng, S.npix)
            g = similar(x)
            f = fn(x, g, S)
            gfd = grad(FDM, y -> fn(y, similar(y), S), x)[1]
            @test isfinite(f) && f > 0
            @test dot(g, gfd)/(norm(g)*norm(gfd)) ≈ 1.0 atol = 1e-9
            @test g ≈ gfd rtol = 1e-6 atol = 1e-10

            # Scale-invariant, like χ² and `"mem"` — and unlike `"tv2"`.
            @test fn(x, similar(x), S) ≈ fn(1000 .* x, similar(x), S) rtol = 1e-10
        end

        # Resolution stability on a smooth map: the area-weighted Sobel forms converge,
        # `"tv2"` does not (it falls by ~×0.4 per nside doubling), which is why a `"tv2"`
        # weight tuned at one resolution does not transfer to another.
        s2 = Float64[]; t2 = Float64[]
        for n in (3, 4)
            Sn = sobel_gradient_healpix(n; T = Float64)
            u, _ = ROTIR.pix2vec_nest(2^n, collect(1:Sn.npix); T = Float64)
            x = 1.0 .+ 0.3 .* (3 .* u[:,3].^2 .- 1) ./ 2 .+ 0.2 .* u[:,1] .* u[:,2]
            push!(s2, spheroid_sobel2_fg(x, similar(x), Sn))
            push!(t2, ROTIR.spheroid_total_variation2_fg(x, similar(x),
                        tv_neighbors_healpix(n); verbose = false))
        end
        @test isapprox(s2[2], s2[1]; rtol = 0.05)      # converged to a few %
        @test t2[2] < 0.6 * t2[1]                       # tv2 has no such limit

        # dispatch wiring, and a typo must still raise
        let x = 4000.0 .+ 400.0 .* randn(rng, S.npix), g = zeros(S.npix)
            f = ROTIR.spheroid_regularization(x, g;
                    regularizers = Any[["sobel2", 1.0, S, 1:S.npix]])
            @test f ≈ spheroid_sobel2_fg(x, similar(x), S) rtol = 1e-12
            @test_throws ErrorException ROTIR.spheroid_regularization(
                x, g; regularizers = Any[["sobell", 1.0, S, 1:S.npix]])
        end
    end

    # A single-patch annulus has no defined variance; it must contribute zero rather than
    # divide by n−1 = 0. nside 2 is where `radflat_bins` warns about exactly this.
    @testset "degenerate bins" begin
        star, _ = rr_star(2)
        bins = @test_logs (:warn,) match_mode = :any radflat_bins(star; nbins = 6)
        x = 4000.0 .+ 400.0 .* randn(MersenneTwister(3), length(bins.idx))
        g = similar(x)
        f = spheroid_radialvar_fg(x, g, bins)
        @test isfinite(f) && all(isfinite, g)
    end
end
