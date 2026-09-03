# The reconstruction's callable surface: the mid-run callback, per-epoch weights, and the
# joint fitter's θ-layout guard.
#
# All three exist because something was silently wrong. The callback is new; the weights
# warned "not implemented" and ignored their argument; and `joint_reconstruct_oi` fell through
# to the SPHERE layout for any surface type it did not recognise, so a Roche run optimised
# parameters that changed nothing and reported a χ² for it.

using Test
using ROTIR
using ROTIR: spheroid_crit_allepochs_fg

@testset "reconstruction hooks" begin
    D = joinpath(pkgdir(ROTIR), "demos", "data")
    file = joinpath(D, "2011Sep02.lam_And_prepped.oifits")
    data = [readoifits(file; verbose = false)[1, 1]]
    tess = tessellation_healpix(2)
    p = default_star_params(0; radius = 1.35, tpole = 5000.0)
    stars = create_star_multiepochs(tess, p, [0.0])
    setup_oi!(data, stars)
    x0 = Float64.(parametric_temperature_map(p, stars[1]))

    @testset "the mid-run callback" begin
        seen = Tuple{Int,Float64}[]
        keep = Vector{Float64}[]
        x = image_reconstruct_oi(copy(x0), data, stars; maxiter = 40, verbose = false,
                                 callback_every = 5,
                                 callback = (xx, n, f) -> (push!(seen, (n, f));
                                                           push!(keep, copy(xx))))
        @test length(seen) >= 3
        @test first(seen)[1] == 1                       # the first evaluation always reports
        @test issorted([s[1] for s in seen])            # counts advance
        @test all(s -> s[1] == 1 || s[1] % 5 == 0, seen)
        # It is descending overall — the whole point of watching it — though not monotonically,
        # since a rejected line-search step evaluates at a worse point and still counts.
        @test last(seen)[2] < first(seen)[2]
        # The map handed over is the live one, and copying it is the caller's job. Each copy
        # must be a real map of the right length, not a view that has since moved.
        @test all(k -> length(k) == length(x0) && all(isfinite, k), keep)
        @test length(x) == length(x0)
    end

    @testset "a callback that throws does not take the run with it" begin
        n = Ref(0)
        bad = (xx, k, f) -> (n[] += 1; error("deliberate"))
        x = @test_logs (:warn, r"callback failed") match_mode = :any (
            image_reconstruct_oi(copy(x0), data, stars; maxiter = 10, verbose = false,
                                 callback_every = 1, callback = bad))
        @test length(x) == length(x0) && all(isfinite, x)
        # Disabled after the first failure rather than throwing once per evaluation.
        @test n[] == 1
    end

    @testset "per-epoch weights" begin
        # Two epochs of the same data: weighting one to zero must halve the criterion of the
        # pair, which is the property that says the weights reach the sum at all.
        data2 = [data[1], deepcopy(data[1])]
        stars2 = create_star_multiepochs(tess, p, [0.0, 0.0])
        setup_oi!(data2, stars2)
        both = spheroid_chi2_allepochs_f(x0, stars2, data2)
        one  = spheroid_chi2_allepochs_f(x0, stars2, data2; epochs_weights = [1.0, 0.0])
        half = spheroid_chi2_allepochs_f(x0, stars2, data2; epochs_weights = [0.5, 0.5])
        @test both isa Real                              # NOT a vector: `f .* w` after the
        @test one isa Real                               #   sum returned one, and every
        @test half isa Real                              #   caller expects a number
        @test one ≈ both / 2  rtol = 1e-10
        @test half ≈ both / 2 rtol = 1e-10

        # And the GRADIENT carries them, or the map descends a different criterion from the
        # one being reported.
        g_un = zeros(length(x0)); g_w = zeros(length(x0))
        f_un = spheroid_crit_allepochs_fg(x0, g_un, stars2, data2)
        f_w  = spheroid_crit_allepochs_fg(x0, g_w, stars2, data2;
                                          epochs_weights = [0.5, 0.5])
        @test f_w ≈ f_un / 2 rtol = 1e-10
        @test g_w ≈ g_un ./ 2 rtol = 1e-10

        # A wrong-length or negative weight is a mistake about WHICH epochs are weighted, so
        # it is rejected rather than broadcast over the wrong nights.
        @test_throws DimensionMismatch spheroid_chi2_allepochs_f(x0, stars2, data2;
                                                                 epochs_weights = [1.0])
        @test_throws ArgumentError spheroid_chi2_allepochs_f(x0, stars2, data2;
                                                             epochs_weights = [1.0, -1.0])
        @test_throws ArgumentError spheroid_chi2_allepochs_f(x0, stars2, data2;
                                                             epochs_weights = [1.0, NaN])
    end

    @testset "joint_reconstruct_oi refuses a layout it does not have" begin
        # surface_type 3 used to fall through to the SPHERE branch: `merge` added a `radius`
        # field the Roche geometry never reads, so every outer iteration rebuilt an identical
        # star and the θ step moved parameters that did nothing.
        pr = default_star_params(3)
        err = try
            joint_reconstruct_oi(copy(x0), [1.0, 90.0, 0.0], data, tess, pr, [0.0];
                                 nouter = 1, maxiter_xmap = 1, maxiter_θ = 1, verbose = false)
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin("no θ layout for surface_type = 3", err.msg)
        @test occursin("no analytic shape gradient", err.msg)

        # The parameter COUNT is checked too: the sphere layout takes three, and four would
        # otherwise be silently truncated to the first three.
        ps = default_star_params(0; radius = 1.35, tpole = 5000.0)
        @test_throws DimensionMismatch joint_reconstruct_oi(
            copy(x0), [1.0, 90.0, 0.0, 0.0], data, tess, ps, [0.0];
            nouter = 1, maxiter_xmap = 1, maxiter_θ = 1, verbose = false)
    end
end
