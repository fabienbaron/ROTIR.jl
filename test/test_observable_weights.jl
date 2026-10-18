# Observable weights `(w_V², w_T3amp, w_T3φ)` in the criteria that take a surface map or a
# shape: the single-star and binary imaging χ² and the shape-fit χ².
#
# Every one of them is a weighted SUM of three χ² terms, so the value and the gradient must
# both be linear in the weights. Checking that against the three unit weights pins the two
# things that can go wrong: a weight applied to the value but not to the gradient (the map is
# then optimised against a different criterion from the one reported), and a weight on the
# wrong observable.

using Test
using ROTIR
using ROTIR: spheroid_chi2_fg, spheroid_chi2_f, binary_chi2_fg, binary_chi2_f,
             binary_phase_shift, shape_chi2_fg!

const UNIT_W = ((1, 0, 0), (0, 1, 0), (0, 0, 1))

@testset "observable weights" begin
    D = joinpath(pkgdir(ROTIR), "demos", "data")
    data = readoifits(joinpath(D, "2007_2012_2015.Spica.oifits"); verbose = false)[1, 1]
    tess = tessellation_healpix(2)
    p1 = default_star_params(3; rpole = 0.447, tpole = 25300.0, q = 0.6188)
    p2 = default_star_params(3; rpole = 0.227, tpole = 20585.0, q = 1 / 0.6188)
    s1 = create_star(tess, p1, 0.0; secondary = false)
    s2 = create_star(tess, p2, 0.0; secondary = true)
    setup_oi!([data], [s1]); setup_oi!([data], [s2])
    x1 = Float64.(parametric_temperature_map(p1, s1; secondary = false))
    x2 = Float64.(parametric_temperature_map(p2, s2; secondary = true))
    w = (2.0, 0.5, 0.0)

    @testset "single-star imaging" begin
        fg(wt) = (g = zeros(length(x1));
                  f = spheroid_chi2_fg(x1, g, s1, data; verbose = false, weights = wt); (f, g))
        parts = [fg(u) for u in UNIT_W]
        f, g = fg(w)
        @test f ≈ sum(w[k] * parts[k][1] for k in 1:3) rtol = 1e-12
        @test g ≈ sum(w[k] .* parts[k][2] for k in 1:3) rtol = 1e-10
        @test f ≈ spheroid_chi2_f(x1, s1, data; weights = w) rtol = 1e-10
        # The unit weights are the three χ² `chi2_breakdown` reports, in that order.
        b = chi2_breakdown(x1, s1, data)
        @test [p[1] for p in parts] ≈ [b.v2, b.t3amp, b.t3phi] rtol = 1e-10
        # The default is all three, unweighted.
        @test fg((1, 1, 1))[1] ≈ spheroid_chi2_fg(x1, zeros(length(x1)), s1, data;
                                                  verbose = false) rtol = 1e-14
    end

    @testset "binary imaging" begin
        ph = binary_phase_shift(data.uv, 0.6, 0.4)
        fg(wt) = (g1 = zeros(length(x1)); g2 = zeros(length(x2));
                  f = binary_chi2_fg(x1, g1, s1, x2, g2, s2, data, ph; weights = wt);
                  (f, g1, g2))
        parts = [fg(u) for u in UNIT_W]
        f, g1, g2 = fg(w)
        @test f ≈ sum(w[k] * parts[k][1] for k in 1:3) rtol = 1e-10
        @test g1 ≈ sum(w[k] .* parts[k][2] for k in 1:3) rtol = 1e-8
        @test g2 ≈ sum(w[k] .* parts[k][3] for k in 1:3) rtol = 1e-8
        @test f ≈ binary_chi2_f(x1, s1, x2, s2, data, ph; weights = w) rtol = 1e-7
    end

    @testset "shape fit" begin
        base = merge(default_star_params(1), (radius_x = 1.3, radius_y = 1.05, radius_z = 0.92))
        θ = [1.3, 1.05, 0.92, 63.0, 27.0]
        tf = tessellation_healpix(2; T = Float64)
        xm = fill(6000.0, tf.npix)
        fg(wt) = (gθ = zeros(5); gx = zeros(tf.npix);
                  f = shape_chi2_fg!(gθ, gx, copy(xm), θ, [data], tf, base, [0.0];
                                     weights = wt); (f, gθ, gx))
        parts = [fg(u) for u in UNIT_W]
        f, gθ, gx = fg(w)
        @test f ≈ sum(w[k] * parts[k][1] for k in 1:3) rtol = 1e-10
        @test gθ ≈ sum(w[k] .* parts[k][2] for k in 1:3) rtol = 1e-8
        @test gx ≈ sum(w[k] .* parts[k][3] for k in 1:3) rtol = 1e-8
    end
end
