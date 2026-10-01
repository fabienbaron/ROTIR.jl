# The binary forward model's two routes, against each other.
#
# `binary_cvis` combines two components as `(F1 + F2·phase)/(flux1 + flux2)`, which needs each
# star's UNNORMALISED transform and its flux separately. That is why it read `star.polyft` and
# `star.polyflux` directly and, for a long time, was the one forward model with no matrix-free
# route: every binary χ² had to run `setup_oi!` twice and build an `nuv × npix` matrix per
# component. MEASURED on 11 Spica epochs at HEALPix 3: 1964 ms that way, 13.4 ms through
# `:nufft` and 2.7 ms through `:turbo`.
#
# `fused_cvis_parts` exposes those two halves — `fused_cvis` is now just their quotient — and
# `binary_cvis` takes the matrix-free route whenever `polyft` is empty, exactly as
# `observables` does for a single star.
#
# THE DENSE ROUTE IS KEPT, and this file is why it is worth keeping: it is the reference the
# fast one is checked against. It is also what the binary IMAGING criterion runs on, since
# imaging fixes the geometry and reuses one matrix over hundreds of iterations — see
# `binary_chi2_fg` and `test_binary_imaging.jl`, which needs `star.polyft` transposed.

using Test
using ROTIR
using ROTIR: binary_cvis, binary_chi2_f, binary_phase_shift, POLYFT_BACKEND

@testset "binary_cvis: dense vs matrix-free" begin
    D = joinpath(pkgdir(ROTIR), "demos", "data")
    # `readoifits` returns `[wavelength_bin, epoch]` and gives ONE epoch here: the file is
    # not split on MJD gaps unless asked, so this is the whole 2007+2012+2015 campaign as a
    # single table. That is what makes it a good test case — plenty of uv coverage — and it
    # is also why the timings quoted above are per call rather than per night.
    data = [readoifits(joinpath(D, "2007_2012_2015.Spica.oifits"); verbose = false)[1, 1]]
    tess = tessellation_healpix(3)
    # Spica's two components, from demos/spica_params.jl.
    p1 = default_star_params(3; rpole = 0.447, tpole = 25300.0, q = 0.6188)
    p2 = default_star_params(3; rpole = 0.227, tpole = 20585.0, q = 1 / 0.6188)
    te = zeros(length(data))
    x1 = nothing; x2 = nothing

    "Both routes for one epoch, returning `(dense, free)`."
    function both(i, backend)
        s1 = create_star_multiepochs(tess, p1, te; secondary = false)
        s2 = create_star_multiepochs(tess, p2, te; secondary = true)
        m1 = parametric_temperature_map(p1, s1[1]; secondary = false)
        m2 = parametric_temperature_map(p2, s2[1]; secondary = true)
        ph = binary_phase_shift(data[i].uv, 0.6, 0.4)
        free = let old = POLYFT_BACKEND[]
            try
                POLYFT_BACKEND[] = backend
                binary_cvis(m1, s1[i], m2, s2[i], ph; data = data[i])
            finally
                POLYFT_BACKEND[] = old
            end
        end
        d1 = create_star_multiepochs(tess, p1, te; secondary = false)
        d2 = create_star_multiepochs(tess, p2, te; secondary = true)
        setup_oi!(data, d1); setup_oi!(data, d2)
        dense = binary_cvis(m1, d1[i], m2, d2[i], ph)
        return dense, free
    end

    @testset "$(backend), epoch $(i)" for backend in (:nufft, :turbo, :scalar), i in eachindex(data)
        backend === :turbo && !ROTIR.turbo_available() && continue
        dense, free = both(i, backend)
        @test length(free) == length(dense)
        @test all(isfinite, free)
        # `:scalar` and `:turbo` are the SAME arithmetic as the dense matrix, just never
        # materialised, so they agree to rounding. `:nufft` is a quadrature and carries its
        # own error — 6.8e-7 at this mesh, per src/polyft_nfft.jl.
        tol = backend === :nufft ? 1e-5 : 1e-5
        @test maximum(abs.(free .- dense)) / maximum(abs, dense) < tol
    end

    @testset "the χ² agrees too" begin
        # What the GUI's binary panel actually reports, so the number a user reads is pinned
        # and not only the visibilities behind it.
        s1 = create_star_multiepochs(tess, p1, te; secondary = false)
        s2 = create_star_multiepochs(tess, p2, te; secondary = true)
        m1 = parametric_temperature_map(p1, s1[1]; secondary = false)
        m2 = parametric_temperature_map(p2, s2[1]; secondary = true)
        d1 = create_star_multiepochs(tess, p1, te; secondary = false)
        d2 = create_star_multiepochs(tess, p2, te; secondary = true)
        setup_oi!(data, d1); setup_oi!(data, d2)
        for i in eachindex(data)
            ph = binary_phase_shift(data[i].uv, 0.6, 0.4)
            cd_ = binary_chi2_f(m1, d1[i], m2, d2[i], data[i], ph)
            cf  = binary_chi2_f(m1, s1[i], m2, s2[i], data[i], ph)
            @test isfinite(cd_) && cd_ > 0
            @test abs(cf - cd_) / cd_ < 1e-4
        end
    end

    @testset "a binary is not the same as one star" begin
        # The guard against the whole path silently collapsing to the primary: the companion
        # contributes flux and a phase term, so the combined visibilities must differ.
        s1 = create_star_multiepochs(tess, p1, te; secondary = false)
        s2 = create_star_multiepochs(tess, p2, te; secondary = true)
        m1 = parametric_temperature_map(p1, s1[1]; secondary = false)
        m2 = parametric_temperature_map(p2, s2[1]; secondary = true)
        ph = binary_phase_shift(data[1].uv, 0.6, 0.4)
        pair = binary_cvis(m1, s1[1], m2, s2[1], ph; data = data[1])
        solo, _ = fused_cvis_parts(m1, s1[1], data[1])
        solo = solo ./ fused_cvis_parts(m1, s1[1], data[1])[2]
        @test maximum(abs.(pair .- solo)) / maximum(abs, solo) > 0.05
    end
end

# =======================================================================================
# The SCALE of the phase shift, which nothing above constrains.
# =======================================================================================
# Every test above passes an arbitrary offset to `binary_phase_shift` and compares two routes
# that both use it, so any factor in its exponent cancels. It shipped with half the right one
# for that reason: `-π` copied from the polygon kernel, whose argument is a SUM OF TWO
# VERTICES and therefore already twice the midpoint. The companion sat at half its separation,
# and no test — and no |V| — could see it.
#
# These two can. Neither compares ROTIR to ROTIR.

@testset "binary_phase_shift: the exponent" begin
    tess = tessellation_healpix(3)
    p = default_star_params(0; radius = 0.3, tpole = 9000.0)
    # uv spanning CHARA: ±3e8 cycles/rad is ~330 m at 1.1 μm.
    nuv = 64
    uv = permutedims(hcat(collect(range(-3e8, 3e8, length = nuv)),
                          collect(range( 2e8, -2e8, length = nuv))))
    d  = (uv = uv,)
    Δx, Δy = 0.37, -0.52     # mas, deliberately asymmetric so a swapped axis shows

    @testset "shift theorem: transforming a moved mesh == moving the transform" begin
        # The definition of a phase shift, tested on the actual transform rather than on a
        # model of it. `proj_west`/`proj_north` are per-component coordinates with the
        # placement held in `center_offsets`, so translating them IS displacing the star.
        s = create_star(tess, p, 0.0)
        m = fill(9000.0, s.npix)
        F0, fl0 = fused_cvis_parts(m, s, d)
        s2 = deepcopy(s)
        s2.proj_west  .= s.proj_west  .+ Δx
        s2.proj_north .= s.proj_north .+ Δy
        F1, fl1 = fused_cvis_parts(m, s2, d)
        # The shoelace area is translation-invariant, so the flux must not move.
        @test abs(fl1 - fl0) / fl0 < 1e-5
        pred = (F0 ./ fl0) .* binary_phase_shift(uv, Δx, Δy)
        # Tolerance is set by `proj_west` being Float32, not by the phase.
        @test maximum(abs.(F1 ./ fl1 .- pred)) < 1e-5
        # And it is genuinely sensitive to the factor that was wrong: half the exponent has
        # to FAIL, or this test would not have caught the bug it exists for.
        half = (F0 ./ fl0) .* cis.(angle.(binary_phase_shift(uv, Δx, Δy)) ./ 2)
        @test maximum(abs.(F1 ./ fl1 .- half)) > 1e-2
    end

    @testset "two unresolved components reduce to the textbook binary" begin
        # For components far smaller than λ/B, each |V| ≈ 1 and the pair must give exactly
        # `(f1 + f2·exp(-2πi(u·Δα + v·Δδ)))/(f1+f2)` with Δα EAST. ROTIR's offset_x is WEST,
        # so Δα = -Δx. This is the check that fixes the sign of both axes as well as the scale.
        tiny1 = default_star_params(0; radius = 0.002, tpole = 9000.0)
        tiny2 = default_star_params(0; radius = 0.002, tpole = 9000.0)
        s1 = create_star(tess, tiny1, 0.0); s2 = create_star(tess, tiny2, 0.0)
        m1 = fill(9000.0, s1.npix);         m2 = fill(9000.0, s2.npix)
        ph = binary_phase_shift(uv, Δx, Δy)
        got = binary_cvis(m1, s1, m2, s2, ph; data = d)
        mas2rad = pi / (180 * 3600000)
        α = -Δx * mas2rad            # East
        δ =  Δy * mas2rad            # North
        _, f1 = fused_cvis_parts(m1, s1, d)
        _, f2 = fused_cvis_parts(m2, s2, d)
        want = (f1 .+ f2 .* cis.(-2π .* (uv[1, :] .* α .+ uv[2, :] .* δ))) ./ (f1 + f2)
        # 2e-3 is the residual resolution of a 0.002 mas disk at 3e8 cycles/rad, not slack.
        @test maximum(abs.(got .- want)) < 2e-3
    end
end
