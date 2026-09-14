# The ELLIPSOID's log-posterior: NUTS on a triaxial von Zeipel star.
#
# The third of three differentiable forward models, and the one where both halves of the model
# respond to θ — the projected GEOMETRY through the three radii and the orientation, the
# temperature MAP through the three radii, β and (under a non-linear intensity law) `tpole`.
#
#     θ = [rx, ry, rz, inc, PA, β, ld1, ld2]      (+ tpole when `tpole_free`)
#
# Neither pullback asks Zygote to trace the mesh: `project_ellipsoid_geometry` wraps
# `projected_vertices_and_derivs` and `ellipsoid_map` wraps
# `temperature_map_vonZeipel_ellipsoid_derivs`, so what is tested here is that those two
# analytic rrules compose into the right total gradient.
#
# NOTE ON THE BUFFER NAMES BELOW. Anywhere an in-place `!` function is called from a closure
# inside a `@testset`, give its scratch arrays names that appear nowhere else in the testset:
# reusing a name aliases the outer array rather than shadowing it, and the finite-difference
# sweep then overwrites the analytic gradient it is supposed to be compared against. That cost
# an hour once already (see test_ellipsoid_map_derivs.jl).

using Test, ROTIR, Zygote, FiniteDifferences, LinearAlgebra

@testset "the ellipsoid's log-posterior" begin
    DATA = joinpath(pkgdir(ROTIR), "demos", "data")
    fs = sort([joinpath(DATA, f) for f in readdir(DATA) if occursin("lam_And", f)])
    data = [readoifits(f)[1, 1] for f in fs[1:2]]
    tepochs = [0.0, 1.0]
    tess = tessellation_healpix(3; T = Float64)

    mkbase(ldtype) = merge(default_star_params(1),
                           (radius_x = 1.30, radius_y = 1.05, radius_z = 0.92, beta = 0.23,
                            tpole = 6200.0, inclination = 63.0, position_angle = 27.0,
                            ldtype = ldtype, ld1 = 0.25,
                            ld2 = ldtype == 2 ? 0.12 : 0.0))

    @testset "the layout" begin
        @test ellipsoid_param_names() ==
              ["radius_x", "radius_y", "radius_z", "inclination", "position_angle",
               "beta", "ld1", "ld2"]
        @test last(ellipsoid_param_names(tpole_free = true)) == "tpole"
        for tpf in (false, true)
            lb, ub = default_ellipsoid_bounds(tpole_free = tpf)
            @test length(lb) == length(ub) == length(ellipsoid_param_names(tpole_free = tpf))
            @test all(lb .< ub)
            @test all(lb[1:3] .> 0)       # a zero radius is 0/0 in the flux normalisation
        end
        # Names in, positions out — and a name the model does not have is an error, not a
        # silent no-op.
        @test ellipsoid_free_indices(["radius_x", "beta"], 3) == [1, 6]
        @test ellipsoid_free_indices(nothing, 3) == collect(1:8)
        @test_throws ErrorException ellipsoid_free_indices(["rpole"], 3)
        # `ld2` is not read by the power law (ldtype 3), so freeing it would add a direction
        # the posterior is exactly flat along — refused rather than sampled.
        @test_throws ErrorException ellipsoid_free_indices(["ld2"], 3)
        @test ellipsoid_free_indices(["ld2"], 2) == [8]        # the quadratic law does read it
    end

    @testset "the gradient against FiniteDifferences" begin
        fdm = central_fdm(5, 1)
        for (ldtype, model, band, tpf) in ((3, :linear, nothing, false),
                                           (2, :linear, nothing, false),
                                           (2, :planck, 1.65e-6, true))
            base = mkbase(ldtype)
            lp = build_ellipsoid_logπ(data, tess, tepochs, base;
                                      intensity_model = model, band = band,
                                      tpole_free = tpf)
            θ = [1.30, 1.05, 0.92, 63.0, 27.0, 0.23, 0.25, ldtype == 2 ? 0.12 : 0.0]
            tpf && push!(θ, 6200.0)
            gz = Zygote.gradient(lp, θ)[1]
            gf = FiniteDifferences.grad(fdm, lp, θ)[1]
            # The tolerance is the finite-difference floor, not the gradient's: logπ is of
            # order 3e6 here, so a difference of two evaluations loses about nine digits.
            @test norm(gz .- gf) / norm(gf) < 1e-4
            # Every entry the model READS must agree tightly and must not be flat.
            live = ldtype == 2 ? eachindex(gz) : [1, 2, 3, 4, 5, 6, 7]
            @test all(abs.((gz[live] .- gf[live]) ./ gf[live]) .< 1e-6)
            @test !any(iszero, gz[live])
            # And the one it does NOT read is a flat direction — the measurement behind
            # refusing to free it. The strong form is on the VALUE: the power law never looks
            # at `ld2`, so moving it changes logπ by exactly nothing. The analytic gradient is
            # then exactly zero as well.
            #
            # The FINITE DIFFERENCE is the one that cannot be held to zero. It came out at
            # exactly 0 here and at 6.0e-9 under `--check-bounds=yes` (which is what
            # `Pkg.test()` runs), from FiniteDifferences' own adaptive step machinery on a
            # function of order 3e6 that is constant along this direction. Relative to the
            # rest of the gradient that is 2.6e-16, so it is bounded rather than asserted zero.
            if ldtype == 3
                @test lp([θ[1:7]; 0.7]) == lp(θ)
                @test gz[8] == 0
                @test abs(gf[8]) < 1e-12 * norm(gf[1:7])
            end
        end
    end

    @testset "tpole is a scale under :linear and identifiable under :planck" begin
        # The claim the `tpole_free` refusal rests on. The map is `tpole·f_i` with `f`
        # independent of it, every step to the visibility is linear, and `cvis = F/flux`
        # divides it out — so under `:linear` the posterior does not move at all. Under
        # `:planck` the map's CONTRAST changes with `tpole` and it carries information.
        one_ep = data[1:1]
        val(tp, model, band) = begin
            b = merge(mkbase(3), (tpole = tp,))
            lp = build_ellipsoid_logπ(one_ep, tess, [0.0], b;
                                      intensity_model = model, band = band)
            lp([1.30, 1.05, 0.92, 63.0, 27.0, 0.23, 0.25, 0.0])
        end
        a = val(6200.0, :linear, nothing); b = val(3 * 6200.0, :linear, nothing)
        @test abs(b - a) / abs(a) < 1e-12                    # measured at 8.8e-16
        p1 = val(6200.0, :planck, 1.65e-6); p3 = val(3 * 6200.0, :planck, 1.65e-6)
        @test abs(p3 - p1) / abs(p1) > 1e-3                  # measured at 1.1e-2

        # So freeing it under `:linear` is refused, with the reason.
        err = try
            build_ellipsoid_logπ(data, tess, tepochs, mkbase(3);
                                 intensity_model = :linear, tpole_free = true)
            ""
        catch e
            sprint(showerror, e)
        end
        @test occursin("pure scale", err) && occursin("planck", err)
    end

    @testset "the two rrules on their own" begin
        # Each primitive against FiniteDifferences separately, so a composition error and a
        # primitive error cannot be confused.
        base = mkbase(2)
        fdm = central_fdm(5, 1)

        # The MAP: five arguments, and the orientation is not among them because the map does
        # not depend on it.
        mfun(v) = sum(ellipsoid_map(v[1], v[2], v[3], v[4], v[5], tess, base))
        v0 = [1.30, 1.05, 0.92, 0.23, 6200.0]
        gm = Zygote.gradient(mfun, v0)[1]
        fm = FiniteDifferences.grad(fdm, mfun, v0)[1]
        @test norm(gm .- fm) / norm(fm) < 1e-7

        # The GEOMETRY: five arguments, contracted to a scalar so FD has something to chew on.
        #
        # A QUADRATIC contraction, and that matters. The obvious `sum(pw) + 2sum(pn) + 3sum(nz)`
        # is a near-exact INVARIANT of a closed symmetric mesh — the projected coordinates and
        # the normal z-components each sum to about zero — so both gradients come out at their
        # noise floors (measured: 1e-13 from the rrule against 1e-8 from FD) and their ratio is
        # meaningless. Squaring removes the cancellation and leaves a derivative worth checking.
        gfun(v) = begin
            pw, pn, nz = project_ellipsoid_geometry(v[1], v[2], v[3], v[4], v[5],
                                                    tess, 0.0, base)
            sum(abs2, pw) + 2 * sum(abs2, pn) + 3 * sum(abs2, nz)
        end
        w0 = [1.30, 1.05, 0.92, 63.0, 27.0]
        gg = Zygote.gradient(gfun, w0)[1]
        fg = FiniteDifferences.grad(fdm, gfun, w0)[1]
        @test norm(gg .- fg) / norm(fg) < 1e-7
    end
end
