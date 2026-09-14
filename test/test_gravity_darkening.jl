# The two gravity-darkening laws, and the registry that makes them selectable.
#
# von Zeipel is derived for a barotropic star: right for a SLOW rotator, and known to
# overestimate the pole-to-equator contrast for a fast one. Espinosa Lara & Rieutord (2011)
# drop barotropy for "the flux is anti-parallel to effective gravity" and add a latitudinal
# flux factor `F_ω = tan²ϑ/tan²θ`, with ϑ from their implicit eq. (24).
#
# WHAT MAKES THIS TESTABLE WITHOUT TRUSTING THE CODE. The paper carries a closed form for the
# pole-to-equator ratio (eq. 32) that follows from eqs. (27), (28) and (30) and touches
# neither the Newton solve nor the map assembly. Agreement between it and `elr_map` exercises
# the whole chain against arithmetic that shares nothing with it.
#
# AND ONE ASSERTION THAT FINITE DIFFERENCES CANNOT MAKE. Eq. (24) is written for θ ∈ (0, π/2].
# The first implementation let the `θ >= π/2 - ε` equator shortcut catch every SOUTHERN
# colatitude too, so the entire southern hemisphere took the equatorial flux factor — 0.6 %
# too hot at θ = 0.2 rad on the wrong side. Every finite-difference check passed, because the
# analytic derivative made the same substitution self-consistently. `T(θ) == T(π - θ)` is what
# catches it, and it is the reason that test is here.

using Test, ROTIR, Zygote, FiniteDifferences, LinearAlgebra

@testset "gravity darkening: von Zeipel and Espinosa Lara-Rieutord" begin
    fdm = central_fdm(5, 1)

    @testset "the law registry" begin
        @test length(GRAVITY_LAWS) == 2
        @test [l.name for l in GRAVITY_LAWS] == [:vonzeipel, :elr]
        @test [l.code for l in GRAVITY_LAWS] == [1, 2]
        # Every way of naming a law reaches the same spec.
        for x in (:elr, "elr", 2, GRAVITY_LAWS[2])
            @test gravity_law_spec(x).name === :elr
            @test gravity_law_code(x) == 2
            @test gravity_law_name(x) === :elr
        end
        @test gravity_law_choices() == [1 => "von Zeipel", 2 => "Espinosa Lara-Rieutord"]
        # A star_params carries it as a field; one WITHOUT the field is von Zeipel, which is
        # what every model written before the laws became selectable means.
        @test gravity_law_name((surface_type = 2, gravity_law = 2)) === :elr
        @test gravity_law_name((surface_type = 2, rpole = 1.0)) === :vonzeipel
        # A wrong name is an error that lists the alternatives, not a silent fallback to
        # von Zeipel — which would quietly fit the wrong law.
        for bad in (:von_zeipel, "ELR2011", 3)
            e = try; gravity_law_spec(bad); ""; catch err; sprint(showerror, err); end
            @test occursin("implemented", e)
        end
    end

    @testset "the schema field" begin
        @test any(ps -> ps.name === :gravity_law && ps.kind === :choice,
                  surface_params(2))
        # Only the rapid rotator has both laws: ELR for a Roche binary is the 2012 paper and
        # is not implemented, so offering the choice there would be a lie.
        for code in (0, 1, 3)
            @test !any(ps -> ps.name === :gravity_law, surface_params(code))
        end
        # ONE source of truth for the labels. The schema spells the choices out because it is
        # included before the registry is, so this is the assertion that keeps them together —
        # a label edited in one place and not the other is a GUI combo that says one thing and
        # a fit that does another.
        gl = only(filter(ps -> ps.name === :gravity_law, surface_params(2)))
        @test gl.choices == gravity_law_choices()
        @test default_star_params(2).gravity_law == 1                 # von Zeipel by default
        @test default_star_params(2; gravity_law = :elr).gravity_law == 2
        @test default_star_params(2; gravity_law = 2).gravity_law == 2
        # Stored as an Int, like `ldtype` and `surface_type`, so nothing iterating the schema
        # reads it as a continuous coordinate.
        @test default_star_params(2; gravity_law = :elr).gravity_law isa Int
        @test isempty(validate_star_params(default_star_params(2; gravity_law = :elr)))
        # An out-of-table code is reported rather than skipped, which is what the validator
        # used to do for every `:choice` field.
        @test any(m -> occursin("gravity_law", m),
                  validate_star_params(merge(default_star_params(2), (gravity_law = 7,))))
    end

    # ---------------------------------------------------------------------------------------
    # The ELR law against the paper
    # ---------------------------------------------------------------------------------------
    @testset "the ω mapping, against eq. (30)" begin
        # ELR eq. (30) at θ = 0 gives r̃_p = R_p/R_e = 2/(2+ω²). ROTIR's own polar-to-equatorial
        # ratio is 1/f(fev). They must agree, or `elr_omega` maps `frac_escapevel` to the wrong
        # rotation rate and everything downstream is evaluated at the wrong ω.
        for fev in (0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99)
            ω = elr_omega(fev)
            @test 0 < ω < 1
            @test abs(2 / (2 + ω^2) - 1 / ROTIR.f_rapid_rot_and_deriv(fev)[1]) < 1e-14
        end
        @test elr_omega(0.0) == 0.0
    end

    @testset "the equator-to-pole ratio, against eq. (32)" begin
        # The independent check. `elr_temperature_ratio` is the paper's closed form; the map is
        # the Newton solve plus the Roche gravity. Both at β = 1/4, which is where eq. (32)
        # lives.
        for fev in (0.1, 0.3, 0.5, 0.7, 0.9, 0.95)
            x = elr_map(1.0, fev, 0.25, 6000.0, [1.0, 0.0], [0.0, 1.0])
            @test abs(x[1] / x[2] - elr_temperature_ratio(elr_omega(fev))) < 1e-12
        end
    end

    @testset "slow rotation reproduces von Zeipel" begin
        # Eq. (33): as ω → 0 the flux factor goes to 1 and eq. (31) collapses to von Zeipel.
        # A LIMIT, so what is asserted is the RATE — a residual that merely stayed small could
        # still be a constant normalisation error.
        #
        # MEASURED, and the rate is O(ω⁴), not the O(ω²) of the flux factor itself: the map
        # divides by the POLAR flux factor, and the leading term of `F_ω(θ)/F_ω(0)` cancels
        # between numerator and denominator.
        #
        #     fev     max |ΔT/T|
        #     1e-1      2.9e-6
        #     3e-2      2.4e-8     <- 3.33x in fev buys 124x ≈ 3.33⁴
        #     1e-2      2.9e-10
        #     3e-3      3.1e-12
        #     1e-3      6.6e-13    <- the Newton solve's own round-off, not the limit any more
        #     1e-5      6.1e-13
        #
        # Below fev ~ 3e-3 the difference is at the floor where `tan²ϑ/tan²θ` can resolve
        # ϑ − θ at all, so the scaling is asserted over the decades that are resolvable and
        # the floor is asserted as a floor.
        θs = collect(range(0.05, π - 0.05, length = 40))
        sl = sin.(θs); cl = cos.(θs)
        err(fev) = (v = vonzeipel_map(1.0, fev, 0.25, 6000.0, sl, cl);
                    e = elr_map(1.0, fev, 0.25, 6000.0, sl, cl);
                    maximum(abs, (e .- v) ./ v))
        e1, e2, e3 = err(1e-1), err(1e-2), err(1e-3)
        @test e1 < 1e-5
        @test e2 < 1e-8
        @test 3e3 < e1 / e2 < 3e4            # one decade in fev, four in the residual
        @test e3 < 1e-11                     # at the solve's floor
        # And exactly von Zeipel at zero rotation, where `elr_flux_factor` short-circuits.
        @test elr_map(1.0, 0.0, 0.25, 6000.0, sl, cl) ==
              vonzeipel_map(1.0, 0.0, 0.25, 6000.0, sl, cl)
    end

    @testset "the southern hemisphere mirrors the northern" begin
        # THE ASSERTION FINITE DIFFERENCES CANNOT MAKE. A rigidly rotating star is symmetric
        # about its equator, so `F_ω(π - θ) = F_ω(θ)` and the map must be too. Eq. (24) is
        # written for the northern hemisphere only; without folding the colatitude the
        # equatorial shortcut swallows the whole southern half.
        θs = [0.02, 0.2, 0.5, 1.0, 1.3, π/2 - 1e-3, π/2]
        for fev in (0.3, 0.9, 0.95)
            n = elr_map(1.0, fev, 0.25, 6000.0, sin.(θs), cos.(θs))
            m = elr_map(1.0, fev, 0.25, 6000.0, sin.(π .- θs), cos.(π .- θs))
            @test maximum(abs, (m .- n) ./ n) < 1e-13
        end
        # And the flux factor itself, which is where the fold lives.
        for θ in (0.2, 0.7, 1.2), q in (0.1, 0.5, 0.9)
            @test elr_flux_factor(q, θ) ≈ elr_flux_factor(q, π - θ) rtol=1e-14
            F1, d1 = elr_flux_factor_and_dq(q, θ)
            F2, d2 = elr_flux_factor_and_dq(q, π - θ)
            @test F1 ≈ F2 rtol=1e-14
            @test d1 ≈ d2 rtol=1e-14
        end
    end

    @testset "the flux factor and its q-derivative" begin
        # `F_ω` is 1 at zero rotation, rises monotonically with q at fixed θ, and its
        # closed-form limits must join the Newton solve continuously.
        @test elr_flux_factor(0.0, 0.7) == 1.0
        for θ in (0.05, 0.3, 0.8, 1.2, π/2 - 1e-5)
            prev = 0.0
            for q in (0.05, 0.2, 0.4, 0.6, 0.8, 0.95)
                F, dF = elr_flux_factor_and_dq(q, θ)
                @test F > prev; prev = F
                @test dF > 0
                @test F == elr_flux_factor(q, θ)          # one solve, two entry points
                fd = fdm(qq -> elr_flux_factor(qq, θ), q)
                @test abs(dF - fd) / abs(fd) < 1e-6
            end
        end
        # The two limits, against eqs. (27) and (28).
        for q in (0.1, 0.5, 0.9)
            @test elr_flux_factor(q, 0.0) ≈ exp(2q/3) rtol=1e-14
            @test elr_flux_factor(q, π/2) ≈ (1 - q)^(-2/3) rtol=1e-14
        end
    end

    @testset "ELR is warmer at the equator than von Zeipel" begin
        # THE PHYSICS, and the reason the law matters for a fit. von Zeipel overestimates the
        # pole-to-equator contrast; ELR's equator is warmer at the same β, by more as the
        # rotation rises. A fit forced to use von Zeipel has to buy the missing equatorial
        # flux somewhere, and limb darkening is what is for sale.
        prev = 0.0
        for fev in (0.3, 0.6, 0.9, 0.95)
            v = vonzeipel_map(1.0, fev, 0.25, 6000.0, [1.0], [0.0])[1]
            e = elr_map(1.0, fev, 0.25, 6000.0, [1.0], [0.0])[1]
            @test e > v
            gain = (e - v) / v
            @test gain > prev; prev = gain           # monotone in rotation
        end
        @test prev > 0.05                            # 9.8 % at fev = 0.95 as measured
        # The POLE is the normalisation and is identical under both laws by construction.
        for fev in (0.3, 0.9)
            @test vonzeipel_map(1.0, fev, 0.25, 6000.0, [0.0], [1.0])[1] ≈
                  elr_map(1.0, fev, 0.25, 6000.0, [0.0], [1.0])[1] rtol=1e-14
        end
    end

    # ---------------------------------------------------------------------------------------
    # The derivatives
    # ---------------------------------------------------------------------------------------
    @testset "elr_map_and_derivs against FiniteDifferences" begin
        θs = collect(range(0.03, π - 0.03, length = 60))
        s = sin.(θs); c = cos.(θs)
        rp, fev, β, tp = 1.37, 0.85, 0.21, 8200.0
        x, drp, dfev, dβ, dtp = elr_map_and_derivs(rp, fev, β, tp, s, c)
        @test x == elr_map(rp, fev, β, tp, s, c)

        # `∂T/∂rpole` IS EXACTLY ZERO, for this law as for von Zeipel: `g_θ/g_pole` is
        # dimensionless in rpole, and `q = (8/27) fev² f(fev sinθ)³` has no rpole in it at all
        # — `f(fev)` cancels out of `ω² r̃³`. So the criterion is ABSOLUTE. A relative one
        # compares two numerical zeros and its ratio means nothing.
        @test maximum(abs, drp) / maximum(x) < 1e-12
        @test maximum(abs, fdm(r -> elr_map(r, fev, β, tp, s, c), rp)) / maximum(x) < 1e-9

        for (label, ana, f) in (("fev", dfev, v -> elr_map(rp, v, β, tp, s, c)),
                                ("beta", dβ,  v -> elr_map(rp, fev, v, tp, s, c)),
                                ("tpole", dtp, v -> elr_map(rp, fev, β, v, s, c)))
            pt = label == "fev" ? fev : (label == "beta" ? β : tp)
            fd = fdm(f, pt)
            @test norm(ana .- fd) / norm(fd) < 1e-7
            @test !all(iszero, ana)
        end

        # The SIGNS, which separate the right formula from a plausible wrong one: a hotter pole
        # scales the whole map, and a larger β pushes every non-polar tessel further below the
        # pole, since `R = (F_ω/F_p)(g/g_pole) < 1` away from the axis.
        @test all(>(0), dtp)
        @test all(<=(0), dβ)
    end

    @testset "the Zygote primitive" begin
        # `elr_map`'s rrule wraps the analytic derivatives, and both laws' rrules have the same
        # shape, which is what lets `gravity_map` swap one for the other inside a closure.
        θs = collect(range(0.03, π - 0.03, length = 30))
        sg = sin.(θs); cg = cos.(θs)
        w = collect(range(0.5, 1.5, length = length(θs)))   # a non-trivial cotangent
        fun(v) = sum(w .* elr_map(v[1], v[2], v[3], v[4], sg, cg))
        v0 = [1.37, 0.85, 0.21, 8200.0]
        gz = Zygote.gradient(fun, v0)[1]
        gf = FiniteDifferences.grad(fdm, fun, v0)[1]
        # Entry 1 is the exact zero above, so it is checked absolutely and the rest
        # relatively.
        @test abs(gz[1]) < 1e-9 * abs(gf[2])
        @test norm(gz[2:end] .- gf[2:end]) / norm(gf[2:end]) < 1e-7
    end

    # ---------------------------------------------------------------------------------------
    # The two laws through the model
    # ---------------------------------------------------------------------------------------
    @testset "the mesh path agrees with the analytic path" begin
        # `temperature_map_rapid_rotator` takes `q` from the tessel's OWN radius while
        # `elr_map` takes it from the Roche shape factor. On a mesh built by `compute_radii`
        # those are the same number, and this is the measurement that says so — it is what
        # keeps the map the GUI draws identical to the map NUTS samples.
        tess = tessellation_healpix(3; T = Float64)
        colat = tess.unit_spherical[:, 5, 2]
        s = sin.(colat); c = cos.(colat)
        for law in (:vonzeipel, :elr)
            p = default_star_params(2; rpole = 1.2, frac_escapevel = 0.88, beta = 0.22,
                                    tpole = 7400.0, gravity_law = law)
            star = create_star(tess, p, 0.0)
            mesh = temperature_map_rapid_rotator(p, star)
            ana  = gravity_map(Val(law), 1.2, 0.88, 0.22, 7400.0, s, c)
            @test maximum(abs, (mesh .- ana) ./ ana) < 1e-11
            # And the same map the χ² path builds, which is the one that must not drift.
            @test parametric_temperature_map(p, star) == mesh
        end
        # A model with no `gravity_law` field is von Zeipel, unchanged from before the laws
        # became selectable.
        p_old = (surface_type = 2, rpole = 1.2, frac_escapevel = 0.88, beta = 0.22,
                 tpole = 7400.0, ldtype = 3, ld1 = 0.2, ld2 = 0.0, inclination = 60.0,
                 position_angle = 0.0, rotation_period = 1.0, B_rot = 0.0)
        star = create_star(tess, p_old, 0.0)
        @test temperature_map_rapid_rotator(p_old, star) ==
              temperature_map_vonZeipel_rapid_rotator(p_old, star)
    end

    @testset "build_parametric_logπ under both laws" begin
        DATA = joinpath(pkgdir(ROTIR), "demos", "data")
        fs = sort([joinpath(DATA, f) for f in readdir(DATA) if occursin("lam_And", f)])
        data = [readoifits(f)[1, 1] for f in fs[1:2]]
        tepochs = [0.0, 1.0]
        tess = tessellation_healpix(3; T = Float64)
        θ = [1.37, 0.85, 78.0, 24.0, 0.21, 0.23, 0.0]

        base_vz = default_star_params(2; ldtype = 3, tpole = 8200.0, gravity_law = :vonzeipel)
        base_el = merge(base_vz, (gravity_law = 2,))
        lp_vz = build_parametric_logπ(data, tess, tepochs, base_vz)
        lp_el = build_parametric_logπ(data, tess, tepochs, base_el)

        # The law is a property of the MODEL: taken from `base_params` unless overridden, so a
        # caller who set it once need not repeat it at every fit call.
        @test lp_el(θ) == build_parametric_logπ(data, tess, tepochs, base_vz;
                                                gravity_law = :elr)(θ)
        @test lp_vz(θ) == build_parametric_logπ(data, tess, tepochs, base_el;
                                                gravity_law = :vonzeipel)(θ)
        # And it MATTERS at this rotation — a law selector that changed nothing would pass
        # every test above and still be broken.
        @test isfinite(lp_vz(θ)) && isfinite(lp_el(θ))
        @test abs(lp_el(θ) - lp_vz(θ)) / abs(lp_vz(θ)) > 1e-4

        # The gradient under ELR, which is what NUTS integrates. The tolerance is the
        # finite-difference floor: logπ is of order 1e6 here, so differencing two
        # evaluations loses about nine digits.
        ge = Zygote.gradient(lp_el, θ)[1]
        gfd = FiniteDifferences.grad(fdm, lp_el, θ)[1]
        @test norm(ge .- gfd) / norm(gfd) < 1e-4
        # β and `fev` are the two the law changes, and they are the ones that must be right:
        # neither flat nor merely close.
        for j in (2, 5)
            @test abs((ge[j] - gfd[j]) / gfd[j]) < 1e-5
            @test !iszero(ge[j])
        end
        # The two laws give DIFFERENT gradients in β and fev, which is the degeneracy the law
        # choice moves.
        gvz = Zygote.gradient(lp_vz, θ)[1]
        @test abs(ge[5] - gvz[5]) / abs(gvz[5]) > 1e-3
    end
end
