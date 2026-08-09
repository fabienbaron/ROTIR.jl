#!/usr/bin/env julia
# Tests for the common binary frame (src/binary_geometry.jl) and the mutual-irradiation
# radiosity solver (src/reflection.jl).
#
# Everything here is analytic or a limiting case — nothing is a regression against a
# previously-recorded number, so a failure means the physics is wrong, not that a
# reference file drifted.
#
# Runs standalone (`julia --project=. test/test_reflection.jl`) or under Pkg.test().

using ROTIR, Test, LinearAlgebra, Printf

const P_ORB = 4.0145
const A_ORB = 1.54
const E_ORB = 0.123
const T0_ORB = 2454189.40
const Q_BIN = 0.6188
const I_ORB = 116.0
const OM_NODE = 309.938
const OM_PERI = 255.0

const R1 = 0.4469        # mas, angular RADIUS (7.40 R☉ at 77 pc)
const R2 = 0.2259        # mas, angular RADIUS (3.74 R☉ at 77 pc)
const T1 = 25300.0
const T2 = 20585.0

sphere_params(r, t) = (surface_type = 0, radius = r, tpole = t, ldtype = 1,
                       ld1 = 0.0, ld2 = 0.0, inclination = 180.0 - I_ORB,
                       position_angle = OM_NODE - 180.0, rotation_period = P_ORB)

const SP1 = sphere_params(R1, T1)
const SP2 = sphere_params(R2, T2)

roche_params(r, t, q) = (surface_type = 3, rpole = r, tpole = t, q = q,
                         ldtype = 3, ld1 = 0.15, ld2 = 0.0, beta = 0.25, d = 77.0,
                         inclination = 180.0 - I_ORB, position_angle = OM_NODE - 180.0,
                         rotation_period = P_ORB,
                         fillout_factor_primary = -1, fillout_factor_secondary = -1,
                         i = I_ORB, Ω = OM_NODE, ω = OM_PERI, P = P_ORB, a = A_ORB,
                         e = E_ORB, T0 = T0_ORB, dP = 0.0, dω = 0.0)

const RP1 = roche_params(R1, T1, Q_BIN)
const RP2 = roche_params(R2, T2, 1.0 / Q_BIN)

const S1P = starparameters(R1, T1, 0.0, 1, 0.0, 0.0, 0.25, 0.0,
                           180.0 - I_ORB, OM_NODE - 180.0, 0.0, P_ORB)
const S2P = starparameters(R2, T2, 0.0, 1, 0.0, 0.0, 0.25, 0.0,
                           180.0 - I_ORB, OM_NODE - 180.0, 0.0, P_ORB)
const BP  = binaryparameters(S1P, S2P, 77.0, I_ORB, OM_NODE, OM_PERI, P_ORB, A_ORB,
                             E_ORB, T0_ORB, Q_BIN, [1.0, 1.0], 0.0, 0.0)

const TEST_T = T0_ORB + 0.31 * P_ORB     # a generic, non-special phase

@testset "binary frame" begin
    R, D, rrel = binary_frame(BP, TEST_T)
    @test size(R) == (3, 3)
    @test norm(R' * R - I) < 1e-12                  # orthonormal
    @test det(R) ≈ 1.0 atol = 1e-12                 # right-handed, no reflection
    @test norm(R * [D, 0.0, 0.0] - rrel) < 1e-10    # x̂ is the line of centres
    @test norm(rrel) ≈ D atol = 1e-10

    # The frame must agree with the offset the rest of ROTIR already uses.
    ow, on = orbit_to_rotir_offset(BP, TEST_T)
    @test rrel[1] ≈ ow atol = 1e-10
    @test rrel[2] ≈ on atol = 1e-10

    # ẑ is the orbital angular momentum: constant in time, ⊥ to the line of centres.
    for t in (T0_ORB, T0_ORB + 0.17P_ORB, T0_ORB + 0.62P_ORB)
        Rt, _, _ = binary_frame(BP, t)
        @test abs(dot(Rt[:, 3], R[:, 3]) - 1) < 1e-10
        @test abs(dot(Rt[:, 1], Rt[:, 3])) < 1e-12
    end

    # Separation matches the independent Kepler route.
    @test D ≈ A_ORB * compute_separation(BP, TEST_T) atol = 1e-8
end

@testset "tidal bulge tracks the companion" begin
    # The Roche body-frame +x̂ (θ=π/2, φ=0) is the sub-companion direction for the primary.
    # After create_binary_star it must project onto the sky along the secondary's offset.
    tes = tessellation_healpix(4)
    for t in (T0_ORB, T0_ORB + 0.25P_ORB, T0_ORB + 0.5P_ORB, T0_ORB + 0.8P_ORB)
        star = create_binary_star(tes, RP1, BP, t)
        # tessel whose body-frame direction is closest to (θ=π/2, φ=0)
        θ = star.vertices_spherical[:, 5, 2]
        φ = star.vertices_spherical[:, 5, 3]
        λ = sin.(θ) .* cos.(φ)                       # direction cosine toward companion
        k = argmax(λ)
        @test λ[k] > 0.99                            # nside 4 resolves the bulge tip
        v = Float64.(star.vertices_xyz[k, 5, :])
        _, _, rrel = binary_frame(BP, t)
        # ORIENTATION: the sub-companion tessel projects along the companion's offset.
        # This part is purely about the frame built in binary_geometry.jl and passes at
        # every phase.
        @test dot(v, rrel) / (norm(v) * norm(rrel)) > 0.99

        # SHAPE: the sub-companion point must also be the LONGEST radius — that is what
        # "tidal bulge" means — at EVERY phase, including apastron.
        #
        # Regression guard for the tidal term of the Roche potential. The frame is centred
        # on the star, which is in free fall toward the companion, so the pseudo-force
        # contributes −q·r·λ/D². Writing it instead as −q·r·λ·D (Aufdenberg+2015 A18)
        # over-weights it by D³ and, past D/a ≈ 1.031 for Spica's q and rpole/a, inverts
        # the bulge — it pointed AWAY from the companion for roughly half of every orbit.
        r = star.vertices_spherical[:, 5, 1]
        @test r[k] ≈ maximum(r) rtol = 1e-3
    end
end

@testset "Roche potential: gravity is the gradient of the potential" begin
    # compute_gravity_* must be exactly |∇Ω| of compute_potential_*. This is what catches
    # a change to one without the matching change to the other — the centre-of-mass term
    # appears in both (as −q·r·λ/D² and as the constant −q/D²).
    q = 0.6188; F = 1.0
    function fd_grad_norm(pf, r, θ, ϕ, D; h = 1e-7)
        λ = cos(ϕ)sin(θ); μ = sin(ϕ)sin(θ); ν = cos(θ)
        p = [r*λ, r*μ, r*ν]
        Ω(v) = begin
            rr = sqrt(sum(abs2, v))
            pf(rr, D, acos(v[3]/rr), atan(v[2], v[1]), q, F)[1]
        end
        g = [(Ω(p .+ h .* (1:3 .== i)) - Ω(p .- h .* (1:3 .== i))) / (2h) for i in 1:3]
        return sqrt(sum(abs2, g))
    end
    for D in (0.85, 1.0, 1.15), (θ, ϕ) in ((0.7, 0.0), (1.4, 0.3), (2.2, 2.0))
        @test compute_gravity_primary(0.30, θ, ϕ, D, q, F) ≈
              fd_grad_norm(ROTIR.compute_potential_primary, 0.30, θ, ϕ, D) rtol = 1e-6
        @test compute_gravity_secondary(0.15, θ, ϕ, D, q, F) ≈
              fd_grad_norm(ROTIR.compute_potential_secondary, 0.15, θ, ϕ, D) rtol = 1e-6
        # and the analytic dΩ/dr must match FD too
        for pf in (ROTIR.compute_potential_primary, ROTIR.compute_potential_secondary)
            r0 = 0.2
            fd = (pf(r0 + 1e-7, D, θ, ϕ, q, F)[1] - pf(r0 - 1e-7, D, θ, ϕ, q, F)[1]) / 2e-7
            @test pf(r0, D, θ, ϕ, q, F)[2] ≈ fd rtol = 1e-6
        end
    end
end

@testset "Roche potential reduces to the D ≡ 1 reference" begin
    # /home/baron/SOFTWARE/roche/RocheLobe.f90 `potential` works in units of the
    # instantaneous separation, so D ≡ 1 there. Rescaling our potential (r = D·s, then
    # multiply by D) must reproduce it exactly — that is the check that pins the
    # D-dependence of the centre-of-mass term.
    q = 0.6188; F = 1.0
    reference(s, θ, ϕ) = 1/s + q/sqrt(1 - 2s*cos(ϕ)sin(θ) + s^2) - q*s*cos(ϕ)sin(θ)
    for D in (0.85, 1.0, 1.15), s in (0.2, 0.3, 0.4), (θ, ϕ) in ((pi/2, 0.0), (pi/2, pi), (1.1, 0.7))
        r = D*s
        # strip the rotational term (Wilson's F convention deliberately carries no D)
        Ω = ROTIR.compute_potential_primary(r, D, θ, ϕ, q, F)[1] -
            F^2*(1+q)*r^2*(1 - cos(θ)^2)/2
        @test D*Ω ≈ reference(s, θ, ϕ) rtol = 1e-12
    end
end

@testset "3-D tessel areas" begin
    tes = tessellation_healpix(4)
    star = create_binary_star(tes, SP1, BP, TEST_T)
    c, A = tessel_centroids_areas(star)
    @test size(c) == (tes.npix, 3)
    @test length(A) == tes.npix
    @test all(A .> 0)
    @test sum(A) ≈ 4pi * R1^2 rtol = 0.01           # flat quads under-cover by <1%
    # centroids sit on the sphere
    @test all(abs.(sqrt.(sum(abs2, c, dims=2)) .- R1) .< 0.02 * R1)

    # `polyflux` is the PROJECTED (2-D shoelace) area used by the Fourier transform, NOT
    # this one — reusing it for view factors would be wrong. Over the visible hemisphere
    # the projected area sums to the disk πR², roughly a quarter of the 3-D total 4πR².
    vis = star.index_quads_visible
    pf = setup_polyflux_single(star.proj_west[vis, :], star.proj_north[vis, :])
    @test sum(abs, pf) ≈ pi * R1^2 rtol = 0.05
    @test sum(abs, pf) < 0.35 * sum(A)
end

@testset "reflection: limiting cases" begin
    tes1 = tessellation_healpix(3)
    tes2 = tessellation_healpix(2)
    s1, s2 = create_binary_geometry(tes1, SP1, tes2, SP2, BP, TEST_T)
    t1 = parametric_temperature_map(SP1, s1)
    t2 = parametric_temperature_map(SP2, s2; secondary = true)
    @test all(t1 .≈ Float32(T1))          # sphere ⇒ uniform
    @test all(t2 .≈ Float32(T2))

    # albedo 0 ⇒ nothing happens at all
    h1, h2 = handle_reflection(s1, t1, s2, t2; albedo1 = 0.0, albedo2 = 0.0)
    @test maximum(abs.(h1 .- t1)) == 0
    @test maximum(abs.(h2 .- t2)) == 0

    # far apart ⇒ no heating
    bp_far = merge(BP, (a = 1000 * A_ORB,))
    f1, f2 = create_binary_geometry(tes1, SP1, tes2, SP2, bp_far, TEST_T)
    g1, g2 = handle_reflection(f1, t1, f2, t2; albedo1 = 1.0, albedo2 = 1.0)
    @test maximum(abs.(g1 .- t1)) < 1e-2
    @test maximum(abs.(g2 .- t2)) < 1e-2

    # heating is one-sided: the far hemisphere of each star is untouched
    hh1, hh2 = handle_reflection(s1, t1, s2, t2; albedo1 = 0.6, albedo2 = 0.6)
    @test maximum(hh2) > maximum(t2)
    @test count(hh2 .≈ t2) > 0.3 * length(t2)     # a substantial night side
    @test all(hh2 .>= t2 .- 1e-3)                 # irradiation never cools

    # refusing to run on stars that were not placed in a common frame
    lone1 = create_star(tes1, SP1, 0.0)
    lone2 = create_star(tes2, SP2, 0.0)
    @test_throws ErrorException handle_reflection(lone1, t1, lone2, t2)
end

@testset "Horvat ≡ Wilson for uniform bolometric LD" begin
    # Algebraic identity: with D ≡ 1 the limb-darkened and Lambert kernels coincide, so
    # Fout (Horvat) and M (Wilson) satisfy the same linear system. Any difference is a
    # bug in one of the two solvers.
    tes1 = tessellation_healpix(3); tes2 = tessellation_healpix(2)
    s1, s2 = create_binary_geometry(tes1, SP1, tes2, SP2, BP, TEST_T)
    t1 = parametric_temperature_map(SP1, s1)
    t2 = parametric_temperature_map(SP2, s2; secondary = true)
    for A in (0.3, 0.6, 1.0)
        a1, a2 = handle_reflection(s1, t1, s2, t2; albedo1 = A, albedo2 = A, method = :horvat)
        b1, b2 = handle_reflection(s1, t1, s2, t2; albedo1 = A, albedo2 = A, method = :wilson)
        @test maximum(abs.(a1 .- b1)) < 1e-4
        @test maximum(abs.(a2 .- b2)) < 1e-4
    end
    # ...and they must genuinely differ once the bolometric LD is not uniform.
    ld = (ldtype = 1, ld1 = 0.6, ld2 = 0.0)
    c1, _ = handle_reflection(s1, t1, s2, t2; albedo1 = 0.9, albedo2 = 0.9,
                              ldbol1 = ld, ldbol2 = ld, method = :horvat)
    d1, _ = handle_reflection(s1, t1, s2, t2; albedo1 = 0.9, albedo2 = 0.9,
                              ldbol1 = ld, ldbol2 = ld, method = :wilson)
    @test maximum(abs.(c1 .- d1)) > 1e-3
end

@testset "bolometric LD normalisation" begin
    # ld_bol_D0 returns D0/π = 2∫₀¹D(μ)μdμ; check against numerical quadrature.
    μ = range(1e-6, 1, length = 20001)
    num(ldt, x, y) = 2 * sum(ld_bol(ldt, x, y, m) * m for m in μ) * step(μ)
    for (ldt, x, y) in ((0, 0.0, 0.0), (1, 0.35, 0.0), (2, 0.35, 0.15), (3, 0.15, 0.0),
                        (3, 0.8, 0.0), (1, 0.9, 0.0))
        @test ld_bol_D0(ldt, x, y) ≈ num(ldt, x, y) rtol = 1e-4
    end
    @test ld_bol_D0(0, 0.0, 0.0) == 1.0                # uniform ⇒ F ≡ F0

    # A normalised law must not change the TOTAL emitted flux, only its angular
    # distribution — so the total heating delivered is nearly LD-independent.
    tes1 = tessellation_healpix(3); tes2 = tessellation_healpix(2)
    s1, s2 = create_binary_geometry(tes1, SP1, tes2, SP2, BP, TEST_T)
    t1 = parametric_temperature_map(SP1, s1)
    t2 = parametric_temperature_map(SP2, s2; secondary = true)
    _, A2 = tessel_centroids_areas(s2)
    tot(ld) = begin
        _, h2 = handle_reflection(s1, t1, s2, t2; albedo1 = 1.0, albedo2 = 1.0,
                                  ldbol1 = ld, ldbol2 = ld)
        sum(Float64.(A2) .* (Float64.(h2) .^ 4 .- Float64.(t2) .^ 4))
    end
    base = tot((ldtype = 0,))
    @test isapprox(tot((ldtype = 1, ld1 = 0.4, ld2 = 0.0)), base, rtol = 0.05)
    @test isapprox(tot((ldtype = 3, ld1 = 0.3, ld2 = 0.0)), base, rtol = 0.05)
end

@testset "reflection: analytic substellar point" begin
    # Two spheres: the flux at a point on the secondary's surface facing the primary is
    # exactly σT₁⁴(R₁/d)² with d = D − R₂, so
    #     T_new = T₂ (1 + A (R₁/d)² (T₁/T₂)⁴)^(1/4).
    # The numerical solver must approach this, and the error must SHRINK with resolution
    # (a fixed offset would mean a systematic error, not discretisation).
    _, D, _ = binary_frame(BP, TEST_T)
    A = 0.6
    pred = T2 * (1 + A * (R1 / (D - R2))^2 * (T1 / T2)^4)^0.25
    errs = Float64[]
    for (n1, n2) in ((3, 2), (4, 3), (5, 4))
        tes1 = tessellation_healpix(n1); tes2 = tessellation_healpix(n2)
        s1, s2 = create_binary_geometry(tes1, SP1, tes2, SP2, BP, TEST_T)
        t1 = parametric_temperature_map(SP1, s1)
        t2 = parametric_temperature_map(SP2, s2; secondary = true)
        _, h2 = handle_reflection(s1, t1, s2, t2; albedo1 = A, albedo2 = A)
        push!(errs, abs(maximum(h2) - pred) / (pred - T2))
    end
    @test errs[1] < 0.10                 # already close at 768/192 tessels
    @test errs[end] < 0.03               # and converging
    @test errs[end] <= errs[1]           # monotone improvement with resolution
end

@testset "energy conservation" begin
    tes1 = tessellation_healpix(3); tes2 = tessellation_healpix(2)
    s1, s2 = create_binary_geometry(tes1, SP1, tes2, SP2, BP, TEST_T)
    t1 = parametric_temperature_map(SP1, s1)
    t2 = parametric_temperature_map(SP2, s2; secondary = true)
    G, L12, L21, A1, A2 = reflection_kernels(s1, s2)
    F0_1 = ROTIR.SIGMA_SB .* Float64.(t1) .^ 4
    F0_2 = ROTIR.SIGMA_SB .* Float64.(t2) .^ 4
    A = 0.6
    o1, o2, in1, in2, nit = solve_radiosity(G, L12, L21, A1, A2, A, A, F0_1, F0_2)
    @test 1 <= nit < 100                                    # converged, not maxed out
    # The extra emitted power is exactly the albedo times the intercepted power.
    @test sum(A1 .* (o1 .- F0_1)) ≈ A * sum(A1 .* in1) rtol = 1e-10
    @test sum(A2 .* (o2 .- F0_2)) ≈ A * sum(A2 .* in2) rtol = 1e-10

    # Total power intercepted by star 2 from star 1, at first order (albedo 0 ⇒ no
    # multiple scattering), equals L₁ × the solid-angle fraction Ω₂/4π seen from star 1.
    p1, p2, _, _, _ = solve_radiosity(G, L12, L21, A1, A2, 0.0, 0.0, F0_1, F0_2)
    @test p1 == F0_1 && p2 == F0_2                          # albedo 0 ⇒ F_out = F_0
    S0_2 = L21' * (A1 .* F0_1)
    L1_tot = sum(A1 .* F0_1)
    _, D, _ = binary_frame(BP, TEST_T)
    frac = (R2 / (2 * D))^2                                  # πR₂²/(4πD²)
    @test sum(A2 .* S0_2) / L1_tot ≈ frac rtol = 0.05
end

@testset "volume-conserving Roche shape" begin
    tes = tessellation_healpix(3)
    Dref, Dper, Dapo = 1.0, 1 - E_ORB, 1 + E_ORB
    Ωref, _ = ROTIR.compute_potential_primary(RP1.rpole / RP1.a, Dref, 0.0, 0.0,
                                              RP1.q, RP1.rotation_period / RP1.P)
    Vref = roche_mesh_volume(tes, RP1, Dref, Ωref)
    for D in (Dper, Dapo)
        Ω = roche_omega_for_volume(tes, RP1, D)
        @test roche_mesh_volume(tes, RP1, D, Ω) ≈ Vref rtol = 1e-4
    end
    # ...whereas the fixed-rpole convention does NOT conserve volume (that is the point)
    Ωfix, _ = ROTIR.compute_potential_primary(RP1.rpole / RP1.a, Dper, 0.0, 0.0,
                                              RP1.q, RP1.rotation_period / RP1.P)
    @test !isapprox(roche_mesh_volume(tes, RP1, Dper, Ωfix), Vref; rtol = 1e-3)

    # Solid angles must tile the sphere.
    @test sum(tessel_solid_angles(tes)) ≈ 4pi rtol = 1e-12
end

@testset "roche_volume / roche_area against the mesh and the sphere limit" begin
    # Sphere limit: q → 0, no rotation, companion far away ⇒ Ω → 1/r exactly.
    # (q·D must also stay negligible, since the centre-of-mass term carries a factor D.)
    for Ω in (2.0, 4.0, 10.0)
        r = 1 / Ω
        @test roche_volume(1e-10, 0.0, 1e4; omega = Ω) ≈ 4pi / 3 * r^3 rtol = 1e-6
        @test roche_area(1e-10, 0.0, 1e4; omega = Ω)   ≈ 4pi * r^2    rtol = 1e-6
    end
    # Distorted lobe: the quadrature must match a direct mesh sum.
    tes = tessellation_healpix(5)
    async = RP1.rotation_period / RP1.P
    for D in (1 - E_ORB, 1.0, 1 + E_ORB)
        Ω, _ = ROTIR.compute_potential_primary(RP1.rpole / RP1.a, D, 0.0, 0.0, RP1.q, async)
        @test roche_volume(RP1.q, async, D; omega = Ω) ≈
              roche_mesh_volume(tes, RP1, D, Ω) rtol = 1e-3
        Rv, Ra = roche_equivalent_radius(RP1.q, async, D; omega = Ω)
        @test Rv > RP1.rpole / RP1.a          # equivalent radius exceeds the polar radius
        @test Rv ≈ Ra rtol = 5e-3             # consistent for a mildly distorted lobe
    end
end

@testset "Planck vs linear intensity" begin
    T = [T1, T2]
    for (λ, expected) in ((1.65e-6, 1.282), (0.65e-6, 1.381), (0.55e-6, 1.415))
        I = intensity(T, :planck, λ)
        @test I[1] / I[2] ≈ expected rtol = 2e-3
    end
    @test intensity(T, :linear, nothing) == T
    @test_throws ErrorException intensity(T, :planck, nothing)
    @test_throws ErrorException intensity(T, :bogus, 1.65e-6)
    # Rayleigh–Jeans: at long wavelength the Planck ratio tends to the linear one.
    I = intensity(T, :planck, 1.0e-2)
    @test I[1] / I[2] ≈ T1 / T2 rtol = 1e-3
end

@testset "occultation diagnostic" begin
    tes1 = tessellation_healpix(3); tes2 = tessellation_healpix(2)
    ts = T0_ORB .+ collect(range(0, P_ORB, length = 60))
    s1, s2 = create_binary_geometry_multiepochs(tes1, SP1, tes2, SP2, BP, ts)
    flags = check_binary_overlap(s1, s2, BP, ts; verbose = false)
    @test length(flags) == length(ts)
    # Spica is a known grazing eclipser (Desmet+2009), so some epochs must flag...
    @test any(flags)
    # ...but most must not.
    @test count(flags) < length(ts) ÷ 3
    @test all(projected_separation.(Ref(BP), ts[flags]) .<
              projected_separation.(Ref(BP), ts[.!flags])[1] + 1)
    # Far apart ⇒ never flagged.
    bp_far = merge(BP, (a = 100 * A_ORB,))
    f1, f2 = create_binary_geometry_multiepochs(tes1, SP1, tes2, SP2, bp_far, ts[1:5])
    @test !any(check_binary_overlap(f1, f2, bp_far, ts[1:5]; verbose = false))
end
