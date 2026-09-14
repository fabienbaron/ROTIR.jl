# The SPHERE's log-posterior: its gradient, and the parameters it is honest about.
#
# `build_sphere_logπ` exists so NUTS can run on a single star. The rapid rotator's
# `build_parametric_logπ` cannot serve: its θ begins `[rpole, frac_escapevel, …]` and it calls
# `vonzeipel_map` directly.
#
# The interesting content of this file is not "the gradient is right" — it is WHICH PARAMETERS
# ARE IN θ, because for a sphere most of them are not identifiable and a sampler handed a flat
# direction wanders over the prior and reports a credible interval that means nothing.

using Test, ROTIR, Zygote, FiniteDifferences, LinearAlgebra

@testset "the sphere's log-posterior" begin
    DATA = joinpath(pkgdir(ROTIR), "demos", "data")
    fs = sort([joinpath(DATA, f) for f in readdir(DATA) if occursin("lam_And", f)])
    data = [readoifits(f)[1, 1] for f in fs[1:2]]
    tepochs = [0.0, 1.0]

    @testset "θ carries only what is identifiable AND differentiable" begin
        # `ld3`/`ld4` are read by Claret's law but `ld_and_derivs` computes no derivative for
        # them — the rrule returns `ZeroTangent()` — so they stay in the forward law and out of
        # θ. Putting them in made them flat directions, and that passes a finite-difference
        # check silently, because FD finds the same zero.
        @test sphere_param_names(1) == ["radius", "ld1"]
        @test sphere_param_names(3) == ["radius", "ld1"]
        @test sphere_param_names(2) == ["radius", "ld1", "ld2"]
        @test sphere_param_names(4) == ["radius", "ld1", "ld2"]   # not ld3, ld4
        for ldtype in (1, 2, 3, 4)
            lb, ub = default_sphere_bounds(ldtype)
            @test length(lb) == length(ub) == length(sphere_param_names(ldtype))
            @test lb[1] > 0            # a zero radius is 0/0 in the flux normalisation
            @test all(lb .< ub)
        end
    end

    @testset "the gradient against FiniteDifferences" begin
        # FD IN FLOAT64. A central difference needs a step well inside the type's precision and
        # Float32 has none to spare, so the meaningful check is at double precision; the Float32
        # mesh is then checked against the Float64 answer instead.
        #
        # THE TOLERANCE IS THE FD FLOOR, NOT THE GRADIENT'S. `logπ` is of order 5e6 here (it is
        # -χ²/2 over two epochs), so a difference of two evaluations loses about nine digits to
        # cancellation, and `∂/∂radius` is ~4e7 — twenty times the limb-darkening components,
        # which is why it dominates the residual.
        #
        # MEASURED, which is what says this is FD and not the gradient: raising the FD order
        # makes it WORSE, because a higher-order rule takes more, finer-spaced samples of the
        # same noisy value. Across the four laws,
        #
        #     central_fdm(5,1):  4.9e-9 … 2.2e-5      <- the best of the three
        #     central_fdm(7,1):  9.9e-3 … 7.1e-2
        #     central_fdm(9,1):  3.4e-3 … 2.9e-1
        #
        # A real gradient error would be order-independent. So the 5-point rule it is, at 1e-4;
        # the LD components, which are not swamped, are held to 1e-5 separately below.
        fdm = central_fdm(5, 1)
        for (ldtype, ld1, ld2) in ((3, 0.25, 0.0), (2, 0.30, 0.15), (1, 0.20, 0.0),
                                   (4, 0.20, 0.10))
            base = merge(default_star_params(0), (ldtype = ldtype, ld1 = ld1, ld2 = ld2))
            nθ = length(sphere_param_names(ldtype))
            θ = [1.05; fill(0.2, nθ - 1)]

            t64 = tessellation_healpix(3; T = Float64)
            l64 = build_sphere_logπ(data, t64, tepochs, base)
            gz = Zygote.gradient(l64, θ)[1]
            gf = FiniteDifferences.grad(fdm, l64, θ)[1]
            @test norm(gz .- gf) / norm(gf) < 1e-4
            # The limb-darkening entries are an order of magnitude smaller than ∂/∂radius and
            # so are not swamped by its cancellation; they get the tight bound.
            @test all(abs.((gz[2:end] .- gf[2:end]) ./ gf[2:end]) .< 1e-5)

            # NO FLAT DIRECTION. This is the assertion that would have caught `ld3`/`ld4`.
            @test !any(iszero, gz)

            # The Float32 mesh, with θ in the mesh's own element type — which is what
            # `_fit_hmc` converts it to.
            t32 = tessellation_healpix(3; T = Float32)
            l32 = build_sphere_logπ(data, t32, tepochs, base)
            g32 = Zygote.gradient(l32, Float32.(θ))[1]
            @test norm(Float64.(g32) .- gz) / norm(gz) < 5e-2
        end
    end

    @testset "what a sphere cannot tell you" begin
        # The reason θ is this short, as a measurement rather than an assertion of principle.
        base = default_star_params(0)
        tess = tessellation_healpix(3)
        function v2_of(; kw...)
            p = merge(base, kw)
            star = create_star(tess, p, 0.0)
            return observables(parametric_temperature_map(p, star), star, data[1])[1]
        end
        ref = v2_of()
        rel(x) = maximum(abs, x .- ref) / maximum(abs, ref)

        # A uniform limb-darkened sphere looks the same from every direction, so its
        # ORIENTATION carries no information: what little there is, is HEALPix faceting.
        d_inc = rel(v2_of(inclination = 5.0))
        d_pa  = rel(v2_of(position_angle = 137.0))
        # And what a sphere DOES measure, an order of magnitude above that floor.
        d_rad = rel(v2_of(radius = base.radius * 1.02))
        d_ld  = rel(v2_of(ld1 = base.ld1 + 0.15))
        @test d_inc < 5e-3 && d_pa < 5e-3
        @test d_rad > 10 * max(d_inc, d_pa)
        @test d_ld  > 10 * max(d_inc, d_pa)
        # The 60 -> 120 degree reflection maps the mesh onto itself, so even the faceting goes.
        @test rel(v2_of(inclination = 180.0 - base.inclination)) < 1e-5

        # `tpole` is a pure multiplicative scale on a uniform map and `cvis = F/flux` divides it
        # out, which is why it is not in θ either.
        @test rel(v2_of(tpole = base.tpole * 3)) < 1e-6
    end
end
