# The ellipsoid's von Zeipel map, differentiated.
#
# This is what a CONSISTENT gradient fit for surface type 1 was missing. `shape_chi2_fg!` has
# analytic derivatives of the projected geometry and the face normals for surface types 0, 1
# and 2, but it holds the TEMPERATURE MAP fixed — which is exact only where the map does not
# depend on the parameters being fitted. For an ellipsoid it does, through `radius_x/y/z`,
# `beta` and `tpole`.
#
# The map is
#
#     g_i = 1/r_i²,  g_pole = 1/rx²,  T_i = tpole·(g_i/g_pole)^β = tpole·(rx/r_i)^{2β}
#
# with `r_i = |(rx·u_x, ry·u_y, rz·u_z)|` for the tessel's body-frame unit direction `u`.
#
# The second result here matters as much as the first: `r_i` is a NORM and rotation is
# orthogonal, so the map does not depend on the ORIENTATION at all. That is what makes a
# shape-only gradient exact for an ellipsoid whose free parameters are just the two angles.

using Test, ROTIR, FiniteDifferences, LinearAlgebra

@testset "the ellipsoid's temperature-map derivatives" begin
    tess = tessellation_healpix(3; T = Float64)
    base = merge(default_star_params(1),
                 (radius_x = 1.30, radius_y = 1.05, radius_z = 0.92,
                  beta = 0.23, tpole = 6200.0,
                  inclination = 63.0, position_angle = 27.0))
    map_of(; kw...) = (p = merge(base, kw);
                       temperature_map_vonZeipel_ellipsoid(p, create_star(tess, p, 0.0)))
    star = create_star(tess, base, 0.0)
    Tm, drx, dry, drz, dβ, dtp = temperature_map_vonZeipel_ellipsoid_derivs(base, star)

    # The derivative routine must return the SAME map the plain one does, or the derivatives
    # belong to a different function than the one being fitted.
    @test maximum(abs, Tm .- map_of()) / maximum(abs, Tm) < 1e-12
    @test length(Tm) == star.npix
    @test all(>(0), Tm)

    fdm = central_fdm(5, 1)
    dmap(key) = FiniteDifferences.jacobian(
        fdm, v -> map_of(; (key => v[1],)...), [Float64(getproperty(base, key))])[1][:, 1]

    for (key, ana) in ((:radius_x, drx), (:radius_y, dry), (:radius_z, drz),
                       (:beta, dβ), (:tpole, dtp))
        fd = dmap(key)
        @test norm(ana .- fd) / norm(fd) < 1e-7
    end

    # THE SIGNS, which is what separates the right formula from the natural mistake.
    #
    # `rx` appears TWICE — as the local radius through `r_i`, and as the polar reference in
    # `g_pole` — so
    #
    #     ∂T/∂rx = 2β·T·( 1/rx − rx·u_x²/r_i² )
    #
    # and since `r_i² = rx²u_x² + ry²u_y² + rz²u_z² ≥ rx²u_x²`, with equality only on the
    # x-axis itself, the bracket is NON-NEGATIVE everywhere: growing the reference axis makes
    # the whole surface hotter relative to its pole. Dropping the `1/rx` term — keeping only
    # the local-radius half — gives `−2β·T·rx·u_x²/r_i²`, which is strictly NEGATIVE, so this
    # assertion is exactly the one that catches it.
    @test all(>=(0), drx)
    @test minimum(drx) < maximum(drx)          # and it is not a constant
    # `ry` and `rz` enter only through `r_i`, so they can only cool their own direction.
    @test all(<=(0), dry) && all(<=(0), drz)
    # A hotter pole scales the whole map, so that derivative is strictly positive.
    @test all(>(0), dtp)

    @testset "the map does not depend on orientation" begin
        # Exactly zero, not small: `r_i` is a norm and rotation is orthogonal. This is the
        # result that lets `gradient_fit_kind` offer the shape gradient for an ellipsoid whose
        # free set is only the two angles.
        scale = maximum(abs, Tm)
        for key in (:inclination, :position_angle)
            @test maximum(abs, dmap(key)) / scale < 1e-11
        end
        # And the map itself is unchanged by a large reorientation.
        @test maximum(abs, map_of(inclination = 12.0, position_angle = 155.0) .- Tm) /
              scale < 1e-12
    end

    @testset "the shape gradient with a parametric map" begin
        # What the derivatives are FOR. `shape_chi2_fg!` held the map fixed while θ moved,
        # which is exact for a sphere and wrong for an ellipsoid; `parametric_map = true`
        # recomputes the map from θ and contracts its derivative into the θ gradient.
        DATA = joinpath(pkgdir(ROTIR), "demos", "data")
        fs = sort([joinpath(DATA, f) for f in readdir(DATA) if occursin("lam_And", f)])
        data = [readoifits(f)[1, 1] for f in fs[1:2]]
        tepochs = [0.0, 1.0]
        npix = tess.npix
        θ0 = [1.30, 1.05, 0.92, 63.0, 27.0]

        # The objective a parametric shape fit actually minimises: χ² with the map following θ.
        #
        # DISTINCT NAMES for the scratch arrays. Named `gθ`/`gx`/`xm` like the ones below, they
        # aliased them inside `@testset` rather than shadowing them, so every finite-difference
        # evaluation wrote its own gradient over the analytic one computed earlier — and the
        # comparison was then the FD gradient against ITSELF at a perturbed point. The analytic
        # gradient was right all along; the test was measuring the wrong thing.
        function chi2_of(θ)
            gbuf = zeros(5); gxbuf = zeros(npix); xbuf = zeros(npix)
            shape_chi2_fg!(gbuf, gxbuf, xbuf, collect(Float64, θ), data, tess, base, tepochs;
                           parametric_map = true)
        end
        gθ = zeros(5); gx = zeros(npix); xm = zeros(npix)
        c = shape_chi2_fg!(gθ, gx, xm, θ0, data, tess, base, tepochs; parametric_map = true)
        @test isfinite(c) && c > 0
        # The map came back in `xmap`, which is an output in this mode.
        @test maximum(abs, xm .- Tm) / maximum(abs, Tm) < 1e-12

        gfd = FiniteDifferences.grad(central_fdm(5, 1), chi2_of, θ0)[1]
        @test norm(gθ .- gfd) / norm(gfd) < 1e-7
        @test all(abs.((gθ .- gfd) ./ gfd) .< 1e-7)

        # AND THE OLD BEHAVIOUR IS MEASURABLY WRONG, which is why type 1 was refused: holding
        # the map fixed misses the radii's effect on the temperature entirely.
        gfix = zeros(5); gxfix = zeros(npix); xfix = copy(xm)
        shape_chi2_fg!(gfix, gxfix, xfix, θ0, data, tess, base, tepochs)
        @test norm(gfix .- gfd) / norm(gfd) > 1e-2          # 6.6 % as measured
        # The orientation components are unaffected — the map never depended on them.
        @test all(abs.((gfix[4:5] .- gfd[4:5]) ./ gfd[4:5]) .< 1e-7)
    end
end
