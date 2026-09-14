@views function temperature_map_vonZeipel_ellipsoid(stellar_parameters, star; offsets = [0.0,0.0,0.0], T=eltype(star))
    p = convert_params(T, stellar_parameters)
    toff= T.(offsets)'
    r_theta = sqrt.(dropdims(sum(abs2, (star.vertices_xyz[:,5,:] .- toff), dims=2), dims=2));
    rpole = p.radius_x # to check, is it really x?
    theta = star.vertices_spherical[:,5,2];
    g_theta = 1 ./(r_theta.^2)
    g_pole = 1/rpole^2
    star_map = p.tpole*(g_theta/g_pole).^p.beta
    return star_map
  end
  
"""
    temperature_map_vonZeipel_ellipsoid_derivs(stellar_parameters, star; offsets, T)
        -> (Tmap, dT_drx, dT_dry, dT_drz, dT_dbeta, dT_dtpole)

The ellipsoid's von Zeipel map AND its analytic derivatives.

This is what a consistent gradient fit for surface type 1 was missing: `shape_chi2_fg!` has the
projected-geometry derivatives but holds the temperature map FIXED, which is exact only when
the map does not depend on the free parameters. For an ellipsoid it does.

THE MAP. `temperature_map_vonZeipel_ellipsoid` computes

    g_i = 1 / r_i²,   g_pole = 1 / rx²,   T_i = tpole · (g_i/g_pole)^β

so, writing `u` for the tessel's body-frame UNIT direction and `r_i = |(rx·u_x, ry·u_y, rz·u_z)|`,

    T_i = tpole · (rx / r_i)^{2β}

Everything below is that expression differentiated, with `∂r_i/∂rx = rx·u_x²/r_i` and likewise
for `ry`, `rz`:

    ∂T/∂rx    = 2β·T·( 1/rx − rx·u_x²/r_i² )      ← two terms: `rx` sets the POLE reference
    ∂T/∂ry    = −2β·T·ry·u_y²/r_i²                    as well as the local radius
    ∂T/∂rz    = −2β·T·rz·u_z²/r_i²
    ∂T/∂β     = 2·T·log(rx/r_i)
    ∂T/∂tpole = T/tpole

**AND NOTHING ELSE.** `r_i` is a NORM, and rotation is orthogonal, so the map does not depend on
`inclination` or `position_angle` at all — `∂T/∂inc = ∂T/∂PA = 0` exactly. That is worth stating
because it is what makes a shape-only gradient CORRECT for an ellipsoid whose free parameters
are just the two orientation angles: the inconsistency `shape_chi2_fg!` warns about is entirely
in the radii, β and tpole.

`rx` appearing twice — once as the local radius through `r_i` and once as the polar reference
in `g_pole` — is why `∂T/∂rx` has a term the other two lack, and it is the term a
finite-difference check catches when it is dropped.
"""
@views function temperature_map_vonZeipel_ellipsoid_derivs(stellar_parameters, star;
                                                           offsets = [0.0, 0.0, 0.0],
                                                           T = eltype(star))
    p = convert_params(T, stellar_parameters)
    toff = T.(offsets)'
    r = sqrt.(dropdims(sum(abs2, (star.vertices_xyz[:, 5, :] .- toff), dims = 2), dims = 2))
    rx, ry, rz = p.radius_x, p.radius_y, p.radius_z
    β, tp = p.beta, p.tpole
    # The BODY-frame unit direction of each tessel centre. `vertices_spherical` carries it as
    # (r, θ, φ); the radius there is the deformed one, so only the angles are wanted.
    θ = star.vertices_spherical[:, 5, 2]
    φ = star.vertices_spherical[:, 5, 3]
    sθ = sin.(θ)
    ux = sθ .* cos.(φ); uy = sθ .* sin.(φ); uz = cos.(θ)

    Tmap = tp .* (rx ./ r) .^ (2β)
    two_β = T(2) * β
    inv_r2 = one(T) ./ (r .^ 2)
    dT_drx = two_β .* Tmap .* (one(T) / rx .- rx .* ux .^ 2 .* inv_r2)
    dT_dry = (-two_β) .* Tmap .* (ry .* uy .^ 2 .* inv_r2)
    dT_drz = (-two_β) .* Tmap .* (rz .* uz .^ 2 .* inv_r2)
    dT_dbeta = T(2) .* Tmap .* log.(rx ./ r)
    dT_dtpole = Tmap ./ tp
    return Tmap, dT_drx, dT_dry, dT_drz, dT_dbeta, dT_dtpole
end

"""
    temperature_map_vonZeipel_ellipsoid_derivs(stellar_parameters, tessels::tessellation; T)

The same map and derivatives from the TESSELLATION alone, with no `stellar_geometry` built.

The map needs only each tessel's body-frame unit direction and the three radii — `r_i` is
`|(rx·u_x, ry·u_y, rz·u_z)|` by construction — so a fit that evaluates it once per iteration
should not have to build and rotate a whole mesh to get it. Identical values to the
`star` method; that is asserted in test/test_ellipsoid_map_derivs.jl.
"""
@views function temperature_map_vonZeipel_ellipsoid_derivs(stellar_parameters,
                                                           tessels::tessellation;
                                                           T = eltype(tessels))
    p = convert_params(T, stellar_parameters)
    rx, ry, rz = p.radius_x, p.radius_y, p.radius_z
    β, tp = p.beta, p.tpole
    ux = tessels.unit_xyz[:, 5, 1]
    uy = tessels.unit_xyz[:, 5, 2]
    uz = tessels.unit_xyz[:, 5, 3]
    r = sqrt.((rx .* ux) .^ 2 .+ (ry .* uy) .^ 2 .+ (rz .* uz) .^ 2)

    Tmap = tp .* (rx ./ r) .^ (2β)
    two_β = T(2) * β
    inv_r2 = one(T) ./ (r .^ 2)
    dT_drx = two_β .* Tmap .* (one(T) / rx .- rx .* ux .^ 2 .* inv_r2)
    dT_dry = (-two_β) .* Tmap .* (ry .* uy .^ 2 .* inv_r2)
    dT_drz = (-two_β) .* Tmap .* (rz .* uz .^ 2 .* inv_r2)
    dT_dbeta = T(2) .* Tmap .* log.(rx ./ r)
    dT_dtpole = Tmap ./ tp
    return Tmap, dT_drx, dT_dry, dT_drz, dT_dbeta, dT_dtpole
end

"""
    shape_map_and_derivs(star_params, tessels, stype) -> (xmap, dxmap_dθ) or nothing

The parametric temperature map and its derivative w.r.t. `shape_chi2_fg!`'s θ, or `nothing`
when the map does not depend on θ at all.

Per surface type, with θ as `shape_chi2_fg!` lays it out:

  * **0, sphere** — `θ = [radius, inc, PA]` and the map is a uniform `tpole`. It depends on
    none of them, so `nothing`: holding the map fixed was already exact here, which is why the
    sphere has always been offered a gradient fit.
  * **1, ellipsoid** — `θ = [rx, ry, rz, inc, PA]`. The map depends on the three radii and NOT
    on the orientation (`r_i` is a norm; see
    [`temperature_map_vonZeipel_ellipsoid_derivs`](@ref)), so the last two columns are zero.
  * **2, rapid rotator** — `nothing`, and deliberately: its map has a Zygote rrule
    (`vonzeipel_map`) and the parametric path uses it. Reaching here would mean the caller
    chose the wrong fitter, so it is an error rather than a silent zero.
"""
function shape_map_and_derivs(star_params, tessels::tessellation{T}, stype::Integer) where {T}
    stype == 0 && return nothing
    stype == 2 && error("shape_map_and_derivs: the rapid rotator's map is differentiated by " *
                        "Zygote through `vonzeipel_map`; use the parametric fit path")
    stype == 1 || error("shape_map_and_derivs: no parametric map for surface_type $(stype)")
    Tmap, drx, dry, drz, _, _ =
        temperature_map_vonZeipel_ellipsoid_derivs(star_params, tessels; T = T)
    d = zeros(T, length(Tmap), 5)
    d[:, 1] .= drx; d[:, 2] .= dry; d[:, 3] .= drz     # columns 4 and 5 (inc, PA) stay zero
    return Tmap, d
end
