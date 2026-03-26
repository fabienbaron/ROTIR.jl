# Analytical shape parameter gradients for joint shape + map optimization.
# Computes ∂χ²/∂θ where θ are shape parameters (surface-type dependent)
# via the chain rule through projected vertex positions and soft visibility.
#
# Supported surface types:
#   0 (Sphere):         θ = [radius, inc, PA]               (3 params)
#   1 (Ellipsoid):      θ = [rx, ry, rz, inc, PA]           (5 params)
#   2 (Rapid Rotator):  θ = [rpole, ω, inc, PA]             (4 params)
#
# Rotation convention: xyz = su * R  (right-multiply), matching rot_vertex in geometry.jl

# ─── Rotation matrix and derivatives ─────────────────────────────────────────
# These are identical to rot_vertex(ψ, inc, PA) from geometry.jl, written with
# named angles for clarity. ψ=rotation phase, inc and PA in radians.

"""
    rotation_matrix(ψ, inc, PA)
Rotation matrix R(ψ, inc, PA) — same as rot_vertex(ψ, inc, PA).
All angles in radians.
"""
function rotation_matrix(ψ, inc, PA)
    c1, s1 = cos(ψ), sin(ψ)
    c2, s2 = cos(inc), sin(inc)
    c3, s3 = cos(PA), sin(PA)
    R = [-s1*c2*s3+c1*c3  s1*c3*c2+c1*s3  -s1*s2;
         -c1*c2*s3-s1*c3  c1*c3*c2-s1*s3  -c1*s2;
                  -s2*s3           s2*c3       c2 ]
    return R
end

"""
    dR_dinc(ψ, inc, PA)
∂R/∂inc — derivative of rotation matrix w.r.t. inclination (in radians).
"""
function dR_dinc(ψ, inc, PA)
    c1, s1 = cos(ψ), sin(ψ)
    c2, s2 = cos(inc), sin(inc)
    c3, s3 = cos(PA), sin(PA)
    dR = [ s1*s2*s3  -s1*s2*c3  -s1*c2;
           c1*s2*s3  -c1*s2*c3  -c1*c2;
              -c2*s3      c2*c3    -s2 ]
    return dR
end

"""
    dR_dPA(ψ, inc, PA)
∂R/∂PA — derivative of rotation matrix w.r.t. position angle (in radians).
"""
function dR_dPA(ψ, inc, PA)
    c1, s1 = cos(ψ), sin(ψ)
    c2, s2 = cos(inc), sin(inc)
    c3, s3 = cos(PA), sin(PA)
    dR = [-s1*c2*c3-c1*s3  -s1*s3*c2+c1*c3  0;
          -c1*c2*c3+s1*s3  -c1*s3*c2-s1*c3  0;
                   -s2*c3          -s2*s3     0 ]
    return dR
end

# ─── Rapid rotator radius and derivative ─────────────────────────────────────

"""
    f_rapid_rot_and_deriv(x)
Compute f(x) = 3cos((π+acos(x))/3)/x and its derivative f'(x).
This is the rapid rotator radius function: r = rpole * f(ω * sin(θ)).

f'(x) = d/dx [3cos(α)/x] where α = (π + acos(x))/3
      = [sin(α)/√(1-x²) - 3cos(α)/x] / x
"""
function f_rapid_rot_and_deriv(x::T) where T
    if abs(x) < T(1e-12)
        return one(T), zero(T)
    end
    α = (T(π) + acos(x)) / 3
    sα, cα = sincos(α)
    f = 3 * cα / x
    onemx2 = one(T) - x^2
    if onemx2 < T(1e-30)
        df = -f / x
    else
        df = (sα / sqrt(onemx2) - f) / x
    end
    return f, df
end

# ─── Projected vertex coordinates and their parameter derivatives ────────────

"""
    projected_vertices_and_derivs(tessels, star_params, t; nparams)

Compute projected quad vertices and their analytical derivatives w.r.t. shape parameters.

Surface type 0 (Sphere):    θ = [radius, inc (deg), PA (deg)]         → nparams=3
Surface type 1 (Ellipsoid): θ = [rx, ry, rz, inc (deg), PA (deg)]    → nparams=5
Surface type 2 (Rapid Rot): θ = [rpole, ω, inc (deg), PA (deg)]      → nparams=4

Returns:
- `projx, projy`: (npix, 4) projected vertex coordinates
- `dprojx_dθ, dprojy_dθ`: (npix, 4, nparams) vertex derivatives
- `normals_z`: (npix,) z-component of face normals
- `dnz_dθ`: (npix, nparams) derivative of normalized normals_z w.r.t. θ
"""
function projected_vertices_and_derivs(tessels::tessellation{T}, star_params, t;
                                        nparams::Int=0) where T
    npix = tessels.npix
    deg2rad = T(π/180)
    stype = star_params.surface_type

    # Determine nparams from surface type if not specified
    if nparams == 0
        nparams = stype == 1 ? 5 : (stype == 2 ? 4 : 3)
    end

    # Indices for inc and PA in θ vector (surface-type dependent)
    if stype == 1       # Ellipsoid: θ = [rx, ry, rz, inc, PA]
        idx_inc = 4; idx_PA = 5
    elseif stype == 2   # Rapid rotator: θ = [rpole, ω, inc, PA]
        idx_inc = 3; idx_PA = 4
    else                # Sphere: θ = [radius, inc, PA]
        idx_inc = 2; idx_PA = 3
    end

    inc_rad = T(star_params.inclination) * deg2rad
    PA_rad  = T(star_params.position_angle) * deg2rad
    ψ = T(2π) * T(t) / T(star_params.rotation_period)

    # Rotation matrix and derivatives (chain rule: inc/PA in degrees, so multiply by deg2rad)
    R    = rotation_matrix(ψ, inc_rad, PA_rad)
    dR_di = dR_dinc(ψ, inc_rad, PA_rad) .* deg2rad
    dR_dp = dR_dPA(ψ, inc_rad, PA_rad)  .* deg2rad

    # For rapid rotator, precompute ω
    if stype == 2
        rpole = T(star_params.rpole)
        ω = T(star_params.frac_escapevel)
    end

    # Compute scaled unit vertices (all 5 vertices, needed for normals)
    su_full = similar(tessels.unit_xyz)  # (npix, 5, 3)
    if stype == 1  # Ellipsoid
        rx = T(star_params.radius_x)
        ry = T(star_params.radius_y)
        rz = T(star_params.radius_z)
        su_full[:, :, 1] = rx .* tessels.unit_xyz[:, :, 1]
        su_full[:, :, 2] = ry .* tessels.unit_xyz[:, :, 2]
        su_full[:, :, 3] = rz .* tessels.unit_xyz[:, :, 3]
    elseif stype == 2  # Rapid rotator
        # Per-vertex radius: r[p,v] = rpole * f(ω * sin(θ_colat[p,v]))
        for p in 1:npix, v in 1:5
            sinθ_v = sin(tessels.unit_spherical[p, v, 2])
            ωsinθ = ω * sinθ_v
            r_pv = abs(ωsinθ) < T(1e-12) ? rpole : rpole * first(f_rapid_rot_and_deriv(ωsinθ))
            su_full[p, v, 1] = r_pv * tessels.unit_xyz[p, v, 1]
            su_full[p, v, 2] = r_pv * tessels.unit_xyz[p, v, 2]
            su_full[p, v, 3] = r_pv * tessels.unit_xyz[p, v, 3]
        end
    else  # Sphere
        rad = T(star_params.radius)
        su_full .= rad .* tessels.unit_xyz
    end

    # Rotate all vertices: xyz = su * R
    xyz_full = reshape(reshape(su_full, npix*5, 3) * R, npix, 5, 3)

    # Projected coordinates (vertices 1:4 only)
    projx = xyz_full[:, 1:4, 1]  # (npix, 4)
    projy = xyz_full[:, 1:4, 2]  # (npix, 4)

    # ─── Normals ─────────────────────────────────────────────────────────────
    vecAC = xyz_full[:, 3, :] - xyz_full[:, 1, :]  # (npix, 3)
    vecBD = xyz_full[:, 4, :] - xyz_full[:, 2, :]  # (npix, 3)
    nx_comp = vecAC[:,2].*vecBD[:,3] - vecAC[:,3].*vecBD[:,2]
    ny_comp = vecAC[:,3].*vecBD[:,1] - vecAC[:,1].*vecBD[:,3]
    nz_comp = vecAC[:,1].*vecBD[:,2] - vecAC[:,2].*vecBD[:,1]
    nn = sqrt.(nx_comp.^2 + ny_comp.^2 + nz_comp.^2)
    normals_z = nz_comp ./ nn
    nn3 = nn.^3

    # ─── Vertex derivatives ──────────────────────────────────────────────────
    u = tessels.unit_xyz[:, 1:4, :]  # (npix, 4, 3)
    dprojx_dθ = zeros(T, npix, 4, nparams)
    dprojy_dθ = zeros(T, npix, 4, nparams)
    dnz_dθ = zeros(T, npix, nparams)

    for j in 1:nparams
        if stype == 1  # Ellipsoid
            if j <= 3
                # ∂proj/∂r_j: xyz = su * R, su_j = r_j * u_j
                # ∂projx/∂r_j = R[j,1] * u_j,  ∂projy/∂r_j = R[j,2] * u_j
                dprojx_dθ[:, :, j] = R[j,1] .* u[:, :, j]
                dprojy_dθ[:, :, j] = R[j,2] .* u[:, :, j]
                # For dnz: dverts_full[:,:,c] = R[j,c] * unit_xyz[:,:,j]
                dverts = zeros(T, npix, 5, 3)
                for c in 1:3
                    dverts[:, :, c] = R[j, c] .* tessels.unit_xyz[:, :, j]
                end
            elseif j == 4  # inc
                dxyz = reshape(reshape(su_full, npix*5, 3) * dR_di, npix, 5, 3)
                dprojx_dθ[:, :, j] = dxyz[:, 1:4, 1]
                dprojy_dθ[:, :, j] = dxyz[:, 1:4, 2]
                dverts = dxyz
            else  # PA
                dxyz = reshape(reshape(su_full, npix*5, 3) * dR_dp, npix, 5, 3)
                dprojx_dθ[:, :, j] = dxyz[:, 1:4, 1]
                dprojy_dθ[:, :, j] = dxyz[:, 1:4, 2]
                dverts = dxyz
            end

        elseif stype == 2  # Rapid rotator: θ = [rpole, ω, inc, PA]
            if j == 1  # rpole
                # ∂r/∂rpole = f(ω sinθ), so ∂su/∂rpole = (r/rpole) * unit_xyz = su/rpole
                dsu = su_full ./ rpole
                dxyz = reshape(reshape(dsu, npix*5, 3) * R, npix, 5, 3)
                dprojx_dθ[:, :, j] = dxyz[:, 1:4, 1]
                dprojy_dθ[:, :, j] = dxyz[:, 1:4, 2]
                dverts = dxyz
            elseif j == 2  # ω
                # ∂r/∂ω = rpole * f'(ωsinθ) * sinθ
                dsu = zeros(T, npix, 5, 3)
                for p in 1:npix, v in 1:5
                    sinθ_v = sin(tessels.unit_spherical[p, v, 2])
                    ωsinθ = ω * sinθ_v
                    if abs(ωsinθ) < T(1e-12)
                        dr_dω = zero(T)
                    else
                        _, f_deriv = f_rapid_rot_and_deriv(ωsinθ)
                        dr_dω = rpole * f_deriv * sinθ_v
                    end
                    dsu[p, v, 1] = dr_dω * tessels.unit_xyz[p, v, 1]
                    dsu[p, v, 2] = dr_dω * tessels.unit_xyz[p, v, 2]
                    dsu[p, v, 3] = dr_dω * tessels.unit_xyz[p, v, 3]
                end
                dxyz = reshape(reshape(dsu, npix*5, 3) * R, npix, 5, 3)
                dprojx_dθ[:, :, j] = dxyz[:, 1:4, 1]
                dprojy_dθ[:, :, j] = dxyz[:, 1:4, 2]
                dverts = dxyz
            elseif j == idx_inc  # inc
                dxyz = reshape(reshape(su_full, npix*5, 3) * dR_di, npix, 5, 3)
                dprojx_dθ[:, :, j] = dxyz[:, 1:4, 1]
                dprojy_dθ[:, :, j] = dxyz[:, 1:4, 2]
                dverts = dxyz
            else  # PA
                dxyz = reshape(reshape(su_full, npix*5, 3) * dR_dp, npix, 5, 3)
                dprojx_dθ[:, :, j] = dxyz[:, 1:4, 1]
                dprojy_dθ[:, :, j] = dxyz[:, 1:4, 2]
                dverts = dxyz
            end

        else  # Sphere: θ = [radius, inc, PA]
            if j == 1  # radius
                dsu = tessels.unit_xyz
                dxyz = reshape(reshape(dsu, npix*5, 3) * R, npix, 5, 3)
                dprojx_dθ[:, :, j] = dxyz[:, 1:4, 1]
                dprojy_dθ[:, :, j] = dxyz[:, 1:4, 2]
                dverts = dxyz
            elseif j == idx_inc
                dxyz = reshape(reshape(su_full, npix*5, 3) * dR_di, npix, 5, 3)
                dprojx_dθ[:, :, j] = dxyz[:, 1:4, 1]
                dprojy_dθ[:, :, j] = dxyz[:, 1:4, 2]
                dverts = dxyz
            else  # PA
                dxyz = reshape(reshape(su_full, npix*5, 3) * dR_dp, npix, 5, 3)
                dprojx_dθ[:, :, j] = dxyz[:, 1:4, 1]
                dprojy_dθ[:, :, j] = dxyz[:, 1:4, 2]
                dverts = dxyz
            end
        end

        # ─── Normal z derivative ─────────────────────────────────────────────
        # ∂(nz/‖n‖)/∂θ = ∂nz/∂θ / ‖n‖ − nz·(n·∂n/∂θ) / ‖n‖³
        dAC = dverts[:, 3, :] - dverts[:, 1, :]
        dBD = dverts[:, 4, :] - dverts[:, 2, :]
        dnx = dAC[:,2].*vecBD[:,3] + vecAC[:,2].*dBD[:,3] - dAC[:,3].*vecBD[:,2] - vecAC[:,3].*dBD[:,2]
        dny = dAC[:,3].*vecBD[:,1] + vecAC[:,3].*dBD[:,1] - dAC[:,1].*vecBD[:,3] - vecAC[:,1].*dBD[:,3]
        dnz = dAC[:,1].*vecBD[:,2] + vecAC[:,1].*dBD[:,2] - dAC[:,2].*vecBD[:,1] - vecAC[:,2].*dBD[:,1]
        n_dot_dn = nx_comp.*dnx + ny_comp.*dny + nz_comp.*dnz
        dnz_dθ[:, j] = dnz ./ nn - nz_comp .* n_dot_dn ./ nn3
    end

    return projx, projy, dprojx_dθ, dprojy_dθ, normals_z, dnz_dθ
end

# ─── Shape chi2 + gradient ────────────────────────────────────────────────────

"""
    shape_chi2_fg!(grad_θ, grad_xmap, xmap, θ, data_epochs, tessels, star_params_base, tepochs;
                    κ=50.0, verbose=false)

Compute chi2 and gradients w.r.t. both shape parameters θ and surface map xmap,
using the fused two-pass approach (no polyft matrix).

θ layout depends on surface_type:
  0 (Sphere):    θ = [radius, inc, PA]
  1 (Ellipsoid): θ = [rx, ry, rz, inc, PA]
  2 (Rapid Rot): θ = [rpole, ω, inc, PA]

Returns chi2 (scalar). Modifies grad_θ and grad_xmap in place.
"""
function shape_chi2_fg!(grad_θ::Vector{T}, grad_xmap::Vector{T},
                         xmap::Vector{T}, θ::Vector{T},
                         data_epochs, tessels::tessellation{T},
                         star_params_base, tepochs;
                         κ::T=T(50), verbose::Bool=false) where T

    nparams = length(θ)
    npix = tessels.npix
    nepochs = length(data_epochs)
    stype = star_params_base.surface_type

    # Build star_params from θ
    if stype == 1
        star_params = merge(star_params_base, (radius_x=θ[1], radius_y=θ[2], radius_z=θ[3],
                                                inclination=θ[4], position_angle=θ[5]))
    elseif stype == 2
        star_params = merge(star_params_base, (rpole=θ[1], frac_escapevel=θ[2],
                                                inclination=θ[3], position_angle=θ[4]))
    else
        star_params = merge(star_params_base, (radius=θ[1],
                                                inclination=θ[2], position_angle=θ[3]))
    end

    grad_θ .= zero(T)
    grad_xmap .= zero(T)
    total_chi2 = zero(T)

    for ep in 1:nepochs
        t = tepochs[ep]
        data = data_epochs[ep]

        # Compute projected vertices + derivatives
        projx, projy, dprojx_dθ, dprojy_dθ, nz, dnz_dθ =
            projected_vertices_and_derivs(tessels, star_params, t, nparams=nparams)

        # Soft visibility
        vis_w, sig_args = soft_visibility(nz, κ=κ)

        # Select visible pixels
        vis_threshold = T(1e-4)
        indx = findall(vis_w .> vis_threshold)
        nvis = length(indx)

        # Weighted pixel values
        xw = xmap[indx] .* vis_w[indx]

        # UV frequencies
        kx = data.uv[1,:] * T(-π / (180*3600000))
        ky = data.uv[2,:] * T( π / (180*3600000))
        nuv = length(kx)
        k2_inv_im = precompute_k2_inv_im(kx, ky)

        # Forward pass
        pjx = projx[indx, :]
        pjy = projy[indx, :]
        F = Vector{Complex{T}}(undef, nuv)
        polyflux_local = Vector{T}(undef, nvis)
        compute_polyflux_and_cvis!(F, polyflux_local, kx, ky, k2_inv_im, pjx, pjy, xw)

        flux = dot(polyflux_local, xw)
        cvis_model = F / flux

        # Observables + chi2
        v2_model = abs2.(cvis_model[data.indx_v2])
        t3_model = cvis_model[data.indx_t3_1] .* cvis_model[data.indx_t3_2] .* cvis_model[data.indx_t3_3]
        t3amp_model = abs.(t3_model)
        t3phi_model = angle.(t3_model) * T(180/π)

        chi2_v2 = sum(abs2, (v2_model - data.v2) ./ data.v2_err)
        chi2_t3amp = sum(abs2, (t3amp_model - data.t3amp) ./ data.t3amp_err)
        chi2_t3phi = sum(abs2, mod360(t3phi_model - data.t3phi) ./ data.t3phi_err)
        total_chi2 += chi2_v2 + chi2_t3amp + chi2_t3phi

        # Adjoint signal: ∂χ²/∂cvis
        adj_cvis = zeros(Complex{T}, nuv)
        for i in eachindex(data.indx_v2)
            k = data.indx_v2[i]
            adj_cvis[k] += 4 * (v2_model[i] - data.v2[i]) / data.v2_err[i]^2 * conj(cvis_model[k])
        end
        t3amp_res = 2 * (t3amp_model - data.t3amp) ./ data.t3amp_err.^2
        for i in eachindex(data.indx_t3_1)
            k1 = data.indx_t3_1[i]; k2 = data.indx_t3_2[i]; k3 = data.indx_t3_3[i]
            c1 = cvis_model[k1]; c2 = cvis_model[k2]; c3 = cvis_model[k3]
            a1 = abs(c1); a2 = abs(c2); a3 = abs(c3)
            adj_cvis[k1] += t3amp_res[i] * conj(c1)/a1 * a2 * a3
            adj_cvis[k2] += t3amp_res[i] * conj(c2)/a2 * a1 * a3
            adj_cvis[k3] += t3amp_res[i] * conj(c3)/a3 * a1 * a2
        end
        t3phi_res = mod360(t3phi_model - data.t3phi) ./ data.t3phi_err.^2
        for i in eachindex(data.indx_t3_1)
            k1 = data.indx_t3_1[i]; k2 = data.indx_t3_2[i]; k3 = data.indx_t3_3[i]
            c1 = cvis_model[k1]; c2 = cvis_model[k2]; c3 = cvis_model[k3]
            t3i = t3_model[i]
            factor = t3phi_res[i] / abs2(t3i) * conj(t3i)
            adj_cvis[k1] -= T(360/π) * im * factor * c2 * c3
            adj_cvis[k2] -= T(360/π) * im * factor * c1 * c3
            adj_cvis[k3] -= T(360/π) * im * factor * c1 * c2
        end

        # Scale adjoint for flux normalization: adj_F = adj_cvis / flux
        adj_F = adj_cvis / flux

        # Adjoint pass for xmap gradient (uses compute_adjoint_cvis!)
        grad_xw = Vector{T}(undef, nvis)
        compute_adjoint_cvis!(grad_xw, adj_F, kx, ky, k2_inv_im, pjx, pjy, polyflux_local)
        flux_adj = -dot(xw, grad_xw) / flux
        grad_xw .+= flux_adj * polyflux_local
        grad_xmap[indx] .+= vis_w[indx] .* grad_xw

        # Adjoint pass for vertex position gradients (non-DC Fourier terms)
        grad_projx = Matrix{T}(undef, nvis, 4)
        grad_projy = Matrix{T}(undef, nvis, 4)
        compute_adjoint_vertices!(grad_projx, grad_projy, adj_F, kx, ky, k2_inv_im, pjx, pjy, xw, polyflux_local)

        # Flux correction for vertex gradients: ∂χ²/∂flux * ∂flux/∂vertex
        # flux = Σ_p xw_p * polyflux_p, where polyflux_p = shoelace area
        # ∂polyflux_p/∂x[j] = 0.5*(y[j+1] - y[j-1])
        # ∂polyflux_p/∂y[j] = 0.5*(x[j-1] - x[j+1])
        @inbounds for (li, pi) in enumerate(indx)
            fa_xw = flux_adj * xw[li]
            for j in 1:4
                jp = mod1(j+1, 4)
                jm = mod1(j-1, 4)
                grad_projx[li, j] += fa_xw * T(0.5) * (pjy[li, jp] - pjy[li, jm])
                grad_projy[li, j] += fa_xw * T(0.5) * (pjx[li, jm] - pjx[li, jp])
            end
        end

        # Chain rule via vertex positions: ∂χ²/∂θ_j = Σ_p,v (∂χ²/∂projx * ∂projx/∂θ + ∂χ²/∂projy * ∂projy/∂θ)
        for j in 1:nparams
            s = zero(T)
            for (li, pi) in enumerate(indx)
                for v in 1:4
                    s += grad_projx[li, v] * dprojx_dθ[pi, v, j]
                    s += grad_projy[li, v] * dprojy_dθ[pi, v, j]
                end
            end
            grad_θ[j] += s
        end

        # Chain rule via soft visibility: ∂χ²/∂θ via w
        # dchi2_dw = dchi2_dxw * xmap, and dw/dθ = dσ/d(κnz) * κ * dnz/dθ
        dchi2_dw_visible = grad_xw .* xmap[indx]  # ∂χ²/∂(w) for visible pixels
        dw_dsigarg = dsigmoid.(sig_args[indx])
        for j in 1:nparams
            s = zero(T)
            for (li, pi) in enumerate(indx)
                s += dchi2_dw_visible[li] * dw_dsigarg[li] * κ * dnz_dθ[pi, j]
            end
            grad_θ[j] += s
        end

        if verbose
            println("Epoch $ep — V2: ", chi2_v2/data.nv2, "\tT3A: ", chi2_t3amp/data.nt3amp,
                    "\tT3P: ", chi2_t3phi/data.nt3phi, "\tFlux: ", flux)
        end
    end

    return total_chi2
end

# ─── Joint reconstruction ────────────────────────────────────────────────────

"""
    joint_reconstruct_oi(xmap_start, θ_start, data, tessels_base, star_params_base, tepochs;
                          maxiter_xmap=200, maxiter_θ=50, nouter=5, reg_weight=1e-5, κ=50.0)

Alternating optimization of surface map and shape parameters.
Outer loop alternates between:
  1. xmap step: optimize surface map with fixed shape (VMLMB)
  2. θ step: optimize shape parameters with fixed map (VMLMB with bounds)

θ layout depends on surface_type (see shape_chi2_fg!).
"""
function joint_reconstruct_oi(xmap_start::Vector{T}, θ_start::Vector{T},
                               data, tessels_base, star_params_base, tepochs;
                               maxiter_xmap::Int=200, maxiter_θ::Int=50,
                               nouter::Int=5, reg_weight::T=T(1e-5),
                               θ_lower::Union{Vector{T},Nothing}=nothing,
                               θ_upper::Union{Vector{T},Nothing}=nothing,
                               κ::T=T(50), verbose::Bool=true) where T

    xmap = copy(xmap_start)
    θ = copy(θ_start)
    npix = length(xmap)
    nparams = length(θ)
    stype = star_params_base.surface_type

    if θ_lower === nothing
        θ_lower = fill(T(-Inf), nparams)
    end
    if θ_upper === nothing
        θ_upper = fill(T(Inf), nparams)
    end

    for outer in 1:nouter
        if verbose
            println("\n=== Outer iteration $outer/$nouter ===")
            println("θ = $θ")
        end

        # Step 1: Rebuild geometry with current θ
        if stype == 1
            star_params = merge(star_params_base, (radius_x=θ[1], radius_y=θ[2], radius_z=θ[3],
                                                    inclination=θ[4], position_angle=θ[5]))
        elseif stype == 2
            star_params = merge(star_params_base, (rpole=θ[1], frac_escapevel=θ[2],
                                                    inclination=θ[3], position_angle=θ[4]))
        else
            star_params = merge(star_params_base, (radius=θ[1],
                                                    inclination=θ[2], position_angle=θ[3]))
        end
        n = npix2n(npix)
        tessels = tessellation_healpix(n)
        stars = create_star_multiepochs(tessels, star_params, tepochs)
        setup_oi!(data, stars)

        # Step 2: Optimize xmap with fixed geometry
        regularizers = [["tv2", reg_weight, tv_neighbours_healpix(n), 1:npix]]
        xmap = image_reconstruct_oi(xmap, data, stars, maxiter=maxiter_xmap,
                                     regularizers=regularizers, verbose=verbose)

        # Step 3: Optimize θ with fixed xmap
        crit_θ = (θ_vec, g_vec) -> begin
            g_xmap_dummy = zeros(T, npix)
            chi2 = shape_chi2_fg!(g_vec, g_xmap_dummy, xmap, θ_vec, data, tessels,
                                   star_params_base, tepochs, κ=κ, verbose=false)
            return chi2
        end
        θ = OptimPackNextGen.vmlmb(crit_θ, θ, verb=verbose, lower=θ_lower, upper=θ_upper,
                                    maxiter=maxiter_θ, blmvm=false, gtol=(0, 1e-6))

        if verbose
            println("θ after step: $θ")
        end
    end

    return xmap, θ
end
