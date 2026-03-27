# Binary star forward model
# Computes combined complex visibilities for two stars using phase-shifted polygon FTs.
# Star 1 is at the origin; star 2 is offset by its orbital position.

"""
    orbit_to_rotir_offset(bparams, tepoch_jd) -> (offset_x, offset_y)

Convert orbital position to ROTIR's projected coordinate frame.
Returns the secondary's offset relative to the primary in mas, in ROTIR's (West, North) frame.

The orbital code (`binary_orbit_abs`) returns positions where x=North, y=East.
ROTIR's projected coordinates have proj_west=West, proj_north=North.
"""
function orbit_to_rotir_offset(bparams, tepoch_jd)
    x1, y1, z1, x2, y2, z2 = binary_orbit_abs(bparams, tepoch_jd)
    dx_north = x2 - x1
    dy_east  = y2 - y1
    offset_x = -dy_east   # West = -East
    offset_y = dx_north   # North
    return offset_x, offset_y
end

"""
    binary_phase_shift(uv, offset_x, offset_y) -> Vector{Complex}

Compute per-baseline phase shift for a star displaced by (offset_x, offset_y) mas
in ROTIR's projected frame (West, North).

Uses the same kx/ky sign convention as the polygon FT in `setup_polyft_single`.
"""
function binary_phase_shift(uv, offset_x, offset_y; T=Float32)
    C = T(180 * 3600000)
    kx = uv[1,:] .* T(-pi / C)
    ky = uv[2,:] .* T( pi / C)
    return cis.(-T(pi) .* (kx .* offset_x .+ ky .* offset_y))
end

"""
    binary_cvis(x1, star1, x2, star2, phase_shift) -> Vector{Complex}

Compute combined complex visibilities for a binary system.
Star 1 is at the origin; star 2's visibilities are multiplied by `phase_shift`.
Both temperature maps are flux-normalized together.
"""
function binary_cvis(x1, star1, x2, star2, phase_shift)
    indx1 = star1.index_quads_visible
    xw1 = x1[indx1] .* star1.vis_weights[indx1]
    flux1 = dot(star1.polyflux, xw1)
    F1 = star1.polyft * xw1

    indx2 = star2.index_quads_visible
    xw2 = x2[indx2] .* star2.vis_weights[indx2]
    flux2 = dot(star2.polyflux, xw2)
    F2 = star2.polyft * xw2

    return (F1 .+ F2 .* phase_shift) ./ (flux1 + flux2)
end

"""
    binary_observables(x1, star1, x2, star2, data, phase_shift) -> (v2, t3amp, t3phi)

Compute model observables (V2, T3amp, T3phi) for a binary system.
Uses `cvis_to_obs` (shared with single-star path) for the cvis→observables step.
"""
function binary_observables(x1, star1, x2, star2, data, phase_shift)
    cvis = binary_cvis(x1, star1, x2, star2, phase_shift)
    return cvis_to_obs(cvis, data)
end

"""
    binary_chi2_f(x1, star1, x2, star2, data, phase_shift; verbose=false) -> Float

Compute chi-squared for a binary model against interferometric data.
"""
function binary_chi2_f(x1, star1, x2, star2, data, phase_shift; verbose::Bool=false)
    v2_model, t3amp_model, t3phi_model = binary_observables(x1, star1, x2, star2, data, phase_shift)
    chi2_v2 = sum(abs2, (v2_model .- data.v2) ./ data.v2_err)
    chi2_t3amp = sum(abs2, (t3amp_model .- data.t3amp) ./ data.t3amp_err)
    chi2_t3phi = sum(abs2, mod360(t3phi_model .- data.t3phi) ./ data.t3phi_err)
    if verbose
        printstyled(@sprintf("V2: %.4f ", chi2_v2/data.nv2), color=:red)
        printstyled(@sprintf("T3A: %.4f ", chi2_t3amp/data.nt3amp), color=:blue)
        printstyled(@sprintf("T3P: %.4f\n", chi2_t3phi/data.nt3phi), color=:green)
    end
    return chi2_v2 + chi2_t3amp + chi2_t3phi
end
