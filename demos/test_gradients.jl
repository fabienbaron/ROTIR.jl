#!/usr/bin/env julia
# Gradient validation for shape_gradient.jl using FiniteDifferences.jl
# Tests: map gradient, shape gradient (sphere/ellipsoid/rapid rotator)

using ROTIR
using FiniteDifferences
using LinearAlgebra

# ─── Load data (single epoch for speed) ─────────────────────────────────────
oifitsfiles = ["./data/2011Sep02.lam_And_prepped.oifits"]
data_all = readoifits_multiepochs(oifitsfiles)
data = data_all[1, :]
nepochs = length(data)
tepochs = [d.mean_mjd for d in data]
tepochs = tepochs .- tepochs[1]

n = 3
tessels = tessellation_healpix(n)
npix = tessels.npix

# ─── Helper: compute chi2 from θ (shape params) and xmap ────────────────────
function chi2_from_θ_and_xmap(θ, xmap, data, tessels, star_params_base, tepochs; κ=50.0)
    grad_θ = zeros(length(θ))
    grad_xmap = zeros(length(xmap))
    chi2 = shape_chi2_fg!(grad_θ, grad_xmap, xmap, θ, data, tessels,
                           star_params_base, tepochs, κ=κ)
    return chi2
end

function gradients_from_θ_and_xmap(θ, xmap, data, tessels, star_params_base, tepochs; κ=50.0)
    grad_θ = zeros(length(θ))
    grad_xmap = zeros(length(xmap))
    chi2 = shape_chi2_fg!(grad_θ, grad_xmap, xmap, θ, data, tessels,
                           star_params_base, tepochs, κ=κ)
    return chi2, grad_θ, grad_xmap
end

# ─── Test function ───────────────────────────────────────────────────────────
function test_gradient(label, star_params_base, θ_start, θ_names; κ=50.0)
    println("\n" * "="^70)
    println("Testing: $label")
    println("="^70)

    stars = create_star_multiepochs(tessels, star_params_base, tepochs)
    tmap = parametric_temperature_map(star_params_base, stars[1])
    setup_oi!(data, stars)

    xmap = copy(tmap)
    θ = copy(θ_start)

    # Analytical gradients
    chi2_a, grad_θ_a, grad_xmap_a = gradients_from_θ_and_xmap(θ, xmap, data, tessels,
                                                                 star_params_base, tepochs, κ=κ)
    println("chi2 = $chi2_a")

    # FD gradient of map (test a random subset of 10 pixels)
    println("\n--- Map gradient (10 random pixels) ---")
    fdm = central_fdm(5, 1)  # 5-point central stencil, 1st derivative
    test_indices = sort(rand(1:npix, 10))

    for idx in test_indices
        fd_grad = fdm(0.0) do δ
            xmap_pert = copy(xmap)
            xmap_pert[idx] += δ
            chi2_from_θ_and_xmap(θ, xmap_pert, data, tessels, star_params_base, tepochs, κ=κ)
        end
        rel = abs(grad_xmap_a[idx]) > 1e-10 ?
              abs(grad_xmap_a[idx] - fd_grad) / abs(grad_xmap_a[idx]) : abs(grad_xmap_a[idx] - fd_grad)
        status = rel < 1e-4 ? "✓" : (rel < 1e-2 ? "~" : "✗")
        println("  pixel $idx: analytical=$(grad_xmap_a[idx])  fd=$fd_grad  rel_err=$rel  $status")
    end

    # FD gradient of shape params
    println("\n--- Shape gradient ---")
    for j in eachindex(θ)
        fd_grad = fdm(0.0) do δ
            θ_pert = copy(θ)
            θ_pert[j] += δ
            chi2_from_θ_and_xmap(θ_pert, xmap, data, tessels, star_params_base, tepochs, κ=κ)
        end
        rel = abs(grad_θ_a[j]) > 1e-10 ?
              abs(grad_θ_a[j] - fd_grad) / abs(grad_θ_a[j]) : abs(grad_θ_a[j] - fd_grad)
        status = rel < 1e-4 ? "✓" : (rel < 1e-2 ? "~" : "✗")
        println("  $(θ_names[j]): analytical=$(grad_θ_a[j])  fd=$fd_grad  rel_err=$rel  $status")
    end
end

# ─── Test 1: Sphere ─────────────────────────────────────────────────────────
star_params_sphere = (
    surface_type    = 0,
    radius          = 1.37,
    tpole           = 4800.0,
    ldtype          = 3,
    ld1             = 0.23,
    ld2             = 0.0,
    inclination     = 78.0,
    position_angle  = 24.0,
    rotation_period = 54.8,
    beta            = 0.08,
    frac_escapevel  = 0.0,
    B_rot           = 0.0
)
test_gradient("Sphere", star_params_sphere,
              [1.37, 78.0, 24.0], ["radius", "inc", "PA"])

# ─── Test 2: Ellipsoid ──────────────────────────────────────────────────────
star_params_ellipsoid = (
    surface_type    = 1,
    radius_x        = 1.5,
    radius_y        = 1.3,
    radius_z        = 1.1,
    tpole           = 4800.0,
    ldtype          = 3,
    ld1             = 0.23,
    ld2             = 0.0,
    inclination     = 78.0,
    position_angle  = 24.0,
    rotation_period = 54.8,
    beta            = 0.08,
    frac_escapevel  = 0.0,
    B_rot           = 0.0
)
test_gradient("Ellipsoid", star_params_ellipsoid,
              [1.5, 1.3, 1.1, 78.0, 24.0], ["rx", "ry", "rz", "inc", "PA"])

# ─── Test 3: Rapid Rotator ──────────────────────────────────────────────────
star_params_rr = (
    surface_type    = 2,
    rpole           = 1.37,
    tpole           = 4800.0,
    ldtype          = 3,
    ld1             = 0.23,
    ld2             = 0.0,
    inclination     = 78.0,
    position_angle  = 24.0,
    rotation_period = 54.8,
    beta            = 0.08,
    frac_escapevel  = 0.9,
    B_rot           = 0.0
)
test_gradient("Rapid Rotator (ω=0.9)", star_params_rr,
              [1.37, 0.9, 78.0, 24.0], ["rpole", "ω", "inc", "PA"])

# ─── Test 4: Rapid Rotator with low ω ───────────────────────────────────────
star_params_rr2 = (
    surface_type    = 2,
    rpole           = 1.37,
    tpole           = 4800.0,
    ldtype          = 3,
    ld1             = 0.23,
    ld2             = 0.0,
    inclination     = 60.0,
    position_angle  = 45.0,
    rotation_period = 54.8,
    beta            = 0.08,
    frac_escapevel  = 0.3,
    B_rot           = 0.0
)
test_gradient("Rapid Rotator (ω=0.3)", star_params_rr2,
              [1.37, 0.3, 60.0, 45.0], ["rpole", "ω", "inc", "PA"])

println("\n\nDone.")
