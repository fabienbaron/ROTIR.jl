using ROTIR
using Statistics
using LinearAlgebra
using Test
using PythonPlot

# ── Inline the two functions under test (di.jl has heavy deps) ────────
function make_circ_spot(temperature_map, star_geometry, spot_radius, lat, long; bright_frac=0.8)
    long = mod(long, 360)
    temperature_map_copy = deepcopy(temperature_map)
    centers = star_geometry.vertices_spherical[:,5,:]
    r_c = @view centers[:, 1]
    θ_c = @view centers[:, 2]
    φ_c = @view centers[:, 3]
    sinθ = sin.(θ_c); cosθ = cos.(θ_c)
    px = r_c .* sinθ .* cos.(φ_c)
    py = r_c .* sinθ .* sin.(φ_c)
    pz = r_c .* cosθ
    spot_colat = (90 - lat) * pi / 180
    spot_lon   = long * pi / 180
    ux = sin(spot_colat) * cos(spot_lon)
    uy = sin(spot_colat) * sin(spot_lon)
    uz = cos(spot_colat)
    r_spot = r_c[argmax(sinθ .* ux .* cos.(φ_c) .+ sinθ .* uy .* sin.(φ_c) .+ cosθ .* uz)]
    d = sqrt.((px .- r_spot * ux).^2 .+ (py .- r_spot * uy).^2 .+ (pz .- r_spot * uz).^2)
    chord_radius = 2 * r_spot * sin(spot_radius * pi / 360)
    spot_mask = findall(d .<= chord_radius)
    temperature_map_copy[spot_mask] .= median(temperature_map_copy) * bright_frac
    return vec(temperature_map_copy)
end

function make_circ_spot_spotfill(star_geometry, spot_radius, lat, long; profile="flat")
    long = mod(long, 360)
    x = zeros(star_geometry.npix)
    centers = star_geometry.vertices_spherical[:,5,:]
    r_c = @view centers[:, 1]
    θ_c = @view centers[:, 2]
    φ_c = @view centers[:, 3]
    sinθ = sin.(θ_c); cosθ = cos.(θ_c)
    px = r_c .* sinθ .* cos.(φ_c)
    py = r_c .* sinθ .* sin.(φ_c)
    pz = r_c .* cosθ
    spot_colat = (90 - lat) * pi / 180
    spot_lon   = long * pi / 180
    ux = sin(spot_colat) * cos(spot_lon)
    uy = sin(spot_colat) * sin(spot_lon)
    uz = cos(spot_colat)
    r_spot = r_c[argmax(sinθ .* ux .* cos.(φ_c) .+ sinθ .* uy .* sin.(φ_c) .+ cosθ .* uz)]
    d = sqrt.((px .- r_spot * ux).^2 .+ (py .- r_spot * uy).^2 .+ (pz .- r_spot * uz).^2)
    chord_radius = 2 * r_spot * sin(spot_radius * pi / 360)
    spot_mask = findall(d .<= chord_radius)
    if profile == "flat"
        x[spot_mask] .= 1.0
    elseif profile == "linear"
        x[spot_mask] .= 1.0 .- d[spot_mask] ./ chord_radius
    end
    return vec(x)
end

# Old haversine method for comparison
function make_circ_spot_haversine(temperature_map, star_geometry, spot_radius, lat, long; bright_frac=0.8)
    long = mod(long, 360)
    temperature_map_copy = deepcopy(temperature_map)
    centers = star_geometry.vertices_spherical[:,5,:]
    θ = centers[:, 3]
    ϕ = -centers[:, 2] .+ pi/2
    da = 2.0*asin.(sqrt.(sin.(0.5*(ϕ .- lat*pi/180)).^2 .+
         cos(lat*pi/180).*cos.(ϕ).*sin.(0.5*(θ .- long*pi/180)).^2))
    spot_mask = findall(da .<= spot_radius*pi/180)
    temperature_map_copy[spot_mask] .= median(temperature_map_copy) * bright_frac
    return vec(temperature_map_copy)
end

# ── Helpers ───────────────────────────────────────────────────────────
spotted_pixels(tmap, tmap_spot) = findall(tmap_spot .!= tmap)

n = 5  # HEALPix nside → 3072 patches
tessels = tessellation_healpix(n)

# ======================================================================
# Test 1: On a sphere, new Euclidean method gives same mask as haversine
# ======================================================================
@testset "Sphere: Euclidean ≡ haversine" begin
    params = (
        surface_type    = 0,
        radius          = 1.0,
        tpole           = 10000.0,
        ldtype          = 3,
        ld1             = 0.3,
        ld2             = 0.0,
        inclination     = 60.0,
        position_angle  = 0.0,
        rotation_period = 1.0,
    )
    star = create_star(tessels, params, 0.0)
    tmap = parametric_temperature_map(params, star)

    for (lat, long, sr) in [(0, 0, 15), (45, 90, 20), (-30, 270, 10), (80, 180, 25)]
        new_mask = Set(spotted_pixels(tmap, make_circ_spot(tmap, star, sr, lat, long)))
        old_mask = Set(spotted_pixels(tmap, make_circ_spot_haversine(tmap, star, sr, lat, long)))
        @test new_mask == old_mask
    end
end

# ======================================================================
# Test 2: On a triaxial ellipsoid, Euclidean method gives uniform
#          boundary distances while haversine does not
# ======================================================================
@testset "Ellipsoid: boundary distance uniformity" begin
    params = (
        surface_type    = 1,
        radius_x        = 2.0,
        radius_y        = 1.5,
        radius_z        = 1.0,
        tpole           = 10000.0,
        ldtype          = 3,
        ld1             = 0.3,
        ld2             = 0.0,
        inclination     = 90.0,
        position_angle  = 0.0,
        rotation_period = 1.0,
        beta            = 0.08,
    )
    star = create_star(tessels, params, 0.0)
    tmap = parametric_temperature_map(params, star)

    lat, lon, sr = 0, 45, 20
    long_m = mod(lon, 360)
    tmap_new = make_circ_spot(tmap, star, sr, lat, lon)
    tmap_old = make_circ_spot_haversine(tmap, star, sr, lat, lon)
    mask_new = spotted_pixels(tmap, tmap_new)
    mask_old = spotted_pixels(tmap, tmap_old)

    # Compute Euclidean distances from spot center for both masks
    centers = star.vertices_spherical[:, 5, :]
    r_c = centers[:, 1]; θ_c = centers[:, 2]; φ_c = centers[:, 3]
    sinθ = sin.(θ_c); cosθ = cos.(θ_c)
    px = r_c .* sinθ .* cos.(φ_c)
    py = r_c .* sinθ .* sin.(φ_c)
    pz = r_c .* cosθ

    spot_colat = (90 - lat) * pi / 180
    spot_lon = long_m * pi / 180
    ux = sin(spot_colat) * cos(spot_lon)
    uy = sin(spot_colat) * sin(spot_lon)
    uz = cos(spot_colat)
    r_spot = r_c[argmax(sinθ .* ux .* cos.(φ_c) .+ sinθ .* uy .* sin.(φ_c) .+ cosθ .* uz)]
    cx, cy, cz = r_spot * ux, r_spot * uy, r_spot * uz
    euclid_dist = sqrt.((px .- cx).^2 .+ (py .- cy).^2 .+ (pz .- cz).^2)

    chord_radius = 2 * r_spot * sin(sr * pi / 360)

    # For the Euclidean method: all spotted patches must be within chord_radius
    @test all(euclid_dist[mask_new] .<= chord_radius + 1e-6)

    # Haversine method should include patches beyond chord_radius
    max_euclid_old = maximum(euclid_dist[mask_old])
    max_euclid_new = maximum(euclid_dist[mask_new])
    println("  Max Euclidean dist — Euclidean method: $(round(max_euclid_new, digits=4)),  chord_radius: $(round(chord_radius, digits=4))")
    println("  Max Euclidean dist — Haversine method: $(round(max_euclid_old, digits=4))")
    @test max_euclid_new <= chord_radius + 1e-6
    @test max_euclid_old > chord_radius

    # The two methods should select different pixels on a non-spherical surface
    @test Set(mask_new) != Set(mask_old)
end

# ======================================================================
# Test 3: Polar spot on oblate star — sphere & ellipsoid agree
# ======================================================================
@testset "Pole spot: sphere vs oblate ellipsoid" begin
    params_sphere = (
        surface_type = 0, radius = 1.0, tpole = 10000.0,
        ldtype = 3, ld1 = 0.3, ld2 = 0.0,
        inclination = 90.0, position_angle = 0.0, rotation_period = 1.0,
    )
    params_oblate = (
        surface_type = 1, radius_x = 1.0, radius_y = 1.0, radius_z = 0.7,
        tpole = 10000.0, ldtype = 3, ld1 = 0.3, ld2 = 0.0,
        inclination = 90.0, position_angle = 0.0, rotation_period = 1.0,
        beta = 0.08,
    )
    star_s = create_star(tessels, params_sphere, 0.0)
    star_e = create_star(tessels, params_oblate, 0.0)
    tmap_s = parametric_temperature_map(params_sphere, star_s)
    tmap_e = parametric_temperature_map(params_oblate, star_e)

    mask_s = Set(spotted_pixels(tmap_s, make_circ_spot(tmap_s, star_s, 10, 90, 0)))
    mask_e = Set(spotted_pixels(tmap_e, make_circ_spot(tmap_e, star_e, 10, 90, 0)))
    @test mask_s == mask_e
end

# ======================================================================
# Test 4: spotfill flat profile matches make_circ_spot mask
# ======================================================================
@testset "spotfill flat ≡ make_circ_spot mask" begin
    params = (
        surface_type = 1, radius_x = 1.8, radius_y = 1.2, radius_z = 1.0,
        tpole = 10000.0, ldtype = 3, ld1 = 0.3, ld2 = 0.0,
        inclination = 60.0, position_angle = 0.0, rotation_period = 1.0,
        beta = 0.08,
    )
    star = create_star(tessels, params, 0.0)
    tmap = parametric_temperature_map(params, star)

    lat, long, sr = 30, 45, 18
    mask_spot = Set(spotted_pixels(tmap, make_circ_spot(tmap, star, sr, lat, long)))
    mask_fill = Set(findall(make_circ_spot_spotfill(star, sr, lat, long; profile="flat") .> 0))
    @test mask_spot == mask_fill
end

# ======================================================================
# Test 5: spotfill linear profile shape
# ======================================================================
@testset "spotfill linear profile" begin
    params = (
        surface_type = 0, radius = 1.0, tpole = 10000.0,
        ldtype = 3, ld1 = 0.3, ld2 = 0.0,
        inclination = 90.0, position_angle = 0.0, rotation_period = 1.0,
    )
    star = create_star(tessels, params, 0.0)
    fill_map = make_circ_spot_spotfill(star, 20, 0, 0; profile="linear")

    inside = findall(fill_map .> 0)
    @test length(inside) > 0
    @test maximum(fill_map) ≈ 1.0 atol=0.1    # center pixel ≈ 1 (limited by patch resolution)
    @test all(fill_map .>= 0)                  # no negatives
    @test all(fill_map .<= 1.0 + eps())        # bounded above
end

# ======================================================================
# Test 6: spot pixel count scales ~ radius²
# ======================================================================
@testset "Spot area ∝ radius²" begin
    params = (
        surface_type = 0, radius = 1.0, tpole = 10000.0,
        ldtype = 3, ld1 = 0.3, ld2 = 0.0,
        inclination = 90.0, position_angle = 0.0, rotation_period = 1.0,
    )
    star = create_star(tessels, params, 0.0)
    tmap = parametric_temperature_map(params, star)

    n10 = length(spotted_pixels(tmap, make_circ_spot(tmap, star, 10, 0, 0)))
    n20 = length(spotted_pixels(tmap, make_circ_spot(tmap, star, 20, 0, 0)))
    ratio = n20 / n10
    @test 2.5 < ratio < 5.5
    println("  Pixel count ratio (r=20/r=10): $(round(ratio, digits=2))  (expect ~4)")
end

println("\nAll spot tests passed.")

# ======================================================================
# Visual comparison: Euclidean vs Haversine on a triaxial ellipsoid
# ======================================================================
outdir = joinpath(@__DIR__, "figures")
mkpath(outdir)

params_ellipsoid = (
    surface_type = 1, radius_x = 2.0, radius_y = 1.5, radius_z = 1.0,
    tpole = 10000.0, ldtype = 3, ld1 = 0.3, ld2 = 0.0,
    inclination = 70.0, position_angle = 0.0, rotation_period = 1.0,
    beta = 0.08,
)

plot_kwargs = (intensity=true, graticules=true,
               inclination=params_ellipsoid.inclination,
               position_angle=params_ellipsoid.position_angle,
               star_params=params_ellipsoid)

# Sub-observer longitude at t=0, PA=0 is ~270° (from rotation convention),
# so place spot at long=270 to face the observer head-on.
spot_long = 270
spot_lat  = 0
spot_rad  = 25

for (phase, label) in [(0.0, "front"), (0.25, "side")]
    star = create_star(tessels, params_ellipsoid, phase)
    tmap = parametric_temperature_map(params_ellipsoid, star)
    tmap_euclid    = make_circ_spot(tmap, star, spot_rad, spot_lat, spot_long; bright_frac=0.5)
    tmap_haversine = make_circ_spot_haversine(tmap, star, spot_rad, spot_lat, spot_long; bright_frac=0.5)

    fig, ax = plot2d(tmap_euclid, star; plot_kwargs..., figtitle="Euclidean (new) — phase=$phase")
    savefig(joinpath(outdir, "spot_euclid_$(label).png"), dpi=150, bbox_inches="tight")
    close(fig)

    fig, ax = plot2d(tmap_haversine, star; plot_kwargs..., figtitle="Haversine (old) — phase=$phase")
    savefig(joinpath(outdir, "spot_haversine_$(label).png"), dpi=150, bbox_inches="tight")
    close(fig)

    println("Saved: spot_euclid_$(label).png, spot_haversine_$(label).png")
end

# Sphere control
params_sphere = (
    surface_type = 0, radius = 1.0, tpole = 10000.0,
    ldtype = 3, ld1 = 0.3, ld2 = 0.0,
    inclination = 70.0, position_angle = 0.0, rotation_period = 1.0,
)
star = create_star(tessels, params_sphere, 0.0)
tmap = parametric_temperature_map(params_sphere, star)

fig, ax = plot2d(make_circ_spot(tmap, star, 25, 0, 270; bright_frac=0.5), star;
    intensity=true, graticules=true,
    inclination=params_sphere.inclination, position_angle=params_sphere.position_angle,
    star_params=params_sphere, figtitle="Sphere (control) — Euclidean")
savefig(joinpath(outdir, "spot_sphere_control.png"), dpi=150, bbox_inches="tight")
close(fig)
println("Saved: spot_sphere_control.png")
