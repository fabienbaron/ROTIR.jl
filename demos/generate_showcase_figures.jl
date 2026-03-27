# Generate showcase figures for ROTIR documentation and README.
#
# This script creates all PNG figures used in the docs and README.
# It is self-contained (no OIFITS data needed) and runs headless.
#
# Usage:
#   julia demos/generate_showcase_figures.jl
#
# Output: docs/src/assets/*.png (15 files)

# Force non-interactive backend BEFORE loading PyPlot (via ROTIR)
ENV["MPLBACKEND"] = "agg"

using ROTIR
using PyPlot

const OUTDIR = joinpath(@__DIR__, "..", "docs", "src", "assets")
mkpath(OUTDIR)
const DPI = 150

# --- Helpers ---

"""Save a figure returned as (fig, ax) by plot2d / plot2d_binary."""
function save_and_close(fig, filename)
    path = joinpath(OUTDIR, filename)
    fig.savefig(path, dpi=DPI, bbox_inches="tight", facecolor="white")
    close(fig)
    println("  Saved $filename")
end

"""Save the current figure (for void-returning plot functions like plot2d_wireframe, plot_mollweide)."""
function save_current(filename)
    fig = gcf()
    save_and_close(fig, filename)
end

# ===========================================================================
# Section 1: Star parameter definitions
# ===========================================================================

sphere_params = (
    surface_type    = 0,
    radius          = 1.0,      # mas
    tpole           = 5800.0,   # K (solar-like)
    ldtype          = 3,        # Hestroffer power law
    ld1             = 0.3,
    ld2             = 0.0,
    inclination     = 60.0,     # degrees
    position_angle  = 30.0,     # degrees
    rotation_period = 1.0,      # days
)

ellipsoid_params = (
    surface_type    = 1,
    radius_x        = 1.5,      # mas
    radius_y        = 1.3,
    radius_z        = 1.1,
    tpole           = 5000.0,   # K
    ldtype          = 3,
    ld1             = 0.23,
    ld2             = 0.0,
    inclination     = 60.0,
    position_angle  = 30.0,
    rotation_period = 10.0,
    beta            = 0.08,
)

rotator_params = (
    surface_type    = 2,
    rpole           = 1.37,     # mas
    tpole           = 4800.0,   # K
    ldtype          = 3,
    ld1             = 0.23,
    ld2             = 0.0,
    inclination     = 70.0,
    position_angle  = 20.0,
    rotation_period = 54.8,     # days
    beta            = 0.08,
    frac_escapevel  = 0.9,
    B_rot           = 0.0,
)

roche_params = (
    surface_type    = 3,
    rpole           = 0.355,    # mas
    tpole           = 4800.0,   # K
    ldtype          = 3,
    ld1             = 0.23,
    ld2             = 0.0,
    inclination     = 60.0,
    position_angle  = 0.0,
    rotation_period = 5.0,      # days
    beta            = 0.08,
    d               = 77.0,     # parsecs
    q               = 1.0,      # mass ratio
    fillout_factor_primary = -1,
    i               = 0.0,
    Ω               = 0.0,
    ω               = 0.0,
    P               = 5.0,      # days
    a               = 1.0,      # mas
    e               = 0.0,
    T0              = 0.0,
    dP              = 0.0,
    dω              = 0.0,
)

# ===========================================================================
# Section 2: Surface type gallery (4 PNGs)
# ===========================================================================
println("Generating surface type gallery...")

for (params, name, show_spin, hp_n) in [
    (sphere_params,    "surface_sphere.png",          false, 4),
    (ellipsoid_params, "surface_ellipsoid.png",       false, 4),
    (rotator_params,   "surface_rapid_rotator.png",   true,  4),
    (roche_params,     "surface_roche.png",           false, 5),  # higher res for L1 region
]
    local tessels = tessellation_healpix(hp_n)
    local star = create_star(tessels, params, 0.0)
    local tmap = parametric_temperature_map(params, star)
    local fig, ax = plot2d(tmap, star;
        intensity      = true,
        graticules     = true,
        compass        = true,
        rotation_axis  = show_spin,
        rotation_arrow = show_spin,
        inclination    = params.inclination,
        position_angle = params.position_angle,
    )
    save_and_close(fig, name)
end

# ===========================================================================
# Section 3: Tessellation comparison (3 PNGs)
# ===========================================================================
println("Generating tessellation comparison figures...")

# HEALPix wireframe
tessels_hp = tessellation_healpix(3)
star_hp = create_star(tessels_hp, sphere_params, 0.0)
plot2d_wireframe(star_hp; compass=true)
save_current("tess_healpix_wireframe.png")

# Lon/lat wireframe
tessels_ll = tessellation_latlong(20, 40)
star_ll = create_star(tessels_ll, sphere_params, 0.0)
plot2d_wireframe(star_ll; compass=true)
save_current("tess_latlong_wireframe.png")

# HEALPix with mesh overlay (temperature + edges visible)
tmap_hp = parametric_temperature_map(sphere_params, star_hp)
fig, ax = plot2d(tmap_hp, star_hp;
    intensity = true,
    plotmesh  = true,
    compass   = true,
)
save_and_close(fig, "tess_healpix_mesh.png")

# Lon/lat with mesh overlay (temperature + edges visible)
tmap_ll = parametric_temperature_map(sphere_params, star_ll)
fig, ax = plot2d(tmap_ll, star_ll;
    intensity = true,
    plotmesh  = true,
    compass   = true,
)
save_and_close(fig, "tess_latlong_mesh.png")

# Lon/lat with a temperature spot (unique lon/lat feature)
ntheta = 20; nphi = 40
tmap_spot = copy(tmap_ll)
tmap_spot = make_circ_spot(tmap_spot, ntheta, nphi, 4, 10, 20; bright_frac=0.5)
fig, ax = plot2d(tmap_spot, star_ll;
    intensity  = true,
    compass    = true,
    graticules = true,
    inclination    = sphere_params.inclination,
    position_angle = sphere_params.position_angle,
)
save_and_close(fig, "tess_latlong_spot.png")

# ===========================================================================
# Section 4: HEALPix resolution progression (3 PNGs)
# ===========================================================================
println("Generating HEALPix resolution progression...")

for n in [2, 3, 4]
    local tessels = tessellation_healpix(n)
    local star = create_star(tessels, rotator_params, 0.0)
    local tmap = parametric_temperature_map(rotator_params, star)
    local fig, ax = plot2d(tmap, star;
        intensity      = true,
        plotmesh       = true,
        compass        = true,
        inclination    = rotator_params.inclination,
        position_angle = rotator_params.position_angle,
    )
    save_and_close(fig, "healpix_n$n.png")
end

# ===========================================================================
# Section 5: Plotting options showcase (4 PNGs)
# ===========================================================================
println("Generating plotting options showcase...")

tessels_hires = tessellation_healpix(4)
star_show = create_star(tessels_hires, rotator_params, 0.0)
tmap_show = parametric_temperature_map(rotator_params, star_show)

# Plain (no decorations)
fig, ax = plot2d(tmap_show, star_show;
    intensity      = false,
    compass        = false,
    graticules     = false,
    rotation_axis  = false,
    rotation_arrow = false,
)
save_and_close(fig, "plot_plain.png")

# Fully decorated
fig, ax = plot2d(tmap_show, star_show;
    intensity      = true,
    compass        = true,
    graticules     = true,
    rotation_axis  = true,
    rotation_arrow = true,
    inclination    = rotator_params.inclination,
    position_angle = rotator_params.position_angle,
)
save_and_close(fig, "plot_decorated.png")

# Intensity (limb-darkened, compass only)
fig, ax = plot2d(tmap_show, star_show;
    intensity = true,
    compass   = true,
)
save_and_close(fig, "plot_intensity.png")

# Mollweide projection
plot_mollweide(tmap_show, star_show;
    colormap = "gist_heat",
    incl     = rotator_params.inclination,
    figtitle = "Mollweide projection",
)
save_current("plot_mollweide.png")

# ===========================================================================
# Section 6: Conventions annotated figure (1 PNG)
# ===========================================================================
println("Generating conventions annotated figure...")

conv_params = (
    surface_type    = 0,
    radius          = 1.0,
    tpole           = 5800.0,
    ldtype          = 3,
    ld1             = 0.3,
    ld2             = 0.0,
    inclination     = 45.0,
    position_angle  = 45.0,
    rotation_period = 1.0,
)
tessels_conv = tessellation_healpix(4)
star_conv = create_star(tessels_conv, conv_params, 0.0)
tmap_conv = parametric_temperature_map(conv_params, star_conv)

fig, ax = plot2d(tmap_conv, star_conv;
    intensity      = true,
    graticules     = true,
    compass        = true,
    rotation_axis  = true,
    rotation_arrow = true,
    inclination    = 45.0,
    position_angle = 45.0,
    figtitle       = "inc = 45°, PA = 45°",
)
save_and_close(fig, "conventions_annotated.png")

# ===========================================================================
# Section 7: Rapid rotator oblateness progression (4 PNGs)
# ===========================================================================
println("Generating rapid rotator oblateness progression...")

for (omega, label) in [(0.0, "00"), (0.5, "50"), (0.9, "90"), (0.99, "99")]
    local params = (
        surface_type    = 2,
        rpole           = 1.37,
        tpole           = 4800.0,
        ldtype          = 3,
        ld1             = 0.23,
        ld2             = 0.0,
        inclination     = 70.0,
        position_angle  = 20.0,
        rotation_period = 54.8,
        beta            = 0.08,
        frac_escapevel  = omega,
        B_rot           = 0.0,
    )
    local tessels = tessellation_healpix(4)
    local star = create_star(tessels, params, 0.0)
    local tmap = parametric_temperature_map(params, star)
    local fig, ax = plot2d(tmap, star;
        intensity      = true,
        plotmesh       = true,
        compass        = true,
        inclination    = params.inclination,
        position_angle = params.position_angle,
    )
    save_and_close(fig, "rotator_omega$label.png")
end

# ===========================================================================
# Section 8: Roche lobe fillout progression (4 PNGs)
# ===========================================================================
println("Generating Roche lobe fillout progression...")

for (fillout, label) in [(0.90, "90"), (0.95, "95"), (0.98, "98"), (0.99, "99")]
    local params = (
        surface_type    = 3,
        rpole           = 0.355,              # overridden by fillout_factor_primary
        tpole           = 4800.0,
        ldtype          = 3,
        ld1             = 0.23,
        ld2             = 0.0,
        inclination     = 60.0,
        position_angle  = 0.0,
        rotation_period = 5.0,
        beta            = 0.08,
        d               = 77.0,
        q               = 1.0,
        fillout_factor_primary = fillout,     # controls lobe size
        i               = 0.0,
        Ω               = 0.0,
        ω               = 0.0,
        P               = 5.0,
        a               = 1.0,
        e               = 0.0,
        T0              = 0.0,
        dP              = 0.0,
        dω              = 0.0,
    )
    local tessels = tessellation_healpix(5)
    local star = create_star(tessels, params, 0.0)
    local tmap = parametric_temperature_map(params, star)
    local fig, ax = plot2d(tmap, star;
        intensity      = true,
        plotmesh       = true,
        compass        = true,
        inclination    = params.inclination,
        position_angle = params.position_angle,
    )
    save_and_close(fig, "roche_fill$label.png")
end

# ===========================================================================
println("\nDone! Generated $(length(readdir(OUTDIR))) files in $OUTDIR")
