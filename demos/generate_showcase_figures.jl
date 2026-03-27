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
using PyPlot, PyCall

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
tmap_spot = make_circ_spot(tmap_spot, ntheta, nphi, 4, 8, 30; bright_frac=0.3)
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
        compass        = true,
        inclination    = params.inclination,
        position_angle = params.position_angle,
    )
    ax.set_xlim(0.75, -0.75)
    ax.set_ylim(-0.75, 0.75)
    save_and_close(fig, "roche_fill$label.png")
end

# ===========================================================================
# Section 9: Binary orbit showcase (Spica-like parameters, 3 PNGs)
# ===========================================================================
println("Generating binary orbit figures...")

# Spica orbital elements (Aufdenberg+2015)
P_orb    = 4.0145       # days
a_orb    = 1.54         # mas
e_orb    = 0.123
T0_orb   = 2454189.40   # JD
q_binary = 0.6188       # M2/M1
i_orb    = 116.0        # degrees (>90 = retrograde)
Omega_orb = 309.938     # degrees
omega_orb = 255.0       # degrees

inc_star = 180.0 - i_orb      # prograde equivalent = 64°
pa_star  = Omega_orb - 180.0  # spin PA on sky = 129.938°

star1p = starparameters(0.93/2, 25300.0, 0.0, 3, 0.15, 0.0, 0.25, 0.0, inc_star, pa_star, 0.0, P_orb)
star2p = starparameters(0.57/2, 20585.0, 0.0, 3, 0.15, 0.0, 0.25, 0.0, inc_star, pa_star, 0.0, P_orb)
bparams = binaryparameters(star1p, star2p, 77.0, i_orb, Omega_orb, omega_orb,
    P_orb, a_orb, e_orb, T0_orb, q_binary, [1.0, 1.0], 0.0, 0.0)

# --- Binary sky-plane image at widest separation ---
star1_nt = (surface_type=0, radius=0.93/2, tpole=25300.0, ldtype=3, ld1=0.15, ld2=0.0,
    inclination=inc_star, position_angle=pa_star, rotation_period=P_orb)
star2_nt = (surface_type=0, radius=0.57/2, tpole=20585.0, ldtype=3, ld1=0.15, ld2=0.0,
    inclination=inc_star, position_angle=pa_star, rotation_period=P_orb)

tessels_b1 = tessellation_healpix(3)
tessels_b2 = tessellation_healpix(2)
star1_geom = create_star(tessels_b1, star1_nt, 0.0)
star2_geom = create_star(tessels_b2, star2_nt, 0.0)
tmap_b1 = parametric_temperature_map(star1_nt, star1_geom)
tmap_b2 = parametric_temperature_map(star2_nt, star2_geom)

# Pick an epoch near widest separation (phase ~0.25 after periastron)
t_wide = T0_orb + 0.25 * P_orb
fig, ax = plot2d_binary(tmap_b1, tmap_b2, star1_geom, star2_geom, bparams, t_wide;
    intensity=true, graticules=true, compass=true,
    inclination1=inc_star, position_angle1=pa_star,
    inclination2=inc_star, position_angle2=pa_star)
save_and_close(fig, "binary_skyplane.png")

# --- Orbital diagram ---
patches_mpl = pyimport("matplotlib.patches")
fig2, ax2 = subplots(1, 1, figsize=(8, 8))
ax2.set_aspect("equal", adjustable="box")

tepochs_orbit = bparams.T0 .+ collect(range(0.0, stop=bparams.P, length=500))
orbit_x = zeros(500)
orbit_y = zeros(500)
for (j, t) in enumerate(tepochs_orbit)
    local x1, y1, z1, x2, y2, z2 = binary_orbit_abs(bparams, t)
    orbit_x[j] = y2 - y1
    orbit_y[j] = x2 - x1
end
ax2.plot(orbit_x, orbit_y, "b-", linewidth=1.5, alpha=0.7)

rpri = bparams.star1.rpole
rsec = bparams.star2.rpole
c_pri = patches_mpl.Circle((0, 0), rpri, facecolor="gold", edgecolor="black", linewidth=0.5, zorder=5)
ax2.add_patch(c_pri)

# Mark a few orbital phases
for (phase, label) in [(0.0, "periastron"), (0.25, ""), (0.5, "apastron"), (0.75, "")]
    local t_jd = T0_orb + phase * P_orb
    local x1, y1, z1, x2, y2, z2 = binary_orbit_abs(bparams, t_jd)
    local east = y2 - y1; local north = x2 - x1
    local c_sec = patches_mpl.Circle((east, north), rsec, facecolor="lightskyblue",
        edgecolor="black", linewidth=0.5, zorder=4, alpha=0.8)
    ax2.add_patch(c_sec)
    if label != ""
        ax2.annotate(label, (east, north), textcoords="offset points", xytext=(8, 8), fontsize=9)
    end
end

ax2.invert_xaxis()
ax2.set_xlabel("East offset (mas)", fontsize=14)
ax2.set_ylabel("North offset (mas)", fontsize=14)
ax2.set_title("Binary orbit (secondary relative to primary)")
ax2.grid(true, alpha=0.3)
save_and_close(fig2, "binary_orbit.png")

# --- Radial velocity curves ---
fig3, ax3 = plot_rv(bparams, K1=123.9, K2=198.8, γ=0.0,
    figtitle="Radial velocities")
save_and_close(fig3, "binary_rv.png")

# ===========================================================================
println("\nDone! Generated $(length(readdir(OUTDIR))) files in $OUTDIR")
