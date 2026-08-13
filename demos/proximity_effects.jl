#!/usr/bin/env julia
# =======================================================================================
# Proximity effects in a close binary: Roche distortion, gravity darkening, irradiation
# =======================================================================================
# A 3 x N panel figure that separates the three effects that only appear when two stars get
# close, by turning them off one at a time while holding everything else fixed.
#
#   COLUMNS  decreasing separation. The stellar polar radii are held constant, so shrinking
#            the orbit is what drives the components deeper into their Roche lobes. The
#            leftmost pair is nearly spherical; the rightmost is strongly distorted.
#
#   ROW 1    everything on:  Roche geometry + gravity darkening + mutual irradiation
#   ROW 2    irradiation OFF (albedo = 0), same geometry and gravity darkening
#   ROW 3    gravity darkening OFF (beta = 0), irradiation back on
#
# Reading it: row 1 minus row 2 isolates irradiation; row 1 minus row 3 isolates gravity
# darkening. They pull in OPPOSITE directions at the same place, which is the whole point —
# the sub-companion point is the tidal bulge tip, so it has the LOWEST local gravity (von
# Zeipel makes it cool) while also receiving the most flux from the companion (irradiation
# makes it hot). Comparing rows is the only way to see which wins.
#
# One shared colour scale across all panels, computed in a first pass. Per-panel
# normalisation would make every panel look the same and destroy the comparison.
#
#   julia --project=demos demos/proximity_effects.jl
#   NSIDE=4 SEPS=3.0,2.0,1.5,1.25 OUT=/tmp/prox.png julia --project=demos demos/proximity_effects.jl

ENV["MPLBACKEND"] = get(ENV, "MPLBACKEND", "Agg")
using ROTIR, PythonPlot, PythonCall, Printf, Statistics

const NSIDE = parse(Int, get(ENV, "NSIDE", "4"))
const OUT   = get(ENV, "OUT", joinpath(@__DIR__, "proximity_effects.png"))
const SEPS  = parse.(Float64, split(get(ENV, "SEPS", "2.60,1.80,1.40,1.20"), ","))
const ALBEDO = parse(Float64, get(ENV, "ALBEDO", "0.6"))
const BETA   = parse(Float64, get(ENV, "BETA", "0.25"))   # von Zeipel, radiative envelope

# --- system -----------------------------------------------------------------------------
# A generic hot detached pair. Radii are POLAR radii in mas and stay fixed across columns;
# only the orbit shrinks.
const RPOLE1, RPOLE2 = 0.45, 0.30      # mas
# Close in temperature ON PURPOSE. With a shared colour scale a 6000 K gap between the
# components eats the whole range and each star's own structure — which is the entire
# subject of this figure — is squeezed into a couple of colour steps.
const TPOLE1, TPOLE2 = 20000.0, 17500.0
const QBIN  = 0.60                     # M2/M1
const PORB  = 4.0                      # d (only sets the tidal locking period)
const INC   = 75.0                     # deg
const OMEGA = 30.0                     # deg
const DIST  = 100.0                    # pc

# NB: a Roche component's params must ALSO carry the orbital elements. `create_binary_star`
# reads `p.a` to reduce the instantaneous separation to units of the semi-major axis, so a
# params NamedTuple without them fails with a bare `getproperty` error.
roche_params(rpole, tpole, q, beta, a) = (
    surface_type              = 3,          # Roche
    rpole                     = rpole,
    tpole                     = tpole,
    ldtype                    = 3,          # Hestroffer power law
    ld1                       = 0.30,
    ld2                       = 0.0,
    inclination               = INC,
    position_angle            = OMEGA,
    rotation_period           = PORB,       # tidally locked
    beta                      = beta,       # 0 turns gravity darkening off
    d                         = DIST,
    q                         = q,
    fillout_factor_primary    = -1,         # rpole defines the potential
    fillout_factor_secondary  = -1,
    i = INC, Ω = OMEGA, ω = 90.0, P = PORB, a = a, e = 0.0, T0 = 0.0,
    dP = 0.0, dω = 0.0,
)

"""
Build both components and their temperature maps for one column (separation `a`) and one
row configuration.

`beta = 0` removes gravity darkening: `T = tpole * (g/g_pole)^beta` becomes uniform, while
the SHAPE stays Roche-distorted. That separation matters — it isolates the temperature
effect from the geometric one.
"""
function make_pair(a_mas; beta = BETA, irradiate = true, nside = NSIDE)
    p1 = roche_params(RPOLE1, TPOLE1, QBIN,     beta, a_mas)
    p2 = roche_params(RPOLE2, TPOLE2, 1 / QBIN, beta, a_mas)   # q is M_companion/M_self
    s1p = starparameters(RPOLE1, TPOLE1, 0.0, 3, 0.30, 0.0, beta, 0.0, INC, OMEGA, 0.0, PORB)
    s2p = starparameters(RPOLE2, TPOLE2, 0.0, 3, 0.30, 0.0, beta, 0.0, INC, OMEGA, 0.0, PORB)
    bp  = binaryparameters(s1p, s2p, DIST, INC, OMEGA, 90.0, PORB, a_mas, 0.0, 0.0,
                           QBIN, [-1.0, -1.0], 0.0, 0.0)
    t1, t2 = tessellation_healpix(nside), tessellation_healpix(nside)
    star1, star2 = create_binary_geometry(t1, p1, t2, p2, bp, 0.0)
    # von Zeipel on the Roche surface: T = tpole*(g/g_pole)^beta, computed from the actual
    # distorted geometry, so beta = 0 gives a uniform map on a still-distorted star.
    m1 = temperature_map_vonZeipel_roche_single(p1, star1, 0.0)
    m2 = temperature_map_vonZeipel_roche_single(p2, star2, 0.0; secondary = true)
    if irradiate
        m1, m2 = handle_reflection(star1, m1, star2, m2;
                                   albedo1 = ALBEDO, albedo2 = ALBEDO)
    end
    return star1, star2, m1, m2, bp
end

const ROWS = [
    (label = "all effects on",           beta = BETA, irradiate = true),
    (label = "irradiation OFF",          beta = BETA, irradiate = false),
    (label = "gravity darkening OFF",    beta = 0.0,  irradiate = true),
]

# One frame for every panel: the columns are only comparable if they share it, and it must
# be tight or the stars are lost in white space.
const AXMAX = 0.62 * (maximum(SEPS) + RPOLE1 + RPOLE2)

println("Building $(length(ROWS)) x $(length(SEPS)) panels at nside = $NSIDE ...")

# --- pass 1: build everything and find the shared colour range ---------------------------
panels = Array{Any}(undef, length(ROWS), length(SEPS))
vmin, vmax = Inf, -Inf
for (r, row) in enumerate(ROWS), (c, a) in enumerate(SEPS)
    panels[r, c] = make_pair(a; beta = row.beta, irradiate = row.irradiate)
    _, _, m1, m2, _ = panels[r, c]
    global vmin = min(vmin, minimum(m1), minimum(m2))
    global vmax = max(vmax, maximum(m1), maximum(m2))
    @printf("  row %d (%-22s) a = %5.2f mas : T = %6.0f .. %6.0f K\n",
            r, row.label, a, min(minimum(m1), minimum(m2)), max(maximum(m1), maximum(m2)))
end
@printf("shared colour scale: %.0f .. %.0f K\n", vmin, vmax)

# --- pass 2: draw ------------------------------------------------------------------------
ROTIR.set_oiplot_defaults()
fig, axs = pyplot.subplots(length(ROWS), length(SEPS),
                           figsize = (3.6 * length(SEPS), 3.6 * length(ROWS)))
for (r, row) in enumerate(ROWS), (c, a) in enumerate(SEPS)
    star1, star2, m1, m2, bp = panels[r, c]
    ax = axs[r - 1][c - 1]          # 0-based: PythonCall hands back a numpy array of axes
    plot2d_binary(m1, m2, star1, star2, bp, 0.0;
                  ax = ax, colorbar_on = false, vmin = vmin, vmax = vmax,
                  axis_max = AXMAX,
                  compass = (r == 1 && c == 1), limb = true, plotmesh = false)
    # annotate what the physics actually did in this panel
    spread = max(maximum(m1) - minimum(m1), maximum(m2) - minimum(m2))
    ax.text(0.03, 0.03, @sprintf("ΔT = %.0f K", spread), transform = ax.transAxes,
            fontsize = 9, color = "0.25", ha = "left", va = "bottom")
    r == 1 && ax.set_title(@sprintf("a = %.2f mas   (r/a = %.2f)", a, RPOLE1/a), fontsize = 11)
    c == 1 && ax.set_ylabel(row.label, fontsize = 10)
    ax.set_xlabel(""); c == 1 || ax.set_ylabel("")
    ax.set_xticks([]); ax.set_yticks([])
end

# one colour bar for the whole figure — the panels are only comparable because they share it
cmap = pyimport("matplotlib").colormaps["gist_heat"]
norm = pyimport("matplotlib").colors.Normalize(vmin = vmin, vmax = vmax)
sm   = pyimport("matplotlib").cm.ScalarMappable(norm = norm, cmap = cmap)
cb   = fig.colorbar(sm, ax = axs, orientation = "vertical", fraction = 0.02, pad = 0.02)
cb.set_label("Temperature (K)")
fig.suptitle("Proximity effects: Roche distortion, gravity darkening, mutual irradiation",
             fontsize = 13)
fig.savefig(OUT, dpi = parse(Int, get(ENV, "DPI", "150")), bbox_inches = "tight",
            facecolor = "white")
pyplot.close("all")
println("wrote $OUT")
