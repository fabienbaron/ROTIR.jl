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

# What to colour. RELATIVE=1 (default) maps ΔT = T - tpole for each component on a SYMMETRIC
# diverging scale: blue is cooler than the unperturbed pole (gravity darkening), red hotter
# (irradiation), white unperturbed. That is the only way to see both effects at once — on an
# absolute scale the two components' different tpole eats the range, and "cooler" and
# "hotter" are not distinguishable from "this star is intrinsically colder".
const RELATIVE = get(ENV, "RELATIVE", "1") == "1"
const CMAP     = get(ENV, "CMAP", RELATIVE ? "RdBu_r" : "gist_heat")

# SIGNED POWER STRETCH. Peak |ΔT| runs from 167 K in the widest column to 3859 K in the
# closest — a 20x spread. On a shared LINEAR scale the strongest panel sets the ends of the
# colormap and the other eleven sit in a narrow band around white, so most of the figure
# carries almost no colour. Mapping sign(x)*|x|^GAMMA before colouring gives the small
# values far more of the colormap while keeping zero at zero (so white still means
# unperturbed) and keeping ONE scale for every panel, which is what makes them comparable.
# GAMMA = 1 is the plain linear scale.
const GAMMA = parse(Float64, get(ENV, "GAMMA", "0.45"))
stretch(x) = sign(x) * abs(x)^GAMMA

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

# VIEWING PHASE. The tidal bulge points along the line of CENTRES, so it foreshortens to
# nothing whenever that line tilts toward us — which is precisely at conjunction, where the
# discs also overlap. At t = T0 the projected separation here is only 26% of a; at
# quadrature it is 100%, i.e. 3.86x wider on the sky, with the line of centres lying in the
# sky plane so the elongation shows at full length.
#
# TEPOCH = "auto" (default) finds that phase by scanning one period. Give a number in days
# to pin it. With e = 0 the quadrature phase does not depend on the separation, so one value
# serves every column — and it MUST, or the columns stop being comparable.
const TEPOCH_ENV = get(ENV, "TEPOCH", "auto")

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
function make_pair(a_mas; beta = BETA, irradiate = true, nside = NSIDE, tepoch = 0.0)
    p1 = roche_params(RPOLE1, TPOLE1, QBIN,     beta, a_mas)
    p2 = roche_params(RPOLE2, TPOLE2, 1 / QBIN, beta, a_mas)   # q is M_companion/M_self
    s1p = starparameters(RPOLE1, TPOLE1, 0.0, 3, 0.30, 0.0, beta, 0.0, INC, OMEGA, 0.0, PORB)
    s2p = starparameters(RPOLE2, TPOLE2, 0.0, 3, 0.30, 0.0, beta, 0.0, INC, OMEGA, 0.0, PORB)
    bp  = binaryparameters(s1p, s2p, DIST, INC, OMEGA, 90.0, PORB, a_mas, 0.0, 0.0,
                           QBIN, [-1.0, -1.0], 0.0, 0.0)
    t1, t2 = tessellation_healpix(nside), tessellation_healpix(nside)
    star1, star2 = create_binary_geometry(t1, p1, t2, p2, bp, tepoch)
    # von Zeipel on the Roche surface: T = tpole*(g/g_pole)^beta, computed from the actual
    # distorted geometry, so beta = 0 gives a uniform map on a still-distorted star.
    m1 = temperature_map_vonZeipel_roche_single(p1, star1, tepoch)
    m2 = temperature_map_vonZeipel_roche_single(p2, star2, tepoch; secondary = true)
    if irradiate
        m1, m2 = handle_reflection(star1, m1, star2, m2;
                                   albedo1 = ALBEDO, albedo2 = ALBEDO)
    end
    return star1, star2, m1, m2, bp
end

"""Epoch of maximum projected separation — one scan over a period, e = 0 so it is phase-only."""
function quadrature_epoch(bp; n = 400)
    best_t, best_s = 0.0, -1.0
    for t in range(0, PORB, length = n)
        ow, on = orbit_to_rotir_offset(bp, t)
        s = hypot(ow, on)
        s > best_s && ((best_t, best_s) = (t, s))
    end
    return best_t, best_s
end

const ROWS = [
    (label = "all effects on",           beta = BETA, irradiate = true),
    (label = "irradiation OFF",          beta = BETA, irradiate = false),
    (label = "gravity darkening OFF",    beta = 0.0,  irradiate = true),
]

# One frame for every panel: the columns are only comparable if they share it, and it must
# be TIGHT or the stars are lost in white space. Sized to the widest column's pair, not to
# some round number.
# Filled in after pass 1 from the real projected extents: the sky-projected separation is
# a*cos-ish of the inclination and the orbital phase, NOT `a`, so sizing the frame from `a`
# leaves most of every panel empty.
AXMAX = 0.0

# GUARD. `fillout_factor = -1` derives the Roche potential from `rpole`, and nothing stops
# rpole exceeding the Roche lobe — the radii then diverge toward L1 and you get a garbage
# shape rather than an error. Measured for q = 0.6: at rpole/a = 0.375 the primary is
# distorted by a healthy 1.18, at 0.429 (just past the lobe) it blows up to 10.6.
const RLOBE = radius_eggleton(QBIN)
for a in SEPS
    RPOLE1 / a < 0.95 * RLOBE || error("""
        separation a = $a mas puts the primary at rpole/a = $(round(RPOLE1/a, digits=3)),
        at or past its Roche lobe (r_L/a = $(round(RLOBE, digits=3)) for q = $QBIN).
        The shape diverges toward L1 rather than failing, so this is refused. Use a larger
        separation, a smaller RPOLE1, or a mass ratio with a wider lobe.""")
end
@printf("mass ratio q = %.2f  ->  Roche lobe r_L/a = %.3f (Eggleton)\n", QBIN, RLOBE)
for a in SEPS
    @printf("   a = %5.2f mas : primary fills %3.0f%% of its lobe\n", a, 100*RPOLE1/(a*RLOBE))
end

let (_, _, _, _, bp0) = make_pair(maximum(SEPS); tepoch = 0.0)
    tq, sq = quadrature_epoch(bp0)
    global TEPOCH = TEPOCH_ENV == "auto" ? tq : parse(Float64, TEPOCH_ENV)
    ow, on = orbit_to_rotir_offset(bp0, TEPOCH)
    @printf("viewing epoch t = %.4f d  ->  projected separation %.0f%% of a  (i = %.0f deg)\n",
            TEPOCH, 100*hypot(ow, on)/maximum(SEPS), INC)
end

println("Building $(length(ROWS)) x $(length(SEPS)) panels at nside = $NSIDE ...")

# --- pass 1: build everything and find the shared colour range ---------------------------
panels = Array{Any}(undef, length(ROWS), length(SEPS))
vmin, vmax = Inf, -Inf
for (r, row) in enumerate(ROWS), (c, a) in enumerate(SEPS)
    panels[r, c] = make_pair(a; beta = row.beta, irradiate = row.irradiate, tepoch = TEPOCH)
    st1, st2, m1, m2, _ = panels[r, c]
    d1, d2 = RELATIVE ? (m1 .- TPOLE1, m2 .- TPOLE2) : (m1, m2)
    # VISIBLE tessels only. Ranging over the whole mesh includes the far hemisphere — whose
    # far-side tidal bulge is the most gravity-darkened point on the star and is never seen —
    # so the scale is set by data that never appears and the ends of the colormap go unused.
    v1, v2 = d1[st1.index_quads_visible], d2[st2.index_quads_visible]
    global vmin = min(vmin, minimum(v1), minimum(v2))
    global vmax = max(vmax, maximum(v1), maximum(v2))
    @printf("  row %d (%-22s) a = %5.2f mas : T = %6.0f .. %6.0f K\n",
            r, row.label, a, min(minimum(m1), minimum(m2)), max(maximum(m1), maximum(m2)))
end
if RELATIVE          # symmetric about zero, or "no change" is not the middle colour
    lim = max(abs(vmin), abs(vmax)); vmin, vmax = -lim, lim
end
# Colour limit from a high PERCENTILE, not the maximum: one extreme panel (the closest
# separation with gravity darkening on) otherwise sets the scale for all twelve and washes
# the rest out. Percentile keeps them comparable while leaving typical structure visible.
if RELATIVE
    # Symmetric about zero so that white always means "unperturbed" — an asymmetric range
    # would move the neutral colour off 0 and the map would stop being readable as
    # cooler/hotter. CLIP < 1 trades a little saturation at the extremes for more contrast
    # everywhere else; CLIP = 1 (default) uses the true visible extremes.
    allmag = Float64[]
    for r in 1:length(ROWS), c in 1:length(SEPS)
        st1, st2, m1, m2, _ = panels[r, c]
        append!(allmag, abs.((m1 .- TPOLE1)[st1.index_quads_visible]))
        append!(allmag, abs.((m2 .- TPOLE2)[st2.index_quads_visible]))
    end
    clip = parse(Float64, get(ENV, "CLIP", "1.0"))
    lim  = clip >= 1 ? maximum(allmag) : quantile(allmag, clip)
    @printf("visible ΔT: %.0f .. %.0f K  ->  symmetric limit %.0f K (gamma = %.2f)\n",
            vmin, vmax, lim, GAMMA)
    global TLIM = lim                       # kept in K, for the colour bar ticks
    vmin, vmax = -stretch(lim), stretch(lim)
end
@printf("shared colour scale: %.0f .. %.0f K%s\n", vmin, vmax,
        RELATIVE ? "  (ΔT from tpole, diverging, clipped at the 99.5th pct)" : "")

# FRAME. `plot2d_binary` centres its axes on star 1, so at quadrature — where the two are
# maximally separated — the frame has to reach a full separation on BOTH sides and half of
# every panel is empty. Recentring on the pair MIDPOINT after the call halves the required
# half-width, i.e. doubles the stars on screen, which is the whole point of viewing at
# quadrature in the first place.
#
# Half-width is shared across panels (so the shrinking orbit stays legible) but each panel is
# centred on its own midpoint.
half = 0.0
for r in 1:length(ROWS), c in 1:length(SEPS)
    star1, star2, _, _, bp = panels[r, c]
    ow, on = orbit_to_rotir_offset(bp, TEPOCH)
    r1 = maximum(hypot.(star1.proj_west, star1.proj_north))
    r2 = maximum(hypot.(star2.proj_west, star2.proj_north))
    global half = max(half, hypot(ow, on)/2 + max(r1, r2))
end
global AXMAX = 1.10 * half
@printf("frame: half-width = %.2f mas, centred per panel on the pair midpoint\n", AXMAX)

# --- pass 2: draw ------------------------------------------------------------------------
ROTIR.set_oiplot_defaults()
fig, axs = pyplot.subplots(length(ROWS), length(SEPS),
                           figsize = (3.6 * length(SEPS), 3.6 * length(ROWS)))
for (r, row) in enumerate(ROWS), (c, a) in enumerate(SEPS)
    star1, star2, m1, m2, bp = panels[r, c]
    d1, d2 = RELATIVE ? (m1 .- TPOLE1, m2 .- TPOLE2) : (m1, m2)
    if RELATIVE; d1 = stretch.(d1); d2 = stretch.(d2); end
    ax = axs[r - 1][c - 1]          # 0-based: PythonCall hands back a numpy array of axes
    plot2d_binary(d1, d2, star1, star2, bp, TEPOCH;
                  ax = ax, colorbar_on = false, vmin = vmin, vmax = vmax,
                  colormap = CMAP, axis_max = AXMAX,
                  compass = (r == length(ROWS) && c == length(SEPS)),
                  limb = true, plotmesh = false)
    # annotate what the physics actually did in this panel
    k1 = (m1 .- TPOLE1)[star1.index_quads_visible]
    k2 = (m2 .- TPOLE2)[star2.index_quads_visible]
    spread = max(maximum(k1) - minimum(k1), maximum(k2) - minimum(k2))
    ax.text(0.03, 0.03, @sprintf("ΔT = %.0f K", spread), transform = ax.transAxes,
            fontsize = 9, color = "0.25", ha = "left", va = "bottom")
    # recentre on the pair midpoint: star 1 sits at the origin and star 2 at (-ow, on)
    ow, on = orbit_to_rotir_offset(bp, TEPOCH)
    cx, cy = -ow/2, on/2
    ax.set_xlim(cx + AXMAX, cx - AXMAX)      # x inverted: East to the left
    ax.set_ylim(cy - AXMAX, cy + AXMAX)
    r == 1 && ax.set_title(@sprintf("a = %.2f mas   (r/a = %.2f)", a, RPOLE1/a), fontsize = 11)
    c == 1 && ax.set_ylabel(row.label, fontsize = 10)
    ax.set_xlabel(""); c == 1 || ax.set_ylabel("")
    ax.set_xticks([]); ax.set_yticks([])
end

# one colour bar for the whole figure — the panels are only comparable because they share it
cmap = pyimport("matplotlib").colormaps[CMAP]
norm = pyimport("matplotlib").colors.Normalize(vmin = vmin, vmax = vmax)
sm   = pyimport("matplotlib").cm.ScalarMappable(norm = norm, cmap = cmap)
cb   = fig.colorbar(sm, ax = axs, orientation = "vertical", fraction = 0.02, pad = 0.02)
if RELATIVE
    # The colour axis is in stretched units; label it with the kelvin values they mean.
    kt = [k for k in (-3000, -1000, -300, -100, 0, 100, 300, 1000, 3000) if abs(k) <= TLIM]
    cb.set_ticks(pylist(stretch.(Float64.(kt))))
    cb.set_ticklabels(pylist(string.(kt)))
end
cb.set_label(RELATIVE ? "ΔT from unperturbed pole (K)" : "Temperature (K)")
# The unperturbed pole temperatures belong in the title: every panel is plotted RELATIVE to
# them, so ΔT has no meaning without knowing what it is a departure from.
fig.suptitle("Proximity effects: Roche distortion, gravity darkening, mutual irradiation\n" *
             "unperturbed pole temperatures $(round(Int, TPOLE1)) K (primary) / " *
             "$(round(Int, TPOLE2)) K (secondary),  mass ratio q = $QBIN",
             fontsize = 12)
fig.savefig(OUT, dpi = parse(Int, get(ENV, "DPI", "150")), bbox_inches = "tight",
            facecolor = "white")
pyplot.close("all")
println("wrote $OUT")
