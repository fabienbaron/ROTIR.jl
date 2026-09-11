#!/usr/bin/env julia
# =======================================================================================
# Spica through one orbit, with the proximity effects left in
# =======================================================================================
# The companion to `proximity_effects.jl`. That figure freezes the system at one phase and
# separates the three close-binary effects by turning them off one at a time. This one does
# the opposite: it holds all three ON and sweeps the epoch over a COMPLETE PERIOD of a real
# system, so what you watch is the thing a static panel cannot show — the effects CHANGING.
#
# Three of them, and they move for different reasons:
#
#   ROCHE DISTORTION  follows the instantaneous separation D(t). Spica's orbit is eccentric
#                     (e = 0.123), so D runs from 1.35 to 1.73 mas and each star is squeezed
#                     and released once per orbit. A circular orbit would show none of this.
#   GRAVITY DARKENING rides on that shape: T = tpole (g/g_pole)^beta is computed on the
#                     distorted surface, so the tidal bulge cools as it grows.
#   IRRADIATION       depends on which face is presented to the companion and how far away
#                     it is, so it waxes and wanes with both the separation and the phase.
#
# And the two geometric things that make it worth animating at all: the components rotate
# (tidally locked, so one face always leads), and the pair swings across the sky, eclipsing
# near conjunction at Spica's inclination of 116 degrees.
#
# VOLUME IS THE INVARIANT, not the polar radius. As D changes the Roche potential changes,
# and holding `rpole` fixed would make the star physically breathe — gaining and losing
# matter once per orbit. `volume_conserving = true` instead solves for the equipotential
# that preserves each component's volume at every separation, which is what a star actually
# does on an orbital timescale. It is a root solve over a quadrature, so `binary_movie`
# tabulates Omega(D) once and interpolates rather than paying it per frame. Set
# VOLUME=0 to see the (wrong, but cheaper) fixed-rpole version.
#
#   julia --project=demos demos/spica_proximity_movie.jl
#   NSIDE=3 NFRAMES=60 julia --project=demos demos/spica_proximity_movie.jl     # a quick look
#   PANELS=compare julia --project=demos demos/spica_proximity_movie.jl          # 3-up diagnostic
#
# PANELS=compare gives the three-panel version — irradiated, intrinsic, and the difference —
# which is the one that isolates what irradiation alone contributes at each phase. The
# default single panel is the system as it would be observed.
#
# Output is an mp4 (`ffmpeg` on PATH; the PNG frames are kept either way, and the encode
# command is printed if ffmpeg is missing).

ENV["MPLBACKEND"] = get(ENV, "MPLBACKEND", "Agg")
using ROTIR, PythonPlot, Printf

const NSIDE   = parse(Int, get(ENV, "NSIDE", "4"))
const NFRAMES = parse(Int, get(ENV, "NFRAMES", "120"))
const FPS     = parse(Int, get(ENV, "FPS", "24"))
const PANELS  = Symbol(get(ENV, "PANELS", "single"))
const ALBEDO  = parse(Float64, get(ENV, "ALBEDO", "0.6"))
const VOLUME  = get(ENV, "VOLUME", "1") == "1"
const OUTDIR  = get(ENV, "OUTDIR", joinpath(@__DIR__, "results", "spica_movie_frames"))
const OUT     = get(ENV, "OUT", joinpath(@__DIR__, "spica_proximity_orbit.mp4"))

# --- Spica ------------------------------------------------------------------------------
# alpha Vir: a close, eccentric, non-synchronous B-type pair, and the best-resolved one
# there is. Elements and component parameters as in `spica_binary_roche.jl`
# (Aufdenberg+2015 angular diameters, Tkachenko+2016 temperatures).
const P_ORB = 4.0145        # d
const A_ORB = 1.54          # mas, semi-major axis of the RELATIVE orbit
const E_ORB = 0.123
const Q_BIN = 0.6188        # M2/M1
const I_ORB = 116.0         # deg; > 90 is retrograde
const OMEGA = 309.938       # deg, ascending node
const OMEGA_P = 255.0       # deg, argument of periapsis of the relative orbit
const DIST  = 77.0          # pc

const RPOLE1, RPOLE2 = 0.93 / 2, 0.57 / 2      # mas, polar radii
const TPOLE1, TPOLE2 = 25300.0, 20585.0        # K
const BETA  = 0.25          # von Zeipel, radiative envelope
const LD1   = 0.15          # Hestroffer power law

# THE TIME ORIGIN IS PERIASTRON, not the published JD. Spica's T0 is 2454189.40, and an
# absolute Julian date is a poor argument to carry through a Float32 mesh: eps(Float32(2.45e6))
# is a QUARTER OF A DAY, which on a 4-day orbit is a sixteenth of a phase. Time DIFFERENCES
# are fine, so the elements are shifted to T0 = 0 and the sweep runs 0 -> P. The system is
# unchanged; only the zero of the clock moves. `spica_binary_roche.jl` sidesteps the same
# trap by passing t = 0 for the geometry and computing the offsets separately.
const T0 = 0.0

# A Roche component's parameters must ALSO carry the orbital elements: `create_binary_star`
# reads `p.a` to express the instantaneous separation in units of the semi-major axis, and
# `compute_separation` reads `P`, `e` and `T0`. A NamedTuple without them fails with a bare
# `getproperty` error several calls down.
#
# `inclination` and `position_angle` follow the ORBIT — a tidally locked component's spin axis
# is the orbital one — which is the `180 - i`, `Omega - 180` convention the rest of ROTIR uses
# for binaries.
roche_params(rpole, tpole, q) = (
    surface_type              = 3,
    rpole                     = rpole,
    tpole                     = tpole,
    ldtype                    = 3,               # Hestroffer power law
    ld1                       = LD1,
    ld2                       = 0.0,
    inclination               = 180.0 - I_ORB,
    position_angle            = OMEGA - 180.0,
    rotation_period           = P_ORB,           # tidally locked
    beta                      = BETA,
    d                         = DIST,
    q                         = q,               # M_companion / M_self
    fillout_factor_primary    = -1,              # rpole defines the equipotential
    fillout_factor_secondary  = -1,
    i = I_ORB, Ω = OMEGA, ω = OMEGA_P, P = P_ORB, a = A_ORB, e = E_ORB, T0 = T0,
    dP = 0.0, dω = 0.0,
)

params1 = roche_params(RPOLE1, TPOLE1, Q_BIN)
params2 = roche_params(RPOLE2, TPOLE2, 1 / Q_BIN)    # inverted for the secondary potential

star1p = starparameters(RPOLE1, TPOLE1, 0.0, 3, LD1, 0.0, BETA, 0.0,
                        180.0 - I_ORB, OMEGA - 180.0, 0.0, P_ORB)
star2p = starparameters(RPOLE2, TPOLE2, 0.0, 3, LD1, 0.0, BETA, 0.0,
                        180.0 - I_ORB, OMEGA - 180.0, 0.0, P_ORB)
bparams = binaryparameters(star1p, star2p, DIST, I_ORB, OMEGA, OMEGA_P,
                           P_ORB, A_ORB, E_ORB, T0, Q_BIN, [-1.0, -1.0], 0.0, 0.0)

# --- how deep in its lobe does each star sit, and where does it get deepest? --------------
# PERIASTRON is the test, not the mean separation: `fillout_factor = -1` derives the potential
# from `rpole`, and nothing stops rpole exceeding the lobe — the radii then diverge toward L1
# and you get a garbage shape rather than an error. On an eccentric orbit the closest approach
# is where that happens first.
const RLOBE = radius_eggleton(Q_BIN)
const DPERI = A_ORB * (1 - E_ORB)
const DAPO  = A_ORB * (1 + E_ORB)
@printf("Spica: P = %.4f d, a = %.2f mas, e = %.3f, q = %.4f, i = %.1f deg\n",
        P_ORB, A_ORB, E_ORB, Q_BIN, I_ORB)
@printf("separation: %.3f mas at periastron .. %.3f mas at apastron  (%.0f%% swing)\n",
        DPERI, DAPO, 100 * (DAPO - DPERI) / A_ORB)
@printf("Roche lobe r_L/a = %.3f (Eggleton, q = %.4f)\n", RLOBE, Q_BIN)
for (nm, rp) in (("primary", RPOLE1), ("secondary", RPOLE2))
    fill_peri = rp / (DPERI * RLOBE)
    @printf("  %-9s rpole = %.3f mas : fills %3.0f%% of its lobe at periastron\n",
            nm, rp, 100 * fill_peri)
    fill_peri < 0.95 || error("""
        the $nm reaches $(round(100*fill_peri))% of its Roche lobe at periastron
        (D = $(round(DPERI, digits=3)) mas). The shape diverges toward L1 rather than
        failing, so this is refused. Reduce its polar radius or widen the orbit.""")
end

@printf("\n%d frames over one period at nside = %d, %s panel%s, volume-conserving = %s\n",
        NFRAMES, NSIDE, PANELS, PANELS === :single ? "" : "s", VOLUME)

tess1 = tessellation_healpix(NSIDE)
tess2 = tessellation_healpix(NSIDE)

# One period, periastron to periastron. `binary_movie` defaults to exactly this — `T0` to
# `T0 + P` — and it is spelled out here because it is the point of the demo.
dir, movie = binary_movie(bparams, tess1, params1, tess2, params2;
                          tstart = T0, tstop = T0 + P_ORB,
                          nframes = NFRAMES, fps = FPS, panels = PANELS,
                          outdir = OUTDIR, prefix = "spica",
                          reflection = true, albedo1 = ALBEDO, albedo2 = ALBEDO,
                          volume_conserving = VOLUME,
                          intensity = true, intensity_model = :planck, band = 1.65e-6,
                          colormap = "gist_heat",
                          # The phase is what a reader needs to locate a frame in the orbit;
                          # the separation says where in the ECCENTRIC orbit it is, which a
                          # phase alone does not on a non-circular one.
                          title_fn = (t, phase, sep) -> @sprintf(
                              "Spica  —  phase %.3f   ρ = %.3f mas   D = %.3f mas",
                              phase, sep, compute_separation(bparams, t) * A_ORB))

if movie !== nothing && movie != OUT
    mv(movie, OUT; force = true)
    println("wrote $OUT")
elseif movie === nothing
    println("frames are in $dir; encode them with the ffmpeg command above")
end
