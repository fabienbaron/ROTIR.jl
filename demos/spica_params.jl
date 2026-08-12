# spica_params.jl — one place where the Spica system is defined.
#
# `include("spica_params.jl")` from the irradiation demos. Nothing here reads data or
# plots; it only builds parameter NamedTuples and prints a units audit.
#
# ---------------------------------------------------------------------------
# UNITS — read this before changing a number
# ---------------------------------------------------------------------------
# ROTIR's `rpole` is an angular **RADIUS** in mas. Published interferometric values are
# very often angular **DIAMETERS**, and published model values are usually physical radii
# in R☉ — so this file takes the physical radii from the literature and converts, rather
# than hard-coding an angular number whose factor-of-two provenance is invisible.
#
# The pre-existing `demos/spica_binary_roche.jl` hard-codes `rpole = 0.57/2 = 0.285 mas`
# for the secondary, which is 4.72 R☉ at 77 pc — 26% larger than the 3.74–3.76 R☉ both
# Aufdenberg+2015 and Tkachenko+2016 report. (Its primary, 0.465 mas = 7.70 R☉, is only 4%
# off 7.40 R☉.) The values below supersede those.
#
# Sources:
#   Aufdenberg et al. 2015 (papers/spica_the_paper_2015.pdf)
#       "The primary polar radius is fixed at 7.40 R☉"; mean equatorial radius of the
#       primary 7.93 R☉ (periastron) → 8.01 R☉ (apastron); mean equatorial radius of the
#       secondary 3.76 R☉, varying by <1% over the orbit.
#       Parallax 12.44 ± 0.86 mas (Hipparcos) / 13.06 ± 0.70 mas (van Leeuwen 2007).
#       Notes a *grazing eclipse* in the MOST photometry (Desmet et al. 2009) — see the
#       overlap audit below, which reproduces it.
#   Tkachenko et al. 2016, MNRAS 458, 1964 (papers/stw255.pdf), Tables 3 & 6
#       R1 = 7.47 ± 0.54 R☉, R2 = 3.74 ± 0.53 R☉;
#       Teff = 25300 ± 500 K / 20585 ± 850 K.
#   Wages & Aufdenberg, "Constraining the Apsidal Constant of Spica by Modelling
#   Ellipsoidal Variations" (papers/Spica_poster8.pdf)
#       e = 0.065 ± 0.015, apsidal period U = 105 ± 2 yr.
#       Explicitly: their model "lacks the reflection effect, which could improve the
#       model fit at phase = 0.56", and one light-curve point WITH reflection "took 127
#       hours with 210 wavelength points" — i.e. the physics is wanted but the existing
#       implementation is too slow to use. That is the gap `src/reflection.jl` fills: the
#       radiosity solve here is a few dense gemv's and runs in milliseconds per epoch.
#
# ---------------------------------------------------------------------------
# ECCENTRICITY IS CONTESTED — and it sets the irradiation modulation
# ---------------------------------------------------------------------------
# Modern values for Spica's e still span nearly a factor of two, and the disagreement is
# acknowledged in the literature rather than resolved:
#     Riddle 2000                    0.067 ± 0.014  (primary alone)
#     Aufdenberg et al. 2015         0.118 – 0.119  (interferometry + RV + spectrophotometry)
#     Robinette & Aufdenberg 2015    0.116 ± 0.004  (MCMC on 1889–2000 RVs; tightest, 4%)
#                                    — explicitly "inconsistent with the value of Riddle (2000)"
#     Wages & Aufdenberg (poster)    0.065 ± 0.015  (ellipsoidal-variation fit)
# The peak-to-peak irradiation contrast over the orbit goes as ((1+e)/(1−e))², i.e.
# 1.30 at e = 0.065 and 1.61 at e = 0.118 — so the *modulation*, which is the cleanest
# signature of the reflection effect, is uncertain at the tens-of-percent level from e
# alone. Override with the ECC environment variable to test sensitivity.
#
# Robinette & Aufdenberg also give q = M2/M1 = 0.622 ± 0.005 and a·sin i = 23.14 ± 0.88 R☉
# (→ a = 25.7 R☉ = 1.555 mas at i = 116°, d = 77 pc), both consistent with the values used
# below, and an apsidal period U = 139 ± 6 yr (P/U = 7.91 ± 0.36 ×10⁻⁵).
# ---------------------------------------------------------------------------

using ROTIR, Printf

# --- unit conversion ---------------------------------------------------------------
const R_SUN_M = 6.957e8          # m
const PC_M    = 3.0856775814e16  # m
const RAD_MAS = 206264806.247    # mas per radian

"Angular RADIUS in mas of a star of physical radius `R` (R☉) at distance `d` (pc)."
ang_radius_mas(R_Rsun, d_pc) = R_Rsun * R_SUN_M / (d_pc * PC_M) * RAD_MAS

"Physical radius in R☉ corresponding to an angular RADIUS `θ` (mas) at distance `d` (pc)."
phys_radius_rsun(θ_mas, d_pc) = θ_mas / RAD_MAS * (d_pc * PC_M) / R_SUN_M

# --- system ---------------------------------------------------------------------------
const D_PC   = 77.0        # pc  (π = 12.99 mas, between the Hipparcos and orbital values)
const P_ORB  = 4.0145      # d
const A_ORB  = 1.54        # mas, semi-major axis of the RELATIVE orbit
const E_ORB  = parse(Float64, get(ENV, "ECC", "0.123"))   # Aufdenberg+2015; see the note above
const T0_ORB = 2454189.40  # JD, periastron
const Q_BIN  = 0.6188      # M2/M1
const I_ORB  = 116.0       # deg (>90 ⇒ retrograde on the sky)
const OMEGA_NODE = 309.938 # deg, longitude of ascending node Ω
const OMEGA_PERI = 255.0   # deg, argument of periapsis ω (relative orbit)

# Physical radii (R☉) → angular polar radii (mas).
# The secondary is spherical to <1%, so its polar radius is taken equal to its mean
# equatorial radius; the primary's polar radius is quoted directly by Aufdenberg+2015.
const R1_RSUN = 7.40
const R2_RSUN = 3.74
const RPOLE1 = ang_radius_mas(R1_RSUN, D_PC)
const RPOLE2 = ang_radius_mas(R2_RSUN, D_PC)

const TPOLE1 = 25300.0     # K
const TPOLE2 = 20585.0     # K

# --- asynchronous rotation --------------------------------------------------------
# Spica's components are NOT tidally locked. From v sin i = 161 ± 2 / 70 ± 5 km/s
# (Smith 1985; Tkachenko+2016) with R = 7.40 / 3.74 R☉ and i = 116°:
#     v_eq(synchronous) = 2πR/P = 93.3 / 47.1 km/s
#     v_eq(observed)    = v sin i / |sin i| = 179.1 / 77.9 km/s
#     F = ω_rot/ω_orb   = 1.92 / 1.65
# The centrifugal term goes as F², so this is a 3.7× / 2.7× change: the primary's
# equatorial bulge is 9.4%, not the 2.1% a synchronous model gives. Aufdenberg+2015 note
# the same ("asynchronous rotation factors").
#
# NB the Roche shape is still taken to be static in the frame where the tidal bulge points
# at the companion — the standard approximation (PHOEBE does the same). That is exact here
# because our surface maps (gravity darkening + irradiation) are symmetric about that
# frame; it would need revisiting for spots, which would rotate at ω_rot.
const F1_SYNC = parse(Float64, get(ENV, "F1", "1.92"))
const F2_SYNC = parse(Float64, get(ENV, "F2", "1.65"))
const USE_ASYNC = get(ENV, "ASYNC", "1") == "1"
_frot(F) = USE_ASYNC ? P_ORB / F : P_ORB          # rotation period in days

# --- apsidal motion ---------------------------------------------------------------
# U = 139 ± 6 yr (Robinette & Aufdenberg 2015, MCMC on 1889-2000 RVs; P/U = 7.91e-5).
# Also 105 ± 2 yr (Wages & Aufdenberg poster) and ~110 yr (Aufdenberg+2015). Over the
# 2007-2015 CHARA span ω advances ~21°, moving the predicted secondary position by up to
# 0.44 mas — ten to forty times the MIRC astrometric precision.
const U_APSIDAL  = parse(Float64, get(ENV, "UAPS", "139.0"))     # yr
const USE_APSIDAL = get(ENV, "APSIDAL", "1") == "1"
const DOMEGA = USE_APSIDAL ? 360.0 / (U_APSIDAL * 365.25) : 0.0   # deg/day
const BETA   = 0.25        # von Zeipel exponent for a radiative envelope.
                           # NB PHOEBE's gravb_bol is 4× this (T = Tpole·(g/gpole)^β here,
                           # T = Tpole·((g/gpole)^gravb_bol)^0.25 there).
const LDTYPE = 3           # Hestroffer power law, μ^ld1
const LD1_H  = 0.15        # H-band; NOTE limb darkening is band-dependent — change this
const LD2_H  = 0.0         # together with `band` if you move to R or V.

# Bolometric albedo for the reflection effect (PHOEBE's irrad_frac_refl_bol default is 0.6;
# 1.0 is the radiative-envelope expectation, so B stars are often modelled with 1.0).
const ALBEDO_DEFAULT = 0.6

# Common Roche/orbit block. Both components are assumed tidally locked (rotation_period =
# P), so async_ratio = 1 and the tidal bulge sits on the line of centres.
const SPICA_BASE = (surface_type = 3,                # Roche
                    ldtype = LDTYPE, ld1 = LD1_H, ld2 = LD2_H,
                    inclination = 180.0 - I_ORB,     # equivalent prograde viewing angle
                    position_angle = OMEGA_NODE - 180.0,
                    rotation_period = P_ORB,         # overridden per component below
                    beta = BETA, d = D_PC,
                    fillout_factor_primary = -1,     # -1 ⇒ rpole defines the potential
                    fillout_factor_secondary = -1,
                    i = I_ORB, Ω = OMEGA_NODE, ω = OMEGA_PERI,
                    P = P_ORB, a = A_ORB, e = E_ORB, T0 = T0_ORB,
                    dP = 0.0, dω = DOMEGA)

# q convention: the primary's potential uses q = M2/M1; the secondary's potential is
# centred on the secondary and needs the INVERTED ratio M1/M2.
const SPICA_P1 = merge(SPICA_BASE, (rpole = RPOLE1, tpole = TPOLE1, q = Q_BIN,
                                    rotation_period = _frot(F1_SYNC)))
const SPICA_P2 = merge(SPICA_BASE, (rpole = RPOLE2, tpole = TPOLE2, q = 1.0 / Q_BIN,
                                    rotation_period = _frot(F2_SYNC)))

const SPICA_S1 = starparameters(RPOLE1, TPOLE1, 0.0, LDTYPE, LD1_H, LD2_H, BETA, 0.0,
                                180.0 - I_ORB, OMEGA_NODE - 180.0, 0.0, _frot(F1_SYNC))
const SPICA_S2 = starparameters(RPOLE2, TPOLE2, 0.0, LDTYPE, LD1_H, LD2_H, BETA, 0.0,
                                180.0 - I_ORB, OMEGA_NODE - 180.0, 0.0, _frot(F2_SYNC))
const SPICA_BP = binaryparameters(SPICA_S1, SPICA_S2, D_PC, I_ORB, OMEGA_NODE, OMEGA_PERI,
                                  P_ORB, A_ORB, E_ORB, T0_ORB, Q_BIN, [1.0, 1.0], 0.0, DOMEGA)

"""
    spica_audit()

Print the system in several unit systems at once. The point is that a radius/diameter slip
is invisible in mas but obvious in R☉ and in R/a, so all three are shown side by side —
this is the check that catches "the stars look way too large".
"""
function spica_audit()
    a_rsun = phys_radius_rsun(A_ORB, D_PC)
    println("─"^78)
    println("Spica — parameter audit  (d = $D_PC pc,  1 R☉ = $(round(ang_radius_mas(1.0, D_PC), sigdigits=5)) mas)")
    println("─"^78)
    @printf("  semi-major axis a      %8.4f mas   %7.3f R☉\n", A_ORB, a_rsun)
    @printf("  primary polar radius   %8.4f mas   %7.3f R☉   R/a = %.4f   (diam %.4f mas)\n",
            RPOLE1, R1_RSUN, RPOLE1 / A_ORB, 2RPOLE1)
    @printf("  secondary polar radius %8.4f mas   %7.3f R☉   R/a = %.4f   (diam %.4f mas)\n",
            RPOLE2, R2_RSUN, RPOLE2 / A_ORB, 2RPOLE2)
    @printf("  (R1+R2)/a              %8.4f\n", (RPOLE1 + RPOLE2) / A_ORB)
    Dper = A_ORB * (1 - E_ORB); Dapo = A_ORB * (1 + E_ORB)
    @printf("  separation D           %8.4f mas (periastron) … %.4f mas (apastron)\n", Dper, Dapo)
    # Minimum sky-projected separation over one orbit.
    ts = T0_ORB .+ range(0, P_ORB, length=2001)
    ρmin = minimum(projected_separation.(Ref(SPICA_BP), collect(ts)))
    @printf("  min projected sep ρ    %8.4f mas  vs  R1+R2 = %.4f mas  →  %s\n",
            ρmin, RPOLE1 + RPOLE2,
            ρmin < RPOLE1 + RPOLE2 ? "GRAZING ECLIPSE (matches Desmet+2009 / Aufdenberg+2015)" :
                                     "no eclipse")
    # Substellar irradiation bump, single-scattering estimate, surface-to-centre distance.
    for (lab, D) in (("periastron", Dper), ("apastron", Dapo))
        d2 = D - RPOLE2; d1 = D - RPOLE1
        f2 = 1 + ALBEDO_DEFAULT * (RPOLE1 / d2)^2 * (TPOLE1 / TPOLE2)^4
        f1 = 1 + ALBEDO_DEFAULT * (RPOLE2 / d1)^2 * (TPOLE2 / TPOLE1)^4
        @printf("  ΔT_substellar (A=%.1f, %-10s)  secondary %+6.0f K   primary %+5.0f K\n",
                ALBEDO_DEFAULT, lab, TPOLE2 * (f2^0.25 - 1), TPOLE1 * (f1^0.25 - 1))
    end
    @printf("  synchronicity F = ω_rot/ω_orb : primary %.2f, secondary %.2f  %s\n",
            synchronicity(SPICA_P1), synchronicity(SPICA_P2),
            USE_ASYNC ? "(from v sin i = 161/70 km/s)" : "(FORCED SYNCHRONOUS, ASYNC=0)")
    @printf("  apsidal motion dω = %.5f deg/day (U = %.0f yr)  %s\n",
            DOMEGA, U_APSIDAL, USE_APSIDAL ? "→ Δω = $(round(DOMEGA*2910, digits=1))° over the data span" : "(DISABLED, APSIDAL=0)")
    @printf("  eccentricity e = %.3f  →  peri/apo irradiation contrast %.2f\n",
            E_ORB, ((1 + E_ORB) / (1 - E_ORB))^2)
    println("     (modern e determinations span 0.065–0.119, a disagreement the")
    println("      literature acknowledges but has not resolved; that alone moves the")
    println("      contrast 1.30–1.61. Set ECC=... to test the sensitivity.)")
    println("  NB the estimates above assume UNIFORM Tpole SPHERES. The full Roche model")
    println("     gives ~900 K on the secondary at periastron, not ~1030 K: gravity")
    println("     darkening cools the primary's tidal bulge — the very face doing the")
    println("     irradiating — which more than offsets the bulge being closer. The")
    println("     primary gets ~160 K instead of ~103 K, since the secondary's own bulge")
    println("     leans toward it. Converged to <1% by nside 4/3.")
    println("─"^78)
end
