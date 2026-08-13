# betlyr_params.jl — β Lyrae (β Lyr A, Sheliak) system definition.
#
# `include("betlyr_params.jl")` from the fitting scripts. Builds parameter NamedTuples and
# prints an audit; reads no data and draws nothing.
#
# ---------------------------------------------------------------------------
# SOURCES (both PDFs are in this folder)
# ---------------------------------------------------------------------------
#   [Z08] Zhao et al. 2008, ApJ 684, L95 — "First resolved images of the eclipsing and
#         interacting binary β Lyrae". CHARA/MIRC H band, 6 epochs 2006-2007. Fits an
#         astrometric orbit with two uniform ellipses.
#   [M18] Mourard, Brož, Nemravová et al. 2018, A&A 618, A112 — "Physical properties of
#         β Lyrae A and its opaque accretion disk". SHELLSPEC radiative transfer, 27 289
#         observations 1970-2017. (The PDF metadata says 616; the running head and CDS say
#         618. 618 is correct.) NB [M18] deliberately EXCLUDED the 2006/2007 MIRC data
#         (their Table 1 note a) because it predates MIRC's photometric channels — so the
#         `older_data/` nights here were never used by them.
#   [Ak07] Ak et al. 2007, the quadratic ephemeris, quoted as Eq. (1) of [M18].
#
# ---------------------------------------------------------------------------
# WHY β LYR IS NOT SPICA
# ---------------------------------------------------------------------------
#   1. CIRCULAR ORBIT. e = 0 exactly ([M18] Table 7; [Z08] §4, both citing Harmanec &
#      Scholz 1993). ω is then undefined — [M18] pins it at 90° purely by convention. With
#      e = 0 the residual degeneracy is (Ω, T0) → (Ω+180°, T0 + P/2), so Ω is restricted to
#      [0,180) and T0 sampled over a full period.
#
#   2. THE PERIOD IS CHANGING, MEASURABLY. Mass transfer at ~2×10⁻⁵ M☉/yr lengthens the
#      orbit. Over the 6.7 yr / 189 orbits these data span, the quadratic ephemeris
#      accumulates ≈0.14 orbits ≈ 50° of phase. dP is a first-class fitted parameter here.
#
#   3. THE GAINER IS HIDDEN IN AN OPAQUE DISC. The donor is a B6-8 II star exactly filling
#      its Roche lobe; the gainer (B0.5 IV-V, R = 6 R☉) is buried in a disc of outer radius
#      30 R☉ — five times the star — which carries most of the H-band flux. Two uniform
#      discs, adequate for Spica, cannot describe this.
#
# ---------------------------------------------------------------------------
# EPHEMERIS: USE A MODERN EPOCH, NOT THE 1881 ONE
# ---------------------------------------------------------------------------
# [Ak07]'s reference epoch is HJD 2408247.968 (1881). Propagating it to our 2013 data means
# E ≈ 189 600 cycles, where the quadratic term reaches 139 189 days — and the published
# uncertainty on c alone, ±0.00369e-6 d/E², maps to ±133 days ≈ ±10.3 ORBITS. The phase is
# completely unconstrained that way.
#
# So the reference epoch is [Z08]'s T0 = JD 2454283.0430 (2007 Jul 1), which sits inside
# the data span. It is exactly E = 3561 on the [Ak07] ephemeris (verified: 2408247.968 +
# 12.913779·3561 + 3.87265e-6·3561² = 2454283.0430), and the instantaneous period there is
# P₀ + 2cE = 12.94136 d, matching [Z08]'s adopted 12.9414 d. A local re-expansion, not a
# different orbit.
#
# ROTIR's dP convention is the same as [Ak07]'s: `compute_eccentric_anomaly` uses
# t_n = T0 + P·n + ½(Ṗ·P)n², so ½ṖP = c ⟹ Ṗ = 2c/P₀ = 5.9977e-7 d/d — reproducing [M18]
# Table 7's Ṗ to five digits. (≈18.93 s/yr, the "≈19 s/yr" of the literature.)
#
# ⚠ Two traps in the sources, verified numerically, recorded so they are not propagated:
#   * [M18] Table 7's T_min = 2408254.4248895 is offset from their own Eq. (1) epoch by
#     EXACTLY P₀/2 (6.4568895 d). Implementing Table 7 literally puts you half a period out.
#   * [Z08] prints "a sin i = 57.87 ± 0.62 M☉" and "0.46 mas (29.4 M☉)" — both should be
#     R☉ (confirmed by unit arithmetic and by [M18] §3.1.2).

using ROTIR, Printf

# --- unit conversion (same helpers as demos/spica_params.jl) --------------------------
const R_SUN_M = 6.957e8
const PC_M    = 3.0856775814e16
const RAD_MAS = 206264806.247
ang_radius_mas(R_Rsun, d_pc) = R_Rsun * R_SUN_M / (d_pc * PC_M) * RAD_MAS
phys_radius_rsun(θ_mas, d_pc) = θ_mas / RAD_MAS * (d_pc * PC_M) / R_SUN_M

# --- ephemeris -----------------------------------------------------------------------
const P_AK   = 12.913779          # d, [Ak07] via [M18] Eq. (1)  (±0.000016)
const C_AK   = 3.87265e-6         # d/E², quadratic term         (±0.00369e-6)
const T0_AK  = 2408247.968        # HJD, epoch of primary minimum (±0.015)
const P_ORB  = parse(Float64, get(ENV, "PORB", "12.9414"))       # [Z08] local period at E=3561
const T0_ORB = parse(Float64, get(ENV, "T0",  "2454283.0430"))   # [Z08] 2007 Jul 1, E=3561
const DP_DD  = parse(Float64, get(ENV, "DP",  string(2C_AK / P_AK)))  # Ṗ, d/d
const DPDT_SY = DP_DD * 365.25 * 86400                            # s/yr, for reporting

# Phase 0 = superior conjunction of the mass DONOR, i.e. primary minimum ([M18] §2).
const E_ORB      = 0.0            # circularised, exactly ([M18] Table 7; [Z08] §4)
const OMEGA_PERI = 90.0           # [M18]'s convention; meaningless at e = 0, pinned

# --- orbit ---------------------------------------------------------------------------
const I_ORB      = parse(Float64, get(ENV, "INCL",   "93.5"))   # deg, [M18] §6 (±1.0)
const OMEGA_NODE = parse(Float64, get(ENV, "OMNODE", "253.7"))  # deg, [M18] §6 (±1.0)
const A_ORB      = parse(Float64, get(ENV, "AMAS",   "0.865"))  # mas, [Z08] Table 3 (±0.048)
const D_PC       = parse(Float64, get(ENV, "DPC",    "319.7"))  # pc, [M18] §6 (±2.7)
# ROTIR's q is M2/M1 with body 1 at the origin. Component 1 here is the DONOR, so
# q = m_gainer/m_donor = 4.50 ([M18] §3.2; = 1/0.222 of [Z08] fn. 8, consistent).
# NOTE it does not enter this fit at all: `orbit_to_rotir_offset` returns the RELATIVE
# separation, which is independent of the mass split. It matters only for Roche geometry.
const Q_BIN      = parse(Float64, get(ENV, "QBIN", "4.50"))

# The orbit is RETROGRADE — i > 90° and position angle decreases with time ([Z08] §4).

# --- components ----------------------------------------------------------------------
# Donor: modelled as a circular uniform disc. [Z08] §3 measured a uniform ELLIPSE of
# 0.62 ± 0.16 × 0.52 ± 0.14 mas (full axes, averaged over epochs) — the elongation is their
# direct detection of Roche-lobe filling. The circular mean is used as a starting value;
# M18's Roche geometry implies ≈0.44 mas, so this is poorly resolved and the fit should be
# allowed to move it.
const UD_DONOR = parse(Float64, get(ENV, "UDDONOR", "0.57"))    # mas

# Gainer + disc: [Z08] §3 measured 1.04 ± 0.11 × 0.63 ± 0.07 mas (full axes). [M18] finds
# R_out = 30.0 ± 1.0 R☉ and semi-thickness H = 6.5 ± 1.0 R☉ — a thick disc that cannot be
# in vertical hydrostatic equilibrium. Modelled here as an inclined elliptical Gaussian.
# Its major axis lies along the line of nodes, so PA ≈ Ω (mod 180).
const DISC_FWHM  = parse(Float64, get(ENV, "DISCFWHM",  "0.90"))            # mas
const DISC_RATIO = parse(Float64, get(ENV, "DISCRATIO", string(0.63/1.04))) # minor/major
const DISC_PA    = parse(Float64, get(ENV, "DISCPA", string(mod(253.7, 180.0))))

# H-band flux. [Z08] §3 adopt donor/disc = 1.24, i.e. donor ≈55%, disc ≈45% out of eclipse.
# Our model parameter is the ratio disc/donor = 1/1.24.
const F_DISC = parse(Float64, get(ENV, "FDISC", string(1/1.24)))

"""
    betlyr_audit()

Print the system in several unit systems at once, and self-check the ephemeris against
[Ak07]. As with Spica, a radius/diameter slip is invisible in mas but obvious in R☉, so
both are always shown.
"""
function betlyr_audit()
    mas_per_rsun = ang_radius_mas(1.0, D_PC)
    println("─"^78)
    println("β Lyrae — parameter audit   (d = $D_PC pc, 1 R☉ = $(round(mas_per_rsun, sigdigits=4)) mas)")
    println("─"^78)
    @printf("  P = %.5f d at T0 = %.4f (JD)   [Z08], E = 3561 on the Ak07 ephemeris\n",
            P_ORB, T0_ORB)
    @printf("  Ṗ = %+.4e d/d = %+.2f s/yr    [Ak07] c = %.5e d/E² ⟹ Ṗ = 2c/P₀\n",
            DP_DD, DPDT_SY, C_AK)
    # Self-check: does the modern epoch sit on the 1881 quadratic ephemeris? The cycle
    # count must come from SOLVING P·E + c·E² = ΔT, not from the linear ΔT/P — over 3561
    # cycles the quadratic term is already 49 d, so the linear estimate is off by 4 cycles.
    E = round((-P_AK + sqrt(P_AK^2 + 4C_AK*(T0_ORB - T0_AK))) / (2C_AK))
    tpred = T0_AK + P_AK*E + C_AK*E^2
    @printf("  ephemeris check: E = %.0f ⟹ T = %.4f, adopted %.4f, Δ = %+.4f d %s\n",
            E, tpred, T0_ORB, tpred - T0_ORB, abs(tpred - T0_ORB) < 0.01 ? "✓" : "✗")
    @printf("  e = %.1f exactly (circularised) ⇒ ω undefined (pinned at %.0f°);\n", E_ORB, OMEGA_PERI)
    println("     residual degeneracy is (Ω, T0) → (Ω+180°, T0 + P/2), handled by Ω ∈ [0,180)")
    @printf("  i = %.2f° (retrograde: i > 90°, PA decreasing),  Ω = %.2f°\n", I_ORB, OMEGA_NODE)
    @printf("  a = %.4f mas = %.2f R☉   [Z08]; [M18] a sin i = 58.19 R☉ ⟹ %.4f mas\n",
            A_ORB, phys_radius_rsun(A_ORB, D_PC), 58.19/sin(deg2rad(I_ORB)) * mas_per_rsun)
    @printf("  q = m_g/m_d = %.2f  (unused: the relative orbit is independent of it)\n", Q_BIN)
    println()
    @printf("  donor  UD  %.4f mas = %.2f R☉ diam   [Z08] ellipse 0.62 × 0.52 mas\n",
            UD_DONOR, phys_radius_rsun(UD_DONOR, D_PC))
    # Sanity: the donor exactly fills its Roche lobe, so its diameter cannot exceed ~2R_L.
    RLd = 2 * 15.1                                   # R☉, Eggleton lobe for q=4.5, a=58.3 R☉
    if phys_radius_rsun(UD_DONOR, D_PC) > 1.15RLd
        @printf("         ⚠ that is %.0f%% larger than the Roche lobe diameter (%.1f R☉);\n",
                100*(phys_radius_rsun(UD_DONOR, D_PC)/RLd - 1), RLd)
        println("           [Z08] note the donor is only partially resolved, so let the fit shrink it")
    end
    @printf("  disc  FWHM %.4f mas = %.2f R☉ , ratio %.2f, PA %.1f°, flux/donor %.2f\n",
            DISC_FWHM, phys_radius_rsun(DISC_FWHM, D_PC), DISC_RATIO, DISC_PA, F_DISC)
    @printf("         [M18] R_out = 30.0 ± 1.0 R☉ (diam %.2f mas), H = 6.5 ± 1.0 R☉\n",
            2*30.0*mas_per_rsun)
    println()
    T = 6.71 * 365.25
    @printf("  over the %.1f yr baseline the quadratic term moves phase by %.3f orbits (%.0f°)\n",
            T/365.25, 0.5*DP_DD*T^2/P_ORB, 360*0.5*DP_DD*T^2/P_ORB)
    println("─"^78)
end
