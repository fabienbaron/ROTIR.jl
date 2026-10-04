#!/usr/bin/env julia
# HOW ROETTENBACHER ET AL. DID IT — their conventions, applied to these data, so that any
# comparison with their published numbers is like for like.
#
#     julia --project=. -t 16 demos/zetand_surfing_conventions.jl
#
# Every choice below is either stated in arXiv:2603.18144 or MEASURED from its Figure 10 by
# `demos/zetand_extract_fig10.py`. Where it was measured, the measurement is given, because
# these are the assumptions a comparison stands or falls on.
#
# 1. TEMPERATURE -> INTENSITY IS LINEAR, I ∝ T.
#    Not stated in the paper; measured from Figure 10, which contains both the temperature map
#    (top) and the H-band intensity images it produces (bottom). De-limb-darkening the intensity
#    panels gives a deepest contrast of 0.881 (1st percentile; min 0.847 over the eight phases),
#    while the temperature map's coldest point is 3742 K against a mean of 4372 K, a ratio of
#    0.856. I ∝ T predicts 0.856. A blackbody at 1.605 um predicts 0.679. The panels say linear.
#
#    This matters: at 1.605 um and 4372 K, hc/(lambda k T) = 2.05, so H band is NOT the
#    Rayleigh-Jeans limit for a K giant and dlnI/dlnT = 2.36, not 1. See `zetand_proper.jl`.
#
# 2. LIMB DARKENING alpha = 0.269, FIXED, not fitted — "after R. M. Roettenbacher et al. (2016)",
#    adopted to avoid "confusion between limb darkening and spots on the limb". Figure 10's own
#    intensity panels give alpha = 0.248 from their centre-to-limb profile, consistent with it.
#
# 3. NO SPATIAL PRIOR. The imaging is MCMC only: each map is "an average of 50 instances
#    (24 walkers over 500 iterations)". Positivity, and smoothness that is emergent from
#    averaging the posterior rather than imposed. See `zetand_ensemble_test.jl` for why
#    averaging OPTIMISATIONS does not reproduce that.
#
# 4. 768 SURFACE ELEMENTS, nside_exp = 3, "matching the CHARA Array's spatial resolution".
#
# 5. ERRORS INFLATED BY 1.65 to reach chi2 ~ 1, so a chi2 computed on the published error bars
#    must be divided by 1.65^2 = 2.72 before it is compared with theirs.
#
# 6. P = 17.769426 d and T0 = MJD 49992.281 (time of nodal passage), both FIXED. The map's
#    longitude zero is HALF A PERIOD from that T0: at their T0 the map fits at chi2/n = 19.36,
#    half a period away at 2.44.
#
# 7. FIGURE, their Table 6, 2026 column: theta_LD = 2.495 +- 0.014 mas, major-to-minor axis
#    ratio 1.063 +- 0.013, inclination 69.3 +- 6.9 deg, pole PA 124.0 +- 2.9 deg.

# 8. THE FIGURE IS AXISYMMETRIC (OBLATE), on the evidence of the ratio itself. Fitting the
#    OBLATE family (radius_x = radius_y) to these data with the size free returns a ratio of
#    1.0625 — their published 1.063 +- 0.013, to three digits. Fitting the TRIAXIAL family
#    instead gives equatorial/polar 1.0546 and chi2/n 2.064 against the oblate 5.812: a
#    Dchi2 of 53292 for ONE extra parameter over 14218 data. So the data reject an axisymmetric
#    figure outright, and they should — an oblate star presents the same outline at every
#    rotation phase and cannot produce an ellipsoidal light curve, which ζ And has.
#
# 9. THE FIGURE DOES NOT TRANSFER DIRECTLY, and this script does not pretend it does.
#    Their theta_LD = 2.495 +- 0.014 mas is a LIMB-DARKENED EQUIVALENT DISC diameter, divided
#    by a 1.0054 wavelength correction; ours is 2*cbrt(rx*ry*rz) = 2.5415 mas for the fitted
#    triaxial figure. The 1.9 % difference is not cosmetic at these baselines — see the table
#    this script prints. Nor is "major-to-minor axis ratio 1.063" unambiguous: read as an
#    oblate figure (rx = ry) it is the shape family that CANNOT represent a tidally distorted
#    star, and it fits far worse than read as prolate.

using ROTIR, Printf, DelimitedFiles, Statistics

const PROT     = 17.769426
const THEIR_T0 = 49992.281
const ALPHA    = 0.269          # fixed, after Roettenbacher et al. 2016
const INFLATE  = 1.65           # their error-bar inflation
const ZET      = joinpath(pkgdir(ROTIR), "demos", "zetAnd")

function main()
    files = sort(filter(f -> endswith(f, ".oifits"), readdir(ZET)))
    data  = [readoifits(joinpath(ZET, f); verbose = false)[1] for f in files]
    mjd   = [minimum(d.v2_mjd) for d in data]
    nd    = sum(d.nv2 + d.nt3amp + d.nt3phi for d in data)
    tess  = tessellation_healpix(3; T = Float64)
    @printf("their T0 = MJD %.3f; our first epoch is at phase %.4f of it\n",
            THEIR_T0, mod((minimum(mjd) - THEIR_T0)/PROT, 1))

    pub  = vec(readdlm(joinpath(ZET, "zetAnd2019A_map.txt")))
    flip = vec(Int.(readdlm(joinpath(ZET, "healpix_flip_3.txt"))))
    # THEIR map is a temperature map normalised by 4600 K — the scale of Figure 10's colourbar,
    # whose 3500-4600 K range it reproduces when multiplied out.
    Tmap = pub[flip] .* 4600.0
    @printf("their map: %.0f - %.0f K (mean %.0f)\n\n", minimum(Tmap), maximum(Tmap), mean(Tmap))

    function chi2(rx, ry, rz, inc, pa, ld, ph; x = Tmap)
        p = (surface_type = 1, radius_x = rx, radius_y = ry, radius_z = rz, tpole = 4372.0,
             ldtype = 3, ld1 = ld, ld2 = 0.0, beta = 0.25, inclination = inc,
             position_angle = pa, rotation_period = PROT)
        stars = create_star_multiepochs(tess, p, (mjd .- minimum(mjd)) .- ph*PROT)
        c = 0.0
        for (i, s) in enumerate(stars)
            # THEIR forward model: intensity proportional to temperature, i.e. ROTIR's :linear
            # handed the temperature map directly. Wrong physically, right for reproducing them.
            v2m, t3am, t3pm = ROTIR.cvis_to_obs(ROTIR.fused_cvis(x, s, data[i];
                                                 intensity_model = :linear), data[i])
            c += sum(abs2, (v2m .- data[i].v2) ./ data[i].v2_err)
            c += sum(abs2, (t3am .- data[i].t3amp) ./ data[i].t3amp_err)
            c += sum(abs2, ROTIR.mod360.(t3pm .- data[i].t3phi) ./ data[i].t3phi_err)
        end
        c/nd
    end

    r, th = 1.063, 2.495
    rzO = th/2/cbrt(r^2); rzP = th/2/cbrt(r)
    @printf("%-32s %10s %10s\n", "figure", "a=0.269", "a=0.145")
    for (tag, rx, ry, rz) in (("their ratio read as OBLATE",  r*rzO, r*rzO, rzO),
                              ("their ratio read as PROLATE", r*rzP, rzP,   rzP),
                              ("our fitted triaxial figure",  1.3332, 1.2545, 1.2269))
        @printf("%-32s %10.3f %10.3f\n", tag,
                chi2(rx, ry, rz, 69.3, 124.0, 0.269, 0.5077),
                chi2(rx, ry, rz, 69.3, 124.0, 0.145, 0.5077))
    end
    @printf("\nmean diameter: theirs %.4f mas (equivalent LD disc), ours %.4f mas\n", th,
            2*cbrt(1.3332*1.2545*1.2269))
    best = chi2(1.3332, 1.2545, 1.2269, 69.3, 124.0, 0.145, 0.5077)
    @printf("\nbest of the above on published errors: chi2/n = %.3f;  on their inflated ones: %.3f\n",
            best, best/INFLATE^2)
    @printf("phase at THEIR T0 instead of half a period away: chi2/n = %.3f\n",
            chi2(1.3332, 1.2545, 1.2269, 69.3, 124.0, 0.145,
                 mod((minimum(mjd) - THEIR_T0)/PROT, 1)))
    @printf("uniform surface, same figure:                    chi2/n = %.3f\n",
            chi2(1.3332, 1.2545, 1.2269, 69.3, 124.0, 0.145, 0.5077; x = fill(4372.0, tess.npix)))
end

main()
