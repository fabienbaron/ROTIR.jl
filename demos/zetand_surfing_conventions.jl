#!/usr/bin/env julia
# HOW ROETTENBACHER ET AL. DID IT — their conventions, applied to these data, so that any
# comparison with their published numbers is like for like.
#
#     julia --project=. -t 16 demos/zetand_surfing_conventions.jl
#
# Every choice below is either stated in arXiv:2603.18144 or MEASURED from its Figure 10 by
# `demos/zetand_extract_fig10.py` and, for the geometry and the orientation, by the pair
# `demos/zetand_fig10_panels.py` + `demos/zetand_fig10_match.jl`. Where it was measured, the
# measurement is given, because these are the assumptions a comparison stands or falls on.
#
# 1. TEMPERATURE -> INTENSITY IS LINEAR, I ∝ T.
#    Not stated in the paper; measured from Figure 10, which contains both the temperature map
#    (top) and the H-band intensity images it produces (bottom). With the projection pinned
#    (`zetand_fig10_match.jl`), every published pixel can be assigned to the facet it shows and
#    the panel value regressed on the map, facet by facet across all eight phases:
#
#        I = T  (linear)           r = 0.9909   rms = 13.0 K
#        I = Planck(T, 1.605 um)   r = 0.9152   rms = 39.1 K
#        I = T^4                   r = 0.7795   rms = 60.7 K
#
#    This matters: at 1.605 um and 4372 K, hc/(lambda k T) = 2.05, so H band is NOT the
#    Rayleigh-Jeans limit for a K giant. See `zetand_proper.jl`.
#
#    CORRECTED 2026-10-04, AND THE CORRECTION IS LARGE. The comparison above uses a BLACKBODY,
#    for which dlnI/dlnT = 2.35 at 4372 K. A real H-band atmosphere is far shallower: Kurucz
#    ATLAS9 at logg 2.5 gives **dlnI/dlnT = 1.51**, because H- opacity falls with temperature
#    and partly cancels the Planck steepening. Repeating the regression against a Kurucz
#    intensity rather than a blackbody:
#
#        I = T       (linear)        r = 0.99091   rms 13.0 K
#        I = Kurucz(T, logg 2.5)     r = 0.97597   rms 21.1 K
#        I = Planck(T, 1.605 um)     r = 0.91519   rms 39.1 K
#
#    Linear still fits their figure best — the best power law is I ~ T^0.96 — but a blackbody
#    is the WRONG thing to call the correct answer, and the size of their error shrinks with it.
#    Reading their map as an intensity and inverting it properly:
#
#        through a Kurucz atmosphere : 3955-4516 K, deficit 417 K, overstated x1.51
#        through a blackbody         : 4098-4463 K, deficit 273 K, overstated x2.31
#
#    So the right statement is "overstated by about 1.5", not 2.3. Their 2016 paper also SAYS
#    they convert with Kurucz models ("for [Fe/H]=-0.25 and appropriate log g"), which the
#    figure's measured exponent of 0.96 does not match — that tension is unresolved and the
#    figure is the only evidence either way.
#
# 2. LIMB DARKENING alpha = 0.269, FIXED, not fitted — "after R. M. Roettenbacher et al. (2016)",
#    adopted to avoid "confusion between limb darkening and spots on the limb". Figure 10's own
#    panels are drawn with alpha = 0.285 best-fitting, and their 0.269 is indistinguishable from
#    it (r 0.99103 vs 0.99091) while 0.145 is excluded (0.9555). So the figure is internally
#    consistent with Table 6 — and the data still want ~0.15 in every shape family, including
#    theirs. See `zetand_fig10_family.jl`: that disagreement is real, not a convention error.
#
# 3. THEY DO HAVE SPATIAL REGULARISERS, and an earlier version of this file said they did not.
#    Roettenbacher et al. 2016's methods section (arXiv:1709.10107 page 17) states: "In addition
#    to minimizing the chi2 statistic, we also incorporated priors on each parameter and could
#    experiment with a variety of imaging regularizers, such at Total Variation (TV) or L2norm
#    of wavelet coefficients." The 2026 paper's page 2 likewise lists "input stellar parameters,
#    regularization, and other" constraints as what reduces the degeneracy.
#
#    WHAT NEITHER PAPER SAYS is which regulariser produced the published maps, or at what
#    weight. "Could experiment with" is the strongest claim available. So a comparison can say
#    that SURFING has TV and a wavelet L2 norm, and cannot say the published images used either.
#
#    The imaging is also an ENSEMBLE: each 2026 map is "an average of 50 instances (each running
#    24 walkers over 500 iterations)", so whatever prior is applied, smoothness also arrives
#    from averaging the posterior. See `zetand_ensemble_test.jl` for why averaging OPTIMISATIONS
#    does not reproduce that.
#
# 4. 768 SURFACE ELEMENTS, nside_exp = 3, "matching the CHARA Array's spatial resolution".
#
# 5. ERRORS INFLATED BY 1.65 to reach chi2 ~ 1, so a chi2 computed on the published error bars
#    must be divided by 1.65^2 = 2.72 before it is compared with theirs.
#
# 6. P = 17.769426 d and T0 = MJD 49992.281 (time of nodal passage, from the ellipsoidal
#    light curve), both FIXED. THEIR T0 IS THE OPPOSITE NODE to ROTIR's phase zero. Passing
#    MJD epochs directly with the surface parameter `t0` (see docs/src/guides/conventions.md):
#
#        t0 = T0          chi2/n = 19.34
#        t0 = T0 + P/2    chi2/n =  2.11      <- and a scan of the whole period agrees,
#                                                finding its optimum at exactly T0 + 0.500 P
#
#    So use `t0 = 49992.281 + 17.769426/2` and no hand-tuned phase. A nodal passage is where
#    the star sits on the line of nodes, i.e. the tidal axis lies across the line of sight —
#    maximum projected elongation, and a light-curve MAXIMUM for an ellipsoidal variable. That
#    is the same kind of instant ROTIR's psi = 0 is; the half period is only WHICH node, which
#    the ephemeris defines and this table measures once.
#
# 7. FIGURE, their Table 6, 2026 column: theta_LD = 2.495 +- 0.014 mas, major-to-minor axis
#    ratio 1.063 +- 0.013, inclination 69.3 +- 6.9 deg, pole PA 124.0 +- 2.9 deg.

# 8. THE FIGURE IS PROLATE, AND THIS IS MEASURED FROM THE FIGURE, NOT INFERRED.
#    The silhouettes of Figure 10's eight panels are PHASE-DEPENDENT — semi-major axis 251.8,
#    245.4, 238.8, 245.2 px at phi = 0.000, 0.125, 0.250, 0.375 with the semi-minor fixed at
#    236.9. A rotationally oblate spheroid projects to the same ellipse at every phase, so the
#    figure is not one. Fitting the orthographic shadow of an ellipsoid to all eight second-
#    moment matrices (`zetand_fig10_panels.py`) reproduces them to 0.1 px out of 240 and gives
#
#        rx : ry : rz = 1.0629 : 1.0008 : 1,   inclination 69.44 deg
#
#    i.e. a tidally elongated prolate spheroid with NO rotational flattening — which is what
#    SURFING's single "major-to-minor axis ratio" can express, and what the paper's own text
#    describes when it says the light curve is "dominated by the star's ellipsoidal shape".
#    zeta And is a synchronised RS CVn binary, so a tidal figure is expected.
#
#    CONFIRMED BY THE SOURCE, after this was first inferred from the figure. Roettenbacher et
#    al. 2016 (arXiv:1709.10107, the Nature paper this work's alpha and geometry come from, and
#    which `demos/zetAnd/` now carries) says it outright on its page 3: "The 2011 and 2013
#    datasets were separately imaged onto a **prolate ellipsoid**". Its Table 1 then labels the
#    parameter "Oblateness (major to polar axis) = 1.060 ± 0.011", which is where the confusion
#    comes from — the word is "oblateness" and the quantity is a prolate figure's major-to-POLAR
#    ratio. Measured off the 2026 figure: rx/rz = 1.0629.
#
#    An earlier version of this file read their ratio as OBLATE. That was wrong, and the body
#    of this script still prints both readings so the difference stays visible.
#
# 9. THEIR theta_LD IS THE MINOR (POLAR) AXIS, and it transfers exactly. The 2016 paper names
#    it in its own Table 1: "Angular **polar** diameter, theta_LD (mas) = 2.502 ± 0.008". The
#    2026 paper relabels the same quantity "Limb-darkened disk diameter", which is what made it
#    look like an equivalent-disc diameter.
#    From Figure 10's own x ticks (7 of them, 0.5 mas apart, 94.917 px) the drawn semi-axes are
#    2*rx = 2.6529, 2*ry = 2.4979, 2*rz = 2.4959 mas. Their Table 6 theta_LD = 2.495 +- 0.014
#    is the MINOR axis to 0.04 %, and the major axis is 1.063 times larger. Fitting their own
#    prolate family to these data returns r_minor = 1.2489 mas against the 1.2479 read off the
#    figure — 0.08 %. Two independent routes to the same ellipsoid.
#
# 10. ORIENTATION AND PHASE, both needed before any overlay.
#    Figure 10's panels are drawn with East LEFT and North DOWN: the y axis runs -1.5 at the top
#    to +1.5 at the bottom, so the vertical axis is inverted from the usual convention and a
#    naive overlay mirrors. Their phase also runs opposite to ROTIR's, and the sub-observer
#    longitude of the published map is 90 deg - 360 deg * phi. The map itself needs no parity
#    surgery: it is used exactly as the Mollweide plots it.

using ROTIR, Printf, DelimitedFiles, Statistics

const PROT     = 17.769426
const THEIR_T0 = 49992.281      # their time of NODAL passage, from the ellipsoidal light curve
# ROTIR's phase zero is the OTHER node: see item 6. Measured, then used as a definition.
const OUR_T0   = THEIR_T0 + PROT/2
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

    # `t0` is a SURFACE PARAMETER now, so the epochs go in as plain MJD and the phase zero is
    # stated rather than tuned. `OUR_T0` is their published nodal passage plus half a period,
    # because their ephemeris refers to the other node — see item 6 above and
    # `docs/src/guides/conventions.md`. The old spelling shifted `tepochs` by a fitted 0.5077.
    function chi2(rx, ry, rz, inc, pa, ld, t0 = OUR_T0; x = Tmap)
        p = (surface_type = 1, radius_x = rx, radius_y = ry, radius_z = rz, tpole = 4372.0,
             ldtype = 3, ld1 = ld, ld2 = 0.0, beta = 0.25, inclination = inc,
             position_angle = pa, rotation_period = PROT, t0 = t0)
        stars = create_star_multiepochs(tess, p, mjd)
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
                chi2(rx, ry, rz, 69.3, 124.0, 0.269),
                chi2(rx, ry, rz, 69.3, 124.0, 0.145))
    end
    @printf("\nmean diameter: theirs %.4f mas (equivalent LD disc), ours %.4f mas\n", th,
            2*cbrt(1.3332*1.2545*1.2269))
    best = chi2(1.3332, 1.2545, 1.2269, 69.3, 124.0, 0.145)
    @printf("\nbest of the above on published errors: chi2/n = %.3f;  on their inflated ones: %.3f\n",
            best, best/INFLATE^2)
    @printf("t0 at THEIR T0 instead of the other node:        chi2/n = %.3f\n",
            chi2(1.3332, 1.2545, 1.2269, 69.3, 124.0, 0.145, THEIR_T0))
    @printf("uniform surface, same figure:                    chi2/n = %.3f\n",
            chi2(1.3332, 1.2545, 1.2269, 69.3, 124.0, 0.145; x = fill(4372.0, tess.npix)))
end

main()
