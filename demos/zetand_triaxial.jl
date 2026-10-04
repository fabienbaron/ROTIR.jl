#!/usr/bin/env julia
# ζ And: the FIGURE — three semi-axes, limb darkening, six epochs over 1.13 rotations.
#
#     julia --project=. demos/zetand_triaxial.jl
#
# WHY TRIAXIAL AND NOT OBLATE. ζ And is a synchronised RS CVn binary, so its distortion is
# TIDAL: the K giant is elongated towards its companion, radius_x > radius_y ≈ radius_z. An
# oblate spheroid (radius_x = radius_y, squashed along the spin axis) cannot represent that at
# any parameter value — asked for one it drives the shape to its bound, which is exactly what
# the previous version of `fit_ellipsoid_ld` did here: a flattening of 0.3997 against a wall at
# 0.4, where a rigid rotator at break-up reaches only 1/3.
#
# THE PERIOD IS ADOPTED, NOT FITTED. P = 17.769426 d, the orbital period, which the system is
# tidally locked to (Roettenbacher et al. 2016). It is required: a triaxial figure presents a
# different projected ellipse at every rotation phase, and that is precisely what makes the
# three radii separately identifiable from six epochs that span 1.13 turns.
#
# WEIGHTS (1, 1, 0) — T3φ excluded, as in `zetand_sphere_ld.jl`: a uniform ellipsoid is still
# centrosymmetric in projection, so it predicts closure phases of 0° or 180° only. The measured
# T3φ is the spot signal and belongs to the imaging step.

using ROTIR, Printf

const PROT = 17.769426        # d, tidally locked
const ZETDIR = joinpath(pkgdir(ROTIR), "demos", "zetAnd")
# Roettenbacher et al.'s single "major-to-minor axis ratio". NOT an equatorial:polar ratio,
# and not a ground truth: it is another group's result under their own model, and ζ And's real
# figure is unknowable. Their model carries ONE ratio, so its two shorter axes are equal — and
# `zetand_fig10_match.jl` confirms that from their own Figure 10, which draws 1.0629 : 1.0008 : 1.
# A ROTIR triaxial fit has three, so `eqpol` below is NOT the same quantity as their 1.063 and
# the two should not be read off against each other. The like-for-like comparison is in
# `zetand_fig10_family.jl`, which fits THEIR family and gets 1.0626 against their 1.063.
const LIT_EQ_POL = 1.06

function zetdata()
    files = sort(filter(f -> endswith(f, ".oifits"), readdir(ZETDIR)))
    d = [readoifits(joinpath(ZETDIR, f); verbose = false)[1] for f in files]
    m = [minimum(x.v2_mjd) for x in d]
    (d, m .- minimum(m))
end

eqpol(r) = (r.radius_x + r.radius_y) / (2r.radius_z)

function show_fit(tag, r, nd, el)
    @printf("%-10s rx=%.4f ry=%.4f rz=%.4f mas   inc=%.1f°  PA=%.1f°  ld1=%.3f\n",
            tag, r.radius_x, r.radius_y, r.radius_z, r.inclination, r.position_angle, r.ld1)
    @printf("%-10s equatorial elongation %.4f (tidal) | polar flattening %.4f (rotational)\n",
            "", r.equatorial_elongation, r.polar_flattening)
    @printf("%-10s equatorial/polar = %.4f   major/minor = %.4f   mean D = %.4f mas\n",
            "", eqpol(r), r.radius_x/min(r.radius_y, r.radius_z), 2r.mean_radius)
    @printf("%-10s chi2/n = %.3f   (sphere %.3f)   %.0f s\n", "", r.chi2/nd, r.chi2_sphere/nd, el)
end

function main()
    data, ts = zetdata()
    nd = sum(d.nv2 + d.nt3amp for d in data)
    @printf("%d epochs, %d V2 + T3amp, phases %s\n\n", length(data), nd,
            join([@sprintf("%.2f", mod(t/PROT, 1)) for t in ts], " "))

    # ---- the measurement, at the two meshes that bracket the resolution ---------------------
    # λ/2B_max = 0.4717 mas against D = 2.56 mas is 5.4 resolution elements across the disc, so
    # nside 3 (768 tessels, 2.9 per element) is the coarsest defensible imaging grid and nside 2
    # (1.4 per element) is below Nyquist — included to show the discretisation, not to be used.
    fits = Dict{Int,Any}()
    for nexp in (2, 3)
        tess = tessellation_healpix(nexp)
        el = @elapsed r = fit_ellipsoid_ld(data, tess; tepochs = ts, rotation_period = PROT,
                                           ldtype = 3, r_bounds = (0.5, 3.0))
        fits[nexp] = r
        show_fit("nside $nexp", r, nd, el)
    end
    @printf("\nmesh-to-mesh difference in equatorial/polar: %.2f %%\n\n",
            100*(eqpol(fits[3])/eqpol(fits[2]) - 1))

    # ---- how much of that shape could the SPOTS be producing? -------------------------------
    # Synthetic data from a KNOWN figure carrying three dark spots at the 20 % contrast of the
    # published 2019A map, on the real uv coverage, noiseless — so the spots are the only
    # perturbation and the bias they induce is measured rather than argued.
    tess = tessellation_healpix(3)
    truth = (surface_type = 1, radius_x = 1.350, radius_y = 1.310, radius_z = 1.255,
             tpole = 5000.0, ldtype = 3, ld1 = 0.20, ld2 = 0.0, beta = 0.25,
             inclination = 65.0, position_angle = 40.0, rotation_period = PROT)
    # [:, 5, :] IS THE TESSEL CENTRE; 1:4 are its corners. Using index 1 puts every spot
    # 12 deg from where this code says it is — about 1.7 facet widths at nside_exp 3.
    ctr = tess.unit_xyz[:, 5, :]
    spot(lon, lat, σ) = begin
        n = [cosd(lat)*cosd(lon), cosd(lat)*sind(lon), sind(lat)]
        exp.(-([acosd(clamp(sum(ctr[i,:] .* n), -1, 1)) for i in 1:size(ctr,1)] ./ σ).^2)
    end
    bright = Dict("uniform" => ones(tess.npix),
                  "20% spots" => 1 .- 0.20 .* (spot(30.0, 40.0, 25.0) .+
                                               spot(160.0, -20.0, 30.0) .+
                                               spot(260.0, 10.0, 20.0)))
    @printf("synthetic truth: equatorial/polar = %.4f\n",
            (truth.radius_x + truth.radius_y)/(2truth.radius_z))
    for tag in ("uniform", "20% spots")
        sdata, sts = zetdata()
        stars = create_star_multiepochs(tess, truth, sts)
        x = truth.tpole .* bright[tag]
        for (i, s) in enumerate(stars)
            v2m, t3am, t3pm = ROTIR.cvis_to_obs(ROTIR.fused_cvis(x, s, sdata[i]), sdata[i])
            sdata[i].v2 .= v2m; sdata[i].t3amp .= t3am; sdata[i].t3phi .= t3pm
        end
        r = fit_ellipsoid_ld(sdata, tess; tepochs = sts, rotation_period = PROT, ldtype = 3,
                             r_bounds = (0.5, 3.0))
        @printf("  surface %-10s -> equatorial/polar = %.4f  (%+.2f %%)   inc=%.1f° PA=%.1f°\n",
                tag, eqpol(r),
                100*(eqpol(r)/((truth.radius_x + truth.radius_y)/(2truth.radius_z)) - 1),
                r.inclination, r.position_angle)
    end

    @printf("\nSo: this analysis gives equatorial/polar %.3f with a spot-induced bias of +1.1 %%\n",
            eqpol(fits[3]))
    @printf("measured on the synthetic star above. Roettenbacher et al.'s %.3f is a MAJOR-TO-MINOR\n",
            LIT_EQ_POL)
    @printf("ratio of a two-axis figure and is not this number; fitting their family instead gives\n")
    @printf("1.0626 against their 1.063 (`zetand_fig10_family.jl`), which is the comparison that\n")
    @printf("means something — and it is not a measurement hitting a truth, which is unknowable\n")
    @printf("here. The residual chi2/n of %.2f is the spot signal, and that is what the imaging\n",
            fits[3].chi2/nd)
    @printf("step is for.\n")
end

main()
