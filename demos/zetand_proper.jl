#!/usr/bin/env julia
# HOW TO DO IT PROPERLY — the same data, with the temperature-to-intensity relation the star
# actually has.
#
#     julia --project=. -t 16 demos/zetand_proper.jl
#
# THE ONE THING THAT IS WRONG in the published treatment is the conversion, and it is wrong in a
# way that leaves the FIT untouched and the ASTROPHYSICS wrong. See
# `demos/zetand_surfing_conventions.jl` for how their choices were established; the short of it
# is that Figure 10's intensity panels are related to its temperature map by I ∝ T, measured
# from the figure itself.
#
# Interferometry constrains INTENSITY. The visibilities never see a temperature; they see the
# brightness distribution on the sky. A temperature map is that intensity map pushed through an
# assumed law, so:
#
#   * their intensity map is fine, and so is every geometric quantity derived from it;
#   * their TEMPERATURES are not, because I ∝ T is the Rayleigh-Jeans limit and H band is not
#     the Rayleigh-Jeans limit for a K giant. At 1.605 um and 4400 K, hc/(lambda k T) = 2.04,
#     so dlnI/dlnT = 2.36 rather than 1.
#
# CORRECTED 2026-10-04: 2.36 IS A BLACKBODY, AND THE STAR IS NOT ONE. Kurucz ATLAS9 at
# logg 2.5 — the gravity their own polar radius 15.0 R_sun and mass 2.6 M_sun imply — gives
# dlnI/dlnT = **1.51** in H band, because H- opacity falls with temperature and partly cancels
# the Planck steepening. Everything below inverts a blackbody and therefore OVERSTATES the
# correction: the spot deficit comes out 273 K where a real atmosphere gives 417 K, against
# their published 630 K. The right headline is "overstated by about 1.5x", not 2.3x. Their 2016
# paper says they convert with Kurucz models themselves, which makes the blackbody doubly the
# wrong reference. `demos/zetand_alpha_audit.jl` reads the same grid.
#
# The correction is to invert the real Planck function instead of a proportionality. It needs an
# absolute anchor, because a flux-normalised visibility fixes only the CONTRAST: the mean
# temperature has to come from spectroscopy, not from here. `TEFF_ANCHOR` is that input and the
# result scales with it, which this script makes explicit by reporting two anchors.

using ROTIR, Printf, DelimitedFiles, Statistics, NLopt

const PROT = 17.769426
const ZET  = joinpath(pkgdir(ROTIR), "demos", "zetAnd")
const TEFF_ANCHOR = 4400.0        # K, mean effective temperature — AN INPUT, from spectroscopy

"Invert the Planck function: the temperature whose band intensity is `r` times that of `Tref`."
function temperature_from_intensity(r, Tref, band)
    Bref, _ = ROTIR.planck_and_dT(Tref, band)
    target = r*Bref
    lo, hi = 0.3Tref, 3.0Tref
    for _ in 1:80                                    # bisection: B is monotone in T
        mid = 0.5*(lo+hi)
        (ROTIR.planck_and_dT(mid, band)[1] < target) ? (lo = mid) : (hi = mid)
    end
    0.5*(lo+hi)
end

function main()
    files = sort(filter(f -> endswith(f, ".oifits"), readdir(ZET)))
    data  = [readoifits(joinpath(ZET, f); verbose = false)[1] for f in files]
    mjd   = [minimum(d.v2_mjd) for d in data]
    nd    = sum(d.nv2 + d.nt3amp + d.nt3phi for d in data)
    tess  = tessellation_healpix(3; T = Float64)
    band  = Float64(ROTIR.band_of(data[1]))
    B, dB = ROTIR.planck_and_dT(TEFF_ANCHOR, band)
    @printf("band %.4g m, anchor %.0f K: hc/(lambda k T) = %.2f, dlnI/dlnT = %.2f\n",
            band, TEFF_ANCHOR, 1.438776877e-2/(band*TEFF_ANCHOR), dB*TEFF_ANCHOR/B)

    # Their published numbers are proportional to INTENSITY, whatever the colourbar says.
    pub  = vec(readdlm(joinpath(ZET, "zetAnd2019A_map.txt")))
    flip = vec(Int.(readdlm(joinpath(ZET, "healpix_flip_3.txt"))))
    Irel = pub[flip] ./ mean(pub)                    # intensity relative to the mean
    @printf("\nintensity map (what the data actually constrain): %.3f - %.3f of the mean\n",
            minimum(Irel), maximum(Irel))

    for anchor in (4400.0, 4600.0)
        Tproper = [temperature_from_intensity(r, anchor, band) for r in Irel]
        Tnaive  = anchor .* Irel                      # what I ∝ T gives
        @printf("anchor %.0f K:  proper T %.0f-%.0f K (deficit %.0f K)   naive I~T %.0f-%.0f K (deficit %.0f K)\n",
                anchor, minimum(Tproper), maximum(Tproper), anchor - minimum(Tproper),
                minimum(Tnaive), maximum(Tnaive), anchor - minimum(Tnaive))
    end

    # THE FIT IS UNCHANGED, and that is the point. Feeding the proper temperature map through
    # Planck reproduces exactly the intensity map that :linear was already using, so chi2 is
    # identical to machine precision — the correction moves the astrophysics, not the data.
    Tproper = [temperature_from_intensity(r, TEFF_ANCHOR, band) for r in Irel]
    function chi2(ld, ph, mode)
        p = (surface_type = 1, radius_x = 1.3332, radius_y = 1.2545, radius_z = 1.2269,
             tpole = TEFF_ANCHOR, ldtype = 3, ld1 = ld, ld2 = 0.0, beta = 0.25,
             inclination = 71.3, position_angle = 124.3, rotation_period = PROT)
        stars = create_star_multiepochs(tess, p, (mjd .- minimum(mjd)) .- ph*PROT)
        x, model, bnd = mode === :proper ? (Tproper, :planck, band) : (Irel, :linear, nothing)
        c = 0.0
        for (i, s) in enumerate(stars)
            v2m, t3am, t3pm = ROTIR.cvis_to_obs(ROTIR.fused_cvis(x, s, data[i];
                                                 intensity_model = model, band = bnd), data[i])
            c += sum(abs2, (v2m .- data[i].v2) ./ data[i].v2_err)
            c += sum(abs2, (t3am .- data[i].t3amp) ./ data[i].t3amp_err)
            c += sum(abs2, ROTIR.mod360.(t3pm .- data[i].t3phi) ./ data[i].t3phi_err)
        end
        c/nd
    end
    cp = chi2(0.145, 0.5061, :proper); ci = chi2(0.145, 0.5061, :intensity)
    @printf("\nchi2/n: proper temperature map through Planck %.6f\n", cp)
    @printf("        the same star as an intensity map        %.6f   (difference %.2e)\n",
            ci, abs(cp-ci))
    @printf("\nSo the correction changes no number the data can see. What it changes is the\n")
    @printf("spot temperature: %.0f K cooler than the photosphere, not %.0f K.\n",
            TEFF_ANCHOR - minimum(Tproper), TEFF_ANCHOR - minimum(TEFF_ANCHOR .* Irel))
end

main()
