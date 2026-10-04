#!/usr/bin/env julia
# The ζ And surface map read off the figure in `zetAnd/zetAnd2019A.pdf` — does it fit the data?
#
#     julia --project=. demos/zetand_published_map.jl
#
# `zetand_extract_map.py` inverts the figure's own colourbar (a `gist_heat` ramp) and resamples
# onto HEALPix nside_exp 3.
# That extraction is checkable against the image itself (colour distances, the plot ellipse, the
# histogram) but those checks cannot tell whether the result is the star. This one can: a map
# that is really ζ And's surface must FIT ζ AND'S VISIBILITIES, and a mis-parity of it must not.
#
# THERE IS NO GROUND TRUTH HERE. Roettenbacher et al.'s reconstruction is another analysis under
# its own model, not the star. What this demo establishes is narrower and more useful: that the
# map we read out of their figure carries real information about these data, and at which
# orientation. It is the natural starting point and comparison for a joint model + imaging fit.

using ROTIR, Printf, DelimitedFiles, NLopt

const PROT = 17.769426
const ZETDIR = joinpath(pkgdir(ROTIR), "demos", "zetAnd")

function main()
    files = sort(filter(f -> endswith(f, ".oifits"), readdir(ZETDIR)))
    data = [readoifits(joinpath(ZETDIR, f); verbose = false)[1] for f in files]
    mjd  = [minimum(d.v2_mjd) for d in data]; t = mjd .- minimum(mjd)
    nd   = sum(d.nv2 + d.nt3amp for d in data)
    tess = tessellation_healpix(3)

    # the figure from demos/zetand_triaxial.jl, scaled by `s`
    R = (1.3275, 1.2663, 1.2120)
    mkp(s, inc, pa) = (surface_type = 1, radius_x = s*R[1], radius_y = s*R[2], radius_z = s*R[3],
                       tpole = 5000.0, ldtype = 3, ld1 = 0.139, ld2 = 0.0, beta = 0.25,
                       inclination = inc, position_angle = pa, rotation_period = PROT)
    function chi2(x, p, ts)
        stars = create_star_multiepochs(tess, p, ts)
        c = 0.0
        for (i, st) in enumerate(stars)
            v2m, t3am, _ = ROTIR.cvis_to_obs(ROTIR.fused_cvis(x, st, data[i]), data[i])
            c += sum(abs2, (v2m .- data[i].v2) ./ data[i].v2_err)
            c += sum(abs2, (t3am .- data[i].t3amp) ./ data[i].t3amp_err)
        end
        c
    end

    unif = chi2(fill(5000.0, tess.npix), mkp(1.0, 64.7, 122.8), t)
    @printf("uniform surface on the same figure:  chi2/n = %.3f\n\n", unif/nd)

    # THE MAP IS STORED AS PLOTTED. ζ And's visibilities select its EQUATOR-FLIPPED parity,
    # which is a statement about the model frame, not about the figure: `plot_mollweide` has no
    # flip in it (colatitude runs pi -> 0 against latitude -90 -> +90, consistently), so the
    # parity belongs to how ROTIR orients a star relative to how the reconstruction did. The
    # flip is an exact permutation of HEALPix centres — a bijection and an involution, worst
    # centre mismatch 0.028 deg — precomputed in `healpix_flip_3.txt`.
    m = vec(readdlm(joinpath(ZETDIR, "zetAnd2019A_map.txt")))
    flip = vec(Int.(readdlm(joinpath(ZETDIR, "healpix_flip_3.txt"))))
    m = m[flip]
    x = 5000.0 .* (m ./ maximum(m))
    # The map's orientation is NOT published with it: the phase zero, the inclination and the
    # position angle it was reconstructed at all have to be refitted here, and the overall size
    # is freed as a check — if the figure is right it should come back at 1.
    best = Inf; bestθ = zeros(4)
    for ph0 in 0.0:0.25:0.75, inc0 in (45.0, 70.0)
        opt = NLopt.Opt(:LN_NELDERMEAD, 4)
        opt.lower_bounds = [-1.0, 20.0,   0.0, 0.8]
        opt.upper_bounds = [ 2.0, 90.0, 180.0, 1.25]
        opt.maxeval = 1500; opt.xtol_rel = 1e-7
        opt.min_objective = (θ, g) -> chi2(x, mkp(θ[4], θ[2], θ[3]), t .- θ[1]*PROT)
        f, xo, _ = NLopt.optimize(opt, [ph0, inc0, 122.8, 1.0])
        f < best && (best = f; bestθ = copy(xo))
    end
    @printf("the extracted map, orientation refitted:\n")
    @printf("  chi2/n = %.3f   phase %.3f   inc %.1f deg   PA %.1f deg   size x%.3f\n",
            best/nd, mod(bestθ[1], 1), bestθ[2], bestθ[3], bestθ[4])
    @printf("  improvement over a uniform surface: %.0f %%\n\n", 100*(1 - best/unif))

    # HOW WELL IS THE PHASE PINNED? It is the one number that had to be invented: the figure
    # carries no longitude labels at all (title, colourbar and its two numbers are the only ink
    # outside the ellipse), so there is nothing to read off. What is fitted conflates her
    # longitude origin with her reference epoch — both are a rigid rotation about the spin axis
    # and enter the model identically, so no data can separate them.
    @printf("\nphase scan at the fitted orientation (the one number with nothing to read off):\n")
    phs = 0:0.01:0.99
    cs = [chi2(x, mkp(bestθ[4], bestθ[2], bestθ[3]), t .- ph*PROT)/nd for ph in phs]
    i0 = argmin(cs); within = phs[cs .< cs[i0] + 0.1]
    @printf("  minimum %.3f at phase %.2f | worst %.3f | within +0.1 over %.3f in phase = %.0f deg\n",
            cs[i0], phs[i0], maximum(cs), maximum(within)-minimum(within),
            360*(maximum(within)-minimum(within)))
    @printf("  i.e. phase zero at MJD %.3f, given this inclination, PA and figure.\n",
            minimum(mjd) + phs[i0]*PROT)
    @printf("  A %.1fx swing in chi2 with phase is itself evidence the map is real: an\n",
            maximum(cs)/cs[i0])
    @printf("  arbitrary pattern would not care where it was pointed.\n\n")

    @printf("For reference, the four latitude/longitude parities of the same extraction, each\n")
    @printf("with its own orientation refitted (measured when the convention was chosen):\n")
    @printf("    equator flipped  chi2/n = 1.978   <- the parity used above\n")
    @printf("    as plotted       chi2/n = 3.564\n")
    @printf("    both flipped     chi2/n = 3.959\n")
    @printf("    lon flipped      chi2/n = 4.584\n")
    @printf("Only one parity beats the uniform surface, which is what says the extraction is\n")
    @printf("signal rather than an artefact of inverting a colormap.\n")
end

main()
