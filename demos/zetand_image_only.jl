#!/usr/bin/env julia
# ζ And: IMAGE-ONLY reconstruction with the figure held fixed.
#
#     julia --project=. -t 16 demos/zetand_image_only.jl
#
# The first half of a joint model + imaging fit: hold the triaxial figure at the value the
# uniform-surface fit found, and solve only for the 768-tessel surface. What this establishes,
# before any joint fit is attempted, is how much of ζ And's χ² the SURFACE alone can explain —
# and whether it beats the externally published map, which reaches χ²/n = 1.98 on this figure.
#
# CLOSURE PHASES ARE IN, and that is the point. `spheroid_chi2_fg` sums V², T3amp and T3φ
# unconditionally, so the imaging path has always used them; it is the SHAPE fitters that
# default to `weights = [1,1,0]`, because a uniform ellipsoid is centrosymmetric and predicts
# T3φ ≡ 0° or 180°. A spotted surface is not, and T3φ is where the spots are.
#
# THE REGULARISER IS `sobel`. Measured on RW Cep at nside 3 (docs/src/guides/reconstruction.md):
# tv2 at its best weight fits WORSE than no prior at all (χ²v2/n 2.29 against 2.00) because its
# k⁴ response flattens structure the closure phases can see, while sobel reaches 1.86 and the
# best χ²t3φ. Its weight is O(1–10²), not tv2's O(1e-5), since the Sobel forms divide by
# mean(x)² and carry the solid angle — so the ladder below is run in that range.
#
# THE OVERALL SCALE OF THE MAP IS NOT DETERMINED, and this is exact rather than a convergence
# issue: the visibilities are flux-normalised, so χ² is invariant under x → c·x, and `sobel`
# divides by mean(x), so it is invariant too. The whole criterion therefore has a flat
# direction. VMLMB simply does not move along it and the scale is inherited from the start, so
# everything reported here is normalised to mean(x) = 1 and only CONTRAST is meaningful.

using ROTIR, Printf, DelimitedFiles, Statistics

const PROT   = 17.769426
# THE MAP'S PHASE ZERO, AND IT IS HALF A PERIOD FROM THE PUBLISHED T0. Roettenbacher et al.
# give a time of nodal passage T0 = MJD 49992.281, and our first epoch sits at phase 0.9923 of
# it — i.e. essentially AT it. Using that directly is catastrophic: chi2/n = 19.36 against 2.44
# half a period away. So the map's longitude zero and the nodal passage differ by half a turn,
# which is a convention difference and not an error in either.
#
# 0.5077 is "half a period after their T0" expressed in this script's frame. An earlier value
# of 0.511 came from fitting the phase WITHOUT closure phases; with T3phi in, the minimum is
# sharp enough (within +0.1 in chi2/n over 0.010 in phase) that the difference matters.
const PHASE0 = 0.5077
const ZET    = joinpath(pkgdir(ROTIR), "demos", "zetAnd")

# the figure from demos/zetand_triaxial.jl, held FIXED here
const FIGURE = (surface_type = 1, radius_x = 1.3275, radius_y = 1.2663, radius_z = 1.2120,
                tpole = 5000.0, ldtype = 3, ld1 = 0.139, ld2 = 0.0, beta = 0.25,
                inclination = 67.8, position_angle = 126.2, rotation_period = PROT)

function setup()
    files = sort(filter(f -> endswith(f, ".oifits"), readdir(ZET)))
    data  = [readoifits(joinpath(ZET, f); verbose = false)[1] for f in files]
    mjd   = [minimum(d.v2_mjd) for d in data]
    ts    = (mjd .- minimum(mjd)) .- PHASE0*PROT      # same frame as the published map
    tess  = tessellation_healpix(3; T = Float64)
    stars = create_star_multiepochs(tess, FIGURE, ts)
    setup_oi!(data, stars)                            # the dense path image_reconstruct_oi uses
    (data, stars, tess)
end

"χ² per datum, split by observable — the split is what says whether the spots are being fitted."
function report(tag, x, stars, data)
    nv2 = sum(d.nv2 for d in data); nt3a = sum(d.nt3amp for d in data)
    nt3p = sum(d.nt3phi for d in data)
    c = [ROTIR.chi2_breakdown(x, stars[i], data[i]) for i in eachindex(data)]
    v2 = sum(ci.v2 for ci in c); t3a = sum(ci.t3amp for ci in c); t3p = sum(ci.t3phi for ci in c)
    n = nv2 + nt3a + nt3p
    @printf("%-26s χ²/n %6.3f | V² %6.3f  T3amp %6.3f  T3φ %6.3f | contrast %.3f–%.3f\n",
            tag, (v2+t3a+t3p)/n, v2/nv2, t3a/nt3a, t3p/nt3p,
            minimum(x)/mean(x), maximum(x)/mean(x))
    (v2+t3a+t3p)/n
end

function main()
    data, stars, tess = setup()
    npix = stars[1].npix
    nd = sum(d.nv2 + d.nt3amp + d.nt3phi for d in data)
    seen = falses(npix); for s in stars; seen[s.index_quads_visible] .= true; end
    @printf("%d epochs, %d data (V²+T3amp+T3φ), %d tessels, %d NEVER visible (%.1f %%)\n\n",
            length(data), nd, npix, count(.!seen), 100*count(.!seen)/npix)

    pub = vec(readdlm(joinpath(ZET, "zetAnd2019A_map.txt")))
    flip = vec(Int.(readdlm(joinpath(ZET, "healpix_flip_3.txt"))))
    xpub = 5000.0 .* (pub[flip] ./ maximum(pub))
    uniform = fill(5000.0, npix)
    report("uniform surface", uniform, stars, data)
    report("published map (as given)", xpub, stars, data)
    println()

    sob = sobel_gradient_healpix(3; T = Float64)
    regs(w) = w == 0 ? Any[] : Any[["sobel", w, sob, 1:npix]]

    # CONVERGENCE FIRST. The repo measures 0.14–0.66 drift in χ²/n between 300 and 3000
    # iterations on RW Cep — larger than the differences a λ ladder is trying to resolve, so a
    # λ comparison made at an arbitrary maxiter compares stopping points, not priors.
    println("convergence ladder at sobel weight 10, cold start:")
    prev = Inf
    for it in (300, 1000, 3000)
        x = image_reconstruct_oi(copy(uniform), data, stars; regularizers = regs(10.0),
                                 maxiter = it, verbose = false)
        c = report(@sprintf("  maxiter %5d", it), x, stars, data)
        @printf("%-26s Δχ²/n vs previous: %+.4f\n", "", c - prev); prev = c
    end
    println()

    # THE LADDER, AND WHY χ² CANNOT PICK THE WEIGHT. χ² falls monotonically as the prior is
    # weakened — at weight 0 it reaches 0.99, BELOW one, with the map running 0 to 8 times its
    # own mean and pixels pinned at the lower bound. 741 visible tessels against 14218 data at
    # 2.9 pixels per resolution element is enough freedom to absorb the noise, so "minimise χ²"
    # selects the least physical map on the ladder.
    #
    # Two criteria that do work are reported beside it: the CONTRAST, which a stellar surface
    # cannot have at 0–8, and the CORRELATION with an independent published reconstruction of
    # the same data. They agree with each other and disagree with χ².
    println("sobel weight ladder at maxiter 1000 (cold start):")
    best = (-Inf, 0.0, nothing)
    for w in (0.0, 10.0, 100.0, 300.0, 1e3, 3e3, 1e4, 3e4, 1e5)
        x = image_reconstruct_oi(copy(uniform), data, stars; regularizers = regs(w),
                                 maxiter = 1000, verbose = false)
        r = cor(x[seen] ./ mean(x[seen]), xpub[seen] ./ mean(xpub[seen]))
        c = report(@sprintf("  weight %7.0f  r=%.3f", w, r), x, stars, data)
        r > best[1] && (best = (r, w, x))
    end
    @printf("\nbest agreement with the published map: r = %.3f at sobel weight %.0f\n\n",
            best[1], best[2])

    # WARM START FROM THE PUBLISHED MAP. Two things this tests: that the solver cannot make the
    # published map worse (if it does, the criterion disagrees with the one that produced it),
    # and whether a cold start finds the same basin.
    @printf("warm start from the published map, sobel weight %.0f:\n", best[2])
    xw = image_reconstruct_oi(copy(xpub), data, stars; regularizers = regs(best[2]),
                              maxiter = 1000, verbose = false)
    cw = report("  warm", xw, stars, data)
    xc = best[3]
    r = cor(xc[seen] ./ mean(xc[seen]), xw[seen] ./ mean(xw[seen]))
    @printf("\ncold vs warm at that weight: correlation over visible tessels %.4f (χ²/n %.3f)\n",
            r, cw)
    @printf("cold vs published map: correlation over visible tessels %.4f\n", best[1])

    out = joinpath(ZET, "zetAnd2019A_image_only.txt")
    writedlm(out, xc ./ mean(xc))      # normalised: the scale is not determined, see the header
    println("wrote ", out)
end

main()
