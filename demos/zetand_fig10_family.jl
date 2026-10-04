#!/usr/bin/env julia
# ζ And: does the limb-darkening disagreement come from the shape model?
#
#     julia --project=. demos/zetand_fig10_family.jl
#
# `demos/zetand_fig10_match.jl` reads Roettenbacher et al.'s geometry straight off their own
# Figure 10 and finds it internally consistent with their Table 6 — a PROLATE spheroid, major
# to minor 1.0629, 2*r_minor = 2.496 mas, i = 69.44 deg, PA = 124.30 deg, drawn with a fixed
# power-law limb darkening alpha = 0.269.
#
# ROTIR's fit to the same interferometry prefers alpha ~ 0.14. Two explanations compete:
#
#   (a) their alpha is simply fixed at the wrong value, or
#   (b) ROTIR's ellipsoid has three independent semi-axes where SURFING's figure has one axis
#       ratio, and that extra freedom is what buys the low alpha — gravity/limb darkening and
#       geometry both remove flux from the limb, which is the degeneracy `gravity_darkening.jl`
#       is about.
#
# So this fits BOTH families to the same data: theirs (r_minor = r_polar exactly, which is what
# a single "major-to-minor axis ratio" means) and ours (three semi-axes). T3phi is excluded,
# weights (1, 1, 0), as in `zetand_triaxial.jl`: a uniform ellipsoid is centrosymmetric in
# projection and predicts closure phases of 0 or 180 deg only, so the measured T3phi is the spot
# signal and belongs to the imaging step.
#
# RESULT. Restricted to their own family the data still want alpha = 0.154, and with their
# geometry held completely fixed, 0.150. So (b) is ruled out: the disagreement is not bought by
# ROTIR's extra shape freedom.
#
# And the geometry agrees to a degree worth stating: fitting THEIR family to the data returns
# rx = 1.3271 and r_minor = 1.2489 mas against the 1.3265 and 1.2490 read off their figure —
# 0.05 % and 0.01 %. Two independent analyses, one from the data and one from a rasterised PDF
# page, landing on the same ellipsoid.
#
# CAVEAT, and it is not small: chi2/n is 3.6 at best, because a uniform ellipsoid cannot
# represent a spotted star. Every formal delta chi2 below inherits that misfit, so they rank
# models but should not be read as significances. The triaxial family's extra 4.5 % polar
# flattening in particular is partly real and partly spot signal — `zetand_triaxial.jl` measures
# a +1.1 % spot-induced bias on the axis ratio from synthetic data on this same uv coverage.

using ROTIR, Printf, NLopt, Statistics

const PROT = 17.769426
const ZET  = joinpath(pkgdir(ROTIR), "demos", "zetAnd")
const W    = (1.0, 1.0, 0.0)
# Figure 10's own geometry, as `zetand_fig10_match.jl` establishes it
const FX, FY, FZ, FI, FP = 2.6529/2, 2.4979/2, 2.4959/2, 69.44, 124.30

function zetdata()
    fs = sort(filter(f -> endswith(f, ".oifits"), readdir(ZET)))
    d  = [readoifits(joinpath(ZET, f); verbose = false)[1] for f in fs]
    m  = [minimum(x.v2_mjd) for x in d]
    (d, m .- minimum(m))
end

mk(rx, ry, rz, inc, pa, ld) =
    (surface_type = 1, radius_x = rx, radius_y = ry, radius_z = rz, tpole = 5000.0,
     ldtype = 3, ld1 = ld, ld2 = 0.0, beta = 0.25, inclination = inc,
     position_angle = pa, rotation_period = PROT)

function main()
    data, ts = zetdata()
    nd   = sum(d.nv2 + d.nt3amp for d in data)
    tess = tessellation_healpix(3)
    chi2(p) = parametric_chi2(p, tess, data, ts; weights = W, uniform = true)
    @printf("%d epochs, %d V2 + T3amp\n\n", length(data), nd)

    println("Figure 10's geometry held FIXED (nothing free but alpha)")
    @printf("  their fixed alpha = 0.269        chi2/n = %.3f\n", chi2(mk(FX,FY,FZ,FI,FP,0.269))/nd)
    best = (Inf, 0.0)
    for al in 0.0:0.005:1.0
        c = chi2(mk(FX,FY,FZ,FI,FP,al)); c < best[1] && (best = (c, al))
    end
    @printf("  alpha free                       chi2/n = %.3f at alpha = %.3f\n\n",
            best[1]/nd, best[2])

    function fit(n, lo, hi, x0, build)
        o = Opt(:LN_NELDERMEAD, n); o.lower_bounds = lo; o.upper_bounds = hi
        o.xtol_rel = 1e-7; o.maxeval = 4000
        o.min_objective = (x, g) -> chi2(build(x))
        bf, bx = Inf, x0
        for s in 1:8
            st = x0 .+ (s == 1 ? zeros(n) : 0.04 .* (hi .- lo) .* (rand(n) .- 0.5))
            f, xo, _ = optimize(o, clamp.(st, lo, hi))
            f < bf && ((bf, bx) = (f, xo))
        end
        (bf, bx)
    end

    println("SURFING's family: prolate spheroid, r_minor = r_polar")
    f1, x1 = fit(4, [1.0,1.0,30.0,90.0], [1.8,1.8,89.0,170.0], [1.33,1.25,69.0,124.0],
                 x -> mk(x[1], x[2], x[2], x[3], x[4], 0.269))
    @printf("  alpha FIXED 0.269  rx=%.4f r_minor=%.4f  ratio=%.4f  inc=%.1f PA=%.1f   chi2/n = %.3f\n",
            x1[1], x1[2], x1[1]/x1[2], x1[3], x1[4], f1/nd)
    f2, x2 = fit(5, [1.0,1.0,30.0,90.0,0.0], [1.8,1.8,89.0,170.0,1.5], [1.33,1.25,69.0,124.0,0.27],
                 x -> mk(x[1], x[2], x[2], x[3], x[4], x[5]))
    @printf("  alpha free         rx=%.4f r_minor=%.4f  ratio=%.4f  inc=%.1f PA=%.1f   chi2/n = %.3f  alpha=%.3f\n",
            x2[1], x2[2], x2[1]/x2[2], x2[3], x2[4], f2/nd, x2[5])
    @printf("  Figure 10 itself   rx=%.4f r_minor=%.4f  ratio=%.4f  inc=%.1f PA=%.1f\n\n",
            FX, FZ, FX/FZ, FI, FP)

    println("ROTIR's family: three independent semi-axes")
    f3, x3 = fit(6, [1.0,1.0,1.0,30.0,90.0,0.0], [1.8,1.8,1.8,89.0,170.0,1.5],
                 [1.33,1.27,1.21,65.0,123.0,0.14], x -> mk(x[1],x[2],x[3],x[4],x[5],x[6]))
    @printf("  alpha free         rx=%.4f ry=%.4f rz=%.4f  ry/rz=%.4f  inc=%.1f PA=%.1f  chi2/n = %.3f  alpha=%.3f\n",
            x3[1], x3[2], x3[3], x3[2]/x3[3], x3[4], x3[5], f3/nd, x3[6])

    println("\nso:")
    @printf("  the data want alpha ~ 0.15 in EVERY family tested, including theirs (%.3f, %.3f, %.3f),\n",
            best[2], x2[5], x3[6])
    println("  so the disagreement with their fixed 0.269 is not an artefact of shape freedom;")
    @printf("  and fitting their family recovers their own figure to %.2f %% in rx and %.2f %% in r_minor.\n",
            100*abs(x2[1]/FX - 1), 100*abs(x2[2]/FZ - 1))
end

main()
