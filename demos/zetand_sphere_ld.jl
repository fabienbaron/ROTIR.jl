#!/usr/bin/env julia
# ζ And: the symmetric model — angular diameter and limb darkening — before any imaging.
#
#     julia --project=. demos/zetand_sphere_ld.jl
#
# WHY THIS FIRST. Before a joint model+imaging fit means anything, the symmetric part has to
# be pinned and checked against the published analysis. This fits ONLY radius and the
# limb-darkening exponent, on the six 2019 MIRC-X epochs.
#
# WEIGHTS (1, 1, 0) — T3φ EXCLUDED, and for ζ And that is not a detail. A centrosymmetric
# model has closure phases identically 0° or 180°, so the measured T3φ cannot be improved by
# any (diameter, α); it IS the spot signal, and on a star with ~20% contrast it is large.
# Including it would drag the diameter to a compromise that fits neither — `fit_sphere_ld`'s
# docstring measures exactly that on RW Cep: adding T3φ moved the diameter 2.93 → 2.56 mas and
# took χ²v2/n from 8.9 to 24.8. Here the T3φ residual is the thing the imaging step then
# reconstructs, which is the whole point of fitting the symmetric model first.
#
# POWER LAW (ldtype = 3, Hestroffer): one coefficient, and sufficient for an H-band disc.
#
# MESH RESOLUTION IS NOT FREE. `fit_sphere_ld` was validated against OITOOLS' analytic
# `visibility_ldpow` on RW Cep: at nside_exp 5 the fitted diameter agrees to 0.003% and α to
# 0.07%; at nside_exp 3, to 0.12% and 1.2%. Both are run here so the discretisation error is
# visible rather than assumed.

using ROTIR, Statistics, Printf

const DIR = joinpath(@__DIR__, "zetAnd")
const FILES = sort(filter(f -> endswith(f, ".oifits"), readdir(DIR, join = true)))

R = String[]; say(s) = (push!(R, s); println(s); flush(stdout); nothing)
sayf(fmt, args...) = say(Printf.format(Printf.Format(fmt), args...))

data = [readoifits(f)[1] for f in FILES]
mjd  = [minimum(d.v2_mjd) for d in data]
say("ζ And — six MIRC-X epochs, 2019 Jul 17 … Aug 06")
sayf("%d epochs, %d data points, %.1f d baseline (P_rot = 17.77 d → %.2f rotations)",
     length(data), sum(d -> length(d.v2) + length(d.t3amp) + length(d.t3phi), data),
     maximum(mjd) - minimum(mjd), (maximum(mjd) - minimum(mjd)) / 17.77)

"One fit, reported with its χ² split so the spot residual is visible."
function run_fit(ds, nexp; label)
    tess = tessellation_healpix(nexp)
    r = fit_sphere_ld(ds, tess; tepochs = zeros(length(ds)), ldtype = 3,
                      radius0 = 1.25, ld0 = 0.3,
                      radius_bounds = (0.3, 4.0), ld_bounds = (0.0, 3.0),
                      weights = [1.0, 1.0, 0.0], method = :neldermead,
                      algorithm = :LN_NELDERMEAD, maxeval = 4000, verbose = false)
    npts = sum(d -> length(d.v2) + length(d.t3amp), ds)
    sayf("%-26s nside_exp=%d  D = %.4f mas   α = %.4f   χ²/n(V²+T3amp) = %.3f",
         label, nexp, 2 * r.radius, r.ld1, r.chi2 / max(npts - 2, 1))
    return (D = 2 * r.radius, α = r.ld1)
end

say("")
say("── all six epochs together ──")
a4 = run_fit(data, 4; label = "all epochs")
a5 = run_fit(data, 5; label = "all epochs")
sayf("mesh sensitivity: ΔD = %.4f mas (%.2f%%), Δα = %.4f (%.2f%%)",
     a5.D - a4.D, 100*(a5.D - a4.D)/a5.D, a5.α - a4.α, 100*(a5.α - a4.α)/max(a5.α, 1e-9))

say("")
say("── per epoch, at nside_exp = 5 ──")
say("The scatter across epochs is an EMPIRICAL error bar, and a more honest one than the")
say("formal σ: each epoch sees a different spot configuration, so the spread folds in the")
say("systematic a symmetric model makes on a spotted star.")
per = [run_fit([data[k]], 5; label = "epoch $(k) (MJD $(round(mjd[k], digits=1)))")
       for k in eachindex(data)]
Ds = [p.D for p in per]; αs = [p.α for p in per]
say("")
sayf("diameter : joint %.4f mas | per-epoch %.4f ± %.4f (rms) | spread %.4f–%.4f",
     a5.D, mean(Ds), std(Ds), minimum(Ds), maximum(Ds))
sayf("α        : joint %.4f     | per-epoch %.4f ± %.4f (rms) | spread %.4f–%.4f",
     a5.α, mean(αs), std(αs), minimum(αs), maximum(αs))
say("")
say("Compare against the published 2019A analysis: diameter and limb darkening are the two")
say("numbers that must agree before a joint model+imaging fit on this star means anything.")
println("\n" * "="^78); for l in R; println(l); end; println("="^78)
