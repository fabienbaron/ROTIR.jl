# Spica (alpha Vir) — global search of the orbit, interferometry + radial velocities
#
#   julia -t 8 --project=demos demos/spica_global_orbit_fit.jl
#
#   NLIVE=2000 ...   # live points (default 1000)
#   NEFF=20000 ...   # effective posterior sample size (default 10000)
#
# The first of two steps. This one samples the whole prior box with Nautilus and writes
# results/spica_orbitfit_joint_rv_best_nautilus.txt; demos/spica_refine_orbit_fit.jl then
# refines that solution with the tessellated Roche model by Levenberg-Marquardt.
#
# ---------------------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------------------
# Interferometry: two power-law limb-darkened discs (fit_orbit, model = :analytic), the
# secondary placed at each uv point's own observation time. Data: 2007_2012_2015.Spica.oifits,
# split into 11 epochs at gaps longer than half a day.
#
# Radial velocities: data/spica_rv.csv with sigma = 7 + 30 exp(-|V2 - V1|/100) km/s
# (spica_rv() in spica_params.jl), the semi-amplitudes following from the parallax and the
# mass ratio (rv_parametrization = :physical).
#
# Fixed:
#   P       4.01454 d, ANOMALISTIC (Shobbrook et al. 1972; Aufdenberg et al. 2015).
#   LD      c1_ld1, c2_ld1: H-band power-law exponents fitted to Kurucz ATLAS9 intensities at
#           each star's polar Teff and logg (spica_kurucz_ld() in spica_params.jl), the values
#           demos/spica_refine_orbit_fit.jl uses. The data do not constrain them: free, the
#           primary's runs to the edge of its box, and the secondary, about one resolution
#           element across at the longest baseline and 15 % of the flux, has no measurable
#           limb darkening.
#
# Free, with boxes set by physics and not by earlier fits:
#   a            0.8-2.5 mas
#   i            0-180 deg (both senses of motion)
#   Omega        0-360 deg, periodic (the velocities break the (Omega, omega) mirror)
#   omega        0-360 deg, periodic
#   e            0-0.4
#   T0           one period
#   domega       0-0.02 deg/day (apsidal period from 49 yr to infinity; apsidal motion from
#                tides, rotation and general relativity is prograde)
#   c1_diameter  0.2-1.6 mas, c2_diameter 0.1-1.0 mas
#   f            0.02-0.8 (flux ratio secondary/primary; below 1, so the brighter star is
#                component 1)
#   plx          5-25 mas, q 0.2-1.5, gamma -50 to 50 km/s
#
# The written values are posterior medians, with the 16 and 84 % quantiles. Their widths take
# the OIFITS errors at face value, and the disc model leaves OI chi2r near 2.

using ROTIR, Nautilus, Printf
using OITOOLS: split_by_time

const HERE = @__DIR__
include(joinpath(HERE, "spica_params.jl"))

const NLIVE = parse(Int, get(ENV, "NLIVE", "1000"))
const NEFF  = parse(Int, get(ENV, "NEFF", "10000"))

data, _ = split_by_time(readoifits(joinpath(HERE, "data", "2007_2012_2015.Spica.oifits"))[1, 1]; gap = 0.5)
rv = spica_rv()

# Starting values only set the RV block's defaults and the centre of the T0 box; the sampler
# draws from the boxes.
elements = (a = 1.54, i = 116.0, Omega = 310.0, omega = 255.0, e = 0.12, P = P_ORB,
            T0 = T0_ORB, dP = 0.0, domega = DOMEGA)
LD1, LD2 = spica_kurucz_ld()
@printf("Kurucz H-band power-law limb darkening: primary %.4f, secondary %.4f\n", LD1, LD2)
star1 = LimbDarkenedDisk(diameter = 2RPOLE1, law = :power, ld1 = LD1)
star2 = LimbDarkenedDisk(diameter = 2RPOLE2, law = :power, ld1 = LD2)

free = [:a, :i, :Omega, :omega, :e, :T0, :domega, :c1_diameter, :c2_diameter, :f,
        :plx, :q, :gamma]
bounds = Dict(:a => (0.8, 2.5), :i => (0.0, 180.0), :Omega => (0.0, 360.0),
              :omega => (0.0, 360.0), :e => (0.0, 0.4), :domega => (0.0, 0.02),
              :c1_diameter => (0.2, 1.6), :c2_diameter => (0.1, 1.0),
              :f => (0.02, 0.8), :plx => (5.0, 25.0), :q => (0.2, 1.5), :gamma => (-50.0, 50.0))

res = fit_orbit(data, star1, star2; elements = elements, flux_ratio = 0.2, free = free,
                bounds = bounds, rv = rv, method = :nautilus, n_live = NLIVE, n_eff = NEFF,
                verbose = true)

# --- report ---------------------------------------------------------------------------------
idx = res.spec.free
cs = res.chi2_split
nv2, nt3a, nt3p = res.data.nv2, res.data.nt3amp, res.data.nt3phi
noi = nv2 + nt3a + nt3p
@printf("\nOI chi2r %.3f over %d points (V2 %.2f, T3amp %.2f, T3phi %.2f); RV chi2 %.1f over %d\n",
        sum(cs[1:3]) / noi, noi, cs[1] / nv2, cs[2] / nt3a, cs[3] / nt3p, cs[4], length(rv.jd))
@printf("logZ = %.3f +- %.3f\n", res.logz, res.logzerr)
println("\nparameter        median          16%             84%")
for (j, k) in enumerate(idx)
    @printf("  %-12s %14.8g  %14.8g  %14.8g\n", res.names[k], res.params[k], res.q16[j], res.q84[j])
end
for (k, v) in pairs(res.rvpar)
    k in (:plx, :q, :gamma, :jit_a, :jit_b) || @printf("  %-12s %14.8g\n", k, v)
end

out = joinpath(HERE, "results", "spica_orbitfit_joint_rv_best_nautilus.txt")
open(out, "w") do io
    println(io, "# demos/spica_global_orbit_fit.jl: Nautilus over the full prior box; power-law discs on")
    println(io, "# 2007_2012_2015.Spica.oifits + spica_rv.csv (sigma = 7 + 30 exp(-|V2-V1|/100) km/s).")
    @printf(io, "# OI chi2r %.3f over %d points, RV chi2 %.1f over %d; logZ %.3f. P fixed at %.5f d.\n",
            sum(cs[1:3]) / noi, noi, cs[4], length(rv.jd), res.logz, P_ORB)
    println(io, "# name\tmedian\tq16\tq84 (q16 = q84 = median when fixed or derived)")
    for (k, n) in enumerate(res.names)
        j = findfirst(==(k), idx)
        lo, hi = j === nothing ? (res.params[k], res.params[k]) : (res.q16[j], res.q84[j])
        println(io, n, '\t', res.params[k], '\t', lo, '\t', hi)
    end
    for (k, v) in pairs(res.rvpar)
        k in res.names || println(io, k, '\t', v, '\t', v, '\t', v)
    end
end
println("saved ", out)
post = joinpath(HERE, "results", "spica_orbitfit_joint_rv_posterior_nautilus.txt")
open(post, "w") do io
    println(io, join(string.(res.names[idx]), '\t'))
    for r in eachrow(res.posterior); println(io, join(r, '\t')); end
end
println("saved ", post)
