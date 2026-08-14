# beta Lyrae — orbit fit through ROTIR's generic path, with the disc PA TIED to the node
#
#   julia --project=demos demos/betlyr/betlyr_orbit_fit_tied.jl
#
#   NLIVE=800 julia --project=demos demos/betlyr/betlyr_orbit_fit_tied.jl   # tighter run
#   TIED=0    julia --project=demos demos/betlyr/betlyr_orbit_fit_tied.jl   # free pa, for comparison
#   METHOD=neldermead ...                                                   # quick smoke test
#
# ---------------------------------------------------------------------------------------
# What is different from betlyr_orbit_fit_ultranest.jl
# ---------------------------------------------------------------------------------------
# That driver fits the bespoke model in `betlyr_model.jl`, in which all ten parameters are
# independent. Here the same system is expressed with `orbit_fit_spec`, which allows a
# parameter to be TIED to another rather than fitted:
#
#     ties = Dict(:c2_pa => "-Omega")
#
# A disc in the orbital plane has its projected major axis along the line of nodes, so its
# position angle is not an independent quantity. Fitting it freely both spends a parameter
# and permits a disc oriented inconsistently with the orbit fitted alongside it. On the
# 2006-2007 data the freely-fitted angle lands 1.9 deg from the tied prediction, so this
# costs almost nothing in fit quality and removes a degree of freedom.
#
# SIGN. `EllipticalGaussian` takes `pa` as the angle whose axis is scaled by `ratio < 1` —
# the MINOR axis — measured from +RA toward +Dec, while Omega is the node measured North
# through East. Major axis on the node line is therefore `-Omega`, not `Omega`. Check this
# before copying the tie to another system: a wrong tie cannot be spotted from the chi2.
#
# NOT A REPLACEMENT for the bespoke driver. This model has no occultation term, and the
# components overlap on the sky at 3 of 7 epochs by half-FWHM (6 of 7 using the Gaussian's
# 2-sigma extent), with the disc in front in October and behind in July. Whatever this fit
# reports, that systematic is still present — see the caveats in docs/src/guides/orbits.md.

using ROTIR, Printf, Statistics, DelimitedFiles

const HERE   = @__DIR__
const OLDDIR = joinpath(HERE, "older_data")     # Zhao et al.'s 2006-2007 MIRC nights
const TIED   = get(ENV, "TIED", "1") == "1"
const METHOD = Symbol(get(ENV, "METHOD", "ultranest"))
const NLIVE  = parse(Int, get(ENV, "NLIVE", "400"))
const OUT    = joinpath(HERE, "results")

files = sort!(filter(f -> endswith(f, ".oifits"), readdir(OLDDIR, join = true)))
isempty(files) && error("no OIFITS under $OLDDIR")
data = [readoifits(f; T = Float64)[1, 1] for f in files]
@printf("%d epochs, %d data points\n", length(data),
        sum(d.nv2 + d.nt3amp + d.nt3phi for d in data))

# Starting point: the converged solution from the bespoke UltraNest run, so this tests the
# tie rather than the ability to find the basin from scratch.
elements = (a = 0.8958, i = 91.244, Omega = 73.529, omega = 90.0, e = 0.0,
            P = 12.9616, T0 = 2454283.0555, dP = 0.0)
donor = UniformDisk(diameter = 0.590)                        # B6-8 II
disc  = EllipticalGaussian(fwhm = 0.5825, ratio = 0.6423, pa = 108.35)

free = Symbol[:a, :i, :Omega, :T0, :f, :c1_diameter, :c2_fwhm, :c2_ratio]
ties = TIED ? Dict(:c2_pa => "-Omega") : Dict{Symbol,String}()
TIED || push!(free, :c2_pa)

@printf("\n%s: %d free parameters%s\n", TIED ? "TIED" : "FREE",
        length(free), TIED ? " + c2_pa tied to -Omega" : "")

res = fit_orbit(data, donor, disc; elements = elements, flux_ratio = 0.8179,
                free = free, ties = ties, method = METHOD,
                min_num_live_points = NLIVE, verbose = true,
                # P is pinned by the eclipse ephemeris, not by 269 d of visibilities; a is
                # held near the published value so the tie is what is being tested.
                bounds = Dict(:a => (0.80, 1.00), :i => (85.0, 95.0),
                              :Omega => (60.0, 90.0), :c2_ratio => (0.2, 1.0)))

@printf("\n%-16s %14s\n", "parameter", "fitted")
for (n, v) in zip(res.names, res.params)
    tag = n in keys(ties) ? "  (tied)" : (n in free ? "" : "  (fixed)")
    @printf("%-16s %14.6f%s\n", n, v, tag)
end
@printf("\nchi2/n = %.4f   (V2 %.2f, T3amp %.2f, T3phi %.2f)\n",
        res.chi2_red, res.chi2_split[1], res.chi2_split[2], res.chi2_split[3])

# The tie is only meaningful if the convention is right: print the major axis in the same
# frame as Omega so the two can be compared directly.
jpa = findfirst(==(:c2_pa), res.names)
@printf("disc major axis astro PA = %.2f deg   line of nodes Omega = %.2f deg\n",
        mod(-res.params[jpa], 180), mod(res.params[3], 180))

mkpath(OUT)
tag = TIED ? "tied" : "free"
writedlm(joinpath(OUT, "orbitfit_$(tag)_best.txt"), hcat(string.(res.names), res.params))
haskey(res, :posterior) &&
    writedlm(joinpath(OUT, "orbitfit_$(tag)_posterior.txt"), res.posterior)
@info "wrote results to $OUT"
