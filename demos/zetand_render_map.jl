#!/usr/bin/env julia
# Render the extracted ζ And map with ROTIR's OWN Mollweide plotter, for comparison with the
# figure it was read out of.
#
#     julia --project=. demos/zetand_render_map.jl        # needs PythonPlot
#
# `mollplot_temperature_healpix` gives the same projection, `ang2pix_nest` on a 2000x1000 grid
# so facets come out sharp and boundary ones follow the ellipse exactly, and `gist_heat` as its
# DEFAULT colormap — which is also the ramp measured off the figure's colourbar (0.0125 rms,
# against 0.233 for `afmhot`). Using it rather than a hand-rolled projection removes the
# convention questions at once: pixel ordering, latitude sense and facet registration are then
# the package's rather than a guess.
#
# WHAT DREW THE ORIGINAL IS NOT KNOWN, and it was probably NOT this. Two details point at
# healpy's `mollview` instead: her colourbar carries only a minimum and a maximum, where this
# function writes seven ticks, and her graticule is unlabelled, where this one labels both axes.
# It does not matter for the extraction, which is driven by the fitted plot geometry rather than
# by whatever drew it — the per-pixel agreement in facet interiors (rms 0.0032, correlation
# 0.9971) is the evidence that the geometry is right. It would matter if healpy's default
# `flip = "astro"` were in play, since that mirrors longitude; it is not, because mirroring
# longitude destroys that agreement and costs a factor two in chi2.
#
# The map is stored AS PLOTTED, so nothing is flipped here. The equator flip that ζ And's
# visibilities select belongs to the forward model and is applied in `zetand_published_map.jl`.

using ROTIR, PythonPlot, DelimitedFiles, Printf

function main()
    zet = joinpath(pkgdir(ROTIR), "demos", "zetAnd")
    m = vec(readdlm(joinpath(zet, "zetAnd2019A_map.txt")))
    @printf("%d tessels, %.4f .. %.4f\n", length(m), minimum(m), maximum(m))
    star = create_star(tessellation_healpix(3),
                       default_star_params(0; radius = 1.0, ldtype = 0), 0.0)
    plot_mollweide(m, star; vmin = 0.8025, vmax = 0.9977, colormap = "gist_heat",
                   figtitle = "zet And 2019A — extracted")
    out = joinpath(zet, "zetAnd2019A_map_rotir.png")
    PythonPlot.savefig(out, dpi = 150, bbox_inches = "tight")
    println("wrote ", out)
end

main()
