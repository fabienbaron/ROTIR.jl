#!/usr/bin/env julia
# Render the image-only reconstruction beside the published map it is being compared against.
#
#     julia --project=<env with PythonPlot> demos/zetand_render_image_only.jl
#
# Both panels on the SAME colour scale (0.80-1.09, the union of the two ranges) and in HER
# plotting frame: the reconstruction is produced in ROTIR's model frame, so it is mapped back
# through the inverse of `healpix_flip_3.txt` before plotting. Comparing two maps on two
# different scales, or in two different frames, is how a disagreement gets hidden.

using ROTIR, PythonPlot, DelimitedFiles, Printf, Statistics
function main()
    zet = joinpath(pkgdir(ROTIR), "demos", "zetAnd")
    pub = vec(readdlm(joinpath(zet, "zetAnd2019A_map.txt")))
    flip = vec(Int.(readdlm(joinpath(zet, "healpix_flip_3.txt"))))
    img = vec(readdlm(joinpath(zet, "zetAnd2019A_image_only.txt")))     # model frame, mean 1
    star = create_star(tessellation_healpix(3),
                       default_star_params(0; radius=1.0, ldtype=0), 0.0)
    # put ours back into HER plotting frame so the two panels are comparable
    inv = zeros(Int, length(flip)); for (i, j) in enumerate(flip); inv[j] = i; end
    ours = img[inv]
    ours ./= mean(ours); p = pub ./ mean(pub)
    @printf("published %.3f-%.3f   ours %.3f-%.3f\n",
            minimum(p), maximum(p), minimum(ours), maximum(ours))
    opts = (colormap="gist_heat", labels=false, lon_color="black", lat_color="black",
            grid_style=":", grid_alpha=0.9, lon_spacing=30, lat_spacing=30,
            axis_lines=true, colorbar_ticks=2, colorbar_label="")
    lo, hi = 0.80, 1.09
    plot_mollweide(ours, star; vmin=lo, vmax=hi, title="ROTIR image-only, sobel 3000", opts...)
    PythonPlot.savefig(joinpath(zet, "zetAnd2019A_image_only.png"), dpi=150, bbox_inches="tight")
    plot_mollweide(p, star; vmin=lo, vmax=hi, title="Roettenbacher et al., extracted", opts...)
    PythonPlot.savefig(joinpath(zet, "zetAnd2019A_published_rescaled.png"), dpi=150, bbox_inches="tight")
    println("wrote both")
end
main()
