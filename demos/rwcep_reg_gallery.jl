#!/usr/bin/env julia
# rwcep_reg_gallery.jl — RW Cep reconstructions side by side: no regularizer, RADFLAT and
# orthoLD, each at several strengths, at two tessellation resolutions.
#
#   julia --project=demos demos/rwcep_reg_gallery.jl
#   NSIDES=3 LD1=1.2 WEIGHTS=1e-1,1e1 julia --project=demos demos/rwcep_reg_gallery.jl
#
# Headless: Agg backend, PNGs into demos/results/, nothing displayed.
#
# ---------------------------------------------------------------------------
# ld1 is held FIXED across every panel, and that is the point
# ---------------------------------------------------------------------------
# The scan in `rwcep_radflat.jl` reconstructs at each of many ld1 and asks where χ² is
# minimised. Plotting the best-fit map of each regularizer would therefore compare images
# made at DIFFERENT ld1, confounding "what the regularizer does to the map" with "where the
# regularizer put the limb darkening". Here ld1 is pinned at the parametric value (α = 1.58
# from the V²+T3amp fit), so every panel differs only by its regularizer, and the visual
# comparison means what it looks like it means.
#
# The colour scale is likewise SHARED across all panels of a figure: per-panel
# normalisation would hide exactly the contrast changes being looked for.

ENV["MPLBACKEND"] = get(ENV, "MPLBACKEND", "Agg")
using ROTIR, PythonPlot, Printf, Statistics

const HERE   = @__DIR__
const OUT    = joinpath(HERE, "results"); mkpath(OUT)
const OIFITS = joinpath(HERE, "MIRCX_L3.022Dec23.RW_Cep.MIRCX_IDL.deepclean.AVG15m.oifits")
const NSIDES = parse.(Int, split(get(ENV, "NSIDES", "3,4"), ','))
const WEIGHTS = parse.(Float64, split(get(ENV, "WEIGHTS", "1e-2,1e0,1e2"), ','))
const MAXIT  = parse(Int, get(ENV, "MAXITER", "20000"))
const NBINS  = parse(Int, get(ENV, "NBINS", "6"))
# TV2 ON. Setting it to zero "isolates the radial priors" in the sense that nothing else
# constrains the map — but the result is not a reconstruction anyone would use: with no
# smoothness prior at all, single tessels reach 15-48x the map mean and the images are
# pixel-noise. Those maps are kept out of the default for that reason. TVW=0 still works if
# you want to see the pathology, but do not read structure off it.
# Smoothing prior. "sobel2" penalises |∇x|² on the tangent plane (k² response, normalised by
# mean(x)² and by solid angle, so one weight transfers across BOTH nside values swept here);
# "tv2" penalises curvature (k⁴) and does neither, so its weight is only valid at one nside.
const SMOOTH = get(ENV, "SMOOTH", "sobel2")
SMOOTH in ("tv2", "sobel", "sobel2") || error("SMOOTH must be tv2, sobel or sobel2")
# Default measured by the weight ladder in demos/rwcep_radflat.jl (ALPHAS=0, TVWS=0,1e1..1e4
# on RW Cep): chi2_v2 runs 2.00 / 1.93 / 2.36 / 3.68 / 5.96, so ~1e1 is the only rung that
# improves on no smoothing at all. Above ~1e3 a sharp chi2(ld1) minimum appears, but it is
# the PRIOR choosing ld1, not the data — chi2 is 3x worse there.
#
# It does NOT transfer blindly to another target: "sobel2" is scale- and resolution-
# normalised, so tpole and nside drop out, but the balance against chi2 still scales with the
# NUMBER of data points. Re-run the ladder per dataset.
const TVW    = parse(Float64, get(ENV, "TVW", SMOOTH == "tv2" ? "1e-4" : "1e1"))
const D_FIT  = parse(Float64, get(ENV, "DFIT", "2.930"))
const LD1    = parse(Float64, get(ENV, "LD1", "1.58"))

data = readoifits_multiepochs([OIFITS], warn = false, verbose = false, T = Float64)[1, :]
d    = data[1]
@printf("RW Cep: %d V², %d T3amp, %d T3φ — D = %.3f mas, ld1 fixed at %.2f, TV2 = %.0e\n",
        d.nv2, d.nt3amp, d.nt3phi, D_FIT, LD1, TVW)

params() = (surface_type = 0, radius = D_FIT/2, tpole = 3900.0, ldtype = 3,
            ld1 = LD1, ld2 = 0.0, inclination = 90.0, position_angle = 0.0,
            rotation_period = 1.0)

"Reconstruct with RADFLAT weight `rf` and orthoLD weight `ol`; returns map, star, metrics."
function reconstruct(nside, rf, ol)
    p     = params()
    stars = create_star_multiepochs(tessellation_healpix(nside), p, [0.0])
    setup_oi!(data, stars)
    x0    = parametric_temperature_map(p, stars[1])
    regs  = Any[]
    TVW > 0 && push!(regs, [SMOOTH, TVW, SMOOTH == "tv2" ? tv_neighbors_healpix(nside) :
                                       sobel_gradient_healpix(nside), 1:length(x0)])
    bins  = radflat_bins(stars[1]; nbins = NBINS)
    rf > 0 && push!(regs, ["radflat", rf, bins, bins.idx])
    if ol > 0
        od = orthold_direction(stars[1], x0, p)
        push!(regs, ["orthold", ol, od, od.idx])
    end
    local x
    for it in (MAXIT, MAXIT ÷ 4, MAXIT ÷ 16, 500, 100)
        try
            x = image_reconstruct_oi(x0, data, stars; maxiter = it,
                                     regularizers = regs, verbose = false); break
        catch err
            err isa AssertionError || rethrow(); it == 100 && rethrow()
        end
    end
    xb = x[bins.idx]; m = mean(xb)
    prof = [mean(xb[bins.bin .== k]) for k in 1:NBINS] ./ m
    within = mean(mean(abs2, xb[bins.bin .== k] .- mean(xb[bins.bin .== k]))
                  for k in 1:NBINS if count(==(k), bins.bin) > 1)
    v2m, t3am, t3pm = cvis_to_obs(ComplexF64.(poly_to_cvis(x, stars[1])), d)
    return x, stars[1],
           (v2   = sum(abs2, (v2m .- d.v2)./d.v2_err)/d.nv2,
            t3a  = sum(abs2, (t3am .- d.t3amp)./d.t3amp_err)/d.nt3amp,
            t3p  = sum(abs2, mod360(t3pm .- d.t3phi)./d.t3phi_err)/d.nt3phi,
            prof = sqrt(mean(abs2, prof .- 1)), az = sqrt(within)/m)
end

# panel list: unregularized control, then each prior at each strength
panels = vcat([("none", 0.0, 0.0)],
              [(@sprintf("RADFLAT %.0e", w), w, 0.0) for w in WEIGHTS],
              [(@sprintf("orthoLD %.0e", w), 0.0, w) for w in WEIGHTS])

for nside in NSIDES
    @printf("\nnside %d (%d tessels)\n", nside, 12*4^nside)
    res = [(nm, reconstruct(nside, rf, ol)...) for (nm, rf, ol) in panels]
    # One shared colour scale over every panel — but on maps NORMALISED TO UNIT MEAN.
    # χ² is invariant under x -> c·x (the visibilities are flux-normalised), so the map's
    # overall level is an unconstrained direction and different regularizers converge to
    # wildly different absolute levels. Sharing an ABSOLUTE scale across panels therefore
    # compresses every panel but the brightest into a single flat colour, which is what a
    # first attempt at this figure did. Normalising per panel compares STRUCTURE, which is
    # the only thing the data constrain and the only thing being compared here.
    maps = [r[2][r[3].index_quads_visible] ./ mean(r[2][r[3].index_quads_visible]) for r in res]
    # ROBUST limits, not min/max. A reconstruction with a weak (or absent) smoothness prior
    # puts a few tessels tens of times above the mean; scaling to the extremes then buries
    # every panel in the bottom of the colormap. The 2nd-98th percentile over the pooled
    # panels keeps the comparison legible while the per-panel range printed on each axes
    # says exactly how much is being clipped.
    pooled = sort(vcat(maps...))
    q(f) = pooled[clamp(round(Int, f*length(pooled)), 1, length(pooled))]
    lo, hi = q(0.02), q(0.98)
    # `rgba`, `get_cmap` and `_padded_norm` are internal to ROTIR (not exported), so
    # qualify them. `_padded_norm` returns a matplotlib Normalize, broadcast below.
    cnorm = ROTIR._padded_norm(lo, hi, hi - lo; cfloor = 0.08)
    cmap  = ROTIR.get_cmap("gist_heat")

    ncol = 4
    nrow = ceil(Int, length(res)/ncol)
    # `squeeze = false` is load-bearing: matplotlib drops a length-1 dimension, so a single
    # row of panels comes back as a 1-D array and the 2-D indexing below raises. With it,
    # `axs` is always 2-D. Indices are 0-based — it is a numpy array, not a Julia Matrix.
    fig, axs = pyplot.subplots(nrow, ncol, figsize = (4.1*ncol, 4.4*nrow), squeeze = false)
    axf = [axs[i-1, j-1] for i in 1:nrow for j in 1:ncol]
    for (k, (nm, x, star, m)) in enumerate(res)
        ax = axf[k]
        add_tessel_collection!(ax, star, ROTIR.rgba(cmap, cnorm.(maps[k])))
        draw_limb!(ax, star; color = "black", linewidth = 0.8)
        r = 1.15 * maximum(abs, star.proj_west)
        ax.set_xlim(r, -r); ax.set_ylim(-r, r); ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(nm, fontsize = 11)
        ax.text(0.02, 0.02,
                @sprintf("χ²ᵥ₂ %.2f  χ²ₜ₃ₚ %.1f\nprof %.4f  az %.3f\nmap %.2f–%.2f × mean",
                         m.v2, m.t3p, m.prof, m.az, minimum(maps[k]), maximum(maps[k])),
                transform = ax.transAxes, fontsize = 7.5, va = "bottom",
                # OPAQUE. At alpha < 1 the tessels show through this box as a washed-out
                # rectangle over the bottom-left of the disk, which reads exactly like a
                # patch of misbehaving tessels.
                bbox = Dict("fc" => "white", "alpha" => 1.0, "ec" => "0.7", "pad" => 2.0))
        @printf("  %-14s χ²ᵥ₂ %.2f  χ²ₜ₃ₐ %.2f  χ²ₜ₃ₚ %6.1f  profrms %.4f  azrms %.4f\n",
                nm, m.v2, m.t3a, m.t3p, m.prof, m.az)
    end
    for k in (length(res)+1):length(axf); axf[k].axis("off"); end
    fig.suptitle(@sprintf("RW Cep — nside %d (%d tessels), ld1 = %.2f fixed, %s = %.0e",
                          nside, 12*4^nside, LD1, SMOOTH, TVW), fontsize = 13)
    fig.tight_layout()
    # The smoothing prior AND its weight go in the FILENAME. They change the answer
    # qualitatively (orthoLD is free at weight 0 and costs 2.7x in chi2 at the tv2 = 1e-4
    # working point), so two runs differing only in SMOOTH or TVW must not land on the same
    # file — an earlier version of this script silently overwrote one with the other.
    f = joinpath(OUT, @sprintf("rwcep_gallery_nside%d_%s%s.png", nside, SMOOTH,
                               TVW == 0 ? "0" : @sprintf("%.0e", TVW)))
    fig.savefig(f, dpi = 130, bbox_inches = "tight")
    pyplot.close(fig)
    @info "wrote $f"
end
