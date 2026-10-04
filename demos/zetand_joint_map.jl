#!/usr/bin/env julia
# The joint fit's surface beside the published one.
#
#     julia --project=. demos/zetand_joint_v0.jl      # first, writes the maps
#     julia --project=. demos/zetand_joint_map.jl     # then this
#
# Both panels are INTENSITY maps on the same scale, which is the only honest way to compare
# them: ROTIR's imaging path reconstructs surface brightness (`intensity_model = :linear`), and
# `zetand_fig10_match.jl` establishes that SURFING's "temperature" map is an intensity map too
# — it converts linearly, so their kelvin axis is a relabelled brightness. Plotting ours as a
# temperature and theirs as published would be comparing two different quantities.
using ROTIR, PythonPlot, Printf, Statistics, DelimitedFiles

const ZET = joinpath(pkgdir(ROTIR), "demos", "zetAnd")
const PROT = 17.769426

ours  = vec(readdlm(joinpath(ZET, "zetAnd_joint_v0.txt")))
theirs = vec(readdlm(joinpath(ZET, "zetAnd_fig10_epochA_T.txt")))
@printf("ours   %.0f - %.0f, mean %.0f, sd %.1f (sd/mean %.4f)\n",
        minimum(ours), maximum(ours), mean(ours), std(ours), std(ours)/mean(ours))
@printf("theirs %.0f - %.0f, mean %.0f, sd %.1f (sd/mean %.4f)\n",
        minimum(theirs), maximum(theirs), mean(theirs), std(theirs), std(theirs)/mean(theirs))

# a COMMON scale, set by the pair, so the eye compares structure and not two normalisations
lo = min(minimum(ours), minimum(theirs)); hi = max(maximum(ours), maximum(theirs))

# `plot_mollweide` takes a star, because it reads the inclination to grey out the polar cap no
# epoch ever sees. The joint fit's own geometry, so the masked region is the one THIS
# reconstruction could not constrain. (`mollplot_temperature_healpix` is the map-only form but
# lives inside ROTIRPythonPlotExt and is not re-exported.)
fit = (surface_type = 1, radius_x = 1.3406, radius_y = 1.2818, radius_z = 1.2271,
       tpole = 4372.0, ldtype = 3, ld1 = 0.22, ld2 = 0.0, beta = 0.25,
       inclination = 62.59, position_angle = 124.08, rotation_period = PROT)
tess = tessellation_healpix(3; T = Float64)
star = create_star(tess, fit, 0.0)

for (nm, m, ttl) in (("ours", ours, "ROTIR joint fit: shape + map + V0"),
                     ("theirs", theirs, "Roettenbacher et al. 2026, Epoch A"))
    # `colorbar_label` OVERRIDDEN, and the point is not cosmetic. The default is
    # "Temperature (K)", and mislabelling an intensity map as a temperature is the exact error
    # `zetand_fig10_match.jl` spent its length establishing in the published figure. Both of
    # these are surface brightness. `incl` is passed too, so the polar cap no epoch ever sees
    # is greyed out rather than shown as if it were constrained.
    plot_mollweide(m, star; vmin = lo, vmax = hi, colormap = "gist_heat", title = ttl,
                   incl = fit.inclination, colorbar_label = "surface brightness (arb.)")
    PythonPlot.savefig(joinpath(ZET, "zetAnd_joint_map_$nm.png"), dpi = 110,
                       bbox_inches = "tight")
    println("wrote zetAnd_joint_map_$nm.png")
end

# ─── comparing the two maps at all requires putting them in the SAME frame ───────────────
#
# A naive facet-by-facet correlation gives −0.21, which says nothing about the star. Our map
# lives in ROTIR's body frame with t = 0 at the first epoch; theirs is referenced to their
# T0 = MJD 49992.281, and `zetand_surfing_conventions.jl` already measures the longitude zero
# as sitting HALF A PERIOD from that T0. So the frames differ by an unknown rotation, and the
# honest comparison scans it rather than assuming it away.
#
# The scan also covers the longitude PARITY, for the same reason `zetand_fig10_match.jl` had
# to: a reflected map has the same contrast and the same statistics, and only the structure
# tells them apart.
c = tess.unit_xyz[:, 5, :]
lon = [mod(atand(c[i,2], c[i,1]), 360) for i in 1:size(c,1)]
lat = [asind(clamp(c[i,3], -1, 1)) for i in 1:size(c,1)]
"The facet nearest (lon + Δ, ±lat-preserving) — a permutation, so no interpolation smooths it."
function rotperm(Δ, flip)
    out = Vector{Int}(undef, length(lon))
    for i in eachindex(lon)
        L = mod((flip ? -lon[i] : lon[i]) + Δ, 360); B = lat[i]
        out[i] = argmax([sind(B)*sind(lat[j]) + cosd(B)*cosd(lat[j])*cosd(L - lon[j])
                         for j in eachindex(lon)])
    end
    out
end
# a Ref, not a bare binding: a top-level `for` opens a soft scope, so an assignment inside it
# would create a NEW local and leave this one untouched.
const BEST = Ref((-2.0, 0.0, false))
for flip in (false, true), Δ in 0:5:355
    r = cor(ours[rotperm(Float64(Δ), flip)], theirs)
    r > BEST[][1] && (BEST[] = (r, Float64(Δ), flip))
end
best = BEST[]
@printf("\nbest alignment: correlation %.4f at a longitude offset of %.0f deg%s\n",
        best[1], best[2], best[3] ? ", with the longitude REFLECTED" : "")
@printf("  (unaligned, for contrast: %.4f — which is why the raw number means nothing)\n",
        cor(ours, theirs))
@printf("contrast ratio (ours/theirs, as sd/mean): %.3f\n",
        (std(ours)/mean(ours)) / (std(theirs)/mean(theirs)))

println("""

READ THIS BEFORE QUOTING THE CORRELATION. The two reconstructions differ in the shape (ours
fitted jointly, theirs held fixed), the limb darkening (0.22 predicted from a Kurucz
atmosphere against their 0.269 fixed from spectral type), the regularisation (a Sobel prior
against an average of 50 MCMC instances), and V0 (fitted here, fitted there, not the same
value). They also differ in what they CAN see: the polar cap below the inclination is
unobservable and is masked in both panels. A high correlation would be remarkable; a moderate
one is what two honest reconstructions of a spotted star from the same visibilities look like.
The contrast ratio is the number that does not depend on any of that.""")
