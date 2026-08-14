# RW Cep — snapshot imaging on a spheroid, with and without the RADFLAT regularizer
#
#   julia --project=demos demos/rwcep_radflat.jl
#   NSIDE=4 ALPHAS=0,1e-3,1e-2,1e-1 julia --project=demos demos/rwcep_radflat.jl
#
# ---------------------------------------------------------------------------
# The workflow: fit the size FIRST, then image
# ---------------------------------------------------------------------------
# A reconstruction on a spheroid needs the star's angular size up front — the tessellation is
# built at a fixed `radius`, and getting it wrong is not a small error. A 55 % radius error on
# this dataset gave χ²/n ≈ 4000. So step 1 is always a parametric fit.
#
# Fit the size and limb darkening on **V² and T3amp only**. Closure phases are excluded on
# purpose: RW Cep's T3φ has an rms of 118°, while any centrosymmetric model can produce only
# 0° or 180°. Including them drags a symmetric model toward a compromise that is worse in the
# amplitudes (χ²v2/n 8.9 → 24.8, χ²t3a/n 5.4 → 9.5) and shrinks the diameter by 13 % — it is
# not fitting the star better, it is contorting to chase a signal the model cannot represent.
# The large residual T3φ is not a failure: it IS the asymmetry, and it is what the imaging
# step reconstructs.
#
# ---------------------------------------------------------------------------
# What RADFLAT is for
# ---------------------------------------------------------------------------
# RADFLAT (J. Monnier, priv. comm. 2026-07-01) breaks a specific degeneracy. On a NON-ROTATING
# star from a single epoch, the surface regularizer is weak enough that the fit can park dark
# or bright spots on the limb, and those trade off almost exactly against the limb-darkening
# coefficient: χ² barely moves while the LD coefficient wanders. Monnier saw exactly this on
# RW Cep with SURFING — α = 1.07 ± 0.19 — and forcing the azimuthally averaged radial profile
# flat gave α = 0.49 ± 0.02, with χ² only slightly worse and the image essentially unchanged.
#
# He got his spread from the scatter over many reconstructions. With one epoch there is no
# such scatter, so this script measures the degeneracy directly instead:
#
#     scan the LD coefficient, reconstruct at each value, and look at χ²(ld1).
#
# A FLAT χ²(ld1) *is* the degeneracy — every LD value fits equally well because the surface
# map absorbs the difference. If RADFLAT works, χ²(ld1) should develop a real minimum, and
# the curvature at that minimum is what turns into an error bar.
#
# NOT a reproduction of Monnier's numbers: SURFING has its own parameterisation, priors and
# error model. The question is whether the degeneracy lifts, not whether 0.49 comes back.
# RADFLAT is deliberately wrong for rotators — several epochs already break this degeneracy,
# and RADFLAT would then suppress real structure.
#
# ---------------------------------------------------------------------------
# ⚠ CONVERGENCE IS PART OF THE MEASUREMENT
# ---------------------------------------------------------------------------
# Every number here is a DIFFERENCE of χ² between reconstructions, so each one has to be at
# its own minimum or the difference is measuring the optimizer instead. At nside 3 the drift
# in χ²/n from 300 to 3000 iterations is 0.14–0.66 depending on ld1 — LARGER than the whole
# spread across the LD grid, which is the quantity being measured. `image_reconstruct_oi`
# uses `gtol = (0, 1e-8)`, tight enough that runs stop on `maxiter` rather than on the
# gradient, so MAXITER is the convergence control here and its default is set from the
# plateau of a convergence ladder, not from what is quick. Raising it is cheap: at nside 3
# a reconstruction is ~1.2 s per 300 iterations.

using ROTIR, Printf, Statistics

const HERE   = @__DIR__
const OIFITS = joinpath(HERE, "MIRCX_L3.022Dec23.RW_Cep.MIRCX_IDL.deepclean.AVG15m.oifits")
const NSIDE  = parse(Int, get(ENV, "NSIDE", "3"))          # 3 matches Monnier's SURFING run
const NBINS  = parse(Int, get(ENV, "NBINS", "6"))          # Monnier uses 6
const MAXIT  = parse(Int, get(ENV, "MAXITER", "20000"))    # convergence control — see above
# Which smoothing prior backs the radial regularizers. "tv2" is ‖Lx‖² (graph Laplacian, so it
# penalises CURVATURE, k⁴ response, and is neither scale- nor resolution-normalised);
# "sobel"/"sobel2" are built on the tangent-plane gradient (k², normalised by mean(x)² and by
# solid angle). Declared before TVWS because the default ladder depends on which family it is.
const SMOOTH = get(ENV, "SMOOTH", "sobel2")
SMOOTH in ("tv2", "sobel", "sobel2") || error("SMOOTH must be tv2, sobel or sobel2")
# Smoothing weights to sweep, INCLUDING ZERO. The smoothing prior pulls the map toward a
# constant on its own, so leaving it on confounds exactly what this script is trying to
# isolate: a χ²(ld1) minimum could be the smoothing flattening the map rather than
# RADFLAT/RADIALVAR doing it. With 2819 data against 428 patches the problem is
# over-determined even with no regularizer at all, so 0 is a legitimate run and is the
# cleanest isolation of the pair.
#
# The default ladder is for SMOOTH = "sobel"/"sobel2", whose weights are O(1–10⁴): they are
# normalised by mean(x)² and by solid angle, so they sit on a completely different scale from
# "tv2" (O(1e-7–1e-4), which rides on the map's absolute temperature). Pass TVWS explicitly
# when switching families — a weight carried across is meaningless.
const TVWS   = parse.(Float64, split(get(ENV, "TVWS",
                          SMOOTH == "tv2" ? "0,1e-7,1e-4" : "0,1e1,1e2,1e3,1e4"), ','))
const ALPHAS = parse.(Float64, split(get(ENV, "ALPHAS", "0,1e-3,1e-2,1e-1,1e0,1e2"), ','))
# Shared weight for RADFLAT/RADIALVAR in the four-way comparison. 100 is the value tuned on
# α Cen A (see the REGW ladder in demos/alphacen_ld.jl): below it the residual ~1 % of map
# structure biases the recovered ld1 high, above ~1e3 the stiff quadratic penalty defeats
# VMLMB. Carried over unchanged so the two benchmarks are directly comparable.
const REGW   = parse(Float64, get(ENV, "REGW", "100"))
const DOPLOT = get(ENV, "PLOT", "1") == "1"
isfile(OIFITS) || error("RW Cep OIFITS not found at $OIFITS")

# --- 1. Data -----------------------------------------------------------------
# readoifits returns Array{OIdata}(nwavbin, ntimebin); defaults give 1x1, hence [1,1].
data_all = readoifits_multiepochs([OIFITS], warn = false, verbose = false, T = Float64)
data     = data_all[1, :]                       # single epoch
d        = data[1]
nd       = d.nv2 + d.nt3amp + d.nt3phi
@printf("RW Cep: %d V², %d T3amp, %d T3φ (n=%d), MJD span %.1f min — single snapshot\n",
        d.nv2, d.nt3amp, d.nt3phi, nd, 1440*(maximum(d.v2_mjd) - minimum(d.v2_mjd)))

# --- 2. Parametric fit: angular size + limb darkening, on V² + T3amp ----------
# Hestroffer power law I(μ) ∝ μ^α, which fits far better here than a linear law — the linear
# coefficient pins at its physical ceiling u = 1 (u > 1 means negative intensity at the limb).
function chi2_v2t3a(p)
    v2m, t3am, _ = cvis_to_obs(ComplexF64.(visibility_ldpow(p, d.uv)), d)
    sum(abs2, (v2m .- d.v2)./d.v2_err) + sum(abs2, (t3am .- d.t3amp)./d.t3amp_err)
end
function fit_size()                                     # coarse grid, then refine
    bD, bA, bc = 0.0, 0.0, Inf
    for D in 1.0:0.05:6.0, α in 0.0:0.10:3.0
        c = chi2_v2t3a([D, α]); c < bc && ((bc, bD, bA) = (c, D, α))
    end
    for D in bD-0.1:0.01:bD+0.1, α in max(0,bA-0.2):0.02:bA+0.2
        c = chi2_v2t3a([D, α]); c < bc && ((bc, bD, bA) = (c, D, α))
    end
    return bD, bA, bc/(d.nv2 + d.nt3amp)
end
D_fit, α_fit, c_fit = fit_size()
@printf("parametric fit (V²+T3amp): D = %.3f mas, α = %.2f, χ²/n = %.2f  →  polar radius %.3f mas\n",
        D_fit, α_fit, c_fit, D_fit/2)
_, _, t3p = cvis_to_obs(ComplexF64.(visibility_ldpow([D_fit, α_fit], d.uv)), d)
@printf("residual χ²t3φ/n = %.0f  (T3φ rms in data %.0f° vs 0/180° for any symmetric model)\n\n",
        sum(abs2, mod360(t3p .- d.t3phi)./d.t3phi_err)/d.nt3phi, sqrt(mean(abs2, mod360(d.t3phi))))

# --- 3. Reconstruction, LD scanned about the fitted value --------------------
tessels = tessellation_healpix(NSIDE)
tvinfo  = tv_neighbors_healpix(NSIDE)
sobelinfo = SMOOTH == "tv2" ? nothing : sobel_gradient_healpix(NSIDE)
smoothinfo = SMOOTH == "tv2" ? tvinfo : sobelinfo
# Wide enough that the minimum can be INTERIOR. A minimum pinned at a grid edge is not a
# measurement of anything — it only says "further that way was never tried".
LDGRID  = parse.(Float64, split(get(ENV, "LDGRID",
              join([@sprintf("%.2f", v) for v in 0.4:0.2:3.0], ',')), ','))

params(ld1) = (surface_type = 0, radius = D_fit/2, tpole = 3900.0,   # RW Cep ≈ M0 Iab
               ldtype = 3, ld1 = ld1, ld2 = 0.0,
               inclination = 90.0, position_angle = 0.0, rotation_period = 1.0)

"""
Regularizer sets. `rf` weights RADFLAT (between-annuli structure), `rv` weights RADIALVAR
(within-annuli). They are the two halves of a variance decomposition of the visible disk.

RADIALVAR IS OFF BY DEFAULT ON THIS STAR, and the reason is measured rather than assumed.
Set `REGSETS=all` to reproduce it: switching RADIALVAR on — alone or inside "both" — sends
χ²ₜ₃ₚ/n from 1.5 to 8.0, while RADFLAT alone leaves it at 1.6. Closure phases are the
signature of genuine asymmetry, so that is RADIALVAR removing real surface structure, not
artefacts. It is the expected outcome once stated plainly: stellar surface features are
AZIMUTHAL (spots at assorted longitudes), which is exactly the quantity RADIALVAR penalises,
whereas the radial profile RADFLAT targets is dominated by limb darkening — a model
degeneracy rather than astrophysics.

That is also why `demos/alphacen_ld.jl` reaches the opposite conclusion without
contradiction: α Cen has no resolved surface structure for RADIALVAR to destroy, so there
"both" is safe and is the only combination that pins ld1. Here it would cost the science.
"""
const REGSETS = get(ENV, "REGSETS", "default") == "all" ?
                [("none",      0.0,  0.0,  0.0),
                 ("radflat",   REGW, 0.0,  0.0),
                 ("radialvar", 0.0,  REGW, 0.0),
                 ("both",      REGW, REGW, 0.0),
                 ("orthold",   0.0,  0.0,  REGW)] :
                # Default: RADFLAT and orthoLD swept over their weight, head to head, with
                # α = 0 as the shared unregularized control. Entries are (name, rf, rv, ol).
                vcat([("none", 0.0, 0.0, 0.0)],
                     [(@sprintf("radflat %.0e", α), α, 0.0, 0.0) for α in ALPHAS if α > 0],
                     [(@sprintf("orthold %.0e", α), 0.0, 0.0, α) for α in ALPHAS if α > 0])

"""
Reconstruct at fixed LD coefficient `ld1`, RADFLAT weight `rf`, RADIALVAR weight `rv` and
smoothing weight `tvw` (applied with whichever prior `SMOOTH` selects).

Returns `(χ²/n, map, star, profrms, azrms, (χ²ᵥ₂/n, χ²ₜ₃ₐ/n, χ²ₜ₃ₚ/n))`.

THE SPLIT χ² IS THE POINT HERE, and it is the diagnostic α Cen could not provide. RW Cep
has closure phases with an rms of 118°, and a centrosymmetric model can only produce 0° or
180° — the parametric fit leaves χ²ₜ₃ₚ/n = 712. So if a regularizer combination annihilates
the surface map, the map becomes centrosymmetric and χ²ₜ₃ₚ must climb back toward that
value. Real structure surviving shows up as χ²ₜ₃ₚ staying low. No eyeballing of images
required: the closure phases measure directly whether genuine asymmetry is being destroyed.
"""
function reconstruct(ld1, rf, rv, tvw, ol = 0.0)
    p     = params(ld1)
    stars = create_star_multiepochs(tessels, p, [0.0])
    setup_oi!(data, stars)
    x0    = parametric_temperature_map(p, stars[1])
    regs  = Any[]
    tvw > 0 && push!(regs, [SMOOTH, tvw, smoothinfo, 1:length(x0)]) # omitted entirely at tvw = 0
    bins  = radflat_bins(stars[1]; nbins = NBINS)        # bins follow the PROJECTED geometry
    rf > 0 && push!(regs, ["radflat",   rf, bins, bins.idx])
    rv > 0 && push!(regs, ["radialvar", rv, bins, bins.idx])
    if ol > 0                                        # orthoLD: rank-1, no binning
        od = orthold_direction(stars[1], x0, p)
        push!(regs, ["orthold", ol, od, od.idx])
    end
    local x
    for it in (MAXIT, MAXIT ÷ 4, MAXIT ÷ 16, 500, 100)
        try
            x = image_reconstruct_oi(x0, data, stars; maxiter = it,
                                     regularizers = regs, verbose = false)
            break
        catch err
            err isa AssertionError || rethrow()   # VMLMB line search at full convergence
            it == 100 && rethrow()
        end
    end
    xb   = x[bins.idx]
    m    = mean(xb)
    prof = [mean(xb[bins.bin .== k]) for k in 1:NBINS] ./ m
    within = mean(mean(abs2, xb[bins.bin .== k] .- mean(xb[bins.bin .== k]))
                  for k in 1:NBINS if count(==(k), bins.bin) > 1)
    cv   = poly_to_cvis(x, stars[1])
    v2m, t3am, t3pm = cvis_to_obs(ComplexF64.(cv), d)
    split = (sum(abs2, (v2m .- d.v2) ./ d.v2_err) / d.nv2,
             sum(abs2, (t3am .- d.t3amp) ./ d.t3amp_err) / d.nt3amp,
             sum(abs2, mod360(t3pm .- d.t3phi) ./ d.t3phi_err) / d.nt3phi)
    return image_reconstruct_oi_chi2(x, data, stars, verbose = false)/nd, x, stars[1],
           sqrt(mean(abs2, prof .- 1)), sqrt(within)/m, split
end

"""
Indices where VMLMB plainly failed to reach the same basin as its neighbours.

Weak TV2 makes the problem badly conditioned and the odd start lands somewhere far worse —
one point at χ²/n = 3.22 against neighbours at 1.58, for instance. The curvature being
measured is worth ~0.1 in χ²/n end to end, so a single stuck run is 15x the signal. Flag
and exclude from the fit, but PRINT it: a stuck run is a fact about the reconstruction,
not noise to be quietly cleaned away.
"""
stuck(c) = findall(c .> minimum(c) + 3*(median(c) - minimum(c)) + 1e-9)

"""
    vertex(grid, c) -> (ld1*, σ_ld1, scatter, a, a/σ_a)

Fit `χ²/n = a·ld1² + b·ld1 + const` over the WHOLE grid by least squares; report the vertex,
the residual scatter, the curvature `a` and its significance `a/σ_a`.

The estimator matters more than it looks. The obvious choice — a parabola through the three
points around the minimum, with σ from Δχ² = 1 — is WRONG here, and confidently so: it
returns σ ≈ 0.02 while the scatter between ADJACENT ld1 values is 0.03–0.07 in χ²/n. Each
reconstruction settles into its own local minimum, so three neighbouring points carry at
least as much reconstruction noise as curvature, and the three-point fit reports that noise
as a precision. A global fit uses every point and takes σ FROM the residual scatter, so
reconstruction noise inflates the error bar instead of shrinking it.

`a/σ_a` is the number that answers the question this script asks. The degeneracy IS the
statement that χ²(ld1) has no minimum, so it lifts exactly when `a` becomes significantly
positive. A huge |ld1*| next to `a/σ_a ≈ 0` is not a measurement — it is a nearly straight
line whose extrapolated vertex ran off the end of the world.
"""
function vertex(grid, c)
    X = hcat(grid.^2, grid, ones(length(grid)))
    p = X \ c
    r = c .- X*p
    s2 = sum(abs2, r) / max(length(grid) - 3, 1)
    C  = s2 * inv(X'X)
    a, b = p[1], p[2]
    abs(a) < 1e-12 && return (NaN, NaN, sqrt(s2), a, 0.0)
    xv = -b/(2a)
    da, db = b/(2a^2), -1/(2a)                                   # ∂xv/∂a, ∂xv/∂b
    var = da^2*C[1,1] + db^2*C[2,2] + 2*da*db*C[1,2]
    return (xv, sqrt(max(var, 0)), sqrt(s2), a, a/sqrt(C[1,1]))
end

for tvw in TVWS
    @printf("\n%s\n%s weight = %.0e%s\n%s\n", "="^96, SMOOTH, tvw,
            tvw == 0 ? "   (OFF — only RADFLAT/RADIALVAR/orthoLD constrain the map)" : "",
            "="^96)
    @printf("%-10s", "regs\\ld1")
    for l in LDGRID; @printf("%7.1f", l); end
    @printf("  %9s %7s %8s %8s %8s   %7s %7s %8s\n",
            "ld1*", "a/sig", "scatter", "profrms", "azrms", "χ²ᵥ₂", "χ²ₜ₃ₐ", "χ²ₜ₃ₚ")
    for (nm, rf, rv, ol) in REGSETS
        c   = [reconstruct(l, rf, rv, tvw, ol)[1] for l in LDGRID]
        bad = stuck(c)
        keep = setdiff(1:length(c), bad)
        ld1s, σ, sc, a, sig = vertex(LDGRID[keep], c[keep])
        at = isfinite(ld1s) && 0 < ld1s < 5 ? ld1s : LDGRID[argmin(c)]
        _, _, _, prms, azrms, sp = reconstruct(at, rf, rv, tvw, ol)
        @printf("%-10s", nm)
        for v in c; @printf("%7.3f", v); end
        @printf("  %9s %7.1f %8.4f %8.4f %8.4f   %7.2f %7.2f %8.1f%s\n",
                sig > 2 ? @sprintf("%.2f±%.2f", ld1s, σ) : "unconstr",
                sig, sc, prms, azrms, sp[1], sp[2], sp[3],
                isempty(bad) ? "" : "  stuck at ld1=" *
                    join((@sprintf("%.2f", LDGRID[j]) for j in bad), ","))
    end
    println("  a/sig > ~2 with a > 0 ⇒ χ²(ld1) really has a minimum ⇒ the degeneracy has lifted.")
    println("  profrms = between-annuli structure (RADFLAT's target); azrms = within-annuli (RADIALVAR's).")
    @printf("  χ²ₜ₃ₚ IS THE TEST: a centrosymmetric map can only give 0/180°, and the parametric\n")
    @printf("  fit leaves χ²ₜ₃ₚ/n = %.0f. A row that drives the map flat must climb back toward that;\n",
            712)
    println("  one that keeps genuine asymmetry stays low. This is what α Cen could not measure.")
end

# --- 4. Images and the radial profile RADFLAT acts on ------------------------
if DOPLOT
    @eval using PythonPlot
    mkpath(joinpath(HERE, "results"))
    tvw = first(TVWS)
    # Use the fitted vertex only when it is actually inside the scanned range. With a flat
    # scan the quadratic is nearly a straight line and its vertex can land at ld1 = 385.
    function best_ld(rf, rv, ol)
        c = [reconstruct(l, rf, rv, tvw, ol)[1] for l in LDGRID]
        v = vertex(LDGRID, c)[1]
        return (isfinite(v) && first(LDGRID) <= v <= last(LDGRID)) ? v : LDGRID[argmin(c)]
    end
    maps = Dict{String,Any}()
    fig, ax = subplots(figsize = (7, 4.5))
    for (nm, rf, rv, ol) in REGSETS
        l = best_ld(rf, rv, ol)
        _, x, s, _, _, sp = reconstruct(l, rf, rv, tvw, ol)
        maps[nm] = (x = x, star = s, ld1 = l, t3p = sp[3])
        fmap, _ = plot2d(x, s, intensity = true, graticules = true, compass = true,
                         figtitle = @sprintf("RW Cep — %s (ld1=%.2f, χ²ₜ₃ₚ/n=%.0f)", nm, l, sp[3]))
        fmap.savefig(joinpath(HERE, "results", "rwcep_$nm.png"), dpi = 130, bbox_inches = "tight")
        pyplot.close(fmap)
        b = radflat_bins(s; nbins = NBINS); xb = x[b.idx]
        ax.plot([(k-0.5)/NBINS for k in 1:NBINS],
                [mean(xb[b.bin .== k]) for k in 1:NBINS] ./ mean(xb), "o-", label = nm)
    end
    ax.axhline(1.0, color = "k", ls = ":", lw = 1); ax.legend(); ax.grid(alpha = 0.3)
    ax.set_xlabel("projected radius ρ / ρ_max"); ax.set_ylabel("mean patch brightness / disk mean")
    ax.set_title("Radial profile — the quantity RADFLAT flattens")
    fig.savefig(joinpath(HERE, "results", "rwcep_radflat_profile.png"), dpi = 130, bbox_inches = "tight")
    pyplot.close("all")
    # How much of the reconstruction each regularizer set destroys, against the unregularized
    # map. Monnier's claim is that RADFLAT costs almost nothing in image terms — quantify it
    # rather than asserting it. Maps normalised to unit mean before comparing.
    ref = maps["none"].x ./ mean(maps["none"].x)
    @printf("\n%-10s %8s %10s %10s %10s\n", "regs", "ld1", "rms(Δ)", "Pearson r", "χ²ₜ₃ₚ/n")
    for (nm, _, _, _) in REGSETS
        y = maps[nm].x ./ mean(maps[nm].x)
        @printf("%-10s %8.2f %10.3f %10.4f %10.1f\n", nm, maps[nm].ld1,
                sqrt(mean(abs2, y .- ref)), cor(vec(ref), vec(y)), maps[nm].t3p)
    end
    println("rms(Δ) and r are against the UNREGULARIZED map: how much structure each set removes.")
    @info "wrote demos/results/rwcep_{none,radflat,radialvar,both,radflat_profile}.png"
end
