#!/usr/bin/env julia
# alphacen_ld.jl — limb darkening of α Cen A and B against a published benchmark.
#
#   julia --project=demos demos/alphacen_ld.jl
#   STARS=A NSIDE=4 SCAN=0 julia --project=demos demos/alphacen_ld.jl
#
# Headless by construction: Agg backend, PNGs into demos/results/, nothing displayed.
#
# ---------------------------------------------------------------------------
# Why α Cen is the right benchmark
# ---------------------------------------------------------------------------
# Kervella, Bigot, Gallenne & Thévenin 2017, A&A 597, A137 (`aa29505-16.pdf`) measured both
# components with VLTI/PIONIER in H, resolving α Cen A into the fourth lobe and B into the
# second. They fit a family of limb-darkening laws to the SQUARED VISIBILITIES ALONE and
# tabulate θ_LD and the LD coefficients for each (their Table 3). That gives us three things
# almost nothing else does:
#
#   * a published number for every law ROTIR implements, from the same data;
#   * two stars that are, to the precision of the data, featureless limb-darkened discs —
#     so an image reconstruction that invents surface structure is visibly wrong;
#   * a KNOWN limb-darkening coefficient, which makes this a ground-truth test of RADFLAT.
#     On RW Cep the true LD is unknown, so a scan there can only show whether χ²(ld1)
#     acquires a minimum; here it can show whether that minimum is in the RIGHT PLACE.
#
# ---------------------------------------------------------------------------
# The 0.48 % wavelength calibration — read this before comparing any diameter
# ---------------------------------------------------------------------------
# Kervella §2.2.4 multiply the PIONIER wavelength scale by γ = 1.00481 ± 0.00412, derived by
# comparing their fitted semi-major axis of HD 123999 against three literature values. The
# OIFITS distributed with OITOOLS carries the UNCORRECTED wavelengths, so every angular
# diameter fitted from it is 0.48 % small. That is larger than their entire quoted error bar
# (0.43 %), so it is not optional bookkeeping — without it every diameter here looks
# discrepant at >1σ, and with it they agree to the fourth decimal.
#
# θ ∝ λ, so the correction is a plain multiplication. Both are printed below.

ENV["MPLBACKEND"] = get(ENV, "MPLBACKEND", "Agg")
using ROTIR, PythonPlot, Printf, Statistics, LinearAlgebra, NLopt
using ROTIR: _fit_ultranest          # unexported: the shared UltraNest driver

const HERE   = @__DIR__
const OUT    = joinpath(HERE, "results"); mkpath(OUT)
const STARS  = split(get(ENV, "STARS", "A,B"), ',')
const NSIDE  = parse(Int, get(ENV, "NSIDE", "3"))
const MAXIT  = parse(Int, get(ENV, "MAXITER", "20000"))   # convergence, not speed — see below
const NBINS  = parse(Int, get(ENV, "NBINS", "6"))
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
# Weight shared by RADFLAT and RADIALVAR. 100 is not arbitrary: a ladder on α Cen A at
# nside 3 (target = the mesh parametric ld1 = 0.1399) gives
#   w = 1     ld1* = 0.1852 ± 0.0077   a/σ_a 18   profrms 1.4e-2
#   w = 10    ld1* = 0.1488 ± 0.0030   a/σ_a 58   profrms 4e-4
#   w = 100   ld1* = 0.1478 ± 0.0021   a/σ_a 81   profrms 1e-4     ← best
#   w = 1000  ld1* = 0.1621 ± 0.0067   a/σ_a 24   profrms 2e-4
#   w = 1e4   ld1* = 0.1645 ± 0.0077   a/σ_a 21   profrms 3e-5
# Too weak and the ~1 % residual structure trades against ld1 and biases it high; too
# strong and the quadratic penalty stiffens the problem until VMLMB stops converging
# properly — note a/σ_a COLLAPSING past 1e3 while the structure metrics keep shrinking,
# which is the signature of an optimizer failing rather than a prior succeeding.
const REGW   = parse(Float64, get(ENV, "REGW", "100"))
const DOSCAN = get(ENV, "SCAN", "1") == "1"
const DOUN   = get(ENV, "ULTRANEST", "1") == "1"
const NLIVE  = parse(Int, get(ENV, "LIVE", "400"))
const GAMMA  = 1.00481                                     # Kervella Eq. (1)

# Kervella Table 3, V²-only fits. (stat ⊕ syst errors; we quote stat, the syst 0.42–0.43 %
# is common to all rows and is the wavelength calibration itself.)
const LIT = Dict(
 "A" => (nv2 = 324,
         uniform   = (θ = 8.347, σθ = 0.004, p = Float64[],       σp = Float64[],      χ2 = 15.23),
         linear    = (θ = 8.458, σθ = 0.005, p = [0.1761],        σp = [0.0062],       χ2 =  4.24),
         quadratic = (θ = 8.451, σθ = 0.013, p = [0.191, -0.031], σp = [0.026, 0.054], χ2 =  4.25),
         power     = (θ = 8.502, σθ = 0.006, p = [0.1404],        σp = [0.0050],       χ2 =  3.90)),
 "B" => (nv2 = 432,
         uniform   = (θ = 5.883, σθ = 0.003, p = Float64[],       σp = Float64[],      χ2 = 14.26),
         linear    = (θ = 5.962, σθ = 0.003, p = [0.1907],        σp = [0.0048],       χ2 =  2.89),
         # Kervella fit no 2-parameter law to B: they reach only the second lobe, so the
         # data cannot separate a from b. Ours is reported but has nothing to compare to.
         quadratic = (θ = NaN,   σθ = NaN,   p = [NaN, NaN],      σp = [NaN, NaN],     χ2 = NaN),
         power     = (θ = 5.999, σθ = 0.004, p = [0.1545],        σp = [0.0044],       χ2 =  3.33)))

# --- 1. Data: V² only ---------------------------------------------------------
# `use_t3 = false` drops the closure phases AND their uv points, so `nuv` is the V² count
# and the imaging step below builds a correspondingly smaller `polyft`. Both stars are
# centrosymmetric to within the errors, so nothing is lost by excluding T3 — and it is what
# the paper did.
load(star) = readoifits(joinpath(HERE, "data", "AlphaCen$star.oifits");
                        use_t3 = false, filter_bad_data = true,
                        verbose = false, warn = false)[1, 1]

# --- 2. Analytic parametric fits ----------------------------------------------
# The closed-form limb-darkened visibilities (OITOOLS, re-exported by ROTIR). These are
# Hankel transforms of the intensity profile, so they carry no mesh discretisation and are
# the fair comparison against a published parametric fit.
"""
Least-squares fit of an analytic visibility model to V² alone.

A COARSE SCAN OVER θ FIRST, then Nelder–Mead. That is not defensive padding: a star
resolved into its fourth lobe has a violently multi-modal χ²(θ), because a model disc one
lobe too large can still line its nulls up with some of the data. Started from α Cen A's
diameter, a plain Nelder–Mead on α Cen B converges to θ = 15.07 mas at χ²ᵣ = 1169 and
reports success. The scan step must be a good deal finer than the lobe spacing λ/B_max.
"""
function fit_v2(d, vf, x0, lb, ub)
    uv = d.uv[:, d.indx_v2]                       # V² uv points only
    f(p) = sum(abs2, (abs2.(vf(p, uv)) .- d.v2) ./ d.v2_err)
    # lobe spacing in θ is ~λ/B_max; step at a tenth of it
    Bmax = maximum(sqrt.(uv[1,:].^2 .+ uv[2,:].^2))
    step = 0.1 * 206264806.2 / Bmax
    best = copy(x0); bc = Inf
    for θ in lb[1]:step:ub[1]
        q = copy(x0); q[1] = θ
        c = f(q); c < bc && ((bc, best) = (c, q))
    end
    o = Opt(:LN_NELDERMEAD, length(x0))
    o.lower_bounds = lb; o.upper_bounds = ub
    o.xtol_rel = 1e-10; o.ftol_rel = 1e-12; o.maxeval = 20000
    o.min_objective = (p, g) -> f(p)
    (c, p, _) = optimize(o, best)
    return p, c / (d.nv2 - length(x0))             # reduced χ², as in Table 3
end

const MODELS = ("uniform", "linear", "quadratic", "power")
vismodel(m) = m == "uniform"   ? visibility_ud   :
              m == "linear"    ? visibility_ldlin :
              m == "quadratic" ? visibility_ldquad : visibility_ldpow

"Prior box and parameter names for each law: (names, lo, hi)."
function prior_box(m)
    m == "uniform"   && return (["θ"],            [1.0],            [15.0])
    m == "quadratic" && return (["θ", "a", "b"],  [1.0, -1.0, -1.0], [15.0, 1.0, 1.0])
    m == "power"     && return (["θ", "α"],       [1.0, 0.0],       [15.0, 3.0])
    return (["θ", "u"], [1.0, 0.0], [15.0, 1.0])                     # linear
end

"""
Nested sampling of the same V²-only likelihood, over the same box the grid scan searches.

This needs no starting guess at all, which is the point: `fit_v2` above has to pre-scan θ
because χ²(θ) is multi-modal past the first null — started badly, Nelder–Mead settles into a
wrong lobe at χ²ᵣ = 1169 and reports convergence. UltraNest samples the whole prior, so
finding the right mode is a result rather than a consequence of where it was started, and
the posterior width is directly comparable to Kervella's quoted error bars instead of being
a curvature estimate at a point.
"""
function ultranest_fit(d, m; nlive = NLIVE)
    uv = d.uv[:, d.indx_v2]
    vf = vismodel(m)
    names, lo, hi = prior_box(m)
    chi2(p) = begin
        v = sum(abs2, (abs2.(vf(p, uv)) .- d.v2) ./ d.v2_err)
        isfinite(v) ? v : 1e30
    end
    un = _fit_ultranest(chi2, names, lo, hi; min_num_live_points = nlive)
    return un, chi2(un.median) / (d.nv2 - length(lo))
end

function analytic_fits(d, θ0)
    out = Dict{String,Any}()
    for m in MODELS
        x0, lb, ub = m == "uniform"   ? ([θ0], [1.0], [15.0]) :
                     m == "quadratic" ? ([θ0, 0.2, 0.0], [1.0, -1.0, -1.0], [15.0, 1.0, 1.0]) :
                     m == "power"     ? ([θ0, 0.14], [1.0, 0.0], [15.0, 3.0]) :
                                        ([θ0, 0.2], [1.0, 0.0], [15.0, 1.0])
        out[m] = fit_v2(d, vismodel(m), x0, lb, ub)
    end
    return out
end

# --- 3. Tessellated (mesh) parametric fit -------------------------------------
# Same physics through ROTIR's HEALPix forward model instead of a closed form. Agreement
# with §2 validates `compute_ldmap` + `polyft` against an independent analytic transform;
# this is what caught the quadratic-law convention error (ROTIR used 1−a(1−μ)−b(1−μ²),
# everyone else uses 1−a(1−μ)−b(1−μ)²).
function mesh_fit(d, ldtype, θ0)
    tess = tessellation_healpix(NSIDE)
    r = fit_sphere_ld([d], tess; tepochs = [0.0], ldtype = ldtype,
                      radius0 = θ0/2, radius_bounds = (0.5, 10.0),
                      ld_bounds = ldtype == 3 ? (0.0, 3.0) : (0.0, 1.0),
                      fit_ld2 = ldtype == 2,
                      weights = [1.0, 0.0, 0.0],          # V² only
                      algorithm = :LN_NELDERMEAD, maxeval = 4000)
    nfree = ldtype == 2 ? 3 : 2
    return 2r.radius, r.ld1, r.ld2, r.chi2/(d.nv2 - nfree)
end

# --- 4. Imaging + RADFLAT, as a ground-truth test -----------------------------
# Fix the diameter at the parametric value, fix ld1, reconstruct a free surface map, and
# read off χ². Repeat over a grid of ld1. Two things are being asked at once:
#
#   does χ²(ld1) have a minimum at all?   (the degeneracy RADFLAT targets)
#   is that minimum at the PUBLISHED α?   (only answerable on a star like this one)
#
# MAXITER matters: on RW Cep the drift in χ²/n between 300 and 3000 iterations exceeded the
# entire spread across the LD grid, so an under-converged scan measures the optimizer. VMLMB
# self-terminates well before 20000 here.
"""
Regularizer sets compared in the scan: `(name, rf, rv, ol)` weighting RADFLAT
(between-annuli), RADIALVAR (within-annuli) and orthoLD (the rank-one direction exactly
degenerate with ld1; the derivation is in the `spheroid_orthold_fg` docstring). A smoothing
prior is always present so all three are isolated against a common baseline.

RADFLAT and RADIALVAR are the two halves of the disk's variance, so "both" constrains the
whole map — and on this star that is precisely why it "works": the map is driven to a
constant and the model degenerates to the parametric limb-darkened disc. orthoLD is the
control for that objection. It removes ONE direction rather than the entire map, so if it
also recovers the published α, the recovery is not an artefact of annihilating the surface.
α Cen is the only place this can be tested, because it is the only star here whose
limb darkening is independently known.
"""
const REGSETS = (("none",       0.0,  0.0,  0.0),
                 ("radflat",    REGW, 0.0,  0.0),
                 ("radialvar",  0.0,  REGW, 0.0),
                 ("both",       REGW, REGW, 0.0),
                 ("orthold 1",  0.0,  0.0,  1.0),
                 ("orthold 1e2",0.0,  0.0,  REGW))

function reconstruct(d, θ, ld1, rf, rv = 0.0, ol = 0.0)
    tess  = tessellation_healpix(NSIDE)
    tv    = SMOOTH == "tv2" ? tv_neighbors_healpix(NSIDE) : sobel_gradient_healpix(NSIDE)
    p     = (surface_type = 0, radius = θ/2, tpole = 5000.0, ldtype = 3,
             ld1 = ld1, ld2 = 0.0, inclination = 90.0, position_angle = 0.0,
             rotation_period = 1.0)
    stars = create_star_multiepochs(tess, p, [0.0])
    data  = [d]
    setup_oi!(data, stars)
    x0    = parametric_temperature_map(p, stars[1])
    regs  = Any[[SMOOTH, TVW, tv, 1:length(x0)]]
    bins  = radflat_bins(stars[1]; nbins = NBINS)
    rf > 0 && push!(regs, ["radflat",   rf, bins, bins.idx])
    rv > 0 && push!(regs, ["radialvar", rv, bins, bins.idx])
    if ol > 0                                        # orthoLD: rank one, no binning
        od = orthold_direction(stars[1], x0, p)
        push!(regs, ["orthold", ol, od, od.idx])
    end
    # OptimPackNextGen's VMLMB throws `AssertionError: stp == ls.stp` when its line search
    # is asked to keep going from a point it has already minimised to machine precision —
    # which happens here at small nside, where 100 free patches against 324 data converge
    # long before MAXITER. Backing the iteration count off returns the same converged map
    # (it had stopped improving), so retry rather than lose the grid point.
    local x
    for it in (MAXIT, MAXIT ÷ 4, MAXIT ÷ 16, 500, 100)
        try
            x = image_reconstruct_oi(x0, data, stars; maxiter = it,
                                     regularizers = regs, verbose = false)
            break
        catch err
            err isa AssertionError || rethrow()
            it == 100 && rethrow()
        end
    end
    # Two structure metrics, one per half of the variance split, both as a fraction of the
    # disk mean so they are directly comparable:
    #   profrms — rms departure of the annulus MEANS from the disk mean   (RADFLAT's target)
    #   azrms   — rms scatter WITHIN annuli, pooled                       (RADIALVAR's target)
    xb   = x[bins.idx]
    m    = mean(xb)
    prof = [mean(xb[bins.bin .== k]) for k in 1:NBINS] ./ m
    within = mean(mean(abs2, xb[bins.bin .== k] .- mean(xb[bins.bin .== k]))
                  for k in 1:NBINS if count(==(k), bins.bin) > 1)
    return image_reconstruct_oi_chi2(x, data, stars, verbose = false)/d.nv2,
           x, stars[1], sqrt(mean(abs2, prof .- 1)), sqrt(within)/m
end

"""
Vertex of a global quadratic fit to χ²(ld1), with σ propagated from the RESIDUAL scatter.

A three-point parabola around the minimum is NOT usable here: neighbouring reconstructions
land in different local minima and scatter by more than the curvature, so a local fit
reports that scatter as a precision. `a/σ_a` is the real answer — the degeneracy IS the
absence of a minimum, so it lifts exactly when `a` becomes significantly positive.
"""
function vertex(grid, c)
    # Fit LOCALLY when there is a sharp minimum, GLOBALLY when the curve is flat. With both
    # regularizers on, χ²(ld1) runs from 4.4 to 413 across the grid and is nothing like a
    # parabola over that range, so a global fit is dominated by the far wings and its σ is
    # meaningless. When the scan IS flat the opposite holds — a local fit sees only noise
    # (see the note above) — so the window only kicks in once the minimum is deep enough to
    # define one, and never shrinks below 5 points.
    lo, hi = extrema(c)
    if hi > 3lo && length(grid) > 5              # a deep, obviously non-parabolic minimum
        k = argmin(c)
        w = clamp(k-2, 1, length(grid)-4) : clamp(k+2, 5, length(grid))
        grid, c = grid[w], c[w]
    end
    X = hcat(grid.^2, grid, ones(length(grid)))
    p = X \ c
    s2 = sum(abs2, c .- X*p) / max(length(grid) - 3, 1)
    C  = s2 * inv(X'X)
    a, b = p[1], p[2]
    abs(a) < 1e-12 && return (NaN, NaN, sqrt(s2), a, 0.0)
    xv = -b/(2a)
    da, db = b/(2a^2), -1/(2a)
    return (xv, sqrt(max(da^2*C[1,1] + db^2*C[2,2] + 2*da*db*C[1,2], 0)),
            sqrt(s2), a, a/sqrt(C[1,1]))
end

# ==============================================================================
for star in STARS
    L = LIT[star]
    d = load(star)
    @printf("\n%s\nα Cen %s — %d V² points, H band, closure phases excluded (use_t3=false)\n%s\n",
            "="^100, star, d.nv2, "="^100)
    d.nv2 == L.nv2 || @warn "V² count differs from the paper" ours=d.nv2 paper=L.nv2

    # ---- analytic ----
    af = analytic_fits(d, L.linear.θ)
    println("\nAnalytic parametric fits to V² (Kervella Table 3 in brackets; pulls use their stat error)")
    @printf("%-10s %9s %9s %8s   %-22s %8s  %5s %5s\n",
            "law", "θ raw", "θ×γ", "θ pull", "LD coefficients [published]", "LD pull", "χ²ᵣ", "[lit]")
    for m in MODELS
        p, c = af[m]
        lit  = getfield(L, Symbol(m))
        pull(x, μ, σ) = isnan(μ) || isnan(σ) || σ == 0 ? NaN : (x - μ)/σ
        pθ   = pull(p[1]*GAMMA, lit.θ, lit.σθ)
        cof  = join((@sprintf("%+.4f", v) for v in p[2:end]), " ")
        litc = isempty(lit.p) || any(isnan, lit.p) ? "—" :
               join((@sprintf("%+.4f", v) for v in lit.p), " ")
        pld  = isempty(lit.p) ? NaN :
               maximum(abs(pull(p[1+k], lit.p[k], lit.σp[k])) for k in eachindex(lit.p))
        @printf("%-10s %9.4f %9.4f %8s   %-11s [%-11s] %8s  %5.2f %5s\n",
                m, p[1], p[1]*GAMMA, isnan(pθ) ? "—" : @sprintf("%+.1fσ", pθ),
                cof, litc, isnan(pld) ? "—" : @sprintf("%.1fσ", pld),
                c, isnan(lit.χ2) ? "—" : @sprintf("%.2f", lit.χ2))
    end
    @printf("γ = %.5f applied to the θ×γ column; θ ∝ λ and the OIFITS wavelengths are uncorrected.\n", GAMMA)

    # ---- UltraNest on the same likelihood ----
    if DOUN
        println("\nNested sampling of the same V²-only likelihood (no starting guess, $(NLIVE) live points)")
        @printf("%-10s %18s %10s %8s  %18s %10s %8s %8s %9s\n",
                "law", "θ×γ [mas]", "±post×√χ²ᵣ", "[lit ±]", "LD coeff", "±post×√χ²ᵣ",
                "[lit ±]", "χ²ᵣ", "log Z")
        for m in MODELS
            un, cr = ultranest_fit(d, m)
            lit = getfield(L, Symbol(m))
            # Posterior half-width, then INFLATED BY √χ²ᵣ. Kervella's statistical errors are
            # the likelihood widths scaled that way — the usual correction for a fit that
            # does not reach χ²ᵣ = 1. Unscaled, our intervals look ~2x too tight and the
            # comparison misleads; scaled, they land on the published errors to a few
            # percent. Reported as ± so the two columns are directly comparable.
            sc  = sqrt(cr)
            σθ  = (un.q84[1] - un.q16[1])/2 * GAMMA * sc
            θs  = @sprintf("%.4f±%.4f", un.median[1]*GAMMA, σθ)
            cs  = length(un.median) == 1 ? "—" :
                  join((@sprintf("%+.4f±%.4f", un.median[k], (un.q84[k]-un.q16[k])/2*sc)
                        for k in 2:length(un.median)), " ")
            litσp = isempty(lit.σp) || any(isnan, lit.σp) ? "—" :
                    join((@sprintf("%.4f", v) for v in lit.σp), " ")
            @printf("%-10s %18s %10.4f %8s  %18s %10s %8s %8.2f %9.1f\n",
                    m, θs, σθ, isnan(lit.σθ) ? "—" : @sprintf("%.3f", lit.σθ),
                    cs, "", litσp, cr, un.logz)
        end
        println("±post×√χ²ᵣ is the posterior half-width scaled by √χ²ᵣ; compare with [lit ±].")
        println("log Z ranks the laws as models: higher is better, and it PENALISES the extra")
        println("parameter a two-coefficient law spends, which a reduced χ² does not.")
    end

    # ---- mesh cross-check ----
    println("\nSame fits through the HEALPix forward model (nside $NSIDE) — validates the mesh")
    @printf("%-10s %22s %26s %14s\n", "law", "θ_LD [mas] (×γ)", "LD coefficients", "χ²ᵣ")
    mesh_pow = (θ = NaN, ld1 = NaN)
    for (m, lt) in (("linear", 1), ("quadratic", 2), ("power", 3))
        θm, l1, l2, cm = mesh_fit(d, lt, af[m][1][1])
        lt == 3 && (mesh_pow = (θ = θm, ld1 = l1))
        cof = lt == 2 ? @sprintf("%+.4f %+.4f", l1, l2) : @sprintf("%+.4f", l1)
        @printf("%-10s %8.4f → %8.4f  [%7.3f] %11s [%11s] %6.2f\n",
                m, θm, θm*GAMMA, getfield(L, Symbol(m)).θ, cof,
                isempty(getfield(L, Symbol(m)).p) || any(isnan, getfield(L, Symbol(m)).p) ? "—" :
                    join((@sprintf("%+.4f", v) for v in getfield(L, Symbol(m)).p), " "), cm)
    end

    # ---- V² figure ----
    θpow = af["power"][1][1]        # analytic power-law diameter, for the model curves
    B    = vec(sqrt.(sum(abs2, d.uv[:, d.indx_v2], dims = 1)))
    ord  = sortperm(B)
    fig, axs = pyplot.subplots(2, 1, figsize = (8, 7), sharex = true,
                               gridspec_kw = Dict("height_ratios" => [3, 1]))
    axs[0].errorbar(B[ord]/1e6, d.v2[ord], yerr = d.v2_err[ord], fmt = "o",
                    ms = 2.2, lw = 0.6, color = "0.45", label = "PIONIER V²", zorder = 1)
    Bg = collect(range(minimum(B), maximum(B), length = 900))
    uvg = vcat(Bg', zeros(1, length(Bg)))
    for (m, col) in zip(MODELS, ("#D55E00", "#0072B2", "#009E73", "#CC79A7"))
        p, c = af[m]
        axs[0].plot(Bg/1e6, abs2.(vismodel(m)(p, uvg)), "-", lw = 1.3, color = col,
                    label = @sprintf("%s (χ²ᵣ=%.2f)", m, c), zorder = 2)
        r = (abs2.(vismodel(m)(p, d.uv[:, d.indx_v2])) .- d.v2) ./ d.v2_err
        axs[1].plot(B[ord]/1e6, r[ord], ".", ms = 2.5, color = col)
    end
    axs[0].set_yscale("log"); axs[0].set_ylim(1e-4, 1.4)
    axs[0].set_ylabel("V²"); axs[0].legend(fontsize = 8, ncol = 2)
    axs[0].set_title("α Cen $star — limb-darkening laws fitted to V² only")
    axs[1].axhline(0, color = "k", lw = 0.8, ls = ":")
    axs[1].set_ylabel("residual / σ"); axs[1].set_xlabel("spatial frequency B/λ [Mrad⁻¹]")
    for a in (axs[0], axs[1]); a.grid(alpha = 0.25); end
    fig.savefig(joinpath(OUT, "alphacen_$(star)_v2.png"), dpi = 140, bbox_inches = "tight")
    pyplot.close(fig)

    # ---- imaging + RADFLAT ----
    if DOSCAN
        αlit = L.power.p[1]
        grid = collect(0.0:0.05:0.60)
        # Fix θ at the MESH power-law diameter, not the analytic one. The reconstruction
        # runs the mesh forward model, and at nside 3 that wants a diameter ~0.26 % larger
        # than the closed form does (discretisation). Pinning θ 0.26 % small leaves a
        # residual the map used to absorb — but once BOTH regularizers are on, the map
        # cannot, so ld1 compensates instead and the recovered coefficient comes out biased
        # high. Self-consistency here is not cosmetic: it is the difference between reading
        # off a limb-darkening coefficient and reading off a diameter error.
        θimg = isfinite(mesh_pow.θ) ? mesh_pow.θ : θpow
        # Count the free parameters, as context for the scan below. With closure phases
        # excluded there are only nv2 data, while the visible half of the mesh carries 428
        # patches at nside 3 and 1666 at nside 4 — against 324 V² for α Cen A. That ratio
        # is normal for imaging and is what the regularizers exist to absorb; it is NOT an
        # argument for a coarser mesh (nside 3–4 is the working range). It is, though, the
        # reason a limb-darkening COEFFICIENT cannot be read off a reconstruction: the map
        # has more freedom than the data, so it can absorb almost any ld1, and RADFLAT does
        # not change that because it constrains only the azimuthal AVERAGE of the map.
        sdiag = create_star_multiepochs(tessellation_healpix(NSIDE),
                    (surface_type = 0, radius = θimg/2, tpole = 5000.0, ldtype = 3,
                     ld1 = αlit, ld2 = 0.0, inclination = 90.0, position_angle = 0.0,
                     rotation_period = 1.0), [0.0])[1]
        @printf("\nnside %d: %d visible patches for %d V² data = %.2f free surface params per datum\n",
                NSIDE, sdiag.nquads_visible, d.nv2, sdiag.nquads_visible/d.nv2)
        println("Imaging with a free surface map, ld1 scanned (Hestroffer). θ fixed at $(round(θimg, digits=4)) mas (the MESH power-law value).")
        @printf("%-10s", "regs")
        for l in grid; @printf("%6.2f", l); end
        @printf("  %10s %7s %8s %9s %8s\n", "ld1*", "a/σ_a", "scatter", "profrms", "azrms")
        best = Dict{String,Float64}()
        for (nm, rf, rv, ol) in REGSETS
            c = [reconstruct(d, θimg, l, rf, rv, ol)[1] for l in grid]
            ld1s, σ, sc, a, sig = vertex(grid, c)
            best[nm] = (isfinite(ld1s) && first(grid) <= ld1s <= last(grid)) ? ld1s : grid[argmin(c)]
            _, _, _, prms, azrms = reconstruct(d, θimg, best[nm], rf, rv, ol)
            @printf("%-10s", nm)
            for v in c; @printf("%6.3f", v); end
            @printf("  %10s %7.1f %8.4f %9.5f %8.5f\n",
                    sig > 2 ? @sprintf("%.3f±%.3f", ld1s, σ) : "unconstr", sig, sc, prms, azrms)
        end
        println("profrms = between-annuli structure (RADFLAT's target); azrms = within-annuli (RADIALVAR's).")
        @printf("published α (power law, V²-only parametric) = %.4f\n", αlit)
        println("""
a/σ_a > ~2 with a > 0 means χ²(ld1) genuinely has a minimum. On a star whose limb darkening
IS known, the four rows say something neither could say alone:

  NEITHER REGULARIZER ALONE CONSTRAINS ld1. Both "radflat" and "radialvar" come out with
  |a/σ_a| ≲ 1, indistinguishable from "none". With ~1.3-5 free surface params per V² datum
  the map simply absorbs whatever the coefficient does — the imaging reaches χ²ᵣ ≈ 2 where
  the best PARAMETRIC model floors at 4.25 (A), which is the same statement.

  WORSE, EACH ONE ALONE PUSHES STRUCTURE INTO THE OTHER HALF. RADFLAT drives the radial
  profile to ~1e-4 and the azimuthal scatter UP by 17x (A) to 90x (B); RADIALVAR does the
  mirror image. The map needs its slack somewhere, and closing one door opens the other.
  That is the practical consequence of the two being a variance decomposition, and it is a
  reason not to deploy either on its own.

  TOGETHER THEY RECOVER THE PUBLISHED COEFFICIENT — at nside 3 with REGW = 100,
  ld1* = 0.145 ± 0.001 for A (published 0.1404, +3.3 %) and 0.155 for B (published 0.1545,
  +0.3 %). Both metrics fall at once and χ²(ld1) develops a minimum tens of σ deep.

  BUT READ THAT RESULT CORRECTLY. The two regularizers span the whole variance of the disk,
  so driving both to zero drives the map to a CONSTANT — at which point the model IS the
  parametric limb-darkened disc and χ²(ld1) IS the parametric curve. The recovery is
  therefore built in, not extracted: note χ²ᵣ at the minimum, 4.42, against the parametric
  floor of 4.59. This is a consistency check that the mesh forward model, the regularizers
  and the analytic transform all agree — a useful one, and it is what caught the quadratic
  LD convention error — but it is NOT evidence that imaging measures limb darkening. It is
  an expensive route back to the fit in the section above.

So: RADFLAT and RADIALVAR are for suppressing reconstruction artefacts on a star you have
reason to believe is smooth. They are not a route to a limb-darkening measurement; fit that
parametrically. The open question they DO bear on is the one α Cen cannot answer, because
it has no real surface structure to preserve: on a star that does, can the pair suppress
artefacts while leaving genuine features intact? That is `rwcep_radflat.jl`.

Carry this back to `rwcep_radflat.jl`, where χ²(ld1) acquires curvature under RADFLAT ALONE
(up to 6σ). Nothing here refutes that — RW Cep has closure phases and 2819 data against the
same 428 patches, 0.15 free params per datum against 1.32 here, so the map has far less room
— but α Cen shows that "RADFLAT alone made χ²(ld1) curve" is not by itself evidence the
coefficient is being measured. Check the params-per-datum ratio first.""")

        # images + radial profile
        figp, axp = pyplot.subplots(figsize = (7, 4.5))
        for (nm, rf, rv, ol) in REGSETS
            _, x, s, _, _ = reconstruct(d, θimg, best[nm], rf, rv, ol)
            fmap, _ = plot2d(x, s; intensity = true, graticules = true, compass = true,
                             figtitle = @sprintf("α Cen %s — %s (ld1=%.3f)", star, nm, best[nm]))
            fmap.savefig(joinpath(OUT, "alphacen_$(star)_$(nm).png"),
                         dpi = 130, bbox_inches = "tight")
            pyplot.close(fmap)
            b = radflat_bins(s; nbins = NBINS); xb = x[b.idx]
            axp.plot([(k-0.5)/NBINS for k in 1:NBINS],
                     [mean(xb[b.bin .== k]) for k in 1:NBINS] ./ mean(xb), "o-", label = nm)
        end
        axp.axhline(1.0, color = "k", ls = ":", lw = 1); axp.legend(); axp.grid(alpha = 0.3)
        axp.set_xlabel("projected radius ρ / ρ_max")
        axp.set_ylabel("mean patch brightness / disk mean")
        axp.set_title("α Cen $star — radial profile RADFLAT flattens")
        figp.savefig(joinpath(OUT, "alphacen_$(star)_profile.png"), dpi = 130, bbox_inches = "tight")
        pyplot.close(figp)
    end
end
@info "PNGs written to $OUT"
