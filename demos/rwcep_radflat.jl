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
#     scan the LD coefficient, reconstruct at each value, and look at χ²(α).
#
# A FLAT χ²(α) *is* the degeneracy — every LD value fits equally well because the surface map
# absorbs the difference. If RADFLAT works, χ²(α) should develop a real minimum.
#
# NOT a reproduction of Monnier's numbers: SURFING has its own parameterisation, priors and
# error model. The question is whether the degeneracy lifts, not whether 0.49 comes back.
# RADFLAT is deliberately wrong for rotators — several epochs already break this degeneracy,
# and RADFLAT would then suppress real structure.

using ROTIR, Printf, Statistics

const HERE   = @__DIR__
const OIFITS = joinpath(HERE, "MIRCX_L3.022Dec23.RW_Cep.MIRCX_IDL.deepclean.AVG15m.oifits")
const NSIDE  = parse(Int, get(ENV, "NSIDE", "4"))          # 3..6 is the useful range
const NBINS  = parse(Int, get(ENV, "NBINS", "6"))          # Monnier uses 6
const MAXIT  = parse(Int, get(ENV, "MAXITER", "300"))
const TVW    = parse(Float64, get(ENV, "TVW", "1e-4"))     # baseline TV2, always on
const ALPHAS = parse.(Float64, split(get(ENV, "ALPHAS", "0,1e-3,1e-2,1e-1,1e0"), ','))
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
LDGRID  = parse.(Float64, split(get(ENV, "LDGRID",
              join([@sprintf("%.2f", α_fit*f) for f in (0.4,0.6,0.8,1.0,1.2,1.4,1.6)], ',')), ','))

params(ld1) = (surface_type = 0, radius = D_fit/2, tpole = 3900.0,   # RW Cep ≈ M0 Iab
               ldtype = 3, ld1 = ld1, ld2 = 0.0,
               inclination = 90.0, position_angle = 0.0, rotation_period = 1.0)

"Reconstruct at fixed LD coefficient `ld1` with RADFLAT weight `α` (α = 0 disables it)."
function reconstruct(ld1, α)
    p     = params(ld1)
    stars = create_star_multiepochs(tessels, p, [0.0])
    setup_oi!(data, stars)
    x0    = parametric_temperature_map(p, stars[1])
    regs  = Any[["tv2", TVW, tvinfo, 1:length(x0)]]      # TV2 always on: isolates RADFLAT
    if α > 0
        bins = radflat_bins(stars[1]; nbins = NBINS)      # bins follow the PROJECTED geometry
        push!(regs, ["radflat", α, bins, bins.idx])       # subset must be the one bins used
    end
    x = image_reconstruct_oi(x0, data, stars; maxiter = MAXIT,
                             regularizers = regs, verbose = false)
    return image_reconstruct_oi_chi2(x, data, stars, verbose = false)/nd, x, stars[1]
end

@printf("%-9s", "α \\ ld1")
for l in LDGRID; @printf(" %8.2f", l); end
@printf(" %9s %10s\n", "min at", "Δχ²/n")
res = Dict{Float64,Vector{Float64}}(); best = Dict{Float64,Float64}()
for α in ALPHAS
    c = [reconstruct(l, α)[1] for l in LDGRID]
    res[α] = c; k = argmin(c); best[α] = LDGRID[k]
    @printf("%-9.0e", α)
    for v in c; @printf(" %8.4f", v); end
    @printf(" %9.2f %10.4f\n", LDGRID[k], maximum(c) - minimum(c))
end

println("\nΔχ²/n is the SPREAD of reduced χ² across the LD range. Small ⇒ every LD value fits")
println("about equally well ⇒ the coefficient is unconstrained. It should GROW with α.")
c0 = res[first(ALPHAS)]; s0 = maximum(c0) - minimum(c0)
for α in ALPHAS
    c = res[α]
    @printf("  α = %-7.0e Δχ²/n = %8.4f (%5.1fx)  best ld1 = %.2f\n",
            α, maximum(c)-minimum(c), (maximum(c)-minimum(c))/max(s0,1e-12), best[α])
end

# --- 4. Images and the radial profile RADFLAT acts on ------------------------
if DOPLOT
    @eval using PythonPlot
    mkpath(joinpath(HERE, "results"))
    αoff, αon = first(ALPHAS), last(ALPHAS)
    _, xoff, s1 = reconstruct(best[αoff], αoff)
    _, xon,  s2 = reconstruct(best[αon],  αon)
    for (x, s, t) in ((xoff, s1, @sprintf("RW Cep — no RADFLAT (ld1=%.2f)", best[αoff])),
                      (xon,  s2, @sprintf("RW Cep — RADFLAT α=%.0e (ld1=%.2f)", αon, best[αon])))
        plot2d(x, s, intensity = true, graticules = true, compass = true, figtitle = t)
    end
    fig, ax = subplots(figsize = (7, 4.5))
    for (x, s, lab) in ((xoff, s1, "no RADFLAT"), (xon, s2, "RADFLAT"))
        b = radflat_bins(s; nbins = NBINS); xb = x[b.idx]
        ax.plot([(k-0.5)/NBINS for k in 1:NBINS],
                [mean(xb[b.bin .== k]) for k in 1:NBINS] ./ mean(xb), "o-", label = lab)
    end
    ax.axhline(1.0, color = "k", ls = ":", lw = 1); ax.legend(); ax.grid(alpha = 0.3)
    ax.set_xlabel("projected radius ρ / ρ_max"); ax.set_ylabel("mean patch brightness / disk mean")
    ax.set_title("Radial profile — the quantity RADFLAT flattens")
    fig.savefig(joinpath(HERE, "results", "rwcep_radflat_profile.png"), dpi = 130, bbox_inches = "tight")
    @info "wrote results/rwcep_radflat_profile.png"
end
