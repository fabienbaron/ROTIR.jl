# beta Lyrae — relative-orbit figure from a STORED fit
#
#   julia --project=demos demos/betlyr/betlyr_orbit_figure.jl
#
# Reads `results/ultranest_best.txt` and `results/ultranest_posterior.txt` and redraws the
# sky-projected relative orbit, so the documentation figure can be regenerated without
# re-running the sampler (the fit itself is hours; this is seconds). The plotting recipe is
# the same one `betlyr_orbit_fit_pigeons.jl` uses inline, kept in sync deliberately.
#
# It also re-evaluates chi2_split at both the fitted and the literature parameters, so the
# numbers quoted in docs/src/guides/orbits.md come from the stored posterior rather than
# from a transcript.
#
# DATASET defaults to "old" to match the stored fit (Zhao et al.'s 2006-2007 MIRC nights).

get!(ENV, "DATASET", "old")
get!(ENV, "MPLBACKEND", "Agg")

using ROTIR, Printf, Statistics, DelimitedFiles, PythonPlot, PythonCall, Dates

const HERE = @__DIR__
include(joinpath(HERE, "betlyr_model.jl"))

const RES = joinpath(HERE, "results")
const OUT = get(ENV, "OUT", joinpath(HERE, "..", "..", "docs", "src", "assets"))
isfile(joinpath(RES, "ultranest_best.txt")) ||
    error("no stored fit at $RES — run betlyr_orbit_fit_ultranest.jl first")

# --- 1. Rebuild the full 11-vector from the 10 sampled parameters -------------
# `dP_dd` (index 5) is held fixed at its literature value and is therefore absent from the
# results file, so the stored names are matched back by NAME rather than by position.
best = readdlm(joinpath(RES, "ultranest_best.txt"))
θ_fit = copy(THETA_LIT)
for r in 1:size(best, 1)
    k = findfirst(==(String(best[r, 1])), PNAMES)
    k === nothing && error("unknown parameter \"$(best[r,1])\" in ultranest_best.txt")
    θ_fit[k] = Float64(best[r, 2])
end
sampled = [findfirst(==(String(best[r, 1])), PNAMES) for r in 1:size(best, 1)]

post = readdlm(joinpath(RES, "ultranest_posterior.txt"))
@printf("stored posterior: %d samples x %d sampled parameters\n", size(post, 1), size(post, 2))

# --- 2. Re-verify the fit quality --------------------------------------------
cl = chi2_split(THETA_LIT)
cf = chi2_split(θ_fit)
lab = ("V²", "T3amp", "T3φ")
println("\nreduced chi2 (n = $NTOT: $NV2 V², $NT3A T3amp, $NT3P T3φ)")
@printf("%-10s %12s %12s\n", "", "literature", "fitted")
for (j, n) in enumerate((NV2, NT3A, NT3P))
    @printf("%-10s %12.2f %12.2f\n", lab[j], cl[j]/n, cf[j]/n)
end
@printf("%-10s %12.2f %12.2f\n", "total", sum(cl)/NTOT, sum(cf)/NTOT)

println("\nelement                literature        fitted     posterior sd")
for k in sampled
    c = findfirst(==(k), sampled)
    @printf("%-16s %14.6f %14.6f %14.6g\n", PNAMES[k], THETA_LIT[k], θ_fit[k], std(post[:, c]))
end

# --- 3. The figure ------------------------------------------------------------
# Relative orbit of the disc about the donor, projected on sky. P is now a FITTED parameter
# (index 11) rather than the fixed P_ORB the Pigeons driver assumed, so it is read from θ.
rd(θ, t) = orbit_to_rotir_offset((i = θ[2], Ω = θ[3], ω = OMEGA_PERI, P = θ[11],
                                  a = θ[1], e = E_ORB, T0 = θ[4], q = Q_BIN,
                                  dP = θ[5], dω = 0.0), t)

# i = 91.2 deg, so the projected ellipse is ~3:1 and a square canvas wastes most of the
# frame. Legend goes BELOW the axes: at this aspect every in-axes corner overlaps the orbit.
fig, ax = subplots(figsize = (9, 4.2)); ax.set_aspect("equal")
tt = θ_fit[4] .+ range(0, θ_fit[11], length = 400)

# Posterior spread first, so it sits UNDER the two curves. With 10 parameters this tight the
# draws overlap almost exactly — that visual degeneracy is the point: it is what a converged
# posterior looks like, and a fan of visibly distinct ellipses would mean the opposite.
ndraw = min(60, size(post, 1))
for r in round.(Int, range(1, size(post, 1), length = ndraw))
    θd = copy(THETA_LIT); θd[sampled] .= post[r, :]
    xy = [rd(θd, t) for t in tt]
    ax.plot([-p[1] for p in xy], [p[2] for p in xy], "-", lw = 0.4, color = "0.75",
            zorder = 1, label = r == 1 ? "posterior draws" : "")
end
for (θ, lb, sty, c) in ((THETA_LIT, "literature (Zhao et al. 2008)", "--", "C0"),
                        (θ_fit,     "fitted (UltraNest)",            "-",  "C1"))
    xy = [rd(θ, t) for t in tt]
    ax.plot([-p[1] for p in xy], [p[2] for p in xy], sty, lw = 1.8, color = c,
            label = lb, zorder = 3)
end
for i in 1:NEP
    p = rd(θ_fit, TMEAN[i])
    ax.plot([-p[1]], [p[2]], "o", ms = 7, mfc = "none", mew = 1.4, color = "C1",
            zorder = 4, label = i == 1 ? "observed epochs" : "")
end
ax.plot([0], [0], "r*", ms = 18, zorder = 5, label = "donor (B6-8 II)")
ax.invert_xaxis()
ax.set_xlabel("ΔRA East (mas)"); ax.set_ylabel("ΔDec North (mas)")
ax.legend(loc = "upper center", bbox_to_anchor = (0.5, -0.22), ncol = 3, fontsize = 9,
          frameon = false)
ax.grid(alpha = 0.3)
ax.set_title(@sprintf("β Lyrae relative orbit — %d epochs, χ²/n %.2f → %.2f",
                      NEP, sum(cl)/NTOT, sum(cf)/NTOT))
mkpath(OUT)
f = joinpath(OUT, "betlyr_orbit.png")
fig.savefig(f, dpi = 130, bbox_inches = "tight")
pyplot.close(fig)
@info "wrote $f"

# --- 4. Per-epoch panels with the components drawn to scale -------------------
# The orbit plot above shows only where the components are, not what they are — and here
# that matters: a = 0.90 mas while the donor is 0.59 mas across, so the two nearly touch.
#
# GEOMETRY, as read off `model_cvis` / `ellgauss_vis` in betlyr_model.jl:
#   * the donor is `ud_donor_mas` in DIAMETER and carries no phase factor, so it sits at the
#     origin; the orbit offset is applied to the disc.
#   * `ellgauss_vis` multiplies the rotated `up` coordinate by `ratio`, which makes the source
#     NARROWER along `pa`. So `disc_pa_deg` is the MINOR axis, the major axis is at pa+90°,
#     and `disc_fwhm_mas` is the MAJOR FWHM (verified numerically against the visibility).
#   * `up = u cos φ + v sin φ` with u↔RA(East) and v↔Dec(North), so φ runs from East toward
#     North and the astronomical PA (North through East) is 90° − φ.
patches = pyimport("matplotlib.patches")

ud, fwj, ar, pa = θ_fit[6], θ_fit[7], θ_fit[8], θ_fit[9]
pa_major_astro = mod(90 - (pa + 90), 180)
@printf("\ndisc major axis: astro PA %.2f deg   |   line of nodes Omega = %.2f deg  (delta %.2f deg)\n",
        pa_major_astro, mod(θ_fit[3], 180), abs(pa_major_astro - mod(θ_fit[3], 180)))
@printf("donor UD diameter %.3f mas, disc FWHM %.3f x %.3f mas, a = %.3f mas\n",
        ud, fwj, fwj*ar, θ_fit[1])

ncol = 4; nrow = ceil(Int, NEP/ncol)
fig2, axs = pyplot.subplots(nrow, ncol, figsize = (3.1*ncol, 3.2*nrow), squeeze = false)
axf = [axs[i-1, j-1] for i in 1:nrow for j in 1:ncol]
lim = 0.62*θ_fit[1] + fwj                        # same scale in every panel
for k in 1:NEP
    ax = axf[k]; ax.set_aspect("equal")
    xy = [rd(θ_fit, t) for t in tt]              # faint orbit for context
    ax.plot([-p[1] for p in xy], [p[2] for p in xy], "-", lw = 0.6, color = "0.8", zorder = 1)
    p = rd(θ_fit, TMEAN[k]); dx, dy = -p[1], p[2]
    # donor: uniform disc, diameter ud, at the origin
    ax.add_patch(patches.Circle((0.0, 0.0), ud/2, fc = "#d62728", ec = "k", lw = 0.8,
                                alpha = 0.85, zorder = 3))
    # disc: FWHM ellipse. matplotlib's `angle` is CCW from +x in DATA coords, and our x is
    # RA(East) and y Dec(North), which is exactly the (φ) frame — so pass pa+90 directly.
    ax.add_patch(patches.Ellipse((dx, dy), fwj, fwj*ar, angle = pa + 90,
                                 fc = "#1f77b4", ec = "k", lw = 0.8, alpha = 0.6, zorder = 2))
    ax.plot([dx], [dy], "k+", ms = 5, zorder = 4)
    ax.set_xlim(lim, -lim); ax.set_ylim(-lim, lim)     # x reversed => East left
    ax.grid(alpha = 0.25)
    ax.set_title(Dates.format(julian2datetime(TMEAN[k]), "yyyy-mm-dd"), fontsize = 10)
    k > (nrow-1)*ncol && ax.set_xlabel("ΔRA East (mas)", fontsize = 8)
    (k-1) % ncol == 0 && ax.set_ylabel("ΔDec North (mas)", fontsize = 8)
    ax.tick_params(labelsize = 7)
end
for k in (NEP+1):length(axf); axf[k].axis("off"); end
fig2.suptitle(@sprintf("β Lyrae components to scale — donor UD %.2f mas (red), disc FWHM %.2f×%.2f mas (blue), a = %.2f mas",
                       ud, fwj, fwj*ar, θ_fit[1]), fontsize = 11)
fig2.tight_layout()
f2 = joinpath(OUT, "betlyr_epochs.png")
fig2.savefig(f2, dpi = 130, bbox_inches = "tight")
pyplot.close(fig2)
@info "wrote $f2"
