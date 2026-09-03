# beta Lyrae — orbit and component figures from a STORED `fit_orbit` result
#
#   julia --project=demos demos/betlyr/betlyr_orbit_fit_figure.jl
#   FIT=free julia --project=demos demos/betlyr/betlyr_orbit_fit_figure.jl   # the TIED=0 run
#
# Reads `results/orbitfit_<tag>_best.txt` (+ posterior) written by
# `betlyr_orbit_fit_tied.jl` and redraws the sky-projected relative orbit and the components
# to scale, so figures can be regenerated in seconds without re-running the sampler.
#
# This is the counterpart of `betlyr_orbit_figure.jl`, which reads the BESPOKE
# `betlyr_model.jl` fit. The two parameter layouts genuinely differ — 11 hand-rolled
# parameters there against `ORBIT_ELEMENTS` + component parameters + `f` here — so they are
# separate scripts rather than one with a branch.

get!(ENV, "MPLBACKEND", "Agg")

using ROTIR, Printf, Statistics, DelimitedFiles, PythonPlot, PythonCall, Dates

const HERE = @__DIR__
# Matches the tag `betlyr_orbit_fit_tied.jl` writes: "<dataset>_tied[_fitP][_fitdP]".
const TAG  = get(ENV, "FIT", "all_tied")
const RES  = joinpath(HERE, "results")
const OUT  = get(ENV, "OUT", joinpath(HERE, "..", "..", "docs", "src", "assets"))
const BEST = joinpath(RES, "orbitfit_$(TAG)_best.txt")
isfile(BEST) || error("no stored fit at $BEST — run betlyr_orbit_fit_tied.jl first")

best = readdlm(BEST)
V    = Dict(Symbol(best[k, 1]) => Float64(best[k, 2]) for k in 1:size(best, 1))
post = let f = joinpath(RES, "orbitfit_$(TAG)_posterior.txt")
    isfile(f) ? readdlm(f) : zeros(0, 0)
end

# Load whichever nights the fit used — the tag's leading field says which. Reading only
# older_data/ would silently plot 7 epochs of an 18-epoch fit.
const DATASET = split(TAG, "_")[1]
DATASET in ("all", "2007", "2013") ||
    error("cannot infer the dataset from FIT=\"$TAG\"; expected it to start all_/2007_/2013_")
_oifits(dir) = sort!(filter(f -> endswith(f, ".oifits"), readdir(dir, join = true)))
files = String[]
DATASET in ("all", "2007") && append!(files, _oifits(joinpath(HERE, "older_data")))
DATASET in ("all", "2013") && append!(files, _oifits(HERE))
data = [readoifits(f; T = Float64, verbose = false)[1, 1] for f in files]
data = [d for d in data if !isempty(d.v2_mjd)]
tmean = [mean(d.v2_mjd) + 2400000.5 for d in data]
ntot  = sum(d.nv2 + d.nt3amp + d.nt3phi for d in data)
@printf("%s fit: %d epochs, %d points\n", uppercase(TAG), length(data), ntot)

# Literature, for the overlay. Ω is quoted as 253.7 deg; the projected orbit is unchanged by
# a half turn, so it is folded into [0,180) to sit on the same branch as the fit.
LIT = (a = 0.865, i = 93.5, Omega = mod(253.7, 180), P = 12.9414, T0 = 2454283.0430,
       ud = 0.57, fwhm = 0.90, ratio = 0.63/1.04)

bp(a, i, Om, P, T0) = (i = i, Ω = Om, ω = 90.0, P = P, a = a, e = 0.0, T0 = T0,
                       q = 1.0, dP = 0.0, dω = 0.0)
rd(el, t) = orbit_to_rotir_offset(bp(el...), t)

fit_el = (V[:a], V[:i], V[:Omega], V[:P], V[:T0])
lit_el = (LIT.a, LIT.i, LIT.Omega, LIT.P, LIT.T0)
tt = V[:T0] .+ range(0, V[:P], length = 400)

# --- 1. Relative orbit --------------------------------------------------------
fig, ax = subplots(figsize = (9, 4.2)); ax.set_aspect("equal")
if !isempty(post)
    # Posterior draws sit UNDER the curves. Columns follow spec.names[spec.free], which
    # `orbit_fit_spec` SORTS by parameter index — not the order `free` was written in the
    # fitting script. P and dP are only present when the tag says they were freed, and they
    # sort into the MIDDLE (indices 6 and 8), so a hardcoded list silently mislabels every
    # column after them.
    cols = Symbol[:a, :i, :Omega]
    occursin("_fitP", TAG)  && push!(cols, :P)
    push!(cols, :T0)
    occursin("_fitdP", TAG) && push!(cols, :dP)
    append!(cols, [:c1_diameter, :c2_fwhm, :c2_ratio])
    occursin("_free", TAG)  && push!(cols, :c2_pa)     # untied runs fit pa as well
    push!(cols, :f)
    size(post, 2) == length(cols) ||
        error("posterior has $(size(post,2)) columns but the tag \"$TAG\" implies " *
              "$(length(cols)): $(join(string.(cols), ", "))")
    ia, ii, iO, iT = indexin([:a, :i, :Omega, :T0], cols)
    for r in round.(Int, range(1, size(post, 1), length = min(60, size(post, 1))))
        xy = [rd((post[r,ia], post[r,ii], post[r,iO], V[:P], post[r,iT]), t) for t in tt]
        ax.plot([-p[1] for p in xy], [p[2] for p in xy], "-", lw = 0.4, color = "0.75",
                zorder = 1, label = r == 1 ? "posterior draws" : "")
    end
end
for (el, lb, sty, c) in ((lit_el, "literature (Z08 / M18)", "--", "C0"),
                         (fit_el, "fitted ($TAG)",          "-",  "C1"))
    xy = [rd(el, t) for t in tt]
    ax.plot([-p[1] for p in xy], [p[2] for p in xy], sty, lw = 1.8, color = c,
            label = lb, zorder = 3)
end
for k in eachindex(tmean)
    p = rd(fit_el, tmean[k])
    ax.plot([-p[1]], [p[2]], "o", ms = 7, mfc = "none", mew = 1.4, color = "C1",
            zorder = 4, label = k == 1 ? "observed epochs" : "")
end
ax.plot([0], [0], "r*", ms = 18, zorder = 5, label = "donor (B6-8 II)")
ax.invert_xaxis()
ax.set_xlabel("ΔRA East (mas)"); ax.set_ylabel("ΔDec North (mas)")
ax.legend(loc = "upper center", bbox_to_anchor = (0.5, -0.22), ncol = 3, fontsize = 9,
          frameon = false)
ax.grid(alpha = 0.3)
ax.set_title(@sprintf("β Lyrae relative orbit — %s fit, %d epochs", TAG, length(data)))
mkpath(OUT)
f1 = joinpath(OUT, "betlyr_orbitfit_$(TAG).png")
fig.savefig(f1, dpi = 130, bbox_inches = "tight"); pyplot.close(fig)
@info "wrote $f1"

# --- 2. Components to scale, per epoch ----------------------------------------
# `EllipticalGaussian` scales the rotated `up` axis by `ratio < 1`, so `pa` is the MINOR
# axis and the major axis lies at pa+90. `up = u cosφ + v sinφ` with u↔RA(East), so φ runs
# East toward North and matplotlib's CCW-from-+x angle is exactly pa+90 in these axes.
patches = pyimport("matplotlib.patches")
ud, fw, ar, pa = V[:c1_diameter], V[:c2_fwhm], V[:c2_ratio], V[:c2_pa]
@printf("donor UD %.3f mas, disc FWHM %.3f x %.3f mas, a = %.3f mas\n", ud, fw, fw*ar, V[:a])
@printf("disc major axis astro PA %.2f deg vs Omega %.2f deg\n",
        mod(-pa, 180), mod(V[:Omega], 180))

ncol = length(data) > 8 ? 6 : 4; nrow = ceil(Int, length(data)/ncol)
fig2, axs = pyplot.subplots(nrow, ncol, figsize = (2.6*ncol, 2.8*nrow), squeeze = false)
axf = [axs[i-1, j-1] for i in 1:nrow for j in 1:ncol]
lim = 0.62*V[:a] + fw
for k in eachindex(tmean)
    ax = axf[k]; ax.set_aspect("equal")
    xy = [rd(fit_el, t) for t in tt]
    ax.plot([-p[1] for p in xy], [p[2] for p in xy], "-", lw = 0.6, color = "0.8", zorder = 1)
    p = rd(fit_el, tmean[k]); dx, dy = -p[1], p[2]
    ax.add_patch(patches.Circle((0.0, 0.0), ud/2, fc = "#d62728", ec = "k", lw = 0.8,
                                alpha = 0.85, zorder = 3))
    ax.add_patch(patches.Ellipse((dx, dy), fw, fw*ar, angle = pa + 90,
                                 fc = "#1f77b4", ec = "k", lw = 0.8, alpha = 0.6, zorder = 2))
    ax.plot([dx], [dy], "k+", ms = 5, zorder = 4)
    ax.set_xlim(lim, -lim); ax.set_ylim(-lim, lim)
    ax.grid(alpha = 0.25)
    ax.set_title(Dates.format(julian2datetime(tmean[k]), "yyyy-mm-dd"), fontsize = 10)
    k > (nrow-1)*ncol && ax.set_xlabel("ΔRA East (mas)", fontsize = 8)
    (k-1) % ncol == 0 && ax.set_ylabel("ΔDec North (mas)", fontsize = 8)
    ax.tick_params(labelsize = 7)
end
for k in (length(data)+1):length(axf); axf[k].axis("off"); end
fig2.suptitle(@sprintf("β Lyrae components to scale (%s fit) — donor UD %.2f mas, disc FWHM %.2f×%.2f mas, a = %.2f mas",
                       TAG, ud, fw, fw*ar, V[:a]), fontsize = 11)
fig2.tight_layout()
f2 = joinpath(OUT, "betlyr_orbitfit_$(TAG)_epochs.png")
fig2.savefig(f2, dpi = 130, bbox_inches = "tight"); pyplot.close(fig2)
@info "wrote $f2"
