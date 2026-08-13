# spica_orbit_fit_ultranest.jl
#
# Fit Spica's orbital elements DIRECTLY to the CHARA visibilities.
#
#     (orbital elements, sizes) → separation at each epoch → visibilities → observables → χ²
#
# One χ² over the whole dataset, sampled in orbital-element space. No intermediate
# per-epoch astrometry and no hand-made error ellipses: the likelihood is the OIFITS error
# bars on V², T3amp and T3φ, which is the honest error model.
#
# Why this is needed at all: with the literature elements the closure phases sit at
# χ²/n ≈ 199 while V² and T3amp are ≈ 15. T3φ is the asymmetry observable, so that is the
# binary's *position* being wrong — and every physics refinement (irradiation, Roche
# distortion, apsidal motion) has been getting judged against a broken baseline.
#
# The surface is multimodal — a scan over (ω, dω) alone runs from χ²/n = 73 to 232 with a
# degeneracy valley that changes sign — so a local optimiser would just report wherever it
# started. Hence nested sampling.
#
# The star model here is two uniform discs: analytic, microseconds per likelihood, against
# ~0.2 s for a tessellated Roche χ². Roche distortion is a few-percent correction to the
# disc sizes and will not move the orbit appreciably; once this converges, hand the orbit
# to the full model and iterate if needed.
#
# Ω is restricted to [0°,180°). Under (Ω,ω) → (Ω+180°,ω+180°) every Thiele–Innes constant
# is unchanged (each term picks up two sign flips), so the predicted separation — and
# therefore every interferometric observable — is *exactly* invariant. Sampling the full
# range would return a perfectly bimodal posterior that is an artefact. Only radial
# velocities break it, and those are deliberately not used yet.
#
# Run:  julia --project=demos demos/spica_orbit_fit_ultranest.jl
# Knobs: LIVE=400  FREEDIAM=1  FITDW=1  OUTDIR=...

ENV["MPLBACKEND"] = get(ENV, "MPLBACKEND", "Agg")
using ROTIR, PythonPlot, PythonCall, Printf, Statistics, LinearAlgebra, DelimitedFiles
using SpecialFunctions: besselj1
include(joinpath(@__DIR__, "spica_params.jl"))

nlive    = parse(Int, get(ENV, "LIVE", "400"))
freediam = get(ENV, "FREEDIAM", "1") == "1"
fitdw    = get(ENV, "FITDW", "1") == "1"
outdir   = get(ENV, "OUTDIR", joinpath(@__DIR__, "results", "spica_orbit"))
mkpath(outdir)

const MAS_PER_RAD = 206264806.2

# =======================================================================================
# Data
# =======================================================================================
spica_audit()
data_all = readoifits(joinpath(@__DIR__, "data", "2007_2012_2015.Spica.oifits"))[1, 1]
mj = sort(data_all.v2_mjd); jm = findall(diff(mj) .> 0.5)
st = mj[[1; jm .+ 1]]; sp = mj[[jm; length(mj)]]
const NEP = length(st)
data = Vector{typeof(data_all)}(undef, NEP); tepochs = zeros(NEP)
for i in 1:NEP
    idx = set_data_filter(data_all; mjd_range = [st[i] - 0.01, sp[i] + 0.01])
    data[i] = filter_data(data_all, idx)
    tepochs[i] = mean(data[i].v2_mjd) + 2400000.5
end
const NDATA = [d.nv2 + d.nt3amp + d.nt3phi for d in data]
const NTOT = sum(NDATA)
@printf("\n%d epochs, %d data points (V² %d, T3amp %d, T3φ %d)\n", NEP, NTOT,
        sum(d.nv2 for d in data), sum(d.nt3amp for d in data), sum(d.nt3phi for d in data))

# Per-epoch caches: everything that does not depend on the fitted parameters.
const UV   = [Float64.(d.uv) for d in data]
const RHO  = [sqrt.(u[1, :].^2 .+ u[2, :].^2) for u in UV]   # uv radius, cycles/rad

# --- per-OBSERVATION times, not one mean time per night ---------------------------------
# Spica's period is 4.01 d and the 2015 nights span ~4.7 h each, over which the secondary
# moves 0.24-0.52 mas and the position angle swings 8-38°. Collapsing a night to its mean
# epoch smears the binary by up to ten times the astrometric precision — and those three
# nights hold 1262 of the 2186 closure phases, so it lands squarely on the observable that
# was mysteriously bad. Each uv point therefore gets its own time.
#
# `uv_mjd` is stored in the data eltype (Float32), which quantises MJD ≈ 57131 to ~5.6 min;
# `v2_mjd` and `t3_mjd` are Float64, so the per-uv times are reconstructed from those via
# the observable→uv index maps.
#
# Only ~430 distinct times exist across all 7526 uv points (each exposure contributes many
# baselines), so the orbit is solved once per unique time and broadcast back — 17.5× fewer
# Kepler solves per likelihood.
function uv_times(d)
    t = zeros(Float64, d.nuv)
    t[d.indx_v2]   .= d.v2_mjd
    t[d.indx_t3_1] .= d.t3_mjd
    t[d.indx_t3_2] .= d.t3_mjd
    t[d.indx_t3_3] .= d.t3_mjd
    return t .+ 2400000.5                       # MJD → JD
end
const TUV  = [uv_times(d) for d in data]
const TUNI = [unique(t) for t in TUV]                      # distinct times per epoch
const TIDX = [[searchsortedfirst(sort(u), t) for t in tv]  # uv point → unique-time slot
              for (tv, u) in zip(TUV, TUNI)]
const TSRT = [sort(u) for u in TUNI]
@printf("per-observation timing: %d uv points, %d distinct times (%.1f× fewer Kepler solves)\n",
        sum(d.nuv for d in data), sum(length, TSRT),
        sum(d.nuv for d in data) / sum(length, TSRT))

# =======================================================================================
# Forward model:  (elements, sizes) → separations → visibilities → observables
# =======================================================================================
# Uniform-disc visibility, same expression as OITOOLS' `visibility_ud` but taking the
# precomputed uv radius so the per-likelihood cost is one jinc per uv point.
@inline ud_vis(diam_mas, ρ) = (t = diam_mas / MAS_PER_RAD * pi .* ρ .+ 1e-10;
                               2 .* besselj1.(t) ./ t)

"""
    model_observables(θ, i) -> (v2, t3amp, t3phi)

θ = [a, incl, Ω, ω, e, T0, dω, ud1, ud2, f]; `i` indexes the epoch.

The whole chain: the orbital elements give the instantaneous separation on the sky, that
becomes a Fourier phase shift between the two discs, and `cvis_to_obs` turns the combined
complex visibility into the observables. `f` is the secondary/primary flux ratio.
"""
# Cache of the two disc visibilities per epoch. They depend only on the diameters, so when
# those are fixed (FREEDIAM=0) this collapses ~15 000 Bessel evaluations per likelihood to
# a one-off. Safe as a mutable cache precisely because the sampler is single-threaded.
const V1C = [Vector{Float64}(undef, length(r)) for r in RHO]
const V2C = [Vector{Float64}(undef, length(r)) for r in RHO]
const UDC = [Ref((NaN, NaN)) for _ in 1:NEP]
@inline function disc_vis!(i, ud1, ud2)
    if UDC[i][] !== (ud1, ud2)
        V1C[i] .= ud_vis(ud1, RHO[i])
        V2C[i] .= ud_vis(ud2, RHO[i])
        UDC[i][] = (ud1, ud2)
    end
    return V1C[i], V2C[i]
end

function model_observables(θ, i)
    a, inc, Ω, ω, e, T0, dω, ud1, ud2, f = θ
    bp = (i = inc, Ω = Ω, ω = ω, P = P_ORB, a = a, e = e, T0 = T0, q = Q_BIN,
          dP = 0.0, dω = dω)
    # solve the orbit once per distinct observation time, then map onto the uv points.
    # The vectorised method is ~6x faster than looping the scalar one (it hoists the
    # orientation cosines out and skips the true-anomaly branch); see src/oichi2_binary.jl.
    ows, decs = orbit_to_rotir_offset(bp, TSRT[i])
    ras = -ows                                     # East = −West (mas)
    uv = UV[i]; ti = TIDX[i]
    v1, v2 = disc_vis!(i, ud1, ud2)
    cvis = Vector{ComplexF64}(undef, size(uv, 2))
    @inbounds for k in axes(uv, 2)
        j = ti[k]
        ph = cis(-2pi / MAS_PER_RAD * (uv[1, k] * ras[j] + uv[2, k] * decs[j]))
        cvis[k] = (v1[k] + f * v2[k] * ph) / (1 + f)
    end
    return cvis_to_obs(cvis, data[i])
end

"Total χ² over all epochs, split by observable."
function chi2_split(θ)
    cv = ca = cp = 0.0
    @inbounds for i in 1:NEP
        d = data[i]
        v2m, t3am, t3pm = model_observables(θ, i)
        cv += sum(abs2, (v2m .- d.v2) ./ d.v2_err)
        ca += sum(abs2, (t3am .- d.t3amp) ./ d.t3amp_err)
        cp += sum(abs2, mod360(t3pm .- d.t3phi) ./ d.t3phi_err)
    end
    return cv, ca, cp
end

# =======================================================================================
# Radial velocities (optional, USERV=1)
# =======================================================================================
# RVs add what interferometry alone cannot supply:
#   * they BREAK the (Ω,ω) → (Ω+180°,ω+180°) mirror, which is exactly degenerate for
#     astrometry, so Ω may then span the full [0,360);
#   * K1, K2 give a·sin i and q, which combined with the interferometric angular a yield
#     dynamical masses and an orbital parallax.
#
# ⚠ The error column in these files is a PLACEHOLDER — every point carries exactly
# 1.00 km/s. Taken literally the RV χ² would be ~400 per point and would bury the 6871
# interferometric points entirely. The residuals are also strongly separation-dependent:
# binned by model line separation |V2−V1| they run 36.7 km/s below 50 km/s down to 7.2 km/s
# above 300, the classic blended-line signature for a v sin i = 161 km/s SB2. So an
# empirical error model calibrated on that trend is used instead, and the blended points
# are down-weighted rather than discarded (cutting them biases K by ~7% and leaves ω
# unconstrained).
const USERV = get(ENV, "USERV", "0") == "1"
const RVD1, RVD2 = if USERV
    (readdlm(joinpath(@__DIR__, "data", "all_rv_1_ORIG.txt")),
     readdlm(joinpath(@__DIR__, "data", "all_rv_2_ORIG.txt")))
else
    (zeros(0, 3), zeros(0, 3))
end
const RVT  = Float64.(RVD1[:, 1])
const RV1  = Float64.(RVD1[:, 2])
const RV2  = Float64.(RVD2[:, 2])
const NRV  = 2 * length(RVT)
"Empirical RV error (km/s) as a function of model line separation — see the note above."
@inline rv_sigma(Δ) = 7.0 + 30.0 * exp(-abs(Δ) / 100.0)

function chi2_rv(θ, K1, K2, γ)
    NRV == 0 && return 0.0
    bp = (i = θ[2], Ω = θ[3], ω = θ[4], P = P_ORB, a = θ[1], e = θ[5], T0 = θ[6],
          q = Q_BIN, dP = 0.0, dω = θ[7])
    m1, _ = binary_RV(bp, RVT; K1 = K1, K2 = K2, γ = γ)
    _, m2 = binary_RV(bp, RVT; K1 = K1, K2 = K2, γ = γ)
    c = 0.0
    @inbounds for k in eachindex(RVT)
        σ = rv_sigma(m2[k] - m1[k])
        c += ((RV1[k] - m1[k])^2 + (RV2[k] - m2[k])^2) / σ^2
    end
    return c
end

function chi2_total(θ)
    (θ[1] > 0 && 0 <= θ[5] < 0.6 && θ[8] > 0.05 && θ[9] > 0.02 && 0.005 < θ[10] < 3) ||
        return 1e12
    c = sum(chi2_split(θ))
    USERV && (c += chi2_rv(θ, θ[11], θ[12], θ[13]))
    return isfinite(c) ? c : 1e12
end

const UD1_0 = 2RPOLE1
const UD2_0 = 2RPOLE2
const F0    = (RPOLE2 / RPOLE1)^2 * (TPOLE2 / TPOLE1)
# Map the literature (Ω, ω) into the Ω ∈ [0,180) branch. The shift must be applied to the
# PAIR: (Ω, ω) → (Ω−180°, ω−180°) is the exact invariance, whereas moving Ω alone gives a
# genuinely different — and much worse — orbit.
const OM_NODE_B, OM_PERI_B = OMEGA_NODE >= 180 ?
    (OMEGA_NODE - 180, mod(OMEGA_PERI - 180, 360)) : (OMEGA_NODE, OMEGA_PERI)
# The (Ω,ω) → (Ω±180°,ω±180°) fold is invariant for ASTROMETRY but NOT for radial
# velocities: an RV depends on ω alone, so shifting ω by 180° inverts the curve. That is
# precisely the degeneracy RVs break — so when they are in play the literature branch is
# used unfolded, and Ω is freed to the full circle instead.
const OM_NODE_S, OM_PERI_S = USERV ? (OMEGA_NODE, OMEGA_PERI) : (OM_NODE_B, OM_PERI_B)
θ_lit = [A_ORB, I_ORB, OM_NODE_S, OM_PERI_S, E_ORB, T0_ORB, DOMEGA, UD1_0, UD2_0, F0,
         123.9, 198.8, 0.0]      # K1, K2, γ — literature semi-amplitudes (km/s)
@printf("literature (Ω, ω) = (%.3f, %.3f) → branch-mapped (%.3f, %.3f)\n",
        OMEGA_NODE, OMEGA_PERI, OM_NODE_B, OM_PERI_B)

println("\n" * "="^78)
println("Starting point: literature elements")
println("="^78)
let c = chi2_split(θ_lit)
    @printf("  χ²ᵥ₂/n = %.1f   χ²ₜ₃ₐ/n = %.1f   χ²ₜ₃ₚ/n = %.1f   total χ²/n = %.1f\n",
            c[1]/sum(d.nv2 for d in data), c[2]/sum(d.nt3amp for d in data),
            c[3]/sum(d.nt3phi for d in data), sum(c)/NTOT)
end
@printf("  one likelihood evaluation: ")
let
    chi2_total(θ_lit)                                   # warm up: exclude JIT
    t = minimum(begin s = time_ns(); chi2_total(θ_lit); (time_ns()-s)/1e6 end for _ in 1:200)
    @printf("%.3f ms\n", t)
end

# BENCH=1 reports where the likelihood time actually goes. The two disc visibilities
# depend only on the diameters, so when those are fixed the cache turns ~15 000 Bessel
# evaluations per call into a one-off; when they are free it is paid every call.
if get(ENV, "BENCH", "0") == "1"
    fixed()   = chi2_total(θ_lit)
    varying() = (t = copy(θ_lit); t[8] += 1e-6rand(); t[9] += 1e-6rand(); chi2_total(t))
    bench(f, n=200) = (f(); minimum(begin s=time_ns(); f(); (time_ns()-s)/1e6 end for _=1:n))
    tf = bench(fixed); tv = bench(varying)
    @printf("\nlikelihood, diameters FIXED   (Bessel cached)   %.3f ms\n", tf)
    @printf("likelihood, diameters VARYING (Bessel each call) %.3f ms\n", tv)
    @printf("=> Bessel recomputation = %.3f ms, %.0f%% of the likelihood\n", tv-tf,
            100*(tv-tf)/max(tv, eps()))
end

# =======================================================================================
# Nested sampling over the elements
# =======================================================================================
println("\n" * "="^78)
println("Direct fit of the orbital elements to the visibilities")
println("="^78)
names = ["a_mas", "incl_deg", "Omega_deg", "omega_deg", "ecc", "T0_JD",
         "domega_degday", "ud1_mas", "ud2_mas", "fratio", "K1_kms", "K2_kms", "gamma_kms"]
lo = [1.0,   90.0,   0.0,   0.0, 0.0,  T0_ORB - P_ORB/2, -3*360/(80*365.25), 0.5, 0.2, 0.05,  80.0, 140.0, -40.0]
hi = [2.2,  150.0, 180.0, 360.0, 0.35, T0_ORB + P_ORB/2,  3*360/(80*365.25), 1.4, 0.9, 0.60, 170.0, 260.0,  40.0]
free = trues(13)
# Without RVs, K1/K2/γ are unconstrained; with them, the (Ω,ω) mirror is broken so Ω is
# freed to the full circle.
if !USERV
    free[11] = free[12] = free[13] = false
else
    hi[3] = 360.0
end
if !fitdw;    free[7] = false; lo[7] = hi[7] = 0.0;        end
if !freediam; free[8] = free[9] = false
              lo[8] = hi[8] = UD1_0; lo[9] = hi[9] = UD2_0 end
idx = findall(free)
@printf("free parameters (%d): %s\n", length(idx), join(names[idx], ", "))
if USERV
    @printf("RVs: %d points (%d per component), empirical σ(Δ) = 7 + 30·exp(−|Δ|/100) km/s\n",
            NRV, length(RVT))
    println("     ⇒ the (Ω,ω) mirror is broken, so Ω spans the full [0,360).")
    @printf("     χ²_RV at the literature start = %.1f (%.1f per point)\n",
            chi2_rv(θ_lit, θ_lit[11], θ_lit[12], θ_lit[13]),
            chi2_rv(θ_lit, θ_lit[11], θ_lit[12], θ_lit[13])/NRV)
else
    println("RVs: NOT used (USERV=1 to include them). Ω restricted to [0,180) because")
    println("     (Ω,ω) → (Ω+180°,ω+180°) is an exact astrometric degeneracy.")
end

# DRYRUN=1 stops here — useful for timing the likelihood or checking a starting point
# without paying for a full sampling run.
if get(ENV, "DRYRUN", "0") == "1"
    @info "DRYRUN=1: stopping before sampling"
    exit(0)
end

println("Ω ∈ [0,180): (Ω,ω) → (Ω+180°,ω+180°) leaves every observable exactly unchanged.")

# UltraNest works on the unit hypercube. Plain (non-vectorized) Julia closures, matching
# the pattern in src/ultranest.jl: `collect(Float64, cube)` is needed because PyCall hands
# the cube over as a PyObject. At ~1 ms per likelihood the per-call overhead is nothing.
# Vectorized UltraNest, following OITOOLS' fit_model_ultranest idiom. The essential
# detail is the DECLARED argument types: PyCall then converts the numpy batch into a Julia
# matrix. Passing the callables through `pyfunction(..., PyObject)` instead suppresses that
# conversion and the transform receives a Vector{Any} of rows, which fails to broadcast.
const Δ = hi[idx] .- lo[idx]
prior_1(u::AbstractVector{<:Real}) = u .* Δ .+ lo[idx]
prior_v(U::AbstractMatrix{<:Real}) = reduce(vcat, (u -> prior_1(u)').(eachrow(U)))

function loglike_1(p::AbstractVector{<:Real})
    θ = copy(θ_lit); θ[idx] .= p
    v = -0.5 * chi2_total(θ)
    return isfinite(v) ? v : -1e30
end
# NOT threaded, deliberately. PyCall is not thread-safe: `PyObject` finalizers call
# Py_DECREF without holding the GIL, and Julia's GC will happily run them on a worker
# thread, which segfaults inside `pydecref_`. It does not matter that `loglike_1` never
# touches Python — merely having Julia threads live while PyObjects exist is enough.
#
# So with an UltraNest (PyCall) sampler the likelihood is single-threaded regardless of
# `-t auto`, and the batching buys only the ~30 μs Python↔Julia crossing per call. A pure
# Julia sampler (Pigeons, see demos/rho_cas_pigeons.jl) has no such restriction and can use
# `multithreaded = true` safely — which is the real argument for switching.
loglike_v(X::AbstractMatrix{<:Real}) = loglike_1.(eachrow(X))

ultranest = pyimport("ultranest")
sampler = ultranest.ReactiveNestedSampler(names[idx], loglike_v;
                                          transform = prior_v, vectorized = true)
res = sampler.run(min_num_live_points = nlive, show_status = true)
sampler.print_results()

post = Array{Float64}(res["samples"])
θ_fit = copy(θ_lit); θ_fit[idx] .= [median(post[:, k]) for k in 1:length(idx)]

println("\n" * "="^78); println("RESULT"); println("="^78)
@printf("%-16s %14s %14s %22s\n", "parameter", "literature", "fitted", "68% interval")
for (k, p) in enumerate(idx)
    @printf("%-16s %14.5f %14.5f  [%.5f, %.5f]\n", names[p], θ_lit[p], θ_fit[p],
            quantile(post[:, k], 0.16), quantile(post[:, k], 0.84))
end
for (lab, θ) in (("literature", θ_lit), ("fitted", θ_fit))
    c = chi2_split(θ)
    @printf("\n%-11s  χ²ᵥ₂/n = %7.2f  χ²ₜ₃ₐ/n = %7.2f  χ²ₜ₃ₚ/n = %7.2f  total = %7.2f",
            lab, c[1]/sum(d.nv2 for d in data), c[2]/sum(d.nt3amp for d in data),
            c[3]/sum(d.nt3phi for d in data), sum(c)/NTOT)
end
println()
if fitdw && abs(θ_fit[7]) > 1e-7
    @printf("\napsidal period U = %.0f yr   (Robinette & Aufdenberg 2015: 139 ± 6)\n",
            360 / (abs(θ_fit[7]) * 365.25))
end
@printf("physical radii at d = %.0f pc: R1 = %.2f R☉, R2 = %.2f R☉  (lit. 7.40 / 3.74)\n",
        D_PC, phys_radius_rsun(θ_fit[8]/2, D_PC), phys_radius_rsun(θ_fit[9]/2, D_PC))
println("\nNOTE: astrometry cannot distinguish (Ω,ω) from (Ω+180°,ω+180°); RVs would.")

# ---- per-epoch breakdown and figure ----------------------------------------------------
println("\nper-epoch reduced χ² :")
@printf("%3s %11s %8s %9s %9s %9s\n", "ep", "JD-2400000", "ρ(mas)", "lit", "fitted", "n")
for i in 1:NEP
    d = data[i]
    c(θ) = begin
        v2m, t3am, t3pm = model_observables(θ, i)
        (sum(abs2, (v2m .- d.v2)./d.v2_err) + sum(abs2, (t3am .- d.t3amp)./d.t3amp_err) +
         sum(abs2, mod360(t3pm .- d.t3phi)./d.t3phi_err)) / NDATA[i]
    end
    ow, on = orbit_to_rotir_offset((i = θ_fit[2], Ω = θ_fit[3], ω = θ_fit[4], P = P_ORB,
                                    a = θ_fit[1], e = θ_fit[5], T0 = θ_fit[6], q = Q_BIN,
                                    dP = 0.0, dω = θ_fit[7]), tepochs[i])
    @printf("%3d %11.4f %8.4f %9.2f %9.2f %9d\n", i, tepochs[i]-2400000,
            hypot(ow, on), c(θ_lit), c(θ_fit), NDATA[i])
end

fig, ax = subplots(figsize = (8, 8)); ax.set_aspect("equal")
tt = T0_ORB .+ range(0, P_ORB, length = 500)
rd(θ, t) = (ow, on) = orbit_to_rotir_offset((i = θ[2], Ω = θ[3], ω = θ[4], P = P_ORB,
                a = θ[1], e = θ[5], T0 = θ[6], q = Q_BIN, dP = 0.0, dω = θ[7]), t)
for (θ, lab, sty) in ((θ_lit, "literature", "--"), (θ_fit, "fitted", "-"))
    xy = [rd(θ, t) for t in tt]
    ax.plot([-p[1] for p in xy], [p[2] for p in xy], sty, lw = 1.5, label = lab)
end
for i in 1:NEP
    for (θ, m) in ((θ_lit, "s"), (θ_fit, "o"))
        p = rd(θ, tepochs[i]); ax.plot([-p[1]], [p[2]], m, ms = 6, mfc = "none",
                                       color = θ === θ_lit ? "C0" : "C1")
    end
end
ax.plot([0], [0], "r*", ms = 18, label = "primary")
ax.invert_xaxis(); ax.set_xlabel("ΔRA East (mas)"); ax.set_ylabel("ΔDec North (mas)")
ax.legend(); ax.grid(alpha = 0.3); ax.set_title("Spica relative orbit — direct fit to visibilities")
fig.savefig(joinpath(outdir, "orbit_fit.png"), dpi = 130, bbox_inches = "tight")
PyPlot.close(fig)
writedlm(joinpath(outdir, "orbit_posterior.txt"), post)
writedlm(joinpath(outdir, "orbit_best.txt"), hcat(names, θ_fit))
@info "results in $outdir"
