# betlyr_model.jl — data, forward model and χ² for the β Lyrae orbit fit.
#
# Shared by both drivers so the model is defined once:
#   betlyr_orbit_fit_ultranest.jl — UltraNest (Python, via PythonCall), single-threaded
#   betlyr_orbit_fit_pigeons.jl   — Pigeons.jl (native Julia), multithreaded
#
# The chain is:
#     (orbital elements, component sizes) → separation at each observation
#         → component visibilities → combined complex visibility → observables → χ²
#
# ---------------------------------------------------------------------------
# Modelling choices, and why they differ from the Spica fit
# ---------------------------------------------------------------------------
#  * CIRCULAR ORBIT (e = 0 exactly). ω is then undefined and pinned; T0 alone carries the
#    phase. The residual degeneracy is (Ω, T0) → (Ω+180°, T0 + P/2), removed by restricting
#    Ω to [0°,180°) while T0 spans a full period.
#
#  * dP IS FITTED. Mass transfer lengthens the period at ≈+19 s/yr; over the 6.7 yr / 189
#    orbits spanned here the quadratic ephemeris accumulates ≈0.14 orbits ≈ 50° of phase.
#
#  * DONOR + DISC, not two stars. The donor is a Roche-filling B6-8 II star (circular
#    uniform disc here); the gainer is buried in an opaque accretion disc of outer radius
#    30 R☉ — five times the star — carrying most of the H-band flux, modelled as an
#    inclined elliptical Gaussian.
#
#  * PER-OBSERVATION TIMES. The orbit is evaluated at each observation's own MJD, solved
#    once per distinct time and broadcast back. `uv_mjd` is stored in the data eltype, so
#    the always-Float64 `v2_mjd`/`t3_mjd` are mapped through the observable→uv index arrays
#    instead — which is also what keeps absolute JD out of a narrower type (see DATA_T).

using ROTIR, Printf, Statistics, LinearAlgebra
using SpecialFunctions: besselj1

# Float64 on purpose, and the model runs Float64 throughout (see DATA_T below), so this is
# consistent. But it is also the exact trap to remember if anyone tries a Float32 pass:
# this constant sits in the hot path (`ud_vis`, `ellgauss_vis`, the phasor), so Float32 data
# multiplied by it promotes straight back to Float64. Measured on the χ² hot loop, one epoch:
#
#   properly typed Float32   primal 1.26x, ForwardDiff grad 1.28x, Zygote grad 1.05x
#   mixed (this constant left Float64)      1.02x,            1.01x,            1.00x
#
# i.e. mixed precision throws away the entire benefit — it is slower than either type alone.
# Converting Float32 would mean typing EVERY literal here (1e-10, 180/pi, 2√(2ln2), …), and
# absolute JD would still have to stay Float64. Reverse mode barely gains regardless: its
# cost is tape and closure overhead, not bandwidth (allocations fall 0.87 → 0.50 MiB for a
# 1.05x speedup).
const MAS_PER_RAD = 206264806.2
const HERE = @__DIR__
isdefined(Main, :P_ORB) || include(joinpath(HERE, "betlyr_params.jl"))

# =======================================================================================
# Data — one file per night
# =======================================================================================
# Which nights to fit, via DATASET=old|new|all (default `old`):
#
#   old  — demos/betlyr/older_data/, the 2006-2007 MIRC nights used by Zhao et al. Fitting
#          these ALONE is what makes the result directly comparable with the published
#          elements, and it converges faster because the baseline coverage is smaller.
#   new  — the 2013 MIRC-X nights sitting in demos/betlyr/.
#   all  — both, which is what `OLDDATA=1` used to mean (it APPENDED old to new rather
#          than selecting it). `OLDDATA=1` is still accepted and now maps to `old`.
const DATASET = lowercase(get(ENV, "DATASET",
                              get(ENV, "OLDDATA", "1") == "1" ? "old" : "new"))
DATASET in ("old", "new", "all") ||
    error("DATASET must be one of old|new|all (got \"$DATASET\")")

# ---------------------------------------------------------------------------------------
# Why Float64 here, when Float32 is the right default for interferometry
# ---------------------------------------------------------------------------------------
# `readoifits` defaults to `T = Float32`, and for the observables themselves that is plainly
# correct: V², T3amp and T3phi carry ~1 % errors, and the whole orbit forward model carries
# Float32 end to end at a cost of ~1e-6 mas in the predicted separation — five orders below
# MIRC's astrometric precision. Nothing in the *model* wants Float64.
#
# The scalar χ² does, and only because of how big it is. Julia's `sum` is pairwise, so the
# accumulation is essentially exact — a Float32 pairwise sum comes back bit-identical to a
# perfect Float64 sum rounded to Float32, and the relative error is flat in n (~4e-8 from
# 1e4 to 1e7 terms). The limit is REPRESENTATION of the result, not accumulation:
#
#     eps(Float32(1.96e6)) = 0.125
#
# and χ² here *is* 1.96e6 (18682 points at χ²_red ≈ 18). MCMC and nested sampling consume
# Δχ² — the 1σ contour is Δχ² = 1 — so a Float32 χ² would quantise the likelihood to ~1/8,
# and measured Δχ² errors were 0.13–0.32 against Float64. That is a real distortion of the
# acceptance ratio.
#
# Rule of thumb: a Float32 χ² is fine while `n × χ²_red` ≲ 1e5–1e6 and degrades linearly
# past that (at 1e7 residuals one ulp is 1.0 and Δχ² = 1 is unresolvable). Two things put
# this fit on the wrong side of that line: 18682 points, and a model that does not yet fit.
# Once χ²_red approaches 1 the χ² drops ~18× and Float32 would be adequate again.
#
# Loading as Float64 is the blunt fix and costs ~5 MB here. The surgical alternative — keep
# Float32 data and widen only the scalar return — is the better long-term answer and belongs
# in OITOOLS' `cvis_to_chi2_f`/`_fg`, not in this script.
const DATA_T = Float64

function load_nights()
    oifits_in(dir) = isdir(dir) ?
        filter(f -> endswith(f, ".oifits"), readdir(dir, join = true)) : String[]
    files = String[]
    DATASET in ("new", "all") && append!(files, oifits_in(HERE))
    DATASET in ("old", "all") && append!(files, oifits_in(joinpath(HERE, "older_data")))
    isempty(files) && error("no .oifits found for DATASET=$DATASET under $HERE")
    sort!(files)
    d = Any[]; nm = String[]
    for f in files
        x = try readoifits(f; T = DATA_T)[1, 1]
            catch err; @warn "skipping $(basename(f)): $err"; continue end
        isempty(x.v2_mjd) && continue
        push!(d, x); push!(nm, basename(f))
    end
    # Narrow `Any[]` to the concrete element type before returning. `OIdata{T<:AbstractFloat}`
    # binds its field types to the struct parameter, so `OIdata{Float32}` is concrete and
    # `d.v2` infers as `Vector{Float32}` — but only if the CONTAINER says so. `DATA` is a
    # const global read inside the χ², and with `Vector{Any}` (or `Vector{OIdata}`, the
    # unbound UnionAll — which is no better) `DATA[i]` gives inference nothing, and
    # model_observables / chi2_split / chi2_total all infer as `Any`.
    #
    # `map(identity, ...)` is the narrowing idiom: it rebuilds the vector at the promoted
    # concrete eltype without naming that type here, so changing DATA_T needs no edit below.
    # The try/catch loop is why it cannot simply be a typed comprehension.
    return map(identity, d), nm
end

const DATA, NIGHTS = load_nights()
const NEP   = length(DATA)
const NDATA = [d.nv2 + d.nt3amp + d.nt3phi for d in DATA]
const NTOT  = sum(NDATA)
const NV2   = sum(d.nv2 for d in DATA)
const NT3A  = sum(d.nt3amp for d in DATA)
const NT3P  = sum(d.nt3phi for d in DATA)
const TMEAN = [mean(d.v2_mjd) + 2400000.5 for d in DATA]

const UV  = [DATA_T.(d.uv) for d in DATA]     # follow the data type, do not hardcode
const UU  = [u[1, :] for u in UV]
const VV  = [u[2, :] for u in UV]
const RHO = [sqrt.(u[1, :].^2 .+ u[2, :].^2) for u in UV]

function uv_times(d)
    t = zeros(Float64, d.nuv)
    t[d.indx_v2]   .= d.v2_mjd
    t[d.indx_t3_1] .= d.t3_mjd
    t[d.indx_t3_2] .= d.t3_mjd
    t[d.indx_t3_3] .= d.t3_mjd
    return t .+ 2400000.5
end
const TUV  = [uv_times(d) for d in DATA]
const TSRT = [sort(unique(t)) for t in TUV]
const TIDX = [[searchsortedfirst(s, t) for t in tv] for (tv, s) in zip(TUV, TSRT)]

# =======================================================================================
# Component visibilities
# =======================================================================================
@inline ud_vis(diam, ρ) = (t = diam / MAS_PER_RAD * pi .* ρ .+ 1e-10;
                           2 .* besselj1.(t) ./ t)

"""
    ellgauss_vis(fwhm, ratio, pa, u, v)

Inclined elliptical Gaussian: `fwhm` the major-axis FWHM (mas), `ratio` the minor/major
axis ratio, `pa` the major-axis position angle in degrees East of North. Stands in for the
gainer's opaque accretion disc — resolved, elongated, and carrying most of the H-band flux.
"""
@inline function ellgauss_vis(fwhm, ratio, pa, u, v)
    φ = deg2rad(pa)
    up =  u .* cos(φ) .+ v .* sin(φ)
    vp = -u .* sin(φ) .+ v .* cos(φ)
    ρe = sqrt.((up .* ratio).^2 .+ vp.^2)
    σ  = fwhm / MAS_PER_RAD / (2 * sqrt(2 * log(2)))
    return exp.(-2 * (pi * σ)^2 .* ρe.^2)
end

# =======================================================================================
# Forward model
# =======================================================================================
# θ = [a, incl, Ω, T0, dP, ud_donor, disc_fwhm, disc_ratio, disc_pa, f_disc, P]
#
# P is APPENDED rather than inserted next to the other orbital elements, so that every
# existing index keeps its meaning (T0 = 4 and disc_pa = 9 are referenced by the wrapped
# parameter list in the sampler driver).
const PNAMES = ["a_mas", "incl_deg", "Omega_deg", "T0_JD", "dP_dd",
                "ud_donor_mas", "disc_fwhm_mas", "disc_ratio", "disc_pa_deg", "f_disc",
                "P_d"]
const NPAR = length(PNAMES)

"""
    model_observables(θ, i) -> (v2, t3amp, t3phi)

No caching: with multiple threads a shared mutable cache would be a data race, and the
Bessel/exp evaluations are only ~30% of the cost anyway. Allocations are local, so this is
safe to call concurrently — which is the whole point of the Pigeons driver.

Written without array mutation (comprehensions and broadcasts, never `x[k] = ...`) so the
whole χ² is reverse-mode differentiable. Zygote refuses `setindex!` on an array it is
taping, and that — not Kepler's equation — was the only thing standing between this model
and a gradient. `orbit_to_rotir_offset` is differentiable because `kepler_E` carries an
analytic `rrule` (see `src/orbits.jl`). `ras[ti]` is a gather, which Zygote handles as a
scatter-add in the pullback.
"""
function model_observables(θ, i)
    a, inc, Ω, T0, dP, ud, fw, ar, pa, f, P = θ
    bp = (i = inc, Ω = Ω, ω = OMEGA_PERI, P = P, a = a, e = E_ORB, T0 = T0,
          q = Q_BIN, dP = dP, dω = 0.0)
    ows, decs = orbit_to_rotir_offset(bp, TSRT[i])      # vectorised over epochs
    ras = -ows                                          # East = −West (mas)
    v1 = ud_vis(ud, RHO[i])
    v2 = ellgauss_vis(fw, ar, pa, UU[i], VV[i])
    uv = UV[i]; ti = TIDX[i]
    ph = cis.((-2pi / MAS_PER_RAD) .*
              (view(uv, 1, :) .* ras[ti] .+ view(uv, 2, :) .* decs[ti]))
    cvis = (v1 .+ f .* v2 .* ph) ./ (1 + f)
    return cvis_to_obs(cvis, DATA[i])
end

"Complex visibilities for epoch `i` — the part of `model_observables` before the observables."
function model_cvis(θ, i)
    a, inc, Ω, T0, dP, ud, fw, ar, pa, f, P = θ
    bp = (i = inc, Ω = Ω, ω = OMEGA_PERI, P = P, a = a, e = E_ORB, T0 = T0,
          q = Q_BIN, dP = dP, dω = 0.0)
    ows, decs = orbit_to_rotir_offset(bp, TSRT[i])
    ras = -ows
    v1 = ud_vis(ud, RHO[i])
    v2 = ellgauss_vis(fw, ar, pa, UU[i], VV[i])
    uv = UV[i]; ti = TIDX[i]
    ph = cis.((-2pi / MAS_PER_RAD) .*
              (view(uv, 1, :) .* ras[ti] .+ view(uv, 2, :) .* decs[ti]))
    return (v1 .+ f .* v2 .* ph) ./ (1 + f)
end

"Per-epoch χ² contributions, split by observable. Not differentiated — see `chi2_ad`."
@inline function chi2_epoch(θ, i)
    d = DATA[i]
    v2m, t3am, t3pm = model_observables(θ, i)
    return (sum(abs2, (v2m  .- d.v2)    ./ d.v2_err),
            sum(abs2, (t3am .- d.t3amp) ./ d.t3amp_err),
            sum(abs2, mod360(t3pm .- d.t3phi) ./ d.t3phi_err))
end

"""
    chi2_ad(θ) -> χ²

Total χ² routed through [`cvis_chi2`](@ref), whose `rrule` calls OITOOLS' hand-written
`cvis_to_chi2_fg` adjoint instead of taping the observable broadcasts. Same value as
`sum(chi2_split(θ))`; use this one when a gradient is wanted.
"""
chi2_ad(θ) = sum(cvis_chi2(model_cvis(θ, i), DATA[i]) for i in 1:NEP)

function chi2_split(θ)
    c = chi2_epoch(θ, 1)
    for i in 2:NEP
        e = chi2_epoch(θ, i)
        c = (c[1] + e[1], c[2] + e[2], c[3] + e[3])
    end
    return c
end

"Bounds check then total χ². Returns a large finite value outside the box."
function chi2_total(θ)
    (θ[1] > 0 && θ[6] > 0.02 && θ[7] > 0.02 && 0.05 < θ[8] <= 1 && 0.002 < θ[10] < 20) ||
        return 1e12
    c = chi2_split(θ)
    return isfinite(sum(c)) ? sum(c) : 1e12
end

# =======================================================================================
# Starting point and prior box
# =======================================================================================
# (Ω, ω) folded into Ω ∈ [0,180). With e = 0, ω is pinned anyway, and the compensating
# shift lives in T0 — which is sampled over a full period, so the branch is complete.
const OM_NODE_B = mod(OMEGA_NODE, 180.0)
const THETA_LIT = [A_ORB, I_ORB, OM_NODE_B, T0_ORB, DP_DD,
                   UD_DONOR, DISC_FWHM, DISC_RATIO, DISC_PA, F_DISC, P_ORB]

# Prior half-width on P, for the FITP consistency check.
#
# Scale it to the physics, not to what the data could in principle detect. beta Lyr's period
# is not constant — it drifts at Pdot = +18.93 s/yr — so `P_ORB` is the LOCAL period at the
# reference epoch. Across this 269 d baseline P itself moves by only 1.6e-4 d (13.9 s), so a
# prior wider than that is asserting an uncertainty the ephemeris flatly contradicts. The
# previous default of 0.05 d was 310x the drift.
#
# 1e-3 d is ~6x the intrinsic drift: generous against the ephemeris without being absurd.
#
# Note what this means for the check. The observable is accumulated phase drift, n*dP/P over
# n orbits; at 20.8 orbits a dP of 1e-3 d moves phase by 0.58 deg, which these data cannot
# see. So the posterior on P WILL come out flat, and that is the honest answer rather than a
# failure: the visibilities carry no information on P over this baseline, and the eclipse
# ephemeris is doing all the work. Widening PMAX to ~1e-2 would make the check discriminating
# but only by adopting a prior the ephemeris rules out.
const PMAX = parse(Float64, get(ENV, "PMAX", "1e-3"))

const DPMAX = 5 * abs(DP_DD)
# f_disc (index 10) is [0.1, 3.0], not [0.02, 10.0]. The literature value is 0.81, so the
# old prior spanned a 500x range and allowed a disc outshining the donor tenfold — not
# physical, and the width is prior volume the sampler has to compress for nothing. 4x either
# side of the published value is still generous.
const LO = [0.4,  60.0,   0.0, T0_ORB - P_ORB/2, -DPMAX, 0.05, 0.05, 0.05,   0.0, 0.10, P_ORB - PMAX]
# disc_pa upper bound is 180, NOT 360: `ellgauss_vis` squares the rotated coordinates, so
# PA and PA+180 give identical visibilities. A [0,360) prior therefore holds two exact copies
# of every solution — the posterior looks unconstrained across the full range when it is in
# fact bimodal by construction, and the sampler wastes half its volume on the duplicate.
const HI = [2.0, 120.0, 180.0, T0_ORB + P_ORB/2,  DPMAX, 1.50, 3.00, 1.00, 180.0,  3.00, P_ORB + PMAX]

function free_indices(; fitdp = get(ENV, "FITDP", "1") == "1",
                        freesize = get(ENV, "FREESIZE", "1") == "1",
                        fitp = get(ENV, "FITP", "0") == "1")
    f = trues(NPAR)
    fitdp    || (f[5] = false)
    freesize || (f[6] = f[7] = f[8] = f[9] = false)
    # P is OFF by default: the eclipse ephemeris determines it far better than 21 orbits of
    # interferometry can. FITP=1 frees it as a consistency check — recovering P_ORB says the
    # ephemeris and the visibilities agree; a significant offset would be a real result.
    fitp     || (f[11] = false)
    return findall(f)
end

function describe_model(idx = free_indices())
    betlyr_audit()
    @printf("\n%d nights, %d data points (V² %d, T3amp %d, T3φ %d)\n", NEP, NTOT, NV2, NT3A, NT3P)
    @printf("baseline %.2f yr = %.0f orbits at P = %.4f d\n",
            (maximum(TMEAN)-minimum(TMEAN))/365.25, (maximum(TMEAN)-minimum(TMEAN))/P_ORB, P_ORB)
    @printf("per-observation timing: %d uv points, %d distinct times (%.1f× fewer Kepler solves)\n",
            sum(d.nuv for d in DATA), sum(length, TSRT),
            sum(d.nuv for d in DATA) / sum(length, TSRT))
    c = chi2_split(THETA_LIT)
    @printf("\nliterature start: χ²ᵥ₂/n = %.1f  χ²ₜ₃ₐ/n = %.1f  χ²ₜ₃ₚ/n = %.1f  total χ²/n = %.1f\n",
            c[1]/NV2, c[2]/NT3A, c[3]/NT3P, sum(c)/NTOT)
    chi2_total(THETA_LIT)
    t = minimum(begin s = time_ns(); chi2_total(THETA_LIT); (time_ns()-s)/1e6 end for _ in 1:100)
    @printf("one likelihood evaluation: %.3f ms\n", t)
    @printf("free parameters (%d): %s\n", length(idx), join(PNAMES[idx], ", "))
    println("e = 0 fixed (circularised) ⇒ ω pinned; T0 alone carries the orbital phase.")
    println("Ω ∈ [0,180): (Ω,T0) → (Ω+180°, T0+P/2) is an exact degeneracy at e = 0.")
end

# =======================================================================================
# Unconstrained parameterisation, for samplers that want R^n
# =======================================================================================
# θ_k = lo_k + (hi_k − lo_k)·σ(z_k). The log-Jacobian keeps the posterior in z proper, so a
# uniform prior on the box really is uniform after the change of variables.
@inline _sigm(x) = x >= 0 ? 1 / (1 + exp(-x)) : (e = exp(x); e / (1 + e))

function z_to_theta(z, idx = free_indices())
    θ = copy(THETA_LIT)
    @inbounds for (k, p) in enumerate(idx)
        θ[p] = LO[p] + (HI[p] - LO[p]) * _sigm(z[k])
    end
    return θ
end

"""
    free_positions(idx) -> Vector{Int}

`pos[p]` is the slot of parameter `p` inside the free-parameter vector, or `0` if `p` is
held fixed. Precompute once and pass to [`z_to_theta_ad`](@ref).
"""
function free_positions(idx = free_indices())
    pos = zeros(Int, NPAR)
    for (k, p) in enumerate(idx); pos[p] = k; end
    return pos
end

"""
    z_to_theta_ad(z, pos) -> θ

Differentiable twin of [`z_to_theta`](@ref).

Two things make the plain version unusable under AD: it writes into the array with
`setindex!`, which Zygote refuses outright, and it starts from `copy(THETA_LIT)`, a
`Vector{Float64}` that cannot hold a `ForwardDiff.Dual`. This builds the vector by
comprehension instead, so the element type promotes on its own and nothing is mutated. The
index map `pos` is taken as an argument rather than recomputed, precisely so that building
*it* stays outside the differentiated code.
"""
z_to_theta_ad(z, pos) =
    [pos[p] == 0 ? convert(eltype(z), THETA_LIT[p]) :
                   LO[p] + (HI[p] - LO[p]) * _sigm(z[pos[p]]) for p in 1:NPAR]

function log_jacobian(z, idx = free_indices())
    s = 0.0
    @inbounds for (k, p) in enumerate(idx)
        u = _sigm(z[k])
        s += log(HI[p] - LO[p]) + log(u) + log1p(-u)
    end
    return s
end

"θ (physical) → z (unconstrained), for building a starting point."
function theta_to_z(θ, idx = free_indices())
    z = Vector{Float64}(undef, length(idx))
    @inbounds for (k, p) in enumerate(idx)
        u = clamp((θ[p] - LO[p]) / (HI[p] - LO[p]), 1e-6, 1 - 1e-6)
        z[k] = log(u / (1 - u))
    end
    return z
end
