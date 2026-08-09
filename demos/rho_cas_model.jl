# Shared ρ Cas model definition for the sampler comparison.
#
# Included by rho_cas_{pigeons,ultranest,pathfinder,bootstrap_mcmc}.jl so that every
# method targets the *same* posterior: identical likelihood, identical uniform-box prior,
# identical parameterisation. Error bars from different samplers are only comparable if
# the posterior is the same one.
#
# Environment knobs (all optional):
#   FREE=...           which parameters to sample. Default is the six identifiable ones,
#                      rpole,omega,inc,PA,beta,ld1. FREE=rpole,ld1 gives the sphere + LD
#                      model (ω frozen at 0). Avoid ld2: ldtype=3 ignores it entirely.
#   NSIDE=3            healpix resolution exponent (3 → 768 tessels)
#   DATAFILE=...       OIFITS path

using ROTIR, Printf

const DATAFILE = get(ENV, "DATAFILE", "./demos/data/rho_Cas_example.oifits")
const NSIDE    = parse(Int, get(ENV, "NSIDE", "3"))

const DATA    = readoifits_multiepochs([DATAFILE]; T = Float64)[1, :]
const TEPOCHS = [0.0]
const TESSELS = tessellation_healpix(NSIDE, T = Float64)

# Values here are only used for parameters that are NOT free (see FREE below). With the
# default six free, only tpole, ld2, rotation_period and B_rot come from this.
# Note: with ω frozen at 0 instead, the rapid-rotator surface is exactly a limb-darkened
# sphere and inc/PA/beta become unobservable — check_identifiability() warns about that.
const BASE = (surface_type = 2, rpole = 1.1, tpole = 4000.0,
              ldtype = 3, ld1 = 0.75, ld2 = 0.0,
              inclination = 25.4, position_angle = 94.6, rotation_period = 1e6,
              beta = 0.08, frac_escapevel = 0.0, B_rot = 0.0)

# Starting point: the oblate + gravity-darkened optimum (χ²ᵣ = 2.82), which is also a
# valid interior start for the sphere case. Every free parameter must start strictly
# inside the box below — the logit transform sends a value sitting on a bound to ±Inf.
const THETA0 = [0.877, 0.940, 25.4, 94.6, 0.086, 0.696, 0.0]

const FREE_NAMES = let s = get(ENV, "FREE", "rpole,omega,inc,PA,beta,ld1")
    s == "all7" ? parametric_param_names() : String.(split(s, ","))
end
const IFREE  = parametric_free_indices(FREE_NAMES)
const LABELS = parametric_param_names()[IFREE]

# ── Prior box (finite: nested sampling needs it, and it keeps every method identical) ──
# inc is capped at 90° and PA at 180° deliberately. An oblate star with a symmetric von
# Zeipel map looks identical under inc → 180−inc and under PA → PA+180, so the full ranges
# would contain four exact copies of every mode. Sampling them costs effort and makes the
# posterior look multimodal for a reason that is pure labelling.
const BOX_LO_FULL = [0.5, 0.0,  0.0,   0.0, 0.0, 0.0, -1.0]
const BOX_HI_FULL = [4.0, 0.99, 90.0, 180.0, 1.0, 2.0,  1.0]
const BOX_LO = BOX_LO_FULL[IFREE]
const BOX_HI = BOX_HI_FULL[IFREE]

const LOGPI = build_parametric_logπ(DATA, TESSELS, TEPOCHS, BASE)

# Full θ from the free subvector, built as θ_frozen + S·v rather than by mutating a copy:
# Zygote cannot differentiate `θ[IFREE] .= v`, and every gradient-based sampler here
# differentiates straight through this function.
const THETA_FROZEN = let θ = copy(THETA0); θ[IFREE] .= 0.0; θ end
const SCATTER = let S = zeros(length(THETA0), length(IFREE))
    for (j, i) in enumerate(IFREE); S[i, j] = 1.0; end
    S
end

"Full θ vector with the free entries replaced (AD-safe)."
theta_full(v) = THETA_FROZEN .+ SCATTER * v

"Likelihood only, on physical parameters (−χ²/2)."
loglike_theta(v) = LOGPI(theta_full(v))

# ── Unconstrained space for the gradient-based samplers ─────────────────────
# θ = lo + (hi−lo)·σ(z) maps ℝⁿ onto the box, so target = loglike + log|dθ/dz| in z-space
# is *exactly* a uniform prior over the box — the same prior UltraNest gets from its unit
# hypercube transform. No sampler can propose outside the box, and rpole can never reach 0
# (where the flux normalisation makes the likelihood NaN).
_logistic(x) = 1 / (1 + exp(-x))

z_to_theta(z) = BOX_LO .+ (BOX_HI .- BOX_LO) .* _logistic.(z)
theta_to_z(v) = log.((v .- BOX_LO) ./ (BOX_HI .- v))
logabsjac(z)  = sum(log.(BOX_HI .- BOX_LO) .+ log.(_logistic.(z)) .+ log.(1 .- _logistic.(z)))

function logpost_z(z)
    any(!isfinite, z) && return -Inf
    v = loglike_theta(z_to_theta(z))
    isfinite(v) || return -Inf
    return v + logabsjac(z)
end

const NDIM = length(IFREE)
const Z0   = theta_to_z(THETA0[IFREE])

"Dispersed starting points in z, for multi-start methods."
function dispersed_starts(n; rng)
    [theta_to_z(BOX_LO .+ (BOX_HI .- BOX_LO) .* (0.05 .+ 0.9 .* rand(rng, NDIM)))
     for _ in 1:n]
end

"""
Warn about free parameters the likelihood cannot see. A flat direction is not merely
wasteful: it returns its prior as if it were a measurement, and for nested sampling it
costs a whole dimension.
"""
function check_identifiability()
    inert = String[]
    # ldtype 3 (Hestroffer, I ∝ μ^ld1) ignores ld2 entirely — see compute_ldmap.
    BASE.ldtype == 3 && "ld2" in LABELS && push!(inert, "ld2 (unused by ldtype=3)")
    BASE.ldtype == 1 && "ld2" in LABELS && push!(inert, "ld2 (unused by ldtype=1)")
    # With ω frozen at 0 the surface is a sphere with a rotationally symmetric intensity
    # map, so orientation is unobservable.
    if !("omega" in LABELS) && THETA0[2] == 0
        for p in ("inc", "PA")
            p in LABELS && push!(inert, "$p (sphere: ω frozen at 0)")
        end
        "beta" in LABELS && push!(inert, "beta (sphere: uniform von Zeipel map)")
    end
    # tpole cancels in the flux normalisation under the linear intensity proxy.
    "tpole" in LABELS && push!(inert, "tpole (degenerate unless intensity_model=:planck)")

    isempty(inert) && return nothing
    @warn """Free parameters the likelihood cannot constrain: $(join(inert, ", ")).
             Their marginals will reproduce the prior, and every extra dimension costs
             sampling effort — nested sampling most of all. Consider
             FREE=rpole,omega,inc,PA,beta,ld1"""
    return nothing
end

function describe_model()
    @printf("ρ Cas: %d points, %d tessels, sampling %d parameter(s): %s\n",
            sum(d -> d.nv2 + d.nt3amp + d.nt3phi, DATA), TESSELS.npix, NDIM,
            join(LABELS, ", "))
    @printf("  box prior: %s\n", join([@sprintf("%s ∈ [%.3g, %.3g]", LABELS[j], BOX_LO[j],
                                                BOX_HI[j]) for j in 1:NDIM], "  "))
    check_identifiability()
end
