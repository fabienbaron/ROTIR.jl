# Shared ρ Cas model definition for the sampler comparison.
#
# Included by rho_cas_{pigeons,ultranest,pathfinder,bootstrap_mcmc}.jl so that every
# method targets the *same* posterior: identical likelihood, identical uniform-box prior,
# identical parameterisation. Error bars from different samplers are only comparable if
# the posterior is the same one.
#
# Environment knobs (all optional):
#   FREE=rpole,ld1     which parameters to sample (default), or FREE=all7
#   NSIDE=3            healpix resolution exponent (3 → 768 tessels)
#   DATAFILE=...       OIFITS path

using ROTIR, Printf

const DATAFILE = get(ENV, "DATAFILE", "./demos/data/rho_Cas_example.oifits")
const NSIDE    = parse(Int, get(ENV, "NSIDE", "3"))

const DATA    = readoifits_multiepochs([DATAFILE]; T = Float64)[1, :]
const TEPOCHS = [0.0]
const TESSELS = tessellation_healpix(NSIDE, T = Float64)

# ω is a free parameter only in the all-7 case; frozen at 0 the rapid-rotator surface is
# exactly a limb-darkened sphere.
const BASE = (surface_type = 2, rpole = 1.1, tpole = 4000.0,
              ldtype = 3, ld1 = 0.75, ld2 = 0.0,
              inclination = 0.0, position_angle = 0.0, rotation_period = 1e6,
              beta = 0.08, frac_escapevel = 0.0, B_rot = 0.0)

# Starting point: the better of the basins found by demos/rho_cas_basins.jl
const THETA0 = [1.105, 0.0, 0.0, 0.0, 0.08, 0.745, 0.0]

const FREE_NAMES = let s = get(ENV, "FREE", "rpole,ld1")
    s == "all7" ? parametric_param_names() : String.(split(s, ","))
end
const IFREE  = parametric_free_indices(FREE_NAMES)
const LABELS = parametric_param_names()[IFREE]

# ── Prior box (finite: nested sampling needs it, and it keeps every method identical) ──
const BOX_LO = [0.5, 0.0,   0.0, -180.0, 0.0, 0.0, -1.0][IFREE]
const BOX_HI = [4.0, 0.99, 180.0, 180.0, 1.0, 2.0,  1.0][IFREE]

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

function describe_model()
    @printf("ρ Cas: %d points, %d tessels, sampling %d parameter(s): %s\n",
            sum(d -> d.nv2 + d.nt3amp + d.nt3phi, DATA), TESSELS.npix, NDIM,
            join(LABELS, ", "))
    @printf("  box prior: %s\n", join([@sprintf("%s ∈ [%.3g, %.3g]", LABELS[j], BOX_LO[j],
                                                BOX_HI[j]) for j in 1:NDIM], "  "))
end
