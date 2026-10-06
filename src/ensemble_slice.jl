# ensemble_slice.jl — Ensemble Slice Sampling (Karamanis & Beutler 2020, the `zeus` algorithm).
#
# A GRADIENT-FREE sampler, complementing the gradient-based backends (`fit_hmc`, Pigeons,
# Nautilus, UltraNest). Validated against zeus 2.5.4 on a correlated Gaussian, Neal's funnel and
# Rosenbrock: at evaluation parity (4.91 / 5.85 / 6.24 log-density calls per step against zeus's
# 4.88 / 5.86 / 6.67) and recovering the funnel's v marginal to 2.96 against a true 3.00.
#
# WHERE IT APPLIES. Chen (2025, arXiv:2505.02987) shows derivative-free affine-invariant ensembles
# need **N = 2d walkers** with O(d) autocorrelation. That makes this the right tool for the
# GEOMETRY fits — 5-6 parameters, so a dozen walkers suffice — and the wrong tool for imaging: a
# 768-tessel surface would need 1536 walkers and tens of hours per independent sample.
#
# WHAT SLICE SAMPLING BUYS, AND WHAT IT DOES NOT. There is no accept/reject: the shrinking
# procedure always returns a point inside the slice, so the acceptance collapse that cripples
# Metropolis proposals in high dimensions does not arise, and the bracket adapts to the target's
# local geometry along the chosen direction. It does NOT remove direction confinement — the
# direction is a walker difference, so the reachable set is the ensemble's affine hull. With K
# walkers that is a (K-1)-dimensional subspace, FOREVER. The confinement is intrinsic to affine
# invariance rather than to any one package, and it is why N > d matters.

"""
    ess_step(logp, x, lp_x, d, logy; max_expand, maxshrink, rng) -> (x, lp, ne, nc, capped)

One slice-sampling step along `d` from `x`, whose log-density `lp_x` is already known. Steps the
bracket out until both ends leave the slice, then shrinks toward the current point until a sample
inside it is found. Always returns a point; there is no accept/reject.

The expansion budget is split RANDOMLY between the two sides, `J ~ U{0..B-1}`, `K = B-1-J`, which
is Neal (2003) Fig. 3 and is what makes a LIMITED stepping-out reversible. A fixed per-side cap
is not reversible when it binds, so it is not used here even though it is simpler.

Returns the accepted point's log-density so the caller can cache it: re-deriving it on the next
sweep costs exactly one evaluation per walker per sweep, 13-17 % of the total.
"""
function ess_step(logp, x::AbstractVector{Float64}, lp_x::Float64, d::AbstractVector{Float64},
                  logy::Float64; max_expand::Int = 10_000, maxshrink::Int = 10_000,
                  rng = Random.GLOBAL_RNG)
    nd = length(x)
    xp = Vector{Float64}(undef, nd)
    J = rand(rng, 0:(max_expand - 1))          # Neal's randomised split of one shared budget
    K = (max_expand - 1) - J
    u = rand(rng)
    L = -u; R = L + 1.0
    ne = 0; capped = false
    while J > 0
        @. xp = x + L * d
        logp(xp) > logy || break
        L -= 1.0; ne += 1; J -= 1
    end
    J == 0 && (capped = true)
    while K > 0
        @. xp = x + R * d
        logp(xp) > logy || break
        R += 1.0; ne += 1; K -= 1
    end
    K == 0 && (capped = true)
    nc = 0
    while nc < maxshrink
        t = L + (R - L) * rand(rng)
        @. xp = x + t * d
        lp = logp(xp)
        if lp > logy
            return copy(xp), lp, ne, nc, capped
        end
        t < 0 ? (L = t) : (R = t)
        nc += 1
    end
    copy(x), lp_x, ne, nc, true                # degenerate: stay put rather than return garbage
end

"""
    ensemble_slice(logp, x0; nsweep, burn, thin, mu0, tune, rng, verb) -> (chain, diagnostics)

`x0` is `(nwalkers, ndim)`. Red-blue split: each half is updated using directions drawn from the
complementary half, which preserves detailed balance under simultaneous updates.

**The half-membership is RESHUFFLED every sweep** (zeus's `shuffle_ensemble=True`). A fixed
partition is the single most dangerous defect this port had: with `nw = 2d` in low dimension the
complementary half can be one pair, confining the direction to ONE FIXED LINE for the whole run.
Measured on 2-D Rosenbrock at nw = 4 that gave sd 1.75-1.99 against a true 3.1623 — 40 % under-
dispersed — while reporting R̂ = 1.0001. A converged-looking, badly wrong chain.

`mu` is tuned toward an equal count of expansions and contractions, and **frozen once the balance
has held `patience` times**, as zeus does. Without a stopping rule `mu` is an undamped
multiplicative random walk with downward drift, so a LONGER burn-in makes the sampler worse:
median final `mu` fell 0.574 -> 0.317 -> 0.190 for burn = 200 -> 2000 -> 20000, costing 2-5x in
autocorrelation.
"""
function ensemble_slice(logp, x0::Matrix{Float64}; nsweep::Int = 5000, burn::Int = 2000,
                        thin::Int = 10, mu0::Float64 = 1.0, tune::Bool = true,
                        max_expand::Int = 10_000, maxshrink::Int = 10_000,
                        patience::Int = 5, tol::Float64 = 0.05,
                        rng = Random.GLOBAL_RNG, verb::Bool = true)
    nw, nd = size(x0)
    isodd(nw) && error("ensemble_slice needs an even number of walkers (red-blue split)")
    nw < 2nd && @warn "nw = $nw is below the N = 2d = $(2nd) that affine-invariant ensembles need \
                       (Chen 2025); directions span at most nw/2 - 1 dimensions per half."
    X = copy(x0); mu = mu0; half = nw ÷ 2
    lpX = [logp(@view X[i, :]) for i in 1:nw]          # cached current log-densities
    nev = Ref(nw)
    wrapped = z -> (nev[] += 1; logp(z))
    kept = Vector{Matrix{Float64}}()
    nexp = 0; ncon = 0; ncap = 0; tuning = tune; nbal = 0; order = collect(1:nw)
    for sweep in 1:(burn + nsweep)
        shuffle!(rng, order)                           # reshuffle the red/blue membership
        Ne = 0; Nc = 0
        for (a, b) in ((1:half, half+1:nw), (half+1:nw, 1:half))
            act = @view order[a]; comp = @view order[b]
            for i in act
                j = rand(rng, comp); k = rand(rng, comp)
                while k == j && length(comp) > 1; k = rand(rng, comp); end
                d = (2.0 * mu) .* (@view(X[j, :]) .- @view(X[k, :]))   # zeus's factor 2
                if !any(!iszero, d)
                    # Degenerate direction: scale the rescue by the COMPLEMENTARY half only, so
                    # the proposal still depends on nothing in the active half.
                    d = [max(std(@view X[comp, q]), 1e-12) * randn(rng) for q in 1:nd]
                end
                logy = lpX[i] - randexp(rng)
                xn, lpn, ne, nc, cap = ess_step(wrapped, @view(X[i, :]), lpX[i], d, logy;
                                                max_expand = max_expand, maxshrink = maxshrink,
                                                rng = rng)
                X[i, :] = xn; lpX[i] = lpn
                Ne += ne; Nc += nc; cap && (ncap += 1)
            end
        end
        nexp += Ne; ncon += Nc
        if tuning && sweep <= burn && (Ne + Nc) > 0
            mu *= 2 * max(Ne, 1) / (max(Ne, 1) + Nc)   # zeus floors Ne, which bounds the factor
            abs(Ne/(Ne + Nc) - 0.5) < tol ? (nbal += 1) : nothing
            nbal > patience && (tuning = false)
        end
        sweep > burn && (sweep - burn) % thin == 0 && push!(kept, copy(X))
        verb && sweep % 500 == 0 &&
            (@printf("  sweep %5d/%d  mu %.4f%s  Ne/Nc %.2f  evals %d\n", sweep, burn + nsweep,
                     mu, tuning ? "" : " (frozen)", Ne/max(Nc,1), nev[]); flush(stdout))
    end
    chain = Array{Float64,3}(undef, length(kept), nw, nd)
    for (si, Xs) in enumerate(kept); chain[si, :, :] = Xs; end
    chain, (mu = mu, tuned = !tuning, n_expand = nexp, n_contract = ncon, n_eval = nev[],
            n_capped = ncap, nkept = length(kept) * nw, nsweep_kept = length(kept), nwalkers = nw)
end
