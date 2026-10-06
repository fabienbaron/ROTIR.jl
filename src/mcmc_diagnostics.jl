# mcmc_diagnostics.jl — convergence and autocorrelation for ensemble chains.
#
# The package has several samplers and had no diagnostics. Both functions take a
# `(nsweep, nwalkers, ndim)` chain, and that shape is not cosmetic: flattening an ensemble chain
# sweep-major gives rows ordered [sweep1 walkers..., sweep2 walkers...], and any estimator that
# then reads contiguous rows as one walker's time series autocorrelates ACROSS walkers and reports
# tau = 1 for everything. That silently invalidated a benchmark during this package's own ESS
# port, so the shape is enforced by the signature.

"""
    rhat(chain) -> Vector

Gelman-Rubin R̂ per dimension from a `(nsweep, nwalkers, ndim)` chain, treating each walker as a
chain. Values above ~1.01 indicate the chain has not converged.

A caveat specific to ENSEMBLE samplers: the walkers are coupled through the shared proposal
directions, so a walker-based R̂ under-reports non-convergence and is blind to an ensemble drifting
as a block. Splitting the chain in time (as `zeus`'s split-R̂ does) is the stricter test.
"""
function rhat(chain::AbstractArray{<:Real,3})
    m, nw, d = size(chain); out = zeros(float(eltype(chain)), d)
    for p in 1:d
        W = mean(var(@view chain[:, w, p]) for w in 1:nw)
        B = m * var([mean(@view chain[:, w, p]) for w in 1:nw])
        out[p] = sqrt(((m - 1)/m * W + B/m) / max(W, eps()))
    end
    out
end

"""
    iat(chain; per_walker=false) -> Float64

Integrated autocorrelation time, averaged over dimensions, from a `(nsweep, nwalkers, ndim)` chain
with Sokal's `t > 5τ` window.

By default the walkers are CONCATENATED into one series per dimension. The `per_walker` variant
subtracts each walker's own mean from a series only tens of τ long and is badly biased low —
measured against AR(1) with known τ it is −24 % at τ = 99 and −40 % at m/τ ≈ 16, where the
concatenated form stays within 5 %.

Both forms are limited by run length, because the window cannot see past the chain: on Neal's
funnel the same estimator gives 98 at 3k sweeps, 377 at 50k and 909 at 1M. **Treat a τ comparable
to `nsweep/50` as a lower bound rather than a measurement.** Neither form models within-sweep
cross-walker correlation, which leaves both roughly a factor 2 low for an ensemble.
"""
function iat(chain::AbstractArray{<:Real,3}; per_walker::Bool = false)
    m, nw, d = size(chain); out = Float64[]
    for p in 1:d
        if per_walker
            acc = zeros(m)
            for w in 1:nw
                y = Float64.(chain[:, w, p]); y .-= mean(y)
                n0 = max(sum(abs2, y)/m, eps())
                acc .+= [sum(y[1:m-t] .* y[1+t:m])/(m-t) for t in 0:m-1] ./ n0
            end
            push!(out, _sokal_window(acc ./ nw))
        else
            y = Float64.(vec(chain[:, :, p])); y .-= mean(y); n = length(y)
            f = [sum(y[1:n-t] .* y[1+t:n]) for t in 0:min(n-1, 10m)] ./ n
            push!(out, _sokal_window(f ./ max(f[1], eps())))
        end
    end
    mean(out)
end

function _sokal_window(acf::AbstractVector)
    τ = 1.0
    for t in 2:length(acf)
        τ += 2 * acf[t]
        t > 5τ && break
    end
    max(τ, 1.0)
end

"Flatten a `(nsweep, nwalkers, ndim)` chain to walker-contiguous rows, for moment estimates."
flatten_chain(c::AbstractArray{<:Real,3}) = reshape(c, size(c,1)*size(c,2), size(c,3))
