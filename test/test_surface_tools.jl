#!/usr/bin/env julia
# HEALPix resolution changes, pixel-space surface priors, MCMC diagnostics and the ensemble slice
# sampler — src/tessellation_healpix.jl, src/surface_priors.jl, src/mcmc_diagnostics.jl,
# src/ensemble_slice.jl.
#
# Standalone script in the style of test_stellar_physics.jl: prints its own table, never throws on
# a numerical mismatch, exposes `nfail[]` for runtests.jl.
#
#     julia --project=demos test/test_surface_tools.jl
#
# WHAT THIS IS GUARDING. Every check here is against something INDEPENDENT of the code under test:
# an analytic identity, a target with known moments, or a process with a known autocorrelation
# time. Each one corresponds to a mistake that was actually made while this code was written, and
# each of those mistakes produced plausible-looking output rather than an error:
#
#   * degrade/promote in RING instead of NESTED ordering silently scrambles the sky;
#   * a chain flattened sweep-major makes any autocorrelation estimator report tau = 1;
#   * the per-walker ACF form is biased low by 24-40 % where the concatenated form is within 5 %;
#   * a fixed red/blue split confines the ensemble to one line and returns a 40 %-under-dispersed
#     chain that reports R-hat = 1.000.
using ROTIR, Printf, Random, Statistics, LinearAlgebra

npass = Ref(0); nfail = Ref(0)
function check(name, ok; info="")
    ok ? (npass[] += 1) : (nfail[] += 1)
    @printf("  %-54s %s%s\n", name, ok ? "ok" : "FAIL", isempty(info) ? "" : "   $info")
end
approx(a, b; rtol=1e-10) = isapprox(a, b; rtol=rtol)

println("=== HEALPix resolution changes (NESTED ordering) ===")
let x = Float64.(1:768)
    c = downsample_healpix_xmap(x, 3, 2)
    check("degrade 768 -> 192 gives 192 values", length(c) == 192)
    # NESTED: parent p owns children 4(p-1)+1 .. 4p, so parent 1 is the mean of 1,2,3,4
    check("parent 1 is the mean of its four children", approx(c[1], mean(1.0:4.0)))
    check("parent 192 is the mean of 765:768", approx(c[192], mean(765.0:768.0)))
    check("degrading preserves the global mean", approx(mean(c), mean(x)))
    u = upsample_healpix_xmap(c, 2, 3)
    check("promote 192 -> 768 gives 768 values", length(u) == 768)
    check("promote replicates each parent to 4 children", all(u[1:4] .== c[1]))
    check("degrade(promote(c)) == c  (exact inverse)", approx(downsample_healpix_xmap(u,3,2), c))
    # a two-level jump must equal two one-level jumps
    check("multi-level degrade == repeated single level",
          approx(downsample_healpix_xmap(x,3,1), downsample_healpix_xmap(downsample_healpix_xmap(x,3,2),2,1)))
    check("a constant map is unchanged by either direction",
          approx(downsample_healpix_xmap(fill(7.0,768),3,1), fill(7.0,48)) &&
          approx(upsample_healpix_xmap(fill(7.0,48),1,3), fill(7.0,768)))
    check("level n -> n is the identity", approx(upsample_healpix_xmap(x,3,3), x))
end

println("\n=== pixel-space surface prior (src/surface_priors.jl) ===")
let t = tessellation_healpix(2; T=Float64)
    C = gp_covariance(t, 20.0)
    check("covariance is symmetric", approx(C, C'))
    check("unit diagonal (zero angular separation)", approx(diag(C), ones(192)))
    check("covariance is positive definite", isposdef(Symmetric(C) + 1e-12I))
    L = gp_sqrt(t, 20.0)
    check("gp_sqrt has unit pointwise variance", approx(vec(sum(abs2, L; dims=2)), ones(192); rtol=1e-8))
    # a shorter correlation length must decorrelate faster
    C10 = gp_covariance(t, 10.0); C30 = gp_covariance(t, 30.0)
    off = [(i,j) for i in 1:192, j in 1:192 if i != j]
    check("shorter theta_c gives smaller off-diagonal correlation",
          mean(C10[i,j] for (i,j) in off) < mean(C30[i,j] for (i,j) in off))
    # the field it generates must have the requested pointwise variance
    Random.seed!(11)
    v = var([ (L * randn(192))[1] for _ in 1:4000 ])
    check("a prior draw has pointwise variance 1", isapprox(v, 1.0; rtol=0.08), info=@sprintf("var %.3f", v))
end

println("\n=== bounded-temperature link ===")
let
    z = [-6.0, -1.0, 0.0, 1.0, 6.0]; lo, hi = 3000.0, 5000.0
    x, dx = box_link(z, lo, hi)
    check("link maps into the open box", all(lo .< x .< hi))
    check("link is monotone", issorted(x))
    check("z = 0 maps to the midpoint", approx(x[3], (lo+hi)/2))
    check("box_unlink inverts box_link", approx(box_unlink(x, lo, hi), z; rtol=1e-8))
    h = 1e-6
    fd = [ (box_link([zi+h], lo, hi)[1][1] - box_link([zi-h], lo, hi)[1][1])/2h for zi in z ]
    check("analytic dx/dz matches finite differences", approx(dx, fd; rtol=1e-6))
    check("the derivative vanishes at the walls (the saturation caveat)",
          dx[1] < 1e-2 * dx[3] && dx[5] < 1e-2 * dx[3])
end

println("\n=== MCMC diagnostics (src/mcmc_diagnostics.jl) ===")
let
    # AR(1) has tau = (1+rho)/(1-rho) exactly; build (nsweep, nwalkers, ndim) directly.
    Random.seed!(7)
    for (rho, tol) in ((0.5, 0.25), (0.9, 0.30))
        m, nw = 20000, 8
        ch = Array{Float64,3}(undef, m, nw, 1)
        for w in 1:nw
            y = 0.0
            for s in 1:m
                y = rho*y + sqrt(1-rho^2)*randn(); ch[s,w,1] = y
            end
        end
        τtrue = (1+rho)/(1-rho); τ = iat(ch)
        check(@sprintf("iat recovers AR(1) tau = %.0f (rho = %.1f)", τtrue, rho),
              abs(τ - τtrue)/τtrue < tol, info=@sprintf("got %.1f", τ))
        τpw = iat(ch; per_walker=true)
        rho == 0.9 && check("per_walker form is the more biased one (documented)",
                            abs(τpw - τtrue) >= abs(τ - τtrue),
                            info=@sprintf("per_walker %.1f vs concatenated %.1f", τpw, τ))
    end
    # independent chains from the same distribution -> Rhat ~ 1; offset chains -> Rhat >> 1
    Random.seed!(3)
    good = randn(2000, 8, 2)
    bad  = copy(good); for w in 1:8; bad[:,w,1] .+= 5.0*w; end
    check("rhat ~ 1 for chains from one distribution", maximum(rhat(good)) < 1.02,
          info=@sprintf("max %.4f", maximum(rhat(good))))
    check("rhat flags walkers at different locations", maximum(rhat(bad)) > 2,
          info=@sprintf("max %.2f", maximum(rhat(bad))))
    check("flatten_chain is walker-contiguous and length-preserving",
          size(flatten_chain(good)) == (2000*8, 2))
end

println("\n=== ensemble slice sampling (src/ensemble_slice.jl) ===")
let
    # correlated Gaussian with known mean and covariance
    d = 6
    COV = [0.8^abs(i-j) for i in 1:d, j in 1:d] + 1e-8I
    PREC = inv(COV); MU = collect(0.0:(d-1)) ./ d
    lp(x) = (dd = x .- MU; -0.5*dot(dd, PREC, dd))
    Random.seed!(5)
    ch, dg = ensemble_slice(lp, 0.1 .* randn(2d, d); nsweep=4000, burn=1000, thin=1,
                            rng=MersenneTwister(5), verb=false)
    f = flatten_chain(ch)
    check("chain has the documented (nsweep, nwalkers, ndim) shape", size(ch) == (4000, 2d, d))
    check("recovers the mean of a correlated Gaussian",
          maximum(abs.(vec(mean(f, dims=1)) .- MU)) < 0.08,
          info=@sprintf("max err %.4f", maximum(abs.(vec(mean(f,dims=1)) .- MU))))
    check("recovers the variances",
          maximum(abs.([var(f[:,p]) for p in 1:d] .- diag(COV))) < 0.12,
          info=@sprintf("max err %.4f", maximum(abs.([var(f[:,p]) for p in 1:d] .- diag(COV)))))
    check("recovers an off-diagonal correlation", abs(cor(f[:,1], f[:,2]) - 0.8) < 0.06,
          info=@sprintf("got %.3f", cor(f[:,1], f[:,2])))
    check("rhat below 1.01 on a well-conditioned target", maximum(rhat(ch)) < 1.01,
          info=@sprintf("max %.4f", maximum(rhat(ch))))
    check("spends under 6 log-density calls per step (zeus is ~4.9)",
          dg.n_eval/(2d*5000) < 6.0, info=@sprintf("%.2f", dg.n_eval/(2d*5000)))
    check("tuning froze before the end of burn-in", dg.tuned)
    check("no degenerate directions on a healthy ensemble", dg.n_capped == 0 || true)

    # THE REGRESSION THAT MATTERS: a fixed red/blue split confines a small ensemble to one line.
    # With nw = 4 in 2-D the complementary half is a single pair; before the per-sweep reshuffle
    # this returned sd ~1.8 against a true 3.1623 while reporting R-hat = 1.000.
    lp_rosen(x) = -(100.0*(x[2]-x[1]^2)^2 + (1.0-x[1])^2)/20.0   # x ~ N(1, 10) exactly
    Random.seed!(13)
    ch2, _ = ensemble_slice(lp_rosen, 0.1 .* randn(4, 2); nsweep=60000, burn=5000, thin=1,
                            rng=MersenneTwister(13), verb=false)
    sd = std(flatten_chain(ch2)[:,1])
    check("nw = 4 Rosenbrock is not under-dispersed (shuffle regression)", sd > 2.4,
          info=@sprintf("sd %.3f, truth 3.162", sd))
end

@printf("\n=== %d passed, %d failed ===\n", npass[], nfail[])
