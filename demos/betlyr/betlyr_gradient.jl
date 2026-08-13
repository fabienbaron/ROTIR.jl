# betlyr_gradient.jl — TRUE ANALYTIC gradient of the β Lyrae χ², no AD anywhere.
#
# Included by betlyr_orbit_fit_pigeons.jl when EXPLORER=mala. Every derivative below is in
# closed form; ForwardDiff appears only in the regression test, as the correctness oracle.
#
#   dχ²/dp = Σ_epochs real( Jᵀ_epoch · g_cvis_epoch )        [OITOOLS convention: TRANSPOSE,
#                                                              no conjugation — see cvis_chi2]
# The pieces:
#   g_cvis          OITOOLS' hand-written cvis_to_chi2_fg (value and adjoint in one sweep)
#   ∂v₁/∂ud         2J₀/t − 4J₁/t², via the recurrence J₂ = 2J₁/t − J₀ so the GENERIC-order
#                   besselj(2,·) is never called (it costs 26.7 µs against 1.9 µs for
#                   besselj1 — a 12x penalty), reusing the J₁ already computed for v₁
#   ∂v₂/∂(fw,ar,pa) closed form for the inclined elliptical Gaussian
#   ∂cvis/∂f        (v₂·ph − v₁)/(1+f)²
#   ∂(ras,decs)/∂q  orbit_jac below — Ω is a pure rotation, a is linear, and T0/Ṗ chain
#                   through Kepler's ∂E/∂M = 1/(1 − e·cos E) (the implicit function theorem,
#                   same identity that backs `kepler_E`'s rrule in src/orbits.jl)
#
# MEASURED: value+gradient in 2.69x the primal, against 8.3x for ForwardDiff and 20x for
# central differences. Matches ForwardDiff to 6.6e-16 on θ and 6.9e-16 in z-space.
#
# MAINTENANCE WARNING. This must be re-derived whenever `ud_vis`, `ellgauss_vis`, or the
# ephemeris changes. It was first committed with a factor-of-2 error in ∂v₁/∂ud that left
# nine of ten parameters exact — only the ForwardDiff cross-check caught it. `test_gradient()`
# at the bottom is that check; run it after any edit to the forward model.

using SpecialFunctions: besselj0, besselj1

const OIT_ = parentmodule(ROTIR.OIdata)      # OITOOLS, for cvis_to_chi2_fg
const K_   = -2pi / MAS_PER_RAD              # phasor scale
const GS_  = 2*sqrt(2*log(2))                # FWHM -> sigma

# d/dt[2 J1(t)/t] = (J0 - J2)/t - 2 J1/t^2, and the recurrence J2 = 2 J1/t - J0 removes the
# GENERIC-ORDER besselj(2,·) — which cost 26.7 us against 1.9 us for besselj1, a 12x penalty
# for the general-order algorithm. Reusing the J1 already computed for v1 leaves:
@inline dud_kernel(J0, J1, t) = 2*J0/t - 4*J1/t^2

# Fused reduction: real(J' * g) with no temporaries. The broadcast form allocated 19 arrays
# (33.6 KiB) per call and ran 24x slower; there are ten of these per epoch.
@inline function dot_r(g, z)
    s = 0.0
    @inbounds @simd for k in eachindex(g)
        s += real(g[k])*real(z[k]) - imag(g[k])*imag(z[k])
    end
    return s
end

function orbit_jac(a, i_deg, Om_deg, T0, dP, ts; P=P_ORB, e=E_ORB, om_deg=OMEGA_PERI)
    Ω = Om_deg*pi/180; inc = i_deg*pi/180; ω = om_deg*pi/180
    β = sqrt(1 - e^2); D2R = pi/180
    cΩ, sΩ, ci, si = cos(Ω), sin(Ω), cos(inc), sin(inc)
    cω, sω = cos(ω), sin(ω)
    L1 =  cΩ*cω - sΩ*sω*ci;  M1 =  sΩ*cω + cΩ*sω*ci
    L2 = -cΩ*sω - sΩ*cω*ci;  M2 = -sΩ*sω + cΩ*cω*ci
    bp = (i=i_deg, Ω=Om_deg, ω=om_deg, P=P, a=a, e=e, T0=T0, q=Q_BIN, dP=dP, dω=0.0)
    E  = compute_eccentric_anomaly(bp, ts)
    cE, sE = cos.(E), sin.(E)
    ras  = a .* (M1 .* cE .+ β*M2 .* sE .- e*M1)
    decs = a .* (L1 .* cE .+ β*L2 .* sE .- e*L1)
    # ∂/∂E
    dras_dE  = a .* (-M1 .* sE .+ β*M2 .* cE)
    ddec_dE  = a .* (-L1 .* sE .+ β*L2 .* cE)
    dE_dM    = 1 ./ (1 .- e .* cE)                       # implicit function theorem
    # ∂M/∂T0, ∂M/∂dP from the quadratic ephemeris
    τ = ts .- T0
    if abs(dP) > 1e-12
        R = 1 .+ 2dP .* τ ./ P;  s = sqrt.(R)
        dM_dτ  = 4pi ./ (P .* (1 .+ s)) .- 4pi .* τ .* dP ./ (P^2 .* s .* (1 .+ s).^2)
        dM_ddP = -4pi .* τ.^2 ./ (P^2 .* s .* (1 .+ s).^2)
    else
        dM_dτ  = fill(2pi/P, length(τ));  dM_ddP = zeros(length(τ))
    end
    dM_dT0 = -dM_dτ
    # ∂/∂i (only ci moves)
    dci = -si*D2R
    dM1_di =  cΩ*sω*dci;  dM2_di =  cΩ*cω*dci
    dL1_di = -sΩ*sω*dci;  dL2_di = -sΩ*cω*dci
    dras_di  = a .* (dM1_di .* cE .+ β*dM2_di .* sE .- e*dM1_di)
    ddec_di  = a .* (dL1_di .* cE .+ β*dL2_di .* sE .- e*dL1_di)
    return (ras=ras, decs=decs,
            dras  = (a=ras./a, i=dras_di, Ω= decs.*D2R, T0=dras_dE.*dE_dM.*dM_dT0, dP=dras_dE.*dE_dM.*dM_ddP),
            ddecs = (a=decs./a, i=ddec_di, Ω=-ras.*D2R, T0=ddec_dE.*dE_dM.*dM_dT0, dP=ddec_dE.*dE_dM.*dM_ddP))
end


function chi2_and_grad(θ)
    # θ carries an 11th entry (P) since the FITP consistency check was added. Julia's
    # destructuring silently DROPS extras, so taking ten here would have used P_ORB in the
    # Jacobian while the forward model used θ[11] — a wrong gradient with no error. Take P
    # explicitly and use it, so the χ² VALUE is always right.
    a, inc, Ω, T0, dP, ud, fw, ar, pa, f, P = θ
    # No ∂χ²/∂P yet. With Ṗ ≠ 0 the mean anomaly is M = 4πτ/(P(1+s)), s = √(1+2Ṗτ/P), so
    # ∂M/∂P is a genuine derivation rather than a rescaling of ∂M/∂T0 (which is what it
    # collapses to only in the Ṗ → 0 limit). Refuse rather than return a gradient that is
    # silently short by one component.
    length(θ) == NPAR || error("chi2_and_grad: expected $NPAR parameters, got $(length(θ))")
    abs(P - P_ORB) < 1e-12 || error("""
        chi2_and_grad has no ∂/∂P: the analytic gradient predates the FITP consistency
        check. Use FITP=0 for any gradient-based sampler (Pigeons/AutoMALA), or run the
        FITP=1 check through UltraNest, which is gradient-free.""")
    g = zeros(10); tot = 0.0
    bp(q) = (i = q[2], Ω = q[3], ω = OMEGA_PERI, P = P, a = q[1], e = E_ORB,
             T0 = q[4], q = Q_BIN, dP = q[5], dω = 0.0)
    qorb = [a, inc, Ω, T0, dP]
    inv1f = 1/(1+f); σ = fw/MAS_PER_RAD/GS_; φp = deg2rad(pa)
    for i in 1:NEP
        ts = TSRT[i]; uv = UV[i]; ti = TIDX[i]; uu = UU[i]; vv = VV[i]; ρ = RHO[i]
        # --- orbit block and its Jacobian, hand-derived (see orbit_jac above) ---
        OJ  = orbit_jac(a, inc, Ω, T0, dP, ts; P = P)
        ras = OJ.ras; decs = OJ.decs
        # --- components ---
        t  = (ud/MAS_PER_RAD*pi) .* ρ .+ 1e-10
        J0 = besselj0.(t); J1 = besselj1.(t)
        v1 = 2 .* J1 ./ t
        up =  uu .* cos(φp) .+ vv .* sin(φp); vp = -uu .* sin(φp) .+ vv .* cos(φp)
        ρe2 = (up .* ar).^2 .+ vp.^2
        v2 = exp.(-2 .* (pi*σ)^2 .* ρe2)
        ph = cis.(K_ .* (view(uv,1,:) .* ras[ti] .+ view(uv,2,:) .* decs[ti]))
        cv = (v1 .+ f .* v2 .* ph) .* inv1f
        # --- adjoint source from OITOOLS ---
        c2, gc = OIT_.cvis_to_chi2_fg(cv, DATA[i]); tot += c2
        # --- J' * g_cvis, real part (transpose, NO conjugation) ---
        g[6]  += dot_r(gc, (dud_kernel.(J0, J1, t) .* (pi/MAS_PER_RAD) .* ρ) .* inv1f)   # ud
        w2     = f .* ph .* inv1f
        g[7]  += dot_r(gc, w2 .* (v2 .* (-4*pi^2*σ*ρe2 ./ (MAS_PER_RAD*GS_))))                # fwhm
        g[8]  += dot_r(gc, w2 .* (v2 .* (-4*(pi*σ)^2 .* up.^2 .* ar)))                       # ratio
        g[9]  += dot_r(gc, w2 .* (v2 .* (-2*(pi*σ)^2 .* 2 .* up .* vp .* (ar^2-1) * (pi/180))))
        g[10] += dot_r(gc, (v2 .* ph .- v1) .* inv1f^2)                                      # f
        dφ     = f .* v2 .* ph .* im .* inv1f
        for k in 1:5                                                                     # orbit
            dras = OJ.dras[k]; ddec = OJ.ddecs[k]
            g[k] += dot_r(gc, dφ .* (K_ .* (view(uv,1,:) .* dras[ti] .+ view(uv,2,:) .* ddec[ti])))
        end
    end
    return tot, g
end



# ---------------------------------------------------------------------------------------
# z-space: θ = LO + (HI−LO)·σ(z), logpost = −0.5χ²(θ) + Σ[log(HI−LO) + log u + log(1−u)]
#   dθ_k/dz_k      = (HI−LO)·u(1−u)
#   d(logjac)/dz_k = u(1−u)·(1/u − 1/(1−u)) = 1 − 2u
# ---------------------------------------------------------------------------------------
_sig(x) = x >= 0 ? 1/(1+exp(-x)) : (e = exp(x); e/(1+e))

"""
    logpost_and_grad_z(z, idx, pos) -> (logpost, gradient)

Log posterior and its **analytic** gradient in the unconstrained sampling space.
"""
function logpost_and_grad_z(z, idx, pos)
    u = _sig.(z)
    θ = copy(THETA_LIT)
    @inbounds for (k,p) in enumerate(idx); θ[p] = LO[p] + (HI[p]-LO[p])*u[k]; end
    c, gθ = chi2_and_grad(θ)
    c >= 1e12 && return (-Inf, zeros(length(z)))
    lj = 0.0
    @inbounds for (k,p) in enumerate(idx); lj += log(HI[p]-LO[p]) + log(u[k]) + log1p(-u[k]); end
    v = -0.5c + lj
    isfinite(v) || return (-Inf, zeros(length(z)))
    g = similar(z)
    @inbounds for (k,p) in enumerate(idx)
        g[k] = -0.5*gθ[p]*(HI[p]-LO[p])*u[k]*(1-u[k]) + (1 - 2u[k])
    end
    return (v, g)
end

"Cross-check the analytic gradient against ForwardDiff. Run after ANY forward-model edit."
function test_gradient(; rtol = 1e-10)
    @eval Main using ForwardDiff
    FD = getfield(Main, :ForwardDiff)
    idx = free_indices(); pos = free_positions(idx)
    z = theta_to_z(THETA_LIT, idx)
    _, ga = logpost_and_grad_z(z, idx, pos)
    gf = FD.gradient(w -> -0.5*chi2_total(z_to_theta_ad(w, pos)) + log_jacobian(w, idx), z)
    err = maximum(abs.(ga .- gf) ./ max.(abs.(gf), 1e-8))
    println(err < rtol ? "analytic gradient OK (max rel err $(err))" :
                         "*** ANALYTIC GRADIENT MISMATCH: max rel err $(err) ***")
    return err < rtol
end
