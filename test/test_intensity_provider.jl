#!/usr/bin/env julia
# The model-atmosphere intensity provider, I(Teff, logg, mu, lambda)
# (src/intensity_provider.jl).
#
# Standalone script in the style of test_parametric_gradient.jl: prints its own table,
# never throws on a numerical mismatch, exposes `nfail[]` for runtests.jl.
#
#     julia --project=demos test/test_intensity_provider.jl
#
# WHY THE TESTS LOOK LIKE THIS. The derivative of a multilinear interpolant is piecewise
# CONSTANT and discontinuous at every cell wall, so an adaptive finite difference whose
# stencil straddles a wall differences across a kink and reports nonsense — a first pass at
# this file "failed" four checks for exactly that reason, with the analytics correct. The
# fix is not a looser tolerance. Sections [1]-[3] use a grid whose values are EXACTLY
# multilinear (a product of four per-axis linear factors), which the interpolant reproduces
# exactly everywhere, so value and partials both have closed forms and no finite
# differencing is needed at all. Finite differences appear only in [4], against Planck,
# which is smooth — and there lambda is differentiated through a DIMENSIONLESS multiplier,
# because a step chosen for an O(1) argument is meaningless against 2e-6.

using ROTIR, LinearAlgebra, Printf

const HAVE_FD = try; @eval using FiniteDifferences; true; catch; false; end
const HAVE_ZYGOTE = try; @eval using Zygote; true; catch; false; end
const FDM = HAVE_FD ? central_fdm(5, 1) : nothing

relerr(a, b; scale = 0.0) = norm(a .- b) / max(norm(a), norm(b), scale, eps())
# The glyph is judged against the SAME tolerance as the verdict. A fixed threshold printed
# a cross next to the O(h^2) interpolation-error checks in [6], which pass at their own
# (deliberately looser) tolerance - a passing check rendered as a failure is how a real
# failure gets overlooked.
status(r, tol) = r < tol ? "\u2713" : (r < 10 * tol ? "~" : "\u2717")
npass = Ref(0); nfail = Ref(0); fd_ran = Ref(false); zygote_ran = Ref(false)
function chk(label, r; tol = 1e-9)
    ok = r < tol
    ok ? (npass[] += 1) : (nfail[] += 1)
    @printf("  %-52s relerr=%.3e  %s\n", label, r, status(r, tol))
    return ok
end
cb(label, ok) = chk(label, ok ? 0.0 : 1.0)

# A grid whose values are EXACTLY multilinear: a product of four per-axis linear factors.
# Multilinear interpolation reproduces such a function exactly at every point, not just on
# nodes, so the interpolant carries NO discretisation error and its partials have closed
# forms. That is what makes sections [1]-[3] a clean test of the derivatives.
aT(t)=1+t/1000; ag(g)=2+g; am(m)=0.5+m; al(l)=1+l*1e6
ml(t,g,m,l)=aT(t)*ag(g)*am(m)*al(l)
dT_(t,g,m,l)=(1/1000)*ag(g)*am(m)*al(l)
dg_(t,g,m,l)=aT(t)*1*am(m)*al(l)
dm_(t,g,m,l)=aT(t)*ag(g)*1*al(l)
dl_(t,g,m,l)=aT(t)*ag(g)*am(m)*1e6

const Ta=collect(range(4000.0,9000.0,length=9)); const ga=collect(range(2.5,4.5,length=7))
const ma=collect(range(1e-3,1.0,length=11));     const la=collect(range(1.5e-6,2.5e-6,length=6))
const gml = RectGrid4(Ta,ga,ma,la,[ml(t,g,m,l) for t in Ta,g in ga,m in ma,l in la])
const pml = TabulatedProvider(gml)

println("\n[1] exactly-multilinear grid: value and all four partials are EXACT")
for pt in ([5123.7,3.31,0.417,1.93e-6], [4001.0,2.51,0.002,1.51e-6], [8999.0,4.49,0.999,2.49e-6])
    v,d1,d2,d3,d4 = interp4_and_grad(gml, pt...)
    chk("value  at $(round.(pt,sigdigits=4))", relerr(v, ml(pt...)); tol=1e-13)
    chk("  d/dTeff",   relerr(d1, dT_(pt...)); tol=1e-11)
    chk("  d/dlogg",   relerr(d2, dg_(pt...)); tol=1e-11)
    chk("  d/dmu",     relerr(d3, dm_(pt...)); tol=1e-11)
    chk("  d/dlambda", relerr(d4, dl_(pt...)); tol=1e-11)
end

if HAVE_ZYGOTE
println("\n[2] Zygote rrule reproduces the analytic partials exactly")
Tv=[5123.7,6210.0,7788.0]; lgv=[3.31,2.9,4.1]; mv=[0.417,0.82,0.11]; lv=[1.93e-6,2.1e-6,1.7e-6]
Ia,dT,dg,dm,dl = provider_intensity_and_derivs(pml,Tv,lgv,mv,lv)
# Loss = sum(w .* I) makes each cotangent exactly w .* dI/dq, so Zygote and the analytic
# derivatives can be compared with NO finite differencing anywhere.
w=[0.3,-1.7,2.1]
gz = Zygote.gradient((a,b,c,d)->dot(w,provider_map(pml,a,b,c,d)), Tv,lgv,mv,lv)
chk("Zygote d/dTeff   == w .* dI/dTeff",   relerr(gz[1], w.*dT); tol=1e-13)
chk("Zygote d/dlogg   == w .* dI/dlogg",   relerr(gz[2], w.*dg); tol=1e-13)
chk("Zygote d/dmu     == w .* dI/dmu",     relerr(gz[3], w.*dm); tol=1e-13)
chk("Zygote d/dlambda == w .* dI/dlambda", relerr(gz[4], w.*dl); tol=1e-13)

println("\n[3] a scalar lambda pulls back as a SUM, not elementwise")
Is,_,_,_,dls = provider_intensity_and_derivs(pml,Tv,lgv,mv,2.0e-6)
gz2 = Zygote.gradient(l->dot(w,provider_map(pml,Tv,lgv,mv,l)), 2.0e-6)
chk("scalar cotangent == dot(w, dI/dlambda)", relerr(gz2[1], dot(w,dls)); tol=1e-13)

zygote_ran[] = true
else
    @warn "Zygote unavailable - sections [2] and [3] skipped"
end

println("\n[4] Planck: dI/dTeff and dI/dlambda against well-scaled FD")
# lambda ~ 2e-6, so FD must differentiate a DIMENSIONLESS multiplier; a step chosen for an
# O(1) argument is meaningless against 2e-6 and was what made the first attempt fail.
Tk=collect(range(4500,8500,length=13)); λ0=1.65e-6
pp=PlanckProvider()
_,dTp,dgp,dmp,dlp = provider_intensity_and_derivs(pp,Tk,nothing,nothing,λ0)
if HAVE_FD
    fdT=[FiniteDifferences.grad(FDM,t->provider_intensity(pp,[t],nothing,nothing,λ0)[1],Tk[i])[1] for i in eachindex(Tk)]
    chk("dI/dTeff vs FD", relerr(dTp,fdT); tol=1e-9)
    fds=FiniteDifferences.grad(FDM,s->sum(provider_intensity(pp,Tk,nothing,nothing,λ0*s)),1.0)[1]
    chk("dI/dlambda vs FD (scaled)", relerr(sum(dlp)*λ0, fds); tol=1e-9)
    fd_ran[] = true
else
    @warn "FiniteDifferences unavailable - the two FD checks in [4] skipped"
end
cb("dI/dlogg, dI/dmu exactly zero", all(iszero,dgp) && all(iszero,dmp))
cb("matches intensity(:planck)", provider_intensity(pp,Tk,nothing,nothing,λ0) ≈ intensity(Tk,:planck,λ0))
# The identity dB/dlambda == dB/dT * (T/lambda), independent of both implementations.
chk("dB/dlam == dB/dT * T/lam", relerr(dlp, dTp .* Tk ./ λ0); tol=1e-13)

println("\n[5] guards and consistency")
cb("mu=0 node refused", try RectGrid4(Ta,ga,[0.0;ma],la,zeros(9,7,12,6)); false catch; true end)
cb("unsorted axis refused", try RectGrid4(reverse(Ta),ga,ma,la,gml.values); false catch; true end)
cb("shape mismatch refused", try RectGrid4(Ta,ga,ma,la,zeros(2,2,2,2)); false catch; true end)
cb("strict rejects out-of-range", try provider_intensity(TabulatedProvider(gml;strict=true),[99000.0],[3.5],[0.5],2e-6); false catch; true end)
cb("non-strict clamps", length(provider_intensity(pml,[99000.0],[3.5],[0.5],2e-6))==1)
sp = default_star_params(2; ldtype=0, d=16.8, rpole=0.849, frac_escapevel=0.92, rotation_period=1/1.12)
cb("ldtype=0 + tabulated clean", isempty(check_provider_consistency(pml,sp)))
cb("ldtype=3 + tabulated REFUSED", length(check_provider_consistency(pml,merge(sp,(ldtype=3,))))==1)
cb("no distance REFUSED", any(m->occursin("logg",m), check_provider_consistency(pml, Base.structdiff(sp,NamedTuple{(:d,)}))))
cb("Planck + ldtype=3 fine", isempty(check_provider_consistency(pp,merge(sp,(ldtype=3,)))))

println("\n[6] smooth non-multilinear grid: interpolation error is O(h^2), as expected")
g2,f2 = analytic_test_grid()
pt=[5123.7,3.31,0.417,1.93e-6]
cb("on-node exact", interp4_and_grad(g2,g2.Teff[3],g2.logg[2],g2.μ[5],g2.λ[4])[1] ≈ f2(g2.Teff[3],g2.logg[2],g2.μ[5],g2.λ[4]))
chk("off-node within O(h^2)", relerr(interp4_and_grad(g2,pt...)[1], f2(pt...)); tol=3e-2)
g3,_ = analytic_test_grid(nT=33,ng=25,nm=41,nl=21)
chk("4x finer grid -> ~16x smaller error", relerr(interp4_and_grad(g3,pt...)[1], f2(pt...)); tol=2e-3)

println("\n[10] Float32: the OUTPUT eltype follows the INPUTS, not the grid")
# THE DEFAULT COMBINATION. ROTIR meshes default to Float32 and `load_intensity_grid` defaults
# to Float64, so a provider returning the GRID's eltype would silently promote every Float32
# model to Float64 — and the whole downstream visibility computation with it. `intensity` in
# src/intensity.jl sets the precedent by returning `similar(x)`.
let g64 = RectGrid4([4000.0,9000.0], [2.0,5.0], [1e-3,1.0], [1.5e-6,2.5e-6],
                    reshape(collect(1.0:16.0), 2,2,2,2)),
    g32 = RectGrid4(Float32[4000,9000], Float32[2,5], Float32[1e-3,1],
                    Float32[1.5e-6,2.5e-6], reshape(collect(1.0f0:16.0f0), 2,2,2,2))
    T32 = Float32[6000, 7000]; L32 = Float32[3, 4]; M32 = Float32[0.5, 0.9]
    for (nm, g) in (("Float64 grid", g64), ("Float32 grid", g32))
        I = provider_intensity(TabulatedProvider(g), T32, L32, M32, 2.0f-6)
        cb("$nm + Float32 inputs -> Float32", eltype(I) === Float32)
        D = provider_intensity_and_derivs(TabulatedProvider(g), T32, L32, M32, 2.0f-6)
        cb("  ... and all four derivatives too", all(d -> eltype(d) === Float32, D))
    end
    # A Float64 model must still stay Float64 — the rule is "follow the inputs", not "narrow".
    I64 = provider_intensity(TabulatedProvider(g32), [6000.0,7000.0], [3.0,4.0], [0.5,0.9], 2.0e-6)
    cb("Float32 grid + Float64 inputs -> Float64", eltype(I64) === Float64)
    # Interpolation still runs in the grid's precision, so narrowing costs only the last step.
    a = provider_intensity(TabulatedProvider(g64), T32, L32, M32, 2.0f-6)
    b = provider_intensity(TabulatedProvider(g64), Float64.(T32), Float64.(L32),
                           Float64.(M32), 2.0e-6)
    chk("Float32 result matches the Float64 one", relerr(Float64.(a), b); tol = 1e-6)
end

@printf("\n=== %d passed, %d failed ===\n", npass[], nfail[])
