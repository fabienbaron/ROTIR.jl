#!/usr/bin/env julia
# =======================================================================================
# Does velocity-resolved data break the rpole/fev degeneracy?
# =======================================================================================
# demos/betcas_fisher_predicted_ld.jl showed that a model-atmosphere intensity removes the
# beta/ld1 degeneracy outright. It left rpole/fev at corr = -0.995 as the strongest remaining
# correlation, and the obvious hope was that differential phase would break it: a differential
# phase measures the VELOCITY field and the rotation axis, which V^2 and closure phase cannot
# see at all.
#
# THE ANSWER IS NO, and the reason is structural rather than a matter of data quality:
#
#     v is proportional to Omega * R_eq,   and   R_eq = rpole * f(fev)
#
# so the velocity field constrains the same PRODUCT that the projected geometry does. Adding
# it tightens everything by 1.4-2.6x but does not rotate the error ellipse. Confirmed against
# an artificially narrow line (best case) as well as the real Br-11, and across every
# inclination from 20 to 90 degrees, where corr(rpole, fev) never leaves [-0.999, -0.96].
#
# WHAT THAT DOES AND DOES NOT MEAN. A correlation of -0.995 is not the same as a poor
# measurement. The sigmas printed below are MARGINALISED — they already account for the
# correlation — and at beta Cas's inclination they are 0.2 % on rpole and 0.27 % on fev. Both
# parameters are well determined. What a high correlation does mean is that the pair is
# vulnerable to a COMMON SYSTEMATIC: an error in the gravity law, in the limb darkening, or in
# the mesh propagates along that one direction and biases both together. Fisher cannot see
# that, and it is precisely why predicting the limb darkening rather than fitting it matters.
#
# Two practical results fall out:
#   * Br-11 works about as well as an ideal narrow line, despite being Stark-broadened to
#     428 km/s against beta Cas's 73 km/s vsini. A broad line spreads the signal over more
#     channels, so more baseline x channel measurements contribute even though each is diluted.
#   * MIRC-X's R = 190 cannot do this at all: its channels are 88 A wide against a 24 A line.
#     GRAVITY's R = 4000 (4.2 A) can.
#
#   julia --project=demos demos/betcas_fisher_diffphase.jl
#
# Needs demos/data/betcas_Br11_korg.fits and betcas_H_korg.fits; both are built by
# demos/betcas_fisher_predicted_ld.jl and by the `build_korg_grid` call recorded there.
# =======================================================================================

using ROTIR, Printf, LinearAlgebra, Statistics
const R = @__DIR__
const GH   = load_intensity_grid(joinpath(R,"data","betcas_H_korg.fits"))
const GBR  = load_intensity_grid(joinpath(@__DIR__,"data","betcas_Br11_korg.fits"))
const PH   = TabulatedProvider(GH;  name="H broadband")
const PBR  = TabulatedProvider(GBR; name="Br-11")
const DATA = readoifits_multiepochs([joinpath(R,"data",
  "MEDIAN5.MIRCX_L2.2025Oct30.HD_432.MIRCX_IDL.bet_Cas.AVG10m.oifits")];T=Float64)[1,1]
const DPC1 = readoifits_multiepochs([joinpath(R,"data",
  "MEDIAN5.MIRCX_L2.2025Oct30.HD_432.MIRCX_IDL.bet_Cas.AVG10m.oifits")];T=Float64,
  polychromatic=true)[1,1]
const TESS = tessellation_healpix(4,T=Float64)
const BAND = band_of(DATA)
# The V2 baselines of one spectral channel: these are the physical baselines a spectrograph
# would deliver a differential phase on, one per channel.
const UVV  = DPC1.uv[:, DPC1.indx_v2]
const LREF = band_of(DPC1)
const NBL  = size(UVV,2)

# --- a NARROW-line twin of the Br-11 grid ----------------------------------------------
# Same Teff/logg/mu structure and the same 20% depth, but 30 km/s FWHM instead of Br-11's
# Stark-broadened 428 km/s. This separates "does velocity-resolved data carry the
# information" from "is the only line available good enough to deliver it".
const λBR = GBR.λ; const λ0 = 16811.0e-10
function narrow_grid(fwhm_kms, depth)
    v = Array{Float64,4}(undef, size(GBR.values))
    σ = fwhm_kms/2.355/299792.458*λ0
    p = [1 - depth*exp(-((l-λ0)/σ)^2/2) for l in λBR]
    nl = length(λBR)
    for i in axes(v,1), j in axes(v,2), k in axes(v,3)
        # Linear continuum through the window edges, so the real slope and LD survive.
        c1 = mean(GBR.values[i,j,k,1:20]); c2 = mean(GBR.values[i,j,k,nl-19:nl])
        cont = c1 .+ (c2-c1) .* (λBR .- λBR[1]) ./ (λBR[end]-λBR[1])
        v[i,j,k,:] = cont .* p
    end
    TabulatedProvider(RectGrid4(GBR.Teff,GBR.logg,GBR.μ,λBR,v); name="narrow $(fwhm_kms) km/s")
end
const PNAR = narrow_grid(30.0, 0.20)

base(; kw...) = default_star_params(2; rpole=0.849,d=16.8,frac_escapevel=0.92,
  rotation_period=1/1.12,tpole=7208.0,inclination=19.9,position_angle=-7.09,
  beta=0.25,gravity_law=2,ldtype=0,kw...)

# cvis at one wavelength with an explicit uv, without building an OIdata.
function cvis_at(st, state, prov, λ, uv)
    I = channel_intensity(prov, state, λ)
    ix = st.index_quads_visible
    xw = I[ix] .* st.vis_weights[ix] .* st.ldmap[ix]
    pjx = Matrix(st.proj_west[ix,:]); pjy = Matrix(st.proj_north[ix,:])
    kx = uv[1,:] .* (-pi/(180*3600000)); ky = uv[2,:] .* (pi/(180*3600000))
    F = ROTIR.type3_cvis(pjx,pjy,xw,kx,ky)
    return F ./ dot(setup_polyflux_single(pjx,pjy), xw)
end

function model(θ; prov=nothing, λchan=nothing, win=nothing)
    sp = base(d=θ[6],rpole=θ[1],frac_escapevel=θ[2],inclination=θ[3],
              position_angle=θ[4],beta=θ[5])
    st = create_star(TESS,sp,0.0); state = surface_state(st,sp)
    # broadband block: real uv, H-band grid
    Ib = channel_intensity(PH, state, BAND)
    v2,_,t3p = cvis_to_obs(ROTIR.fused_cvis(Ib,st,DATA;intensity_model=:linear), DATA)
    prov === nothing && return (v2, t3p, nothing)
    # spectroscopic block: uv scaled to each channel, since uv = B/lambda
    cvs = [cvis_at(st,state,prov,λc, UVV .* (LREF/λc)) for λc in λchan]
    _, dp = differential_observables(cvs, λchan, win)
    inl = findall(win[1] .<= λchan .<= win[2])
    return (v2, t3p, vec(dp[:, inl]))
end

whiten(o; σdp=0.3) = o[3]===nothing ? vcat(o[1]./DATA.v2_err, o[2]./DATA.t3phi_err) :
    vcat(o[1]./DATA.v2_err, o[2]./DATA.t3phi_err, o[3]./σdp)
function dwhiten(p,m; σdp=0.3)
    a = vcat((p[1].-m[1])./DATA.v2_err, mod360.(p[2].-m[2])./DATA.t3phi_err)
    p[3]===nothing ? a : vcat(a, mod360.(p[3].-m[3])./σdp)
end
function fisher(f,θ; frac=1e-6, σdp=0.3)
    n=length(θ); J=Matrix{Float64}(undef,length(whiten(f(θ);σdp=σdp)),n)
    for j in 1:n
        h=frac*max(abs(θ[j]),1e-3); tp=copy(θ);tp[j]+=h; tm=copy(θ);tm[j]-=h
        J[:,j]=dwhiten(f(tp),f(tm);σdp=σdp)./(2h)
    end
    J'*J
end
cm(C)=[C[i,j]/sqrt(C[i,i]*C[j,j]) for i in axes(C,1), j in axes(C,2)]
const PRI=Diagonal([0.0,0,0,0,0,1/0.1^2])
const NM=["rpole","fev","inc","PA","beta","d"]
const θ0=[0.849,0.92,19.9,-7.09,0.25,16.8]

# GRAVITY R=4000 channels across Br-11, plus continuum on both sides.
Rres=4000.0; dλ=λ0/Rres
λchan=collect(range(λ0-70e-10, λ0+70e-10, step=dλ))
win=(λ0-30e-10, λ0+30e-10)
@printf("%d baselines, %d channels (%.1f A each, R=%.0f), %d in-line -> %d dphase obs\n",
        NBL, length(λchan), dλ*1e10, Rres, count(win[1].<=λchan.<=win[2]),
        NBL*count(win[1].<=λchan.<=win[2]))

function show(title,F)
    C=inv(F); Rm=cm(C)
    @printf("\n%s\n  %-6s %11s   %s\n", title, "param","sigma", join(rpad.(NM,8)))
    for i in 1:6
        @printf("  %-6s %11.5g   %s\n", NM[i], sqrt(C[i,i]),
                join([rpad(@sprintf("%+.3f",Rm[i,j]),8) for j in 1:6]))
    end
    C,Rm
end
println("="^92); println("Does differential phase break the rpole/fev degeneracy?"); println("="^92)
F0=fisher(θ->model(θ), θ0)
C0,R0=show("V2 + T3phi only (broadband)", F0+PRI)
Fb=fisher(θ->model(θ;prov=PBR,λchan=λchan,win=win), θ0)
Cb,Rb=show("+ differential phase across Br-11 (FWHM 428 km/s, real line)", Fb+PRI)
Fn=fisher(θ->model(θ;prov=PNAR,λchan=λchan,win=win), θ0)
Cn,Rn=show("+ differential phase across a NARROW 30 km/s line (best case)", Fn+PRI)

println("\n"*"="^92); println("Verdict"); println("="^92)
@printf("  %-22s %12s %12s %12s\n","quantity","broadband","+Br-11","+narrow")
@printf("  corr(rpole,fev)        %+12.4f %+12.4f %+12.4f\n", R0[1,2],Rb[1,2],Rn[1,2])
for (i,lab) in enumerate(("sigma(rpole)","sigma(fev)","sigma(inc)","sigma(PA)","sigma(beta)"))
    k=(1,2,3,4,5)[i]
    @printf("  %-22s %12.5g %12.5g %12.5g   (%.2fx, %.2fx)\n", lab,
            sqrt(C0[k,k]),sqrt(Cb[k,k]),sqrt(Cn[k,k]),
            sqrt(C0[k,k])/sqrt(Cb[k,k]), sqrt(C0[k,k])/sqrt(Cn[k,k]))
end

println("\n"*"="^92)
println("Is rpole/fev an INCLINATION problem? Scan, broadband and +narrow-line dphase")
println("="^92)
@printf("%-6s %-10s | %-32s | %-32s\n","inc","sin(i)",
        "broadband: corr  sig(rp)  sig(fev)","+dphase:   corr  sig(rp)  sig(fev)")
for inc in (19.9, 30.0, 45.0, 60.0, 75.0, 90.0)
    t=copy(θ0); t[3]=inc
    C1=inv(fisher(θ->model(θ), t)+PRI); R1=cm(C1)
    C2=inv(fisher(θ->model(θ;prov=PNAR,λchan=λchan,win=win), t)+PRI); R2=cm(C2)
    @printf("%-6.1f %-10.3f | %+7.4f %9.5f %9.5f | %+7.4f %9.5f %9.5f\n",
            inc, sind(inc), R1[1,2], sqrt(C1[1,1]), sqrt(C1[2,2]),
            R2[1,2], sqrt(C2[1,1]), sqrt(C2[2,2]))
end
