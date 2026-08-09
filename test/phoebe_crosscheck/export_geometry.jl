ENV["MPLBACKEND"] = "Agg"
using ROTIR, DelimitedFiles, LinearAlgebra, Printf

OUT = joinpath(@__DIR__, "data")
mkpath(OUT)

P=4.0145; a=1.54; e=0.123; T0=2454189.40; q=0.6188; i_orb=116.0; Om=309.938; om=255.0
R1, R2, T1, T2 = 0.4469, 0.2259, 25300.0, 20585.0
sph(r,t) = (surface_type=0, radius=r, tpole=t, ldtype=1, ld1=0.0, ld2=0.0,
            inclination=180.0-i_orb, position_angle=Om-180.0, rotation_period=P)
s1p=starparameters(R1,T1,0.,1,0.,0.,.25,0.,180.0-i_orb,Om-180.,0.,P)
s2p=starparameters(R2,T2,0.,1,0.,0.,.25,0.,180.0-i_orb,Om-180.,0.,P)
bp=binaryparameters(s1p,s2p,77.,i_orb,Om,om,P,a,e,T0,q,[1.,1.],0.,0.)
t = T0 + 0.31P

st1, st2 = create_binary_geometry(tessellation_healpix(3), sph(R1,T1),
                                  tessellation_healpix(2), sph(R2,T2), bp, t)

"""Split each HEALPix quad into 2 triangles. Returns global vertex list, connectivity
(1-based), per-triangle centroid, outward unit normal and area — all in the COMMON frame
(the star's own offset already added)."""
function triangulate(star)
    off = Float64.(star.center_offsets)
    n = star.npix
    V = Array{Float64}(undef, 4n, 3)
    Tr = Array{Int}(undef, 2n, 3)
    C = Array{Float64}(undef, 2n, 3)
    N = Array{Float64}(undef, 2n, 3)
    A = Array{Float64}(undef, 2n)
    for i in 1:n
        for k in 1:4
            V[4(i-1)+k, :] = Float64.(star.vertices_xyz[i, k, :]) .+ off
        end
        b = 4(i-1)
        for (j, (p,r,s)) in enumerate(((1,2,3), (1,3,4)))
            ti = 2(i-1)+j
            Tr[ti, :] = [b+p, b+r, b+s]
            v1 = V[b+p,:]; v2 = V[b+r,:]; v3 = V[b+s,:]
            C[ti, :] = (v1 .+ v2 .+ v3) ./ 3
            cr = cross(v2 .- v1, v3 .- v1)
            A[ti] = norm(cr) / 2
            nu = cr ./ norm(cr)
            # orient outward (away from the body centre)
            if dot(nu, C[ti,:] .- off) < 0; nu = -nu; end
            N[ti, :] = nu
        end
    end
    return V, Tr, C, N, A
end

V1,Tr1,C1,N1,A1 = triangulate(st1)
V2,Tr2,C2,N2,A2 = triangulate(st2)
@printf("triangles: star1 %d, star2 %d;  ΣA1=%.6f (4πR²=%.6f)  ΣA2=%.6f (4πR²=%.6f)\n",
        size(Tr1,1), size(Tr2,1), sum(A1), 4pi*R1^2, sum(A2), 4pi*R2^2)

const SB = ROTIR.SIGMA_SB
F0_1 = fill(SB*T1^4, size(C1,1))
F0_2 = fill(SB*T2^4, size(C2,1))
alb = 0.6

for (tag, ld) in (("uniform", (ldtype=0,)), ("linear", (ldtype=1, ld1=0.6, ld2=0.0)))
    G, L12, L21 = crossbody_kernels(C1, N1, ld, C2, N2, ld)
    for meth in (:horvat, :wilson)
        o1, o2, _, _, nit = solve_radiosity(G, L12, L21, A1, A2, alb, alb, F0_1, F0_2;
                                            method=meth, epsF=1e-12, maxiter=100)
        writedlm(joinpath(OUT, "jl_$(tag)_$(meth)_1.txt"), o1)
        writedlm(joinpath(OUT, "jl_$(tag)_$(meth)_2.txt"), o2)
        @printf("  %-8s %-7s: %2d sweeps, max Fout/F0 = %.6f, %.6f\n",
                tag, meth, nit, maximum(o1./F0_1), maximum(o2./F0_2))
    end
end

# geometry for PHOEBE (0-based connectivity)
for (k,(V,Tr,N,A,F0)) in enumerate(((V1,Tr1,N1,A1,F0_1), (V2,Tr2,N2,A2,F0_2)))
    writedlm(joinpath(OUT,"V$k.txt"), V);  writedlm(joinpath(OUT,"Tr$k.txt"), Tr .- 1)
    writedlm(joinpath(OUT,"N$k.txt"), N);  writedlm(joinpath(OUT,"A$k.txt"), A)
    writedlm(joinpath(OUT,"F0_$k.txt"), F0); writedlm(joinpath(OUT,"R$k.txt"), fill(alb, length(A)))
end
println("exported to $OUT")
