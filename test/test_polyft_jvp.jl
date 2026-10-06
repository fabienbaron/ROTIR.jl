# Verify the forward-mode JVP of the exact polygon Fourier transform three ways.
#
# The decisive one is the ADJOINT IDENTITY: <J v, w> = <v, J' w> for random tangents v and
# cotangents w, using the EXISTING, independently finite-difference-checked VJP. It ties the two
# directions together, so a sign or index error in either shows up — which finite differences on
# the JVP alone would not catch if the forward model itself were misread.
using ROTIR, LinearAlgebra, Printf, Random, Test

Random.seed!(20260205)
npix, nuv = 60, 90
T = Float64
pw = randn(T, npix, 4) .* 0.4
pn = randn(T, npix, 4) .* 0.4
xw = abs.(randn(T, npix)) .+ T(0.1)
kx = randn(T, nuv) .* 2.0
ky = randn(T, nuv) .* 2.0
k2 = ROTIR.precompute_k2_inv_im(kx, ky)

"the forward transform this JVP differentiates: F[k] = sum_p xw[p] * sum_e E(p,e,k)"
function fwd(pw, pn, xw)
    F = zeros(Complex{T}, nuv)
    for k in 1:nuv, p in 1:npix, e in 1:4
        j1 = e; j2 = mod1(e+1, 4)
        dx = pw[p,j2]-pw[p,j1]; dy = pn[p,j2]-pn[p,j1]
        cx = pw[p,j2]+pw[p,j1]; cy = pn[p,j2]+pn[p,j1]
        kdd = kx[k]*dx + ky[k]*dy; kdc = kx[k]*cx + ky[k]*cy
        cr  = ky[k]*dx - kx[k]*dy
        F[k] += xw[p] * k2[k] * sinc(kdd) * cis(-π*kdc) * cr
    end
    F
end

nfail = Ref(0)
chk(lab, got, want; tol=1e-7) = begin
    rel = abs(got-want)/max(abs(want), 1e-12)
    rel > tol && (nfail[] += 1)
    @printf("  %-52s %10.3e %s\n", lab, rel, rel <= tol ? "✓" : "*** FAIL ***")
end

println("\n[1] JVP against central finite differences of the forward transform")
for (lab, dxw, dpw, dpn) in (
        ("map tangent only",      randn(T,npix),  zeros(T,npix,4), zeros(T,npix,4)),
        ("geometry tangent only", zeros(T,npix),  randn(T,npix,4), randn(T,npix,4)),
        ("both together",         randn(T,npix),  randn(T,npix,4), randn(T,npix,4)))
    Fd = zeros(Complex{T}, nuv)
    compute_jvp_cvis!(Fd, dxw, dpw, dpn, kx, ky, k2, pw, pn, xw)
    h = 1e-6
    fd = (fwd(pw .+ h.*dpw, pn .+ h.*dpn, xw .+ h.*dxw) .-
          fwd(pw .- h.*dpw, pn .- h.*dpn, xw .- h.*dxw)) ./ (2h)
    @printf("  %-52s %10.3e %s\n", lab, norm(Fd .- fd)/norm(fd),
            norm(Fd .- fd)/norm(fd) <= 1e-6 ? "✓" : "*** FAIL ***")
    norm(Fd .- fd)/norm(fd) > 1e-6 && (nfail[] += 1)
end

println("\n[2] ADJOINT IDENTITY  <J v, w> = <v, J' w>   (against the existing VJP)")
for trial in 1:3
    dxw = randn(T, npix); dpw = randn(T, npix, 4); dpn = randn(T, npix, 4)
    w   = randn(Complex{T}, nuv)
    Fd  = zeros(Complex{T}, nuv)
    compute_jvp_cvis!(Fd, dxw, dpw, dpn, kx, ky, k2, pw, pn, xw)
    # PAIRING CONVENTION: ROTIR's adjoints return real(sum(adj[k]*E[k])) with NO conjugation,
    # so the inner product they are adjoint to is Re(w' F) in the TRANSPOSE sense, not Julia's
    # `dot` (which conjugates its first argument). Using `dot` here fails by O(1) — not a bug in
    # either direction, but exactly the mismatch that would silently corrupt a VarInf wiring.
    lhs = real(sum(w .* Fd))                     # <J v, w>
    gxw = zeros(T, npix); gpw = zeros(T, npix, 4); gpn = zeros(T, npix, 4)
    ROTIR.compute_adjoint_cvis!(gxw, w, kx, ky, k2, pw, pn, xw)
    ROTIR.compute_adjoint_vertices!(gpw, gpn, w, kx, ky, k2, pw, pn, xw, zeros(T, npix))
    rhs = dot(gxw, dxw) + dot(gpw, dpw) + dot(gpn, dpn)    # <v, J' w>
    chk("trial $trial:  <Jv,w> = $(round(lhs,digits=6))", lhs, rhs; tol=1e-9)
end

println("\n[3] linearity in the tangent")
dxw = randn(T,npix); dpw = randn(T,npix,4); dpn = randn(T,npix,4)
A = zeros(Complex{T}, nuv); B = zeros(Complex{T}, nuv)
compute_jvp_cvis!(A, dxw, dpw, dpn, kx, ky, k2, pw, pn, xw)
compute_jvp_cvis!(B, 3 .* dxw, 3 .* dpw, 3 .* dpn, kx, ky, k2, pw, pn, xw)
@printf("  %-52s %10.3e %s\n", "J(3v) = 3 J(v)", norm(B .- 3 .*A)/norm(B),
        norm(B .- 3 .*A)/norm(B) <= 1e-12 ? "✓" : "*** FAIL ***")
norm(B .- 3 .*A)/norm(B) > 1e-12 && (nfail[] += 1)

@printf("\n=== %s ===\n", nfail[] == 0 ? "all JVP checks passed" : "$(nfail[]) FAILURES")
nfail[] == 0 || error("JVP verification failed")
