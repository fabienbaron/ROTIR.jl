# surface_priors.jl — priors and links for a tessellated surface, in PIXEL space.
#
# WHY PIXEL SPACE. A global angular power spectrum `C_l` is a statement about SCALES, not about
# PLACES, so it cannot express "these tessels were never observed" — and on an interferometric
# surface map some always are. Measured on zeta And: 54 of 768 tessels have |dchi2/dx| < 1e-14,
# and the latitude band -69..-45 deg (12 % of the visible area) absorbed **83 %** of a
# spherical-harmonic reconstruction's variance against 12 % of the published map's. Replacing the
# spectrum with the isotropic pixel-space prior built here took the correlation with the published
# map from 0.46 to 0.78 on the same data, geometry and sampler. See Optimization.md §26.

"""
    gp_covariance(tess, theta_c; T=Float64) -> Matrix

Isotropic Gaussian-process covariance on a tessellation's facet centres,
`C_ij = exp(-θ_ij² / 2θ_c²)` with `θ_ij` the great-circle angle between centres and `θ_c` in
DEGREES.

`θ_c` is a physical statement about feature size rather than a spectral index, which is the point:
a spot of angular radius `θ` has a power spectrum that is flat below `ℓ ≈ π/θ` and falls steeply
above it, so a measured spectral knee can be read straight off as a correlation length.
"""
function gp_covariance(tess, theta_c::Real; T = Float64)
    U = tess.unit_xyz[:, 5, :]
    n = size(U, 1)
    s = 2 * deg2rad(T(theta_c))^2
    C = Matrix{T}(undef, n, n)
    @inbounds for j in 1:n, i in 1:n
        c = clamp(T(U[i,1]*U[j,1] + U[i,2]*U[j,2] + U[i,3]*U[j,3]), -one(T), one(T))
        C[i,j] = exp(-acos(c)^2 / s)
    end
    C
end

"""
    gp_sqrt(tess, theta_c; jitter=1e-10, T=Float64) -> L

`C^{1/2}` for [`gp_covariance`](@ref), normalised so `diag(L*L') == 1` — the prior then has unit
pointwise variance and its amplitude is carried by one scalar the caller supplies, which keeps the
two knobs (how big, how smooth) from trading against each other.

One symmetric eigendecomposition; at HEALPix level 3 that is a 768×768 matrix, 0.6 MB and ~0.1 s,
so it is affordable precisely because surface meshes are small. `jitter` floors the eigenvalues
against round-off.

A field is then `σ .* (L * z)` for `z ~ N(0, I)`, which is the NON-centred form: cheap, and what
VarInf's standardised-latent interface expects. Note that giving `σ` its own latent turns this
into Neal's funnel (Optimization.md §25.8), so prefer to fix it or to centre the parameterisation.
"""
function gp_sqrt(tess, theta_c::Real; jitter::Real = 1e-10, T = Float64)
    C = gp_covariance(tess, theta_c; T = T)
    F = eigen(Symmetric(C))
    L = F.vectors * Diagonal(sqrt.(max.(F.values, T(jitter)))) * transpose(F.vectors)
    L ./ sqrt.(max.(sum(abs2, L; dims = 2), eps(T)))
end

"""
    box_link(z, tmin, tmax) -> (x, dx_dz)

Logistic link from an unconstrained latent to a bounded temperature,
`x = tmin + (tmax - tmin)·σ(z)`, with its derivative.

This is the link half of the logit-normal of §15.2. Used WITHOUT a prior on `z` it is purely a
reparameterisation imposing the box — legitimate, because a likelihood's argmax is invariant under
reparameterisation — and then the box itself is the only regularisation. Measured on zeta And at
HEALPix level 2 that is enough to reach r = 0.80 against the published map, but only because the
box binds: widen it to ±30 % about the mean and the correlation falls to 0.45, narrow it to ±5 %
and the fit starves. The range is the regulariser, not a sanity check.

Caveat: the logistic saturates, so `dx/dz → 0` at the walls and an optimiser pinned against them
makes very slow progress. With a prior `z ~ N(0, s²)` this becomes the full logit-normal.
"""
function box_link(z::AbstractVector, tmin::Real, tmax::Real)
    s = @. 1 / (1 + exp(-z))
    x = @. tmin + (tmax - tmin) * s
    dx = @. (tmax - tmin) * s * (1 - s)
    x, dx
end

"Inverse of [`box_link`](@ref): the latent that produces temperatures `x`."
box_unlink(x::AbstractVector, tmin::Real, tmax::Real) =
    (u = clamp.((x .- tmin) ./ (tmax - tmin), 1e-12, 1 - 1e-12); log.(u ./ (1 .- u)))
