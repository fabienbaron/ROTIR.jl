# reflection.jl
# ---------------------------------------------------------------------------
# Mutual irradiation ("the reflection effect") between the two components of a
# detached binary, as a radiosity problem.
#
# Port of PHOEBE 2.5 (`phoebe/lib/reflection.h`, M. Horvat), specialised to the
# two-convex-body case. References: Wilson 1990 ApJ 356, 613; Horvat et al. 2019;
# Budaj 2011 AJ 141, 59; Kallrath & Milone 2009.
#
# ---------------------------------------------------------------------------
# The maths
# ---------------------------------------------------------------------------
# For a pair of surface elements i (on one star) and j (on the other), separated by
# a = r_j − r_i with s = |a|, and cosθ_i = n̂_i·â, cosθ_j = −n̂_j·â:
#
#     G[i,j]  = cosθ_i cosθ_j / (π s²)                              geometric kernel
#     F0[i←j] = A_j G[i,j]                                          Lambert view factor
#     F [i←j] = A_j G[i,j] · π D_j(cosθ_j) / D0_j                   limb-darkened
#     D0      = 2π ∫₀¹ D(μ) μ dμ
#
# The π·D/D0 factor is what conserves the emitted energy: D0 already contains the π, so
# F ≡ F0 exactly for a uniform (Lambert) law. Omitting it — as a naive `F = F0·ld(μ)`
# would — silently changes the total flux leaving each element.
#
#     Horvat:  S₀   = L_LD·F₀
#              F_in ← S₀ + L₀·diag(R)·F_in            (Lambertian re-emission)
#              F_out = F₀ + diag(R)·F_in
#     Wilson:  M    ← M₀ + diag(R)·L_LD·M             (LD at every scattering order)
#     both:    T_new = T·(F_out/F₀)^(1/4),   F₀ = σT⁴
#
# Note the albedo index: EMITTER (R_j) in Horvat, RECEIVER (R_i) in Wilson. With a
# uniform bolometric LD law the two schemes are algebraically identical — `test_reflection.jl`
# asserts that, which pins both implementations at once.
#
# ---------------------------------------------------------------------------
# Simplifications relative to PHOEBE (all deliberate)
# ---------------------------------------------------------------------------
# * Two CONVEX bodies only. A convex body cannot see itself, so there is no
#   self-irradiation block, and for exactly two disjoint convex bodies the double
#   cosθ > 0 test is provably sufficient — PHOEBE's own 2-body fast path
#   (reflection.h:1166-1230) likewise runs with no occlusion test at all. Contact /
#   common-envelope systems are out of scope.
# * Per-element (centroid) support, matching ROTIR's polygon-FT element model. PHOEBE
#   defaults to per-vertex but implements both.
# * Bolometric only. ROTIR's von Zeipel maps are already Teff, so the σT⁴ ↔ Teff round
#   trip is exactly PHOEBE's.
# * No redistribution of absorbed flux (PHOEBE 2.5 does not wire it up either:
#   `redistribution.h` exists but is never called from Python). The non-reflected
#   fraction 1 − albedo is simply lost.
# * Dense n₁×n₂ kernels rather than PHOEBE's sparse pair list: for 3072×768 a pair list
#   is ~2.4 M tuples (>100 MB) whereas the three dense blocks are 57 MB in Float64 (28 MB in
#   Float32) and every solver sweep becomes a BLAS gemv.
#
# Forward model only — no ChainRules rrules. The temperature bump is a post-processing
# step on the intrinsic (gravity-darkened) maps, applied before the flux operators.
# ---------------------------------------------------------------------------

const SIGMA_SB = 5.670374419e-8   # W m⁻² K⁻⁴
const EPS_COSTHETA = 0.00872654   # cos(89.5°) — PHOEBE's grazing-angle cutoff

"""
    tessel_centroids_areas(g::stellar_geometry) -> (centroids, areas)

Per-tessel centroid `(npix, 3)` and 3-D area `(npix,)`, in the frame of `g.vertices_xyz`
(mas, so areas are in mas²).

The centroid is the mean of the four quad corners `vertices_xyz[i, 1:4, :]`; the area is the
sum of triangles (1,2,3) and (1,3,4). This is the surface area, not the projected area
`polyflux`.
"""
function tessel_centroids_areas(g)
    T = eltype(g.vertices_xyz)
    n = g.npix
    c = Array{T}(undef, n, 3)
    A = Array{T}(undef, n)
    v = g.vertices_xyz
    @inbounds for i in 1:n
        c1 = (v[i,1,1] + v[i,2,1] + v[i,3,1] + v[i,4,1]) / 4
        c2 = (v[i,1,2] + v[i,2,2] + v[i,3,2] + v[i,4,2]) / 4
        c3 = (v[i,1,3] + v[i,2,3] + v[i,3,3] + v[i,4,3]) / 4
        c[i,1] = c1; c[i,2] = c2; c[i,3] = c3
        e1 = (v[i,2,1]-v[i,1,1], v[i,2,2]-v[i,1,2], v[i,2,3]-v[i,1,3])
        e2 = (v[i,3,1]-v[i,1,1], v[i,3,2]-v[i,1,2], v[i,3,3]-v[i,1,3])
        e3 = (v[i,4,1]-v[i,1,1], v[i,4,2]-v[i,1,2], v[i,4,3]-v[i,1,3])
        A[i] = (_cross_norm(e1, e2) + _cross_norm(e2, e3)) / 2
    end
    return c, A
end

@inline function _cross_norm(u, w)
    x = u[2]*w[3] - u[3]*w[2]
    y = u[3]*w[1] - u[1]*w[3]
    z = u[1]*w[2] - u[2]*w[1]
    return sqrt(x*x + y*y + z*z)
end

"""
    ld_bol_D0(ldtype, ld1, ld2) -> D0/π

Hemispheric normalisation `D0 = 2π ∫₀¹ D(μ) μ dμ` of a bolometric limb-darkening law,
divided by π (a uniform law gives 1).

| ldtype | law | D0/π |
|---|---|---|
| 0 | uniform / Lambert, `D ≡ 1` | 1 |
| 1 | linear, `1 − x(1−μ)` | `1 − x/3` |
| 2 | quadratic, `1 − x(1−μ) − y(1−μ)²` | `1 − x/3 − y/6` |
| 3 | Hestroffer power, `μ^x` | `2/(x+2)` |

These are bolometric coefficients for the irradiation, separate from the passband
`ld1`/`ld2` that set `ldmap`. Other `ldtype` values raise an error.
"""
function ld_bol_D0(ldtype::Integer, ld1::Real, ld2::Real)
    if ldtype == 0
        return 1.0
    elseif ldtype == 1
        return 1.0 - ld1 / 3
    elseif ldtype == 2
        # 2∫₀¹ (1 − x(1−μ) − y(1−μ)²) μ dμ = 2(½ − x/6 − y/12) = 1 − x/3 − y/6
        return 1.0 - ld1 / 3 - ld2 / 6
    elseif ldtype == 3
        return 2.0 / (ld1 + 2)
    else
        error("ld_bol_D0: unknown bolometric ldtype $(ldtype) (use 0, 1, 2 or 3)")
    end
end

"""
    ld_bol(ldtype, ld1, ld2, μ) -> D(μ)

Un-normalised bolometric limb-darkening profile; see [`ld_bol_D0`](@ref) for the laws.
"""
@inline function ld_bol(ldtype::Integer, ld1::Real, ld2::Real, μ::Real)
    if ldtype == 0
        return one(μ)
    elseif ldtype == 1
        return 1 - ld1 * (1 - μ)
    elseif ldtype == 2                       # quadratic: 1 − a(1−μ) − b(1−μ)²
        return 1 - ld1 * (1 - μ) - ld2 * (1 - μ)^2
    else
        return μ > 0 ? μ^ld1 : zero(μ)
    end
end

_ldspec(p) = (ldtype = Int(get(p, :ldtype, 0)),
              ld1    = Float64(get(p, :ld1, 0.0)),
              ld2    = Float64(get(p, :ld2, 0.0)))

"""
    crossbody_kernels(c1, n1, ld1, c2, n2, ld2; epsC=EPS_COSTHETA, T=Float64)
        -> (G, L12, L21)

Dense `n₁ × n₂` irradiation kernels between two convex bodies, from element centroids `c1`,
`c2` (in the same frame) and outward unit normals `n1`, `n2`. `ld1`, `ld2` are the
bolometric laws `(ldtype=, ld1=, ld2=)`.

* `G[i,j] = cosθ_i cosθ_j / (π s²)`, the Lambert kernel; zero unless both cosines exceed
  `epsC` (`cos 89.5°`, PHOEBE's grazing cutoff).
* `L12[i,j] = G[i,j] · D₂(cosθ_j)/(D0₂/π)`: star 2 emitting toward star 1.
* `L21[i,j] = G[i,j] · D₁(cosθ_i)/(D0₁/π)`: star 1 emitting toward star 2 (used transposed).

`L12 = L21 = G` for uniform laws. Emitter areas are applied in [`solve_radiosity`](@ref).

Memory is `3·n₁·n₂` elements of type `T`: 57 MB at 3072×768 in `Float64`, half in
`Float32`, which is sufficient. The loop is threaded over `j`. At 3072×768 it takes 6.8 ms
on 1 thread, 2.1 ms on 32 threads and 1.1 ms on 32 threads in `Float32` (memory-bandwidth
bound). A full frame (geometry, two temperature maps, kernels, solve) takes about 5 ms.
"""
function crossbody_kernels(c1::AbstractMatrix, n1::AbstractMatrix, ld1,
                           c2::AbstractMatrix, n2::AbstractMatrix, ld2;
                           epsC::Float64 = EPS_COSTHETA, T::Type = Float64)
    m = size(c1, 1); n = size(c2, 1)
    G   = zeros(T, m, n)
    L12 = zeros(T, m, n)
    L21 = zeros(T, m, n)
    l1 = _ldspec(ld1); l2 = _ldspec(ld2)
    # ld_bol_D0 already returns D0/π, so the emission weight π·D(μ)/D0 is D(μ)/ld_bol_D0
    # — do NOT multiply by π again (that would inflate every direct term by π and, worse,
    # break the Horvat ≡ Wilson identity, since only Horvat's higher orders use bare G).
    invD01 = 1.0 / ld_bol_D0(l1.ldtype, l1.ld1, l1.ld2)
    invD02 = 1.0 / ld_bol_D0(l2.ldtype, l2.ld1, l2.ld2)
    invπ = 1 / pi
    Threads.@threads for j in 1:n
        @inbounds begin
        c2x = Float64(c2[j,1]); c2y = Float64(c2[j,2]); c2z = Float64(c2[j,3])
        n2x = Float64(n2[j,1]); n2y = Float64(n2[j,2]); n2z = Float64(n2[j,3])
        for i in 1:m
            ax = c2x - Float64(c1[i,1])
            ay = c2y - Float64(c1[i,2])
            az = c2z - Float64(c1[i,3])
            s2 = ax*ax + ay*ay + az*az
            s2 > 0 || continue
            s = sqrt(s2)
            hi =  Float64(n1[i,1])*ax + Float64(n1[i,2])*ay + Float64(n1[i,3])*az
            hj = -(n2x*ax + n2y*ay + n2z*az)
            (hi > epsC * s && hj > epsC * s) || continue   # both must face the other
            μi = hi / s; μj = hj / s
            g = μi * μj * invπ / s2
            G[i,j]   = g
            L12[i,j] = g * ld_bol(l2.ldtype, l2.ld1, l2.ld2, μj) * invD02
            L21[i,j] = g * ld_bol(l1.ldtype, l1.ld1, l1.ld2, μi) * invD01
        end
        end # @inbounds
    end
    return G, L12, L21
end

"""
    solve_radiosity(G, L12, L21, A1, A2, R1, R2, F0_1, F0_2;
                    method=:horvat, epsF=1e-9, maxiter=100)
        -> (Fout_1, Fout_2, Fin_1, Fin_2, niter)

Solve the two-body radiosity problem by fixed-point iteration. `G`, `L12`, `L21` come from
[`crossbody_kernels`](@ref); `A1`, `A2` are 3-D element areas, `R1`, `R2` the bolometric
albedos (scalars or per-element vectors) and `F0 = σT⁴` the intrinsic exitances.

`method = :horvat` (PHOEBE's default) applies the emitter's limb darkening to the direct
irradiation and re-emits the reflected light as a Lambert source; `:wilson` uses the
limb-darkened kernel at every order. Each sweep is four matrix–vector products; convergence
takes 3–10 sweeps at typical albedos.

`epsF` is a relative L∞ tolerance (PHOEBE uses 1e-12).
"""
function solve_radiosity(G::AbstractMatrix, L12::AbstractMatrix, L21::AbstractMatrix,
                         A1::AbstractVector, A2::AbstractVector,
                         R1, R2,
                         F0_1::AbstractVector, F0_2::AbstractVector;
                         method::Symbol = :horvat, epsF::Float64 = 1e-9,
                         maxiter::Int = 100)
    method in (:horvat, :wilson) ||
        error("solve_radiosity: unknown method $(method) (use :horvat or :wilson)")
    m, n = size(G)
    (length(A1) == m && length(A2) == n) ||
        error("solve_radiosity: area vectors do not match the kernel size")
    r1 = R1 isa AbstractVector ? Float64.(R1) : fill(Float64(R1), m)
    r2 = R2 isa AbstractVector ? Float64.(R2) : fill(Float64(R2), n)

    f1 = Float64.(F0_1); f2 = Float64.(F0_2)
    a1 = Float64.(A1);   a2 = Float64.(A2)

    if method === :horvat
        # S0 = L_LD F0 ; F_in = S0 + L0 diag(R) F_in ; F_out = F0 + diag(R) F_in
        S0_1 = L12  * (a2 .* f2)
        S0_2 = L21' * (a1 .* f1)
        in1 = copy(S0_1); in2 = copy(S0_2)
        niter = 0
        for it in 1:maxiter
            niter = it
            new1 = S0_1 .+ G  * (a2 .* r2 .* in2)
            new2 = S0_2 .+ G' * (a1 .* r1 .* in1)
            d = max(_linf(new1, in1), _linf(new2, in2))
            scale = max(maximum(abs, new1; init=0.0), maximum(abs, new2; init=0.0), eps())
            in1, in2 = new1, new2
            d <= epsF * scale && break
        end
        return f1 .+ r1 .* in1, f2 .+ r2 .* in2, in1, in2, niter
    else
        # (1 − diag(R) L_LD) M = M0, albedo at the RECEIVER
        m1 = copy(f1); m2 = copy(f2)
        niter = 0
        for it in 1:maxiter
            niter = it
            new1 = f1 .+ r1 .* (L12  * (a2 .* m2))
            new2 = f2 .+ r2 .* (L21' * (a1 .* m1))
            d = max(_linf(new1, m1), _linf(new2, m2))
            scale = max(maximum(abs, new1; init=0.0), maximum(abs, new2; init=0.0), eps())
            m1, m2 = new1, new2
            d <= epsF * scale && break
        end
        # Report the reflected part so callers can check energy the same way as Horvat.
        return m1, m2, (m1 .- f1) ./ max.(r1, eps()), (m2 .- f2) ./ max.(r2, eps()), niter
    end
end

@inline _linf(a, b) = isempty(a) ? 0.0 : maximum(abs, a .- b)

"""
    reflection_kernels(star1, star2; ldbol1=(ldtype=0,), ldbol2=(ldtype=0,),
                       kernel_eltype=Float64) -> (G, L12, L21, A1, A2)

View-factor kernels and 3-D element areas for two components placed by `center_offsets`.

They depend only on the geometry and the bolometric laws. Pass the result to
[`handle_reflection`](@ref) as `kernels` to reuse it across an albedo scan. Raises an error
if both components have the same centre.
"""
function reflection_kernels(star1, star2; ldbol1 = (ldtype = 0,), ldbol2 = (ldtype = 0,),
                            kernel_eltype::Type = Float64)
    off1 = Float64.(star1.center_offsets)
    off2 = Float64.(star2.center_offsets)
    if !(norm(off2 .- off1) > 0)
        error("reflection_kernels: the two components sit at the same centre. " *
              "Build them with create_binary_geometry, which records the orbital offset " *
              "in `center_offsets`.")
    end
    c1, A1 = tessel_centroids_areas(star1)
    c2, A2 = tessel_centroids_areas(star2)
    c1 = Float64.(c1) .+ reshape(off1, 1, 3)
    c2 = Float64.(c2) .+ reshape(off2, 1, 3)
    G, L12, L21 = crossbody_kernels(c1, star1.normals, ldbol1, c2, star2.normals, ldbol2;
                                    T=kernel_eltype)
    return G, L12, L21, Float64.(A1), Float64.(A2)
end

"""
    handle_reflection(star1, tmap1, star2, tmap2;
                      albedo1=0.6, albedo2=0.6,
                      ldbol1=(ldtype=0,), ldbol2=(ldtype=0,),
                      method=:horvat, epsF=1e-9, maxiter=100, verbose=false,
                      kernel_eltype=Float64, kernels=nothing)
        -> (tmap1_heated, tmap2_heated)

Heat both components' intrinsic (gravity-darkened) temperature maps by their mutual
irradiation, as PHOEBE's `System.handle_reflection` does.

The meshes are placed using `center_offsets`, which [`create_binary_geometry`](@ref) sets.
Stars built with `create_star` (both offsets `[0,0,0]`) raise an error.

- `albedo1`, `albedo2`: bolometric albedos (PHOEBE's `irrad_frac_refl_bol`); the rest of the
  absorbed flux is lost.
- `ldbol1`, `ldbol2`: bolometric limb-darkening law of each irradiator, `(ldtype=, ld1=, ld2=)`;
  `ldtype = 0` is Lambert (PHOEBE's fallback). Separate from the passband `ld1`/`ld2`.
- `method`, `epsF`, `maxiter`: see [`solve_radiosity`](@ref). Warns if `maxiter` is reached.
- `kernel_eltype`: element type of the kernels.
- `kernels`: output of [`reflection_kernels`](@ref), to reuse across an albedo scan.

Returns new arrays; the inputs are not modified.
"""
function handle_reflection(star1, tmap1, star2, tmap2;
                           albedo1 = 0.6, albedo2 = 0.6,
                           ldbol1 = (ldtype = 0,), ldbol2 = (ldtype = 0,),
                           method::Symbol = :horvat, epsF::Float64 = 1e-9,
                           maxiter::Int = 100, verbose::Bool = false,
                           kernel_eltype::Type = Float64, kernels = nothing)
    G, L12, L21, A1, A2 = kernels === nothing ?
        reflection_kernels(star1, star2; ldbol1=ldbol1, ldbol2=ldbol2,
                           kernel_eltype=kernel_eltype) : kernels

    F0_1 = SIGMA_SB .* Float64.(tmap1) .^ 4
    F0_2 = SIGMA_SB .* Float64.(tmap2) .^ 4

    Fout1, Fout2, _, _, niter = solve_radiosity(G, L12, L21, A1, A2, albedo1, albedo2,
                                                F0_1, F0_2;
                                                method=method, epsF=epsF, maxiter=maxiter)
    if verbose
        @info "handle_reflection ($(method)): converged in $niter sweeps; " *
              "max ΔT = ($(round(maximum(Float64.(tmap1) .* ((Fout1 ./ F0_1) .^ 0.25 .- 1)), digits=1)), " *
              "$(round(maximum(Float64.(tmap2) .* ((Fout2 ./ F0_2) .^ 0.25 .- 1)), digits=1))) K"
    end
    if niter == maxiter
        @warn "handle_reflection: radiosity iteration hit maxiter = $maxiter without " *
              "reaching epsF = $epsF; the heated maps may be under-converged."
    end

    T = eltype(tmap1)
    return T.(Float64.(tmap1) .* (Fout1 ./ F0_1) .^ 0.25),
           eltype(tmap2).(Float64.(tmap2) .* (Fout2 ./ F0_2) .^ 0.25)
end
