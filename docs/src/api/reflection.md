# Mutual irradiation (the reflection effect)

`src/reflection.jl` computes how each component of a close binary heats the facing side of
the other, as a radiosity problem. It follows PHOEBE 2.5 (`phoebe/lib/reflection.h`,
M. Horvat), restricted to two detached convex bodies.

It is a forward-model step: it takes both meshes and their intrinsic (gravity-darkened)
temperature maps and returns heated maps. It is not part of the imaging optimiser and has no
`ChainRulesCore` rules.

## The model

For a pair of surface elements `i` (on one star) and `j` (on the other), separated by
`a = r_j − r_i` with `s = |a|`, `cosθ_i = n̂_i·â`, `cosθ_j = −n̂_j·â`:

```
G[i,j]  = cosθ_i cosθ_j / (π s²)                         geometric (Lambert) kernel
F0[i←j] = A_j G[i,j]                                     Lambert view factor
F [i←j] = A_j G[i,j] · D_j(cosθ_j) / (D0_j/π)            limb-darkened view factor
D0/π    = 2 ∫₀¹ D(μ) μ dμ                                ⇒ F ≡ F0 for a uniform law
```

The `D/(D0/π)` factor conserves the energy emitted by each element; only its angular
distribution changes.

Two solvers, selected by `method`:

```
Horvat:  S₀    = L_LD·F₀
         F_in  ← S₀ + L₀·diag(R)·F_in          (reflected light re-emitted Lambertian)
         F_out = F₀ + diag(R)·F_in
Wilson:  M     ← M₀ + diag(R)·L_LD·M           (limb-darkened at every order)
both:    T_new = T·(F_out/F₀)^(1/4),   F₀ = σT⁴
```

The albedo is indexed at the emitter (`R_j`) in Horvat and at the receiver (`R_i`) in
Wilson. With a uniform bolometric limb-darkening law the two give the same result.

## Functions

| Function | Description |
|----------|-------------|
| `tessel_centroids_areas(g)` | Per-tessel centroid `(npix,3)` and 3-D area `(npix,)` |
| `ld_bol_D0(ldtype, ld1, ld2)` | Normalisation `D0/π = 2∫₀¹D(μ)μdμ`, closed form |
| `ld_bol(ldtype, ld1, ld2, μ)` | Un-normalised bolometric limb-darkening profile |
| `crossbody_kernels(c1, n1, ld1, c2, n2, ld2; epsC, T)` | Dense `n₁×n₂` kernels `(G, L12, L21)` |
| `reflection_kernels(star1, star2; ldbol1, ldbol2, kernel_eltype)` | Kernels and 3-D areas for a placed pair |
| `solve_radiosity(G, L12, L21, A1, A2, R1, R2, F0_1, F0_2; method, epsF, maxiter)` | Fixed-point solve |
| `handle_reflection(star1, tmap1, star2, tmap2; ...)` | Intrinsic maps → heated maps |

`handle_reflection` keywords and defaults: `albedo1 = albedo2 = 0.6` (bolometric albedo,
PHOEBE's `irrad_frac_refl_bol`), `ldbol1 = ldbol2 = (ldtype = 0,)`, `method = :horvat`
(or `:wilson`), `epsF = 1e-9` (relative L∞ tolerance), `maxiter = 100`,
`kernel_eltype = Float64`, `kernels = nothing`, `verbose = false`. It warns when `maxiter`
is reached.

### Bolometric limb darkening

`ldbol1` / `ldbol2` are NamedTuples `(ldtype=, ld1=, ld2=)` giving the bolometric law of
each irradiator. They are separate from the passband `ld1`/`ld2` that set `ldmap`.

| `ldtype` | law | `D0/π` |
|---|---|---|
| 0 | uniform (Lambert), `D ≡ 1` (default) | 1 |
| 1 | linear, `1 − x(1−μ)` | `1 − x/3` |
| 2 | quadratic, `1 − x(1−μ) − y(1−μ)²` | `1 − x/3 − y/6` |
| 3 | Hestroffer power, `μ^x` | `2/(x+2)` |

The quadratic law is the one used by PHOEBE, Claret & Bloemen, and Kervella et al. (2017)
Eq. 5, so `ld2` corresponds to the published `b`.

## Usage

Both meshes must be in a common frame. [`create_binary_geometry`](binary_geometry.md)
records the secondary's sky-frame offset in `center_offsets`, which `handle_reflection`
uses to place them. Stars built with `create_star` (both offsets `[0,0,0]`) raise an error.

```julia
star1, star2 = create_binary_geometry(tes1, params1, tes2, params2, bparams, tepoch_jd)
tmap1 = parametric_temperature_map(params1, star1)
tmap2 = parametric_temperature_map(params2, star2; secondary=true)
heated1, heated2 = handle_reflection(star1, tmap1, star2, tmap2;
                                     albedo1=0.6, albedo2=0.6, method=:horvat)
```

The kernels depend only on the geometry and the bolometric laws. For an albedo scan at a
fixed epoch, build them once:

```julia
K = reflection_kernels(star1, star2)
for A in albedos
    h1, h2 = handle_reflection(star1, tmap1, star2, tmap2;
                               albedo1=A, albedo2=A, kernels=K)
end
```

## Performance

The kernels are three dense `n₁×n₂` matrices: 57 MB at 3072×768 in `Float64`, 28 MB with
`kernel_eltype = Float32` (sufficient; the mesh discretisation error is larger). Each solver
sweep is a matrix–vector product; convergence takes 3–10 sweeps at typical albedos.
`crossbody_kernels` is threaded; at 3072×768 it takes 6.8 ms on 1 thread, 2.1 ms on 32
threads, 1.1 ms on 32 threads in `Float32`.

## Differences from PHOEBE

* Two convex bodies only: no self-irradiation and no occlusion test (PHOEBE's two-body convex
  path, `reflection.h:1166-1230`, does the same). Contact systems are not supported.
* Per-element (centroid) quantities, as in ROTIR's polygon Fourier transform. PHOEBE
  defaults to per-vertex.
* Bolometric only; ROTIR's temperature maps are already `Teff`.
* No redistribution of absorbed flux (PHOEBE 2.5 does not use it either). The fraction
  `1 − albedo` is lost.
* Mutual occultation is handled separately, by the `occultation` keyword of `binary_cvis`;
  see `check_binary_overlap` and `occultation_weights`.

## Validation

`test/test_reflection.jl` covers the analytic limits (albedo 0, `D → ∞`, the substellar
point for two spheres and its convergence with resolution), energy conservation, the `D0`
normalisation against numerical quadrature, and Horvat = Wilson for uniform limb darkening.

`test/phoebe_crosscheck/` (not part of `Pkg.test()`) runs the same triangulation through
this solver and through PHOEBE's C++ `libphoebe.mesh_radiosity_problem_nbody_convex`. On
1536 + 384 triangles at albedo 0.6, `max|Δ(reflected)|/max` is between 3×10⁻¹⁵ and
2.4×10⁻¹³ for both methods with uniform and linear bolometric limb darkening.

## Gravity-darkening exponents

PHOEBE's `gravb_bol` is 4× ROTIR's `beta`. PHOEBE writes
`T = T_pole·((g/g_pole)^gravb_bol)^0.25` (von Zeipel: `gravb_bol = 1`); ROTIR writes
`T = tpole·(g/g_pole)^beta` (von Zeipel: `beta = 0.25`).

## References

- Wilson, R. E. 1990, ApJ 356, 613
- Horvat, M. et al. 2019 (PHOEBE 2 radiosity)
- Budaj, J. 2011, AJ 141, 59
- Kallrath, J. & Milone, E. F. 2009, *Eclipsing Binary Stars: Modeling and Analysis*
- PHOEBE 2 source: `phoebe/lib/reflection.h`, `phoebe/backend/universe.py`
