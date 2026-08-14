# Mutual irradiation (the reflection effect)

Each component of a close binary heats the face of the other. `src/reflection.jl` solves
that as a radiosity problem, following PHOEBE 2.5 (`phoebe/lib/reflection.h`, M. Horvat)
specialised to the detached two-convex-body case.

This is a **forward-model / geometry step**: it takes both meshes and their intrinsic
(gravity-darkened) temperature maps and returns *heated* maps. It is not part of the
inverse imaging optimiser, and it defines no `ChainRulesCore` rules.

## The model

For a pair of surface elements `i` (on one star) and `j` (on the other), separated by
`a = r_j − r_i` with `s = |a|`, `cosθ_i = n̂_i·â`, `cosθ_j = −n̂_j·â`:

```
G[i,j]  = cosθ_i cosθ_j / (π s²)                         geometric (Lambert) kernel
F0[i←j] = A_j G[i,j]                                     Lambert view factor
F [i←j] = A_j G[i,j] · D_j(cosθ_j) / (D0_j/π)            limb-darkened view factor
D0/π    = 2 ∫₀¹ D(μ) μ dμ                                ⇒ F ≡ F0 for a uniform law
```

The `D/(D0/π)` factor is what conserves the emitted energy. A naive `F = F0 · ld(μ)`
would change the total flux leaving each element, not just its angular distribution.

Two solvers, selected by `method`:

```
Horvat:  S₀    = L_LD·F₀
         F_in  ← S₀ + L₀·diag(R)·F_in          (reflected light re-emitted Lambertian)
         F_out = F₀ + diag(R)·F_in
Wilson:  M     ← M₀ + diag(R)·L_LD·M           (limb-darkened at every order)
both:    T_new = T·(F_out/F₀)^(1/4),   F₀ = σT⁴
```

The albedo is indexed at the **emitter** (`R_j`) in Horvat and at the **receiver** (`R_i`)
in Wilson. With a uniform bolometric limb-darkening law the two are algebraically
identical — `test/test_reflection.jl` asserts that, which pins both implementations at
once.

## Functions

| Function | Description |
|----------|-------------|
| `tessel_centroids_areas(g)` | Per-tessel centroid `(npix,3)` and **true 3-D** area `(npix,)` |
| `ld_bol_D0(ldtype, ld1, ld2)` | Hemispheric normalisation `D0/π = 2∫₀¹D(μ)μdμ`, closed form |
| `ld_bol(ldtype, ld1, ld2, μ)` | Un-normalised bolometric LD profile |
| `crossbody_kernels(c1, n1, ld1, c2, n2, ld2; epsC, T)` | Dense `n₁×n₂` kernels `(G, L12, L21)` |
| `reflection_kernels(star1, star2; ldbol1, ldbol2, kernel_eltype)` | Kernels + 3-D areas for a placed pair |
| `solve_radiosity(G, L12, L21, A1, A2, R1, R2, F0_1, F0_2; method, epsF, maxiter)` | Fixed-point solve |
| `handle_reflection(star1, tmap1, star2, tmap2; ...)` | Top-level driver: intrinsic maps → heated maps |

### Bolometric limb darkening

`ldbol1` / `ldbol2` are NamedTuples `(ldtype=, ld1=, ld2=)` describing the *bolometric*
law of each irradiator. These are **not** the passband `ld1`/`ld2` that shape `ldmap`.

| `ldtype` | law | `D0/π` |
|---|---|---|
| 0 | uniform / Lambert, `D ≡ 1` (default) | 1 |
| 1 | linear, `1 − x(1−μ)` | `1 − x/3` |
| 2 | quadratic, `1 − x(1−μ) − y(1−μ)²` | `1 − x/3 − y/6` |
| 3 | Hestroffer power, `μ^x` | `2/(x+2)` |

ROTIR's quadratic law is the standard one, `1 − x(1−μ) − y(1−μ)²`, as used by PHOEBE,
Claret & Bloemen, and Kervella et al. (2017) Eq. 5, so `ld2` is directly comparable to the
published `b`. Note `(1−μ)²`, not `(1−μ²)`: the two differ by a factor `(1+μ)`.

## Usage

The two meshes must be in a common frame. That is what
[`create_binary_geometry`](binary_geometry.md) is for — it records the secondary's
sky-frame offset in `center_offsets`, which is the only channel `handle_reflection` uses.
Calling it on stars built by plain `create_star` (both offsets `[0,0,0]`) is an error.

```julia
star1, star2 = create_binary_geometry(tes1, params1, tes2, params2, bparams, tepoch_jd)
tmap1 = parametric_temperature_map(params1, star1)
tmap2 = parametric_temperature_map(params2, star2; secondary=true)
heated1, heated2 = handle_reflection(star1, tmap1, star2, tmap2;
                                     albedo1=0.6, albedo2=0.6, method=:horvat)
```

Scanning albedo at a fixed epoch? Build the geometry once:

```julia
K = reflection_kernels(star1, star2)
for A in albedos
    h1, h2 = handle_reflection(star1, tmap1, star2, tmap2;
                               albedo1=A, albedo2=A, kernels=K)
end
```

## Performance

The kernels are dense `n₁×n₂` matrices rather than PHOEBE's sparse pair list: at
3072×768 a pair list is ~2.4 M tuples (>100 MB) whereas three dense blocks are 57 MB in
`Float64` (or 28 MB with `kernel_eltype=Float32`, which is ample — the mesh truncation
error dominates), and every solver sweep is a BLAS `gemv`. Convergence is typically 3–10
sweeps at realistic albedos.

## Deliberate simplifications vs PHOEBE

* **Two convex bodies only.** A convex body cannot see itself, so there is no
  self-irradiation block; and for two disjoint convex bodies the double `cosθ > 0` test is
  provably sufficient, so no occlusion test is needed at all — PHOEBE's own 2-body convex
  fast path (`reflection.h:1166-1230`) does the same. Contact systems are out of scope.
* **Per-element (centroid) support**, matching ROTIR's polygon-FT element model. PHOEBE
  defaults to per-vertex but implements both.
* **Bolometric only.** ROTIR's von Zeipel maps are already `Teff`, so the `σT⁴ ↔ Teff`
  round trip is exactly PHOEBE's.
* **No redistribution** of absorbed flux — PHOEBE 2.5 does not wire it up either
  (`redistribution.h` exists but is never called from Python). The non-reflected fraction
  `1 − albedo` is simply lost.
* **No mutual occultation** in the forward model. See `check_binary_overlap`.

## Validation

`test/test_reflection.jl` covers the analytic limits (albedo 0, `D → ∞`, the substellar
point for two spheres and its convergence with mesh resolution), energy conservation, the
`D0` normalisation against numerical quadrature, and the Horvat ≡ Wilson identity.

`test/phoebe_crosscheck/` additionally runs **the same triangulation** through this solver
and through PHOEBE's own C++ `libphoebe.mesh_radiosity_problem_nbody_convex`. They agree
to machine precision — `max|Δ(reflected)|/max ≈ 10⁻¹⁴` on 1536 + 384 triangles, for both
methods and for both uniform and linear bolometric limb darkening. The non-uniform case is
the meaningful one: it is what exercises the `D0` normalisation.

## Convention note

PHOEBE's `gravb_bol` is **4×** ROTIR's `beta`: PHOEBE writes
`T = T_pole·((g/g_pole)^gravb_bol)^0.25` (so von Zeipel is `gravb_bol = 1`), ROTIR writes
`T = tpole·(g/g_pole)^beta` (von Zeipel is `beta = 0.25`). Remember this when transcribing
parameters between the two.

## References

- Wilson, R. E. 1990, ApJ 356, 613 — *Accuracy and efficiency in the binary star
  reflection effect*
- Horvat, M. et al. 2019 — the PHOEBE 2 radiosity scheme
- Budaj, J. 2011, AJ 141, 59
- Kallrath, J. & Milone, E. F. 2009, *Eclipsing Binary Stars — Modeling and Analysis*
- PHOEBE 2 source: `phoebe/lib/reflection.h`, `phoebe/backend/universe.py`
