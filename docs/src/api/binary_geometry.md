# Binary geometry: two stars in one frame

`create_star` orients a star with `rotate_star`: a spin at `rotation_period` about a fixed
sky axis. `src/binary_geometry.jl` instead derives one body→sky rotation from the orbit and
applies it to both components, so a Roche lobe's tidal bulge follows the companion and each
surface element's position relative to the companion is known (needed for mutual
irradiation).

## Frames

| Frame | Axes |
|---|---|
| ROTIR sky | (1) West, (2) North, (3) toward the observer; right-handed |
| `src/orbits.jl` | x = North, y = East, z = receding |
| Binary co-rotating | x̂ = primary→secondary, ẑ = orbital angular momentum, ŷ = ẑ × x̂ |

`sky_of_orbit(x, y, z) = (-y, x, -z)` maps the second frame to the first. It is
`orbit_to_rotir_offset` plus the line-of-sight component.

### Both lobes share one rotation

`compute_potential_primary` puts the companion at `+D x̂`; `compute_potential_secondary` puts
it at `−D x̂` (`solve_R_L1` uses `θ = −π/2` for the secondary). The primary's body `+x̂`
therefore points at the secondary and the secondary's body `+x̂` points away from the
primary, which is the same direction in space. Both spin about the orbital angular momentum,
so one 3×3 matrix orients both components.

## Functions

| Function | Description |
|----------|-------------|
| `sky_of_orbit(x, y, z)` | Orbit frame → ROTIR sky frame |
| `binary_frame(bparams, tepoch_jd)` | `(R, D_mas, r_rel_sky)`; the columns of `R` are `x̂, ŷ, ẑ` in the sky frame, so `v_sky = R * v_body` |
| `create_binary_star(tessels, star_params, bparams, t; secondary, volume_conserving, omega)` | One component, deformed by the shared `D(t)` and oriented by `binary_frame` |
| `create_binary_geometry(tes1, p1, tes2, p2, bparams, t; ...)` | Both components; sets `center_offsets` |
| `create_binary_geometry_multiepochs(tes1, p1, tes2, p2, bparams, tepochs; ...)` | The above over a vector of epochs |
| `projected_separation(bparams, t)` | Sky-projected centre-to-centre separation (mas) |
| `check_binary_overlap(star1, star2, bparams, tepochs)` | One `Bool` per epoch, `true` where the disks overlap; warns |

Times are JD. ROTIR stores vertices as row vectors, so the rotation is applied as
`xyz * transpose(R)`.

`create_binary_star` differs from `create_star` in two places: `D` comes from the shared
orbit, and the orientation is `R` instead of the spin rotation. The rest (normals, soft
visibility, projection, `ldmap`, `setup_oi!`, `polyft`) is the shared `finish_star`.

### `center_offsets`

`create_binary_geometry` stores the secondary's sky-frame offset (West, North,
toward observer; mas) in `star2.center_offsets`, and sets `star1.center_offsets = [0,0,0]`.
`handle_reflection` and the animation code use this field to place the two meshes. The
vertex arrays stay centred on their own component, so `polyft` and `binary_phase_shift` are
unaffected.

## Volume-conserving Roche shape

| Function | Description |
|----------|-------------|
| `tessel_solid_angles(tessels)` | Solid angle of each tessel on the unit sphere, normalised so `Σ = 4π` |
| `roche_mesh_volume(tessels, star_params, D, Ω; secondary)` | `V = (1/3) Σ rᵢ³ ΔΩᵢ` on the mesh, in units of a³ |
| `roche_omega_for_volume(tessels, star_params, D; secondary, V_target, D_ref)` | Ω giving the same volume at `D` as at `D_ref` |
| `roche_omega_table(tessels, star_params, bparams; secondary, npoints)` | Linear interpolant `Ω(D)` over `[1−e, 1+e]` |

By default `update_roche_radii` holds `rpole` fixed and lets Ω follow `D`, so the volume
changes over an eccentric orbit. `volume_conserving = true` holds the volume fixed instead
(PHOEBE re-solves Ω from a fixed equivalent-volume radius). The difference between the two
is the ellipsoidal variation. Pass `volume_conserving = true`, or precompute a table and
pass `omega`:

```julia
Ωtab = roche_omega_table(tes1, params1, bparams; secondary=false)
D = compute_separation(bparams, t)
star = create_binary_star(tes1, params1, bparams, t; omega=Ωtab(D))
```

`roche_mesh_volume` is the volume of the mesh ROTIR integrates over, not the analytic
volume, so both sides of `V(Ω) = V_target` use the same discretisation.

Reference phase: ROTIR and ELISa conserve the volume the star has at `D = a` (`D_ref = 1`);
Wilson–Devinney and PHOEBE legacy use periastron (`D = 1 − e`); PHOEBE 2's `requiv` is
phase-independent. The choice changes the absolute size slightly, not the modulation. Pass
`D_ref = 1 - e` to `roche_omega_for_volume` for the WD convention.

`demos/spica_binary_roche.jl` holds the polar radius fixed (`fillout_factor_primary = -1`);
`demos/spica_irradiation_movie.jl` and `demos/spica_proximity_movie.jl` use
`volume_conserving = true`.

## Mutual occultation

| Function | Description |
|----------|-------------|
| `occultation_weights(star1, star2; method, softness, nbins)` | Per-tessel occultation factors `(w1, w2, occulting)` |
| `convex_hull_2d(xs, ys)` | Counter-clockwise convex hull (Andrew's monotone chain) |
| `polygon_convex_clip_area(...)` | Area of a polygon clipped by a convex polygon (in `rasterize.jl`) |
| `projected_limb_profile(star; nbins=72)` | Maximum projected radius of a component's outline, per position-angle bin |
| `limb_radius(profile, φ)` | Limb radius at position angle `φ`, interpolated |
| `silhouette_polygon(star; nbins=72)` | `(sx, sy, profile)`: that profile as a closed counter-clockwise polygon |
| `check_binary_overlap(...)` | Which epochs overlap |

`binary_cvis` sums the two components without occultation by default. Pass
`occultation = true` (or `:exact`, `:soft`) to apply it, or a precomputed `(w1, w2)` to
reuse it across a scan:

```julia
v2, t3a, t3p = binary_observables(tm1, star1, tm2, star2, data, phase; occultation=true)
```

The components are treated as disjoint convex bodies, as in the radiosity solver: the one
whose centre is nearer the observer is in front everywhere, and only the far star is
masked. A tessel of the far star is hidden where it falls inside the near star's projected
outline, which costs O(n₁+n₂).

Two methods:

* `method = :exact` (default): the fraction of each far tessel's projected quad outside the
  near component's convex hull, by Sutherland–Hodgman clipping (`convex_hull_2d` →
  `polygon_convex_clip_area`). No free parameters. The clip requires a convex polygon; the
  projection of a convex body is convex.
* `method = :soft`: the tessel centre is tested against the azimuthal
  `projected_limb_profile`, smoothed by `wᵢ = σ((dᵢ − R_limb(φᵢ))/(width·R_limb(φᵢ)))`, with
  `width = softness` (default `√(4π/npix)/2`). Needs no convexity and is C∞, for use
  under differentiation.

Occulted fraction of the primary's flux at Spica's closest approach (ν + ω ≈ 270°,
secondary in front), Roche lobes, converged value 2.33 %:

| nside₁/₂ | npix₂ | `:soft` | `:exact` |
|---|---|---|---|
| 2/1 | 48 | 3.59 (+54 %) | 1.83 (−21 %) |
| 3/2 | 192 | 2.08 (−11 %) | 2.18 (−6 %) |
| 4/3 | 768 | 2.25 (−3.2 %) | 2.29 (−1.4 %) |
| 5/4 | 3072 | 2.36 (+1.3 %) | 2.32 (−0.3 %) |
| 6/5 | 12288 | 2.33 | 2.33 |

At nside 4/3, `:exact` takes 1.15 ms and `:soft` 0.22 ms. Only tessels straddling the limb
are clipped; the others are classified by a circle test against the hull's inscribed and
circumscribed radii. A bounding-circle test on the two components returns
`(nothing, nothing, 0)` without allocating when the disks do not overlap.

For Spica (R₁ = 7.40 R☉, R₂ = 3.74 R☉, a = 25.5 R☉) the minimum projected separation is
0.594 mas against R₁+R₂ = 0.673 mas: a grazing eclipse, consistent with the MOST photometry
(Desmet et al. 2009, reported in Aufdenberg et al. 2015). At closest approach the secondary
hides about 2.3 % of the primary's projected flux, 1.9 % of the system total.

## Apsidal motion and period change

`omega_at(bparams, tepoch)` returns the argument of periapsis at `tepoch`,
`ω(t) = ω₀ + dω·(t − T0)`, with `dω` in degrees/day (positive in the direction of orbital
motion). `binary_orbit_rel`, `binary_orbit_abs`, `binary_RV`, `binary_proj_plane` and
`binary_frame` use it. For Spica (U ≈ 139 yr, dω/dt ≈ 2.6°/yr) ω advances about 21° across
the 2007–2015 CHARA observations, which moves the predicted secondary position by up to
0.44 mas (MIRC precision: 0.01–0.05 mas).

`compute_eccentric_anomaly` applies `dP` (Ṗ, days/day; positive = lengthening) through the
quadratic ephemeris `t_n = T0 + P·n + ½·Ṗ·P·n²`, for either sign of Ṗ. At β Lyrae's rate,
neglecting Ṗ is worth 65° of orbital phase after 10⁴ d. If the extrapolated period would
reach zero within the epoch range, it warns and uses a constant period.

## Synchronicity

`synchronicity(p) = p.P / p.rotation_period` is `F = ω_rot/ω_orb`, the ratio of angular
rates. It is PHOEBE's `syncpar`, ELISa's `synchronicity` and `P` in `RocheLobe.f90`, and enters
the centrifugal term as `½F²(1+q)r²(1−ν²)`.

For Spica's primary (v sin i = 161 km/s, F ≈ 1.92) the reciprocal `rotation_period/P` would
scale the centrifugal term by F⁴ ≈ 13.6 and give 0.5 % rotational flattening instead of
9.4 %. Checked against `libphoebe.roche_Omega(q, F, d, [x,y,z])` at F = 0.521, 1.0 and 1.92.

The Spica demos set `rotation_period = P_orb` (synchronous). Spica's measured rotation
corresponds to `rotation_period = P/1.92` for the primary and `P/1.65` for the secondary.

## The tidal term for an eccentric orbit

The frame is centred on the star. That origin accelerates toward the companion at
`G·M₂/D²`, so the frame carries a uniform pseudo-force with potential `−G·M₂·x/D²`; in
units of `G·M₁/a`:

```
    −q·x/D²  =  −q·r·λ/D²
```

This term is translational and assumes nothing about the frame's rotation. The rotational
term is the centrifugal potential about the star's own spin axis,
`½F²(1+q)r²(1−ν²)`, with `F` referenced to the mean orbital rate, so it has no `D`
dependence. Sepinsky, Willems & Kalogera (2007, ApJ 660, 1624) derive this frame for
nonsynchronous eccentric binaries.

A centre-of-mass origin, with the linear term written `ω²·x_cm·x` and `x_cm = D·q/(1+q)`,
would also require the rotational coefficient `ω² = (1+q)/D³`, which is `D`-dependent;
neither ROTIR nor the reference codes do this. (`(1+q)/D³` is the circular rate at
separation `D`; the instantaneous orbital rate is larger by `1 + e·cos ν`.)

Aufdenberg et al. 2015 (eqs. A18/A27/A30) give `−q·r·λ·D`, which differs by `D³`. PHOEBE 2
and ELISa both use `−q·r·λ/D²` and a `D`-independent rotational term:

```
PHOEBE 2  phoebe/lib/gen_roche.h
    Ω = 1/ρ + q[(δ² + ρ² − 2ρλδ)^(−1/2) − ρλ/δ²] + ½F²(1+q)ρ²(1−ν²)

ELISa     elisa/binary_system/model.py, pre_calculate_for_potential_value_primary
    b = d²;  c = 2·d·cs;  d_coef = q·cs/b;  e = ½·F²·(1+q)·sin²θ
    Psi1 = 1/r + q/√(b + r² − c·r) − d_coef·r + e·r²          (cs = λ, d = D)
```

With the `·D` form, the term exceeds the tidal attraction beyond a threshold separation and
points the bulge away from the companion. For Spica (`q = 0.619`, `rpole/a = 0.290`) the
threshold is `D/a = 1.031`; the published eccentricities (0.065, 0.118, 0.123) all put
apastron beyond it. Radii toward (λ = +1) and away from (λ = −1) the companion:

```
  D/a     front(λ=+1)  back(λ=−1)   difference
                                                 −q·r·λ·D:
  0.90     0.31953      0.29892      +0.0206   toward
  1.00     0.30661      0.30268      +0.0039   toward
  1.04     0.30302      0.30413      −0.0011   away
  1.12     0.29734      0.30692      −0.0096   away
                                                 −q·r·λ/D² (ROTIR):
  0.90     0.31173      0.30502      +0.0067   toward
  1.12     0.30318      0.30088      +0.0023   toward
```

The form also sets the amplitude of the ellipsoidal variation, which is the observable
fitted for the apsidal constant. For a circular orbit (`D ≡ 1`) the two forms agree.

Tests in `test/test_reflection.jl`:

* "tidal bulge tracks the companion": the bulge is the longest radius at every orbital
  phase;
* "Roche potential: gravity is the gradient of the potential": `compute_gravity_*` equals
  `|∇Ω|` of `compute_potential_*` by finite differences;
* "Roche potential reduces to the D ≡ 1 reference": rescaling to units of the instantaneous
  separation (`r = D·s`, then multiplying by `D`) reproduces `RocheLobe.f90`'s `potential`.

### `roche_volume` and `roche_area`

Both integrate over the equipotential in (μ, φ) by Romberg quadrature, using the μ- and
φ-symmetries, and are tested against a direct mesh sum and the sphere limit.

* The radius in each direction is a Halley solve seeded from `roche_polar_radius`, the
  minimum radius on the equipotential, so the iteration converges to the lobe surface and
  not to a root beyond the L1 saddle.
* `roche_area` includes the projection factor `dA = r² dΩ / cos γ`, with
  `cos γ = |∂Ω/∂r| / |∇Ω|`.
