# Binary geometry: two stars in one frame

`create_star` builds a star at its own origin and orients it with
`rotate_star(2πt/P_rot, inclination, position_angle)` — a spin about a fixed sky axis.
That is right for a single rotating star, but it leaves nothing tying a Roche lobe's
tidal bulge to the line of centres, so the bulge does not follow the companion around
the orbit and there is no way to ask "where is the companion, relative to this surface
element?" — which is exactly what mutual irradiation needs.

`src/binary_geometry.jl` derives one body→sky rotation from the orbit itself and applies
it to both components.

## Frames

| Frame | Axes |
|---|---|
| ROTIR sky | (1) West, (2) North, (3) toward the observer — **right-handed** |
| `src/orbits.jl` | x = North, y = East, z = +receding |
| Binary co-rotating | x̂ = primary→secondary, ẑ = orbital angular momentum, ŷ = ẑ × x̂ |

`sky_of_orbit(x, y, z) = (-y, x, -z)` maps the second to the first. It is
`orbit_to_rotir_offset` plus the line-of-sight component.

East-left / North-up on the sky makes `(East, North, toward-obs)` *left*-handed, which is
why ROTIR's internal frame uses West.

### Both lobes share one rotation

Reading the potentials in `src/geometry_rochelobe.jl`:

* `compute_potential_primary` puts the companion at `+D x̂`;
* `compute_potential_secondary` puts the companion at `−D x̂` (confirmed by `solve_R_L1`
  using `θ = −π/2` for the secondary).

So the primary's body `+x̂` points **at** the secondary and the secondary's body `+x̂`
points **away from** the primary — the same direction in space. Both also spin about the
orbital angular momentum. One 3×3 matrix therefore orients both components.

## Functions

| Function | Description |
|----------|-------------|
| `sky_of_orbit(x, y, z)` | Orbit frame → ROTIR sky frame |
| `binary_frame(bparams, tepoch_jd)` | `(R, D_mas, r_rel_sky)`; `R` columns are `x̂, ŷ, ẑ` in the sky frame, so `v_sky = R * v_body` |
| `create_binary_star(tessels, star_params, bparams, t; secondary, volume_conserving, omega)` | One component, deformed by the shared `D(t)` and oriented by `binary_frame` |
| `create_binary_geometry(tes1, p1, tes2, p2, bparams, t; ...)` | Both components; fills `center_offsets` |
| `create_binary_geometry_multiepochs(tes1, p1, tes2, p2, bparams, tepochs; ...)` | The above over a vector of epochs |
| `projected_separation(bparams, t)` | Sky-projected centre-to-centre separation (mas) |
| `check_binary_overlap(star1, star2, bparams, tepochs)` | Warn where the disks overlap |

ROTIR stores vertices as row vectors, so the right-multiply form is `xyz * transpose(R)`.

`create_binary_star` is `create_star` with exactly two substitutions — `D` from the shared
orbit instead of this component's own `t`, and `R` instead of `rot_vertex`. Everything
downstream (normals, soft visibility, projection, `ldmap`, `setup_oi!`, `polyft`) is
unchanged, because both share the tail helper `finish_star`.

### `center_offsets`

`create_binary_geometry` records the secondary's sky-frame offset (West, North,
toward-observer; mas) in the previously-unused `center_offsets` field of
`stellar_geometry`, with `star1.center_offsets = [0,0,0]`. That single field is the whole
interface used by `handle_reflection` and the animation code to place the two meshes
relative to each other. The vertex arrays themselves stay centred on their own component,
so `polyft` and `binary_phase_shift` are unaffected.

## Volume-conserving Roche shape

| Function | Description |
|----------|-------------|
| `tessel_solid_angles(tessels)` | Solid angle of each tessel on the unit sphere (normalised so `Σ = 4π`) |
| `roche_mesh_volume(tessels, star_params, D, Ω; secondary)` | `V = (1/3) Σ rᵢ³ ΔΩᵢ` on the mesh |
| `roche_omega_for_volume(tessels, star_params, D; secondary, V_target, D_ref)` | Ω giving the same volume at `D` as at `D_ref` |
| `roche_omega_table(tessels, star_params, bparams; secondary, npoints)` | Interpolating `Ω(D)` over `[1−e, 1+e]` |

`update_roche_radii` normally holds `rpole` fixed and lets Ω follow `D`, which makes the
star's *volume* breathe over an eccentric orbit. Physically the volume is the invariant on
an orbital timescale — PHOEBE re-solves Ω from a fixed equivalent-volume radius — and the
difference between the two conventions *is* the ellipsoidal variation. Pass
`volume_conserving=true`, or precompute a table and pass `omega`:

```julia
Ωtab = roche_omega_table(tes1, params1, bparams; secondary=false)
D = compute_separation(bparams, t)
star = create_binary_star(tes1, params1, bparams, t; omega=Ωtab(D))
```

`roche_mesh_volume` is deliberately the volume of the mesh ROTIR actually integrates over,
not the exact analytic volume: the same discretisation then appears on both sides of
`V(Ω) = V_target`, so the conserved quantity is exactly what the forward model sees.

Volume conservation is what every modern code does for eccentric orbits — none holds the
polar radius, `Ω`, or the fillout factor fixed. They differ only in the **reference
phase**: ROTIR and ELISa reference `D = a`, while Wilson–Devinney and PHOEBE legacy
reference periastron (`D = 1 − e`); PHOEBE 2's `requiv` is phase-independent. The choice
shifts the absolute size slightly but not the modulation. Pass `D_ref = 1 - e` to
`roche_omega_for_volume` for the WD convention.

Note that the pre-existing `demos/spica_binary_roche.jl` predates this and uses the
fixed-polar-radius convention (`fillout_factor_primary = -1`); the irradiation demos
default to `volume_conserving = true`.

## Mutual occultation is NOT modelled

`binary_cvis` sums the two components' visibilities unconditionally — the hidden part of
the far star still contributes flux. That is harmless while the disks are well separated
and wrong once they overlap. `check_binary_overlap` flags the epochs where
`ρ < R₁ + R₂`; it does not change the model.

For Spica this is not academic: with R₁ = 7.40 R☉, R₂ = 3.74 R☉ and a = 25.5 R☉ the
minimum projected separation is 0.594 mas against R₁+R₂ = 0.673 mas, i.e. a grazing
eclipse — independently consistent with the MOST photometry (Desmet et al. 2009, reported
in Aufdenberg et al. 2015).

## The tidal term for an eccentric orbit

The frame is centred on the star itself, not on the centre of mass. That origin is in
**free fall**, accelerating toward the companion at `G·M₂/D²`, so the non-inertial frame
carries a uniform pseudo-force with potential `−G·M₂·x/D²`, i.e. in units of `G·M₁/a`

```
    −q·x/D²  =  −q·r·λ/D²
```

This is purely translational and tidal — it assumes nothing about the frame's rotation
rate. The rotational term is separately the centrifugal potential about the star's *own*
spin axis, `½F²(1+q)r²(1−ν²)` with `F = ω_rot/ω_orb` referenced to the **mean** orbital
rate, and therefore carries no `D`. (Sepinsky, Willems & Kalogera 2007, ApJ 660, 1624,
derive exactly this frame for nonsynchronous eccentric binaries.)

The framing matters, because the tempting alternative — origin at the centre of mass,
linear term written as `ω²·x_cm·x` with `x_cm = D·q/(1+q)` — is **not** self-consistent
with the code: it would also force the rotational coefficient to be `ω² = (1+q)/D³`,
i.e. `D`-dependent, which is neither what is implemented nor what any reference code does.
(`(1+q)/D³` is in any case the *circular* rate at separation `D`; the true instantaneous
orbital rate is larger by `1 + e·cos ν`.)

This term was previously `−q·r·λ·D`, matching Aufdenberg et al. 2015 eqs. A18/A27/A30 —
internally self-consistent with his A1, but wrong by `D³` relative to the above. Two
independent modern codes write the correct form explicitly:

```
PHOEBE 2  phoebe/lib/gen_roche.h
    Ω = 1/ρ + q[(δ² + ρ² − 2ρλδ)^(−1/2) − ρλ/δ²] + ½F²(1+q)ρ²(1−ν²)

ELISa     elisa/binary_system/model.py, pre_calculate_for_potential_value_primary
    b = d²;  c = 2·d·cs;  d_coef = q·cs/b;  e = ½·F²·(1+q)·sin²θ
    Psi1 = 1/r + q/√(b + r² − c·r) − d_coef·r + e·r²          (cs = λ, d = D)
```

Both give `−q·r·λ/D²` for the linear term and a rotational term with **no** `D`
dependence, matching what is implemented here.

The old `·D` form over-weighted the term by `D³`, and past a threshold it overwhelms the
tidal attraction and points the bulge *away* from the companion.
For Spica (`q = 0.619`, `rpole/a = 0.290`) that threshold is `D/a = 1.031`, and every
modern eccentricity determination (0.065, 0.118, 0.123) puts apastron beyond it — so the
bulge was inverted for roughly half of every orbit:

```
  D/a     front(λ=+1)  back(λ=−1)   difference
  0.90     0.31953      0.29892      +0.0206   forward
  1.00     0.30661      0.30268      +0.0039   forward     mean-rate (old):
  1.04     0.30302      0.30413      −0.0011   INVERTED    sign-flipping and ~3× larger
  1.12     0.29734      0.30692      −0.0096   INVERTED

  0.90     0.31173      0.30502      +0.0067   forward     instantaneous (current):
  1.12     0.30318      0.30088      +0.0023   forward     smooth, always forward
```

The choice also sets the *amplitude* of the ellipsoidal variation, which is the observable
being fitted for the apsidal constant — so this is not a cosmetic difference.

Circular-orbit results are unchanged either way (`D ≡ 1`).

The rotational term keeps Wilson's convention: `async_ratio = ω_rot/ω_orb` is referenced
to the **mean** orbital rate, so it carries no `D` dependence. PHOEBE and
`RocheLobe.f90` both do this.

Three tests in `test/test_reflection.jl` lock this down:

* the tidal bulge is the longest radius at every orbital phase;
* `compute_gravity_*` equals `|∇Ω|` of `compute_potential_*` by finite differences
  (they share the centre-of-mass term, so one cannot be changed without the other);
* rescaling the potential to units of the instantaneous separation (`r = D·s`, then
  multiply by `D`) reproduces `RocheLobe.f90`'s `potential` exactly — this is what pins
  the `D`-dependence.

### Related fixes in `roche_volume` / `roche_area`

Verified against a direct mesh sum and the exact sphere limit:

* both quadratures carried a **spurious factor of 2** (the μ- and φ-symmetry factors were
  applied on top of an integrand that already included them), so volumes came out exactly
  2× too large and `R_vol` 2^(1/3) ≈ 1.26× too large;
* the Halley solve was seeded with the L1 radius, which lies *outside* an under-filled
  lobe and let the iteration converge to a spurious root past the saddle — returning a
  radius roughly twice too large near the companion. Now seeded from
  `roche_polar_radius`, the minimum radius on the equipotential and hence always inside;
* `roche_area` now includes the projection factor `1/cos γ` (with
  `cos γ = |∂Ω/∂r| / |∇Ω|`) that the original omitted. `dA = r² dΩ / cos γ`, so dropping
  it under-estimates the area of a distorted lobe.
