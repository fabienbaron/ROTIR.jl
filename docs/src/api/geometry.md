# Geometry & surfaces

## Star creation

| Function | Description |
|----------|-------------|
| `create_star(tessels, star_params, t; secondary=false, T=eltype(tessels), κ=50)` | `stellar_geometry` for one epoch at time `t` |
| `create_star_multiepochs(tessels, star_params, tepochs; kwargs...)` | `Vector{stellar_geometry}`, one per epoch |
| `create_binary_geometry(tes1, p1, tes2, p2, bparams, t; ...)` | Both components of a binary in one frame; see [Binary geometry](binary_geometry.md) |

## Star parameters

Parameters are a `NamedTuple`; the required fields depend on `surface_type`. The binary code
builds the same fields with `starparameters()`; see [Orbits](orbits.md).

### Common fields

| Field | Type | Description |
|-------|------|-------------|
| `surface_type` | `Int` | 0 = sphere, 1 = ellipsoid, 2 = rapid rotator, 3 = Roche |
| `tpole` | `Real` | Polar temperature (K) |
| `ldtype` | `Int` | Limb-darkening law: 0 = none, 1 = linear, 2 = quadratic, 3 = power law (Hestroffer), 4 = Claret |
| `ld1`, `ld2` | `Real` | Limb-darkening coefficients (`ld3`, `ld4` for `ldtype = 4`) |
| `inclination` | `Real` | Spin-axis inclination (degrees) |
| `position_angle` | `Real` | Spin-axis position angle (degrees) |
| `rotation_period` | `Real` | Days |
| `t0` | `Real` | Rotation phase zero, in the epochs' units (optional, default 0) |

### Sphere (`surface_type = 0`)

| Field | Description |
|-------|-------------|
| `radius` | Angular radius (mas) |

### Ellipsoid (`surface_type = 1`)

| Field | Description |
|-------|-------------|
| `radius_x`, `radius_y`, `radius_z` | Semi-axes (mas) |
| `beta` | Gravity-darkening exponent |

### Rapid rotator (`surface_type = 2`)

| Field | Description |
|-------|-------------|
| `rpole` | Polar radius (mas) |
| `frac_escapevel` | Ω/Ω_crit, from 0 to 1 |
| `beta` | Gravity-darkening exponent |
| `gravity_law` | `1` or `:vonzeipel` (default), `2` or `:elr` (Espinosa Lara & Rieutord 2011); optional |
| `B_rot` | Differential rotation coefficient |

### Roche lobe (`surface_type = 3`)

| Field | Description |
|-------|-------------|
| `rpole` | Polar radius (mas) |
| `beta` | Gravity-darkening exponent |
| `q` | Mass ratio M_companion/M_self: M₂/M₁ on the primary, M₁/M₂ on the secondary |
| `fillout_factor_primary`, `fillout_factor_secondary` | Fill-out factor; `-1` uses `rpole` |
| `P` | Orbital period (days) |
| `a` | Semi-major axis (mas) |
| `e` | Eccentricity |
| `T0` | Time of periastron, in the epochs' time system |
| `dP` | Period derivative (days/day); required by `compute_separation` |
| `i`, `Ω`, `ω` | Orbital inclination, ascending node, argument of periapsis (degrees); used by `binary_frame` |
| `d` | Distance (pc), for physical-unit conversions |

`create_star` orients a Roche component with `inclination`, `position_angle` and
`rotation_period`, like a single star; `create_binary_star` orients it from the orbit.
`Ω` and `ω` are Unicode field names.

The Roche potentials, gravities, volumes and fill-out conversions (`compute_potential_*`,
`compute_gravity_*`, `roche_volume`, `fillout_to_rpole`, `rpole_to_fillout`) take the binary
mass ratio M₂/M₁ for both components, with `secondary = true` for the secondary.
`roche_q(params; secondary)` converts a component's `q` to it; `update_roche_radii`, `rl1`,
`max_rpole`, `temperature_map_vonZeipel_roche_single`, `roche_mesh_volume` and
`roche_omega_for_volume` apply it themselves.

## Temperature maps

| Function | Description |
|----------|-------------|
| `parametric_temperature_map(star_params, star)` | Map for the star's `surface_type`: uniform `tpole` (0), ellipsoid (1), `temperature_map_rapid_rotator` (2), Roche (3) |
| `temperature_map_vonZeipel_ellipsoid(params, star)` | `T = tpole·(radius_x/r)^(2β)` |
| `temperature_map_rapid_rotator(params, star)` | Von Zeipel or ELR, from `gravity_law` |
| `temperature_map_vonZeipel_rapid_rotator(params, star)` | Von Zeipel on the rotating Roche model |
| `temperature_map_ELR_rapid_rotator(params, star)` | Espinosa Lara & Rieutord (2011) |
| `temperature_map_vonZeipel_roche_single(params, star, t; secondary)` | Von Zeipel with the Roche-potential gravity |

## Rapid rotator helpers

| Function | Description |
|----------|-------------|
| `oblate_const(star_params)` | Oblate-spheroid approximation, `(a, b, c) = (R_eq, R_eq, rpole)` |
| `calc_omega(rpole, oblateness)` | `(ω, rpole, R_eq)` from `oblateness = R_eq/rpole − 1` |
| `calc_rotspin(rpole, R_equ, omega, Mass)` | `(rotational_vel, rotation_period)`; radii in R☉, mass in M☉, period in days |

## Roche lobe functions

Radii are in units of the semi-major axis `a` unless noted.

| Function | Description |
|----------|-------------|
| `update_roche_radii(tessels, params, D; use_fillout_factor, secondary, omega)` | Radii (mas) at all vertices from the surface potential |
| `get_surface_potential(rpole_a, D, q, async_ratio, fillout_factor; secondary, use_fillout_factor)` | `(potential, initial radius)` from `rpole/a` or from the fill-out factor |
| `compute_potential_primary(r, D, θ, ϕ, q, async_ratio)` | Roche potential and its first two radial derivatives, primary |
| `compute_potential_secondary(r, D, θ, ϕ, q, async_ratio)` | The same, secondary |
| `solve_radius(r0, pot, D, θ, ϕ, q, async_ratio, potential_function)` | `r` where the potential equals `pot` (Halley's method) |
| `solve_R_L1(r0, D, q, async_ratio, potential_function; secondary)` | Distance to the L1 point (Newton) |
| `compute_gravity_primary(r, θ, ϕ, D, q, async_ratio)` | Magnitude of the local gravity, primary |
| `compute_gravity_secondary(r, θ, ϕ, D, q, async_ratio)` | The same, secondary |
| `radius_eggleton(q)` | Eggleton (1983) Roche-lobe radius |
| `radius_leahy(q)` | Leahy & Leahy (2015) approximation of the L1 radius |
| `ROTIR.fillout_to_rpole(fillout, D, q, async_ratio)` | Fill-out factor → `rpole/a` |
| `rpole_to_fillout(rpole, D, q, async_ratio; secondary)` | `rpole/a` → fill-out factor |
| `max_rpole(D, params)` | Polar radius (mas) at fill-out 1 |

`async_ratio` is the synchronicity `F = P/rotation_period`; see [`synchronicity`](binary_geometry.md#Synchronicity).

## Rotation

| Function | Description |
|----------|-------------|
| `ROTIR.rotate_star(xyz, star_params, t)` | Rotates vertices by phase `2π(t − t0)/rotation_period`, `inclination` and `position_angle` |
| `compute_separation(bparams, tepoch)` | Binary separation in units of `a` |

## Visibility analysis

| Function | Description |
|----------|-------------|
| `sometimes_visible(stars)` | Pixels visible in at least one epoch |
| `never_visible(stars)` | Pixels never visible |
| `invisible_neighbors(n, stars)` | Never-visible pixels adjacent to a sometimes-visible one (HEALPix level `n`) |
| `with_invisible_neighbors(n, stars)` | Sometimes-visible pixels adjacent to an `invisible_neighbors` pixel |
| `without_invisible_neighbors(n, stars)` | Pixels with at least one neighbour in `invisible_neighbors` |
