# Geometry & surfaces

## Star creation

| Function | Description |
|----------|-------------|
| `create_star(tessels, star_params, t; secondary=false, T=Float64, kappa=50)` | Create `stellar_geometry` for one epoch at time `t` |
| `create_star_multiepochs(tessels, star_params, tepochs; kwargs...)` | Create geometries for all epochs, returns `Vector{stellar_geometry}` |
| `create_binary(star1_tessels, star2_tessels, binary_params, t)` | Create binary system geometry |

## Star parameters

Parameters are passed as a `NamedTuple`. Required fields depend on `surface_type`:

### Common fields

| Field | Type | Description |
|-------|------|-------------|
| `surface_type` | `Int` | 0=sphere, 1=ellipsoid, 2=rapid rotator, 3=Roche |
| `tpole` | `Real` | Polar temperature (K) |
| `ldtype` | `Int` | Limb-darkening law: 1=linear, 2=quadratic, 3=Hestroffer |
| `ld1`, `ld2` | `Real` | Limb-darkening coefficients |
| `inclination` | `Real` | Degrees |
| `position_angle` | `Real` | Degrees |
| `rotation_period` | `Real` | Days |

### Sphere (`surface_type = 0`)

| Field | Description |
|-------|-------------|
| `radius` | Angular radius (mas) |

### Ellipsoid (`surface_type = 1`)

| Field | Description |
|-------|-------------|
| `radius_x`, `radius_y`, `radius_z` | Semi-axes (mas) |
| `beta` | von Zeipel exponent |

### Rapid rotator (`surface_type = 2`)

| Field | Description |
|-------|-------------|
| `rpole` | Polar radius (mas) |
| `frac_escapevel` | Fractional rotational velocity omega (0 to 1) |
| `beta` | von Zeipel exponent |
| `B_rot` | Differential rotation coefficient |

### Roche lobe (`surface_type = 3`)

| Field | Description |
|-------|-------------|
| `rpole` | Polar radius (mas) |
| `beta` | von Zeipel exponent |
| `q` | Mass ratio M2/M1 |
| `d` | Distance (pc) |
| `fillout_factor_primary` | Fillout factor (negative = use rpole) |
| `i`, `Omega`, `omega` | Orbital angles (degrees) |
| `P` | Orbital period (days) |
| `a` | Semi-major axis (mas) |
| `e` | Eccentricity |
| `T0` | Time of periastron (JD) |

## Temperature maps

| Function | Description |
|----------|-------------|
| `parametric_temperature_map(star_params, star)` | Dispatcher: returns von Zeipel map for any surface type |
| `temperature_map_vonZeipel_ellipsoid(params, star)` | Von Zeipel gravity darkening for ellipsoid |
| `temperature_map_vonZeipel_rapid_rotator(params, star)` | Von Zeipel with centrifugal distortion |
| `temperature_map_vonZeipel_roche_single(params, star, t)` | Von Zeipel with Roche potential gravity |

## Rapid rotator helpers

| Function | Description |
|----------|-------------|
| `oblate_const(star_params)` | Approximate rapid rotator as oblate spheroid, returns (a, b, c) |
| `calc_omega(rpole, oblateness)` | Convert oblateness to fractional angular velocity |
| `calc_rotspin(rpole, R_equ, omega, Mass)` | Compute velocity (km/s), period (days), angular velocity |

## Roche lobe functions

| Function | Description |
|----------|-------------|
| `update_roche_radii(tessels, params, D)` | Solve Roche potential for radii at all vertices |
| `get_surface_potential(rpole, D, q, async_ratio)` | Potential at the pole |
| `update_roche_geom(star, params, D)` | Full Roche geometry update |
| `compute_potential_primary(r, D, theta, phi, q, async_ratio)` | Roche potential and derivatives for primary |
| `compute_potential_secondary(r, D, theta, phi, q, async_ratio)` | Roche potential and derivatives for secondary |
| `solve_radius(r0, pot, D, theta, phi, q, async_ratio, func)` | Find r such that potential = pot (Halley's method) |
| `solve_R_L1(r0, D, q, async_ratio, func)` | Find L1 Lagrange point |
| `compute_gravity_primary(r, theta, phi, D, q, async_ratio)` | Local gravity vector for primary |
| `compute_gravity_secondary(r, theta, phi, D, q, async_ratio)` | Local gravity vector for secondary |
| `radius_equivalent_eggleton(q)` | Eggleton (1983) Roche lobe radius approximation |
| `radius_leahy(q)` | Leahy & Leahy (2015) Roche lobe radius |
| `fillout_to_rpole(fillout, D, q, async_ratio)` | Convert fillout factor to polar radius |
| `rpole_to_fillout(rpole, D, q, async_ratio)` | Convert polar radius to fillout factor |
| `max_rpole(D, params)` | Maximum polar radius before L1 overflow |

## Rotation

| Function | Description |
|----------|-------------|
| `rotate_single_star(geom, params, orbit_incl, long_node, rot_angle)` | Rotate star for binary orbital motion |
| `compute_separation(binary_params, tepoch)` | Dimensionless binary separation at epoch |

## Visibility analysis

| Function | Description |
|----------|-------------|
| `sometimes_visible(stars)` | Pixels visible in at least one epoch |
| `never_visible(stars)` | Pixels never visible across all epochs |
| `invisible_neighbors(n, stars)` | Invisible pixels adjacent to visible ones |
| `with_invisible_neighbors(n, stars)` | Visible pixels that border invisible |
| `without_invisible_neighbors(n, stars)` | Interior visible pixels |
