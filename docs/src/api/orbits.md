# Orbits

Kepler orbits for binary stars: positions, separations and radial velocities. Times are JD.
`bparams` is a NamedTuple with the fields of `binaryparameters` below; a Roche
`star_params` carries the same orbital fields and works too.

## Kepler's equation

| Function | Description |
|----------|-------------|
| `compute_E_NR(M, e)` | Solves `E − e·sin(E) = M` to machine precision, scalar or array `M` (`T` keyword accepted and ignored) |
| `compute_eccentric_anomaly(bparams, tepoch)` | Eccentric anomaly, with the period derivative `dP` |
| `compute_true_anomaly(bparams, tepoch)` | True anomaly, by the Broucke–Cefola formula |

## Orbital positions

| Function | Description |
|----------|-------------|
| `compute_coeff(Ω, i, ω)` | The six direction cosines `(L1, M1, N1, L2, M2, N2)`; angles in radians |
| `compute_xyz_rel(a, β, e, L1, M1, N1, L2, M2, N2, cosE, sinE)` | Secondary relative to primary, `β = √(1−e²)` |
| `binary_orbit_rel(bparams, tepoch)` | Relative orbit: `(0, 0, 0, x, y, z)` |
| `binary_orbit_abs(bparams, tepoch)` | Both components about the centre of mass: `(x1, y1, z1, x2, y2, z2)` |
| `binary_proj_plane(bparams, tepochs)` | Relative orbit on the sky: `(x, y, ρ, θ)`, `θ` in radians |
| `omega_at(bparams, tepoch)` | Argument of periapsis at `tepoch` with apsidal motion `dω` (radians) |

Positions are in mas, with x = North, y = East, z = receding. `ω` is the argument of
periapsis of the relative orbit; the primary's is `ω + π`.

## Coordinate conversion

| Function | Description |
|----------|-------------|
| `orbit_to_rotir_offset(bparams, tepoch)` | Secondary's offset from the primary in ROTIR's (West, North) frame, mas; scalar or vector of epochs |

## Separation

| Function | Description |
|----------|-------------|
| `compute_separation(bparams, tepoch)` | Separation in units of `a`, `D = 1 − e·cos(E)` |

## Radial velocities

| Function | Description |
|----------|-------------|
| `binary_RV(bparams, tepoch; K1, K2, γ)` | `(Vrad1, Vrad2)`; positive is receding; units of `K1`, `K2`, `γ` |

## Stellar parameters

`starparameters(rpole, tpole, frac_escapevel, ldtype, ld1, ld2, beta, B_rot, inclination,
position_angle, rotation_offset, rotation_period)` returns a NamedTuple with these fields:

| Field | Description |
|-------|-------------|
| `rpole` | Polar radius (mas) |
| `tpole` | Polar temperature (K) |
| `frac_escapevel` | Ω/Ω_crit |
| `ldtype` | Limb-darkening law (0 = none, 1 = linear, 2 = quadratic, 3 = power law, 4 = Claret) |
| `ld1`, `ld2` | Limb-darkening coefficients |
| `beta` | Gravity-darkening exponent |
| `B_rot` | Differential rotation coefficient |
| `inclination` | Spin-axis inclination (degrees) |
| `position_angle` | Spin-axis position angle (degrees) |
| `rotation_offset` | Fixed rotation offset (degrees) |
| `rotation_period` | Rotation period (days) |

## Binary parameters

`binaryparameters(star1, star2, d, i, Ω, ω, P, a, e, T0, q, fillout_factor, dP, dω)`
returns a NamedTuple with these fields:

| Field | Description |
|-------|-------------|
| `star1`, `star2` | NamedTuple from `starparameters()` for each component |
| `d` | Distance (pc) |
| `i` | Orbital inclination (degrees) |
| `Ω` | Longitude of the ascending node (degrees) |
| `ω` | Argument of periapsis of the relative orbit (degrees) |
| `P` | Orbital period (days) |
| `a` | Semi-major axis of the relative orbit (mas) |
| `e` | Eccentricity |
| `T0` | Time of periastron (JD) |
| `q` | Mass ratio M₂/M₁ |
| `fillout_factor` | Roche-lobe fill-out factors `[primary, secondary]` |
| `dP` | Period derivative (days/day) |
| `dω` | Apsidal motion (degrees/day) |

`Ω`, `ω` and `dω` are Unicode field names; ASCII `Omega`, `omega` or `domega` are not read.
