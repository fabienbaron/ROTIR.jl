# Orbits

Orbital mechanics for binary star systems. Solves Kepler's equation and
computes positions, separations, and radial velocities.

## Kepler's equation

| Function | Description |
|----------|-------------|
| `compute_E_NR(M, e; T=Float64)` | Solve `E - e*sin(E) = M` via Newton-Raphson (adaptive for high eccentricity) |
| `compute_eccentric_anomaly(bparams, tepoch)` | Eccentric anomaly from orbital parameters, handles period derivative `dP` |
| `compute_true_anomaly(bparams, tepoch)` | True anomaly via Broucke-Cefola formula (numerically stable) |

## Orbital positions

| Function | Description |
|----------|-------------|
| `compute_coeff(Omega, i, omega)` | 6 Laplace coefficients for 3D position computation |
| `compute_xyz_rel(a, beta, e, L1, M1, N1, L2, M2, N2, cosE, sinE)` | Relative position (secondary w.r.t. primary) in observer frame |
| `binary_orbit_rel(bparams, tepoch)` | Relative orbit: returns `(0, 0, 0, x, y, z)` |
| `binary_orbit_abs(bparams, tepoch)` | Absolute orbit of both components w.r.t. center of mass: returns `(x1, y1, z1, x2, y2, z2)` |
| `binary_proj_plane(bparams, tepochs)` | Project orbit into observer plane: returns `(x, y, rho, theta)` |

## Coordinate conversion

| Function | Description |
|----------|-------------|
| `orbit_to_rotir_offset(bparams, tepoch)` | Convert orbital position (North, East) to ROTIR's (West, North) projected frame; returns `(offset_x, offset_y)` in mas |

## Separation

| Function | Description |
|----------|-------------|
| `compute_separation(bparams, tepoch)` | Dimensionless separation `D = 1 - e*cos(E)`. Multiply by semi-major axis `a` for physical distance |

## Radial velocities

| Function | Description |
|----------|-------------|
| `binary_RV(bparams, tepoch; K1, K2, gamma)` | Radial velocities of both stars given semi-amplitudes K1, K2 and systemic velocity gamma |

## Stellar parameters

`starparameters(...)` is a convenience constructor that returns a NamedTuple:

| Field | Description |
|-------|-------------|
| `rpole` | Polar radius (mas) |
| `tpole` | Polar temperature (K) |
| `frac_escapevel` | Fractional rotational velocity (0–1) |
| `ldtype` | Limb-darkening law (1=linear, 2=quadratic, 3=Hestroffer) |
| `ld1`, `ld2` | Limb-darkening coefficients |
| `beta_vZ` | von Zeipel gravity-darkening exponent |
| `B_rot` | Differential rotation coefficient |
| `inclination` | Spin-axis inclination (degrees) |
| `position_angle` | Spin-axis position angle (degrees) |
| `rotation_offset` | Fixed rotation offset (degrees) |
| `rotation_period` | Rotation period (days) |

## Binary parameters

`binaryparameters(...)` is a convenience constructor that returns a NamedTuple:

| Field | Description |
|-------|-------------|
| `star1`, `star2` | NamedTuple from `starparameters()` for each component |
| `d` | Distance (pc) |
| `i` | Orbital inclination (degrees) |
| `Omega` | Longitude of ascending node (degrees) |
| `omega` | Argument of periapsis (degrees) |
| `P` | Orbital period (days) |
| `a` | Semi-major axis (mas) |
| `e` | Eccentricity |
| `T0` | Time of periastron passage (JD) |
| `q` | Mass ratio M2/M1 |
| `fillout_factor` | Roche lobe fillout factors |
| `dP` | Period derivative (days/day) |
| `domega` | Periapsis precession (degrees/day) |
