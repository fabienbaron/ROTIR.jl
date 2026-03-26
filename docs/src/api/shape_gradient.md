# Shape gradients

Analytical gradients of the chi-squared with respect to shape parameters
(radii, inclination, position angle). Enables joint optimization of the surface
map and stellar geometry.

## Rotation matrix

| Function | Description |
|----------|-------------|
| `rotation_matrix(psi, inc, PA)` | 3x3 rotation matrix R(psi, inc, PA), all angles in radians. Same as `rot_vertex` in geometry.jl |
| `dR_dinc(psi, inc, PA)` | Derivative of R w.r.t. inclination |
| `dR_dPA(psi, inc, PA)` | Derivative of R w.r.t. position angle |

Convention: `xyz = su * R` (right-multiply), where `su` are the scaled
unit-sphere vertices.

## Rapid rotator derivative

| Function | Description |
|----------|-------------|
| `f_rapid_rot_and_deriv(x)` | Returns `(f(x), f'(x))` where `f(x) = 3*cos((pi+acos(x))/3)/x`. Near x=0, returns `(1, 0)` |

## Projected vertices and derivatives

| Function | Description |
|----------|-------------|
| `projected_vertices_and_derivs(tessels, star_params, t; nparams)` | Compute projected vertices and their analytical derivatives w.r.t. shape parameters |

Returns `(projx, projy, dprojx_dtheta, dprojy_dtheta, normals_z, dnz_dtheta)`:

| Output | Shape | Description |
|--------|-------|-------------|
| `projx`, `projy` | (npix, 4) | Projected quad vertex coordinates |
| `dprojx_dtheta`, `dprojy_dtheta` | (npix, 4, nparams) | Vertex derivatives w.r.t. shape parameters |
| `normals_z` | (npix,) | Normalized z-component of face normals |
| `dnz_dtheta` | (npix, nparams) | Normal z-component derivatives |

The derivative computation is surface-type-specific:

| Surface | Parameter | Derivative |
|---------|-----------|------------|
| Sphere | radius | `d(su)/d(r) = unit_xyz` |
| Ellipsoid | rx, ry, rz | `d(su_j)/d(r_j) = unit_xyz_j` (axis-aligned) |
| Rapid rotator | rpole | `d(su)/d(rpole) = su / rpole` |
| Rapid rotator | omega | `d(r)/d(omega) = rpole * f'(omega*sin(theta)) * sin(theta)` |
| All | inc | `d(xyz)/d(inc) = su * dR_dinc` |
| All | PA | `d(xyz)/d(PA) = su * dR_dPA` |

## Shape chi2 + gradient

| Function | Description |
|----------|-------------|
| `shape_chi2_fg!(grad_theta, grad_xmap, xmap, theta, data, tessels, star_params_base, tepochs; kappa, verbose)` | Chi2 and gradients w.r.t. both shape and map parameters |

The gradient chain has three components:

1. **Vertex positions**: `d(chi2)/d(theta) = sum_p sum_v (d(chi2)/d(projx) * d(projx)/d(theta) + ...)`
2. **Flux normalization**: shoelace area derivatives for the flux correction term
3. **Soft visibility**: `d(chi2)/d(theta) += d(chi2)/d(w) * d(sigmoid)/d(kappa*nz) * kappa * d(nz)/d(theta)`

### Shape parameter vector layout

| Surface type | theta vector | nparams |
|-------------|-------------|---------|
| 0 (Sphere) | `[radius, inclination, PA]` | 3 |
| 1 (Ellipsoid) | `[rx, ry, rz, inclination, PA]` | 5 |
| 2 (Rapid Rotator) | `[rpole, omega, inclination, PA]` | 4 |

Inclination and PA are in **degrees** (converted internally).

## Joint reconstruction

| Function | Description |
|----------|-------------|
| `joint_reconstruct_oi(xmap, theta, data, tessels, star_params, tepochs; kwargs...)` | Alternating optimization of map and shape parameters |

### Keywords

| Keyword | Default | Description |
|---------|---------|-------------|
| `maxiter_xmap` | `200` | Max iterations for map optimization |
| `maxiter_theta` | `50` | Max iterations for shape optimization |
| `nouter` | `5` | Number of alternating cycles |
| `reg_weight` | `1e-5` | TV2 regularization weight for map step |
| `kappa` | `50.0` | Sigmoid sharpness for soft visibility |
| `theta_lower` | `nothing` | Lower bounds on shape parameters |
| `theta_upper` | `nothing` | Upper bounds on shape parameters |
| `verbose` | `true` | Print diagnostics |
