# Shape gradients

Analytic gradients of χ² with respect to the shape parameters (radii, inclination, position
angle), for joint fits of the surface map and the geometry.

## Rotation matrix

| Function | Description |
|----------|-------------|
| `rotation_matrix(psi, inc, PA)` | 3×3 rotation matrix R(psi, inc, PA), angles in radians; same as `rot_vertex` |
| `dR_dinc(psi, inc, PA)` | ∂R/∂inc |
| `dR_dPA(psi, inc, PA)` | ∂R/∂PA |

`xyz = su * R` (right multiplication), where `su` are the scaled unit-sphere vertices.

## Rapid rotator derivative

| Function | Description |
|----------|-------------|
| `f_rapid_rot_and_deriv(x)` | `(f(x), f'(x))` with `f(x) = 3cos((π+acos x)/3)/x = 3sin(asin(x)/3)/x`; below `|x| = eps^(1/4)` the series `(1 + 4x²/27, 8x/27)` |

## Projected vertices and derivatives

| Function | Description |
|----------|-------------|
| `projected_vertices_and_derivs(tessels, star_params, t; nparams, derivs=true)` | projected vertices and their derivatives with respect to the shape parameters |

Returns `(proj_west, proj_north, dproj_west_dtheta, dproj_north_dtheta, normals_z, dnz_dtheta)`:

| Output | Shape | Description |
|--------|-------|-------------|
| `proj_west`, `proj_north` | (npix, 4) | projected quad vertices |
| `dproj_west_dtheta`, `dproj_north_dtheta` | (npix, 4, nparams) | vertex derivatives |
| `normals_z` | (npix,) | z-component of the unit face normals |
| `dnz_dtheta` | (npix, nparams) | its derivatives |

Per surface type:

| Surface | Parameter | Derivative |
|---------|-----------|------------|
| Sphere | radius | `d(su)/d(r) = unit_xyz` |
| Ellipsoid | rx, ry, rz | `d(su_j)/d(r_j) = unit_xyz_j` (axis-aligned) |
| Rapid rotator | rpole | `d(su)/d(rpole) = su / rpole` |
| Rapid rotator | omega | `d(r)/d(omega) = rpole * f'(omega*sin(theta)) * sin(theta)` |
| All | inc | `d(xyz)/d(inc) = su * dR_dinc` |
| All | PA | `d(xyz)/d(PA) = su * dR_dPA` |

## Shape χ² and gradient

`shape_chi2_fg!(grad_θ, grad_xmap, xmap, θ, data, tessels, star_params_base, tepochs; κ=50,
verbose=false, parametric_map=false, fit_v0=false, weights=(1,1,1))` returns χ² and fills
the gradients with respect to the shape parameters and the map.

- `parametric_map = true` recomputes the map from θ (overwriting `xmap`) and adds its
  derivative to `grad_θ`. Surface types 0 and 1 only.
- `fit_v0 = true` reads a visibility scale `V0` from the last entry of θ (see
  `interferometric_chi2` in [Parametric model and AD](parametric_gradient.md)).
- `weights = (w_V², w_T3amp, w_T3φ)` scale each χ² term and its gradients.

The gradient has three parts:

1. **Vertex positions**: `d(chi2)/d(theta) = sum_p sum_v (d(chi2)/d(proj_west) * d(proj_west)/d(theta) + ...)`
2. **Flux normalization**: shoelace area derivatives for the flux correction term
3. **Per-tessel weight** `w = sigmoid(kappa*nz) * ld(nz)`:
   `d(chi2)/d(theta) += d(chi2)/d(w) * d(w)/d(nz) * d(nz)/d(theta)`, with
   `d(w)/d(nz) = d(sigmoid)/d(kappa*nz) * kappa * ld + sigmoid * d(ld)/d(nz)`

Limb darkening uses `ldtype`, `ld1`, `ld2` from `star_params_base` and `mu = max(nz, 0)`, as
`create_star`/`compute_ldmap` do. The LD coefficients are held fixed: with a free map they
are largely degenerate with it. The [parametric model](parametric_gradient.md) fits them.

### Shape parameter vector layout

| Surface type | theta vector | nparams |
|-------------|-------------|---------|
| 0 (Sphere) | `[radius, inclination, PA]` | 3 |
| 1 (Ellipsoid) | `[rx, ry, rz, inclination, PA]` | 5 |
| 2 (Rapid Rotator) | `[rpole, frac_escapevel, inclination, PA]` | 4 |

Inclination and PA are in degrees. With `fit_v0 = true`, `V0` is appended. Surface type 3
(Roche) has no shape gradient.

## Joint reconstruction

| Function | Description |
|----------|-------------|
| `joint_reconstruct_oi(xmap, θ, data, tessels, star_params, tepochs; kwargs...)` | returns `(xmap, θ)`; alternates VMLMB on the map (shape fixed) and on θ (map fixed, bounded) |

### Keywords

| Keyword | Default | Description |
|---------|---------|-------------|
| `maxiter_xmap` | `200` | iterations per map step |
| `maxiter_θ` | `50` | iterations per shape step |
| `nouter` | `5` | alternating cycles |
| `regularizer` | `"sobel"` | regulariser of the map step |
| `reg_weight` | `10` | its weight |
| `κ` | `50` | soft-visibility sharpness |
| `θ_lower`, `θ_upper` | `nothing` | bounds on the shape parameters |
| `fit_v0` | `false` | also fit the visibility scale `V0` (last entry of θ) |
| `verbose` | `true` | print progress |

Surface type 3 (Roche) is refused; reconstruct at a fixed Roche geometry with
`image_reconstruct_oi`.
