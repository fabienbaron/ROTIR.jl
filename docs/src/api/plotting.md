# Plotting

All plotting functions use [PythonPlot](https://github.com/JuliaPy/PythonPlot.jl)
(Matplotlib). Import ROTIR to get access to these functions — the style is applied
automatically and matches OITOOLS, so figures from the two packages sit together in one
document without a font or tick mismatch.

## 2D projections

| Function | Description |
|----------|-------------|
| `plot2d(tmap, star; kwargs...)` | Plot temperature map in observer's sky plane for one epoch |
| `plot2d_wireframe(star; kwargs...)` | Tessellation only, no colour fill |
| `plot2d_allepochs(tmap, stars; kwargs...)` | Multi-epoch subplot grid with shared colour scale |

### `plot2d` keywords

| Keyword | Default | Description |
|---------|---------|-------------|
| `intensity` | `false` | Plot surface *brightness* (map × limb darkening) rather than temperature |
| `intensity_model` | `:linear` | `:linear` or `:planck`. Only applies when `intensity = true`, since a temperature map is band-independent |
| `band` | `nothing` | Wavelength in metres, required by `intensity_model = :planck` |
| `plotmesh` | `false` | Show tessel edges |
| `colormap` | `"gist_heat"` | Matplotlib colormap name |
| `figtitle` | `""` | Figure title |
| `flipx` | `false` | Flip East-West |
| `pad` | `0.5` | Padding around star (mas) |
| `background` | `"white"` | Background colour |
| `xlim`, `ylim` | `Float64[]` | Axis limits (auto if empty) |
| `vmin`, `vmax` | `nothing` | Pin the colour scale (auto if `nothing`) |
| `contours` | `Float64[]` | Contour levels to overlay |
| `contour_color` | `"gray"` | Contour line colour |
| `contour_labels` | `true` | Label the contours |
| `contour_fontsize` | `10` | Contour label size |

### `plot2d_wireframe` keywords

| Keyword | Default | Description |
|---------|---------|-------------|
| `hidden` | `true` | Also draw the far-side tessels, behind the near side |
| `front_color` | `"black"` | Edge colour for the near hemisphere |
| `hidden_color` | `"lightgrey"` | Edge colour for the far hemisphere |
| `linewidth` | `0.5` | Edge width |
| `compass` | `true` | Draw N/E compass arrows |
| `rotation_axis` | `false` | Draw the rotation axis |

`hidden = true` needs the near side unfilled, so the mesh reads as a transparent globe and
the two hemispheres are distinguished by edge weight. `hidden = false` restores the opaque
near-side-only view.

### `plot2d_allepochs` keywords

| Keyword | Default | Description |
|---------|---------|-------------|
| `plotmesh` | `false` | Show tessel edges |
| `tepochs` | `[]` | Epoch labels |
| `colormap` | `"gist_heat"` | Matplotlib colormap name |
| `ncols` | `nothing` | Pin the number of columns (auto: up to 4 across) |
| `arr_box` | `nothing` | Legacy two-digit `<rows><cols>` layout, e.g. `23` |
| `compass` | `true` | Draw the compass on the first panel |

The grid and figure size follow the number of epochs, so three epochs give a 1×3 row rather
than a half-empty 2×3 block.

### Decorations

| Keyword | Default | Description |
|---------|---------|-------------|
| `compass` | `true` | Draw N/E compass arrows |
| `limb` | `true` | Outline the projected stellar limb |
| `limb_color` | `"black"` | Limb colour |
| `limb_linewidth` | `0.8` | Limb width |
| `graticules` | `false` | Draw latitude/longitude grid lines on the surface |
| `rotation_axis` | `false` | Draw the rotation axis (dashed line through the poles) |
| `rotation_arrow` | `false` | Draw the spin-direction arrow at the north pole |
| `inclination`, `position_angle` | `NaN` | Override the orientation used for the decorations |
| `star_params` | `nothing` | Star parameters NamedTuple, for closed-form graticule *geometry* (sphere, triaxial ellipsoid, rapid rotator) |
| `graticule_kwargs` | `(;)` | Style overrides passed to `draw_graticules` (see below) |

!!! note "Where the decoration orientation comes from"
    With `inclination`/`position_angle` given, the decorations use those angles. Otherwise
    the orientation is recovered from the **mesh itself**, which is correct whatever built
    the star and also captures the rotation phase. `star_params` supplies the surface
    *shape*, not the orientation.

!!! note "Graticules on a Roche lobe"
    Sphere, triaxial ellipsoid and rapid rotator have closed forms, and `star_params`
    selects them. **Everything else — a Roche lobe (`surface_type = 3`), or any star drawn
    without `star_params` — reads the shape off the mesh's own `r(θ, φ)`**, interpolated
    between tessel centres.

    That is not a fallback so much as the general case. A Roche lobe is elongated toward
    its companion and has no axisymmetric description, so the curves have to come from the
    surface that was actually built. Visibility comes from the mesh normal rather than
    `z > 0`, which is what makes the curves stop exactly at the drawn limb instead of at
    the terminator of an imaginary sphere.

    Nothing extra is required to get this: pass `star_params` for a Roche star as you would
    for any other, or omit it entirely.

The limb outline matters more than it sounds: a surface that maps to the pale end of a
colormap is indistinguishable from a white background, so without it the disk can vanish
and take its decorations with it.

### Graticule style (`graticule_kwargs`)

| Key | Default | Description |
|-----|---------|-------------|
| `nlat` | `5` | Number of latitude circles |
| `nlon` | `8` | Number of longitude lines |
| `color` | `"black"` | Line colour |
| `linewidth` | `0.8` | Line width |
| `alpha` | `0.5` | Opacity |
| `npoints` | `200` | Points per curve (resolution) |
| `limb` | `true` | Include the projected limb with the graticule |

### Standalone decoration helpers

These take an existing axes, so they can be added to a figure you build yourself:

| Function | Description |
|----------|-------------|
| `draw_compass(ax, axis_max; ...)` | E/N compass arrows, sized to stay legible in small subplots |
| `draw_limb!(ax, star; ...)` | Projected limb (convex hull of the visible vertices) |
| `draw_graticules(ax, star; ...)` | Latitude/longitude grid on the surface |
| `draw_rotation_axis(ax, star; ...)` | Pole-to-pole line with an arrow at the north pole |
| `draw_rotation_arrow(ax, star; ...)` | Curved arrow showing the sense of rotation |
| `add_tessel_collection!(ax, star, colours; ...)` | The tessellated surface as one `PolyCollection` |

## Binary sky-plane

| Function | Description |
|----------|-------------|
| `plot2d_binary(tmap1, tmap2, star1, star2, bparams, tepoch; kwargs...)` | Plot a binary with depth ordering and the orbital offset |

`plot2d_binary` accepts the same decoration keywords as `plot2d` plus:

| Keyword | Default | Description |
|---------|---------|-------------|
| `star_params1`, `star_params2` | `nothing` | Star parameters for exact graticule geometry on each component |
| `inclination1`, `position_angle1` | `NaN` | Override star 1 orientation for decorations |
| `inclination2`, `position_angle2` | `NaN` | Override star 2 orientation for decorations |
| `vmin`, `vmax` | `nothing` | Pin the shared colour scale — essential across a frame sequence |
| `colorbar_on` | `true` | Draw the colour bar |
| `ax` | `nothing` | Draw into an existing axes instead of creating a figure |
| `axis_max` | `nothing` | Pin the plot half-width (mas) |
| `graticule_kwargs` | `(;)` | Graticule style overrides (shared by both components) |

!!! warning "`inclination1` / `position_angle1` are rarely what you want"
    `create_binary_geometry` orients **both** components by the shared
    `binary_frame`, built from the orbital elements `(i, Ω, ω)`. A component's own
    `inclination` / `position_angle` play no part. Passing the single-star angles here
    therefore decorates a differently-oriented star — for a Spica-like system the two
    answers differ by over 100°. Leave them at `NaN` and the mesh orientation is used.

## Radial velocity

| Function | Description |
|----------|-------------|
| `plot_rv(bparams; K1, K2, γ=0.0, rv_data1=nothing, rv_data2=nothing, figtitle="Radial Velocities")` | RV model curves vs orbital phase |

`rv_data1` / `rv_data2` should be N×3 matrices with columns `[JD, RV_km/s, error_km/s]`.
With only N×2, the points are drawn without error bars.

## 3D surface

| Function | Description |
|----------|-------------|
| `plot3d(tmap, star)` | 3D coloured surface (`Poly3DCollection`) |

## Mollweide projection

| Function | Description |
|----------|-------------|
| `plot_mollweide(tmap, star; kwargs...)` | Full-surface Mollweide projection (auto-selects HEALPix or lat/long) |

### `plot_mollweide` keywords

| Keyword | Default | Description |
|---------|---------|-------------|
| `visible_pixels` | `[]` | Pixel subset that was actually observed |
| `mask_unobserved` | `true` | Blank everything outside `visible_pixels` |
| `bad_color` | `"lightgray"` | Colour for masked pixels (needs `visible_pixels` to have any effect) |
| `vmin` | `-Inf` | Colour scale minimum (auto if `-Inf`) |
| `vmax` | `Inf` | Colour scale maximum (auto if `Inf`) |
| `colormap` | `"gist_heat"` | Matplotlib colormap name |
| `incl` | `90.0` | Draw the observability limit at latitude `-incl` |
| `lon_color` | `"white"` | Colour of the meridians and their tick labels |
| `lat_color` | `"black"` | Colour of the parallels and their tick labels |
| `figtitle` | `"Mollweide"` | Figure title |
| `ntheta`, `nphi` | `nothing` | Lat/long grid dimensions; recovered from the mesh when omitted |

`ntheta` / `nphi` apply only to a lat/long tessellation. `stellar_geometry` does not store
them, so they are recovered from the mesh unless you pin them.
