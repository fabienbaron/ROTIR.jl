# Plotting

Two plotting layers, each a package extension:

- **PythonPlot** (matplotlib): `plot2d`, `plot2d_wireframe`, `plot2d_allepochs`,
  `plot2d_binary`, `plot3d`, `plot_rv`, `plot_mollweide` and the `draw_*` helpers. They use
  OITOOLS' matplotlib style (`set_oiplot_defaults`).
- **Makie** (any backend): the same plots with a `_makie` suffix, each returning
  `(fig, ax)` (see [Makie equivalents](@ref)).

Plot coordinates: x = East = `-proj_west`, y = North, in mas; East is to the left.

## 2D projections

| Function | Description |
|----------|-------------|
| `plot2d(tmap, star; kwargs...)` | Map on the sky plane at one epoch; returns `(fig, ax)` |
| `plot2d_wireframe(star; kwargs...)` | Tessellation edges only |
| `plot2d_allepochs(tmap, stars; kwargs...)` | One panel per epoch, shared colour scale |

### `plot2d` keywords

| Keyword | Default | Description |
|---------|---------|-------------|
| `intensity` | `false` | Colour by surface brightness `I(T) × ldmap` instead of temperature |
| `intensity_model` | `:linear` | `I(T)`: `:linear` (`I = T`) or `:planck`; used only with `intensity = true` |
| `band` | `nothing` | Wavelength in metres, required by `:planck` |
| `plotmesh` | `false` | Draw tessel edges |
| `colormap` | `"gist_heat"` | Matplotlib colormap name |
| `figtitle` | `""` | Axes title |
| `flipx` | `false` | East to the right |
| `pad` | `0.5` | Margin around the star (mas) |
| `background` | `"white"` | `"white"` or `"black"` |
| `vmin`, `vmax` | `nothing` | Fix the colour scale (automatic if `nothing`) |
| `contours` | `Float64[]` | Contour levels |
| `contour_color` | `"gray"` | Contour colour |
| `contour_labels` | `true` | Label the contours (in K) |
| `contour_fontsize` | `10` | Contour label size |

Fix `vmin`/`vmax` when rendering a sequence of frames, so the colour scale is the same in all.

### `plot2d_wireframe` keywords

| Keyword | Default | Description |
|---------|---------|-------------|
| `hidden` | `true` | Also draw the far-side edges, behind an unfilled near side |
| `front_color` | `"black"` | Near-side edge colour |
| `hidden_color` | `"lightgrey"` | Far-side edge colour |
| `linewidth` | `0.5` | Edge width |
| `compass` | `true` | N/E compass |
| `rotation_axis` | `false` | Rotation axis |

`hidden = false` draws the near side only, filled white.

### `plot2d_allepochs` keywords

| Keyword | Default | Description |
|---------|---------|-------------|
| `plotmesh` | `false` | Draw tessel edges |
| `tepochs` | `[]` | Epoch labels for the panel titles |
| `colormap` | `"gist_heat"` | Matplotlib colormap name |
| `ncols` | `nothing` | Number of columns (default: up to 4) |
| `arr_box` | `nothing` | Layout as a two-digit `<rows><cols>` number, e.g. `23` |
| `compass` | `true` | Compass on the first panel |

The figure size follows the number of rows and columns.

### Decorations

Keywords of `plot2d` (and `plot2d_binary`):

| Keyword | Default | Description |
|---------|---------|-------------|
| `compass` | `true` | N/E compass arrows |
| `limb` | `true` | Outline of the projected limb |
| `limb_color` | `"black"` | Limb colour (also the graticule colour) |
| `limb_linewidth` | `0.8` | Limb width |
| `graticules` | `false` | Latitude/longitude lines on the surface |
| `rotation_axis` | `false` | Rotation axis, dashed, with an arrow at the north pole |
| `rotation_arrow` | `false` | Arrow around the north pole showing the sense of rotation |
| `inclination`, `position_angle` | `NaN` | Orientation for the decorations (degrees); from the mesh if `NaN` |
| `star_params` | `nothing` | Star parameters, for the closed-form surface shape of the graticules |
| `graticule_kwargs` | `(;)` | Passed to `draw_graticules` (see below) |

The orientation comes from `inclination`/`position_angle` when given, otherwise from the mesh,
which includes the rotation phase. `star_params` gives the surface shape only.

Graticules use a closed-form shape for a sphere, triaxial ellipsoid or rapid rotator when
`star_params` is given. Otherwise, including for a Roche lobe (`surface_type = 3`), the shape
is interpolated from the mesh's `r(θ, φ)`. Curves end at the limb, where the mesh normal turns
away from the observer.

The limb outline keeps the disk visible when the surface maps to the pale end of the colormap
on a white background.

### Graticule style (`graticule_kwargs`)

| Key | Default | Description |
|-----|---------|-------------|
| `nlat` | `5` | Number of latitude circles |
| `nlon` | `8` | Number of meridians |
| `color` | `"black"` | Line colour |
| `linewidth` | `0.8` | Line width |
| `alpha` | `0.5` | Opacity |
| `npoints` | `200` | Points per curve |
| `limb` | `true` | Also draw the limb |

### Decoration helpers

Each draws into an existing matplotlib axes:

| Function | Description |
|----------|-------------|
| `draw_compass(ax, axis_max; size_frac=0.12, fontsize=nothing, color="black")` | N/E arrows in the lower-right corner; `fontsize = nothing` scales with the axes size |
| `draw_limb!(ax, star; ...)` | Projected limb (convex hull of the visible vertices) |
| `draw_graticules(ax, star; ...)` | Latitude/longitude lines |
| `draw_rotation_axis(ax, star; ...)` | Pole-to-pole line with an arrow at the north pole |
| `draw_rotation_arrow(ax, star; ...)` | Arrow around the pole showing the sense of rotation |
| `add_tessel_collection!(ax, star, colours; ...)` | The visible tessels as one `PolyCollection` |

## Binary sky-plane

| Function | Description |
|----------|-------------|
| `plot2d_binary(tmap1, tmap2, star1, star2, bparams, tepoch; kwargs...)` | Both components at their orbital positions, nearer one in front; returns `(fig, ax)` |

`tepoch` is in the time system of `bparams.T0`. `plot2d_binary` takes the decoration keywords of
`plot2d`, `intensity`, `intensity_model`, `band`, `plotmesh`, `colormap`, `background`,
`figtitle`, `pad` (default `1.0`), and:

| Keyword | Default | Description |
|---------|---------|-------------|
| `star_params1`, `star_params2` | `nothing` | Star parameters for each component's graticules |
| `inclination1`, `position_angle1` | `NaN` | Orientation of star 1 for the decorations |
| `inclination2`, `position_angle2` | `NaN` | Orientation of star 2 for the decorations |
| `vmin`, `vmax` | `nothing` | Fix the shared colour scale |
| `colorbar_on` | `true` | Draw the colour bar |
| `ax` | `nothing` | Draw into an existing axes |
| `axis_max` | `nothing` | Plot half-width (mas) |
| `graticule_kwargs` | `(;)` | Graticule style, both components |

!!! warning "`inclination1` / `position_angle1`"
    `create_binary_geometry` orients both components from the orbital elements `(i, Ω, ω)`
    through `binary_frame`; a component's own `inclination` and `position_angle` are not used.
    Passing them here draws the decorations for a different orientation (over 100° apart for
    a Spica-like system). Leave them at `NaN` to use the mesh orientation.

## Radial velocity

| Function | Description |
|----------|-------------|
| `plot_rv(bparams; K1, K2, γ=0.0, rv_data1=nothing, rv_data2=nothing, figtitle="Radial Velocities")` | Model RV curves against orbital phase; returns `(fig, ax)` |

`rv_data1` / `rv_data2` are N×3 matrices `[time, RV, error]` (km/s), with time in the system of
`bparams.T0`. N×2 matrices are drawn without error bars.

## 3D surface

| Function | Description |
|----------|-------------|
| `plot3d(tmap, star)` | Coloured 3D surface (`Poly3DCollection`) |

## Mollweide projection

| Function | Description |
|----------|-------------|
| `plot_mollweide(tmap, star; kwargs...)` | Whole surface in Mollweide projection, HEALPix or lat/long mesh |

### `plot_mollweide` keywords

| Keyword | Default | Description |
|---------|---------|-------------|
| `visible_pixels` | `[]` | Observed pixels; sets the colour range |
| `mask_unobserved` | `true` | Paint pixels outside `visible_pixels` in `bad_color` |
| `bad_color` | `"lightgray"` | Colour of masked pixels |
| `vmin`, `vmax` | `-Inf`, `Inf` | Colour scale (automatic if infinite) |
| `colormap` | `"gist_heat"` | Matplotlib colormap name |
| `incl` | `90.0` | Draw a dash-dot line at latitude `-incl` (none at 90) |
| `lon_color` | `"white"` | Colour of the meridians and their labels |
| `lat_color` | `"black"` | Colour of the parallels and their labels |
| `figtitle` | `"Mollweide"` | Figure window name; also the axes title on a lat/long mesh |
| `title` | `nothing` | Axes title |
| `ntheta`, `nphi` | `nothing` | Lat/long grid size; read from the mesh if omitted |

### Controlling the furniture

| Keyword | Default | Description |
|---------|---------|-------------|
| `graticule` | `true` | Draw meridians and parallels |
| `lon_spacing` | `20` (HEALPix), `30` (lat/long) | Meridian spacing (degrees) |
| `lat_spacing` | `15` (HEALPix), `30` (lat/long) | Parallel spacing (degrees); must divide 90 |
| `labels` | `true` | Degree labels |
| `grid_style` | `"-"` | Graticule linestyle, e.g. `":"` |
| `grid_alpha` | `0.6` | Graticule opacity |
| `axis_lines` | `false` | Draw the equator and central meridian as solid lines |
| `axis_color`, `axis_style`, `axis_width` | `"black"`, `"-"`, `1.0` | Style of those two lines |
| `colorbar` | `true` | Draw the colour bar |
| `colorbar_label` | `"Temperature (K)"` | Colour bar label; `""` for none |
| `colorbar_ticks` | `7` | Number of colour bar ticks |

A `lat_spacing` that does not divide 90 raises an error: matplotlib labels the parallels in
multiples of the spacing, so other values put the labels on the wrong lines.

The settings below match healpy's `mollview` (dotted graticule, solid equator and central
meridian, no degree labels, colour bar with its two end values). `demos/zetand_render_map.jl`
uses them.

```julia
plot_mollweide(map, star; vmin = 0.8025, vmax = 0.9977, colormap = "gist_heat",
               title = "ζ And 2019A",
               labels = false, lon_color = "black", lat_color = "black",
               grid_style = ":", grid_alpha = 0.9,
               lon_spacing = 30, lat_spacing = 30, axis_lines = true,
               colorbar_ticks = 2, colorbar_label = "")
```

## Makie equivalents

Available once a Makie backend is loaded (`CairoMakie` for files, `GLMakie` for windows). Each
returns `(fig, ax)`; `plot2d_allepochs_makie` returns `fig`.

| Function | Matplotlib counterpart |
|----------|------------------------|
| `plot2d_makie(tmap, star; ...)` | `plot2d`; adds `colorbar` and `figsize` |
| `plot2d_wireframe_makie(star; ...)` | `plot2d_wireframe` (see note below) |
| `plot2d_allepochs_makie(tmap, stars; ...)` | `plot2d_allepochs`; adds `pad`, `panelsize` |
| `plot2d_binary_makie(tmap1, tmap2, star1, star2, bparams, tepoch; ...)` | `plot2d_binary` |
| `plot3d_makie(tmap, star; ...)` | `plot3d`; rotatable, orthographic, `axistype = :lscene` or `:axis3` |
| `plot3d_binary_makie(tmap1, tmap2, star1, star2, bparams, tepoch; ...)` | both components in 3D, with the relative orbit (`orbit = true`) |
| `plot_rv_makie(bparams; K1, K2, ...)` | `plot_rv` |
| `plot_mollweide_makie(tmap, star; ...)` | `plot_mollweide` |
| `draw_compass_makie!`, `draw_limb_makie!`, `draw_graticules_makie!`, `draw_rotation_axis_makie!`, `draw_rotation_arrow_makie!`, `add_tessel_collection_makie!` | the `draw_*` helpers |

`plot_mollweide_makie` takes `visible_pixels`, `mask_unobserved`, `bad_color`, `vmin`, `vmax`,
`colormap`, `figtitle`, `incl`, `lon_color`, `lat_color`, `labels`, `colorbar`, `nlon`/`nlat`
(resampling grid, default 720×360) and `figsize`; not the styling keywords of
[Controlling the furniture](@ref).

!!! note "`hidden` in the wireframe plots"
    `plot2d_wireframe(star; hidden = true)` draws the far side; `plot2d_wireframe_makie(star;
    hidden = true)` does not, and draws it with `hidden = false`.
