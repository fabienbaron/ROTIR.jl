# Plotting

All plotting functions use PyPlot (Matplotlib). Import ROTIR to get access to
these functions.

## 2D projections

| Function | Description |
|----------|-------------|
| `plot2d(tmap, star; kwargs...)` | Plot temperature map in observer's sky plane for one epoch |
| `plot2d_wireframe(star)` | Overlay projected pixel edges (no color fill) |
| `plot2d_allepochs(tmap, stars; kwargs...)` | Multi-epoch subplot grid with shared color scale |

### `plot2d` keywords

| Keyword | Default | Description |
|---------|---------|-------------|
| `intensity` | `false` | Multiply map by limb-darkening weights |
| `plotmesh` | `false` | Show pixel edges |
| `colormap` | `"gist_heat"` | Matplotlib colormap name |
| `figtitle` | `""` | Figure title |
| `flipx` | `false` | Flip East-West |
| `pad` | `0.5` | Padding around star (mas) |
| `background` | `"black"` | Background color |
| `xlim`, `ylim` | `[]` | Axis limits (auto if empty) |

### `plot2d_allepochs` keywords

| Keyword | Default | Description |
|---------|---------|-------------|
| `plotmesh` | `false` | Show pixel edges |
| `tepochs` | `[]` | Epoch labels |
| `colormap` | `"gist_heat"` | Matplotlib colormap name |
| `arr_box` | auto | Subplot layout (e.g., `23` = 2 rows, 3 columns) |

### Decorations

| Keyword | Default | Description |
|---------|---------|-------------|
| `compass` | `false` | Draw N/E compass arrows |
| `graticules` | `false` | Draw latitude/longitude grid lines on the surface |
| `rotation_axis` | `false` | Draw the rotation axis (dashed line through poles) |
| `rotation_arrow` | `false` | Draw spin direction arrow at the north pole |

## Binary sky-plane

| Function | Description |
|----------|-------------|
| `plot2d_binary(tmap1, tmap2, star1, star2, bparams, tepoch; kwargs...)` | Plot binary system with correct occlusion and orbital offset |

`plot2d_binary` accepts the same decoration keywords as `plot2d` plus:

| Keyword | Default | Description |
|---------|---------|-------------|
| `inclination1`, `position_angle1` | `NaN` | Override star 1 orientation for decorations |
| `inclination2`, `position_angle2` | `NaN` | Override star 2 orientation for decorations |

## Radial velocity

| Function | Description |
|----------|-------------|
| `plot_rv(bparams; K1, K2, γ=0.0, rv_data1=nothing, rv_data2=nothing, figtitle="Radial Velocities")` | Plot RV model curves vs orbital phase |

`rv_data1` / `rv_data2` should be Nx3 matrices with columns `[JD, RV_km/s, error_km/s]`.
If only Nx2, data points are plotted without error bars.

## 3D surface

| Function | Description |
|----------|-------------|
| `plot3d(tmap, star)` | 3D colored surface (Poly3DCollection) |
| `plot3d_vertices(star)` | Debug: show quad vertices (blue) and centers (red) |

## Mollweide projection

| Function | Description |
|----------|-------------|
| `plot_mollweide(tmap, star; kwargs...)` | Full-surface Mollweide projection (auto-selects HEALPix or lon/lat) |

### `plot_mollweide` keywords

| Keyword | Default | Description |
|---------|---------|-------------|
| `visible_pixels` | `[]` | Highlight visible pixel subset |
| `vmin` | `-Inf` | Color scale minimum (auto if -Inf) |
| `vmax` | `Inf` | Color scale maximum (auto if Inf) |
| `colormap` | `"gist_heat"` | Matplotlib colormap name |
| `incl` | `90.0` | Draw inclination line at this angle |
| `figtitle` | `"Mollweide"` | Figure title |
