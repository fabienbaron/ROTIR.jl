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
