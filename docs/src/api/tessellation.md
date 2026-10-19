# Tessellation

## Data structures

| Type | Description |
|------|-------------|
| `tessellation{T}` | unit-sphere tessellation: `tessellation_type` (0 HEALPix, 1 lon/lat), `npix`, `unit_xyz` and `unit_spherical` (npix × 5 × 3: four corners and the centre; `(x, y, z)` or `(r, colatitude, longitude)`) |
| `stellar_geometry{T}` | one epoch's surface: `surface_type`, `vertices_xyz`, `vertices_spherical`, `normals`, `index_quads_visible`, `nquads_visible`, `proj_west`, `proj_north`, `ldmap`, `vis_weights`, `polyflux`, `polyft`, epoch time `t` |

## HEALPix

Nested ordering: pixel `p` at level `n` has children `4(p−1)+1 … 4p` at level `n+1`.

| Function | Description |
|----------|-------------|
| `tessellation_healpix(n; T=Float32)` | HEALPix tessellation, nside = 2ⁿ, npix = 12·nside² |
| `nside2npix(nside)` | number of pixels |
| `npix2n(npix)` | level `n` from the number of pixels |
| `sobel_gradient_healpix(n; T=Float32, power=2)` | tangent-plane gradient `(Gx, Gy)` and squared-difference operator `D`, for `"sobel"`/`"sobel2"` |
| `laplace_beltrami_healpix(n; T=Float32)` | cotangent Laplace–Beltrami operator `D`, for `"laplace_beltrami"` |
| `tv_neighbors_healpix(n; T=Float32)` | graph-Laplacian difference matrix and its Hessian, for `"tv"`/`"tv2"` |
| `tv_neighbors_healpix_visible(n, stars; T=Float32)` | the same restricted to pixels visible at some epoch; use with pixel subset `sometimes_visible(stars)` |
| `upsample_map_stars(tmap, stars, star_params, tepochs)` | level `n+1`: each pixel copied to its four children; returns `(tmap, stars)` |
| `downsample_map_stars(tmap, stars, star_params, tepochs)` | level `n−1`: mean of each group of four; returns `(tmap, stars)` |
| `upsample_healpix_xmap(x, n_coarse, n_fine)` | any number of levels up |
| `downsample_healpix_xmap(x, n_fine, n_coarse)` | any number of levels down |

## Spots and bands

These take a `stellar_geometry` and work on either tessellation. Angles are in degrees:
latitude in [−90, 90], longitude in [−180, 180] or [0, 360).

| Function | Description |
|----------|-------------|
| `make_circ_spot(tmap, star, radius, lat, lon; bright_frac=0.8)` | sets the tessels within angular `radius` of (lat, lon) to `bright_frac` × median of the map |
| `ROTIR.make_circ_spot_spotfill(star, radius, lat, lon; profile="flat")` | spot filling factor per tessel: 1 inside (`"flat"`) or falling linearly to 0 at the edge (`"linear"`) |
| `ROTIR.make_meridional_band(tmap, star, width, lon; bright_frac=0.8)` | multiplies a band of longitude `width` centred on `lon` by `bright_frac` |
| `ROTIR.make_equatorial_band(tmap, star, width, lat; bright_frac=0.8)` | multiplies a band of latitude `width` centred on `lat` by `bright_frac` |

Spot membership uses the 3-D chord distance from the spot centre, so a spot stays circular on
a non-spherical surface.

## Longitude/latitude

Pixel `(i, j)`, ring `i` from the north pole and column `j` in longitude, has index
`(i−1)·nphi + j`.

| Function | Description |
|----------|-------------|
| `tessellation_latlong(ntheta, nphi; T=Float32)` | `ntheta` rings of `nphi` pixels |
| `tv_neighbors_longlat(ntheta, nphi)` | south and west neighbour lists |
| `make_circ_spot(tmap, ntheta, nphi, radius, lat, lon; bright_frac=0.8)` | spot on the pixel grid: `radius`, `lat`, `lon` are pixel indices; fills a 2·radius square with the map mean, and the disk inside it with `bright_frac` × mean |
| `make_spot_move(tmap, ntheta, nphi, period, B_rot, tepochs=[0.0, 1.0])` | shifts each ring in longitude by its rotation (differential rotation `B_rot`) between `tepochs[1]` and `tepochs[2]` |

## Internal HEALPix functions

| Function | Description |
|----------|-------------|
| `ROTIR.pix2ang_nest(nside, ipix)` | pixel index → (colatitude, longitude) |
| `ROTIR.pix2vec_nest(nside, ipix)` | pixel index → unit vectors of the centre and four corners |
| `ROTIR.ang2pix_nest(nside, theta, phi)` | (colatitude, longitude) → pixel index |
| `ROTIR.vec2ang(vector)` | unit vector → (colatitude, longitude) |
| `ROTIR.neighbors_nest(nside, ipix)` | the 7 or 8 neighbours of a pixel |
| `ROTIR.all_neighbors_nest(n)` | neighbours of every pixel |
