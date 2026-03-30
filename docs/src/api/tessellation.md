# Tessellation

## Data structures

| Type | Description |
|------|-------------|
| `tessellation{T}` | Base tessellation on unit sphere: `npix`, `unit_xyz` (npix,5,3), `unit_spherical` (npix,5,3), `tessellation_type` (0=HEALPix, 1=lon/lat) |
| `stellar_geometry{T}` | Rotated/projected geometry for one epoch: vertices, normals, visibility mask, proj_west/proj_north, polyflux, polyft, epoch time `t` |

## HEALPix

| Function | Description |
|----------|-------------|
| `tessellation_healpix(n; T=Float32)` | Create HEALPix tessellation with nside=2^n, npix=12*nside^2 |
| `nside2npix(nside)` | Convert nside to number of pixels |
| `npix2n(npix)` | Convert npix to resolution parameter n |
| `tv_neighbors_healpix(n; T=Float32)` | Sparse difference matrix and Hessian for TV regularization |
| `tv_neighbors_healpix_visible(n, stars; T=Float32)` | TV neighbors excluding always-invisible pixels |
| `upsample_map_stars(tmap, stars, star_params, tepochs)` | Double resolution (each pixel splits into 4) |
| `downsample_map_stars(tmap, stars, star_params, tepochs)` | Halve resolution (average groups of 4) |

## Longitude/latitude

| Function | Description |
|----------|-------------|
| `tessellation_latlong(ntheta, nphi; T=Float32)` | Create regular lat/lon grid with ntheta rings of nphi pixels |
| `tv_neighbors_longlat(ntheta, nphi)` | TV neighbor structure for lat/lon grid |
| `make_circ_spot(tmap, ntheta, nphi, radius, lat, lon; bright_frac)` | Create circular temperature spot |
| `make_spot_move(tmap, ntheta, nphi, period, B_rot, tepochs)` | Move spot in longitude with differential rotation |

## Internal HEALPix functions

| Function | Description |
|----------|-------------|
| `pix2ang_nest(nside, ipix)` | Pixel index to (colatitude, longitude) |
| `pix2vec_nest(nside, ipix)` | Pixel index to 3D unit vectors (center + 4 vertices) |
| `ang2pix_nest(nside, theta, phi)` | (colatitude, longitude) to pixel index |
| `vec2ang(vector)` | 3D vector to (colatitude, longitude) |
| `neighbors_nest(nside, ipix)` | Find 8 neighbors of a pixel in nested scheme |
| `all_neighbors_nest(n)` | Neighbors for all pixels |
