
"""
Healpix
"""
#=neighbors,south_neighbors,west_neighbors,south_neighbors_reverse,west_neighbors_reverse = tv_neighbours_healpix(n)

# spot 1
x = 5200*ones(star_epoch_geom[1].npix);

ipix = 230;
x[ipix] = 4200.;
x[neighbors[ipix]] = 4200;
x[vcat(neighbors[neighbors[ipix]]...)] = 4200.;
x[vcat(neighbors[vcat(neighbors[neighbors[230]]...)]...)]=4200.;
x[vcat(neighbors[vcat(neighbors[vcat(neighbors[neighbors[230]]...)]...)]...)]=4200.;

# spot 2
ipix = 367;
x[ipix] = 4500.;
x[neighbors[ipix]] = 4500;
x[vcat(neighbors[neighbors[ipix]]...)] = 4500.;
x[vcat(neighbors[vcat(neighbors[neighbors[230]]...)]...)]=4500.;=#

"""
Longitude/Latitude pixel scheme
"""
neighbors,south_neighbors,west_neighbors,south_neighbors_reverse,west_neighbors_reverse = tv_neighbors_longlat(ntheta,nphi)

# spot 1
x = 5200*ones(star_epoch_geom[1].npix);

ipix1 = 80*10+15;
x[vcat(neighbors[vcat(neighbors[vcat(neighbors[neighbors[ipix1]]...)]...)]...)]=4000.;
x[vcat(neighbors[vcat(neighbors[neighbors[ipix1]]...)]...)]=4000.;
x[vcat(neighbors[neighbors[ipix1]]...)] = 4200.;
x[neighbors[ipix1]] = 4200;
x[ipix] = 4200.;

# spot 2
ipix2 = 367;
x[vcat(neighbors[vcat(neighbors[neighbors[ipix2]]...)]...)]=4300.;
x[vcat(neighbors[neighbors[ipix2]]...)] = 4300.;
x[neighbors[ipix2]] = 4500;
x[ipix] = 4500.;

x_true = deepcopy(x);
