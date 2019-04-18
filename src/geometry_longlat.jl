
function latlong_ellipsoid_star(ntheta,nphi,a,b,c)
  npix = ntheta*nphi;
  vertices_spherical = zeros(Float64, npix, 5, 3); # r, theta, phi
  vertices_xyz = zeros(Float64, npix, 5, 3); # x, y, z
  dphi = 2*pi/nphi;
  dtheta = pi/ntheta;

  # get corners
  # calculates radius
  vertices_spherical[:,:, 1] .= 1.0;
  # calculates theta and phi
  theta = collect(range(0,stop=pi-dtheta,length=ntheta));
  phi = collect(range(0,stop=2*pi-dphi,length=nphi));

  # TODO: redo this more cleanly ?
  for i = 1:ntheta
    ilong_range = Int((i-1)*nphi+1):(i*nphi);
    vertices_spherical[ilong_range,1,2] .= theta[i]
    vertices_spherical[ilong_range,1,3] .= phi;
    vertices_spherical[ilong_range,2,2] .= theta[i] + dtheta
    vertices_spherical[ilong_range,2,3] .= phi;
    vertices_spherical[ilong_range,3,2] .= theta[i] + dtheta;
    vertices_spherical[ilong_range,3,3] .= phi .+ dphi;
    vertices_spherical[ilong_range,4,2] .= theta[i]
    vertices_spherical[ilong_range,4,3] .= phi .+ dphi;
    vertices_spherical[ilong_range,5,2] .= theta[i] + dtheta*0.5
    vertices_spherical[ilong_range,5,3] .= phi .+ dphi*0.5;
  end

  vertices_xyz[:,:,1] = a*sin.(vertices_spherical[:,:,2]).*cos.(vertices_spherical[:,:,3]); # X
  vertices_xyz[:,:,2] = b*sin.(vertices_spherical[:,:,2]).*sin.(vertices_spherical[:,:,3]); # Y
  vertices_xyz[:,:,3] = c*cos.(vertices_spherical[:,:,2]); # Z
  vertices_spherical[:,:,1] = sqrt.(vertices_xyz[:,:,1].^2 + vertices_xyz[:,:,2].^2 + vertices_xyz[:,:,3].^2); # radius, TBD replace by norm()

  star_base_geom = base_geometry(npix,vertices_xyz, vertices_spherical);
end

function latlong_rapidrot_star(ntheta,nphi,stellar_parameters)
  npix = ntheta*nphi;
  vertices_spherical = zeros(Float64, npix, 5, 3); # r, theta, phi
  vertices_xyz = zeros(Float64, npix, 5, 3); # x, y, z
  dphi = 2*pi/nphi;
  dtheta = pi/ntheta;

  # get corners
  # calculates theta and phi
  theta = collect(linspace(0,pi-dtheta,ntheta));
  phi = collect(linspace(0,2*pi-dphi,nphi));

  # make vertices -- go counterclockwise
  for i = 1:ntheta
    ilong_range = Int((i-1)*nphi+1):(i*nphi);

    vertices_spherical[ilong_range,1,2], vertices_spherical[ilong_range,1,3] = theta[i], phi;
    vertices_spherical[ilong_range,2,2], vertices_spherical[ilong_range,2,3] = theta[i] + dtheta, phi;
    vertices_spherical[ilong_range,3,2], vertices_spherical[ilong_range,3,3] = theta[i] + dtheta, phi + dphi;
    vertices_spherical[ilong_range,4,2], vertices_spherical[ilong_range,4,3] = theta[i], phi + dphi;
    theta_avg, phi_avg = theta[i] + dtheta*0.5, phi + dphi*0.5;
    vertices_spherical[ilong_range,5,2], vertices_spherical[ilong_range,5,3] = theta_avg, phi_avg;
  end

  R_pole = stellar_parameters.radius;
  ω = stellar_parameters.frac_escapevel;
  vertices_spherical[:,:,1] = 3.0*R_pole.*cos.((pi + acos.(ω*sin.(vertices_spherical[:,:,2])))/3.0)./(ω*sin.(vertices_spherical[:,:,2]));
  # Rewrite pole radius values
  # top of star
  vertices_spherical[1:nphi,1,1] = R_pole;
  vertices_spherical[1:nphi,4, 1] = R_pole;
  # bottom of star
  vertices_spherical[(end-nphi+1):end,2,1] = R_pole;
  vertices_spherical[(end-nphi+1):end,3,2] = R_pole;
  #vertices_spherical[find(vertices_spherical .== Inf)] = R_pole;

  vertices_xyz[:,:,1] = vertices_spherical[:,:,1].*sin.(vertices_spherical[:,:,2]).*cos.(vertices_spherical[:,:,3]); # X
  vertices_xyz[:,:,2] = vertices_spherical[:,:,1].*sin.(vertices_spherical[:,:,2]).*sin.(vertices_spherical[:,:,3]); # Y
  vertices_xyz[:,:,3] = vertices_spherical[:,:,1].*cos.(vertices_spherical[:,:,2]); # Z

  star_base_geom = base_geometry(npix,vertices_xyz, vertices_spherical);
end

function longlat_ang2pix(ntheta, nphi, THETA, PHI)
    vert_len = size(THETA,1);
    horiz_len = size(THETA,2);

    grid_pix = zeros(Int,vert_len,horiz_len);
    npix = 1;

    for i = 1:ntheta
        vert_range = ((i-1)*Int(vert_len/ntheta)+1):(i*Int(vert_len/ntheta));
        for j = 1:nphi
            horiz_range = ((j-1)*Int(horiz_len/nphi)+1):(j*Int(horiz_len/nphi));
            grid_pix[vert_range,horiz_range] .= npix;
            npix += 1;
        end
    end
    return grid_pix
end

function make_circ_spot(temperature_map,ntheta,nphi,spot_radius,lat,long;bright_frac=0.8)
    temperature_map = reshape(temperature_map,nphi,ntheta);
    ilat = collect((lat-spot_radius):(lat+spot_radius));
    ilong = collect((long-spot_radius):(long+spot_radius));
    indx_arr = Array{UnitRange{Int64},2}(length(ilat),2);
    for i = 1:length(ilat)
        indx_long_range = ilong[find(spot_radius^2 .>= ((ilat[i] - lat)^2 + (ilong - long).^2))];
        indx_arr[i,1] = indx_long_range[1]:indx_long_range[end];
        indx_arr[i,2] = ilat[i]:ilat[i];
    end
    for i = 1:length(ilat)
        temperature_map[indx_arr[i,1],indx_arr[i,2]] *= bright_frac;
    end
    temperature_map = vec(temperature_map);
    return temperature_map
end

function latlong_harmon(ntheta,nphi,a,b,c)
dtheta = pi/ntheta;
theta = collect(linspace(0,pi-dtheta,ntheta));
starting = dtheta/2
#ending = pi/2 -starting;
#Int.(ceil.(90*cos.(pi/2 - collect(linspace(starting,ending,ntheta/2)))))

# Close to pixelation scheme from LCI? (Used in Harmon & Crews 2000)
arr_long = Int64[];
arr_long = vcat(arr_long,starting)
for i = 2:ntheta
    arr_long = vcat(arr_long,arr_long[i-1]+dtheta)
end
indx_long = Int.(ceil.(nphi*cos.(pi/2 - arr_long)))
npix = sum(indx_long);
vertices_spherical = zeros(Float64, npix, 5, 3); # r, theta, phi
vertices_xyz = zeros(Float64, npix, 5, 3); # x, y, z
vertices_spherical[:,:,1] = 1.;

first=0
last=1
  for i = 1:ntheta
    nphi = indx_long[i];
    last = sum(indx_long[1:i])
    ilong_range=first+1:last
    dphi= 2*pi/nphi;
    phi = collect(linspace(0,2*pi-dphi,nphi));
    vertices_spherical[ilong_range,1,2], vertices_spherical[ilong_range,1,3] = theta[i], phi;
    vertices_spherical[ilong_range,2,2], vertices_spherical[ilong_range,2,3] = theta[i] + dtheta, phi;
    vertices_spherical[ilong_range,3,2], vertices_spherical[ilong_range,3,3] = theta[i] + dtheta, phi + dphi;
    vertices_spherical[ilong_range,4,2], vertices_spherical[ilong_range,4,3] = theta[i], phi + dphi;
    theta_avg, phi_avg = theta[i] + dtheta*0.5, phi + dphi*0.5;
    vertices_spherical[ilong_range,5,2], vertices_spherical[ilong_range,5,3] = theta_avg, phi_avg;
    first=deepcopy(last);
  end

  vertices_xyz[:,:,1] = a*sin.(vertices_spherical[:,:,2]).*cos.(vertices_spherical[:,:,3]); # X
  vertices_xyz[:,:,2] = b*sin.(vertices_spherical[:,:,2]).*sin.(vertices_spherical[:,:,3]); # Y
  vertices_xyz[:,:,3] = c*cos.(vertices_spherical[:,:,2]); # Z
  vertices_spherical[:,:,1] = sqrt.(vertices_xyz[:,:,1].^2 + vertices_xyz[:,:,2].^2 + vertices_xyz[:,:,3].^2); # radius, TBD replace by norm()

star_base_geom = base_geometry(npix,vertices_xyz, vertices_spherical);
end


function neighbors_longlat(ntheta,nphi)
    npix = ntheta*nphi;

    neighbors = Array{Array{Int64,1}}(undef,npix);

    for i=1:npix
        west_neighbor = i-1;
        west_neighbor_reverse = i+1;

        # correct if west pixel is on a different latitude
        if (i%nphi == 1)
            west_neighbor = i-1+nphi;
        elseif (i%nphi == 0)
            west_neighbor_reverse = i+1-nphi;
        end

        south_neighbor = i+nphi;
        south_neighbor_reverse = i-nphi

        # correct for poles
        if ((i <= nphi) & (i <= div(nphi,2)))
            south_neighbor_reverse = i+div(nphi,2);
        elseif (i <= nphi) & (i > div(nphi,2))
            south_neighbor_reverse = i-div(nphi,2);
        elseif ((i > (npix - nphi)) && (i <= (npix - div(nphi,2))))
            south_neighbor = i+div(nphi,2);
            #println("i = $i \t south = $south_neighbors")
        elseif (i > (npix - div(nphi,2)))
            south_neighbor = i-div(nphi,2);
            #println("i = $i \t south = $south_neighbors")
        end

        neighbors[i] = [south_neighbor,south_neighbor_reverse,west_neighbor,west_neighbor_reverse];
    end
    return neighbors
end

function tv_neighbors_longlat(ntheta,nphi)
    # Complete Neighbor setup (longlat)
    neighbors = neighbors_longlat(ntheta,nphi); #neighbors[ipix] will give the list of all neighbors of pixel [ipix]
    south_neighbors = Array{Int64}(undef,length(neighbors));
    west_neighbors = Array{Int64}(undef,length(neighbors));
    south_neighbors_reverse=Array{Array{Int64}}(undef,length(neighbors));
    west_neighbors_reverse=Array{Array{Int64}}(undef,length(neighbors));
    for i=1:length(neighbors) # TBD: find better expression with map or broadcast
        south_neighbors[i]=(neighbors[i])[1];
        west_neighbors[i]=(neighbors[i])[3];
    end

    for i=1:length(neighbors)
        #south_neighbors_reverse[i]=find(south_neighbors.==i);
        #west_neighbors_reverse[i]=find(west_neighbors.==i);
        south_neighbors_reverse[i]=[(neighbors[i])[2]];
        west_neighbors_reverse[i]=[(neighbors[i])[4]];
    end
    return neighbors,south_neighbors,west_neighbors,south_neighbors_reverse,west_neighbors_reverse
end
