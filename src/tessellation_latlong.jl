@views function tessellation_latlong(ntheta,nphi; T=Float32)
    npix = ntheta*nphi;
    unit_spherical = zeros(T, npix, 5, 3); # r, theta, phi
    unit_xyz = zeros(T, npix, 5, 3); # x, y, z
    dphi   = 2pi/nphi;
    dtheta = pi/ntheta;
    unit_spherical[:,:, 1] .= T(1.0);
    # calculates theta and phi
    theta = collect(range(0,stop=pi-dtheta,length=ntheta));
    phi = collect(range(0,stop=2*pi-dphi,length=nphi));
    for i = 1:ntheta
      ilong_range = Int((i-1)*nphi+1):(i*nphi);
      unit_spherical[ilong_range,1,2] .= theta[i]
      unit_spherical[ilong_range,1,3] .= phi;
      unit_spherical[ilong_range,2,2] .= theta[i] + dtheta
      unit_spherical[ilong_range,2,3] .= phi;
      unit_spherical[ilong_range,3,2] .= theta[i] + dtheta;
      unit_spherical[ilong_range,3,3] .= phi .+ dphi;
      unit_spherical[ilong_range,4,2] .= theta[i]
      unit_spherical[ilong_range,4,3] .= phi .+ dphi;
      unit_spherical[ilong_range,5,2] .= theta[i] + dtheta/2
      unit_spherical[ilong_range,5,3] .= phi .+ dphi/2;
    end
    # Setup unit vectors
    unit_xyz[:,:,1] = sin.(unit_spherical[:,:,2]).*cos.(unit_spherical[:,:,3]); # X
    unit_xyz[:,:,2] = sin.(unit_spherical[:,:,2]).*sin.(unit_spherical[:,:,3]); # Y
    unit_xyz[:,:,3] = cos.(unit_spherical[:,:,2]); # Z
    # Setup radius
    unit_spherical[:,:,1] .= T(1.0);# sqrt.(unit_xyz[:,:,1].^2 + unit_xyz[:,:,2].^2 + unit_xyz[:,:,3].^2); # radius, TBD replace by norm()
    return tessellation(1, npix, unit_xyz, unit_spherical);
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

    # make a circle on a square grid
    x = (repeat(collect(1:2*spot_radius),1,2*spot_radius) .- (2*spot_radius + 1)/2).^2
    spot = sqrt.(x + x')
    mask = findall(spot.>spot_radius)
    spot[mask] .= 0.0; # made spot

    # redefine temperature values to match temperature map
    spot_loc = findall(spot.>0.0)
    background = findall(spot .< 1e-5)
    spot[spot_loc]   .= mean(temperature_map).*bright_frac
    spot[background] .= mean(temperature_map)

    temperature_map[long-spot_radius:long+spot_radius-1,lat-spot_radius:lat+spot_radius-1] = spot
    return vec(temperature_map)
end

# maybe this isn't necessary?
function make_spot_move(temperature_map,ntheta,nphi,period,B,tepoch=[0.0,1.0])
    temperature_map = reshape(temperature_map,nphi,ntheta);
    old_temp_map = deepcopy(temperature_map);

    # calculate, based on number of ntheta, degree of central pixel for a given latitude
    central_pix = collect(range(180.0/ntheta/2.0,stop=180.0-180.0/ntheta/2.0,length=ntheta))*pi/180.0; # needs to be in radians

    long_width = 360.0/nphi; # degrees/pixel

    # Take averages of two different cells based on how much the spot has moved
    for i = 1:ntheta
        # do shift in longitude
        omega = omega_rotation(360.0/period, B, central_pix[i]);
        delta_long = omega*(tepoch[2]-tepoch[1]); # degrees
        first_pix = (delta_long - floor(delta_long/long_width)*long_width)/long_width; # percentage of first pixel covered
        second_pix = 1 - first_pix; # percentage of second pixel covered
        pix_num1 = Int(floor(delta_long/long_width)); # starting pixel that's getting changed
        pix_num2 = Int(ceil(delta_long/long_width)); # second pixel that's getting changed

        # do averages
        for j = 1:nphi
            # use circshift later?
            if (j + pix_num1 <= nphi)
                temperature_map[j+pix_num1,i] = old_temp_map[j,i]*first_pix + old_temp_map[j+1,i]*second_pix;
                #temperature_map[j+pix_num2-1,i] = old_temp_map[j,i]*second_pix + old_temp_map[j+1,i]*first_pix;
            elseif ((j + pix_num1 > nphi) & (j+1 <= nphi))
                temperature_map[j+pix_num1-nphi,i] = old_temp_map[j,i]*first_pix + old_temp_map[j+1,i]*second_pix;
                #temperature_map[j+pix_num2-nphi-1,i] = old_temp_map[j,i]*second_pix + old_temp_map[j+1,i]*first_pix;
            #elseif ((j + second_pix > nphi) & (j + first_pix > nphi) & (j+1 <= nphi))
                #println("third")
                #temperature_map[j+pix_num1-nphi,i] = old_temp_map[j,i]*first_pix + old_temp_map[j+1,i]*second_pix;
                #temperature_map[j+pix_num2-nphi-1,i] = old_temp_map[j,i]*second_pix + old_temp_map[j+1,i]*first_pix;
            else
                temperature_map[j+pix_num1-nphi,i] = old_temp_map[j,i]*first_pix + old_temp_map[j+1-nphi,i]*second_pix;
                #temperature_map[j+pix_num2-nphi,i] = old_temp_map[j,i]*second_pix + old_temp_map[j+1-nphi,i]*first_pix;
            end
            #println(j+pix_num1, " ", j)
        end
    end
    return vec(temperature_map)
end

function latlong_harmon(ntheta,nphi,a,b,c)
dtheta = pi/ntheta;
theta = collect(range(0,stop=pi-dtheta,length=ntheta));
starting = dtheta/2
#ending = pi/2 -starting;
#Int.(ceil.(90*cos.(pi/2 - collect(linspace(starting,ending,ntheta/2)))))

# Close to pixelation scheme from LCI? (Used in Harmon & Crews 2000)
arr_long = Int64[];
arr_long = vcat(arr_long,starting)
for i = 2:ntheta
    arr_long = vcat(arr_long,arr_long[i-1]+dtheta)
end
indx_long = Int.(ceil.(nphi*cos.(pi/2 .- arr_long)))
npix = sum(indx_long);
unit_spherical = zeros(Float64, npix, 5, 3); # r, theta, phi
unit_xyz = zeros(Float64, npix, 5, 3); # x, y, z
unit_spherical[:,:,1] .= 1.;

first=0
last=1
  for i = 1:ntheta
    nphi = indx_long[i];
    last = sum(indx_long[1:i])
    ilong_range=first+1:last
    dphi= 2*pi/nphi;
    phi = collect(range(0,stop=2*pi-dphi,length=nphi));
    unit_spherical[ilong_range,1,2], unit_spherical[ilong_range,1,3] .= theta[i], phi;
    unit_spherical[ilong_range,2,2], unit_spherical[ilong_range,2,3] .= theta[i] + dtheta, phi;
    unit_spherical[ilong_range,3,2], unit_spherical[ilong_range,3,3] .= theta[i] + dtheta, phi + dphi;
    unit_spherical[ilong_range,4,2], unit_spherical[ilong_range,4,3] .= theta[i], phi + dphi;
    theta_avg, phi_avg = theta[i] + dtheta*0.5, phi + dphi*0.5;
    unit_spherical[ilong_range,5,2], unit_spherical[ilong_range,5,3] .= theta_avg, phi_avg;
    first=deepcopy(last);
  end

  unit_xyz[:,:,1] = a*sin.(unit_spherical[:,:,2]).*cos.(unit_spherical[:,:,3]); # X
  unit_xyz[:,:,2] = b*sin.(unit_spherical[:,:,2]).*sin.(unit_spherical[:,:,3]); # Y
  unit_xyz[:,:,3] = c*cos.(unit_spherical[:,:,2]); # Z
  unit_spherical[:,:,1] = sqrt.(unit_xyz[:,:,1].^2 + unit_xyz[:,:,2].^2 + unit_xyz[:,:,3].^2); # radius, TBD replace by norm()

star_tessellation = tessellation(npix,unit_xyz, unit_spherical);
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
