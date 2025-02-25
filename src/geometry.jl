mutable struct tessellation{T}
  tessellation_type::Int64 # 0: Healpix, 1: Longitude/Latitude
  npix::Int64
  unit_xyz::Array{T,3}
  unit_spherical::Array{T,3}
end

mutable struct stellar_geometry{T} # typically one per epoch, rotation and projection of the base geometry
  surface_type::Int64
  tessellation_type::Int64 # 0: Healpix, 1: Longitude/Latitude # inherited from setup
  npix::Int64
  vertices_xyz::Array{T,3}
  vertices_spherical::Array{T,3}
  normals::Array{T,2}
  index_quads_visible::Array{Int64,1}
  nquads_visible::Int64
  projx::Array{T,2}
  projy::Array{T,2}
  ldmap::Array{T,1}   # Limb-darkening map
  center_offsets::Array{T,1} # Center of mass within star
  polyflux::Array{T,1}
  polyft::Matrix{Complex{T}}
  t::T # epoch time
end

function Base.display(x::tessellation)
  if x.tessellation_type==0
    println("Tessellation type: Healpix")
  elseif x.tessellation_type==1
    println("Tessellation type: Latitude-Longitude")
  else 
    println("Unknown tessellation type");
  end
  println("Number of tessels = $(x.npix)")
end

function Base.display(x::stellar_geometry)
  if x.surface_type==0
    println("Surface type: Sphere")
  elseif x.surface_type==1
    println("Surface type: Ellipsoid")
  elseif x.surface_type==2
    println("Surface type: Rapid Rotator")
  elseif x.surface_type==3
    println("Surface type: Roche surface")
  else
    println("Unknown Surface type");
  end
  if x.tessellation_type==0
    println("Tessellation type: Healpix")
  elseif x.tessellation_type==1
    println("Tessellation type: Latitude-Longitude")
  else 
    println("Unknown tessellation type");
  end
  println("Number of tessels = $(x.npix)")
  println("nquads_visible = $(x.nquads_visible)")
  println("Other fields:")
  println("--------------------------------------------------")
  println("index_quads_visible   : list of the visible tessels")
  println("vertices_xyz          : (x,y,z) coordinates of the vertices")
  println("vertices_spherical    : (r,θ,ϕ) coordinates of the vertices")
  println("normals               : coordinates of the vertex normals")
  println("projx                 : projected x vertex coordinates")
  println("projy                 : projected y vertex coordinates") 
  println("ldmap                 : limb-darkening map")
  if x.polyft == []
  println("polyflux              : temperature to flux vector (not defined yet)")
  println("polyft                : temperature to visibility matrix (not defined yet)")
  else
  println("polyflux              : temperature to flux vector (set)")
  println("polyft                : temperature to visibility matrix (set)")
  end
  println("Epoch (time)          : time corresponding to this stellar shape")
end


function Base.display(x::Array{stellar_geometry,1})
  println("Array of Stellar geometries - Healpix")
  println("Number of epochs = $(length(x))")
  println("npix = $(x[1].npix)")
end

# mutable struct tessellation  # this is only affected by stellar parameters (e.g. Roche parameters)
#   npix::Int64
#   # Healpix vertices
#   vertices_xyz::Array{Float64,3}
#   vertices_spherical::Array{Float64,3}
# end

# Base.iterate(star_tessellation::tessellation, state = 1) = state <= 1 ? (star_tessellation, state+1) : nothing
# Base.length(star_tessellation::tessellation) = min(length(star_tessellation.npix), length(star_tessellation.vertices_xyz), length(star_tessellation.vertices_spherical))

# mutable struct stellar_geometry # typically one per epoch, rotation and projection of the base geometry
#   npix::Int64
#   # Healpix vertices
#   vertices_xyz::Array{Float64,3}
#   vertices_spherical::Array{Float64,3}
#   normals::Array{Float64,2}
#   # Healpix projection onto the 2D imaging plane
#  # quads_visible::Array{Bool,1} not necessary anymore
#   index_quads_visible::Array{Int64,1}
#   nquads_visible::Int64
#   projx::Array{Float64,2}
#   projy::Array{Float64,2}
#   # Limb-darkening map
#   ldmap::Array{Float64,1}
#   # Center of mass within star
#   offsets::Array{Float64,1}
# end

# Base.iterate(star_epoch_geom::stellar_geometry, state = 1) = state <= 1 ? (star_epoch_geom, state+1) : nothing
# Base.length(star_epoch_geom::stellar_geometry) = 1

include("tessellation_healpix.jl");
include("tessellation_latlong.jl");
include("geometry_rochelobe.jl");
include("geometry_rapidrotator.jl")

function rot_vertex(angle_r1, angle_r2, angle_r3) # new rotational matrix
  c1 = cos(angle_r1)
  s1 = sin(angle_r1)
  c2 = cos(angle_r2)
  s2 = sin(angle_r2)
  c3 = cos(angle_r3)
  s3 = sin(angle_r3)
  dcm = [-s1*c2*s3+c1*c3  s1*c3*c2+c1*s3 -s1*s2;
         -c1*c2*s3-s1*c3  c1*c3*c2-s1*s3 -c1*s2 ;
                  -s2*s3           s2*c3     c2 ];
  return dcm
end
  

function compute_radii(tessels::tessellation, star_params, t; T=Float32) 
  npix = tessels.npix
  xyz = [];
  r = [];
  # compute radii and xyz based on stellar parameters
  if star_params.surface_type  == 0  # Spherical:0, , Rapid Rotator:2, Roche: 3
    xyz = star_params.radius*tessels.unit_xyz;
    r = repeat([star_params.radius],npix, 5);
  elseif star_params.surface_type == 1 # Ellipsoid: 1
    xyz = reshape(reshape(tessels.unit_xyz, (npix*5,3)).*[star_params.radius_x star_params.radius_y star_params.radius_z], (npix,5,3));;
    # TODO: just repeat of radius_x, radius_y and radius_z then block multiply
    r = sqrt.(sum(xyz.^2, dims=3));
  elseif star_params.surface_type == 2 # Rapid rotator
    r = update_radii_rapidrot(tessels, star_params);
    xyz = r.*tessels.unit_xyz;
  elseif star_params.surface_type == 3
    # Star params are actually binary parameters
    D = T(compute_separation(star_params, t))
    r = update_roche_radii(tessels, star_params, D, use_fillout_factor = star_params.fillout_factor_primary>-1) 
    xyz = r.*tessels.unit_xyz;
  end
  return r, xyz
end

function rotate_star(xyz, star_params, t; T=Float32)
  # TODO: reimplement differential rotation for compatible surfaces (see old ROTIR)
  npix = size(xyz,1)
  compound_rotation = rot_vertex( T(2pi)*t/star_params.rotation_period, star_params.inclination*T(pi/180.), star_params.position_angle*T(pi/180.));
  return reshape(reshape(xyz,(npix*5,3))*compound_rotation, (npix,5,3));
end


function compute_ldmap(μ,star_params; T=Float32)
  # Limb-darkening map
  if (star_params.ldtype == 1) # 1: quadratic
    ldmap = T(1.0) .- star_params.ld1*(T(1.0) .-μ) 
  elseif (star_params.ldtype == 2) # 1: quadratic
    ldmap = T(1.0) .- star_params.ld1*(T(1.0) .-μ) - star_params.ld2*(T(1.0).-μ.^2)  
  elseif (star_params.ldtype == 3)  # 3; Hestroffer
    ldmap = μ.^star_params.ld1
  end
end

# Generate geometry and ld map from tesselation and stellar parameters
@views function create_star(tessels::tessellation, star_params, t; secondary=false, T=Float32)
  npix = tessels.npix;
  # Compute radii 
  r, xyz = compute_radii(tessels, star_params, t);
  
  # Compute rotation
  xyz = rotate_star(xyz, star_params, t); 

  # Determine visible pixels based on whether their normals point toward us (z>0)
  vecAC = xyz[:, 3, :]-xyz[:, 1, :];
  vecBD = xyz[:, 4, :]-xyz[:, 2, :];
  normals_tmp = [ vecAC[:,2].*vecBD[:,3] - vecAC[:,3].*vecBD[:,2] vecAC[:,3].*vecBD[:,1] - vecAC[:,1].*vecBD[:,3] vecAC[:,1].*vecBD[:,2] - vecAC[:,2].*vecBD[:,1]];  
  normals = normals_tmp./sqrt.(sum(abs2, normals_tmp, dims=2))
  index_quads_visible = findall(normals[:,3].>0);
  nquads_visible = length(index_quads_visible);
  # projx = vertices_xyz[:, 1:4, 1]; 
  # projy = vertices_xyz[:, 1:4, 2];
  # μ = normals[:,3].*max.(normals[:,3], 0.)
  # 2D the projection onto the (x,y) observing plane
  projx = xyz[index_quads_visible, 1:4, 1];
  projy = xyz[index_quads_visible, 1:4, 2];
  # Limb-darkening map
  μ = abs.(normals[:,3].*max.(normals[:,3], 0))
  ldmap = compute_ldmap(μ,star_params)
  spherical = copy(tessels.unit_spherical);
  spherical[:,:,1] = r
  # Single star
  center = T.([0.0,0.0,0.0]);
  return stellar_geometry{T}(star_params.surface_type, tessels.tessellation_type, npix, xyz, spherical, normals, index_quads_visible,  nquads_visible, projx,  projy, ldmap, center, T[], zeros(Complex{T}, 0, 0), t);
end

function create_binary(star1::tessellation, star2::tessellation, binary_params::binaryparameters, t)
# Update a binary system
# Step 1: compute location of components
# Step 2: compute Roche lobes for both components
 
## TEST star1=prim_base; star2=sec_base;t = tepochs[i]
  x1, y1, z1, x2, y2, z2 = binary_orbit_abs(binary_params,t);
  D = sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2)/binary_params.a # reduced separation

  # Update Roche geometry
  star1_roche_geom, star2_roche_geom =  update_roche_radii(star1, star2, binary_params, D); # updates both primary and secondary

  # Determine rotation angles
  ζ1 = 360.0/binary_params.star1.rotation_period*(t-binary_params.T0) + 180.0 
  ζ2 = 360.0/binary_params.star2.rotation_period*(t-binary_params.T0)
  Ω = binary_params.Ω*pi/180.0; # longitude of ascending node
  i = binary_params.i*pi/180.0;

  # Rotate and translate both stars
  star_epoch_geom1 = rotate_single_star(star1_roche_geom, binary_params.star1, i*pi/180.0, Ω*pi/180.0, ζ1*pi/180.0, offsets = [x1, z1, y1])
  star_epoch_geom2 = rotate_single_star(star2_roche_geom, binary_params.star2, i*pi/180.0, Ω*pi/180.0, ζ2*pi/180.0, offsets = [x2, z2, y2], recenter=[-D, 0, 0])
  return star_epoch_geom1, star_epoch_geom2
end

# The following function is used for generating stars for a binary
function rotate_single_star(base_geom, star_params, orbit_incl, long_ascending_node, rot_angle; offsets = [], recenter=[])
  npix = base_geom.npix;
  compound_rotation = rot_vertex(orbit_incl*pi/180.0, long_ascending_node*pi/180.0, rot_angle*pi/180.0);
  vertices_xyz_rot = Array{Float64}(undef, npix, 5, 3)
  if recenter ==[]
    for ii=1:npix
      for jj = 1:5
       vertices_xyz_rot[ii, jj, :] = compound_rotation * base_geom.vertices_xyz[ii, jj, :] 
      end
    end
  else 
    for ii=1:npix
      for jj = 1:5
        vertices_xyz_rot[ii, jj, :] = compound_rotation * (base_geom.vertices_xyz[ii, jj, :] + recenter) 
      end
    end
  end
  # Offsets applied after rotation
  if offsets !=[]
    vertices_xyz_rot[:,:,1] .+= offsets[1]
    vertices_xyz_rot[:,:,2] .+= offsets[2]
    vertices_xyz_rot[:,:,3] .+= offsets[3]
  end

  # TODO: Compute the rest as usual including normals, mu, ldmap ,visible parts
  vertices_spherical = copy(base_geom.vertices_spherical)
  return stellar_geometry(npix, vertices_xyz_rot, vertices_spherical, normals, index_quads_visible,  nquads_visible, projx,  projy, ldmap, offsets);
 end

function never_visible(star_epoch_geom)
  # return the list of pixels which are never visible at any epochs
  hidden = Array{Bool}(undef, star_epoch_geom[1].npix)
  hidden[:] .= true;
  for t=1:length(star_epoch_geom)
  hidden[star_epoch_geom[t].index_quads_visible] .= false;
  end
  return findall(hidden.==true)
end

function sometimes_visible(star_epoch_geom)
  # return the list of pixels which are at least visible during one epoch
  sometimes = Array{Bool}(undef, star_epoch_geom[1].npix)
  sometimes[:] .= false;
  for t=1:length(star_epoch_geom)
  sometimes[star_epoch_geom[t].index_quads_visible] .= true;
  end
  return findall(sometimes.==true)
end

function create_star_multiepochs(tessels::tessellation, star_params, tepochs; kwargs...)
nepochs = length(tepochs);
npix = tessels.npix
star_epoch_geom = Array{stellar_geometry}(undef, nepochs);
#println("Creating geometry for $(nepochs) epochs x $(npix) tessels");
for i=1:nepochs
  star_epoch_geom[i] = create_star(tessels, star_params, tepochs[i]; kwargs...);
end
return star_epoch_geom
end

function create_binary_geometry(tessels, binary_parameters, tepochs; kwargs...)
  nepochs = length(tepochs);
  binary_epoch_geom = Array{stellar_geometry}(undef, nepochs);
  #offsets = Array{Float64}(undef,nepochs,3);
  for i=1:nepochs
    binary_epoch_geom[i] = update_binary(tessels, binary_parameters,tepochs[i]; kwargs...);
  end
  return binary_epoch_geom
end
