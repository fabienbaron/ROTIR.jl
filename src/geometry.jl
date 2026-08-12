# =====================================================================
# ROTIR Coordinate Conventions
# =====================================================================
# After projection via rot_vertex and create_star:
#   proj_west  = projected coordinate pointing West  (mas)
#   proj_north = projected coordinate pointing North (mas)
#
# The interferometric sky frame has u=East, v=North (OIFITS standard).
# The polygon FT (setup_polyft_single) uses:
#   kx = u * (-π/C)    ← minus compensates West vs East
#   ky = v * ( π/C)
# so the net FT matches the standard exp(-2πi(u·East + v·North)).
#
# For plotting: East = -proj_west, North = proj_north
# =====================================================================

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
  proj_west::Array{T,2}   # Projected West coordinate of quad vertices (mas)
  proj_north::Array{T,2}  # Projected North coordinate of quad vertices (mas)
  ldmap::Array{T,1}   # Limb-darkening map
  vis_weights::Array{T,1}  # Soft visibility weights (sigmoid of normals_z), length npix
  sig_args::Array{T,1}     # Sigmoid arguments (κ * normals_z), for gradient computation
  center_offsets::Array{T,1} # Center of mass within star
  polyflux::Array{T,1}
  polyft::Matrix{Complex{T}}
  t::T # epoch time
end

# =====================================================================
# PRECISION CONVENTION
# =====================================================================
# Element types FOLLOW THE INPUTS. Never hardcode `Float32`/`Float64` in a keyword default
# or an array allocation: `tessellation_healpix(n; T=Float64)` is the root of the type
# chain, and a hardcoded `T=Float32` downstream silently narrows the whole mesh back at the
# first step. Write `T = eltype(tessels)` / `float(real(eltype(x)))` instead.
#
# Float32 is the right default for interferometry — the observables carry ~1 % errors, and
# the orbit forward model costs ~1e-6 mas in predicted separation at Float32, five orders
# below CHARA astrometry. There are exactly four documented exceptions:
#
#  1. ABSOLUTE MJD/JD. `eps(Float32(2.45e6)) = 0.25 d`. Epoch arguments are annotated
#     `::Float64` (`binary_orbit_abs`, `binary_RV`) or forced with `Float64(tepoch_jd)`.
#     Time DIFFERENCES are fine in Float32 (21 s over a 2450 d span ⇒ 1e-4 mas).
#  2. SUMMED SCALARS. Per-element values take the data's type; a scalar accumulated over the
#     whole mesh or dataset takes the wider one. Julia's `sum` is pairwise so the summation
#     itself is near-exact — the issue is REPRESENTING a large result: `eps(Float32(2e6))`
#     is 0.125, against a Δχ² = 1 confidence contour. See `betlyr/betlyr_model.jl`.
#  3. GEOMETRIC PREDICATES. Sutherland–Hodgman clipping decides inside/outside from the
#     sign of a cross product; a Float32 sign flip changes vertex COMBINATORICS, not just
#     the area. Exact predicates on inexact input — see `occultation_weights`.
#  4. ROOT SOLVES / QUADRATURE with a nested Halley stack (`roche_volume`, `roche_area`,
#     `romberg_integrate`, `solve_R_L2/L3`), which need the headroom. Their docstrings
#     say so explicitly.
#
# Element type accessors, so callers write `T = eltype(tessels)` rather than a literal.
Base.eltype(::tessellation{T}) where {T} = T
Base.eltype(::Type{tessellation{T}}) where {T} = T
Base.eltype(::stellar_geometry{T}) where {T} = T
Base.eltype(::Type{stellar_geometry{T}}) where {T} = T

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
  println("proj_west             : projected West vertex coordinates (mas)")
  println("proj_north            : projected North vertex coordinates (mas)")
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
#   proj_west::Array{Float64,2}
#   proj_north::Array{Float64,2}
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
include("geometry_ellipsoid.jl");
# ZXZ Euler rotation matrix: R = R_z(-a1) * R_x(a2) * R_z(-a3)
# For rotate_star (right multiply: xyz * R):
#   angle_r1 = ψ (rotation phase), angle_r2 = inclination, angle_r3 = PA
#   Spin axis = R[3,:] = [-sin(inc)*sin(PA), sin(inc)*cos(PA), cos(inc)]
#   i.e. the projected pole is at position angle PA from North toward East.
function rot_vertex(angle_r1, angle_r2, angle_r3)
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
  

"""
    compute_radii(tessels, star_params, t; secondary=false, T=eltype(tessels), D=nothing, omega=nothing)

Vertex radii and Cartesian positions in the star's *body* frame.

For a Roche surface (`surface_type == 3`) the shape depends on the instantaneous
separation, normally `D = compute_separation(star_params, t)` in units of the semi-major
axis. Pass `D` explicitly to override that — `create_binary_star` does, because in a binary
the separation comes from the shared orbit rather than from this component's own `t`.
Pass `omega` to fix the surface equipotential directly instead of deriving it from `rpole`
(used by the volume-conserving path, see `roche_omega_for_volume`).
"""
function compute_radii(tessels::tessellation, star_params, t; secondary=false, T=eltype(tessels),
                       D=nothing, omega=nothing)
  p = convert_params(T, star_params)
  npix = tessels.npix
  xyz = [];
  r = [];
  # compute radii and xyz based on stellar parameters
  if p.surface_type  == 0  # Spherical:0, , Rapid Rotator:2, Roche: 3
    xyz = p.radius*tessels.unit_xyz;
    r = repeat([p.radius],npix, 5);
  elseif p.surface_type == 1 # Ellipsoid: 1
    xyz = reshape(reshape(tessels.unit_xyz, (npix*5,3)).*T.([p.radius_x p.radius_y p.radius_z]), (npix,5,3));;
    # TODO: just repeat of radius_x, radius_y and radius_z then block multiply
    r = sqrt.(sum(xyz.^2, dims=3));
  elseif p.surface_type == 2 # Rapid rotator
    r = update_radii_rapidrot(tessels, p);
    xyz = r.*tessels.unit_xyz;
  elseif p.surface_type == 3
    # Star params are actually binary parameters
    Duse = D === nothing ? T(compute_separation(p, t)) : T(D)
    ff = secondary ? p.fillout_factor_secondary : p.fillout_factor_primary
    r = update_roche_radii(tessels, p, Duse, use_fillout_factor = ff>-1, secondary=secondary,
                           omega=omega, T=T)
    xyz = r.*tessels.unit_xyz;
  end
  return r, xyz
end

function rotate_star(xyz, star_params, t; T = float(real(eltype(xyz))))
  # TODO: reimplement differential rotation for compatible surfaces (see old ROTIR)
  # Right multiply: xyz * R where R = rot_vertex(ψ, inc, PA)
  #   ψ  = rotation phase (spinning surface features)
  #   inc = inclination (tilt spin axis from line of sight)
  #   PA  = position angle (orient projected spin axis on sky, N through E)
  # Spin axis direction: R[3,:] = [-sin(inc)*sin(PA), sin(inc)*cos(PA), cos(inc)]
  p = convert_params(T, star_params)
  npix = size(xyz,1)
  compound_rotation = rot_vertex(T(2pi)*t/p.rotation_period, p.inclination*T(pi/180.), p.position_angle*T(pi/180.));
  return reshape(reshape(xyz,(npix*5,3))*compound_rotation, (npix,5,3));
end


# Limb cosine μ from the line-of-sight normal component nz (single source of truth,
# shared by the forward model in create_star and by the LD gradient rules).
# The face normal is unit-length (see create_star), so nz is exactly cos(emergent
# angle) = the limb cosine; visible-side μ = max(nz, 0), back-facing tessels clamp to 0.
@inline limb_mu(nz::T) where {T} = nz > zero(T) ? nz : zero(T)
@inline mu_and_dmu(nz::T) where {T} = nz > zero(T) ? (nz, one(T)) : (zero(T), zero(T))

function compute_ldmap(μ, star_params; T = float(real(eltype(μ))))
  p = convert_params(T, star_params)
  # Limb-darkening map
  if (p.ldtype == 1) # 1: quadratic
    ldmap = T(1.0) .- p.ld1*(T(1.0) .-μ)
  elseif (p.ldtype == 2) # 2: quadratic
    ldmap = T(1.0) .- p.ld1*(T(1.0) .-μ) - p.ld2*(T(1.0).-μ.^2)
  elseif (p.ldtype == 3)  # 3; Hestroffer
    ldmap = μ.^p.ld1
  end
end

# Generate geometry and ld map from tesselation and stellar parameters
@views function create_star(tessels::tessellation, star_params, t; secondary=false, T=eltype(tessels), κ=T(50))
  # Compute radii
  r, xyz = compute_radii(tessels, star_params, t; secondary=secondary);

  # Compute rotation
  xyz = rotate_star(xyz, star_params, t);

  return finish_star(xyz, r, tessels, star_params, t; T=T, κ=κ)
end

"""
    finish_star(xyz, r, tessels, star_params, t; T=eltype(tessels), κ=T(50)) -> stellar_geometry

Common tail of the geometry pipeline: everything that follows the choice of vertex
positions. Given already-deformed, already-rotated sky-frame vertices `xyz` (npix,5,3)
and their radii `r` (npix,5), compute face normals, soft visibility, the orthographic
projection and the limb-darkening map, and package them into a `stellar_geometry`.

Split out of [`create_star`](@ref) so that alternative orientations — notably
`create_binary_star`, which rotates into the *orbital* frame instead of spinning about
a fixed spin axis — reuse this verbatim rather than duplicating it.
"""
@views function finish_star(xyz, r, tessels::tessellation, star_params, t; T=eltype(tessels), κ=T(50))
  npix = tessels.npix;

  # Determine normals via cross product of quad diagonals
  vecAC = xyz[:, 3, :]-xyz[:, 1, :];
  vecBD = xyz[:, 4, :]-xyz[:, 2, :];
  normals_tmp = [ vecAC[:,2].*vecBD[:,3] - vecAC[:,3].*vecBD[:,2] vecAC[:,3].*vecBD[:,1] - vecAC[:,1].*vecBD[:,3] vecAC[:,1].*vecBD[:,2] - vecAC[:,2].*vecBD[:,1]];
  normals = normals_tmp./sqrt.(sum(abs2, normals_tmp, dims=2))

  # Soft visibility: sigmoid weights replace hard mask
  vis_weights, sig_args = soft_visibility(normals[:,3], κ=κ)

  # Hard index kept for backward compat (pixels with significant visibility)
  vis_threshold = T(0.01)
  index_quads_visible = findall(vis_weights .> vis_threshold);
  nquads_visible = length(index_quads_visible);

  # Project ALL pixels onto the observing plane (needed for soft visibility gradients)
  # proj_west = xyz component 1 (West on sky), proj_north = xyz component 2 (North on sky)
  proj_west = xyz[:, 1:4, 1];
  proj_north = xyz[:, 1:4, 2];

  # Limb-darkening map. μ = limb cosine max(nz,0) on the visible side (see limb_mu).
  μ = limb_mu.(normals[:,3])
  ldmap = compute_ldmap(μ,star_params)
  spherical = copy(tessels.unit_spherical);
  spherical[:,:,1] = r
  # Single star
  center = T.([0.0,0.0,0.0]);
  return stellar_geometry{T}(star_params.surface_type, tessels.tessellation_type, npix, xyz, spherical, normals, index_quads_visible, nquads_visible, proj_west, proj_north, ldmap, vis_weights, sig_args, center, T[], zeros(Complex{T}, 0, 0), t);
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

function invisible_neighbors(n, stars)
  # all invisible tessels that neighbor visible tessels
  # this should be all tessels just beyond the limb 
  # they will affect TV values!
  nlist = tv_neighbors_healpix(n)[1]
  return intersect(unique(vcat(nlist[sometimes_visible(stars)]...)), never_visible(stars))
end

function with_invisible_neighbors(n, stars)
  # all visible tessels who have at least an invisible neighbor
  # they should be the limb tessels
  nlist = tv_neighbors_healpix(n)[1]
  return intersect(unique(vcat(nlist[invisible_neighbors(n, stars)]...)),sometimes_visible(stars) )
  return 
end

function without_invisible_neighbors(n, stars)
  # all visible tessels who have at least an invisible neighbor
  # they should be the limb tessels
  nlist = tv_neighbors_healpix(n)[1]
  npix= length(nlist)
  invisible = invisible_neighbors(n, stars)
  return findall(.!prod.([nlist[i] .∉ Ref(invisible) for i=1:npix]))
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

