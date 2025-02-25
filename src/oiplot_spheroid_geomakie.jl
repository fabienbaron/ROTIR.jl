using GLMakie, GeometryBasics


# #using GeoMakie, Makie.GeometryBasics
# #using Healpix
# function mollweide_makie(tmap, star)
#     T = eltype(star.vertices_xyz)
#     nside = 2^npix2n(star.npix)
#     m = HealpixMap{T, NestedOrder}(nside)
#     m.pixels[:] = tmap
#     img, _, _ = Healpix.equirectangular(m);
#     meshimage(-180..180, -90..90, transpose(img); npoints = 200, axis = (; type = GeoAxis, dest = "+proj=moll"), )
# end

function plot3d_temperature_makie(tmap, star)
npix = star.npix
x= star.vertices_xyz[:,1:4,1]'
y= star.vertices_xyz[:,1:4,2]'
z= star.vertices_xyz[:,1:4,3]'
ps = vec(Point3.(Float64.(x),Float64.(y),Float64.(z)))
fs = [QuadFace{Int64}(1+4*i,2+4*i,3+4*i,4+4*i) for i=0:size(x,2)-1]
ns = face_normals(ps, fs)
map = (tmap.-minimum(tmap))/(maximum(tmap)-minimum(tmap))
FT = eltype(fs)
cs = FaceView(RGBf(1.0,1.0,1.0).*map, [FT(i) for i in 1:npix])
m = GeometryBasics.mesh(ps, fs, normal=ns, color=cs)
mesh(m, shading=NoShading)
end

# using CairoMakie
# function plot2d_temperature_makie(tmap, star)
#     npix = star.npix
#     x = star.projx
#     y = star.projy
#     f = Figure();
#     Axis(f[1, 1])
#     ps = [ Polygon(Point2f.(x, y)[i,:]) for i=1:size(x,1)]
#     map = (tmap.-minimum(tmap))/(maximum(tmap)-minimum(tmap))
#     map = map[star.index_quads_visible]
#     poly!(ps, color = RGBf(1.0).*map)
#     return f
# end