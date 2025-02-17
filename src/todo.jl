
# # experimental stuff
# using QHull
# function hull(star_epoch_geom)
# println("Computing convex hull");
# ch = chull(hcat(vec(star_epoch_geom.projx),vec(star_epoch_geom.projy)))
# hull=[ch.vertices;ch.vertices[1]]; #convex hull, closed
# scatter(ch.points[:,1], ch.points[:,2])
# plot(ch.points[hull,1], ch.points[hull,2])
#
# bounding_box = [minimum(star_epoch_geom.projx), minimum(star_epoch_geom.projy), maximum(star_epoch_geom.projx), maximum(star_epoch_geom.projy)]
# end

# function triangle_orientation(a,b,c)
# #CrossProductZ(a,b) = a[1] * b[2] - a[2] * b[1]
# #CrossProductZ(b,c) = b[1] * c[2] - b[2] * c[1]
# #CrossProductZ(c,a) = c[1] * a[2] - c[2] * a[1]
# #  return CrossProductZ(a, b) + CrossProductZ(b, c) + CrossProductZ(c, a)
# return a[1] * b[2] - a[2] * b[1] + b[1] * c[2] - b[2] * c[1] + c[1] * a[2] - c[2] * a[1]
# end

# function sort2D_quad_counter(quad)
# #given four points, return them sorted counterclockwise
# #initial order
# a = quad[1,:];
# b = quad[2,:];
# c = quad[3,:];
# d = quad[4,:];

# if triangle_orientation(a, b, c) > 0.0
#         # Triangle abc is already clockwise.  Where does d fit?
#         if triangle_orientation(a, c, d) > 0.0
#             return quad;
#         elseif triangle_orientation(a, b, d) > 0.0
#             return quad[[1,2,4,3],:];
#         else
#           return quad[[4,2,3,1],:];
#         end
# elseif triangle_orientation(a, c, d) > 0.0
#         # Triangle abc is counterclockwise, i.e. acb is clockwise.
#         # Also, acd is clockwise.
#         if triangle_orientation(a, b, d) > 0.0
#           return quad[[1,3,2,4],:];
#         else
#           return quad[[2,1,3,4],:];
#         end
# else
#   # Triangle abc is counterclockwise, and acd is counterclockwise.
#   # Therefore, abcd is counterclockwise.
#     return quad[[3,2,1,4],:];
# end

# end
