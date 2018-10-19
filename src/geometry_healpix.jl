# Basic Healpix functions, painstakingly adapted from IDL
function nside2npix(nside)
return 12 * nside^2
end

function npix2nside(npix::Integer)
    (npix % 12 == 0) || throw("Invalid number of pixels")
    square_root = sqrt(npix / 12)
    (square_root^2 == npix / 12) || throw("Invalid number of pixels")
    return Int(round(square_root))
end

function sub_compute_vertices(z, z_nv, z_sv, phi, phi_nv, phi_sv, hdelta_phi)
# function called by pix2vec_nest to compute the
# 3D vectors pointing toward the vertices of a pixel from their angles
sth = sqrt.((1.0.-z).*(1.0.+z))
sth_nv = sqrt.((1.0.-z_nv).*(1.0.+z_nv))
phi_wv = phi - hdelta_phi
sth_sv = sqrt.((1.0.-z_sv).*(1.0.+z_sv))
phi_ev = phi + hdelta_phi

vertex = Array{Float64}(undef, length(z), 3, 4);
vertex[:,:,1] = hcat(sth_nv.*cos.(phi_nv),sth_nv.*sin.(phi_nv),z_nv); # north vertex
vertex[:,:,2] = hcat(sth.*cos.(phi_wv),sth.*sin.(phi_wv),z); # west vertex
vertex[:,:,3] = hcat(sth_sv.*cos.(phi_sv),sth_sv.*sin.(phi_sv),z_sv); # south vertex
vertex[:,:,4] = hcat(sth.*cos.(phi_ev),sth.*sin.(phi_ev),z); # east vertex
#
return vertex
end


function vec2ang(vector::Array{Float64,2};astro = false)
#        * if ASTRO is NOT set (default) : geometric coordinates
#       Theta : colatitude in RADIAN measured Southward from North pole
#       Phi   : longitude in RADIAN, increasing Eastward
#        * if ASTRO is set : astronomic coordinates
#       Theta : latitude in DEGREE measured Northward from Equator
#       Phi   : longitude in DEGREE, increasing Eastward

# vector is N x 3
theta_rad = atan.( sqrt.(vector[:,1].^2+vector[:,2].^2) , vector[:,3]) # in [0,Pi]
phi_rad   = atan.( vector[:,2], vector[:,1] )  # in [-Pi,Pi]
phi_rad = phi_rad + (2*pi).* (phi_rad .< 0.)
if astro == true
    theta = 90.0 .- theta_rad .*(180.0/pi)
    phi   = phi_rad .* (180.0/pi)
else
    theta = theta_rad
    phi = phi_rad
end
return theta, phi
end


function swapLSBMSB(i)
# Returns i with even and odd bit positions interchanged.
oddbits  = 89478485         # 2^0 + 2^2 + 2^4+..+2^26
evenbits = 178956970         # 2^1 + 2^3 + 2^4+..+2^27
li = Int64(i)
swapLSBMSB = div(li & evenbits, 2) + (li & oddbits)*2
return swapLSBMSB
end

function invLSBMSB(i)
# Returns NOT(i)
invLSBMSB = ~(Int64(i))
return invLSBMSB
end

function invswapLSBMSB(i)
# Returns NOT(i) with even and odd bit positions interchanged.
return ~swapLSBMSB(i)
end

function invLSB(i)
# Returns i with odd (1,3,5,...) bits inverted.
oddbits = 89478485
invLSB = xor(Int64(i), oddbits);
return invLSB
end

function invMSB(i)
# Returns i with even (0,2,4,...) bits inverted.
evenbits=178956970
invMSB = xor(Int64(i), evenbits)
return invMSB
end

#first setup arrays
const x2pix = [
       0,     1,     4,     5,    16,    17,    20,    21,    64,    65,
      68,    69,    80,    81,    84,    85,   256,   257,   260,   261,
     272,   273,   276,   277,   320,   321,   324,   325,   336,   337,
     340,   341,  1024,  1025,  1028,  1029,  1040,  1041,  1044,  1045,
    1088,  1089,  1092,  1093,  1104,  1105,  1108,  1109,  1280,  1281,
    1284,  1285,  1296,  1297,  1300,  1301,  1344,  1345,  1348,  1349,
    1360,  1361,  1364,  1365,  4096,  4097,  4100,  4101,  4112,  4113,
    4116,  4117,  4160,  4161,  4164,  4165,  4176,  4177,  4180,  4181,
    4352,  4353,  4356,  4357,  4368,  4369,  4372,  4373,  4416,  4417,
    4420,  4421,  4432,  4433,  4436,  4437,  5120,  5121,  5124,  5125,
    5136,  5137,  5140,  5141,  5184,  5185,  5188,  5189,  5200,  5201,
    5204,  5205,  5376,  5377,  5380,  5381,  5392,  5393,  5396,  5397,
    5440,  5441,  5444,  5445,  5456,  5457,  5460,  5461 ];

################################################################################

const y2pix = [
        0,     2,     8,    10,    32,    34,    40,    42,   128,   130,
      136,   138,   160,   162,   168,   170,   512,   514,   520,   522,
      544,   546,   552,   554,   640,   642,   648,   650,   672,   674,
      680,   682,  2048,  2050,  2056,  2058,  2080,  2082,  2088,  2090,
     2176,  2178,  2184,  2186,  2208,  2210,  2216,  2218,  2560,  2562,
     2568,  2570,  2592,  2594,  2600,  2602,  2688,  2690,  2696,  2698,
     2720,  2722,  2728,  2730,  8192,  8194,  8200,  8202,  8224,  8226,
     8232,  8234,  8320,  8322,  8328,  8330,  8352,  8354,  8360,  8362,
     8704,  8706,  8712,  8714,  8736,  8738,  8744,  8746,  8832,  8834,
     8840,  8842,  8864,  8866,  8872,  8874, 10240, 10242, 10248, 10250,
    10272, 10274, 10280, 10282, 10368, 10370, 10376, 10378, 10400, 10402,
    10408, 10410, 10752, 10754, 10760, 10762, 10784, 10786, 10792, 10794,
    10880, 10882, 10888, 10890, 10912, 10914, 10920, 10922 ];

################################################################################

const pix2x = [
  0,  1,  0,  1,  2,  3,  2,  3,  0,  1,  0,  1,  2,  3,  2,  3,
  4,  5,  4,  5,  6,  7,  6,  7,  4,  5,  4,  5,  6,  7,  6,  7,
  0,  1,  0,  1,  2,  3,  2,  3,  0,  1,  0,  1,  2,  3,  2,  3,
  4,  5,  4,  5,  6,  7,  6,  7,  4,  5,  4,  5,  6,  7,  6,  7,
  8,  9,  8,  9, 10, 11, 10, 11,  8,  9,  8,  9, 10, 11, 10, 11,
 12, 13, 12, 13, 14, 15, 14, 15, 12, 13, 12, 13, 14, 15, 14, 15,
  8,  9,  8,  9, 10, 11, 10, 11,  8,  9,  8,  9, 10, 11, 10, 11,
 12, 13, 12, 13, 14, 15, 14, 15, 12, 13, 12, 13, 14, 15, 14, 15,
  0,  1,  0,  1,  2,  3,  2,  3,  0,  1,  0,  1,  2,  3,  2,  3,
  4,  5,  4,  5,  6,  7,  6,  7,  4,  5,  4,  5,  6,  7,  6,  7,
  0,  1,  0,  1,  2,  3,  2,  3,  0,  1,  0,  1,  2,  3,  2,  3,
  4,  5,  4,  5,  6,  7,  6,  7,  4,  5,  4,  5,  6,  7,  6,  7,
  8,  9,  8,  9, 10, 11, 10, 11,  8,  9,  8,  9, 10, 11, 10, 11,
 12, 13, 12, 13, 14, 15, 14, 15, 12, 13, 12, 13, 14, 15, 14, 15,
  8,  9,  8,  9, 10, 11, 10, 11,  8,  9,  8,  9, 10, 11, 10, 11,
 12, 13, 12, 13, 14, 15, 14, 15, 12, 13, 12, 13, 14, 15, 14, 15,
 16, 17, 16, 17, 18, 19, 18, 19, 16, 17, 16, 17, 18, 19, 18, 19,
 20, 21, 20, 21, 22, 23, 22, 23, 20, 21, 20, 21, 22, 23, 22, 23,
 16, 17, 16, 17, 18, 19, 18, 19, 16, 17, 16, 17, 18, 19, 18, 19,
 20, 21, 20, 21, 22, 23, 22, 23, 20, 21, 20, 21, 22, 23, 22, 23,
 24, 25, 24, 25, 26, 27, 26, 27, 24, 25, 24, 25, 26, 27, 26, 27,
 28, 29, 28, 29, 30, 31, 30, 31, 28, 29, 28, 29, 30, 31, 30, 31,
 24, 25, 24, 25, 26, 27, 26, 27, 24, 25, 24, 25, 26, 27, 26, 27,
 28, 29, 28, 29, 30, 31, 30, 31, 28, 29, 28, 29, 30, 31, 30, 31,
 16, 17, 16, 17, 18, 19, 18, 19, 16, 17, 16, 17, 18, 19, 18, 19,
 20, 21, 20, 21, 22, 23, 22, 23, 20, 21, 20, 21, 22, 23, 22, 23,
 16, 17, 16, 17, 18, 19, 18, 19, 16, 17, 16, 17, 18, 19, 18, 19,
 20, 21, 20, 21, 22, 23, 22, 23, 20, 21, 20, 21, 22, 23, 22, 23,
 24, 25, 24, 25, 26, 27, 26, 27, 24, 25, 24, 25, 26, 27, 26, 27,
 28, 29, 28, 29, 30, 31, 30, 31, 28, 29, 28, 29, 30, 31, 30, 31,
 24, 25, 24, 25, 26, 27, 26, 27, 24, 25, 24, 25, 26, 27, 26, 27,
 28, 29, 28, 29, 30, 31, 30, 31, 28, 29, 28, 29, 30, 31, 30, 31,
  0,  1,  0,  1,  2,  3,  2,  3,  0,  1,  0,  1,  2,  3,  2,  3,
  4,  5,  4,  5,  6,  7,  6,  7,  4,  5,  4,  5,  6,  7,  6,  7,
  0,  1,  0,  1,  2,  3,  2,  3,  0,  1,  0,  1,  2,  3,  2,  3,
  4,  5,  4,  5,  6,  7,  6,  7,  4,  5,  4,  5,  6,  7,  6,  7,
  8,  9,  8,  9, 10, 11, 10, 11,  8,  9,  8,  9, 10, 11, 10, 11,
 12, 13, 12, 13, 14, 15, 14, 15, 12, 13, 12, 13, 14, 15, 14, 15,
  8,  9,  8,  9, 10, 11, 10, 11,  8,  9,  8,  9, 10, 11, 10, 11,
 12, 13, 12, 13, 14, 15, 14, 15, 12, 13, 12, 13, 14, 15, 14, 15,
  0,  1,  0,  1,  2,  3,  2,  3,  0,  1,  0,  1,  2,  3,  2,  3,
  4,  5,  4,  5,  6,  7,  6,  7,  4,  5,  4,  5,  6,  7,  6,  7,
  0,  1,  0,  1,  2,  3,  2,  3,  0,  1,  0,  1,  2,  3,  2,  3,
  4,  5,  4,  5,  6,  7,  6,  7,  4,  5,  4,  5,  6,  7,  6,  7,
  8,  9,  8,  9, 10, 11, 10, 11,  8,  9,  8,  9, 10, 11, 10, 11,
 12, 13, 12, 13, 14, 15, 14, 15, 12, 13, 12, 13, 14, 15, 14, 15,
  8,  9,  8,  9, 10, 11, 10, 11,  8,  9,  8,  9, 10, 11, 10, 11,
 12, 13, 12, 13, 14, 15, 14, 15, 12, 13, 12, 13, 14, 15, 14, 15,
 16, 17, 16, 17, 18, 19, 18, 19, 16, 17, 16, 17, 18, 19, 18, 19,
 20, 21, 20, 21, 22, 23, 22, 23, 20, 21, 20, 21, 22, 23, 22, 23,
 16, 17, 16, 17, 18, 19, 18, 19, 16, 17, 16, 17, 18, 19, 18, 19,
 20, 21, 20, 21, 22, 23, 22, 23, 20, 21, 20, 21, 22, 23, 22, 23,
 24, 25, 24, 25, 26, 27, 26, 27, 24, 25, 24, 25, 26, 27, 26, 27,
 28, 29, 28, 29, 30, 31, 30, 31, 28, 29, 28, 29, 30, 31, 30, 31,
 24, 25, 24, 25, 26, 27, 26, 27, 24, 25, 24, 25, 26, 27, 26, 27,
 28, 29, 28, 29, 30, 31, 30, 31, 28, 29, 28, 29, 30, 31, 30, 31,
 16, 17, 16, 17, 18, 19, 18, 19, 16, 17, 16, 17, 18, 19, 18, 19,
 20, 21, 20, 21, 22, 23, 22, 23, 20, 21, 20, 21, 22, 23, 22, 23,
 16, 17, 16, 17, 18, 19, 18, 19, 16, 17, 16, 17, 18, 19, 18, 19,
 20, 21, 20, 21, 22, 23, 22, 23, 20, 21, 20, 21, 22, 23, 22, 23,
 24, 25, 24, 25, 26, 27, 26, 27, 24, 25, 24, 25, 26, 27, 26, 27,
 28, 29, 28, 29, 30, 31, 30, 31, 28, 29, 28, 29, 30, 31, 30, 31,
 24, 25, 24, 25, 26, 27, 26, 27, 24, 25, 24, 25, 26, 27, 26, 27,
 28, 29, 28, 29, 30, 31, 30, 31, 28, 29, 28, 29, 30, 31, 30, 31 ];

const pix2y = [
  0,  0,  1,  1,  0,  0,  1,  1,  2,  2,  3,  3,  2,  2,  3,  3,
  0,  0,  1,  1,  0,  0,  1,  1,  2,  2,  3,  3,  2,  2,  3,  3,
  4,  4,  5,  5,  4,  4,  5,  5,  6,  6,  7,  7,  6,  6,  7,  7,
  4,  4,  5,  5,  4,  4,  5,  5,  6,  6,  7,  7,  6,  6,  7,  7,
  0,  0,  1,  1,  0,  0,  1,  1,  2,  2,  3,  3,  2,  2,  3,  3,
  0,  0,  1,  1,  0,  0,  1,  1,  2,  2,  3,  3,  2,  2,  3,  3,
  4,  4,  5,  5,  4,  4,  5,  5,  6,  6,  7,  7,  6,  6,  7,  7,
  4,  4,  5,  5,  4,  4,  5,  5,  6,  6,  7,  7,  6,  6,  7,  7,
  8,  8,  9,  9,  8,  8,  9,  9, 10, 10, 11, 11, 10, 10, 11, 11,
  8,  8,  9,  9,  8,  8,  9,  9, 10, 10, 11, 11, 10, 10, 11, 11,
 12, 12, 13, 13, 12, 12, 13, 13, 14, 14, 15, 15, 14, 14, 15, 15,
 12, 12, 13, 13, 12, 12, 13, 13, 14, 14, 15, 15, 14, 14, 15, 15,
  8,  8,  9,  9,  8,  8,  9,  9, 10, 10, 11, 11, 10, 10, 11, 11,
  8,  8,  9,  9,  8,  8,  9,  9, 10, 10, 11, 11, 10, 10, 11, 11,
 12, 12, 13, 13, 12, 12, 13, 13, 14, 14, 15, 15, 14, 14, 15, 15,
 12, 12, 13, 13, 12, 12, 13, 13, 14, 14, 15, 15, 14, 14, 15, 15,
  0,  0,  1,  1,  0,  0,  1,  1,  2,  2,  3,  3,  2,  2,  3,  3,
  0,  0,  1,  1,  0,  0,  1,  1,  2,  2,  3,  3,  2,  2,  3,  3,
  4,  4,  5,  5,  4,  4,  5,  5,  6,  6,  7,  7,  6,  6,  7,  7,
  4,  4,  5,  5,  4,  4,  5,  5,  6,  6,  7,  7,  6,  6,  7,  7,
  0,  0,  1,  1,  0,  0,  1,  1,  2,  2,  3,  3,  2,  2,  3,  3,
  0,  0,  1,  1,  0,  0,  1,  1,  2,  2,  3,  3,  2,  2,  3,  3,
  4,  4,  5,  5,  4,  4,  5,  5,  6,  6,  7,  7,  6,  6,  7,  7,
  4,  4,  5,  5,  4,  4,  5,  5,  6,  6,  7,  7,  6,  6,  7,  7,
  8,  8,  9,  9,  8,  8,  9,  9, 10, 10, 11, 11, 10, 10, 11, 11,
  8,  8,  9,  9,  8,  8,  9,  9, 10, 10, 11, 11, 10, 10, 11, 11,
 12, 12, 13, 13, 12, 12, 13, 13, 14, 14, 15, 15, 14, 14, 15, 15,
 12, 12, 13, 13, 12, 12, 13, 13, 14, 14, 15, 15, 14, 14, 15, 15,
  8,  8,  9,  9,  8,  8,  9,  9, 10, 10, 11, 11, 10, 10, 11, 11,
  8,  8,  9,  9,  8,  8,  9,  9, 10, 10, 11, 11, 10, 10, 11, 11,
 12, 12, 13, 13, 12, 12, 13, 13, 14, 14, 15, 15, 14, 14, 15, 15,
 12, 12, 13, 13, 12, 12, 13, 13, 14, 14, 15, 15, 14, 14, 15, 15,
 16, 16, 17, 17, 16, 16, 17, 17, 18, 18, 19, 19, 18, 18, 19, 19,
 16, 16, 17, 17, 16, 16, 17, 17, 18, 18, 19, 19, 18, 18, 19, 19,
 20, 20, 21, 21, 20, 20, 21, 21, 22, 22, 23, 23, 22, 22, 23, 23,
 20, 20, 21, 21, 20, 20, 21, 21, 22, 22, 23, 23, 22, 22, 23, 23,
 16, 16, 17, 17, 16, 16, 17, 17, 18, 18, 19, 19, 18, 18, 19, 19,
 16, 16, 17, 17, 16, 16, 17, 17, 18, 18, 19, 19, 18, 18, 19, 19,
 20, 20, 21, 21, 20, 20, 21, 21, 22, 22, 23, 23, 22, 22, 23, 23,
 20, 20, 21, 21, 20, 20, 21, 21, 22, 22, 23, 23, 22, 22, 23, 23,
 24, 24, 25, 25, 24, 24, 25, 25, 26, 26, 27, 27, 26, 26, 27, 27,
 24, 24, 25, 25, 24, 24, 25, 25, 26, 26, 27, 27, 26, 26, 27, 27,
 28, 28, 29, 29, 28, 28, 29, 29, 30, 30, 31, 31, 30, 30, 31, 31,
 28, 28, 29, 29, 28, 28, 29, 29, 30, 30, 31, 31, 30, 30, 31, 31,
 24, 24, 25, 25, 24, 24, 25, 25, 26, 26, 27, 27, 26, 26, 27, 27,
 24, 24, 25, 25, 24, 24, 25, 25, 26, 26, 27, 27, 26, 26, 27, 27,
 28, 28, 29, 29, 28, 28, 29, 29, 30, 30, 31, 31, 30, 30, 31, 31,
 28, 28, 29, 29, 28, 28, 29, 29, 30, 30, 31, 31, 30, 30, 31, 31,
 16, 16, 17, 17, 16, 16, 17, 17, 18, 18, 19, 19, 18, 18, 19, 19,
 16, 16, 17, 17, 16, 16, 17, 17, 18, 18, 19, 19, 18, 18, 19, 19,
 20, 20, 21, 21, 20, 20, 21, 21, 22, 22, 23, 23, 22, 22, 23, 23,
 20, 20, 21, 21, 20, 20, 21, 21, 22, 22, 23, 23, 22, 22, 23, 23,
 16, 16, 17, 17, 16, 16, 17, 17, 18, 18, 19, 19, 18, 18, 19, 19,
 16, 16, 17, 17, 16, 16, 17, 17, 18, 18, 19, 19, 18, 18, 19, 19,
 20, 20, 21, 21, 20, 20, 21, 21, 22, 22, 23, 23, 22, 22, 23, 23,
 20, 20, 21, 21, 20, 20, 21, 21, 22, 22, 23, 23, 22, 22, 23, 23,
 24, 24, 25, 25, 24, 24, 25, 25, 26, 26, 27, 27, 26, 26, 27, 27,
 24, 24, 25, 25, 24, 24, 25, 25, 26, 26, 27, 27, 26, 26, 27, 27,
 28, 28, 29, 29, 28, 28, 29, 29, 30, 30, 31, 31, 30, 30, 31, 31,
 28, 28, 29, 29, 28, 28, 29, 29, 30, 30, 31, 31, 30, 30, 31, 31,
 24, 24, 25, 25, 24, 24, 25, 25, 26, 26, 27, 27, 26, 26, 27, 27,
 24, 24, 25, 25, 24, 24, 25, 25, 26, 26, 27, 27, 26, 26, 27, 27,
 28, 28, 29, 29, 28, 28, 29, 29, 30, 30, 31, 31, 30, 30, 31, 31,
28, 28, 29, 29, 28, 28, 29, 29, 30, 30, 31, 31, 30, 30, 31, 31 ];
# All the following translated from the not-so-great original IDL Healpix code

function ang2pix_nest(nside::Int64, theta::Array{Float64,1}, phi::Array{Float64,1})
  ns_max = nside #8192 # replace by nside ?
  npix = nside2npix(nside)
  np = length(theta)
  ix = Array{Int64}(undef,np)
  iy = Array{Int64}(undef,np)
  face_num = Array{Int8}(undef,np)
  z = cos.(theta)
  z_abs = abs.(z)
  scaled_phi = mod2pi.(phi) / (π / 2)

  pix_eqt = findall(z_abs .<= 2/3) # equatorial strip
  n_eqt = length(pix_eqt)
  if (n_eqt > 0)
#     (the index of edge lines increase when the longitude=phi goes up)
     jp = floor.(Integer, ns_max * (0.5 .+ scaled_phi[pix_eqt] - z[pix_eqt] * 0.75))
     jm = floor.(Integer, ns_max * (0.5 .+ scaled_phi[pix_eqt] + z[pix_eqt] * 0.75))

     #     finds the face
     face_n = Array{Int8}(undef,n_eqt)
     ifp = div.(jp, ns_max)   # in {0,4}
     ifm = div.(jm, ns_max)
     p_np = findall(ifp .== ifm);
     n_np = length(p_np)
     p_eq = findall(ifp .< ifm)
     n_eq = length(p_eq)
     p_sp = findall(ifp .> ifm)
     n_sp = length(p_sp)
     if (n_np > 0)
          face_n[p_np] = ifp[p_np] .%4 .+ 4
      end
      if (n_eq > 0)
           face_n[p_eq] = ifp[p_eq] .% 4
      end
      if (n_sp > 0)
          face_n[p_sp] = ifm[p_sp] .% 4 .+ 8
      end
      face_num[pix_eqt] = face_n
      ix[pix_eqt] = mod.(jm, ns_max)
      iy[pix_eqt] = ns_max - 1 .- mod.(jp , ns_max)
  end

  pix_pol = findall(z_abs .> 2/3)
  n_pol = length(pix_pol)
  if (n_pol > 0)
      tmp = sqrt(6.0) * sin.(0.5*min.(theta[pix_pol], π .- theta[pix_pol]))
      ntt = floor.(Integer, scaled_phi[pix_pol])
      ntt[findall(ntt.>3)] .= 3
      tp = scaled_phi[pix_pol] - ntt
      jp = floor.(Integer, ns_max * tp .* tmp)
      jm = floor.(Integer, ns_max * (1 .- tp) .* tmp)
      jp = min.(jp, ns_max - 1)
      jm = min.(jm, ns_max - 1)

#     finds the face and pixel's (x,y)
      p_sp = findall(z[pix_pol] .< 0.)
      n_sp = length(p_sp)
      p_np = findall(z[pix_pol] .>= 0.)
      n_np = length(p_np)
      if (n_np > 0)
          face_num[pix_pol[p_np]] = ntt[p_np]
          ix[pix_pol[p_np]] = ns_max - 1 .- jm[p_np]
          iy[pix_pol[p_np]] = ns_max - 1 .- jp[p_np]
      end
      if (n_sp > 0)
          face_num[pix_pol[p_sp]] = ntt[p_sp] .+ 8
          ix[pix_pol[p_sp]] = jp[p_sp]
          iy[pix_pol[p_sp]] = jm[p_sp]
      end
  end
  ix_hi = div.(ix, 128)
  ix_low = mod.(ix, 128)
  iy_hi = div.(iy, 128)
  iy_low = mod.(iy, 128)
  ipf = ((x2pix[ix_hi .+ 1] + y2pix[iy_hi .+ 1]) * (128^2) + (x2pix[ix_low .+ 1] + y2pix[iy_low .+ 1]))
  ipf = floor.(Integer, ipf / ((ns_max / nside) ^ 2))
  return ipf + face_num * nside * nside .+ 1
end

function pix2vec_nest(nside::Int64, ipix::Array{Int64,1})
#   coordinate of the lowest corner of each face
  jrll = [2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4] # in unit of nside
  jpll = [1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7] # in unit of nside/2

  npface = nside * nside
  np = length(ipix)
  fact1 = 1.0/(3.0*nside*nside)
  fact2 = 2.0/(3.0*nside)
 # nr     = Array{Int64}(undef, np)
 # kshift = Array{Int8}(undef,np)
  nr     = Array{Int64}(undef, np)
  kshift = Array{Int8}(undef, np)
  face_num = div.(ipix .- 1, nside*nside).+1
  ipf = mod.(ipix .- 1, nside*nside)
  ip_trunc = div.(ipf, 1024)
  ip_low = mod.(ipf, 1024)
  ip_hi = div.(ip_trunc, 1024)
  ip_med = mod.(ip_trunc, 1024)
  ix = 1024 * pix2x[ip_hi.+1] + 32 * pix2x[ip_med.+1] + pix2x[ip_low.+1]
  iy = 1024 * pix2y[ip_hi.+1] + 32 * pix2y[ip_med.+1] + pix2y[ip_low.+1]
  # Transforms this in (horizontal, vertical) coordinates
  jrt = ix + iy # 'vertical' in {0,2*(nside-1)}
  jpt = ix - iy # 'horizontal' in {-nside+1,nside-1}
  jr = jrll[face_num] * nside - jrt .- 1 #ring number in {1,4*nside-1}

  pix_eqt = findall( (jr .>= nside) .& (jr .<= 3*nside) )
  n_eqt  = length(pix_eqt)
  nr[pix_eqt]     .= nside
  kshift[pix_eqt]  = mod.(jr[pix_eqt] .- nside, 2)

  pix_npl = findall( jr .< nside )
  n_npl = length(pix_npl)
  nr[pix_npl]      = jr[pix_npl]
  kshift[pix_npl] .= 0

  pix_spl = findall( jr .> 3*nside )
  n_spl = length(pix_spl)
  nr[pix_spl]      = 4*nside .- jr[pix_spl]
  kshift[pix_spl] .= 0

#     computes the phi coordinate on the sphere, in [0,2Pi]
  jp = div.(jpll[face_num].*nr .+ jpt .+ 1 + kshift, 2) # 'phi' number in the ring in {1,4*nr}
  #jpt = 0
  #face_num = 0         # free memory
  jp = jp - 4*nside * (jp .> 4*nside)
  jp = jp + 4*nside * (jp .< 1)

# vertex
  iphi_mod = mod.(jp.-1, nr) # in {0,1,... nr-1}
  iphi_rat = div.(jp.-1, nr)      # in {0,1,2,3}
  phi_up = pi / 2.0 * (iphi_rat +  iphi_mod   ./(max.(nr.-1,1.0)))
  phi_dn = pi / 2.0 * (iphi_rat + (iphi_mod.+1)./(nr.+1))

  vertex = Array{Float64}(undef, np, 3, 4);
  vec_out = Array{Float64}(undef, np, 3);

  if (n_eqt > 0)
      phi = (jp[pix_eqt] - (kshift[pix_eqt].+1)*0.50) .* (pi / 2.0 ./ nr[pix_eqt])
      z           = (2*nside .- jr[pix_eqt])*fact2
      sz          = sqrt.(1.0 .- z.^2)
      vec_out[pix_eqt,1] = sz.*cos.(phi)
      vec_out[pix_eqt,2] = sz.*sin.(phi)
      vec_out[pix_eqt,3] = z
      z_nv = (2*nside+1 .-jr[pix_eqt])*fact2
      z_sv = (2*nside-1 .-jr[pix_eqt])*fact2
      k1 = findall(jr[pix_eqt] .==   nside) # northern transition
      nk1 = length(k1)
      k3 = findall(jr[pix_eqt] .== 3*nside) # southern transition
      nk3 = length(k3)
      hdelta_phi = pi./(4.00*nr[pix_eqt])
      vertex[pix_eqt,:,:] = sub_compute_vertices(z, z_nv, z_sv, phi, phi, phi, hdelta_phi)
      if (nk1 > 0)
          z_nv[k1] .=  1.00 - (nside-1.0)^2 * fact1
          # phi_nv = phi_up
          vertex[pix_eqt[k1],:,:] = sub_compute_vertices(z[k1],z_nv[k1],z_sv[k1],phi[k1],phi_up[pix_eqt[k1]],phi[k1], hdelta_phi[k1])
      end
      if (nk3 > 0)
          z_sv[k3] .= -1.00 + (nside-1.0)^2 * fact1
          # phi_sv = phi_up
          vertex[pix_eqt[k3],:,:] = sub_compute_vertices(z[k3],z_nv[k3],z_sv[k3],phi[k3],phi[k3],phi_up[pix_eqt[k3]], hdelta_phi[k3])
      end
  end

  if (n_npl > 0)
      phi = (jp[pix_npl] - (kshift[pix_npl].+1)*0.50) .* (pi / 2.0 ./ nr[pix_npl])
      z           = 1.0 .- (nr[pix_npl]).^2 * fact1
      sz          = nr[pix_npl] .* sqrt.( fact1 * (1.0.+z) )
      vec_out[pix_npl,1] = sz.*cos.(phi)
      vec_out[pix_npl,2] = sz.*sin.(phi)
      vec_out[pix_npl,3] = z
      z_nv = 1.00 .- (nr[pix_npl].-1.0).^2*fact1
      z_sv = 1.00 .- (nr[pix_npl].+1.0).^2*fact1
      vertex[pix_npl,:,:] = sub_compute_vertices(z, z_nv, z_sv, phi, phi_up[pix_npl], phi_dn[pix_npl], pi./(4.00*nr[pix_npl]))
  end

  if (n_spl > 0)
      phi = (jp[pix_spl] - (kshift[pix_spl].+1)*0.50) .* (pi / 2.0 ./ nr[pix_spl])
      z           = -1.0 .+ (nr[pix_spl]).^2 * fact1
      # sz          = SQRT(1.0 - z^2)
      sz          = nr[pix_spl] .* sqrt.( fact1 * (1.0.-z) )
      vec_out[pix_spl,1] = sz.*cos.(phi)
      vec_out[pix_spl,2] = sz.*sin.(phi)
      vec_out[pix_spl,3] = z
      z_nv = - 1.00 .+ (nr[pix_spl].+1.0).^2*fact1
      z_sv = - 1.00 .+ (nr[pix_spl].-1.0).^2*fact1
      vertex[pix_spl,:,:] = sub_compute_vertices(z, z_nv, z_sv, phi, phi_dn[pix_spl], phi_up[pix_spl], pi./(4.00*nr[pix_spl]))
   end

  return vec_out, vertex
end # pix2vec_nest


function pix2ang_nest(nside::Int64, ipix::Array{Int64,1})
    jrll = [ 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 ]
    jpll = [ 1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7 ]

    fact2 = 2.0 / (3.0 * nside)
    sfact = 1.0 / (sqrt(6.0) * nside)
    # face number in {1,12} and pixel number within the face
    face_num = div.(ipix .- 1, nside*nside).+1
    ipf = mod.(ipix .- 1, nside*nside)
    ip_trunc = div.(ipf, 1024)
    ip_low = mod.(ipf, 1024)
    ip_hi = div.(ip_trunc, 1024)
    ip_med = mod.(ip_trunc, 1024)
    ix = 1024 * pix2x[ip_hi.+1] + 32 * pix2x[ip_med.+1] + pix2x[ip_low.+1]
    iy = 1024 * pix2y[ip_hi.+1] + 32 * pix2y[ip_med.+1] + pix2y[ip_low.+1]
    # Transforms this in (horizontal, vertical) coordinates
    jrt = ix + iy # 'vertical' in {0,2*(nside-1)}
    jpt = ix - iy # 'horizontal' in {-nside+1,nside-1}
    jr = jrll[face_num] * nside - jrt .- 1

#    nr     = Array{Int64}(undef, np)
#    theta  = Array{Float64}(undef, np)
    np = length(ipix)
    nr     = Array{Int64}(undef, np)
    theta  = Array{Float64}(undef, np)

    # equatorial region: [nside,3*nside]
    pix_eqt = findall( (jr .>= nside) .& (jr .< 3*nside) )
    nr[pix_eqt]     .= nside   # equatorial region
    theta[pix_eqt]  = acos.(  (2*nside .- jr[pix_eqt])*fact2  )

    # north pole: [1,nside]
    pix_npl = findall( jr .< nside )
    nr[pix_npl]     = jr[pix_npl]
    theta[pix_npl] = 2.0 * asin.( nr[pix_npl] * sfact )

    # south pole: [3*nside,4*nside-1]
    pix_spl = findall( jr .>= 3*nside )
    nr[pix_spl]     = 4*nside .- jr[pix_spl]
    theta[pix_spl] = pi .- 2.0 * asin.( nr[pix_spl] * sfact )

 #  computes the phi coordinate on the sphere, in [0,2Pi]
    jp = jpll[face_num].*nr + jpt
    jp += 8*nside .* (jp .< 0)    # in [0,8Nside-1] because of 1/2 step stagger
    phi = (pi/4) * (jp ./ nr)
  return theta, phi
end

function pix2xy_nest(nside, ipf)
ip_trunc = div.(ipf, 1024)
ip_low = mod.(ipf, 1024)
ip_hi = div.(ip_trunc, 1024)
ip_med = mod.(ip_trunc, 1024)
ix = 1024 * pix2x[ip_hi.+1] + 32 * pix2x[ip_med.+1] + pix2x[ip_low.+1]
iy = 1024 * pix2y[ip_hi.+1] + 32 * pix2y[ip_med.+1] + pix2y[ip_low.+1]
return ix, iy
end


function xy2pix_nest(nside, ix, iy, face_num::Int64) # ix and iy may be vectors
ix_low = ix .& 127     # last 7 bits
iy_low = iy .& 127     # last 7 bits
ix_hi = div.(ix, 128) ;              # truncate out last 7 bits
iy_hi = div.(iy, 128);
ipnest = (x2pix[ix_low.+1] + y2pix[iy_low.+1]) + (x2pix[ix_hi.+1] + y2pix[iy_hi.+1]) * 16384 .+ face_num *nside^2
return ipnest
end

function neighbours_nest(nside::Int64, ipix::Int64)
  #       Find nearest neighbors pixels in Nested scheme
  #       This is a simplified version -- nside > 1 and nside < 2^13
  #       The neighbours are ordered in the following way:
  #       First pixel is the one to the south (the one west of the south
  #        direction is taken for the pixels which don't have a southern neighbour). From
  #       then on the neighbours are ordered in the clockwise direction
  #       (from outside the sphere) about the pixel with number ipix.
  #       number of neighbours: usually 8, sometimes 7 (for 8 pixels) or 6 (for Nside=1)
  npix = nside2npix(nside);
  ipix0 = ipix-1; # we remove one for julia <-> IDL
  nsidesq = div(npix, 12);
  face_num= div(ipix0,nsidesq);
  local_magic1=div((nsidesq-1),3);
  local_magic2=2*local_magic1;
  ipf = mod(ipix0,nsidesq);
  ix, iy = pix2xy_nest(nside,ipf);
  ixm=ix-1;
  ixp=ix+1;
  iym=iy-1;
  iyp=iy+1;
  nneigh=8;                  #Except in special cases below
  icase = 0;                 #general case
  #     Exclude corners
  if (ipf == local_magic2)
    icase = 5 #WestCorner
  elseif (ipf == (nsidesq-1))
    icase = 6 #NorthCorner
  elseif (ipf == 0)
    icase = 7 #SouthCorner
  elseif (ipf == local_magic1)
    icase = 8 #EastCorner
  end


  if (icase == 0)
    #     Detect edges
    if ((ipf & local_magic1) == local_magic1)
      icase = 1 #NorthEast
    end
    if ((ipf & local_magic1) == 0)
      icase = 2 #SouthWest
    end
    if ((ipf & local_magic2) == local_magic2)
      icase = 3 #NorthWest
    end
    if ((ipf & local_magic2) == 0)
      icase = 4 #SouthEast
    end
    if icase == 0
      # inside a face
      vx = [ixm, ixm, ixm, ix , ixp, ixp, ixp, ix];
      vy = [iym, iy , iyp, iyp, iyp, iy , iym, iym];
      list = xy2pix_nest(nside, vx, vy, face_num).+1;
      return list
    end
  end

  ia=  div(face_num, 4)                  # in {0,2}
  ib=  mod(face_num, 4)          # in {0,3}
  ibp= mod((ib+1), 4)
  ibm= mod((ib+4-1), 4)
  ib2= mod((ib+2), 4)

  if ia == 0          #North Pole region
    if icase ==1   #NorthEast edge
      other_face=0+ibp
      n_8 = xy2pix_nest(nside, ix , iym, face_num);
      n_1 = xy2pix_nest(nside, ixm, iym, face_num);
      n_2 = xy2pix_nest(nside, ixm, iy , face_num);
      n_3 = xy2pix_nest(nside, ixm, iyp, face_num);
      n_4 = xy2pix_nest(nside, ix , iyp, face_num);
      ipo = mod(swapLSBMSB(ipf), nsidesq); #East-West flip
      ixo,iyo = pix2xy_nest(nside,ipo);
      n_5 = xy2pix_nest( nside, ixo+1 , iyo, other_face);
      n_6 = other_face*nsidesq+ipo;
      n_7 = xy2pix_nest( nside, ixo-1, iyo, other_face);
    elseif icase == 2 #SouthWest edge
      other_face=4+ib
      ipo=mod(invLSB(ipf), nsidesq); #SW-NE flip
      ixo,iyo = pix2xy_nest(nside,ipo);
      n_1 = xy2pix_nest(nside, ixo, iyo-1, other_face);
      n_2 = other_face*nsidesq+ipo;
      n_3 = xy2pix_nest(nside, ixo, iyo+1, other_face);
      n_8 = xy2pix_nest(nside, ix , iym, face_num);
      n_4 = xy2pix_nest(nside, ix , iyp, face_num);
      n_7 = xy2pix_nest(nside, ixp, iym, face_num);
      n_6 = xy2pix_nest(nside, ixp, iy , face_num);
      n_5 = xy2pix_nest(nside, ixp, iyp, face_num);
    elseif icase==3 #NorthWest edge
      other_face=0+ibm
      ipo=mod(swapLSBMSB(ipf), nsidesq); #East-West flip
      ixo,iyo = pix2xy_nest(nside,ipo);
      n_3 = xy2pix_nest( nside, ixo, iyo-1, other_face);
      n_4 = other_face*nsidesq+ipo;
      n_5 = xy2pix_nest(nside, ixo, iyo+1, other_face);
      n_1 = xy2pix_nest(nside, ixm, iym, face_num);
      n_2 = xy2pix_nest(nside, ixm, iy , face_num);
      n_8 = xy2pix_nest(nside, ix , iym, face_num);
      n_7 = xy2pix_nest(nside, ixp, iym, face_num);
      n_6 = xy2pix_nest(nside, ixp, iy , face_num);
    elseif icase == 4 #SouthEast edge
      other_face=4+ibp
      n_2 = xy2pix_nest(nside, ixm, iy , face_num);
      n_3 = xy2pix_nest(nside, ixm, iyp, face_num);
      n_4 = xy2pix_nest(nside, ix , iyp, face_num);
      n_5 = xy2pix_nest(nside, ixp, iyp, face_num);
      n_6 = xy2pix_nest(nside, ixp, iy , face_num);
      ipo=mod(invMSB(ipf), nsidesq) #SE-NW flip
      ixo,iyo = pix2xy_nest(nside,ipo);
      n_7 = xy2pix_nest(nside, ixo+1, iyo, other_face);
      n_8 = other_face*nsidesq+ipo;
      n_1 = xy2pix_nest(nside, ixo-1, iyo, other_face);
    elseif icase==5 #West corner
      nneigh=7
      other_face=4+ib
      n_2=other_face*nsidesq+nsidesq-1
      n_1=n_2-2
      other_face=0+ibm
      n_3=other_face*nsidesq+local_magic1
      n_4=n_3+2
      n_5=ipix0+1
      n_6=ipix0-1
      n_7=ipix0-2
    elseif icase==6 #North corner
      n_1=ipix0-3
      n_2=ipix0-1
      n_8=ipix0-2
      other_face=0+ibm
      n_4=other_face*nsidesq+nsidesq-1
      n_3=n_4-2
      other_face=0+ib2
      n_5=other_face*nsidesq+nsidesq-1
      other_face=0+ibp
      n_6=other_face*nsidesq+nsidesq-1
      n_7=n_6-1
    elseif icase==7 #South corner
      other_face=8+ib
      n_1=other_face*nsidesq+nsidesq-1
      other_face=4+ib
      n_2=other_face*nsidesq+local_magic1
      n_3=n_2+2
      n_4=ipix0+2
      n_5=ipix0+3
      n_6=ipix0+1
      other_face=4+ibp
      n_8=other_face*nsidesq+local_magic2
      n_7=n_8+1
    elseif icase==8 #East corner
      nneigh=7
      n_2=ipix0-1
      n_3=ipix0+1
      n_4=ipix0+2
      other_face=0+ibp
      n_6=other_face*nsidesq+local_magic2
      n_5=n_6+1
      other_face=4+ibp
      n_7=other_face*nsidesq+nsidesq-1
      n_1=n_7-1
    end
    # north
  end

  if ia == 1 #Equatorial region
    if icase ==1 #NorthEast edge
      other_face=0+ib
      n_8 = xy2pix_nest(nside, ix , iym, face_num);
      n_1 = xy2pix_nest(nside, ixm, iym, face_num);
      n_2 = xy2pix_nest(nside, ixm, iy , face_num);
      n_3 = xy2pix_nest(nside, ixm, iyp, face_num);
      n_4 = xy2pix_nest(nside, ix , iyp, face_num);
      ipo = mod(invLSB(ipf), nsidesq) #NE-SW flip
      ixo,iyo = pix2xy_nest(nside,ipo);
      n_5 = xy2pix_nest(nside, ixo , iyo+1, other_face);
      n_6 = other_face*nsidesq+ipo;
      n_7 = xy2pix_nest(nside, ixo, iyo-1, other_face);
    elseif icase == 2 #SouthWest edge
      other_face=8+ibm;
      ipo = mod(invLSB(ipf), nsidesq); #SW-NE flip
      ixo,iyo = pix2xy_nest(nside,ipo);
      n_1 = xy2pix_nest(nside, ixo, iyo-1, other_face);
      n_2 = other_face*nsidesq+ipo;
      n_3 = xy2pix_nest(nside, ixo, iyo+1, other_face);
      n_8 = xy2pix_nest(nside, ix , iym, face_num);
      n_4 = xy2pix_nest(nside, ix , iyp, face_num);
      n_7 = xy2pix_nest(nside, ixp, iym, face_num);
      n_6 = xy2pix_nest(nside, ixp, iy , face_num);
      n_5 = xy2pix_nest(nside, ixp, iyp, face_num);
    elseif icase== 3 #NorthWest edge
      other_face=0+ibm;
      ipo=mod(invMSB(ipf), nsidesq); #NW-SE flip
      ixo,iyo = pix2xy_nest(nside,ipo);
      n_3 = xy2pix_nest(nside, ixo-1, iyo, other_face);
      n_4 = other_face*nsidesq+ipo;
      n_5 = xy2pix_nest(nside, ixo+1, iyo, other_face);
      n_1 = xy2pix_nest(nside, ixm, iym, face_num);
      n_2 = xy2pix_nest(nside, ixm, iy , face_num);
      n_8 = xy2pix_nest(nside, ix , iym, face_num);
      n_7 = xy2pix_nest(nside, ixp, iym, face_num);
      n_6 = xy2pix_nest(nside, ixp, iy , face_num);
    elseif icase ==4 #SouthEast edge
      other_face=8+ib;
      n_2 = xy2pix_nest(nside, ixm, iy , face_num);
      n_3 = xy2pix_nest(nside, ixm, iyp, face_num);
      n_4 = xy2pix_nest(nside, ix , iyp, face_num);
      n_5 = xy2pix_nest(nside, ixp, iyp, face_num);
      n_6 = xy2pix_nest(nside, ixp, iy , face_num);
      ipo = mod(invMSB(ipf), nsidesq);#SE-NW flip
      ixo,iyo = pix2xy_nest(nside,ipo);
      n_7 = xy2pix_nest(nside, ixo+1, iyo, other_face);
      n_8 = other_face*nsidesq+ipo;
      n_1 = xy2pix_nest(nside, ixo-1, iyo, other_face);
    elseif icase ==5 #West corner
      other_face=8+ibm;
      n_2=other_face*nsidesq+nsidesq-1;
      n_1=n_2-2;
      other_face=4+ibm;
      n_3=other_face*nsidesq+local_magic1;
      other_face=0+ibm;
      n_4=other_face*nsidesq;
      n_5=n_4+1;
      n_6=ipix0+1;
      n_7=ipix0-1;
      n_8=ipix0-2;
    elseif icase ==6 #North corner
      nneigh=7;
      n_1=ipix0-3;
      n_2=ipix0-1;
      other_face=0+ibm;
      n_4=other_face*nsidesq+local_magic1;
      n_3=n_4-1;
      other_face=0+ib;
      n_5=other_face*nsidesq+local_magic2;
      n_6=n_5-2;
      n_7=ipix0-2;
    elseif  icase == 7 #South corner
      nneigh=7;
      other_face=8+ibm;
      n_1=other_face*nsidesq+local_magic1;
      n_2=n_1+2;
      n_3=ipix0+2;
      n_4=ipix0+3;
      n_5=ipix0+1;
      other_face=8+ib;
      n_7=other_face*nsidesq+local_magic2;
      n_6=n_7+1;
    elseif icase==8 #East corner
      other_face=8+ib;
      n_8=other_face*nsidesq+nsidesq-1;
      n_1=n_8-1;
      n_2=ipix0-1;
      n_3=ipix0+1;
      n_4=ipix0+2;
      other_face=0+ib;
      n_6=other_face*nsidesq;
      n_5=n_6+2;
      other_face=4+ibp;
      n_7=other_face*nsidesq+local_magic2;
    end
    # equator
  end

  if ia == 2
    if icase == 1 #NorthEast edge
      other_face=4+ibp
      n_8 = xy2pix_nest(nside, ix , iym, face_num);
      n_1 = xy2pix_nest(nside, ixm, iym, face_num);
      n_2 = xy2pix_nest(nside, ixm, iy , face_num);
      n_3 = xy2pix_nest(nside, ixm, iyp, face_num);
      n_4 = xy2pix_nest(nside, ix , iyp, face_num);
      ipo = mod(invLSB(ipf), nsidesq); #NE-SW flip
      ixo,iyo =  pix2xy_nest(nside,ipo);
      n_5 = xy2pix_nest(nside, ixo , iyo+1, other_face);
      n_6 = other_face*nsidesq+ipo;
      n_7 = xy2pix_nest(nside, ixo, iyo-1, other_face);
    elseif icase == 2 #SouthWest edge
      other_face=8+ibm
      ipo=mod(swapLSBMSB(ipf), nsidesq); #W-E flip
      ixo,iyo = pix2xy_nest(nside,ipo);
      n_1 = xy2pix_nest(nside, ixo-1, iyo, other_face);
      n_2 = other_face*nsidesq+ipo
      n_3 = xy2pix_nest(nside, ixo+1, iyo, other_face);
      n_8 = xy2pix_nest(nside, ix , iym, face_num);
      n_4 = xy2pix_nest(nside, ix , iyp, face_num);
      n_7 = xy2pix_nest(nside, ixp, iym, face_num);
      n_6 = xy2pix_nest(nside, ixp, iy , face_num);
      n_5 = xy2pix_nest(nside, ixp, iyp, face_num);
    elseif icase == 3 #NorthWest edge
      other_face=4+ib;
      ipo=mod(invMSB(ipf) , nsidesq); #NW-SE flip
      ixo,iyo = pix2xy_nest(nside,ipo);
      n_3 = xy2pix_nest(nside, ixo-1, iyo, other_face);
      n_4 = other_face*nsidesq+ipo;
      n_5 = xy2pix_nest(nside, ixo+1, iyo, other_face);
      n_1 = xy2pix_nest(nside, ixm, iym, face_num);
      n_2 = xy2pix_nest(nside, ixm, iy , face_num);
      n_8 = xy2pix_nest(nside, ix , iym, face_num);
      n_7 = xy2pix_nest(nside, ixp, iym, face_num);
      n_6 = xy2pix_nest(nside, ixp, iy , face_num);
    elseif icase == 4 #SouthEast edge
      other_face=8+ibp
      n_2 = xy2pix_nest(nside, ixm, iy , face_num);
      n_3 = xy2pix_nest(nside, ixm, iyp, face_num);
      n_4 = xy2pix_nest(nside, ix , iyp, face_num);
      n_5 = xy2pix_nest(nside, ixp, iyp, face_num);
      n_6 = xy2pix_nest(nside, ixp, iy , face_num);
      ipo = mod(swapLSBMSB(ipf), nsidesq); #E-W flip
      ixo,iyo = pix2xy_nest(nside,ipo);
      n_7 = xy2pix_nest(nside, ixo, iyo+1, other_face);
      n_8 = other_face*nsidesq+ipo
      n_1 = xy2pix_nest(nside, ixo, iyo-1, other_face);
    elseif icase == 5 #West corner
      nneigh=7;
      other_face=8+ibm;
      n_2=other_face*nsidesq+local_magic1;
      n_1=n_2-1;
      other_face=4+ib;
      n_3=other_face*nsidesq;
      n_4=n_3+1;
      n_5=ipix0+1;
      n_6=ipix0-1;
      n_7=ipix0-2;
    elseif icase == 6 #North corner
      n_1=ipix0-3;
      n_2=ipix0-1;
      other_face=4+ib;
      n_4=other_face*nsidesq+local_magic1;
      n_3=n_4-1;
      other_face=0+ib;
      n_5=other_face*nsidesq;
      other_face=4+ibp;
      n_6=other_face*nsidesq+local_magic2;
      n_7=n_6-2;
      n_8=ipix0-2;
    elseif icase == 7 #South corner
      other_face=8+ib2;
      n_1=other_face*nsidesq;
      other_face=8+ibm;
      n_2=other_face*nsidesq;
      n_3=n_2+1;
      n_4=ipix0+2;
      n_5=ipix0+3;
      n_6=ipix0+1;
      other_face=8+ibp;
      n_8=other_face*nsidesq;
      n_7=n_8+2;
    elseif  icase==8 #East corner
      nneigh=7;
      other_face=8+ibp;
      n_7=other_face*nsidesq+local_magic2;
      n_1=n_7-2;
      n_2=ipix0-1;
      n_3=ipix0+1;
      n_4=ipix0+2;
      other_face=4+ibp;
      n_6=other_face*nsidesq;
      n_5=n_6+2;
    end
  end                           # south

  if nneigh == 8
    return  [n_1, n_2, n_3, n_4, n_5, n_6, n_7, n_8].+1
  else
    return [n_1, n_2, n_3, n_4, n_5, n_6, n_7].+1
  end

end

function healpix_round_star(n::Int64;radius=1.0) # round, radius = 1
# instantiate the base geometry
nside = 2^n # note: this is correct, n is the map resolution parameter
npix = nside2npix(nside)
vertices_xyz = Array{Float64}(undef, npix, 3, 5);
vertices_spherical = Array{Float64}(undef , npix, 3, 5);

# Explanation:
# * we need to use angles to compute the vertices of centers and corners
# * right now we get the xyz of centers and corners, convert to angles, then scale by radius and convert back to xyz

#cartesian coordinates for center & corners
vertices_xyz[:,:,5],vertices_xyz[:,:,1:4]=pix2vec_nest(nside,collect(1:npix));
vertices_spherical[:,1,:] .= radius;
# setup spherical (theta,phi) for center
vertices_spherical[:,2,5],vertices_spherical[:,3,5] = pix2ang_nest(nside, collect(1:npix));
# and now for corners
for j = 1:4
	vertices_spherical[:,2,j], vertices_spherical[:,3,j] = vec2ang(vertices_xyz[:,:,j]);
end

vertices_xyz[:,1,:] = radius*sin.(vertices_spherical[:,2,:]).*cos.(vertices_spherical[:,3,:]); # X
vertices_xyz[:,2,:] = radius*sin.(vertices_spherical[:,2,:]).*sin.(vertices_spherical[:,3,:]); # Y
vertices_xyz[:,3,:] = radius*cos.(vertices_spherical[:,2,:]); # Z

star_base_geom = base_geometry(npix, vertices_xyz, vertices_spherical);
end



function healpix_ellipsoid_star(n, a, b, c) # ellipsoid
# instantiate the base geometry
nside = 2^n
npix = nside2npix(nside)
vertices_xyz = Array{Float64}(undef, npix, 3, 5);
vertices_spherical = Array{Float64}(undef , npix, 3, 5);
#cartesian coordinates for center & corners
vertices_xyz[:,:,5],vertices_xyz[:,:,1:4]=pix2vec_nest(nside,collect(1:npix));
vertices_spherical[:,1,:] .= radius;
# setup spherical (theta,phi) for center
vertices_spherical[:,2,5],vertices_spherical[:,3,5] = pix2ang_nest(nside, collect(1:npix));
# and now for corners
for j = 1:4
	vertices_spherical[:,2,j], vertices_spherical[:,3,j] = vec2ang(vertices_xyz[:,:,j]);
end

# but now we change the radius
vertices_xyz[:,1,:] = a*sin.(vertices_spherical[:,2,:]).*cos.(vertices_spherical[:,3,:]); # X
vertices_xyz[:,2,:] = b*sin.(vertices_spherical[:,2,:]).*sin.(vertices_spherical[:,3,:]); # Y
vertices_xyz[:,3,:] = c*cos.(vertices_spherical[:,2,:]); # Z

vertices_spherical[:,1,:] = sqrt.(vertices_xyz[:,1,:].^2 + vertices_xyz[:,2,:].^2 + vertices_xyz[:,3,:].^2); # radius, TBD replace by norm()

star_base_geom = base_geometry(npix, vertices_xyz, vertices_spherical);
end

# Figure out where pole radius is in terms of pixels
function healpix_rapidrot_star(n, stellar_parameters) # ellipsoid
nside = 2^n
npix = nside2npix(nside)
vertices_xyz = Array{Float64}(undef, npix, 3, 5);
vertices_spherical = Array{Float64}(undef , npix, 3, 5);
#cartesian coordinates for center & corners
vertices_xyz[:,:,5],vertices_xyz[:,:,1:4]=pix2vec_nest(nside,collect(1:npix));
vertices_spherical[:,1,:] .= radius;
# setup spherical (theta,phi) for center
vertices_spherical[:,2,5],vertices_spherical[:,3,5] = pix2ang_nest(nside, collect(1:npix));
# and now for corners
for j = 1:4
	vertices_spherical[:,2,j], vertices_spherical[:,3,j] = vec2ang(vertices_xyz[:,:,j]);
end

R_pole = stellar_parameters.radius;
ω = stellar_parameters.frac_escapevel;
vertices_spherical[:,1,:] = 3.0*R_pole.*cos.((pi + acos.(ω*sin.(vertices_spherical[:,2,:])))/3.)./(ω*sin.(vertices_spherical[:,2,:]));

# Rewrite pole radius values
vertices_spherical[:,1,:][find(vertices_spherical[:,1,:] .== Inf)] = R_pole;

# but now we change the radius
for i=1:npix
  vertices_xyz[i,1,:] = vertices_spherical[i,1,:]*sin.(vertices_spherical[i,2,:]).*cos.(vertices_spherical[i,3,:]); # X
  vertices_xyz[i,2,:] = vertices_spherical[i,1,:]*sin.(vertices_spherical[i,2,:]).*sin.(vertices_spherical[i,3,:]); # Y
  vertices_xyz[i,3,:] = vertices_spherical[i,1,:]*cos.(vertices_spherical[i,2,:]); # Z
end

star_base_geom = base_geometry(npix, vertices_xyz, vertices_spherical);
end

function all_neighbours_nest(n)
nside = 2^n;
neighbors = [neighbours_nest(nside, i) for i=1:nside2npix(nside)];
return neighbors
end

function tv_neighbours_healpix(n)
# Complete Neighbor setup (healpix)
neighbors = all_neighbours_nest(n); #neighbors[ipix] will give the list of all neighbors of pixel [ipix]
south_neighbors=[(neighbors[i])[1] for i=1:length(neighbors)]
west_neighbors=[(neighbors[i])[3] for i=1:length(neighbors)]
south_neighbors_reverse=[findall(south_neighbors.==i) for i=1:length(neighbors)]
west_neighbors_reverse=[findall(west_neighbors.==i) for i=1:length(neighbors)]
return neighbors,south_neighbors,west_neighbors,south_neighbors_reverse,west_neighbors_reverse
end
