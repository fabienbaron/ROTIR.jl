using FFTW, LinearAlgebra, Statistics, OptimPackNextGen, JLD2# PyPlot
cd("/home/baron/SOFTWARE/rotir2/demos")
include("../src/ROTIR.jl"); using Main.ROTIR;
include("io-utils.jl");
files = "./IODATA/".*["Kraken_Io_ 2025-02-23 02_44_54.271_I_R_V.fits",
"Kraken_Io_2025-02-23 02_56_05.527_I_R_V.fits",
"Kraken_Io_2025-02-23 03_17_34.8_I_R_V.fits",
"Kraken_Io_2025-02-23 03_24_30.641_I_R_V.fits"]

 tepochs=[sum([02 44 54.271].*[1/24 1/(24*60) 1/(24*3600)]), 
 sum([02 56 05.527].*[1/24 1/(24*60) 1/(24*3600)]),
 sum([03 17 34.8].*[1/24 1/(24*60) 1/(24*3600)]),
 sum([03 24 30.641].*[1/24 1/(24*60) 1/(24*3600)])]  # expressed in days

nepochs = length(files)
N = 256
img = zeros(Float32, N, N, 3, nepochs)
for i=1:nepochs
    for j=1:3
       img[:,:,j,i]=rotl90(readfits(files[i])[128:128+N-1,128:128+N-1,j])
   end
end
# Let's try working on I band only
img = img[:,:,1,:]
mask = sum(img, dims=3); mask = (mask .> maximum(mask)/10)
img = img .* mask
data = [img[:,:,i] for i=1:nepochs]#[vec(rfft(fftshift(img[:,:,i]))) for i=1:nepochs]

nx = 256
pixsize = 0.011
oversample = 1
nx_samples = nx*oversample
pixsize_samples = pixsize/oversample
freqs_u = Float32.((fftfreq(nx_samples))/pixsize_samples*2.0626480624709636e8)
freqs_v = Float32.((rfftfreq(nx_samples))/pixsize_samples*2.0626480624709636e8)
const uv = hcat([ [i,j] for i in freqs_u for j in freqs_v]...)

# Io is an ellipsoid, dimensions~ 3,660.0 × 3,637.4 × 3,630.6 km
star_params = (
              surface_type   = 1      ,   # Round:0, Ellipsoid: 1, Rapid Rotator:2, Roche: 3
              radius_x              = 1.0,   # milliarcseconds (radius at pole)
              radius_y              = 1.0,
              radius_z              = 1.0, 
              ldtype         =      1,   # LD type  1: Linear 2: quadratic 3: power (Hestroffer)
              ld1            =      0, # limb darkening,first coefficient is for LD law type, then LD coefficients
              ld2            =    0.0,   # second ld coeff, used if needed
              inclination    = 86,  # degrees; inclination
              position_angle =      0,  # degrees; position_angle
              rotation_period=    1.769138,  # rotation period in days
              )

# Size of polyft = length(uv)*nside2npix(2^n)*64 bytes = length(uv)*nside2npix(2^n)*8/(1e9) GB
nn =   [1,2,3,4,5,6,7]
tmap = [zeros(Float32, nside2npix(2^n)) for n in nn]
model = zeros(Float32, N,N,length(nn))
ftol= (0,1f-12); xtol=(0,1f-12); gtol=(0,1f-12);
i=5
#while (i<7)
GC.gc()
stars = create_star_multiepochs(tessellation_healpix(nn[i]), star_params, tepochs);
println("Tesselation Healpix level: $(nn[i])")
setup_fullft!(uv, stars)

x_sol = tmap[i][stars[1].index_quads_visible] # temperature map including only the visible pixels
psfs = zeros(Float32, nx, nx, nepochs) # Placeholder for PSF, at the moment just flux
psfs[div(nx,2)+1, div(nx,2)+1,:] .= 1.0

# Object step
#
crit_obj = (x,g)->crit_obj_fg(x,g,stars,psfs,data)
x_sol = OptimPackNextGen.vmlmb(crit_obj, x_sol, verb=true,lower=0f0, maxiter=200, blmvm=false, xtol = xtol, ftol = ftol, gtol=gtol);
#for z=1:3
x_sol = OptimPackNextGen.vmlmb(crit_obj, x_sol+2000*rand(Float32, length(x_sol)), verb=true,lower=0f0, maxiter=200, blmvm=false, xtol = xtol, ftol = ftol, gtol=gtol);
#end

# PSF/Flux step
tmap_visible = copy(x_sol)
crit_psfs = (x,g)->crit_psf_fg(x,g,stars,tmap_visible,data)
psfs = OptimPackNextGen.vmlmb(crit_psfs, psfs, verb=true,lower=0f0, maxiter=200, blmvm=false, xtol = xtol, ftol = ftol, gtol=gtol);


# Plot result
tmap[i][stars[1].index_quads_visible] = x_sol;
plot2d(tmap[i], stars[1])




# clf(); imshow(model[:,:,i])    
# tmap_new, _ = upsample_map_stars(tmap[i], stars, star_params, tepochs)
# if i<length(nn)
#     tmap[i+1][:] = tmap_new
#     end
#    i+=1
#end 
