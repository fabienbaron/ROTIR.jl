include("../src/ROTIR.jl"); using Main.ROTIR
include("../src/oitools-ext.jl")

using BenchmarkTools; 
# LOAD DATA
oifitsfiles=["./data/SU_Lyn.oifits"]
nepochs, tepochs, data = readoifits_multiepochs(oifitsfiles,filter_bad_data=true);
tepochs = tepochs .- tepochs[1]; # First epoch set as t=0
data = dataF32.(data)
tepochs = Float32.(tepochs)

# To use a Healpix scheme
n=3; 
tessels = tessellation_healpix(n)

## Rapid rotator
star_params = (
surface_type   = 1f0      ,   # Round:0, Ellipsoid: 1, Rapid Rotator:2, Roche: 3
radius_x              = 1.75f0,   # milliarcseconds (radius at pole)
radius_y              = 1.75f0,
radius_z              =  1.6f0, 
tpole          =   3500.0f0, #  # Kelvin (at pole)
ldtype         =      3f0,   # LD type  1: Linear 2: quadratic 3: power (Hestroffer)
ld1            =  0.26f0, # limb darkening,first coefficient is for LD law type, then LD coefficients
ld2            =    0.0f0,   # second ld coeff, used if needed
inclination    = 76f0,  # degrees; inclination
position_angle =      49f0,  # degrees; position_angle
rotation_period=    254.8f0,  # rotation period in days
beta  =0.08f0
)

stars = create_star_multiepochs(tessels, star_params, tepochs);
tmap_start = 1000*ones(Float32, stars[1].npix)
setup_oi!(data, stars)

# # SETUP REGULARIZATION
# regularizers_1 = [["tv", 0.1, tv_neighbours_healpix(n),1:length(tmap_start)]];
# tmap_1 =  image_reconstruct_oi(tmap_start, data, stars, maxiter = 500, lower=1000, regularizers = regularizers_1, verbose = true);

# regularizers_2 = [["tv", 0.1, tv_neighbours_healpix_visible(n, stars),1:length(tmap_start)]];
# tmap_2 =  image_reconstruct_oi(tmap_start, data, stars, maxiter = 500, lower=1000, regularizers = regularizers_2, verbose = true);

# plot2d(tmap_1, stars[1])

f = x->chi2s(x, stars[1], data[1], verbose=false)
x = rand(Float32, tessels.npix)
@benchmark f(x)


@views function setup_polyft_single(data, star_epoch_geom; T=Float32)
  # Polyflux is the weight of each pixel, proportional to the surface
  pjx = star_epoch_geom.projx;
  pjy = star_epoch_geom.projy;
  kx =  data.uv[1,:] * T(-pi / (180*3600000));
  ky =  data.uv[2,:] * T( pi / (180*3600000));
  # note: check definition sinc(x) = sin(pi*x)/(pi*x)
  polyft = -im*T(1/(2pi))*(
    (sinc.( (kx*transpose(pjx[:,2]-pjx[:,1])) + (ky*transpose(pjy[:,2]-pjy[:,1])) ).*cis.(-T(pi)*( (kx*transpose(pjx[:,2]+pjx[:,1])) +  (ky*transpose(pjy[:,2]+pjy[:,1])) )).* ( (ky*transpose(pjx[:,2]-pjx[:,1]))  - (kx*transpose(pjy[:,2]-pjy[:,1])) ) 
   + sinc.( (kx*transpose(pjx[:,3]-pjx[:,2])) + (ky*transpose(pjy[:,3]-pjy[:,2])) ).*cis.(-T(pi)*( (kx*transpose(pjx[:,3]+pjx[:,2])) +  (ky*transpose(pjy[:,3]+pjy[:,2])) )).* ( (ky*transpose(pjx[:,3]-pjx[:,2]))  - (kx*transpose(pjy[:,3]-pjy[:,2])) )
   + sinc.( (kx*transpose(pjx[:,4]-pjx[:,3])) + (ky*transpose(pjy[:,4]-pjy[:,3])) ).*cis.(-T(pi)*( (kx*transpose(pjx[:,4]+pjx[:,3])) +  (ky*transpose(pjy[:,4]+pjy[:,3])) )).* ( (ky*transpose(pjx[:,4]-pjx[:,3]))  - (kx*transpose(pjy[:,4]-pjy[:,3])) ) 
   + sinc.( (kx*transpose(pjx[:,1]-pjx[:,4])) + (ky*transpose(pjy[:,1]-pjy[:,4])) ).*cis.(-T(pi)*( (kx*transpose(pjx[:,1]+pjx[:,4])) +  (ky*transpose(pjy[:,1]+pjy[:,4])) )).* ( (ky*transpose(pjx[:,1]-pjx[:,4]))  - (kx*transpose(pjy[:,1]-pjy[:,4])) )
   )./((kx.*kx+ky.*ky)));
  return polyft;
end

setup_polyft_single(data[1], stars[1])




function stcis(x1,x2,y1,y2,kx,ky)
    return sinc.(kx*(x2-x1) + ky*(y2-y1)).*cis.(-T(pi)*(kx*(x2+x1)+ ky*(y2+y1))).*(ky*(x2-x1)-kx*(y2-y1))
end

@views function setup_polyft_single_alt(data, star_epoch_geom; T=Float32)
pjx = star_epoch_geom.projx;
pjy = star_epoch_geom.projy;
kx =  data.uv[1,:] * T(-pi / (180*3600000));
ky =  data.uv[2,:] * T( pi / (180*3600000));
x1=Array(pjx[:,1]'); x2=Array(pjx[:,2]'); 
x3=Array(pjx[:,3]'); x4=Array(pjx[:,4]'); 
y1=Array(pjy[:,1]'); y2=Array(pjy[:,2]');
y3=Array(pjy[:,3]'); y4=Array(pjy[:,4]');
term1 = stcis(x1,x2,y1,y2,kx,ky);
term2 = stcis(x2,x3,y2,y3,kx,ky);
term3 = stcis(x3,x4,y3,y4,kx,ky);
term4 = stcis(x4,x1,y4,y1,kx,ky);
factor = -im*T(1/(2pi))./(kx.*kx+ky.*ky)
polyft = factor.*(term1+term2+term3+term4)
return polyft;
end

@views function setup_polyft_single_alt2(data, star_epoch_geom; T=Float32)
    pjx = star_epoch_geom.projx
    pjy = star_epoch_geom.projy
    kx =  data.uv[1,:] * T(-pi / (180*3600000))
    ky =  data.uv[2,:] * T( pi / (180*3600000))
    x = [Array(pjx[:,i])' for i in 1:4]
    y = [Array(pjy[:,i])' for i in 1:4]

    term1 = Threads.@spawn stcis(x[1], x[2], y[1], y[2], kx, ky)
    term2 = Threads.@spawn stcis(x[2], x[3], y[2], y[3], kx, ky)
    term3 = Threads.@spawn stcis(x[3], x[4], y[3], y[4], kx, ky)
    term4 = Threads.@spawn stcis(x[4], x[1], y[4], y[1], kx, ky)

    factor = -im * T(1/(2pi)) ./ (kx .* kx + ky .* ky)
    
    # Wait for the tasks to complete and get their results
    term1 = fetch(term1)
    term2 = fetch(term2)
    term3 = fetch(term3)
    term4 = fetch(term4)

    polyft = factor .* (term1 + term2 + term3 + term4)
    return polyft
end
