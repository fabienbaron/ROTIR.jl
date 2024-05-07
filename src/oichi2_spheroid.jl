using OptimPackNextGen

# The Fourier transform for polygonal surfaces
# f0 = unresolved flux
function poly_to_cvis(x, polyflux, polyft; f0=0.0)
  flux = sum(polyflux.*x); # star flux
  return (1.0-f0)*(polyft * x / flux);
end

function setupft_single(data, star_epoch_geom; ld = true)
  # Polyflux is the weight of each pixel, proportional to the surface
  polyflux = zeros(Float64, star_epoch_geom.npix);
  polyflux[star_epoch_geom.index_quads_visible] =0.5*(star_epoch_geom.projx[:,1].*star_epoch_geom.projy[:,2]
  - star_epoch_geom.projx[:,2].*star_epoch_geom.projy[:,1]
  + star_epoch_geom.projx[:,2].*star_epoch_geom.projy[:,3]
  - star_epoch_geom.projx[:,3].*star_epoch_geom.projy[:,2]
  + star_epoch_geom.projx[:,3].*star_epoch_geom.projy[:,4]
  - star_epoch_geom.projx[:,4].*star_epoch_geom.projy[:,3]
  + star_epoch_geom.projx[:,4].*star_epoch_geom.projy[:,1]
  - star_epoch_geom.projx[:,1].*star_epoch_geom.projy[:,4]);

  if ld == true
  # Take into account limb-darkening
  polyflux = polyflux.*star_epoch_geom.ldmap;
  end

  # and now the core of the quad/polygon FT
  # this assumes the (x,y) coordinates of the quads are in milliarcseconds
  polyft = zeros(Complex{Float64}, data.nuv, star_epoch_geom.npix); #note: size = npix, but we will fill only the nquads_visible ones
  #polyft = SharedMatrix{Complex{Float64}}(data.nuv, star_epoch_geom.npix);
  kx = data.uv[1,:] * (pi / 180.0) / 3600000.0;
  ky = data.uv[2,:] * (pi / 180.0) / 3600000.0;
  # note: danger, check definition sinc(x) = sin(pi*x)/(pi*x)
  polyft[:,star_epoch_geom.index_quads_visible] = (sinc.( (kx*transpose(star_epoch_geom.projx[:,2]-star_epoch_geom.projx[:,1])) + (ky*transpose(star_epoch_geom.projy[:,2]-star_epoch_geom.projy[:,1])) ).*cis.(-pi*( (kx*transpose(star_epoch_geom.projx[:,2]+star_epoch_geom.projx[:,1])) +  (ky*transpose(star_epoch_geom.projy[:,2]+star_epoch_geom.projy[:,1])) )).* ( (ky*transpose(star_epoch_geom.projx[:,2]-star_epoch_geom.projx[:,1]))  - (kx*transpose(star_epoch_geom.projy[:,2]-star_epoch_geom.projy[:,1])) )
    + sinc.( (kx*transpose(star_epoch_geom.projx[:,3]-star_epoch_geom.projx[:,2])) + (ky*transpose(star_epoch_geom.projy[:,3]-star_epoch_geom.projy[:,2])) ).*cis.(-pi*( (kx*transpose(star_epoch_geom.projx[:,3]+star_epoch_geom.projx[:,2])) +  (ky*transpose(star_epoch_geom.projy[:,3]+star_epoch_geom.projy[:,2])) )).* ( (ky*transpose(star_epoch_geom.projx[:,3]-star_epoch_geom.projx[:,2]))  - (kx*transpose(star_epoch_geom.projy[:,3]-star_epoch_geom.projy[:,2])) )
    + sinc.( (kx*transpose(star_epoch_geom.projx[:,4]-star_epoch_geom.projx[:,3])) + (ky*transpose(star_epoch_geom.projy[:,4]-star_epoch_geom.projy[:,3])) ).*cis.(-pi*( (kx*transpose(star_epoch_geom.projx[:,4]+star_epoch_geom.projx[:,3])) +  (ky*transpose(star_epoch_geom.projy[:,4]+star_epoch_geom.projy[:,3])) )).* ( (ky*transpose(star_epoch_geom.projx[:,4]-star_epoch_geom.projx[:,3]))  - (kx*transpose(star_epoch_geom.projy[:,4]-star_epoch_geom.projy[:,3])) )
    + sinc.( (kx*transpose(star_epoch_geom.projx[:,1]-star_epoch_geom.projx[:,4])) + (ky*transpose(star_epoch_geom.projy[:,1]-star_epoch_geom.projy[:,4])) ).*cis.(-pi*( (kx*transpose(star_epoch_geom.projx[:,1]+star_epoch_geom.projx[:,4])) +  (ky*transpose(star_epoch_geom.projy[:,1]+star_epoch_geom.projy[:,4])) )).* ( (ky*transpose(star_epoch_geom.projx[:,1]-star_epoch_geom.projx[:,4]))  - (kx*transpose(star_epoch_geom.projy[:,1]-star_epoch_geom.projy[:,4])) ));
  polyft[:,star_epoch_geom.index_quads_visible] = broadcast(/,polyft[:,star_epoch_geom.index_quads_visible],im*2*pi*(kx.*kx+ky.*ky));
  # Add limb-darkening
  if ld == true
    polyft = broadcast(*,polyft',star_epoch_geom.ldmap)';
  end
  return polyflux, polyft
end


function setup_polygon_ft(data, star_epoch_geom; ld = true)
nepochs = size(data,1);
polyflux = Array{Array{Float64,1}}(undef,nepochs);
polyft = Array{Array{Complex{Float64},2}}(undef,nepochs);
for i=1:nepochs
    polyflux[i], polyft[i] = setupft_single(data[i], star_epoch_geom[i], ld=ld);
end
return polyflux, polyft
end

function mod360(x)
  mod.(mod.(x.+180.0,360.0).+360.0, 360.0) .- 180.0
end

function cvis_to_v2(cvis, indx)
  v2_model = abs2.(cvis[indx]);
end

function cvis_to_t3(cvis, indx1, indx2, indx3)
  t3 = cvis[indx1].*cvis[indx2].*cvis[indx3];
  t3amp = abs.(t3);
  t3phi = angle.(t3)*180.0/pi;
  return t3, t3amp, t3phi
end

function cvis_to_t4(cvis, indx1, indx2, indx3, indx4)
  t4 = cvis[indx1].*cvis[indx2]./(cvis[indx3]*conj(cvis[indx4]));
  t4amp = abs.(t4);
  t4phi = angle.(t4)*180.0/pi;
  return t4, t4amp, t4phi
end

function observables(xx, polyflux, polyft, data)
  x  = xx[2:end]; #temperature map
  f0 = xx[1];   # unresolved flux
  cvis_model = poly_to_cvis(x, polyflux, polyft, f0=f0);
  # compute observables from all cvis
  v2_model = cvis_to_v2(cvis_model, data.indx_v2);
  ~, t3amp_model, t3phi_model = cvis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);
  return v2_model, t3amp_model, t3phi_model
  end
  

function chi2s(xx, polyflux, polyft, data; verbose = true)
  v2_model, t3amp_model, t3phi_model = observables(xx, polyflux, polyft, data);
  chi2_v2 = sum( ((v2_model - data.v2)./data.v2_err).^2);
  chi2_t3amp = sum( ((t3amp_model - data.t3amp)./data.t3amp_err).^2);
  chi2_t3phi = sum( (mod360(t3phi_model - data.t3phi)./data.t3phi_err).^2);
  if verbose == true
    println("V2: ", chi2_v2/data.nv2, " T3A: ", chi2_t3amp/data.nt3amp, " T3P: ", chi2_t3phi/data.nt3phi)
  end
  return chi2_v2, chi2_t3amp, chi2_t3phi
end
    


function spheroid_chi2_f(xx, polyflux, polyft, data; verbose::Bool = true)
x  = xx[2:end]; #temperature map
f0 = xx[1];   # unresolved flux
cvis_model = poly_to_cvis(x, polyflux, polyft, f0=f0);
# compute observables from all cvis
v2_model = cvis_to_v2(cvis_model, data.indx_v2);
t3_model, t3amp_model, t3phi_model = cvis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);
chi2_v2 = sum( ((v2_model - data.v2)./data.v2_err).^2);
chi2_t3amp = sum( ((t3amp_model - data.t3amp)./data.t3amp_err).^2);
chi2_t3phi = sum( (mod360(t3phi_model - data.t3phi)./data.t3phi_err).^2);
if verbose == true
  flux_star = sum(polyflux .*x);
  println("V2: ", chi2_v2/data.nv2, " T3A: ", chi2_t3amp/data.nt3amp, " T3P: ", chi2_t3phi/data.nt3phi," Raw Star Flux: ", flux_star, " Unresolved flux fraction: ", f0)
end
return chi2_v2 + chi2_t3amp + chi2_t3phi
end

@views function spheroid_chi2_fg(xx, g, polyflux, polyft, data; verbose::Bool = true) 
  x  = xx[2:end]; #temperature map
  f0 = xx[1];   # unresolved flux
  cvis_model = poly_to_cvis(x, polyflux, polyft, f0=f0); # cplx vis including unresolved flux
  v2_model = cvis_to_v2(cvis_model, data.indx_v2);
  t3_model, t3amp_model, t3phi_model = cvis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);
  chi2_v2 = sum( ((v2_model - data.v2)./data.v2_err).^2);
  chi2_t3amp = sum( ((t3amp_model - data.t3amp)./data.t3amp_err).^2);
  chi2_t3phi = sum( (mod360(t3phi_model - data.t3phi)./data.t3phi_err).^2);
  # # # Grad chi2 w/r temperatute map x
  #  g_v2 = real(transpose(polyft[data.indx_v2,:])*(4*((v2_model-data.v2)./data.v2_err.^2).*conj(cvis_model[data.indx_v2])));
  #  g_t3amp = real(transpose(polyft[data.indx_t3_1,:]) *(2.0*((t3amp_model-data.t3amp)./data.t3amp_err.^2).*conj(cvis_model[data.indx_t3_1])./abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_2]).*abs.(cvis_model[data.indx_t3_3]) ))+real(transpose(polyft[data.indx_t3_2,:])*(2.0*((t3amp_model-data.t3amp)./data.t3amp_err.^2).*conj(cvis_model[data.indx_t3_2])./abs.(cvis_model[data.indx_t3_2]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_3]) ))+real(transpose(polyft[data.indx_t3_3,:])*(2.0*((t3amp_model-data.t3amp)./data.t3amp_err.^2).*conj(cvis_model[data.indx_t3_3])./abs.(cvis_model[data.indx_t3_3]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_2]) ));
  #  g_t3phi = 360.0/pi*imag(transpose(polyft[data.indx_t3_1,:])*(((mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2)./abs2.(t3_model)).*cvis_model[data.indx_t3_2].*cvis_model[data.indx_t3_3].*conj(t3_model))+transpose(polyft[data.indx_t3_2,:])*(((mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2)./abs2.(t3_model)).*cvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_3].*conj(t3_model))+transpose(polyft[data.indx_t3_3,:])*(((mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2)./abs2.(t3_model)).*cvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_2].*conj(t3_model)));
  #  g_temp = g_v2 + g_t3amp + g_t3phi;
  #  flux = sum(polyflux .*x); # star flux
  #  g1 = (g_temp .- sum(x.*g_temp)*polyflux / flux ) / flux; # gradient correction to take into account the non-normalization

  # S = polyflux'
  # V = Hx./ S x 
  # dV/dx =  ( Sx .* H' - Hx .* S' ) ./ (S x)^2
  #       =  H' / Sx - Hx S' / flux^2 
  #       = transpose(polyft)/flux - (polyft*x)*polyflux/flux^2
  # Could use V 
  flux_star = sum(polyflux .*x); # star flux
  dV = transpose((polyft - (polyft*x/flux_star).*transpose(polyflux))/flux_star)
  g_v2 = real(dV[:,data.indx_v2]*(4*((v2_model-data.v2)./data.v2_err.^2).*conj(cvis_model[data.indx_v2])));
  g_t3amp = real(  dV[:,data.indx_t3_1] *(2.0*((t3amp_model-data.t3amp)./data.t3amp_err.^2).*conj(cvis_model[data.indx_t3_1])./abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_2]).*abs.(cvis_model[data.indx_t3_3]) ))+real( dV[:,data.indx_t3_2]*(2.0*((t3amp_model-data.t3amp)./data.t3amp_err.^2).*conj(cvis_model[data.indx_t3_2])./abs.(cvis_model[data.indx_t3_2]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_3]) ))+real( dV[:,data.indx_t3_3]*(2.0*((t3amp_model-data.t3amp)./data.t3amp_err.^2).*conj(cvis_model[data.indx_t3_3])./abs.(cvis_model[data.indx_t3_3]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_2]) ));
  g_t3phi = 360.0/pi*imag(dV[:,data.indx_t3_1] *(((mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2)./abs2.(t3_model)).*cvis_model[data.indx_t3_2].*cvis_model[data.indx_t3_3].*conj(t3_model))
                         +dV[:,data.indx_t3_2]*(((mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2)./abs2.(t3_model)).*cvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_3].*conj(t3_model))
                         +dV[:,data.indx_t3_3]*(((mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2)./abs2.(t3_model)).*cvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_2].*conj(t3_model)));
  g[2:end] = (1-f0)*(g_v2 + g_t3amp + g_t3phi);
  g[1] = 0;
  if verbose == true
    println("V2: ", chi2_v2/data.nv2, " T3A: ", chi2_t3amp/data.nt3amp, " T3P: ", chi2_t3phi/data.nt3phi," Raw Star Flux: ", flux_star, " Unresolved flux fraction: ", f0)
  end
  return chi2_v2 + chi2_t3amp + chi2_t3phi
end


  
# dV = cvis_model/(f0-1)
# g_v2 = sum(real(dV[data.indx_v2].*(4*((v2_model-data.v2)./data.v2_err.^2).*conj(cvis_model[data.indx_v2]))));
# g_t3amp = sum(real(dV[data.indx_t3_1].*(2.0*((t3amp_model-data.t3amp)./data.t3amp_err.^2).*conj(cvis_model[data.indx_t3_1])./abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_2]).*abs.(cvis_model[data.indx_t3_3]) ))+real(dV[data.indx_t3_2].*(2.0*((t3amp_model-data.t3amp)./data.t3amp_err.^2).*conj(cvis_model[data.indx_t3_2])./abs.(cvis_model[data.indx_t3_2]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_3]) ))+real(dV[data.indx_t3_3].*(2.0*((t3amp_model-data.t3amp)./data.t3amp_err.^2).*conj(cvis_model[data.indx_t3_3])./abs.(cvis_model[data.indx_t3_3]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_2]))));
# g_t3phi = sum(360.0/pi*imag(dV[data.indx_t3_1].*(((mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2)./abs2.(t3_model)).*cvis_model[data.indx_t3_2].*cvis_model[data.indx_t3_3].*conj(t3_model))
#                        +dV[data.indx_t3_2].*(((mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2)./abs2.(t3_model)).*cvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_3].*conj(t3_model))
#                        +dV[data.indx_t3_3].*(((mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2)./abs2.(t3_model)).*cvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_2].*conj(t3_model))));
# g[1] = g_v2 + g_t3amp + g_t3phi; # gradient w/r unresolved flux


function proj_positivity(ztilde)
  z = copy(ztilde)
  z[ztilde.>0]=0
return z
end


function spheroid_crit_multiepochs_fg(xx, gg, polyflux, polyft, data; regularizers=[], epochs_weights = [], verb=true )
chi2_f = 0.0;
gg[:] .= 0.0;
singleepoch_g = zeros(Float64, length(xx));
nepochs = length(data)
if epochs_weights == []
    epochs_weights = ones(Float64, nepochs)/nepochs
end
for i=1:nepochs # weighted sum -- should probably do the computation in parallel
  chi2_f += epochs_weights[i]*spheroid_chi2_fg(xx, singleepoch_g, polyflux[i], polyft[i], data[i], verbose=verb);
  gg[:] += epochs_weights[i]*singleepoch_g;
end
  println("All epochs, weighted chi2: ", chi2_f);
# Map regularization
reg_g = zeros(Float64, length(xx)-1);
reg_f = spheroid_regularization(xx[2:end], reg_g, regularizers=regularizers, verb = verb); #note: adds to reg_g
gg[2:end] += reg_g
return chi2_f + reg_f;
end


function spheroid_total_variation_fg(x, tv_g, tvinfo; ϵ = 1e-13, verb = true)
  # Add total variation regularization
  #tvinfo: neighbors,south_neighbors,west_neighbors,south_neighbors_reverse,west_neighbors_reverse
  npix = length(x)
  xs = x[tvinfo[2]];
  xw = x[tvinfo[3]];
  tv_f = sum(sqrt.( (x-xs).^2 + (x-xw).^2) .+ ϵ )
  tv_g[1:npix] = (2*x-xs-xw)./(sqrt.( (x-xs).^2 + (x-xw).^2) .+ ϵ)
  for j=1:length(x)
    k = tvinfo[4][j]
    l = tvinfo[5][j]
    if length(k)>0
      kk=k[1]
      tv_g[j] += (x[j]-x[kk])/sqrt((x[kk]-x[j])^2+(x[kk]-x[tvinfo[3][kk]])^2 +ϵ)
    end
    if length(l)>0
      ll=l[1]
      tv_g[j] += (x[j]-x[ll])/sqrt((x[ll]-x[tvinfo[2][ll]])^2+(x[ll]-x[j])^2 +ϵ)
    end
   end
if verb == true
      println("TV: ", tv_f);
end
  return tv_f
end


function spheroid_total_variation_f(x, tvinfo; ϵ = 1e-13, verb = true)
  # Add total variation regularization
  #tvinfo: neighbors,south_neighbors,west_neighbors,south_neighbors_reverse,west_neighbors_reverse
  xs = x[tvinfo[2]];
  xw = x[tvinfo[3]];
  tv_f = sum(sqrt.( (x-xs).^2 + (x-xw).^2) .+ ϵ )
  if verb == true
      println("TV: ", tv_f);
end
  return tv_f
end

function spheroid_total_variation2_fg(x, tv_g, tvinfo; ϵ = 1e-13, verb = true)
  npix = length(x)
  tv_f = norm(tvinfo[6]*x)^2
  tv_g[:] = 2*(tvinfo[7]*x)
  if verb == true
      println("TV2: ", tv_f);
  end
  return tv_f
end

function spheroid_l2_fg(x, g; verb = true)
l2f = sum(abs.(x .-sum(x)/length(x)))
g[:] = sign.(x .-sum(x)/length(x))
if verb == true
println(" L2: ", l2f);
end
return l2f;
end

function spheroid_harmon_bias_fg(x, g, B::Float64; verb = true)
n = length(x);
avg_x = mean(x);
bcorr = (B-1.0)*Int.((x.-avg_x).>0).+1.0
#reg_f = sum(bcorr.*(x.-avg_x).^2);
reg_f = sum(bcorr.*(x.-avg_x).^2)/n;
#reg_g = 2*(n-1)/n*bcorr.*(x.-avg_x)
reg_g = 2/n*bcorr.*(x.-avg_x)
g[:] = reg_g .- mean(reg_g); # necessary ?

if verb == true
println(" Bias: ", reg_f);
end
return reg_f;
end

function spheroid_regularization(x,g; printcolor = :black, regularizers=[], verb=false)
    reg_f = 0.0;
    for ireg =1:length(regularizers)
        x_sub = x[regularizers[ireg][4]] # take the pixel subset if needed (example: only regularize visible pixels)
        temp_g = similar(x_sub);
        if regularizers[ireg][1] == "tv"
            reg_f += regularizers[ireg][2]*spheroid_total_variation_fg(x_sub, temp_g, regularizers[ireg][3], verb = verb);
        elseif regularizers[ireg][1] == "tv2"
            reg_f += regularizers[ireg][2]*spheroid_total_variation2_fg(x_sub, temp_g, regularizers[ireg][3], verb = verb);
        elseif regularizers[ireg][1] == "l2"
            reg_f += regularizers[ireg][2]*spheroid_l2_fg(x_sub, temp_g, verb = verb);
        elseif regularizers[ireg][1] == "bias"
            reg_f += regularizers[ireg][2]*spheroid_harmon_bias_fg(x_sub, temp_g, regularizers[ireg][3], verb = verb );
        end
        g[regularizers[ireg][4]] += regularizers[ireg][2]*temp_g
    end
    if verb == true
        println("\n");
    end
  return reg_f
end

using NLopt
function optimize_bg_flux(xx,polyflux, polyft, data)
  x=xx[2:end]
  f0_start = xx[1]
  f = (f0,dummy)->spheroid_chi2_f(vcat(f0,x), polyflux[1], polyft[1], data[1], verbose = true);
  optimizer = Opt(:LN_NELDERMEAD, 1);
  min_objective!(optimizer, f);
  lower_bounds!(optimizer, [0.0]);
  upper_bounds!(optimizer, [0.2]);
  minchi2,params_opt,ret = optimize(optimizer, [f0_start]);
  println("Unresolved flux -- start=", f0_start, " optimized= ", params_opt[1])
  return vcat(params_opt[1], x)
end



using OptimPackNextGen
function spheroid_oi_reconstruct(xx_start::Array{Float64,1}, data::Array{OIdata,1}, polyflux::Array{Array{Float64,1},1}, polyft::Array{Array{Complex{Float64},2},1}; lower=0.0, upper=1e99 , epochs_weights =[], printcolor= [], verb = true, maxiter = 100, regularizers =[], bgflux=false )
  xx_sol = [];
  crit_imaging = (x,g)->spheroid_crit_multiepochs_fg(x, g, polyflux, polyft, data, regularizers=regularizers, epochs_weights=  epochs_weights, verb=verb);
  xx_sol = OptimPackNextGen.vmlmb(crit_imaging, xx_start, verb=verb, lower=lower, upper=upper,  maxiter=maxiter, blmvm=false, gtol=(0,1e-8));  
    return xx_sol
end
