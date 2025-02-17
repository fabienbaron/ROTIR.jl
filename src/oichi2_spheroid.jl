using OptimPackNextGen
# The Fourier transform for polygonal surfaces
function poly_to_cvis(x, polyflux, polyft)
  flux = sum(polyflux.*x); # get the total flux
  cvis_model = polyft * x / flux;
end

function setupft_single(data, star_epoch_geom, ld = true)
# Polyflux is the weight of each Healpixel, proportional to the surface
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
for uu=1:data.nuv
  kx = -data.uv[1,uu] * (pi / 180.0) / 3600000.0;
  ky = data.uv[2,uu] * (pi / 180.0) / 3600000.0;
  for i=1:4 # note: danger, check definition sinc(x) = sin(pi*x)/(pi*x)
    polyft[uu,star_epoch_geom.index_quads_visible] += sinc.( (star_epoch_geom.projx[:,mod(i,4)+1]-star_epoch_geom.projx[:,i]).*kx + (star_epoch_geom.projy[:,mod(i,4)+1]-star_epoch_geom.projy[:,i]).*ky ).*cis.(-pi*( (star_epoch_geom.projx[:,mod(i,4)+1]+star_epoch_geom.projx[:,i]).*kx +  (star_epoch_geom.projy[:,mod(i,4)+1]+star_epoch_geom.projy[:,i]).*ky )).* ( (star_epoch_geom.projx[:,mod(i,4)+1]-star_epoch_geom.projx[:,i]).*ky  - (star_epoch_geom.projy[:,mod(i,4)+1]-star_epoch_geom.projy[:,i]).*kx )
  end
  polyft[uu,star_epoch_geom.index_quads_visible] /= im*2*pi*(kx*kx+ky*ky);
  # Add limb-darkening
  if ld == true
    polyft[uu,:] = polyft[uu,:].*star_epoch_geom.ldmap
  end
end
return polyflux, polyft
end


function setupft_single_alt(data, star_epoch_geom, ld = true)
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
  kx = -data.uv[1,:] * (pi / 180.0) / 3600000.0;
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


function setup_polygon_ft(data, star_epoch_geom, ld = true)
nepochs = size(data,1);
polyflux = Array{Array{Float64,1}}(undef,nepochs);
polyft = Array{Array{Complex{Float64},2}}(undef,nepochs);
for i=1:nepochs
    polyflux[i], polyft[i] = setupft_single_alt(data[i], star_epoch_geom[i], ld);
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
function spheroid_chi2_f(x, polyflux, polyft, data, verbose = true)
cvis_model = poly_to_cvis(x, polyflux, polyft);
# compute observables from all cvis
v2_model = cvis_to_v2(cvis_model, data.indx_v2);
t3_model, t3amp_model, t3phi_model = cvis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);
chi2_v2 = sum( ((v2_model - data.v2)./data.v2_err).^2);
chi2_t3amp = sum( ((t3amp_model - data.t3amp)./data.t3amp_err).^2);
chi2_t3phi = sum( (mod360(t3phi_model - data.t3phi)./data.t3phi_err).^2);
if verbose == true
  flux = sum(polyflux .*x);
  println("V2: ", chi2_v2/data.nv2, " T3A: ", chi2_t3amp/data.nt3amp, " T3P: ", chi2_t3phi/data.nt3phi," Flux: ", flux)
end
return chi2_v2 + chi2_t3amp + chi2_t3phi
end

function spheroid_chi2_fg(x, g, polyflux, polyft, data, verbose = true) # criterion function plus its gradient w/r x
nx2 = length(x);
cvis_model = poly_to_cvis(x, polyflux, polyft);
v2_model = cvis_to_v2(cvis_model, data.indx_v2);
t3_model, t3amp_model, t3phi_model = cvis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);
chi2_v2 = sum( ((v2_model - data.v2)./data.v2_err).^2);
chi2_t3amp = sum( ((t3amp_model - data.t3amp)./data.t3amp_err).^2);
chi2_t3phi = sum( (mod360(t3phi_model - data.t3phi)./data.t3phi_err).^2);
g_v2 = 4.0*sum(((v2_model-data.v2)./data.v2_err.^2).*real(conj(cvis_model[data.indx_v2]).*polyft[data.indx_v2,:]),1);
g_t3amp = 2.0*sum(((t3amp_model-data.t3amp)./data.t3amp_err.^2).*
                  (   real( conj(cvis_model[data.indx_t3_1]./abs.(cvis_model[data.indx_t3_1])).*polyft[data.indx_t3_1,:]).*abs.(cvis_model[data.indx_t3_2]).*abs.(cvis_model[data.indx_t3_3])       + real( conj(cvis_model[data.indx_t3_2]./abs.(cvis_model[data.indx_t3_2])).*polyft[data.indx_t3_2,:]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_3])+ real( conj(cvis_model[data.indx_t3_3]./abs.(cvis_model[data.indx_t3_3])).*polyft[data.indx_t3_3,:]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_2])),1);

t3model_der = polyft[data.indx_t3_1,:].*cvis_model[data.indx_t3_2].*cvis_model[data.indx_t3_3] + polyft[data.indx_t3_2,:].*cvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_3] + polyft[data.indx_t3_3,:].*cvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_2];
g_t3phi=360.0/pi*sum(((mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2)./abs2.(t3_model)).*(-imag(t3_model).*real(t3model_der)+real(t3_model).*imag(t3model_der)),1);
g[1:nx2] = g_v2 + g_t3amp + g_t3phi;
flux = sum(polyflux .*x);
g[1:nx2] = (g - sum(x.*g)*polyflux / flux ) / flux; # gradient correction to take into account the non-normalization
if verbose == true
  println("V2: ", chi2_v2/data.nv2, " T3A: ", chi2_t3amp/data.nt3amp, " T3P: ", chi2_t3phi/data.nt3phi," Flux: ", flux)
end
return chi2_v2 + chi2_t3amp + chi2_t3phi
end

function spheroid_chi2_fg_alt(x, g, polyflux, polyft, data; verbose::Bool = true) # criterion function plus its gradient w/r x, alternative version
  # this may be faster or slower depending on number of data points vs number of pixels
  nx2 = length(x);
  cvis_model = poly_to_cvis(x, polyflux, polyft);
  v2_model = cvis_to_v2(cvis_model, data.indx_v2);
  t3_model, t3amp_model, t3phi_model = cvis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);
  chi2_v2 = sum( ((v2_model - data.v2)./data.v2_err).^2);
  chi2_t3amp = sum( ((t3amp_model - data.t3amp)./data.t3amp_err).^2);
  chi2_t3phi = sum( (mod360(t3phi_model - data.t3phi)./data.t3phi_err).^2);
  g_v2 = real(transpose(polyft[data.indx_v2,:])*(4*((v2_model-data.v2)./data.v2_err.^2).*conj(cvis_model[data.indx_v2])));
  g_t3amp = real(transpose(polyft[data.indx_t3_1,:]) *(2.0*((t3amp_model-data.t3amp)./data.t3amp_err.^2).*conj(cvis_model[data.indx_t3_1])./abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_2]).*abs.(cvis_model[data.indx_t3_3]) ))+real(transpose(polyft[data.indx_t3_2,:])*(2.0*((t3amp_model-data.t3amp)./data.t3amp_err.^2).*conj(cvis_model[data.indx_t3_2])./abs.(cvis_model[data.indx_t3_2]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_3]) ))+real(transpose(polyft[data.indx_t3_3,:])*(2.0*((t3amp_model-data.t3amp)./data.t3amp_err.^2).*conj(cvis_model[data.indx_t3_3])./abs.(cvis_model[data.indx_t3_3]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_2]) ));
  g_t3phi = 360.0/pi*imag(transpose(polyft[data.indx_t3_1,:])*(((mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2)./abs2.(t3_model)).*cvis_model[data.indx_t3_2].*cvis_model[data.indx_t3_3].*conj(t3_model))+transpose(polyft[data.indx_t3_2,:])*(((mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2)./abs2.(t3_model)).*cvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_3].*conj(t3_model))+transpose(polyft[data.indx_t3_3,:])*(((mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2)./abs2.(t3_model)).*cvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_2].*conj(t3_model)));
  g[1:nx2] = g_v2 + g_t3amp + g_t3phi;
  flux = sum(polyflux .*x);
  g[1:nx2] = (g .- sum(x.*g)*polyflux / flux ) / flux; # gradient correction to take into account the non-normalization
  if verbose == true
    println("V2: ", chi2_v2/data.nv2, " T3A: ", chi2_t3amp/data.nt3amp, " T3P: ", chi2_t3phi/data.nt3phi," Flux: ", flux)
  end
  return chi2_v2 + chi2_t3amp + chi2_t3phi
end


function spheroid_chi2_allepochs_fg(x, g, epochs_weights, polyflux, polyft, data)
f = 0;
g[:] = 0;
npix = size(x);
singleepoch_g = zeros(Float64, npix);
for i=1:nepochs # weighted sum -- should probably do the computation in parallel
  f += epochs_weights[i]*spheroid_chi2_fg_alt(x, singleepoch_g, polyflux[i], polyft[i], data[i], verbose=true);
  g[:] += epochs_weights[i]*singleepoch_g;
end
println("All epochs, weighted chi2: ", f, "\n");
return f;
end

function spheroid_chi2_allepochs_f(x, epochs_weights, polyflux, polyft, data)
f = 0;
npix = size(x);
for i=1:nepochs # weighted sum -- should probably do the computation in parallel
  f += epochs_weights[i]*spheroid_chi2_f(x, polyflux[i], polyft[i], data[i], true);
end
println("All epochs, weighted chi2: ", f, "\n");
return f;
end


function proj_positivity(ztilde)
z = copy(ztilde)
z[ztilde.>0]=0
return z
end


function spheroid_crit_multiepochs_fg(x, g, polyflux, polyft, data; regularizers=[], epochs_weights = [], verb=true)
chi2_f = 0.0;
g[:] .= 0.0;
npix = size(x);
singleepoch_g = zeros(Float64, npix);
nepochs = length(data)
if epochs_weights == []
    epochs_weights = ones(Float64, nepochs)/nepochs
end
for i=1:nepochs # weighted sum -- should probably do the computation in parallel
  chi2_f += epochs_weights[i]*spheroid_chi2_fg_alt(x, singleepoch_g, polyflux[i], polyft[i], data[i], verbose=verb);
  g[:] += epochs_weights[i]*singleepoch_g;
end
reg_f = spheroid_regularization(x, g, regularizers=regularizers, verb = verb); #note: this adds to g
if (verb==true)
    println("All epochs, weighted chi2: ", chi2_f, " crit: ", chi2_f + reg_f);
end
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

function spheroid_total_variation_latitude_fg(x, tv_g, tvinfo; ϵ = 1e-13, verb = true)
  # Add total variation regularization only in the latitudinal direction
  #tvinfo: neighbors,south_neighbors,west_neighbors,south_neighbors_reverse,west_neighbors_reverse
  npix = length(x)
  xs = x[tvinfo[2]];
  xw = x[tvinfo[3]];
  tv_f = sum(sqrt.( (x-xw).^2) .+ ϵ )
  tv_g[1:npix] = (x-xw)./(sqrt.((x-xs).^2 + (x-xw).^2) .+ ϵ) # or just 1.0?
  for j=1:length(x)
    l = tvinfo[5][j]
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

function spheroid_total_variation_latitude_f(x, tvinfo; ϵ = 1e-13, verb = true)
  # Add total variation regularization
  #tvinfo: neighbors,south_neighbors,west_neighbors,south_neighbors_reverse,west_neighbors_reverse
  xw = x[tvinfo[3]];
  tv_f = sum(sqrt.( (x-xw).^2) .+ ϵ )
  if verb == true
      println("TV: ", tv_f);
end
  return tv_f
end


function spheroid_l2_fg(x, g; verb = true)
l2f = sum(abs.(x.-sum(x)/length(x)))
g[:] = sign.(x.-sum(x)/length(x))
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

function max_entropy_fg(x, g; verb=true, ϵ=1e-9)
  # mmap = sum(x) / length(x)
  # nmap = x ./ mmap
  # reg_f = sum(nmap .* log.(nmap))
  # g[:] = sum(-nmap ./ sum(x) .* (log.(nmap) .+ 1.0)) .+ (log.(nmap) .+ 1)./mmap
  xm = x ./ (mean(x) + ϵ)
  reg_f = sum(xm .* log.(xm .+ ϵ))
  g[:] = ((mean(x) .- (x ./ length(x))) ./ (mean(x)^2 + ϵ)) .* (log.(xm .+ ϵ) .+ 1)
  if verb == true
    println(" MEM: ", reg_f)
  end
  return reg_f
end

function spheroid_regularization(x,g; printcolor = :black, regularizers=[], verb=false)
    reg_f = 0.0;
    for ireg =1:length(regularizers)
      x_sub = x[regularizers[ireg][4]] # take the pixel subset if needed (example: only regularize visible pixels)
      temp_g = similar(x_sub);
      if regularizers[ireg][1] == "tv"
        reg_f += regularizers[ireg][2]*spheroid_total_variation_fg(x_sub, temp_g, regularizers[ireg][3], verb = verb);
      elseif regularizers[ireg][1] == "tv_lat"
        reg_f += regularizers[ireg][2]*spheroid_total_variation_latitude_fg(x_sub, temp_g, regularizers[ireg][3], verb = verb);
      elseif regularizers[ireg][1] == "l2"
        reg_f += regularizers[ireg][2]*spheroid_l2_fg(x_sub, temp_g, verb = verb);
      elseif regularizers[ireg][1] == "bias"
        reg_f += regularizers[ireg][2]*spheroid_harmon_bias_fg(x_sub, temp_g, regularizers[ireg][3], verb = verb );
      elseif regularizers[ireg][1] == "mem"
        reg_f += regularizers[ireg][2]*max_entropy_fg(x_sub, temp_g, verb = verb);
      end
      g[regularizers[ireg][4]] .+= regularizers[ireg][2].*temp_g
    end
    if verb == true
      println("\n");
    end
  return reg_f
end


using OptimPackNextGen
function spheroid_oi_reconstruct(x_start::Array{Float64,1}, data::Array{OIdata,1}, polyflux::Array{Array{Float64,1},1}, polyft::Array{Array{Complex{Float64},2},1}; epochs_weights =[], printcolor= [], verb = true, lower=0, upper=0, maxiter = 100, regularizers =[])
  x_sol = [];
  crit_imaging = (x,g)->spheroid_crit_multiepochs_fg(x, g, polyflux, polyft, data, regularizers=regularizers, epochs_weights=  epochs_weights, verb = verb);
  # if (upper == 0)
  #   x_sol = OptimPackNextGen.vmlmb(crit_imaging, x_start, verb=verb, lower=0, maxiter=maxiter, blmvm=false, gtol=(0,1e-8));
  # elseif (upper > 0)
  #   x_sol = OptimPackNextGen.vmlmb(crit_imaging, x_start, verb=verb, lower=0, upper=upper, maxiter=maxiter, blmvm=false, gtol=(0,1e-8));
  # end
  x_sol = OptimPackNextGen.vmlmb(crit_imaging, x_start, verb=verb, lower=0, maxiter=maxiter, blmvm=false, gtol=(0,1e-8));
  dummy = similar(x_sol);
  crit_opt = crit_imaging(x_sol,dummy);
  return x_sol, crit_opt
end

function oi_reconstruct_mutitemporal(x_start::Array{Float64,1}, data::Array{OIdata,1}, polyflux::Array{Array{Float64,1},1}, polyft::Array{Array{Complex{Float64},2},1}; epochs_weights =[], printcolor= [], verb = true, maxiter = 100, regularizers =[])
  crit_imaging = (x,g)->oi_multitemporal_fg(x, g, polyflux, polyft, data, regularizers=regularizers, epochs_weights=  epochs_weights);
  x_sol = OptimPackNextGen.vmlmb(crit_imaging, x_start, verb=verb, lower=0, maxiter=maxiter, blmvm=false, gtol=(0,1e-8));
  return reshape(x_sol, size(polyflux[1],1), length(data))
end

function oi_multitemporal_fg(x, g, polyflux, polyft, data; regularizers=[], epochs_weights = [], verb=true)
  # Explanation of the following: optimpack optimizes a vector, but we want an array of images
  npix = size(polyflux[1],1);
  chi2_f = 0.0;
  g[:] .= 0.0;
  #temp_g = deepcopy(g[:]);
  #npix = size(x);
  nepochs = length(data);
  if epochs_weights == []
    epochs_weights = ones(Float64, nepochs)/nepochs;
  end

  #singleepoch_g = zeros(Float64, npix);
  #for i=1:nepochs # weighted sum -- should probably do the computation in parallel
  #  chi2_f += epochs_weights[i]*spheroid_chi2_fg_alt(x[1:npix], singleepoch_g, polyflux[i], polyft[i], data[i], verbose=verb);
  #  g[1:npix] += epochs_weights[i]*singleepoch_g;
  #end


  for i=1:nepochs # weighted sum -- in the future, do the computation in parallel
    tslice = 1+(i-1)*npix:i*npix; # temporal slice
    subg = zeros(Float64, npix);
    chi2_f += epochs_weights[i]*spheroid_chi2_fg_alt(x[tslice], subg, polyflux[i], polyft[i], data[i], verbose=verb);
    g[tslice] = epochs_weights[i]*subg;
    #x[tslice] = x[1:npix];
    #g[tslice] = g[1:npix];
  end

  # cross temporal regularization -- weight needs to be defined in the "regularizers" variable
  if length(regularizers)>nepochs
    if (regularizers[nepochs+1][1][1] == "temporal_tvsq")  & (regularizers[nepochs+1][1][2] > 0.0) & (nepochs>1)
      y = reshape(x,(npix,nepochs))
      temporalf = sum( (y[:,2:end]-y[:,1:end-1]).^2 )
      tv_g = Array{Float64}(undef, npix,nepochs)
      if nepochs>2
         tv_g[:,1] = 2*(y[:,1] - y[:,2])
         tv_g[:,2:end-1] = 4*y[:,2:end-1]-2*(y[:,1:end-2]+y[:,3:end])
         tv_g[:,end] = 2*(y[:,end] - y[:,end-1])
      else
         tv_g[:,1] = 2*(y[:,1]-y[:,2]);
         tv_g[:,2] = 2*(y[:,2]-y[:,1]);
      end
      chi2_f += regularizers[nepochs+1][1][2]*temporalf
      g[:] += regularizers[nepochs+1][1][2]*vec(tv_g);
      if verb == true
           printstyled("Temporal regularization: $temporalf\n", color=:yellow)
      end
    end
  end
  reg_f = spheroid_regularization(x, g, regularizers=regularizers[1:nepochs], verb = verb); # Note: adds to g
  return chi2_f + reg_f;
end
