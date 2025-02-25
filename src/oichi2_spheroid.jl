using OptimPackNextGen
# The Fourier transform for polygonal surfaces
function poly_to_cvis(x, polyflux, polyft)
  flux = sum(polyflux.*x); # get the total flux
  cvis_model = polyft * x / flux;
end

function poly_to_cvis(x, star)
   flux = dot(star.polyflux, x[star.index_quads_visible]); # get the total flux
   return star.polyft * x[star.index_quads_visible] / flux;
end

function poly_to_flux(x, star)
   flux = dot(star.polyflux,x[star.index_quads_visible]); # get the total flux
   return flux;
end


@views function setup_oi!(data, stars)
  nepochs = size(data,1);
  T = eltype(stars[1].vertices_xyz);
  #polyflux = Array{Array{T,1}}(undef,nepochs);
  #polyft = Array{Array{Complex{T},2}}(undef,nepochs);
  Threads.@threads for i=1:nepochs
       stars[i].polyflux = setup_polyflux_single(stars[i])
       stars[i].polyft = setup_polyft_single(data[i], stars[i]);
     end
end

@views function setup_polygon_ft(data, star_epoch_geom)
  nepochs = size(data,1);
  T = eltype(star_epoch_geom[1].vertices_xyz);
  polyflux = Array{Array{T,1}}(undef,nepochs);
  polyft = Array{Array{Complex{T},2}}(undef,nepochs);
  Threads.@threads for i=1:nepochs
       polyflux[i] = setup_polyflux_single(star_epoch_geom[i])
       polyft[i] = setup_polyft_single(data[i], star_epoch_geom[i]);
     end
  return polyflux, polyft
end

@views function setup_polyflux_single(star_epoch_geom)
  # Polyflux is the weight of each pixel, proportional to the surface
  pjx = star_epoch_geom.projx;
  pjy = star_epoch_geom.projy;
  T = eltype(star_epoch_geom.vertices_xyz);
  polyflux =T(0.5)*(
    pjx[:,1].*pjy[:,2]
  - pjx[:,2].*pjy[:,1]
  + pjx[:,2].*pjy[:,3]
  - pjx[:,3].*pjy[:,2]
  + pjx[:,3].*pjy[:,4]
  - pjx[:,4].*pjy[:,3]
  + pjx[:,4].*pjy[:,1]
  - pjx[:,1].*pjy[:,4]);
  return polyflux;
end

@views function setup_polyft_single(data, star_epoch_geom; T=Float32)
  # Polyflux is the weight of each pixel, proportional to the surface
  pjx = star_epoch_geom.projx;
  pjy = star_epoch_geom.projy;
  kx =  data.uv[1,:] * T(-pi / (180*3600000));
  ky =  data.uv[2,:] * T( pi / (180*3600000));
  # note: check definition sinc(x) = sin(pi*x)/(pi*x)
  polyft = -im*T(1/(2pi))*(((sinc.( (kx*transpose(pjx[:,2]-pjx[:,1])) + (ky*transpose(pjy[:,2]-pjy[:,1])) ).*cis.(-T(pi)*( (kx*transpose(pjx[:,2]+pjx[:,1])) +  (ky*transpose(pjy[:,2]+pjy[:,1])) )).* ( (ky*transpose(pjx[:,2]-pjx[:,1]))  - (kx*transpose(pjy[:,2]-pjy[:,1])) ) + sinc.( (kx*transpose(pjx[:,3]-pjx[:,2])) + (ky*transpose(pjy[:,3]-pjy[:,2])) ).*cis.(-T(pi)*( (kx*transpose(pjx[:,3]+pjx[:,2])) +  (ky*transpose(pjy[:,3]+pjy[:,2])) )).* ( (ky*transpose(pjx[:,3]-pjx[:,2]))  - (kx*transpose(pjy[:,3]-pjy[:,2])) )+ sinc.( (kx*transpose(pjx[:,4]-pjx[:,3])) + (ky*transpose(pjy[:,4]-pjy[:,3])) ).*cis.(-T(pi)*( (kx*transpose(pjx[:,4]+pjx[:,3])) +  (ky*transpose(pjy[:,4]+pjy[:,3])) )).* ( (ky*transpose(pjx[:,4]-pjx[:,3]))  - (kx*transpose(pjy[:,4]-pjy[:,3])) ) + sinc.( (kx*transpose(pjx[:,1]-pjx[:,4])) + (ky*transpose(pjy[:,1]-pjy[:,4])) ).*cis.(-T(pi)*( (kx*transpose(pjx[:,1]+pjx[:,4])) +  (ky*transpose(pjy[:,1]+pjy[:,4])) )).* ( (ky*transpose(pjx[:,1]-pjx[:,4]))  - (kx*transpose(pjy[:,1]-pjy[:,4])) )))./((kx.*kx+ky.*ky)));
  return polyft;
end


function mod360(x)
  return mod.(mod.(x.+180,360).+360, 360) .- 180
end

function cvis_to_v2(cvis, indx)
  v2_model = abs2.(cvis[indx]);
end

function cvis_to_t3(cvis, indx1, indx2, indx3; T=Float32)
  t3 = cvis[indx1].*cvis[indx2].*cvis[indx3];
  t3amp = abs.(t3);
  t3phi = angle.(t3)*T(180/pi);
  return t3, t3amp, t3phi
end

# function cvis_to_t4(cvis, indx1, indx2, indx3, indx4)
#   t4 = cvis[indx1].*cvis[indx2]./(cvis[indx3]*conj(cvis[indx4]));
#   t4amp = abs.(t4);
#   t4phi = angle.(t4)*180.0/pi;
#   return t4, t4amp, t4phi
# end

function observables(x, star, data)
  cvis_model = poly_to_cvis(x, star);
  # compute observables from all cvis
  v2_model = cvis_to_v2(cvis_model, data.indx_v2);
  _, t3amp_model, t3phi_model = cvis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);
  return v2_model, t3amp_model, t3phi_model
end

function chi2s(x, star, data; verbose::Bool = true)
  v2_model, t3amp_model, t3phi_model = observables(x, star, data);
  chi2_v2 = sum(abs2, (v2_model - data.v2)./data.v2_err);
  chi2_t3amp = sum(abs2, (t3amp_model - data.t3amp)./data.t3amp_err);
  chi2_t3phi = sum(abs2, mod360(t3phi_model - data.t3phi)./data.t3phi_err);
  if verbose == true
    println("V2: ", chi2_v2/data.nv2, " T3A: ", chi2_t3amp/data.nt3amp, " T3P: ", chi2_t3phi/data.nt3phi)
  end
  return chi2_v2, chi2_t3amp, chi2_t3phi
end

function chi2s2(x, star, data; verbose::Bool = true)
  v2_model, t3amp_model, t3phi_model = observables(x, star, data);
  chi2_v2 = n2((v2_model - data.v2)./data.v2_err);
  chi2_t3amp = n2((t3amp_model - data.t3amp)./data.t3amp_err);
  chi2_t3phi = n2( mod360(t3phi_model - data.t3phi)./data.t3phi_err);
  if verbose == true
    println("V2: ", chi2_v2/data.nv2, " T3A: ", chi2_t3amp/data.nt3amp, " T3P: ", chi2_t3phi/data.nt3phi)
  end
  return chi2_v2, chi2_t3amp, chi2_t3phi
end

function spheroid_chi2_f(x, star, data; verbose::Bool = false) 
  v2_model, t3amp_model, t3phi_model = observables(x, star, data);
  chi2_v2 = sum(abs2, (v2_model - data.v2)./data.v2_err);
  chi2_t3amp = sum(abs2, (t3amp_model - data.t3amp)./data.t3amp_err);
  chi2_t3phi = sum(abs2, mod360(t3phi_model - data.t3phi)./data.t3phi_err);
  if verbose == true
    println("V2: ", chi2_v2/data.nv2, " T3A: ", chi2_t3amp/data.nt3amp, " T3P: ", chi2_t3phi/data.nt3phi)
  end
  return chi2_v2 + chi2_t3amp + chi2_t3phi
end

# function n2(A) # sum of squares. only makes a difference for thousands of uv points
#   x = zero(eltype(A))
#   @inbounds  @simd for v in A
#     @fastmath x += v * v
#   end
#   return x
# end

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

function spheroid_chi2_allepochs_f(x, stars, data; epochs_weights=[], verbose=false)
nepochs = length(data)
chi2_t = zeros(eltype(x), nepochs);
Threads.@threads for i=1:nepochs # weighted sum -- should probably do the computation in parallel
  chi2_t[i] = spheroid_chi2_f(x, stars[i], data[i], verbose=verbose);
end
f = sum(chi2_t)
if epochs_weights!=[]
  f = f.*epochs_weights
end
if verbose == true 
  println("All epochs, weighted chi2: ", f, "\n");
end
return f;
end

function parametric_temperature_map(parameters, star) # dispatches parametric 
  if star.surface_type == 3
    return temperature_map_vonZeipel_roche_single(parameters,star, star.t);
  elseif star.surface_type == 2
    return temperature_map_vonZeipel_rapid_rotator(parameters,star);
  else
    println("Unimplemented parametric von Zeipel function")
  end
  return
end

function spheroid_parametric_f(parameters, tessels, data, tepochs) 
  stars = create_star_multiepochs(tessels, parameters, tepochs);
  setup_oi!(data, stars) 
  x = parametric_temperature_map(parameters, stars[1]);
  return spheroid_chi2_allepochs_f(x, stars, data)
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
