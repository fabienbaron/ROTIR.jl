using LinearAlgebra, OptimPackNextGen
# The Fourier transform for polygonal surfaces
function poly_to_cvis(x, polyflux, polyft)
  flux = sum(polyflux.*x); # get the total flux
  cvis_model = polyft * x / flux;
end

function poly_to_cvis(x, star)
   indx = star.index_quads_visible
   xw = x[indx] .* star.vis_weights[indx]
   flux = dot(star.polyflux, xw); # get the total flux
   return star.polyft * xw / flux;
end

function poly_to_flux(x, star)
   indx = star.index_quads_visible
   xw = x[indx] .* star.vis_weights[indx]
   flux = dot(star.polyflux, xw); # get the total flux
   return flux;
end

@views function setup_oi!(data, stars)
  nepochs = size(data,1);
  T = eltype(stars[1].vertices_xyz);
  if nepochs>1
  Threads.@threads for i=1:nepochs
       indx = stars[i].index_quads_visible
       stars[i].polyflux = setup_polyflux_single(stars[i].projx[indx,:], stars[i].projy[indx,:])
       stars[i].polyft = setup_polyft_single(data[i].uv, stars[i].projx[indx,:], stars[i].projy[indx,:]);
     end
  else # single epoch, thread over calculation
    indx = stars[1].index_quads_visible
    stars[1].polyflux = setup_polyflux_single(stars[1].projx[indx,:], stars[1].projy[indx,:])
    stars[1].polyft = setup_polyft_single_alt(data[1].uv, stars[1].projx[indx,:], stars[1].projy[indx,:]);
  end
end

@views function setup_polygon_ft(data, star_epoch_geom)
  nepochs = size(data,1);
  T = eltype(star_epoch_geom[1].vertices_xyz);
  polyflux = Array{Array{T,1}}(undef,nepochs);
  polyft = Array{Array{Complex{T},2}}(undef,nepochs);
  Threads.@threads for i=1:nepochs
       indx = star_epoch_geom[i].index_quads_visible
       polyflux[i] = setup_polyflux_single(star_epoch_geom[i].projx[indx,:], star_epoch_geom[i].projy[indx,:])
       polyft[i] = setup_polyft_single(data[i].uv, star_epoch_geom[i].projx[indx,:], star_epoch_geom[i].projy[indx,:]);
     end
  return polyflux, polyft
end

@views function setup_polyflux_single(pjx, pjy)
  # Polyflux is the projected area of each pixel (shoelace formula)
  T = eltype(pjx);
  polyflux = T(0.5)*(
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

function stcis(x1,x2,y1,y2,kx,ky)
  return sinc.(kx*(x2-x1) + ky*(y2-y1)).*cis.(-π*(kx*(x2+x1)+ ky*(y2+y1))).*(ky*(x2-x1)-kx*(y2-y1))
end

@views function setup_polyft_single_alt(uv, pjx, pjy; T=Float64)
  kx =  uv[1,:] * T(-pi / (180*3600000))
  ky =  uv[2,:] * T( pi / (180*3600000))
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

@views function setup_polyft_single(uv, pjx, pjy; T=Float64)
  kx =  uv[1,:] * T(-pi / (180*3600000));
  ky =  uv[2,:] * T( pi / (180*3600000));
  # note: check definition sinc(x) = sin(pi*x)/(pi*x)
  polyft = -im*T(1/(2pi))*(((sinc.( (kx*transpose(pjx[:,2]-pjx[:,1])) + (ky*transpose(pjy[:,2]-pjy[:,1])) ).*cis.(-T(pi)*( (kx*transpose(pjx[:,2]+pjx[:,1])) +  (ky*transpose(pjy[:,2]+pjy[:,1])) )).* ( (ky*transpose(pjx[:,2]-pjx[:,1]))  - (kx*transpose(pjy[:,2]-pjy[:,1])) ) + sinc.( (kx*transpose(pjx[:,3]-pjx[:,2])) + (ky*transpose(pjy[:,3]-pjy[:,2])) ).*cis.(-T(pi)*( (kx*transpose(pjx[:,3]+pjx[:,2])) +  (ky*transpose(pjy[:,3]+pjy[:,2])) )).* ( (ky*transpose(pjx[:,3]-pjx[:,2]))  - (kx*transpose(pjy[:,3]-pjy[:,2])) )+ sinc.( (kx*transpose(pjx[:,4]-pjx[:,3])) + (ky*transpose(pjy[:,4]-pjy[:,3])) ).*cis.(-T(pi)*( (kx*transpose(pjx[:,4]+pjx[:,3])) +  (ky*transpose(pjy[:,4]+pjy[:,3])) )).* ( (ky*transpose(pjx[:,4]-pjx[:,3]))  - (kx*transpose(pjy[:,4]-pjy[:,3])) ) + sinc.( (kx*transpose(pjx[:,1]-pjx[:,4])) + (ky*transpose(pjy[:,1]-pjy[:,4])) ).*cis.(-T(pi)*( (kx*transpose(pjx[:,1]+pjx[:,4])) +  (ky*transpose(pjy[:,1]+pjy[:,4])) )).* ( (ky*transpose(pjx[:,1]-pjx[:,4]))  - (kx*transpose(pjy[:,1]-pjy[:,4])) )))./((kx.*kx+ky.*ky)));
  return polyft;
end




@views function setup_polyft_single(kx, ky, pjx, pjy; T=Float64)
  polyft = -im*T(1/(2pi))*(((sinc.( (kx*transpose(pjx[:,2]-pjx[:,1])) + (ky*transpose(pjy[:,2]-pjy[:,1])) ).*cis.(-T(pi)*( (kx*transpose(pjx[:,2]+pjx[:,1])) +  (ky*transpose(pjy[:,2]+pjy[:,1])) )).* ( (ky*transpose(pjx[:,2]-pjx[:,1]))  - (kx*transpose(pjy[:,2]-pjy[:,1])) ) + sinc.( (kx*transpose(pjx[:,3]-pjx[:,2])) + (ky*transpose(pjy[:,3]-pjy[:,2])) ).*cis.(-T(pi)*( (kx*transpose(pjx[:,3]+pjx[:,2])) +  (ky*transpose(pjy[:,3]+pjy[:,2])) )).* ( (ky*transpose(pjx[:,3]-pjx[:,2]))  - (kx*transpose(pjy[:,3]-pjy[:,2])) )+ sinc.( (kx*transpose(pjx[:,4]-pjx[:,3])) + (ky*transpose(pjy[:,4]-pjy[:,3])) ).*cis.(-T(pi)*( (kx*transpose(pjx[:,4]+pjx[:,3])) +  (ky*transpose(pjy[:,4]+pjy[:,3])) )).* ( (ky*transpose(pjx[:,4]-pjx[:,3]))  - (kx*transpose(pjy[:,4]-pjy[:,3])) ) + sinc.( (kx*transpose(pjx[:,1]-pjx[:,4])) + (ky*transpose(pjy[:,1]-pjy[:,4])) ).*cis.(-T(pi)*( (kx*transpose(pjx[:,1]+pjx[:,4])) +  (ky*transpose(pjy[:,1]+pjy[:,4])) )).* ( (ky*transpose(pjx[:,1]-pjx[:,4]))  - (kx*transpose(pjy[:,1]-pjy[:,4])) )))./((kx.*kx+ky.*ky)));
  return polyft;
end





function mod360(x)
  return mod.(mod.(x.+180,360).+360, 360) .- 180
end

function cvis_to_v2(cvis, indx)
  v2_model = abs2.(cvis[indx]);
end

function cvis_to_t3(cvis, indx1, indx2, indx3; T=Float64)
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

@views function spheroid_chi2_fg(x, g, star, data; verbose::Bool = true)
  npix = star.npix;
  T = eltype(x);
  indx = star.index_quads_visible
  w = star.vis_weights[indx]  # soft visibility weights for visible pixels
  xw = x[indx] .* w           # weighted pixel values
  cvis_model = poly_to_cvis(x, star);
  v2_model = cvis_to_v2(cvis_model, data.indx_v2);
  t3_model, t3amp_model, t3phi_model = cvis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);
  chi2_v2 = sum(abs2, (v2_model - data.v2)./data.v2_err);
  chi2_t3amp = sum(abs2, (t3amp_model - data.t3amp)./data.t3amp_err);
  chi2_t3phi = sum(abs2, mod360(t3phi_model - data.t3phi)./data.t3phi_err);
  # Gradient w.r.t. weighted pixels (xw), computed via polyft^T @ adjoint_signal
  g_v2 = real(transpose(star.polyft[data.indx_v2,:])*(4*((v2_model-data.v2)./data.v2_err.^2).*conj(cvis_model[data.indx_v2])));
  g_t3amp = real(transpose(star.polyft[data.indx_t3_1,:]) *(2*((t3amp_model-data.t3amp)./data.t3amp_err.^2).*conj(cvis_model[data.indx_t3_1])./abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_2]).*abs.(cvis_model[data.indx_t3_3]) ))+real(transpose(star.polyft[data.indx_t3_2,:])*(2*((t3amp_model-data.t3amp)./data.t3amp_err.^2).*conj(cvis_model[data.indx_t3_2])./abs.(cvis_model[data.indx_t3_2]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_3]) ))+real(transpose(star.polyft[data.indx_t3_3,:])*(2*((t3amp_model-data.t3amp)./data.t3amp_err.^2).*conj(cvis_model[data.indx_t3_3])./abs.(cvis_model[data.indx_t3_3]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_2]) ));
  g_t3phi = T(360/pi)*imag(transpose(star.polyft[data.indx_t3_1,:])*(((mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2)./abs2.(t3_model)).*cvis_model[data.indx_t3_2].*cvis_model[data.indx_t3_3].*conj(t3_model))+transpose(star.polyft[data.indx_t3_2,:])*(((mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2)./abs2.(t3_model)).*cvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_3].*conj(t3_model))+transpose(star.polyft[data.indx_t3_3,:])*(((mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2)./abs2.(t3_model)).*cvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_2].*conj(t3_model)));
  gsum = g_v2 + g_t3amp + g_t3phi;
  flux = poly_to_flux(x, star);
  # Gradient w.r.t. x: chain rule through soft visibility weights
  # gsum is ∂χ²/∂(xw), so ∂χ²/∂x = w .* ∂χ²/∂(xw) (after flux normalization)
  g_normalized = (gsum .- dot(xw, gsum) * star.polyflux / flux) / flux;
  g[indx] = w .* g_normalized;
  if verbose == true
    println("V2: ", chi2_v2/data.nv2, "\tT3A: ", chi2_t3amp/data.nt3amp, "\tT3P: ", chi2_t3phi/data.nt3phi,"\tFlux: ", flux)
  end
  return chi2_v2 + chi2_t3amp + chi2_t3phi
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

function spheroid_crit_allepochs_fg(x, g, stars, data; regularizers=[], epochs_weights=[], verbose=false, T=Float64)
  T = eltype(x)
  #g[:] .= T(0);
  nepochs = length(data)
  chi2_t = zeros(T, nepochs);
  #npix = stars[1].npix
  singleepoch_g = [zeros(T, length(x)) for i=1:nepochs];
  Threads.@threads for i=1:nepochs # weighted sum -- should probably do the computation in parallel
    chi2_t[i] = spheroid_chi2_fg(x, singleepoch_g[i], stars[i], data[i], verbose=verbose);
  end
  f = sum(chi2_t)
  if epochs_weights!=[]
  #  f = f.*epochs_weights
    @warn "Epoch weights not implemented"
   end 
  g[:] .= sum(singleepoch_g)
  if verbose == true
    println("Total weighted chi2: $f");
  end
  # Map regularization
  if regularizers!=[]
    reg_g = zeros(T, length(x));
    f += spheroid_regularization(x, reg_g, regularizers=regularizers, verbose = verbose);
    g[:] += reg_g
  end
  return f;
end

function parametric_temperature_map(parameters, star) # dispatches parametric 
  if star.surface_type == 3
    return temperature_map_vonZeipel_roche_single(parameters,star, star.t);
  elseif star.surface_type == 2
    return temperature_map_vonZeipel_rapid_rotator(parameters,star);
  elseif star.surface_type == 1
    return temperature_map_vonZeipel_ellipsoid(parameters,star);
  elseif star.surface_type == 0 # sphere
    return parameters.tpole*ones(eltype(star.vertices_spherical), star.npix );  
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


# function spheroid_total_variation_fg(x, tv_g, tvinfo; ϵ = 1e-13, T=Float32, verbose = true)
#   # Add total variation regularization
#   #tvinfo: neighbors,south_neighbors,west_neighbors,south_neighbors_reverse,west_neighbors_reverse
#   ϵ = T(ϵ)
#   npix = length(x)
#   xs = x[tvinfo[2]];
#   xw = x[tvinfo[3]];
#   tv_f = sum(sqrt.( (x-xs).^2 + (x-xw).^2) .+ ϵ )
#   tv_g[1:npix] = (2*x-xs-xw)./(sqrt.( (x-xs).^2 + (x-xw).^2) .+ ϵ)
#   for j=1:length(x)
#     k = tvinfo[4][j]
#     l = tvinfo[5][j]
#     if length(k)>0
#       kk=k[1]
#       tv_g[j] += (x[j]-x[kk])/sqrt((x[kk]-x[j])^2+(x[kk]-x[tvinfo[3][kk]])^2 +ϵ)
#     end
#     if length(l)>0
#       ll=l[1]
#       tv_g[j] += (x[j]-x[ll])/sqrt((x[ll]-x[tvinfo[2][ll]])^2+(x[ll]-x[j])^2 +ϵ)
#     end
#    end
# if verbose == true
#       println("TV: ", tv_f);
# end
#   return tv_f
# end


function spheroid_total_variation2_fg(x, tv_g, tvinfo; verbose = true)
  npix = length(x)
  tv_f = norm(tvinfo[6]*x)^2
  tv_g[:] = 2*(tvinfo[7]*x)
  if verbose == true
      println("TV2: ", tv_f);
  end
  return tv_f
end

function spheroid_total_variation_fg(x, tv_g, tvinfo; verbose = true)
  npix = length(x)
  tv_f = norm(tvinfo[6]*x)
  if tv_f>0 
  tv_g[:] = (tvinfo[7]*x)/tv_f
  else
    tv_g[:] .= 0
  end
  if verbose == true
      println("TV: ", tv_f);
  end
  return tv_f
end

function spheroid_mean_fg(x, g; verbose = true)
f = sum(abs.(x.-sum(x)/length(x)))
g[:] = sign.(x.-sum(x)/length(x))
if verbose == true
println(" MeanReg: ", f);
end
return f;
end

function spheroid_harmon_bias_fg(x, g, B::Float64; verbose = true)
n = length(x);
avg_x = mean(x);
bcorr = (B-1.0)*Int.((x.-avg_x).>0).+1.0
#reg_f = sum(bcorr.*(x.-avg_x).^2);
reg_f = sum(bcorr.*(x.-avg_x).^2)/n;
#reg_g = 2*(n-1)/n*bcorr.*(x.-avg_x)
reg_g = 2/n*bcorr.*(x.-avg_x)
g[:] = reg_g .- mean(reg_g); # necessary ?

if verbose == true
println(" Bias: ", reg_f);
end
return reg_f;
end

function max_entropy_fg(x, g; verbose=false, ϵ=1e-9)
  # mmap = sum(x) / length(x)
  # nmap = x ./ mmap
  # reg_f = sum(nmap .* log.(nmap))
  # g[:] = sum(-nmap ./ sum(x) .* (log.(nmap) .+ 1.0)) .+ (log.(nmap) .+ 1)./mmap
  xm = x ./ (mean(x) + ϵ)
  reg_f = sum(xm .* log.(xm .+ ϵ))
  g[:] = ((mean(x) .- (x ./ length(x))) ./ (mean(x)^2 + ϵ)) .* (log.(xm .+ ϵ) .+ 1)
  if verbose == true
    println(" MEM: ", reg_f)
  end
  return reg_f
end

function spheroid_regularization(x,g; printcolor = :black, regularizers=[], verbose=false)
  reg_f = 0.0;
  for ireg =1:length(regularizers)
      x_sub = x[regularizers[ireg][4]] # take the pixel subset if needed (example: only regularize visible pixels)
      temp_g = similar(x_sub);
      if regularizers[ireg][1] == "mem"
          reg_f += regularizers[ireg][2]*max_entropy_fg(x_sub, temp_g, verbose = verbose);
      elseif regularizers[ireg][1] == "tv2"
          reg_f += regularizers[ireg][2]*spheroid_total_variation2_fg(x_sub, temp_g, regularizers[ireg][3], verbose = verbose);
      elseif regularizers[ireg][1] == "tv"
          reg_f += regularizers[ireg][2]*spheroid_total_variation_fg(x_sub, temp_g, regularizers[ireg][3], verbose = verbose);
      elseif regularizers[ireg][1] == "mean"
          reg_f += regularizers[ireg][2]*spheroid_mean_fg(x_sub, temp_g, verbose = verbose);
      elseif regularizers[ireg][1] == "bias"
          reg_f += regularizers[ireg][2]*spheroid_harmon_bias_fg(x_sub, temp_g, regularizers[ireg][3], verbose = verbose );
      end
      g[regularizers[ireg][4]] += regularizers[ireg][2]*temp_g
  end
  if verbose == true
      println("\n");
  end
return reg_f
end


using OptimPackNextGen

function image_reconstruct_oi(x_start, data, stars; epochs_weights =[], printcolor= [], verbose = true, lower=0, upper=Inf, maxiter = 100, regularizers =[])
  x_sol = [];
  crit_imaging = (x,g)->spheroid_crit_allepochs_fg(x, g, stars, data, regularizers=regularizers, epochs_weights=  epochs_weights, verbose = verbose);
  x_sol = OptimPackNextGen.vmlmb(crit_imaging, x_start, verb=verbose, lower=lower, upper=upper, maxiter=maxiter, blmvm=false, gtol=(0,1e-8));
  dummy = similar(x_sol);
  crit_opt = crit_imaging(x_sol,dummy);
  return x_sol
end

function image_reconstruct_oi_crit(x, data, stars; regularizers =[],  verbose = verbose)
  g = similar(x);
  crit = spheroid_crit_allepochs_fg(x, g, stars, data, regularizers=regularizers,epochs_weights=[],  verbose = verbose);
  return crit
end

function image_reconstruct_oi_chi2(x, data, stars;  verbose = verbose)
  g = similar(x);
  crit = spheroid_crit_allepochs_fg(x, g, stars, data, regularizers=[], epochs_weights= [], verbose = verbose);
  return crit
end

function image_reconstruct_oi_chi2_fg(x, data, stars;  verbose = verbose)
  g = similar(x);
  crit = spheroid_crit_allepochs_fg(x, g, stars, data, regularizers=[], epochs_weights= [], verbose = verbose);
  return crit,g
end

function multires_reconstruct_oi(data, star_params, tepochs; n_start=2, n_end=4, maxiter=500, reg_weight=1e-5, reg_type="tv2", verbose=true, kwargs...)
  tmap = nothing
  stars = nothing
  for n in n_start:n_end
    tessels = tessellation_healpix(n)
    stars = create_star_multiepochs(tessels, star_params, tepochs)
    if tmap === nothing
      tmap = parametric_temperature_map(star_params, stars[1])
    else
      tmap = vec(repeat(tmap, 1, 4)')  # upsample from previous level
    end
    setup_oi!(data, stars)
    regularizers = [[reg_type, reg_weight, tv_neighbours_healpix(n), 1:length(tmap)]]
    if verbose
      println("Multi-resolution: HEALPix level n=$n, npix=$(nside2npix(2^n))")
    end
    tmap = image_reconstruct_oi(tmap, data, stars; maxiter=maxiter, regularizers=regularizers, verbose=verbose, kwargs...)
  end
  return tmap, stars
end





# "Multitemporal" = dynamical reconstruction

# function oi_reconstruct_mutitemporal(x_start::Array{Float64,1}, data::Array{OIdata,1}, polyflux::Array{Array{Float64,1},1}, polyft::Array{Array{Complex{Float64},2},1}; epochs_weights =[], printcolor= [], verbose = true, maxiter = 100, regularizers =[])
#   crit_imaging = (x,g)->oi_multitemporal_fg(x, g, polyflux, polyft, data, regularizers=regularizers, epochs_weights=  epochs_weights);
#   x_sol = OptimPackNextGen.vmlmb(crit_imaging, x_start, verbose=verbose, lower=0, maxiter=maxiter, blmvm=false, gtol=(0,1e-8));
#   return reshape(x_sol, size(polyflux[1],1), length(data))
# end

# function oi_multitemporal_fg(x, g, polyflux, polyft, data; regularizers=[], epochs_weights = [], verbose=true)
#   # Explanation of the following: optimpack optimizes a vector, but we want an array of images
#   npix = size(polyflux[1],1);
#   chi2_f = 0.0;
#   g[:] .= 0.0;
#   #temp_g = deepcopy(g[:]);
#   #npix = size(x);
#   nepochs = length(data);
#   if epochs_weights == []
#     epochs_weights = ones(Float64, nepochs)/nepochs;
#   end

#   #singleepoch_g = zeros(Float64, npix);
#   #for i=1:nepochs # weighted sum -- should probably do the computation in parallel
#   #  chi2_f += epochs_weights[i]*spheroid_chi2_fg_alt(x[1:npix], singleepoch_g, polyflux[i], polyft[i], data[i], verbose=verbose);
#   #  g[1:npix] += epochs_weights[i]*singleepoch_g;
#   #end


#   for i=1:nepochs # weighted sum -- in the future, do the computation in parallel
#     tslice = 1+(i-1)*npix:i*npix; # temporal slice
#     subg = zeros(Float64, npix);
#     chi2_f += epochs_weights[i]*spheroid_chi2_fg_alt(x[tslice], subg, polyflux[i], polyft[i], data[i], verbose=verbose);
#     g[tslice] = epochs_weights[i]*subg;
#     #x[tslice] = x[1:npix];
#     #g[tslice] = g[1:npix];
#   end

#   # cross temporal regularization -- weight needs to be defined in the "regularizers" variable
#   if length(regularizers)>nepochs
#     if (regularizers[nepochs+1][1][1] == "temporal_tvsq")  & (regularizers[nepochs+1][1][2] > 0.0) & (nepochs>1)
#       y = reshape(x,(npix,nepochs))
#       temporalf = sum( (y[:,2:end]-y[:,1:end-1]).^2 )
#       tv_g = Array{Float64}(undef, npix,nepochs)
#       if nepochs>2
#          tv_g[:,1] = 2*(y[:,1] - y[:,2])
#          tv_g[:,2:end-1] = 4*y[:,2:end-1]-2*(y[:,1:end-2]+y[:,3:end])
#          tv_g[:,end] = 2*(y[:,end] - y[:,end-1])
#       else
#          tv_g[:,1] = 2*(y[:,1]-y[:,2]);
#          tv_g[:,2] = 2*(y[:,2]-y[:,1]);
#       end
#       chi2_f += regularizers[nepochs+1][1][2]*temporalf
#       g[:] += regularizers[nepochs+1][1][2]*vec(tv_g);
#       if verbose == true
#            printstyled("Temporal regularization: $temporalf\n", color=:yellow)
#       end
#     end
#   end
#   reg_f = spheroid_regularization(x, g, regularizers=regularizers[1:nepochs], verbose = verbose); # Note: adds to g
#   return chi2_f + reg_f;
# end
