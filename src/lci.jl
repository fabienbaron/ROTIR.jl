using DelimitedFiles
using LinearAlgebra
using Printf

mutable struct LCI
  mjd::Array{Float64,1}
  phase::Array{Float64,1}
  flux::Array{Float64,1}
  fluxerr::Array{Float64,1}
  nepochs::Int64
end

function setup_stellar_parameters_lci(starparams, lcidata)
stellar_parameters = Array{starparameters}(undef, lcidata.nepochs);
for i=1:lcidata.nepochs
    stellar_parameters[i]=starparameters(starparams[1],starparams[2],starparams[3],starparams[4],starparams[5],starparams[6],starparams[7],starparams[8],mod(360.0/starparams[10].*lcidata.mjd[i],360),starparams[10]);
end
return stellar_parameters
end


function split_lcidata_by_period(lcidata,polyflux, P)
    mjd_intervals = [range(minimum(lcidata.mjd), maximum(lcidata.mjd), step=P);Inf];
    nframes = length(mjd_intervals)-1;
    Y = Array{Array{Float64,1}}(undef,nframes);
    E = Array{Array{Float64,1}}(undef,nframes);
    W = Array{Diagonal{Float64,Array{Float64,1}}}(undef,nframes);
    H = Array{Array{Float64,2}}(undef,nframes);
    for i=1:nframes
        epoch_filter = findall((lcidata.mjd.>=mjd_intervals[i]) .& (lcidata.mjd.<mjd_intervals[i+1]))
        Y[i] = lcidata.flux[epoch_filter];
        E[i] = lcidata.fluxerr[epoch_filter];
        W[i]=Diagonal(1 ./lcidata.fluxerr[epoch_filter].^2);
        H[i] = polyflux[epoch_filter,:]
    end
    return nframes, Y, E, W, H
end

function setup_regularization_matrices(star_epoch_geom)
visible_pixels = sometimes_visible(star_epoch_geom); # lists all pixels at least visible once
# Setup regularization matrices
npix = star_epoch_geom[1].npix;
A=zeros(npix)';
A[visible_pixels].=1.0/length(visible_pixels); # A*x    <=> mean(x[visible_pixels])
C = (Diagonal(ones(npix)).-A); ;# C*x    <=> x.-mean(x[visible_pixels])
# Total variation
neighbors,south_neighbors,west_neighbors,south_neighbors_reverse,west_neighbors_reverse = tv_neighbours_healpix(npix2n(npix));
∇s = sparse(1:npix, 1:npix, 1.0) + sparse(1:npix, south_neighbors, -1.0);
∇w = sparse(1:npix, 1:npix, 1.0) + sparse(1:npix, west_neighbors, -1.0);
return C, ∇s, ∇w
end

function read_lci_relative(file)
lc = readdlm(file, skipstart=6)
phase = vec(lc[:,1]);
flux = vec(lc[:,2]);
fluxerr = ones(size(flux))
return LCI([],phase, flux, fluxerr, length(phase))
end

function read_lci_absolute(file; headerlength = 6, use_mjd = true, use_phase = false, normalized = false) # if mjd = false, column 1 is the phase
if use_phase == true
  use_mjd = false
end
lc = readdlm(file, skipstart=headerlength)
mjd = [];
phase = [];
if use_mjd == true
  mjd = vec(lc[:,1]);
else
  phase = vec(lc[:,1]);
end
flux = vec(lc[:,2]);
fluxerr = vec(lc[:,3]);
return LCI(mjd,phase, flux, fluxerr, length(flux))
end

function write_lci(file, lcidata)
open(file, "w") do f
           write(f, "; Fake Kepler light curve\n");
           write(f, "; No Quarter ! Yarr\n");
           write(f, "; CBVs have been removed.\n");
           write(f, "; Column 1:  Time (MBJD)\n");
           write(f, "; Column 2:  Normalized flux with CBVs removed\n");
           write(f, "; Column 3:  Errors on the flux data\n");
           for i in 1:length(lcidata.phase)
              write(f, "$(lcidata.phase[i]) $(lcidata.flux[i]) $(lcidata.fluxerr[i])\n")
           end
       end
close(file)
end


function modelflux_lci(x, polyflux)
  flux = polyflux*x;
  return flux;
end

function modelflux_lci_rel(x, polyflux)
  flux = polyflux*x;
  flux /= mean(flux);
  return flux;
end

function setup_lci(star_epoch_geom; ld = true)
# Polyflux is the weight of each Healpixel, proportional to the surface
nepochs = length(star_epoch_geom);
polyflux = zeros(Float64, nepochs, star_epoch_geom[1].npix);
for i=1:nepochs
polyflux[i,star_epoch_geom[i].index_quads_visible] =0.5*(star_epoch_geom[i].projx[:,1].*star_epoch_geom[i].projy[:,2]
- star_epoch_geom[i].projx[:,2].*star_epoch_geom[i].projy[:,1]
+ star_epoch_geom[i].projx[:,2].*star_epoch_geom[i].projy[:,3]
- star_epoch_geom[i].projx[:,3].*star_epoch_geom[i].projy[:,2]
+ star_epoch_geom[i].projx[:,3].*star_epoch_geom[i].projy[:,4]
- star_epoch_geom[i].projx[:,4].*star_epoch_geom[i].projy[:,3]
+ star_epoch_geom[i].projx[:,4].*star_epoch_geom[i].projy[:,1]
- star_epoch_geom[i].projx[:,1].*star_epoch_geom[i].projy[:,4]);
if ld == true
  polyflux[i,:] .*= star_epoch_geom[i].ldmap;
end
end
visible_pixels = sometimes_visible(star_epoch_geom); # lists all pixels at least visible once
hidden_pixels =  never_visible(star_epoch_geom);
return polyflux, visible_pixels, hidden_pixels
end


function chi2_lci_f(x, polyflux, lcidata)
flux = polyflux*x;
res = (flux-lcidata.flux)./lcidata.fluxerr;
chi2_f = norm(res,2)^2 ;  # norm((modelflux_lci(x,polyflux)-lcidata.flux)./lcidata.fluxerr,2)^2
return chi2_f
end

function lci_absolute_crit_fg(x::Array{Float64,1}, g::Array{Float64,1},  dataflux::Array{Float64,1}, datafluxerr::Array{Float64,1}, polyflux::Array{Float64,2}, visible_pixels::Array{Int64,1}; regularizers=[], verb = true )
flux = polyflux*x;
res = (flux-dataflux)./datafluxerr;
chi2_f = norm(res,2)^2 ;  # norm((modelflux_lci(x,polyflux)-lcidata.flux)./lcidata.fluxerr,2)^2
g[:] = 2.0*polyflux'*(res./datafluxerr);
reg_f = spheroid_regularization(x, g, regularizers=regularizers, verb = false); #note: this adds to g
if verb == true
  y= x[visible_pixels]
  low_x = minimum(y);
  mean_x = mean(y)
  ratio = abs((mean_x-low_x)/mean_x)
  @printf("Crit: %8.1f chi2r: %8.2f λ*reg:%8.2f minT:%8.1f meanT:%8.1f ratio:%8.2f\n", chi2_f+reg_f, chi2_f/length(flux), reg_f, low_x, mean_x, ratio);
end
return chi2_f + reg_f
end


function lci_relative_crit_fg(x, g, dataflux, datafluxerr, polyflux, visible_pixels; regularizers=[] , verb = true)
flux = polyflux*x;
flux_avg = mean(flux);
res = (flux/flux_avg-dataflux)./datafluxerr;
chi2_f = norm(res,2)^2 ;  # norm((modelflux_lci(x,polyflux)-lcidata.flux)./lcidata.fluxerr,2)^2
g[:] = 2.0/flux_avg*polyflux'*(res./datafluxerr)-2.0/flux_avg^2*sum(flux.*(res./datafluxerr))*vec(mean(polyflux,dims=1));
reg_f = spheroid_regularization(x, g, regularizers=regularizers, verb = false); #note: this adds to g
if verb == true
    v = x[visible_pixels]
    low_x = minimum(v);
    hi_x = maximum(v);
    Δx = hi_x - low_x
    println("Crit: $(chi2_f) \t chi2r: $(chi2_f/length(flux)) \t λ*reg: $(reg_f) \t minT: $(low_x) \t maxT: $(hi_x) \t ΔT: $(Δx)");
end
return chi2_f + reg_f
end



function lci_multitemporal_absolute_fg(x::Array{Float64,1}, g::Array{Float64,1}, F, E, H, visible_pixels; printcolor= [], verb = true, regularizers =[])
  # Explanation of the following: optimpack optimizes a vector, but we want an array of images
  # So we need to unpack/repack it
  # H: polyflux
  # F: flux data, E: error on flux data
  npix = size(H[1],2)
  nframes = length(F)
  f = 0.0;
  for i=1:nframes # weighted sum -- in the future, do the computation in parallel
    tslice = 1+(i-1)*npix:i*npix; # temporal slice
    subg = Array{Float64}(undef, npix);
     if verb == true
     printstyled("Frame $i ",color=:black);
    end
    f+= lci_absolute_crit_fg(x[tslice], subg,  F[i], E[i], H[i], visible_pixels, regularizers=regularizers[i], verb = true )
    g[tslice] = subg
  end

  # cross temporal regularization -- weight needs to be defined in the "regularizers" variable
   if length(regularizers)>nframes
     if (regularizers[nframes+1][1][1] == "temporal_tvsq")  & (regularizers[nframes+1][1][2] > 0.0) & (nframes>1)
      y = reshape(x,(npix,nframes))
      temporalf = sum( (y[:,2:end]-y[:,1:end-1]).^2 )
      tv_g = Array{Float64}(undef, npix,nframes)
      if nframes>2
         tv_g[:,1] = 2*(y[:,1] - y[:,2])
         tv_g[:,2:end-1] = 4*y[:,2:end-1]-2*(y[:,1:end-2]+y[:,3:end])
         tv_g[:,end] = 2*(y[:,end] - y[:,end-1])
      else
         tv_g[:,1] = 2*(y[:,1]-y[:,2]);
         tv_g[:,2] = 2*(y[:,2]-y[:,1]);
      end
      f+= regularizers[nframes+1][1][2]*temporalf
      g[:] += regularizers[nframes+1][1][2]*vec(tv_g);
      if verb == true
           printstyled("Temporal regularization: $temporalf\n", color=:yellow)
      end
     end
    end
   return f;
end


function lci_reconstruct(x_start::Array{Float64,1}, lcidata::LCI, polyflux::Array{Float64,2}, visible_pixels; lowtemp = 0, relative = false, printcolor= [], verb = true, maxiter = 100, regularizers =[])
  x_sol = [];
  crit_imaging=x->x; #TODO: proper dummy init
  if relative == true
  crit_imaging = (x,g)->lci_relative_crit_fg(x, g, lcidata.flux, lcidata.fluxerr, polyflux, visible_pixels, regularizers=regularizers, verb = verb);
  else
    crit_imaging = (x,g)->lci_absolute_crit_fg(x, g, lcidata.flux, lcidata.fluxerr, polyflux, visible_pixels, regularizers=regularizers, verb = verb);
  end
  x_sol = OptimPackNextGen.vmlmb(crit_imaging, x_start, verb=false, lower=lowtemp, maxiter=maxiter, blmvm=false, gtol=(0,1e-8));
  return x_sol
end

function lci_reconstruct_mutitemporal(x_start::Array{Float64,1}, F, E, H, visible_pixels; lowtemp = 0, relative = false, printcolor= [], verb = true, maxiter = 100, regularizers =[])
  crit_imaging = (x,g)->lci_multitemporal_absolute_fg(x, g, F, E, H, visible_pixels, verb = true, regularizers = regularizers);
  x_sol = OptimPackNextGen.vmlmb(crit_imaging, x_start, verb=false, lower=lowtemp, maxiter=maxiter, blmvm=false, gtol=(0,1e-8));
  return reshape(x_sol, (size(H[1],2), length(F)))
end

function lci_linear_inversion_frame(xstart, H, W, Y, C, ∇s, ∇w, B, λ , μ; tol = 1e-4, maxiter = 50, verbose = true)
ϵ = 1e99
iter = 1
x = deepcopy(xstart)
xprev = deepcopy(xstart);
npix = size(H,2)
# Note: equation is:  x = (H'*W*H+λ*C'*(V*C)+μ*(∇s'∇s + ∇w'∇w))\(H'*W*F);
left = H'*W*H+μ*(∇s'∇s + ∇w'∇w)
right = H'*W*Y
while (ϵ>tol)&(iter<=maxiter)
    #V = Diagonal(((B-1.0)*Int.(C*x.>0).+1.0)/nx)
    V = spdiagm(0=>((B-1.0)*Int.(C*x.>0).+1.0)/npix)
    x = (left+λ*C'*(V*C))\right;
    x .*= (x.>0)
    if iter>1
        ϵ = norm(x - xprev)
    end
    if verbose == true
        println("Inner iteration: $iter $ϵ")
    end
    xprev = deepcopy(x);
    iter +=1;
end
return x;
end

function rescale_temperature(x, Tphot_wanted, visible_pixels)
nframes = size(x,2)
for i=1:nframes
    x[:,i] .*= Tphot_wanted/mean(x[visible_pixels,i])
end
return x
end


#
# function chi2_lci_bias_fg(x, g, polyflux, lcidata, B, λ, somvis; verbose = true)
# flux = polyflux*x;
# flux_avg = mean(flux);
# res = (flux/flux_avg-lcidata.flux)./lcidata.fluxerr;
# chi2_f = norm(res,2)^2 ;  # norm((modelflux_lci(x,polyflux)-lcidata.flux)./lcidata.fluxerr,2)^2
# chi2_g = 2.0/flux_avg*polyflux'*(res./lcidata.fluxerr)-2.0/flux_avg^2*sum(flux.*(res./lcidata.fluxerr))*vec(mean(polyflux,dims=1));
# # Harmon regularization
# xvis = x[somvis];
# avg_xvis = mean(xvis);
# n = length(xvis);
# bcorr = (B-1.0)*Int.(xvis.>avg_xvis)+1.0
# reg_f = sum(bcorr.*(xvis-avg_xvis).^2)/n;
# reg_g = 2*(xvis-avg_xvis).*bcorr/n;
# Δx = maximum(x)-minimum(x)
# if verbose == true
#     low_x = minimum(x[somvis]);
#     hi_x = maximum(x[somvis]);
#     println("Crit: $(chi2_f + 2*λ*reg_f) \t chi2r: $(chi2_f/length(flux)) \t λ*reg: $(λ*reg_f) \t min(x): $(low_x) \t max(x): $(hi_x) \t Δx: $(Δx)");
# end
# g[:] =  chi2_g;
# g[somvis] += 2*λ*reg_g;
# return chi2_f + 2*λ*reg_f
# end
