using DelimitedFiles
using LinearAlgebra

struct LCI
  mjd::Array{Float64,1}
  phase::Array{Float64,1}
  flux::Array{Float64,1}
  fluxerr::Array{Float64,1}
  nepochs::Int64
end

function read_lci_relative(file)
lc = readdlm(file, skipstart=6)
phase = vec(lc[:,1]);
flux = vec(lc[:,2]);
fluxerr = ones(size(flux))
return LCI([],phase, flux, fluxerr, length(phase))
end

function read_lci_absolute(file)
lc = readdlm(file, skipstart=6)
phase = vec(lc[:,1]);
flux = vec(lc[:,2]);
fluxerr = vec(lc[:,3]);
return LCI([],phase, flux, fluxerr, length(phase))
end

function read_lci_absolute_mjd(file)
lc = readdlm(file, skipstart=6)
mjd = vec(lc[:,1]);
flux = vec(lc[:,2]);
fluxerr = vec(lc[:,3]);
return LCI(mjd, [], flux, fluxerr, length(mjd))
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
return polyflux
end


function chi2_lci_f(x, polyflux, lcidata)
flux = polyflux*x;
res = (flux-lcidata.flux)./lcidata.fluxerr;
chi2_f = norm(res,2)^2 ;  # norm((modelflux_lci(x,polyflux)-lcidata.flux)./lcidata.fluxerr,2)^2
return chi2_f
end

function lci_absolute_crit_fg(x, g,  lcidata, polyflux, visible_pixels; regularizers=[], verb = true )
flux = polyflux*x;
res = (flux-lcidata.flux)./lcidata.fluxerr;
chi2_f = norm(res,2)^2 ;  # norm((modelflux_lci(x,polyflux)-lcidata.flux)./lcidata.fluxerr,2)^2
g[:] = 2.0*polyflux'*(res./lcidata.fluxerr);
reg_f = spheroid_regularization(x, g, regularizers=regularizers, verb = false); #note: this adds to g
if verb == true
  y= x[visible_pixels]
  low_x = minimum(y);
  hi_x = maximum(y);
  Δx = hi_x - low_x
  print("Crit: $(chi2_f+reg_f) \t chi2r: $(chi2_f/length(flux)) \t λ*reg: $(reg_f) \t min(x[vis]): $(low_x) \t max(x[vis]): $(hi_x) \t Δx: $(Δx) ");
end
return chi2_f + reg_f
end


function lci_relative_crit_fg(x, g,  lcidata, polyflux, visible_pixels; regularizers=[] , verb = true)
flux = polyflux*x;
flux_avg = mean(flux);
res = (flux/flux_avg-lcidata.flux)./lcidata.fluxerr;
chi2_f = norm(res,2)^2 ;  # norm((modelflux_lci(x,polyflux)-lcidata.flux)./lcidata.fluxerr,2)^2
g[:] = 2.0/flux_avg*polyflux'*(res./lcidata.fluxerr)-2.0/flux_avg^2*sum(flux.*(res./lcidata.fluxerr))*vec(mean(polyflux,dims=1));
reg_f = spheroid_regularization(x, g, regularizers=regularizers, verb = false); #note: this adds to g
if verb == true
y= x[visible_pixels]
low_x = minimum(y);
hi_x = maximum(y);
Δx = hi_x - low_x
println("Crit: $(chi2_f) \t chi2r: $(chi2_f/length(flux)) \t λ*reg: $(reg_f) \t min(x[vis]): $(low_x) \t max(x[vis]): $(hi_x) \t Δx: $(Δx)");
end
return chi2_f + reg_f
end


function lci_reconstruct(x_start::Array{Float64,1}, lcidata::LCI, polyflux::Array{Float64,2}, visible_pixels; lowtemp = 1000, relative = true, printcolor= [], verb = true, maxiter = 100, regularizers =[])
  x_sol = [];
  crit_imaging=x->x; #TODO: proper dummy init
  if relative == true
  crit_imaging = (x,g)->lci_relative_crit_fg(x, g, lcidata, polyflux, visible_pixels, regularizers=regularizers, verb = verb);
  else
    crit_imaging = (x,g)->lci_absolute_crit_fg(x, g, lcidata, polyflux, visible_pixels, regularizers=regularizers, verb = verb);
  end
  x_sol = OptimPackNextGen.vmlmb(crit_imaging, x_start, verb=verb, lower=lowtemp, maxiter=maxiter, blmvm=false, gtol=(0,1e-8));
  return x_sol
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
