#
# Simple light curve inversion demo code
#

include("lci.jl");
include("lciplot.jl");
include("oiplot.jl");
include("oistars.jl");
include("geometry.jl");

using OptimPack

# Import light curve
#lcifile = "5110407_intphaselc_pcs_0001.txt";

#lcifile = "kplr005110407_LC_CBV_Q04.txt"
#lcidata = read_lci_absolute_mjd(lcifile);
ΔT_target = 1000;
Tphot_target = 5200;

#Debug
 lcifile = "KIC_flux.txt"
 lcidata = read_lci_absolute(lcifile);
lcidata.mjd = lcidata.phase* 3.4693
 #lcidata.mjd = vcat(lcidata.phase* 3.4693, 3.4693+lcidata.phase* 3.4693, 6.9386 +lcidata.phase* 3.4693)  ;
 #lcidata.flux = vcat(lcidata.flux, lcidata.flux, lcidata.flux)
 #lcidata.fluxerr = vcat(lcidata.fluxerr, lcidata.fluxerr, lcidata.fluxerr)
 #lcidata.nepochs = length(lcidata.flux)



# Define stellar parameters for each epoch
stellar_parameters = Array{starparameters}(lcidata.nepochs);
params = [1.,   # radius (fixed for LCI)
60.0, # inclination angle
0.0, # position angle
0,    # rotation angle (for LCI we use phase instead)
3.4693, # rotation period
[3,0.2,0.0],#[2, 0.7248, 0.1941],    # Limb darkening law 1: quadratic, 2: logarithmic, 3; Hestroffer, then coefficients
0,    # differential rotation coefficient 1
0     # differential rotation coefficient 2
];
for i=1:lcidata.nepochs
          stellar_parameters[i]=starparameters(params[1],params[2],params[3],360.0/params[5].*lcidata.mjd[i], 0*params[5], params[6], params[7], params[8]);
end

# Create 3D geometry from parameters
const nhealpix = 4; #Healpix tesselation level
star_epoch_geom = create_geometry( healpix_round_star(nhealpix), stellar_parameters);

# Initial temperature maps
const npix = star_epoch_geom[1].npix
const ntimes = 3
xranges = Array{UnitRange{Int64},1}(ntimes);
iepochranges = Int64.(round.(collect(linspace(0, length(lcidata.mjd), ntimes+1))));
for i=1:ntimes
          xranges[i]= iepochranges[i]+1:iepochranges[i+1]
end

# Surface precalc
polyflux = setup_lci(star_epoch_geom);
somvis = sometimes_visible(star_epoch_geom); # lists all pixels at least visible once
nevervis=never_visible(star_epoch_geom);
#x_init = 1000.*ones(ntimes*npix); #start from lowest possible temperature
x_init = vec(repmat(polyflux'*lcidata.flux, ntimes));

function chi2_lci_f(x, xranges, polyflux, lcidata)
          const nepochs = size(polyflux,1)
          const npix = size(polyflux,2)
          flux = Array{Float64}(nepochs);
          for i=1:length(xranges)
                    flux[xranges[i]] = polyflux[xranges[i],:]*x[npix*(i-1)+1:npix*i];
          end
          res = (flux-lcidata.flux)./lcidata.fluxerr;
          chi2_f = norm(res,2)^2 ;  # norm((modelflux_lci(x,polyflux)-lcidata.flux)./lcidata.fluxerr,2)^2
          return chi2_f
end

function modelflux_lci(x, xranges, polyflux)
          const nepochs = size(polyflux,1)
          const npix = size(polyflux,2)
          flux = Array{Float64}(nepochs);
          for i=1:length(xranges)
                    flux[xranges[i]] = polyflux[xranges[i],:]*x[npix*(i-1)+1:npix*i];
          end
          return flux
end

function chi2_lci_bias_fg(x, g, polyflux, lcidata, B, λ1, λ2, somvis; verbose = true)
          # Chi2 and gradient
          #g=zeros(size(x_init));x=deepcopy(x_init);
          const nepochs = size(polyflux,1)
          const npix = size(polyflux,2)
          const ntimes = length(xranges)
          const nvis = length(somvis);

          flux = Array{Float64}(nepochs);
          for i=1:ntimes
                    flux[xranges[i]] = polyflux[xranges[i],:]*x[npix*(i-1)+1:npix*i];
          end
          res = (flux-lcidata.flux)./lcidata.fluxerr;
          chi2_f = norm(res,2)^2 ;
          #println("Chi2: $(chi2_f/length(lcidata.flux))")
          # Gradient
          chi2_g = Array{Float64}(length(x))
          for i=1:ntimes
                    chi2_g[npix*(i-1)+1:npix*i] = 2.0*polyflux[xranges[i],:]'*((res./lcidata.fluxerr)[xranges[i]])
          end
          g[:] = chi2_g;

          # Harmon regularization (independent for each epoch)
          reg_f = 0.0;
          for i=1:ntimes
                    xvis = x[npix*(i-1)+somvis];
                    avg_xvis = mean(xvis);
                    bcorr = (B-1.0)*Int.(xvis.>avg_xvis)+1.0
                    reg_f += sum(bcorr.*(xvis-avg_xvis).^2)/nvis;
                    reg_g = 2*(xvis-avg_xvis).*bcorr/nvis;
                    reg_g = reg_g - mean(reg_g);
                    g[npix*(i-1)+somvis] += λ1*reg_g;
          end
          #
          # Group normalization
          # tv_f = 0;
          # if(ntimes>1)
          #           y = reshape(x,(npix,ntimes))
          #           yy = sqrt.(sum(y.^2,2))
          #           tv_f = sum( yy );
          #           tv_g = vec(y./yy)
          #           g[:] += λ2*tv_g;
          # end

          # cross-temporal L2 diff
          tv_f = 0;
          if(ntimes>1)
                  y = reshape(x,(npix,ntimes))
                  tv_f = sum( (y[:,2:end]-y[:,1:end-1]).^2 )
                  tv_g = Array{Float64}(npix,ntimes)
                  if(ntimes >2)
                    tv_g[:,1] = 2*(y[:,1] - y[:,2])
                    tv_g[:,2:end-1] = 4*y[:,2:end-1]-2*(y[:,1:end-2]+y[:,3:end])
                    tv_g[:,end] = 2*(y[:,end] - y[:,end-1])
                  else
                    tv_g[:,1] = 2*(y[:,1]-y[:,2]);
                    tv_g[:,2] = 2*(y[:,2]-y[:,1]);
                  end
                  g[:] += λ2*vec(tv_g);
          end

          return chi2_f + λ1*reg_f+ λ2*tv_f
end

function chi2_lci_f(x, polyflux, lcidata; verbose = true)
          const npix = size(polyflux,2)
          const ntimes = length(xranges)
          const nepochs = size(polyflux,1)
          flux = Array{Float64}(nepochs);
          for i=1:ntimes
                    flux[xranges[i]] = polyflux[xranges[i],:]*x[npix*(i-1)+1:npix*i];
          end
          res = (flux-lcidata.flux)./lcidata.fluxerr;
          chi2_f = norm(res,2)^2 ;
          return chi2_f
end


#n=30;x=rand(n);y=deepcopy(x);bcorr = (B-1.0)*Int.(x.>mean(x))+1.0;y[2]+=1e-7;g=2*bcorr.*(x.-mean(x))/n;g=g-mean(g);g[2]-(sum(bcorr.*(y-mean(y)).^2)/n-sum(bcorr.*(x-mean(x)).^2)/n)/1e-7
#debug initial situation
λ1=0.0002; # bias regularization hyperparameter
λ2= 1e-2; # transpectral regularization hyperparameter
B = 6000;  # bias factor from Harmonλ=0.0002

crit_imaging = (x,g)->chi2_lci_bias_fg(x, g, polyflux, lcidata, B, λ1, λ2, somvis);
# Optimize temperature_map
x_opt = OptimPack.vmlmb(crit_imaging, x_init, lower=0.1, verb=true, maxiter=1000, blmvm=false); # note: lower temp needs to be < min

#x_opt[never_visible(star_epoch_geom)] = mean(x_opt[sometimes_visible(star_epoch_geom)]);
lciplot_vs_model_mjd(lcidata, modelflux_lci(x_opt, xranges, polyflux));

x_opt = reshape(x_opt,(npix,ntimes));

# overwrite unseen flux with average then rescale; for display only
#x_opt[nevervis,:]=repmat(mean(x_opt[somvis,:],1),length(nevervis))
#x_opt .*= Tphot_target/mean(x_opt[somvis,:])
mollplot_temperature(x_opt[:,1]);


#include("lci_fakemap.jl")

dummy=zeros(size(x_init))
for B in 10.^linspace(-2,10,12)
          println("B=$(B) \n");
          for λ1 in 10.^linspace(-10,4,15)
                    x_opt = OptimPack.vmlmb(crit_imaging, x_init,lower=0.1, verb=false, maxiter=500, blmvm=false);
                    chi2_r = chi2_lci_f(x_opt, polyflux, lcidata)/length(lcidata.flux)
                    x_opt .*= Tphot_target/mean(x_opt[somvis,:])
                    minT= minimum(x_opt[somvis,:]);
                    maxT= maximum(x_opt[somvis,:]);
                    ΔT_model = Tphot_target - minT
                    #if abs(chi2_r - 1)<6./sqrt(length(lcidata.flux))
                              println("B=$(B) \t λ1=$(λ1) \t λ2=$(λ2)\t chi2_r =$(chi2_r) Crit=$(crit_imaging(x_opt,dummy)) \t minT=$(minimum(x_opt[somvis,:])) maxT=$(maximum(x_opt[somvis,:])) \t δ=$(δ)");
                              if abs((ΔT_model-ΔT_target)/ΔT_target)<.2
                                        x_opt = reshape(x_opt,(npix,ntimes));
                                        # overwrite unseen flux with average then rescale; for display only
                                        x_opt[nevervis,:]=repmat(mean(x_opt[somvis,:],1),length(nevervis))
                                        x_opt .*= Tphot_target/mean(x_opt[somvis,:])
                                        mollplot_temperature(x_opt[:,1]);
                                        print_with_color(:red, "Good solution B=$(B) \t λ1=$(λ1)\n")
                                        # print_with_color(:red, "Good solution B=$(B) \t λ1=$(λ1) dist ", norm(x_opt[somvis,1]-x_true[somvis],2),"\n")
                              end
                    #end
          end
end
