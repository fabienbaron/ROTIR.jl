
#
# Try grid search of lambda & B
#

include("lci_fakemap.jl"); # this sets xtrue to be the true map
crit_imaging = (x,g)->chi2_lci_bias_fg(x, g, polyflux, lcidata, B, λ, somvis, verbose = false); # turn off verbose
ΔT_target = 1000;
minerr = 1e99;
dummy=zeros(star_epoch_geom[1].npix);
for B in 10.^linspace(2,10,10)
    println("B=$(B) \n");
    for λ in 10.^linspace(-8,4,15)
        x_opt = OptimPack.vmlmb(crit_imaging, temperature_map,lower=3000.0, verb=false, maxiter=400, blmvm=false);
        # for niter=1:3
        #     x_opt = OptimPack.vmlmb(crit_imaging, x_opt, verb=false, maxiter=300, blmvm=false);
        # end
        lowflux = norm(Int.(x_opt.<(0.9*mean(x_opt))),0);
        ΔT=maximum(x_opt[somvis])-minimum(x_opt[somvis]);
        chi2_r = chi2_lci_f(x_opt, polyflux, lcidata)/length(lcidata.flux)
        #if abs(chi2_r - 1)<3./sqrt(length(lcidata.flux))
            println("B=$(B) \t λ=$(λ) \t chi2_r =$(chi2_r) Crit=$(crit_imaging(x_opt,dummy)) \t minT=$(minimum(x_opt[somvis])) maxT=$(maximum(x_opt[somvis])) \t ΔT=$(ΔT)\n");
            if abs((ΔT-ΔT_target)/ΔT_target)<.3
            #        mollplot_temp_healpix(x_opt);
                    mollplot_temp_longlat(x_opt,ntheta,nphi);
                    print_with_color(:red, "Good solution B=$(B) \t λ=$(λ)\n")
            end
        #end
#        err = norm(x_opt[somvis]-x_true[somvis], 2);

#        if err<minerr
#            minerr = deepcopy(err)

        end
end




crit_imaging = (x,g)->chi2_lci_tv_fg(x, g, polyflux, lcidata, λ);
somvis = sometimes_visible(star_epoch_geom); # lists all pixels at least visible once
λvals = 10.^linspace(-8,6,40);
chi2_opt_r = zeros(100);
ΔT= zeros(100);
#x_optim = deepcopy(x_opt);
for i=1:40
    λ=λvals[i];
# Optimize temperature_map
    x_optim = OptimPack.vmlmb(crit_imaging, temperature_map, verb=false, maxiter=100, blmvm=false);
    for niter=1:5 # error with OptimPack here...
        x_optim = OptimPack.vmlmb(crit_imaging, x_optim, verb=false, maxiter=100, blmvm=false);
    end
     if typeof(x_optim) ==Void
         x_optim =  deepcopy(x_opt);
     end
    chi2_opt_r[i] = norm((modelflux_lci(x_optim, polyflux)-lcidata.flux)./lcidata.fluxerr)^2 /length(lcidata.flux);
    ΔT[i] = maximum(x_opt[somvis])-minimum(x_opt[somvis]);
    err = norm(x_optim[somvis] - x_true[somvis],1)/length(somvis);
    println("i=$i λ=$λ \t chi2r=$(chi2_opt_r[i]) ΔT =$(ΔT[i])  err =$err\n");
end


lciplot_vs_model(lcidata, modelflux_lci(x_optim, polyflux));
mollplot_temperature_healpix(x_optim);
plot2d_temperature(x_optim, star_epoch_geom[1]);
plot(log.(λvals),log.(chi2_opt_r[1:40]))
