#include("di.jl")
#include("lci.jl")
#include("oichi2_spheroid.jl")
using Crayons

function reconstruct(x_start::Vector{Float64}; diparams=[], lciparams=[], oiparams=[], w=[1.0, 1.0, 1.0], lowtemp=0, hightemp=10000, printcolor=[], verb=true, maxiter=200, regularizers=[])
    ## Use VMLMB from OptimPackNextGen to find the temperature map that minimizes the χ^2 between the observed profiles and the models
    x_sol = [];
    crit_imaging = x -> x; #TODO: proper dummy init
    crit_imaging = (x, g) -> crit_fg(x, g, diparams=diparams, lciparams=lciparams, oiparams=oiparams, w=w, regularizers=regularizers);
    # global N = 1
    x_sol = OptimPackNextGen.vmlmb(crit_imaging, x_start, verb=false, lower=lowtemp, upper=hightemp, maxiter=maxiter, blmvm=false, gtol=(0, 1e-8));
    # x_sol = OptimPack.nlcg(crit_imaging, x_start, verb=false, maxiter=maxiter)
    return x_sol
end

function crit_fg(x, g; diparams=[], lciparams=[], oiparams=[], w=Dict("OI"=>1.0, "DI"=>1.0, "LCI"=>1.0), verb=true, regularizers=[])    
    
    g_oi = zeros(length(x))
    if (oiparams != []) && (w["OI"] != 0.0)
        f_oi = spheroid_crit_multiepochs_fg(x, g_oi, oiparams["polyflux"], oiparams["polyft"], oiparams["data"]; regularizers=[], verb=false)
        f_oi_red = f_oi / length(x)
    else
        f_oi = f_oi_red = 0.0
    end
    
    g_di = zeros(length(x))
    if (diparams != []) && (w["DI"] != 0.0)
        f_di = di_crit_fg(x, g_di, diparams["didata"], diparams["polyflux"], diparams["star_epoch_geom"], diparams["star_epoch_μ"], diparams["star_epoch_v"], diparams["stellar_parameters"], diparams["model_grid"], diparams["vsini"], diparams["visible_pixels"]; regularizers=[], verb = false, wavelength_regions=diparams["wavelength_regions"], normalization_region=diparams["normalization_region"])
        f_di_red = f_di / diparams["didata"].nepochs / diparams["didata"].nwvls
    else
        f_di = f_di_red = 0.0
    end

    g_lci = zeros(length(x))
    if (lciparams != []) && (w["LCI"] != 0.0)
        relative = lciparams["relative"]
        if relative == true
            f_lci = lci_relative_crit_fg(x, g_lci, lciparams["dataflux"], lciparams["datafluxerr"], lciparams["polyflux"], lciparams["visible_pixels"]; regularizers=[], verb = false)
        else
            f_lci = lci_absolute_crit_fg(x, g_lci, lciparams["dataflux"], lciparams["datafluxerr"], lciparams["polyflux"], lciparams["visible_pixels"]; regularizers=[], verb = false)
        end
        f_lci_red = f_lci / length(x)
    else
        f_lci = f_lci_red = 0.0
    end

    f = f_oi*w["OI"] + f_di*w["DI"] + f_lci*w["LCI"]
    g .= g_oi.*w["OI"] .+ g_di.*w["DI"] .+ g_lci.*w["LCI"]
    
    f_reg = spheroid_regularization(x, g, regularizers=regularizers, verb = false)
    f += f_reg

    if verb == true
        if oiparams != []
            y = x[oiparams["visible_pixels"]]
        elseif diparams != []
            y = x[diparams["visible_pixels"]]
        elseif lciparams != []
            y = x[lciparams["visible_pixels"]]
        end
        low_x = minimum(y);
        high_x = maximum(y);
        mean_x = mean(y);
        printinfo(f, w, f_oi, f_oi_red, f_di, f_di_red, f_lci, f_lci_red, f_reg, low_x, high_x, mean_x)
        # @printf("f_tot: %8.2f\tf_oi: %8.2f (%2.2f)\tf_di: %8.2f (%2.2f)\tf_lci: %8.2f (%2.2f)\tf_reg: %8.2f\tminT: %5.1f\tmaxT: %5.1f\tmeanT: %5.1f\n", f, w["OI"]*f_oi, f_oi_red, w["DI"]*f_di, f_di_red, w["LCI"]*f_lci, f_lci_red, f_reg, low_x, high_x, mean_x);
    end

    return f
end

function printinfo(f, w, f_oi, f_oi_red, f_di, f_di_red, f_lci, f_lci_red, f_reg, low_x, high_x, mean_x)
    s1 = @sprintf("f_tot: %8.2f\t", f)
    print(crayon"214", s1)
    s2 = @sprintf("f_oi: %8.2f (%2.2f)\t", w["OI"]*f_oi, f_oi_red)
    print(crayon"red", s2)
    s3 = @sprintf("f_di: %8.2f (%2.2f)\t", w["DI"]*f_di, f_di_red)
    print(crayon"green", s3)
    s4 = @sprintf("f_lci: %8.2f (%2.2f)\t", w["LCI"]*f_lci, f_lci_red)
    print(crayon"blue", s4)
    s5 = @sprintf("f_reg: %8.2f\tminT: %5.1f\tmaxT: %5.1f\tmeanT: %5.1f\n", f_reg, low_x, high_x, mean_x)
    print(crayon"white", s5)
end