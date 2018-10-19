using PyCall
using PyPlot

function lciplot_phase(lcidata)
errorbar(lcidata.phase, lcidata.flux,yerr=lcidata.fluxerr,fmt="o", color="Black")
end
function lciplot_mjd(lcidata)
errorbar(lcidata.phase, lcidata.flux,yerr=lcidata.fluxerr,fmt="o", color="Black")
end
function lciplot_vs_model_phase(lcidata, modelflux)
fig = figure("Light curve - Data vs model",figsize=(15,10),facecolor="White")
errorbar(lcidata.phase, lcidata.flux,yerr=lcidata.fluxerr,fmt="o", color="Black")
plot(lcidata.phase, modelflux)
end

function lciplot_vs_model_mjd(lcidata, modelflux)
fig = figure("Light curve - Data vs model",figsize=(15,10),facecolor="White")
errorbar(lcidata.mjd, lcidata.flux,yerr=lcidata.fluxerr,fmt="o", color="Black")
plot(lcidata.mjd, modelflux)
end
