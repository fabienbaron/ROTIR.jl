using PythonCall
using PythonPlot

# Prefer an interactive backend, but only when the user has not already chosen one.
# An unconditional assignment here runs at `using ROTIR` time and would force Qt onto
# headless sessions (CI, batch fits, the test suite), where importing it fails or hangs.
get!(ENV, "MPLBACKEND", "qt5agg")
rcParams = pyimport("matplotlib").rcParams
#rcParams["font.family"] = "serif"
#rcParams["font.serif"] = "Times New Roman"
# rcParams["font.size"] = 15
rcParams["xtick.top"] = true
rcParams["xtick.direction"] = "in"
rcParams["xtick.minor.visible"] = true
rcParams["ytick.right"] = true
rcParams["ytick.direction"] = "in"
rcParams["ytick.minor.visible"] = true

function lciplot_phase(lcidata)
errorbar(lcidata.phase, lcidata.flux,yerr=lcidata.fluxerr,fmt="o", color="Black")
end
function lciplot_mjd(lcidata)
errorbar(lcidata.mjd, lcidata.flux,yerr=lcidata.fluxerr,fmt="o", color="Black")
end
function lciplot_vs_model_phase(lcidata, modelflux)
fig = figure("Light curve - Data vs model",figsize=(15,10),facecolor="White")
errorbar(lcidata.phase, lcidata.flux,yerr=lcidata.fluxerr,fmt="o", color="Black")
plot(lcidata.phase, modelflux)
end

function lciplot_vs_model_mjd(lcidata, modelflux)
fig = figure("Light curve - Data vs model",figsize=(15,10),facecolor="White")
clf();
errorbar(lcidata.mjd, lcidata.flux,yerr=lcidata.fluxerr,fmt="o", color="Black")
plot(lcidata.mjd, modelflux)
end

function lciplot_vs_model_phase_residuals(lcidata, modelflux, relative=true)
    # fig = figure("Light curve - Data vs model",figsize=(15,10),facecolor="White")
    fig, ax = pyplot.subplots(2, 1, figsize=(15,10), sharex=true, gridspec_kw=Dict("height_ratios"=>[1.0, 0.3], "hspace"=>0.0))
    # 0-BASED indexing: `subplots` returns a Py numpy array of axes, not a Julia Matrix.
    ax[0].errorbar(lcidata.phase, lcidata.flux, yerr=lcidata.fluxerr, fmt="o", color="C0", label="Observed", zorder=1, mfc="none")
    ax[0].plot(lcidata.phase, modelflux, label="Inferred", zorder=2, c="C1")
    ax[1].axhline(0.0, c="k", ls="--", zorder=1)
    ax[1].plot(lcidata.phase, lcidata.flux .- modelflux, zorder=2)
    ax[1].set_ylabel("O-C")
    ax[1].set_xlabel("MBJD [days]")
    if relative==true
        ax[0].set_ylabel("Normalized Flux")
    else
        ax[0].set_ylabel("Flux [counts]")
    end
    ax[0].legend()
end 
