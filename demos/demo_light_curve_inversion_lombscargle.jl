#
# Finding a period in a light curve using Lomb Scargle
#
include("../src/ROTIR.jl"); using Main.ROTIR;
# Import light curve
lcifile = "./data/kplr005110407_LC_CBV_Q02.txt"
lcidata = read_lci_absolute(lcifile);
#lciplot_mjd(lcidata)

# Determine period (basic)
using LombScargle
pgram = lombscargle(lcidata.mjd, lcidata.flux, lcidata.fluxerr, maximum_frequency=1.0)
period = findmaxperiod(pgram)
using PyPlot
plot(freqpower(pgram)...)
println("Period determined by Lomb Scargle: $(period) days.")
