# include("di.jl");
using FFTW
using PyPlot
using FITSIO
using Dierckx
using Statistics
using DelimitedFiles
# using LsqFit

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.family"] = "serif"
#rcParams["font.serif"] = "Times New Roman"
rcParams["xtick.top"] = true
rcParams["xtick.direction"] = "in"
rcParams["xtick.minor.visible"] = true
rcParams["ytick.right"] = true
rcParams["ytick.direction"] = "in"
rcParams["ytick.minor.visible"] = true

function doppler_shift_spectrum(wavelength::Vector{Float64}, intensity::Vector{Float64}, v::Float64)
    return Spline1D(wavelength .* (1.0 + v/3E5), intensity, k=1, bc="nearest", s=0.0)(wavelength)
end

function normalize_UVES(file, rv, output; regions=[], write=false, plot=false)
    println("Normalizing $file")

    f = FITS(file, "r+")
    mjd = read_key(f[2], "TMID")[1]
    vhelio = read_key(f[1], "ESO QC VRAD HELICOR")[1]
    vbary  = read_key(f[1], "ESO QC VRAD BARYCOR")[1]
    SNR = read_key(f[1], "SNR")[1]
    R = read_key(f[1], "SPEC_RES")[1]
    wavelength = vec(read(f[2], "WAVE"))
    intensity = convert(Vector{Float64}, vec(read(f[2], "FLUX")[:, 1]))
    close(f)
    
    intensity_cut = intensity[findall(intensity .> 0)]
    wavelength_cut = wavelength[findall(intensity .> 0)]

    intensity_shifted = doppler_shift_spectrum(wavelength_cut, intensity_cut, vbary - rv)
    continuum = fit_continuum(wavelength_cut, intensity_shifted, mode="spline", regions=regions, sorder=9, Nσ_low=0.5)
    intensity_normalized = intensity_shifted ./ continuum

    fig, ax = subplots(2, 1, sharex=true, gridspec_kw=Dict("height_ratios" => [3.0, 1.0], "hspace" => 0.0))
    ax[1].plot(wavelength_cut, intensity_shifted)
    ax[1].plot(wavelength_cut, continuum)
    ax[2].plot(wavelength_cut, intensity_normalized)
    ax[2].axhline(1.0, c="k", ls="--")

    ax[1].set_ylabel(raw"F$_\lambda$ [erg/s/cm$^2$/Å]")
    ax[2].set_ylabel("Normalized")
    ax[2].set_xlabel("Wavelength [Å]")

    if plot == true
        PyPlot.show()
    end

    header = "# MJD:\t$(mjd)\tSNR:\t$(SNR)\tR:\t$(R)\n"
    if write == true
        if plot == true
            print("Accept? [y]/n: ")
            accept = readline()
        else
            accept = "y"
        end

        if accept == ""
            accept = "y"
        end

        if accept == "y"
            open(output, "w") do o
                Base.write(o, header)
                writedlm(o, [wavelength_cut intensity_normalized], "\t")
            end
            savefig(replace("$(output)", ".dat"=>".png"), dpi=300, bbox_inches="tight")
            println("Written to $output")
        end
    end

    return wavelength_cut, intensity_normalized, header
end

function fit_continuum(wavelength, intensity; porder=3, sorder=5, Nσ_high=99.0, Nσ_low=2.0, niter=10, mode="spline", regions=[], mflx=[])
    if regions != []
        ixs = Array{Any, 1}(undef, length(regions))
        for i=1:length(regions)
            ix = findall(regions[i][1] .< wavelength .< regions[i][2])
            ixs[i] = ix
        end
        ixs = vcat(ixs...)
    else
        ixs = 1:length(wavelength)
    end

    if mode == "spline"
        ixknots = Int.(floor.(range(1, stop=length(wavelength[ixs]), length=sorder+1)))[2:end-1]
        xknots = wavelength[ixs][ixknots]

        wavelength_sub = deepcopy(wavelength[ixs])
        intensity_sub = deepcopy(intensity[ixs])
        for n=1:niter
            global f
            f = Spline1D(wavelength_sub, intensity_sub, xknots, k=porder)
            continuum = f(wavelength_sub)
            res = intensity_sub .- continuum
            σ = std(res)
            ixrej = findall(-Nσ_low*σ .< res .< Nσ_high*σ)
            wavelength_sub = wavelength_sub[ixrej]
            intensity_sub = intensity_sub[ixrej]
            # plot(wavelength, intensity, "-")
            # plot(wavelength_sub, f(wavelength_sub), "-")
            # show()
        end
        f = Spline1D(wavelength_sub, intensity_sub, xknots, k=porder)
        return f(wavelength)
    elseif mode == "median"
        Δλ = 0.01
        bin_width = 100.0
        nixs = bin_width / Δλ  # Å
        nbins = Int(floor((maximum(wavelength) - minimum(wavelength)) / bin_width))
        median_bins = zeros(nbins+2)
        wavelength_bins = zeros(nbins+2)
        wavelength_bins[1] = minimum(wavelength); median_bins[1] = median(intensity[1:Int(nixs÷2)])
        wavelength_bins[end] = maximum(wavelength); median_bins[end] = median(intensity[end-Int(nixs÷2):end])
        for i=1:nbins
            ix1 = Int((i-1)*nixs + 1)
            ix2 = Int(i*nixs + 1)
            wavelength_bins[i+1] = mean(wavelength[ix1:ix2])
            median_bins[i+1] = median(intensity[ix1:ix2])
        end
        f = Spline1D(wavelength_bins, median_bins)
        return f(wavelength)
    #= elseif mode == "polynomial"
        wavelength_sub = deepcopy(wavelength[ixs])
        intensity_sub = deepcopy(intensity[ixs])
        for n=1:niter
            global f
            f = polyfit(wavelength_sub, intensity_sub, porder)
            continuum = f(wavelength_sub)
            res = intensity_sub .- continuum
            σ = std(res)
            ixrej = findall(-Nσ_low*σ .< res .< Nσ_high*σ)
            wavelength_sub = wavelength_sub[ixrej]
            intensity_sub = intensity_sub[ixrej]
            # plot(wavelength, intensity, "-")
            # plot(wavelength_sub, f(wavelength_sub), "-")
            # show()
        end
        f = polyfit(wavelength_sub, intensity_sub, porder)
        return f(wavelength) =#
    end

end

function fit_rv(wvl, flx, mwvl, mflx; grow=20)
    ccf = real.(fftshift(ifft(fft(flx) .* conj.(fft(mflx)))))
    ccf .-= mean(ccf)
    ccf ./= maximum(ccf)
    maxval, maxix = findmax(ccf)
    lnλ0 = log(mean(wvl))
    Δlnλ = log.(wvl) .- lnλ0
    Δv = Δlnλ .* 3e5

    gaussian(x, p) = exp.( -0.5 .* ((x .- p[1]) ./ p[2]).^2 )

    fit = curve_fit(gaussian, Δv[maxix-grow:maxix+grow], ccf[maxix-grow:maxix+grow], [Δv[maxix], 1.0])
    rv = coef(fit)[1]
    return rv
end