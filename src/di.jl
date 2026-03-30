using FFTW
using FITSIO
using Printf
using Dierckx
using AstroLib
using PrettyTables
using LinearAlgebra
using Interpolations
using OptimPackNextGen

using BenchmarkTools

ENV["MPLBACKEND"] = "qt5agg"
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
# rcParams["font.family"] = "serif"
# rcParams["font.serif"] = "Times New Roman"
# rcParams["font.size"] = 15
rcParams["xtick.top"] = true
rcParams["xtick.direction"] = "in"
rcParams["xtick.minor.visible"] = true
rcParams["ytick.right"] = true
rcParams["ytick.direction"] = "in"
rcParams["ytick.minor.visible"] = true

mutable struct modelGrid
    ## Structure to hold model grid
    temperatures::Array{Float64, 1}
    logg::Float64
    mus::Array{Float64, 1}
    nwvls::Int64
    wavelength::Vector{Float64}
    intensities::Array{Vector{Float64}, 3}
    continuums::Array{Vector{Float64}, 2}
    intensities_itp::Interpolations.GriddedInterpolation
    continuums_itp::Interpolations.GriddedInterpolation
    intensities_gradient_itp
    continuums_gradient_itp
    dRdT::Interpolations.GriddedInterpolation
end

mutable struct globalLineProfiles
    ## Structure to hold model line profiles
    npix::Int64
    nepochs::Int64
    nwvls::Int64
    wavelength::Vector{Float64}
    intensities::Array{Vector{Float64}, 1}
end

mutable struct observedLineProfiles
    ## Structure to hold data line profiles
    fname::Array{String, 1}
    nepochs::Int64
    nwvls::Int64
    snr::Array{Float64, 1}
    mjd::Array{Float64, 1}
    phase::Array{Float64, 1}
    wavelength::Vector{Float64}
    intensities::Array{Vector{Float64}, 1}
end

function Base.display(profiles::observedLineProfiles)
    #= @printf("Nepochs: %3d \tNwvls: %5d\n", profiles.nepochs, profiles.nwvls)
    println();
    println("MJD\tPhase\tSNR\tλmin\tλmax\tΔλ");
    println("---------------------------------------------------------------------");
    for i=1:profiles.nepochs
        @printf("%5.4f\t%0.4f\t%3.1f\t%4.2f\t%4.2f\t%1.2f\n", profiles.mjd[i], profiles.phase[i], profiles.snr[i], minimum(profiles.wavelength), maximum(profiles.wavelength), (maximum(profiles.wavelength) - minimum(profiles.wavelength))/profiles.nwvls)
    end =#
    #= ixs = sortperm(profiles.phase)
    profiles.fname .= profiles.fname[ixs]
    profiles.mjd .= profiles.mjd[ixs]
    profiles.snr .= profiles.snr[ixs]
    profiles.intensities .= profiles.intensities[ixs]
    profiles.phase .= profiles.phase[ixs] =#

    folder = splitdir(profiles.fname[1])[1]
    fnames = [splitdir(profiles.fname[i])[end] for i=1:profiles.nepochs]

    minλ = [minimum(profiles.wavelength) for i=1:profiles.nepochs]
    maxλ = [maximum(profiles.wavelength) for i=1:profiles.nepochs]
    Δλ = [(maximum(profiles.wavelength) - minimum(profiles.wavelength))/profiles.nwvls for i=1:profiles.nepochs]
    data = hcat(fnames, profiles.mjd, profiles.phase, profiles.snr, minλ, maxλ, Δλ)
    header = (["Filename", "MJD", "Phase", "SNR", "λmin", "λmax", "Δλ"],
                ["$(folder)", "", "", "", "[Å]", "[Å]", "[Å]"])
    pretty_table(data, header=header, header_crayon=crayon"#ffaf00", title="Observed Spectra", title_alignment=:c)
    println()
end

function setup_di(star_epoch_geom; ld = true)
    # Polyflux is the weight of each Healpixel, proportional to the surface
    nepochs = length(star_epoch_geom);
    polyflux = zeros(Float64, nepochs, star_epoch_geom[1].npix);
    for i=1:nepochs
        polyflux[i,star_epoch_geom[i].index_quads_visible] =0.5*(star_epoch_geom[i].proj_west[:,1].*star_epoch_geom[i].proj_north[:,2]
        - star_epoch_geom[i].proj_west[:,2].*star_epoch_geom[i].proj_north[:,1]
        + star_epoch_geom[i].proj_west[:,2].*star_epoch_geom[i].proj_north[:,3]
        - star_epoch_geom[i].proj_west[:,3].*star_epoch_geom[i].proj_north[:,2]
        + star_epoch_geom[i].proj_west[:,3].*star_epoch_geom[i].proj_north[:,4]
        - star_epoch_geom[i].proj_west[:,4].*star_epoch_geom[i].proj_north[:,3]
        + star_epoch_geom[i].proj_west[:,4].*star_epoch_geom[i].proj_north[:,1]
        - star_epoch_geom[i].proj_west[:,1].*star_epoch_geom[i].proj_north[:,4]);
        if ld == true
          polyflux[i,:] .*= star_epoch_geom[i].ldmap;
        end
    end
    visible_pixels = sometimes_visible(star_epoch_geom); # lists all pixels at least visible once
    hidden_pixels =  never_visible(star_epoch_geom);
    return polyflux, visible_pixels, hidden_pixels
end

function make_circ_spot_DI(temperature_map,star_geometry,spot_radius,lat,long;bright_frac=0.8)
    ## Makes a circular spot at a given (lat, long)
    ##  spot_radius ∈ [0, 360]
    ##  lat ∈ [-90, 90]
    ##  long ∈ [-180, 180]
    ##  bright_frac ∈ [0, ∞)

    long = mod(long, 360)

    temperature_map_copy = deepcopy(temperature_map)
    centers = star_geometry.vertices_spherical[:,5,:]
    θ = centers[:, 3]
    ϕ = -centers[:, 2] .+ pi/2

    da = 2.0*asin.(sqrt.(sin.(0.5*(ϕ .- lat*pi/180)).^2+cos(lat*pi/180).*cos.(ϕ).*sin.(0.5*(θ .- long*pi/180)).^2))

    spot_mask = findall(da .<= spot_radius*pi/180)
    temperature_map_copy[spot_mask] .= median(temperature_map_copy) * bright_frac
    return vec(temperature_map_copy)
end

function make_meridional_band_DI(temperature_map,star_geometry,band_width,long;bright_frac=0.8)
    ## Makes a meridional band at a given long
    ##  spot_radius ∈ [0, 360]
    ##  long ∈ [-180, 180]
    ##  bright_frac ∈ [0, ∞)

    long = mod(long, 360)

    temperature_map_copy = deepcopy(temperature_map)
    centers = star_geometry.vertices_spherical[:,5,:]
    θ = centers[:, 3]
    
    ix = findall(long - band_width/2 .<= θ .* 180/pi .<= long + band_width/2)
    temperature_map_copy[ix] .*= bright_frac

    return vec(temperature_map_copy)
end

function make_equatorial_band_DI(temperature_map,star_geometry,band_width,lat;bright_frac=0.8)
    ## Makes an equatorial band at a given lat
    ##  spot_radius ∈ [0, 360]
    ##  lat ∈ [-90, 90]
    ##  bright_frac ∈ [0, ∞)

    temperature_map_copy = deepcopy(temperature_map)
    centers = star_geometry.vertices_spherical[:,5,:]
    ϕ = -centers[:, 2] .+ pi/2

    ix = findall(lat - band_width/2 .<= ϕ .* 180/pi .<= lat + band_width/2)
    temperature_map_copy[ix] .*= bright_frac

    return vec(temperature_map_copy)
end

function compute_grid_atlas9_SPECTRUM(temperatures, logg, mus, wavelength_min, wavelength_max, dwvl, grid_directory; supergrid="ap00k2odfnew.dat", linelist="luke.lst", vmicro=0.8, vmacro=2.0)
    ## Computes a grid of spectra using ATLAS9 atmospheres and the SPECTRUM radiative transfer code
    cwd = pwd()
    cd(grid_directory)
    for i=eachindex(temperatures)
        T = Int(temperatures[i])
        atmosphere_outname = string("at", T, "g", replace(string(logg), "."=>""),".dat")
        # println(outname)
        redirect_stdout(devnull) do
            run(`selectmod $supergrid $atmosphere_outname $T $logg`)
        end

        atmosphere_inname = atmosphere_outname
        for j=eachindex(mus)
            mu = mus[j]
            println("Temp: $T, μ: $mu")

            spectrum_outname = string("iabst", T, "g", replace(string(logg), "."=>""), "mu", mu, ".spc")
            continuum_outname = string("cabst", T, "g", replace(string(logg), "."=>""), "mu", mu, ".spc")

            open("iresponse", "w") do f  # Input file for specific intensity computation
                write(f, string(atmosphere_inname, "\n"))
                write(f, string(linelist, "\n"))
                write(f, string(spectrum_outname, "\n"))
                write(f, string(vmicro, "\n"))
                write(f, string(mu, "\n"))
                write(f, string(wavelength_min, ",", wavelength_max, "\n"))
                write(f, string(dwvl, "\r"))
            end
            open("cresponse", "w") do f  # Input file for continuum computation
                write(f, string(atmosphere_inname, "\n"))
                write(f, string(linelist, "\n"))
                write(f, string(continuum_outname, "\n"))
                write(f, string(vmicro, "\n"))
                write(f, string(mu, "\n"))
                write(f, string(wavelength_min, ",", wavelength_max, "\n"))
                write(f, string(dwvl, "\r"))
            end

            open("run.sh", "w") do f
                write(f, "spectrum nM < iresponse && spectrum nMc < cresponse")
            end

            vmacro_outname = replace(spectrum_outname, ".spc" => "_turb.spc")

            redirect_stdout(devnull) do  # Run SPECTRUM and apply macroturbulence in SPECTRUM
                run(`sh run.sh`)
                run(`macturb $spectrum_outname $vmacro_outname $dwvl $vmacro`)
            end
        end
    end
    cd(cwd)
end

function read_spectrum_grid(grid_directory, temperatures, mus, logg)
    ## Reads in model grid files
    cwd = pwd()
    cd(grid_directory)
    wavelengths = Array{Vector{Float64}}(undef, length(temperatures), length(mus))
    intensities = Array{Vector{Float64}}(undef, length(temperatures), length(mus))
    continuums = Array{Vector{Float64}}(undef, length(temperatures), length(mus))
    @views for i=eachindex(temperatures)
        T = temperatures[i]
        for j=eachindex(mus)
            μ = mus[j]
            spectrum_filename = string("iabst$T","g",replace(string(logg),"."=>""),"mu$(μ)_turb.spc")
            continuum_filename = string("cabst$T","g",replace(string(logg),"."=>""),"mu$(μ).spc")

            spectrum = readdlm(spectrum_filename)
            continuum = readdlm(continuum_filename)
            wavelengths[i, j] = vec(spectrum[:, 1])
            if μ == 0.0
                intensities[i, j] = zeros(length(wavelengths[i, j]))
                continuums[i, j] = zeros(length(wavelengths[i, j]))
            else
                intensities[i, j] = vec(spectrum[:, 2])
                continuums[i, j] = vec(continuum[:, 2])
            end
        end
    end
    cd(cwd)
    return wavelengths[1, 1], intensities, continuums
end

function gaussian_instrumental_kernel(wavelength, resolution)
    ## Creates a gaussian broadening kernel to match models to observed spectral resolution
    # μ = mean(wavelength)
    μ = wavelength[length(wavelength)÷2 + (1 - length(wavelength)%2)]
    # println(μ, " ", x[length(x)÷2 - (length(x)%2)])
    FWHM = μ / resolution
    σ = FWHM / (2.0 * sqrt(2.0 * log(2.0)))
    kernel = 1.0/sqrt(2.0*pi*σ^2) .* exp.(-1.0/2.0 .* (μ .- wavelength).^2 ./ σ^2)
    return kernel ./ sum(kernel)
end

function instrumental_broadening(wavelength, intensity, resolution; kernel=[], pad=50)
    ## Applies gaussian broadening kernel to match models to observed spectral resolution
    ## Pads spectrum in zeros to minimize edge effects
    if kernel == []
        kernel = gaussian_instrumental_kernel(wavelength, resolution)
    end

    kernel_pad = cat(zeros(pad), kernel, zeros(pad), dims=1)
    intensity_pad = cat(zeros(pad), intensity, zeros(pad), dims=1)
    intensity_broad = real(fftshift(ifft(fft(intensity_pad) .* fft(kernel_pad))))
    intensity_broad_unpad = intensity_broad[pad+1:end-pad]
    return intensity_broad_unpad
end

function setup_model_grid(temperatures, logg, mus, vs, wavelength, intensities, continuums; resolution=0, wavelength_bins=[], wavelength_regions=[], verb=true)
    if verb==true
        # println("log(g) \t $(logg)")
        data = ["log(g)" "$(logg)"; "T[K]" "$(temperatures)"; "μ" "$(mus)"; "v[km/s]" "$(vs)"]
        pretty_table(data, noheader=true, title="Model Grid", title_alignment=:c)
        println()
    end
    
    ## Creates model grid structure and computes gradients
    #= if wavelength_regions == []
        ixs = collect(1:length(wavelength))
    else
        ixs = vcat([findall(bounds[1] .< wavelength .< bounds[2]) for bounds in wavelength_regions]...)
        sort!(ixs)
    end
    wavelength_cut = wavelength[ixs] =#
    nwvls = length(wavelength_bins) 

    kernel = gaussian_instrumental_kernel(wavelength, resolution)
    intensities_cut = Array{Vector{Float64}}(undef, length(temperatures), length(mus), length(vs))
    continuums_cut = Array{Vector{Float64}}(undef, length(temperatures), length(mus))
    for i=eachindex(temperatures)
        for j=eachindex(mus)
            for k=eachindex(vs)
                if resolution != 0
                    ## if I interpolate and shift in lnλ the shift will be linear!!
                    intensities_cut[i, j, k] = doppler_shift_spectrum(wavelength, instrumental_broadening(wavelength, intensities[i, j], resolution, kernel=kernel), vs[k])
                    Iitp = Spline1D(wavelength, intensities_cut[i, j, k], k=1, bc="nearest", s=0.0)
                    intensities_cut[i, j, k] = Iitp(wavelength_bins)
                else
                    intensities_cut[i, j, k] = doppler_shift_spectrum(wavelength, intensities[i, j], vs[k])
                    Iitp = Spline1D(wavelength, intensities_cut[i, j, k], k=1, bc="nearest", s=0.0)
                    intensities_cut[i, j, k] = Iitp(wavelength_bins)
                end
            end
            continuums_cut[i, j] = continuums[i, j]
            Icitp = Spline1D(wavelength, continuums_cut[i, j], k=1, bc="nearest", s=0.0)
            continuums_cut[i, j] = Icitp(wavelength_bins)
        end
    end
    intensities_itp = Interpolations.interpolate((temperatures, mus, vs), intensities_cut, Gridded(Linear()))
    continuums_itp = Interpolations.interpolate((temperatures, mus), continuums_cut, Gridded(Linear()))
    intensities_gradient_itp = (T, μ, v) -> Interpolations.gradient(intensities_itp, T, μ, v)[1]
    continuums_gradient_itp = (T, μ) -> Interpolations.gradient(continuums_itp, T, μ)[1]

    I = Array{Float64}(undef, nwvls)
    Ic = Array{Float64}(undef, nwvls)
    dIdT = Array{Float64}(undef, nwvls)
    dIcdT = Array{Float64}(undef, nwvls)
    dRdT = Array{Vector{Float64}}(undef, length(temperatures), length(mus), length(vs))
    for i=eachindex(temperatures)
        for j=eachindex(mus)
            for k=eachindex(vs)
                if mus[j] > 0.0
                    I .= intensities_itp(temperatures[i], mus[j], vs[k])
                    Ic .= continuums_itp(temperatures[i], mus[j])
                    dIdT .= intensities_gradient_itp(temperatures[i], mus[j], vs[k])
                    dIcdT .= continuums_gradient_itp(temperatures[i], mus[j])
                    dRdT[i, j, k] = (dIdT .* Ic .- I .* dIcdT) ./ Ic.^2
                else
                    dRdT[i, j, k] = zeros(nwvls)
                end
            end
        end
    end
    dRdT_itp = interpolate((temperatures, mus, vs), dRdT, Gridded(Linear()))

    return modelGrid(temperatures, logg, mus, nwvls, wavelength_bins, intensities_cut, continuums_cut, intensities_itp, continuums_itp, intensities_gradient_itp, continuums_gradient_itp, dRdT_itp)
end

function doppler_shift_spectrum(wavelength::Vector{Float64}, intensity::Vector{Float64}, v::Float64)
    ## Doppler shift a spectrum and interpolate to original wavelength axis
    # return LinearInterpolation(wavelength .* (1.0 + v/3E5), intensity, extrapolation_bc=Flat())(wavelength)
    return Spline1D(wavelength .* (1.0 + v/3E5), intensity, k=1, bc="nearest", s=0.0)(wavelength)
end

function calculate_global_profiles(temperature_map::Vector{Vector{Float64}}, stellar_geometry::Array{stellar_geometry, 1}, star_epoch_μ::Array{Float64, 2}, star_epoch_v::Array{Float64, 2}, stellar_parameters::Vector, polyflux::Matrix{Float64}, model_grid::modelGrid, vsini::Float64)
# function calculate_global_profiles(temperature_map, stellar_geometry, star_epoch_μ, star_epoch_v, stellar_parameters, polyflux, model_grid, vsini)
    global_intensities = Array{Vector{Float64}}(undef, length(stellar_geometry))
    # fill!(global_intensities, zeros(model_grid.nwvls))

    I = zeros(model_grid.nwvls)
    Ic = zeros(model_grid.nwvls)
    F = zeros(model_grid.nwvls)  # disk-integrated flux per pixel
    Fc = zeros(model_grid.nwvls)  # disk-integrated flux intensity per pixel
    
    @views for i = eachindex(stellar_geometry)  # loop over all phases
        if length(temperature_map) > 1  ## Use multiple-dispatch for this later on
            n = i
        else
            n = 1
        end

        for j in findall(star_epoch_μ[:, i] .> 0.0)  # loop over all visible pixels, i.e. where μ > 0
            I .= model_grid.intensities_itp(temperature_map[n][j], star_epoch_μ[j, i], star_epoch_v[j, i])# .* polyflux[i, j]
            Ic .= model_grid.continuums_itp(temperature_map[n][j], star_epoch_μ[j, i])# .* polyflux[i, j]  # Grab the model continuum intensity for the pixel's temperature and μ
            F .+= I .* polyflux[i, j]
            Fc .+= Ic .* polyflux[i, j]
        end
        # ixs = findall(star_epoch_μ[:, i] .> 0.0)
        # F .= sum(model_grid.intensities_itp.(temperature_map[ixs], star_epoch_μ[ixs, i], star_epoch_v[ixs, i]) .* polyflux[i, ixs])
        # Fc .= sum(model_grid.continuums_itp.(temperature_map[ixs], star_epoch_μ[ixs, i]) .* polyflux[i, ixs])
        
        @inbounds global_intensities[i] = F ./ Fc  # global profile is the total disk flux divided by the total disk intensity
        fill!(F, 0.0)
        fill!(Fc, 0.0)
    end

    return globalLineProfiles(length(temperature_map[1]), length(stellar_geometry), model_grid.nwvls, model_grid.wavelength, global_intensities)
end

function calculate_global_profiles_unnormalized(temperature_map::Vector{Vector{Float64}}, stellar_geometry::Array{stellar_geometry, 1}, star_epoch_μ::Array{Float64, 2}, star_epoch_v::Array{Float64, 2}, stellar_parameters::Vector, polyflux::Matrix{Float64}, model_grid::modelGrid, vsini::Float64)
# function calculate_global_profiles(temperature_map, stellar_geometry, star_epoch_μ, star_epoch_v, stellar_parameters, polyflux, model_grid, vsini)
    global_intensities = Array{Vector{Float64}}(undef, length(stellar_geometry))
    # fill!(global_intensities, zeros(model_grid.nwvls))

    I = zeros(model_grid.nwvls)
    F = zeros(model_grid.nwvls)  # disk-integrated flux per pixel
    
    @views for i = eachindex(stellar_geometry)  # loop over all phases
        if length(temperature_map) > 1  ## Use multiple-dispatch for this later on
            n = i
        else
            n = 1
        end

        for j in findall(star_epoch_μ[:, i] .> 0.0)  # loop over all visible pixels, i.e. where μ > 0
            I .= model_grid.intensities_itp(temperature_map[n][j], star_epoch_μ[j, i], star_epoch_v[j, i])# .* polyflux[i, j]
            F .+= I .* polyflux[i, j]
        end
        # ixs = findall(star_epoch_μ[:, i] .> 0.0)
        # F .= sum(model_grid.intensities_itp.(temperature_map[ixs], star_epoch_μ[ixs, i], star_epoch_v[ixs, i]) .* polyflux[i, ixs])
        # Fc .= sum(model_grid.continuums_itp.(temperature_map[ixs], star_epoch_μ[ixs, i]) .* polyflux[i, ixs])

        @inbounds global_intensities[i] = deepcopy(F)  # global profile is the total disk flux divided by the total disk intensity
        fill!(F, 0.0)
    end

    return globalLineProfiles(length(temperature_map[1]), length(stellar_geometry), model_grid.nwvls, model_grid.wavelength, global_intensities)
end

function calculate_μv(npix::Int64, stellar_geometry::Array{stellar_geometry, 1}, stellar_parameters::Vector, vsini::Float64)
    nepochs = length(stellar_geometry)
    μ = zeros(npix, nepochs)  # μ angle per pixel
    xcosPA = zeros(npix)
    ysinPA = zeros(npix)
    v = zeros(npix, nepochs)  # velocity per pixel

    PA = -stellar_parameters[1].position_angle
    cosPA = cos(PA * pi/180)
    sinPA = sin(PA * pi/180)

    @views for i=1:nepochs
        μ[:, i] .= stellar_geometry[i].normals[:, 3]  # μ angle is just the z-coordinate of each pixel
        xcosPA .= stellar_geometry[i].normals[:, 1] .* cosPA
        ysinPA .= stellar_geometry[i].normals[:, 2] .* sinPA
        v[:, i] .= vsini .* (xcosPA .- ysinPA) # velocity of each pixel is vsini*(xcos(PA) - ysin(PA)), from rotation about z-axis
    end

    return μ, v
end

# function calculate_global_profiles(temperature_map::Vector{Vector{Float64}}, stellar_geometry::Array{stellar_geometry}, stellar_parameters::Vector, polyflux::Matrix{Float64}, model_grid::modelGrid, vsini::Float64)
#= function calculate_global_profiles(temperature_map::Vector{Vector{Float64}}, stellar_geometry, stellar_parameters, polyflux, model_grid, vsini)
    global_intensities = Array{Vector{Float64}}(undef, length(stellar_geometry))

    μ = Array{Float64}(undef, length(temperature_map[1]))  # μ angle per pixel
    xcosPA = Array{Float64}(undef, length(temperature_map[1]))
    ysinPA = Array{Float64}(undef, length(temperature_map[1]))
    v = Array{Float64}(undef, length(temperature_map[1]))  # velocity per pixel

    PA = -stellar_parameters[1].position_angle
    cosPA = cos(PA * pi/180)
    sinPA = sin(PA * pi/180)

    F = zeros(model_grid.nwvls)  # disk-integrated flux per pixel
    Fc = zeros(model_grid.nwvls)  # disk-integrated flux intensity per pixel
    I = Array{Float64}(undef, model_grid.nwvls)  # specific intensity per pixel
    Ic = Array{Float64}(undef, model_grid.nwvls)  # continuum specific intensity per pixel
    @views for i=1:length(stellar_geometry)  # loop over all phases
        μ .= stellar_geometry[i].normals[:, 3]  # μ angle is just the z-coordinate of each pixel
        xcosPA .= stellar_geometry[i].normals[:, 1] .* cosPA
        ysinPA .= stellar_geometry[i].normals[:, 2] .* sinPA
        v .= vsini .* (xcosPA .- ysinPA) # velocity of each pixel is vsini*(xcos(PA) - ysin(PA)), from rotation about z-axis
        
        for j in findall(μ .> 0.0)  # loop over all visible pixels, i.e. where μ > 0
            I .= model_grid.intensities_itp(temperature_map[i][j], μ[j], v[j]) .* polyflux[i, j]  # Grab the model intensity for the pixel's temperature and μ
            Ic .= model_grid.continuums_itp(temperature_map[i][j], μ[j]) .* polyflux[i, j]  # Grab the model continuum intensity for the pixel's temperature and μ
            F .+= I  # non-normalized flux is the sum of all pixel fluxes
            Fc .+= Ic  # total continuum is the sum of all pixel continuum fluxes
        end
        global_intensities[i] = F ./ Fc  # global profile is the total disk flux divided by the total disk intensity
        fill!(F, 0.0)
        fill!(Fc, 0.0)
    end

    return globalLineProfiles(length(temperature_map[1]), length(stellar_geometry), model_grid.nwvls, model_grid.wavelength, global_intensities)
end =#

function bin_spectrum(wavelength, intensity, bins)
    ixs = searchsortedfirst.(Ref(bins), wavelength)
    intensity_binned = vec([mean(intensity[ixs .== i]) for i=1:length(bins)])
    return intensity_binned
end

function read_di_profiles(files; P=1.0, T0=0.0, wavelength_bins=[], verb=true)
    ## Setup data structure and pull relevant keywords from fits header
    nepochs = length(files)
    intensities = Array{Vector{Float64}}(undef, nepochs)
    mjd = Array{Float64}(undef, nepochs)
    snr = Array{Float64}(undef, nepochs)
    res = Array{Float64}(undef, nepochs)
    phase = Array{Float64}(undef, nepochs)

    @views for i=1:nepochs
        open(files[i], "r") do f
            data = readdlm(files[i], '\t', header=true)
            wavelength = data[1][:, 1]
            intensity = data[1][:, 2]
            if length(wavelength_bins) > length(wavelength)
                println("WARNING: Upsampling wavelength grid...")
                f = Spline1D(wavelength, intensity, k=1, s=0.0, bc="nearest")
                intensities[i] = f(wavelength_bins)
            else
                ## bin the spectrum
                intensities[i] = bin_spectrum(wavelength, intensity, wavelength_bins)
            end

            ## cut the wavelengths needed
            # ixs = findall(wavelength_bins[1] .<= wavelength .<= wavelength_bins[end])
            # f = Spline1D(wavelength[ixs], intensity[ixs], k=1, s=0.0, bc="nearest")
            # wavelength = wavelength[ixs]
            # Wintensity = f(wavelength)
            
            mjd[i] = parse(Float64, data[2][2])
            snr[i] = parse(Float64, data[2][4])
            res[i] = parse(Float64, data[2][6])
        end
    end

    phase = jd2phase(mjd, P, T0)
    nwvls = length(wavelength_bins)
    observed_profiles = observedLineProfiles(files, nepochs, nwvls, snr, mjd, phase, wavelength_bins, intensities)
    if verb == true
        display(observed_profiles)
    end
    return observed_profiles
end

function jd2phase(jd::Vector{Float64}, P::Float64, T0::Float64)
    ## Compute phase from time of reference (e.g. time of minimum, periastron passage, etc.) and period
    return (jd .- T0) ./ P .- floor.((jd .- T0) ./ P)
end

function setup_stellar_parameters_di(starparams::Vector{Any}, didata::observedLineProfiles; mjd_long0 = 0.0)
    stellar_parameters = Vector(undef, didata.nepochs);
    for i=1:didata.nepochs
        stellar_parameters[i]=starparameters(starparams[1],starparams[2],starparams[3],starparams[4],starparams[5],starparams[6],starparams[7],starparams[8],mod(360.0/starparams[10].*(didata.mjd[i]-mjd_long0),360),starparams[10]);
    end
    return stellar_parameters
end

function di_reconstruct(x_start::Vector{Float64}, didata::observedLineProfiles, polyflux::Matrix{Float64}, star_epoch_geom::Array{stellar_geometry}, star_epoch_μ::Array{Float64, 2}, star_epoch_v::Array{Float64, 2}, stellar_parameters::Vector, model_grid::modelGrid, vsini::Float64, visible_pixels::Vector{Int64}; lowtemp = 0, hightemp=10000, relative = false, printcolor= [], verb = true, maxiter = 200, regularizers =[], wavelength_regions=[], normalization_region=[])
    ## Use VMLMB from OptimPackNextGen to find the temperature map that minimizes the χ^2 between the observed profiles and the models
    x_sol = [];
    crit_imaging = x -> x; #TODO: proper dummy init
    crit_imaging = (x, g) -> di_crit_fg(x, g, didata, polyflux, star_epoch_geom,  star_epoch_μ, star_epoch_v, stellar_parameters, model_grid, vsini, visible_pixels, regularizers=regularizers, verb = verb, wavelength_regions=wavelength_regions, normalization_region=normalization_region);
    x_sol = OptimPackNextGen.vmlmb(crit_imaging, x_start, verb=false, lower=lowtemp, upper=hightemp, maxiter=maxiter, blmvm=false, gtol=(0, 1e-8), maxeval=100000);#, fmin=0.99*length(x_start));
    return x_sol
end

function di_crit_fg(x::Vector{Float64}, g::Vector{Float64}, didata::observedLineProfiles, polyflux::Matrix{Float64}, star_epoch_geom::Array{stellar_geometry, 1}, star_epoch_μ::Array{Float64, 2}, star_epoch_v::Array{Float64, 2}, stellar_parameters::Vector, model_grid::modelGrid, vsini::Float64, visible_pixels::Vector{Int64}; regularizers=[], verb = true, wavelength_regions=[], normalization_region=[])
# function di_crit_fg(x, g, didata, polyflux, star_epoch_geom, star_epoch_μ, star_epoch_v, stellar_parameters, model_grid, vsini, visible_pixels; regularizers=[], verb = true, wavelength_regions=[])
## Compute χ^2 and ∇χ^2
    # x .= max.(x, 3500.0)
    # x .= min.(x, 5500.0)

    if wavelength_regions == []  # Specific spectral lines are selected here
        ixs = collect(1:length(model_grid.wavelength))
    else
        ixs = vcat([findall(bounds[1] .< model_grid.wavelength .< bounds[2]) for bounds in wavelength_regions]...)
        sort!(ixs)
    end

    global_line_profiles::globalLineProfiles = calculate_global_profiles([x], star_epoch_geom, star_epoch_μ, star_epoch_v, stellar_parameters, polyflux, model_grid, vsini)
    # global_line_profiles::globalLineProfiles = calculate_global_profiles_unnormalized([x], star_epoch_geom, star_epoch_μ, star_epoch_v, stellar_parameters, polyflux, model_grid, vsini)
    res_di = Array{Vector{Float64}, 1}(undef, length(star_epoch_geom))
    scale = Vector{Float64}(undef, didata.nepochs)
    @views for i=1:length(star_epoch_geom)
        @inbounds global_line_profiles.intensities[i] .= rescale(didata.wavelength, global_line_profiles.intensities[i], didata.intensities[i], normalization_region)
        # @inbounds didata.intensities[i] ./= fit_continuum(didata.wavelength, didata.intensities[i], Nσ_low=0.5, porder=1)
        # scale[i] = mean(didata.intensities[i][ixs]) / mean(global_line_profiles.intensities[i][ixs])
        # @inbounds global_line_profiles.intensities[i] .*= scale[i]
        @inbounds res_di[i] = didata.intensities[i][ixs] .- global_line_profiles.intensities[i][ixs]
    end

    chi2_f_di::Float64 = calculate_chi2_f(res_di, didata.snr)  # Calculate χ^2, weighted by SNR
    g_di::Vector{Float64} = zeros(length(g))
    g_di .= calculate_chi2_g(x, g_di, star_epoch_geom, star_epoch_μ, star_epoch_v, stellar_parameters, model_grid, polyflux, vsini, res_di, didata.snr, wavelength_regions=wavelength_regions)  # calculate ∇χ^2
    # g_di .= calculate_chi2_g_unnormalized(x, g_di, star_epoch_geom, star_epoch_μ, star_epoch_v, stellar_parameters, model_grid, polyflux, vsini, res_di, didata.snr, wavelength_regions=wavelength_regions, scale=scale)  # calculate ∇χ^2

    # g_numerical::Vector{Float64} = zeros(length(g))
    # g_di .= calculate_chi2_g_numerical(x, g_di, star_epoch_geom, star_epoch_μ, star_epoch_v, stellar_parameters, model_grid, polyflux, vsini, res_di, didata.snr, wavelength_regions=wavelength_regions)  # calculate ∇χ^2

    reg_f::Float64 = spheroid_regularization(x, g, regularizers=regularizers, verb = false);  #Apply regularizers; note: this adds to g
    chi2_f::Float64 = chi2_f_di
    g .= g_di
    if verb == true
        y = x[visible_pixels]
        low_x = minimum(y);
        high_x = maximum(y);
        mean_x = mean(y);
        @printf("ftot: %8.2f\tχ^2: %8.2f\tχ^2r: %4.2f\tμ*freg:%8.2f\tminT: %5.1f\tmaxT: %5.1f\tmeanT: %5.1f\n", chi2_f+reg_f, chi2_f_di, chi2_f_di/length(ixs)/length(res_di), reg_f, low_x, high_x, mean_x);
    end

    return chi2_f + reg_f
end

function rescale(wavelength, global_line_profile, didata_intensity, normalization_region)
    if normalization_region == []
        ixs = Int.(collect(range(1, stop=length(wavelength), length=length(wavelength))))
    else
        ixs = findall(normalization_region[1] .<= wavelength .<= normalization_region[2])
    end
    return global_line_profile .* mean(didata_intensity[ixs]) / mean(global_line_profile[ixs])  # renormalize
end

function calculate_chi2_f(res::Vector{Vector{Float64}}, SNR::Vector{Float64})
    chi2_f::Float64 = 0.0
    for i=eachindex(res)  # Here σ = 1/SNR
        chi2_f += sum(res[i].^2) / (1/SNR[i]^2)
    end

    return chi2_f
end

function calculate_chi2_g(x::Vector{Float64}, g::Vector{Float64}, star_epoch_geom::Vector{stellar_geometry}, star_epoch_μ::Array{Float64, 2}, star_epoch_v::Array{Float64, 2}, stellar_parameters::Vector, model_grid::modelGrid, polyflux::Matrix{Float64}, vsini::Float64, res::Vector{Vector{Float64}}, SNR::Vector{Float64}; wavelength_regions=[])
# function calculate_chi2_g(x, g, star_epoch_geom, star_epoch_μ, star_epoch_v, stellar_parameters, model_grid, polyflux, vsini, res, SNR; wavelength_regions=[])
    if wavelength_regions == [] # Specific spectral lines are selected here
        ixs = collect(1:length(model_grid.wavelength))
    else
        ixs = vcat([findall(bounds[1] .< model_grid.wavelength .< bounds[2]) for bounds in wavelength_regions]...)
        sort!(ixs)
    end

    @views for i=eachindex(star_epoch_geom)
        for j in findall(star_epoch_μ[:, i] .> 0.0)
            @inbounds g[j] += -2.0 * polyflux[i, j] * sum(model_grid.dRdT(x[j], star_epoch_μ[j, i], star_epoch_v[j, i])[ixs] .* res[i]) / (1/SNR[i])^2
        end
    end

    return g
end

function calculate_chi2_g_unnormalized(x::Vector{Float64}, g::Vector{Float64}, star_epoch_geom::Vector{stellar_geometry}, star_epoch_μ::Array{Float64, 2}, star_epoch_v::Array{Float64, 2}, stellar_parameters::Vector, model_grid::modelGrid, polyflux::Matrix{Float64}, vsini::Float64, res::Vector{Vector{Float64}}, SNR::Vector{Float64}; wavelength_regions=[], scale=1.0)
# function calculate_chi2_g(x, g, star_epoch_geom, star_epoch_μ, star_epoch_v, stellar_parameters, model_grid, polyflux, vsini, res, SNR; wavelength_regions=[])
    if wavelength_regions == [] # Specific spectral lines are selected here
        ixs = collect(1:length(model_grid.wavelength))
    else
        ixs = vcat([findall(bounds[1] .< model_grid.wavelength .< bounds[2]) for bounds in wavelength_regions]...)
        sort!(ixs)
    end

    @views for i=eachindex(star_epoch_geom)
        for j in findall(star_epoch_μ[:, i] .> 0.0)
            @inbounds g[j] += -2.0 * polyflux[i, j] * scale[j] * sum(model_grid.intensities_gradient_itp(x[j], star_epoch_μ[j, i], star_epoch_v[j, i])[ixs] .* res[i]) / (1/SNR[i])^2
        end
    end

    return g
end

function calculate_chi2_g_numerical(x::Vector{Float64}, g::Vector{Float64}, star_epoch_geom::Vector{stellar_geometry}, star_epoch_μ::Array{Float64, 2}, star_epoch_v::Array{Float64, 2}, stellar_parameters::Vector, model_grid::modelGrid, polyflux::Matrix{Float64}, vsini::Float64, res::Vector{Vector{Float64}}, SNR::Vector{Float64}; wavelength_regions=[])
    # function calculate_chi2_g(x, g, star_epoch_geom, star_epoch_μ, star_epoch_v, stellar_parameters, model_grid, polyflux, vsini, res, SNR; wavelength_regions=[])
    if wavelength_regions == [] # Specific spectral lines are selected here
        ixs = collect(1:length(model_grid.wavelength))
    else
        ixs = vcat([findall(bounds[1] .< model_grid.wavelength .< bounds[2]) for bounds in wavelength_regions]...)
        sort!(ixs)
    end

    ΔT = 10.0
    I = zeros(model_grid.nwvls)
    I_plus_ΔT = zeros(model_grid.nwvls)
    Ic = zeros(model_grid.nwvls)
    Ic_plus_ΔT = zeros(model_grid.nwvls)
    @views for i=eachindex(star_epoch_geom)
        for j in findall(star_epoch_μ[:, i] .> 0.0)
            I .= model_grid.intensities_itp(x[j], star_epoch_μ[j, i], star_epoch_v[j, i])
            I_plus_ΔT .= model_grid.intensities_itp(x[j].+ΔT, star_epoch_μ[j, i], star_epoch_v[j, i])
            Ic .= model_grid.continuums_itp(x[j], star_epoch_μ[j, i])
            Ic_plus_ΔT .= model_grid.continuums_itp(x[j].+ΔT, star_epoch_μ[j, i])       
            @inbounds g[j] += -2.0 * polyflux[i, j] * sum((((I_plus_ΔT./Ic_plus_ΔT) .- (I./Ic)) ./ ΔT)[ixs] .* res[i]) / (1/SNR[i])^2
        end
    end

    return g
end

function plot_profiles(line_profiles; didata=[], line_centers=[], line_labels=[], wavelength_regions=[], vsini=0.0, grow=1.0)
    ## Plots the observed and model profiles
    ## Can either plot the whole profile or just the relevant lines
    
    if line_centers != []
        fig, ax = PyPlot.subplots(1, length(line_centers), figsize=(14, 8), sharey=true, gridspec_kw=Dict("wspace" => 0.05))
        # subplots_adjust(wspace=0.0, hspace=0.0)
        # fig.text(0.5, 0.04, "Velocity [km/s]", ha="center", fontweight="bold", fontsize=15)
        # fig.text(0.04, 0.5, "Normalized Flux", va="center", rotation="vertical", fontweight="bold", fontsize=15)
        ax[1].set_ylabel("Normalized Flux", fontweight="bold", fontsize=15)
        if (length(ax) > 1) && (length(ax) % 2 != 0)
            mid = Int(median(collect(1:length(ax))))
            ax[mid].set_xlabel(raw"Velocity [km s$^{-1}$]", fontweight="bold", fontsize=15)
        else
            fig.text(0.5, 0.04, raw"Velocity [km s$^{-1}$]", ha="center", fontweight="bold", fontsize=15)
        end

        if length(line_centers) == 1; ax = [ax]; end
        for i=eachindex(line_centers)
            # wvlmin = line_centers[i] * (1 - 1.25*vsini / 3E5)
            # wvlmax = line_centers[i] * (1 + 1.25*vsini / 3E5)
            wvlmin = wavelength_regions[i][1]
            wvlmax = wavelength_regions[i][2]
            vmin = (wvlmin .- line_centers[i]) ./ line_centers[i] .* 3E5
            vmax = (wvlmax .- line_centers[i]) ./ line_centers[i] .* 3E5
            ax[i].set_title(line_labels[i], fontweight="bold", fontsize=15)
            for j=1:line_profiles.nepochs
                if j == line_profiles.nepochs
                    mlabel = "Inferred"
                    dlabel = "Observed"
                else
                    mlabel = dlabel = "_nolabel_"
                end
                v = (line_profiles.wavelength .- line_centers[i]) ./ line_centers[i] .* 3E5
                ax[i].plot(v, line_profiles.intensities[j].-0.05*(j-1), "C1--", lw=1.5, label=mlabel, zorder=10)
                if didata != []
                    ax[i].plot(v, didata.intensities[j].-0.05*(j-1), "C0-", lw=1.5, alpha=0.75, label=dlabel)
                end
            end
            ax[i].set_xlim((vmin, vmax))
            ax[i].set_xticks(range(-vsini, stop=vsini, step=vsini/2))
            ax[i].set_xticklabels(range(-vsini, stop=vsini, step=vsini/2), fontsize=8)
        end
        # ax[3].legend(loc="center left", bbox_to_anchor=(1, 0.5), markerscale=3, fontsize=14, frameon=false)

        for j=1:line_profiles.nepochs
            vels = 3e5 .* ((didata.wavelength .- line_centers[end]) ./ line_centers[end])
            ixs = findall(-vsini*grow .<= vels .<= vsini*grow)
            xdata = vels[ixs]; ydata = didata.intensities[j][ixs]
            yphi = mean(ydata[end-5:end])
            phi = @sprintf("%1.3f", didata.phase[j])
            ax[end].annotate(string(" - ", phi), xy=(xdata[end], yphi - 0.05*(j-1)), ha="center", va="center", xytext=(20, 0), textcoords="offset points", fontsize=12, annotation_clip=false)
            if j == 1
                ax[end].annotate(raw"  $\phi$", xy=(xdata[end], yphi - 0.05*(j-1)), ha="center", va="center", xytext=(20, 15), textcoords="offset points", fontsize=12)
            end
        end
        ax1_x1coord = ax[1].get_position().x0
        axend_x2coord = ax[end].get_position().x1
        xmiddle = (axend_x2coord + ax1_x1coord) / 2
        ybottom = ax[1].get_position().y0
        # ax[end].legend(loc="upper center", bbox_to_anchor=(xdata[end], yphi - 0.05*j), bbox_transform=ax[end].transData, fontsize=14)
        ax[end].legend(loc="lower center", bbox_to_anchor=(xmiddle, ybottom+0.01), fontsize=14, bbox_transform=fig.transFigure, framealpha=0.95)

    else
        fig, ax = PyPlot.subplots(1, 1, figsize=(14, 8))
        # subplots_adjust(wspace=0.0, hspace=0.0)
        ax.set_xlabel(raw"Velocity [km s$^{-1}$]", fontweight="bold", fontsize=15)
        ax.set_ylabel("Normalized Flux", fontweight="bold", fontsize=15)
        for j=1:line_profiles.nepochs
            xmodel = line_profiles.wavelength; ymodel = line_profiles.intensities[j]
            xdata = didata.wavelength; ydata = didata.intensities[j]
            yphi = mean(sort(ydata)[end-(length(ydata)÷2):end])
            if j == 1
                mlabel = "Inferred"
                dlabel = "Observed"
            else
                mlabel = "_nolabel_"
                dlabel = "_nolabel_"
            end
            ax.plot(xmodel, ymodel .-0.05*(j-1), "C1--", lw=1.5, zorder=10, label=mlabel)
            if didata != []
                ax.plot(xdata, ydata .-0.05*(j-1), "C0-", lw=1.5, label=dlabel)
                ax.annotate(string("- ", round(didata.phase[j], digits=3)), xy=(xdata[end], yphi - 0.05*(j-1)), ha="center", va="center", xytext=(50, 0), textcoords="offset pixels", fontsize=12)
                if j == 1
                    ax.annotate(raw"$\phi$", xy=(xdata[end], yphi - 0.05*(j-1)), ha="center", va="center", xytext=(50, 75), textcoords="offset pixels", fontsize=12)
                end
            end
        end
        ax.legend(fontsize=14, markerscale=3)
    end
    PyPlot.draw()
end

function di_reconstruct_spotfill(x_start::Vector{Float64}, didata::observedLineProfiles, polyflux::Matrix{Float64}, star_epoch_geom::Array{stellar_geometry}, stellar_parameters::Vector, model_grid::modelGrid, vsini::Float64, visible_pixels::Vector{Int64}, Tphot::Float64, Tspot::Float64; printcolor= [], verb = true, maxiter = 200, regularizers =[], wavelength_regions=[])
    ## Use VMLMB from OptimPackNextGen to find the temperature map that minimizes the χ^2 between the observed profiles and the models, where the fitted parameter is the spot-filling fctor per pixel
    x_sol = [];
    crit_imaging = x -> x; #TODO: proper dummy init
    crit_imaging = (x, g) -> di_crit_fg_spotfill(x, g, didata, polyflux, star_epoch_geom, stellar_parameters, model_grid, vsini, visible_pixels, Tphot, Tspot, regularizers=regularizers, verb = verb, wavelength_regions=wavelength_regions);
    x_sol = OptimPackNextGen.vmlmb(crit_imaging, x_start, verb=false, lower=0.0, upper=1.0, maxiter=maxiter, blmvm=false, gtol=(0, 1e-8));#, fmin=0.99*length(x_start));
    return x_sol
end

function di_crit_fg_spotfill(x::Vector{Float64}, g::Vector{Float64}, didata::observedLineProfiles, polyflux::Matrix{Float64}, star_epoch_geom::Array{stellar_geometry}, stellar_parameters::Vector, model_grid::modelGrid, vsini::Float64, visible_pixels::Vector{Int64}, Tphot::Float64, Tspot::Float64; regularizers=[], verb = true, wavelength_regions=[])
    ## Compute χ^2 and ∇χ^2
    if wavelength_regions == []  # Specific spectral lines are selected here
        ixs = collect(1:length(model_grid.wavelength))
    else
        ixs = vcat([findall(bounds[1] .< model_grid.wavelength .< bounds[2]) for bounds in wavelength_regions]...)
        sort!(ixs)
    end

    global_line_profiles = calculate_global_profiles_spotfill(x, star_epoch_geom, stellar_parameters, polyflux, model_grid, vsini, Tphot, Tspot)
    res = Array{Vector{Float64}}(undef, length(star_epoch_geom))
    for i=1:length(star_epoch_geom)
        global_line_profiles.intensities[i] ./= (mean(global_line_profiles.intensities[i][ixs]) / mean(didata.intensities[i][ixs]))  # renormalize
        res[i] = didata.intensities[i][ixs] .- global_line_profiles.intensities[i][ixs]
    end
    chi2_f = calculate_chi2_f(res, didata.snr)  # Calculate χ^2, weighted by SNR

    g .= zeros(length(g))
    g .= calculate_chi2_g_spotfill(x, g, star_epoch_geom, stellar_parameters, model_grid, polyflux, vsini, res, didata.snr, Tphot, Tspot, wavelength_regions=wavelength_regions)  # calculate ∇χ^2

    reg_f = spheroid_regularization(x, g, regularizers=regularizers, verb = false);  #Apply regularizers; note: this adds to g
    if verb == true
        y = x[visible_pixels]
        low_x = minimum(y);
        high_x = maximum(y);
        mean_x = mean(y);
        @printf("ftot: %8.2f\tχ^2: %8.2f\tχ^2r: %4.2f\tμ*freg: %8.2f\t∑f/Npix: %2.4f\n", chi2_f+reg_f, chi2_f, chi2_f/length(ixs)/length(res), reg_f, sum(x)/length(x));
    end

    return chi2_f + reg_f
end

function calculate_global_profiles_spotfill(f_map::Vector{Float64}, stellar_geometry::Array{stellar_geometry}, stellar_parameters::Vector, polyflux::Matrix{Float64}, model_grid::modelGrid, vsini::Float64, Tphot::Float64, Tspot::Float64)
    # local_line_profiles = calculate_local_profiles(temperature_map, stellar_geometry, model_grid, vsini)
    global_intensities = Array{Vector{Float64}}(undef, length(stellar_geometry))
    fill!(global_intensities, zeros(model_grid.nwvls))

    μ = Array{Float64}(undef, length(f_map))  # μ angle per pixel
    xcosPA = Array{Float64}(undef, length(f_map))  # μ angle per pixel
    ysinPA = Array{Float64}(undef, length(f_map))  # μ angle per pixel
    v = Array{Float64}(undef, length(f_map))  # velocity per pixel

    PA = -stellar_parameters[1].position_angle
    cosPA = cos(PA * pi/180)
    sinPA = sin(PA * pi/180)

    intensity = zeros(model_grid.nwvls)
    cintensity = zeros(model_grid.nwvls)
    F = zeros(model_grid.nwvls)  # disk-integrated flux per pixel
    Fc = zeros(model_grid.nwvls)  # disk-integrated flux intensity per pixel
    @views for i=1:length(stellar_geometry)  # loop over all phases
        μ .= stellar_geometry[i].normals[:, 3]  # μ angle is just the z-coordinate of each pixel
        xcosPA .= stellar_geometry[i].normals[:, 1] .* cosPA
        ysinPA .= stellar_geometry[i].normals[:, 2] .* sinPA
        v .= vsini .* (xcosPA .- ysinPA)
        for j in findall(μ .> 0.0)  # loop over all visible pixels, i.e. where μ > 0
            intensity .= f_map[j] .* model_grid.intensities_itp(Tspot, μ[j], v[j]) .+ (1 - f_map[j]) .* model_grid.intensities_itp(Tphot, μ[j], v[j])
            # intensity .= doppler_shift_spectrum(model_grid.wavelength, intensity, v[j])
            F .+= intensity .* polyflux[i, j]
            cintensity .= (f_map[j] .* model_grid.continuums_itp(Tspot, μ[j]) .+ (1 - f_map[j]) .* model_grid.continuums_itp(Tphot, μ[j]))
            Fc .+= cintensity .* polyflux[i, j]
        end
        global_intensities[i] = F ./ Fc  # global profile is the total disk flux divided by the total disk continuum flux
        fill!(F, 0.0)
        fill!(Fc, 0.0)
    end

    return globalLineProfiles(length(f_map), length(stellar_geometry), model_grid.nwvls, model_grid.wavelength, global_intensities)
end

function calculate_chi2_g_spotfill(x::Vector{Float64}, g::Vector{Float64}, star_epoch_geom::Vector{stellar_geometry}, stellar_parameters::Vector, model_grid::modelGrid, polyflux::Matrix{Float64}, vsini::Float64, res::Vector{Vector{Float64}}, SNR::Vector{Float64}, Tphot::Float64, Tspot::Float64; wavelength_regions=[])
    if wavelength_regions == [] # Specific spectral lines are selected here
        ixs = collect(1:length(model_grid.wavelength))
    else
        ixs = vcat([findall(bounds[1] .< model_grid.wavelength .< bounds[2]) for bounds in wavelength_regions]...)
        sort!(ixs)
    end

    μ = Array{Float64}(undef, length(x))  # μ angle per pixel
    xcosPA = Array{Float64}(undef, length(x))  # μ angle per pixel
    ysinPA = Array{Float64}(undef, length(x))  # μ angle per pixel
    v = Array{Float64}(undef, length(x))  # velocity per pixel

    PA = -stellar_parameters[1].position_angle
    cosPA = cos(PA * pi/180)
    sinPA = sin(PA * pi/180)

    F = zeros(model_grid.nwvls)  # disk-integrated flux per pixel
    Fc = zeros(model_grid.nwvls)  # disk-integrated flux intensity per pixel
    Iphot = zeros(model_grid.nwvls)  # photosphere specific intensity per pixel
    Ispot = zeros(model_grid.nwvls)  # spot specific intensity per pixel
    Icphot = zeros(model_grid.nwvls)  # photosphere continuum intensity per pixel
    Icspot = zeros(model_grid.nwvls)  # spot continuum intensity per pixel
    dRdf = zeros(model_grid.nwvls)  # spot continuum intensity per pixel
    @views for i=1:length(star_epoch_geom)
        μ .= star_epoch_geom[i].normals[:, 3]
        xcosPA .= star_epoch_geom[i].normals[:, 1] .* cosPA
        ysinPA .= star_epoch_geom[i].normals[:, 2] .* sinPA
        v .= vsini .* (xcosPA .- ysinPA)
        for j in findall(μ .> 0.0)
            Ispot .= model_grid.intensities_itp(Tspot, μ[j], v[j])
            Iphot .= model_grid.intensities_itp(Tphot, μ[j], v[j])
            Icspot .= model_grid.continuums_itp(Tspot, μ[j])
            Icphot .= model_grid.continuums_itp(Tphot, μ[j])
            F .= x[j].*Ispot .+ (1-x[j]).*Iphot
            Fc .= x[j].*Icspot .+ (1-x[j]).*Icphot
            dRdf .= (Fc.*(Ispot .- Iphot) - F.*(Icspot .- Icphot)) ./ (Fc .* Fc)
            # dRdf .= doppler_shift_spectrum(model_grid.wavelength, dRdf, v[j])
            @inbounds g[j] += -2.0 * polyflux[i, j] * sum(dRdf[ixs] .* res[i]) / (1/SNR[i])^2  # Calculate ∇χ^2
        end
    end

    return g
end

function make_circ_spot_DI_spotfill(star_geometry,spot_radius,lat,long; profile="flat")
    ## Makes a circular spot at a given (lat, long)
    ##  spot_radius ∈ [0, 360]
    ##  lat ∈ [-90, 90]
    ##  long ∈ [-180, 180]
    ##  bright_frac ∈ [0, ∞)

    long = mod(long, 360)

    x = zeros(star_geometry.npix)
    centers = star_geometry.vertices_spherical[:,5,:]
    θ = centers[:, 3]
    ϕ = -centers[:, 2] .+ pi/2

    da = 2.0*asin.(sqrt.(sin.(0.5*(ϕ .- lat*pi/180)).^2+cos(lat*pi/180).*cos.(ϕ).*sin.(0.5*(θ .- long*pi/180)).^2))

    spot_mask = findall(da .<= spot_radius*pi/180)
    if profile == "flat"
        x[spot_mask] .= 1.0
    elseif profile == "linear"
        x[spot_mask] .= 1.0 .- da[spot_mask] ./ (spot_radius*pi/180)
    end

    return vec(x)
end

function mollplot_healpix(image; visible_pixels = [], vmin = -Inf, vmax = Inf, incl=90.0, colormap="gist_heat", figtitle="Mollweide", clabel="Temperature [K]")
    xsize = 2000
    ysize = div(xsize,2)
    theta = collect(range(pi, stop=0.0, length=ysize))
    phi   = collect(range(-pi, stop=pi, length=xsize))
    longitude = collect(range(-180, stop=180, length=xsize))/180*pi
    latitude = collect(range(-90, stop=90, length=ysize))/180*pi
    # project the map to a rectangular matrix xsize x ysize
    nside = npix2nside(length(image))
    PHI = [i for j in theta, i in phi]
    THETA = [j for j in theta, i in phi]
    grid_pix = reshape(ang2pix_nest(nside, vec(THETA), vec(PHI)), size(PHI))
    grid_map = image[grid_pix]
    fig = figure(figtitle, figsize=(10, 7))
    clf();
    ax = subplot(111,projection="mollweide")
    # rasterized makes the map bitmap while the labels remain vectorial
    # flip longitude to the astro convention
    if visible_pixels == []
      if (vmin == -Inf)
        vmin = minimum(image);
      end
      if (vmax == Inf)
        vmax = maximum(image);
      end
    else
      if (vmin == -Inf)
        vmin = minimum(image[visible_pixels]);
      end
      if (vmax == Inf)
        vmax = maximum(image[visible_pixels]);
      end
    end
    moll = pcolormesh(longitude, latitude, grid_map, vmin=vmin, vmax=vmax, rasterized=true, cmap=colormap)
    # graticule
    ax.set_longitude_grid(30)
    ax.set_latitude_grid(30)
    ax.set_longitude_grid_ends(90)
    spacing = 0.04
    subplots_adjust(bottom=spacing, top=1-spacing, left=spacing, right=1-spacing)
    grid(true)
    if incl != 90.0
        ax.axhline(-incl * pi/180, c=:black, ls="-.")
    end
  # class ThetaFormatterShiftPi(GeoAxes.ThetaFormatter):
  #     """Shifts labelling by pi
  #     Shifts labelling from -180,180 to 0-360"""
  #     def __call__(self, x, pos=None):
  #         if x != 0:
  #             x *= -1
  #         if x < 0:
  #             x += 2*np.pi
  #         return GeoAxes.ThetaFormatter.__call__(self, x, pos)
  #
  #
  # ax[:xaxis][:set_major_formatter](ThetaFormatterShiftPi(60))
  # colorbar
  if vmin != vmax
    ticks = collect(range(vmin, stop=vmax, length=7));
    cb = colorbar(moll, orientation="horizontal", shrink=.6, pad=0.05, ticks=ticks)
    cb.ax.xaxis.labelpad=5
    cb.ax.xaxis.set_label_text(clabel)

    # workaround for issue with viewers, see colorbar docstring
    cb.solids.set_edgecolor("face")
  end
  ax.tick_params(axis="x", labelsize=15)
  ax.tick_params(axis="y", labelsize=15)
end

function mollplot_longlat(image, ntheta, nphi; visible_pixels = [], vmin = -Inf, vmax = Inf, colormap="gist_heat", figtitle="Mollweide", clabel="Temperature [K]")
    xsize = 2000
    ysize = div(xsize,2)

    ## DJ: had to set xsize=ntheta and ysize=nphi to get plotting to work for large
    ##  number of latitudes and longitudes
    # xsize = ntheta
    # ysize = nphi

    theta = collect(range(pi, stop=0, length=ysize))
    phi   = collect(range(-pi, stop=pi, length=xsize))
    longitude = collect(range(-180.0, stop=180.0, length=xsize))/180.0*pi
    latitude = collect(range(90.0, stop=-90.0, length=ysize))/180.0*pi
    # project the map to a rectangular matrix xsize x ysize
    PHI = [i for j in theta, i in phi]
    THETA = [j for j in theta, i in phi]
    grid_pix = longlat_ang2pix(ntheta, nphi, THETA, PHI); # for long lat scheme
    #grid_map = image[grid_pix]
    grid_map = image[circshift(grid_pix,(0,Int(xsize/2)))]
    fig = figure(figtitle, figsize=(10, 7))
    clf();
    ax = subplot(111,projection="mollweide",title=title)
    # rasterized makes the map bitmap while the labels remain vectorial
    # flip longitude to the astro convention
    # if visible_pixels == []
    #     vmin = minimum(image);
    #     vmax = maximum(image)*1.01;
    # else
    #    vmin = minimum(image[visible_pixels]);
    #    vmax = maximum(image[visible_pixels])*1.01;
    # end
    if visible_pixels == []
        if (vmin == -Inf)
        vmin = minimum(image);
        end
        if (vmax == Inf)
        vmax = maximum(image);
        end
    else
        if (vmin == -Inf)
        vmin = minimum(image[visible_pixels]);
        end
        if (vmax == Inf)
        vmax = maximum(image[visible_pixels]);
        end
    end
    moll = pcolormesh(longitude, latitude, grid_map, vmin=vmin, vmax=vmax, rasterized=true, cmap=colormap)
    # graticule
    ax.set_longitude_grid(30)
    ax.set_latitude_grid(30)
    ax.set_longitude_grid_ends(90)
    spacing = 0.04
    subplots_adjust(bottom=spacing, top=1-spacing, left=spacing, right=1-spacing)
    grid(true)

    ticks = collect(range(vmin, stop=vmax, length=7));
    cb = colorbar(moll, orientation="horizontal", shrink=.6, pad=0.05, ticks=ticks)
    cb.ax.xaxis.labelpad=5
    cb.ax.xaxis.set_label_text(clabel)
    # workaround for issue with viewers, see colorbar docstring
    cb.solids.set_edgecolor("face")
    ax.tick_params(axis="x", labelsize=15)
    ax.tick_params(axis="y", labelsize=15)
end
