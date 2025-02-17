using OptimPackNextGen

function split_didata_by_period(didata, P)
    mjd_intervals = [range(minimum(didata.mjd), maximum(didata.mjd), step=P);Inf];
    nframes = length(mjd_intervals)-1;
    frames_mjd = Vector{Float64}(undef, nframes)
    didata_split = Array{observedLineProfiles}(undef, nframes)
    epoch_ixs = Array{Vector{Int}}(undef, nframes)
    for i=1:nframes
        epoch_filter = findall((didata.mjd .>= mjd_intervals[i]) .& (didata.mjd .< mjd_intervals[i+1]))
        didata_split[i] = observedLineProfiles(didata.fname[epoch_filter], length(epoch_filter), didata.nwvls, didata.snr[epoch_filter], didata.mjd[epoch_filter], didata.phase[epoch_filter], didata.wavelength, didata.intensities[epoch_filter])
        frames_mjd[i] = mean(didata.mjd[epoch_filter])
        epoch_ixs[i] = vec(epoch_filter)
    end
    return nframes, didata_split, epoch_ixs, frames_mjd
end

function di_reconstruct_multitemporal(x_start, didata_split, epoch_ixs, nframes, polyflux, star_epoch_geom, star_epoch_μ, star_epoch_v, stellar_parameters, model_grid, vsini, visible_pixels, lowtemp, hightemp; verb = true, maxiter = 200, regularizers = [], wavelength_regions = [])
    crit_imaging = (x,g) -> di_multitemporal_fg(x, g, didata_split, epoch_ixs, nframes, polyflux, star_epoch_geom, star_epoch_μ, star_epoch_v, stellar_parameters, model_grid, vsini, visible_pixels, wavelength_regions=wavelength_regions, verb = verb, regularizers = regularizers, lowtemp=lowtemp, hightemp=hightemp);
    x_sol = OptimPackNextGen.vmlmb(crit_imaging, x_start, verb=false, maxiter=maxiter, blmvm=false, gtol=(0,1e-8), lower=lowtemp, upper=hightemp);
    return x_sol
end

function di_multitemporal_fg(x, g, didata_split, epoch_ixs, nframes, polyflux, star_epoch_geom, star_epoch_μ, star_epoch_v, stellar_parameters, model_grid, vsini, visible_pixels; regularizers=[], verb = true, wavelength_regions=[], lowtemp=0, hightemp=Inf)
    g .= zeros(length(g))
    npix = star_epoch_geom[1].npix

    f = 0.0;
    printstrings = Array{String}(undef, nframes)
    @views @inbounds for i=1:nframes
        tslice = 1+(i-1)*npix:i*npix; # temporal slice
        subg = zeros(npix)
        (critf, printstring) = di_crit_multitemporal_fg(x[tslice], subg, didata_split[i], polyflux[epoch_ixs[i], :], star_epoch_geom[epoch_ixs[i]], star_epoch_μ[:, epoch_ixs[i]], star_epoch_v[:, epoch_ixs[i]], stellar_parameters, model_grid, vsini, visible_pixels, regularizers=regularizers[i], verb = verb)
        f += critf
        g[tslice] .= subg
        printstrings[i] = printstring
    end
    if verb==true
        for j=1:nframes
            printstyled("Frame $j ",color=128);  # purple
            print(printstrings[j])
        end
    end

    if length(regularizers) > nframes
        if (regularizers[nframes+1][1][1] == "Δt") & (regularizers[nframes+1][1][2] > 0.0) & (nframes>1)
                tvinfo = regularizers[nframes+1][1][3]
                y = reshape(x, (npix, nframes))
    
                Δtg = zeros(size(y))
                Δtf = calculate_Δt_fg(y, Δtg, nframes)
                f += regularizers[nframes+1][1][2] * Δtf
                g .+= regularizers[nframes+1][1][2] .* vec(Δtg)
    
                if verb == true
                    tempstring = @sprintf("Temporal regularization: %1.4e\n", Δtf)
                    printstyled(tempstring, color=214)  # orange
                    fstring = @sprintf("Total f: %8.2f\n", f)
                    printstyled(fstring, color=214)  # orange
                end
        else
            printstyled("-------------------\n", color=214)
        end
    end

    return f
end

function di_crit_multitemporal_fg(x, g, didata, polyflux, star_epoch_geom, star_epoch_μ, star_epoch_v, stellar_parameters, model_grid, vsini, visible_pixels; regularizers=[], verb = true, wavelength_regions=[])
    ## Compute χ^2 and ∇χ^2
    if wavelength_regions == []  # Specific spectral lines are selected here
        ixs = collect(1:length(model_grid.wavelength))
    else
        ixs = vcat([findall(bounds[1] .< model_grid.wavelength .< bounds[2]) for bounds in wavelength_regions]...)
        sort!(ixs)
    end

    global_line_profiles = calculate_global_profiles([x], star_epoch_geom, star_epoch_μ, star_epoch_v, stellar_parameters, polyflux, model_grid, vsini)
    res = Array{Vector{Float64}}(undef, length(star_epoch_geom))
    @inbounds for i=eachindex(star_epoch_geom)
        global_line_profiles.intensities[i] ./= (mean(global_line_profiles.intensities[i][ixs]) / mean(didata.intensities[i][ixs]))  # renormalize
        res[i] = didata.intensities[i][ixs] .- global_line_profiles.intensities[i][ixs]
    end
    
    chi2_f = calculate_chi2_f(res, didata.snr)  # Calculate χ^2, weighted by SNR
    g .= calculate_chi2_g(x, g, star_epoch_geom, stellar_parameters, model_grid, polyflux, vsini, res, didata.snr, wavelength_regions=wavelength_regions)  # calculate ∇χ^2
    reg_f = spheroid_regularization(x, g, regularizers=regularizers, verb = false);  #Apply regularizers; note: this adds to g
    if verb == true
        y = x[visible_pixels]
        low_x = minimum(y);
        high_x = maximum(y);
        mean_x = mean(y);
        printstring = @sprintf("ftot: %8.2f\tχ^2: %8.2f\tχ^2r: %4.2f\tμ*freg:%8.2f\tminT: %5.1f\tmaxT: %5.1f\tmeanT: %5.1f\n", chi2_f+reg_f, chi2_f, chi2_f/length(ixs)/length(res), reg_f, low_x, high_x, mean_x);
    end

    return chi2_f + reg_f, printstring
end

#= 
function di_reconstruct_multitemporal(x_start, Y, E, H, star_epoch_geom, stellar_parameters, model_grid, vsini, visible_pixels, lowtemp, hightemp, mjd, ntheta, nphi; relative = false, printcolor= [], verb = true, maxiter = 200, regularizers =[], wavelength_regions=[])
    crit_imaging = (x,g) -> di_multitemporal_fg(x, g, Y, E, H, star_epoch_geom, stellar_parameters, model_grid, vsini, visible_pixels, mjd, ntheta, nphi, wavelength_regions=wavelength_regions, verb = verb, regularizers = regularizers, lowtemp=lowtemp, hightemp=hightemp);
    x_sol = OptimPackNextGen.vmlmb(crit_imaging, x_start, verb=false, maxiter=maxiter, blmvm=false, gtol=(0,1e-8), lower=0.0, mem=20);
    g_dummy = zeros(size(x_start));
    return x_sol
end

function di_multitemporal_fg(x, g, I_j, SNR_j, W_j, star_epoch_geom, stellar_parameters, model_grid, vsini, visible_pixels, mjd, ntheta, nphi; regularizers=[], verb = true, wavelength_regions=[], lowtemp=0, hightemp=Inf)
    g .= zeros(length(g))
    npix = size(W_j[1], 2)
    δθ = 180 / ntheta
    δϕ = 180 / nphi
    nframes = length(I_j)
    # δθ = δϕ = 360 / sqrt(npix)

    x[1:npix*nframes] .= max.(x[1:npix*nframes], lowtemp)
    x[1:npix*nframes] .= min.(x[1:npix*nframes], hightemp)
    # x[end-npix+1:end] .= min.(x[end-npix+1:end], 0.0)

    printstrings = Array{String}(undef, nframes)

    f = 0.0;
    for i=1:nframes
        nphase = length(I_j[i])
        tslice = 1+(i-1)*npix:i*npix; # temporal slice
        tslice_geom = 1+(i-1)*nphase:i*nphase
        subg = zeros(npix)
        (critf, printstring) = di_crit_fg(x[tslice], subg, I_j[i], SNR_j[i], W_j[i], star_epoch_geom[tslice_geom], stellar_parameters, model_grid, vsini, visible_pixels, regularizers=regularizers[i], verb = verb)
        f += critf
        g[tslice] .= subg
        printstrings[i] = printstring
    end
    if verb==true
        for j=1:nframes
            printstyled("Frame $j ",color=128);  # purple
            print(printstrings[j])
        end
    end
    # f = sum(f_threads)

    if length(regularizers) > nframes
        if (regularizers[nframes+1][1][1] == "flow")  & (regularizers[nframes+1][1][2] > 0.0) & (nframes>1)
            tvinfo = regularizers[nframes+1][1][3]
            y = reshape(x, (npix, nframes+1))
            # y = reshape(x[1:end-nphi], (npix, nframes))
            # vθ = x[end-nphi+1:end]
            # vθ = vec(ones((nphi, ntheta)) .* vθ)

            # dRdI = similar(y)
            # fill!(dRdI, 0.0)
            # dRdv = similar(vθ)
            # fill!(dRdv, 0.0)

            flowg = zeros(size(y))
            flowf = calculate_flow_fg(y, flowg, mjd, nframes, npix, ntheta, nphi, tvinfo, δθ, δϕ, star_epoch_geom)
            f += regularizers[nframes+1][1][2] * flowf
            g .+= regularizers[nframes+1][1][2] .* vec(flowg)

            flowtvgθ = zeros(npix)
            fvθtv = 0.0
            fvθtv = spheroid_regularization(y[:, end], flowtvgθ, regularizers=regularizers[nframes+2], verb=false)
            f += regularizers[nframes+2][1][2] * fvθtv# + fvϕtv)
            g[end-npix+1:end] .+= regularizers[nframes+2][1][2] .* flowtvgθ
            
            if verb == true
                tempstring = @sprintf("Temporal regularization: %1.4e\tFlow TV: %6.2f\n", flowf, fvθtv)#+fvϕtv)
                printstyled(tempstring, color=214)  # orange
                fstring = @sprintf("Total f: %8.2f\n", f)
                printstyled(fstring, color=214)  # orange
            end

        elseif (regularizers[nframes+1][1][1] == "Δt") & (regularizers[nframes+1][1][2] > 0.0) & (nframes>1)
                tvinfo = regularizers[nframes+1][1][3]
                y = reshape(x, (npix, nframes))
    
                Δtg = zeros(size(y))
                Δtf = calculate_Δt_fg(y, Δtg, nframes)
                f += regularizers[nframes+1][1][2] * Δtf
                g .+= regularizers[nframes+1][1][2] .* vec(Δtg)
    
                if verb == true
                    tempstring = @sprintf("Temporal regularization: %1.4e\n", Δtf)
                    printstyled(tempstring, color=214)  # orange
                    fstring = @sprintf("Total f: %8.2f\n", f)
                    printstyled(fstring, color=214)  # orange
                end
        else
            printstyled("-------------------\n", color=214)
        end
    end

    return f
end

function calculate_flow_fg(y, g, mjd, nframes, npix, ntheta, nphi, tvinfo, δθ, δϕ, star_epoch_geom)
    fflow = 0.0
    mθ = zeros(npix);
    mϕ = zeros(npix)
    vθ = y[:, nframes+1]
    # mollplot_temperature_longlat(tvinfo[3], 40, 80)
    # show()
    
    # vθ = vec((x[end-ntheta+1:end] * ones(nphi)')')
    # y = reshape(x[1:npix*nframes], npix, nframes)
    # g_reshape = reshape(g[1:npix*nframes], npix, nframes)
    # dRdv = g[end-ntheta+1:end]

    for j=1:nframes-1
        δt = mjd[j+1] - mjd[j]
        mθ .= (δt / δθ) .* vθ

        F_Ij = evolve_frame(y[:, j], mθ, mϕ, δt, δθ, δϕ, tvinfo, star_epoch_geom[j])
        # mollplot_temperature_longlat(F_Ij, 40, 80)
        # show()F
        diff = y[:, j+1] .- F_Ij
        # mollplot_temperature_longlat(diff, 40, 80)
        # show()

        fflow += norm(diff, 2)^2

        ## Calculate dRflow / dI
        if j > 1
            F_Ijm1 = evolve_frame(y[:, j-1], mθ, mϕ, mjd[j] - mjd[j-1], δθ, δϕ, tvinfo, star_epoch_geom[j-1])
            g[:, j] .+= 2.0 .* (y[:, j] - F_Ijm1)  # 2(I_j - F_flow ⋅ I_j-1)
        end
        if j < nframes
            g[:, j] .+= -2.0 .* (diff .+ mθ .* ∇θ(diff, tvinfo))
        end

        ## Calculate dRflow / dm
        # g[:, nframes+1] .+= (δt / δθ) .* (-2.0 .* y[:, j] .* ∇θ(diff, tvinfo))
        g[:, nframes+1] .+= (δt / δθ) .* ((2.0 .* diff .* ∇θ(y[:, j], tvinfo)) .- ∇θ(2.0 .* diff .* y[:, j], tvinfo))
    end

    δt = mjd[nframes] - mjd[nframes-1]
    mθ .= (δt / δθ) .* vθ
    F_Ijm1 = evolve_frame(y[:, nframes-1], mθ, mϕ, δt, δθ, δϕ, tvinfo, star_epoch_geom[nframes-1])
    g[:, nframes] .+= 2.0 .* (y[:, nframes] - F_Ijm1)

    return fflow
end

function calculate_Δt_fg(y, g, nframes)
    Δtf = 0.0
    for i=1:nframes-1
        Δtf += norm(y[:, i] .- y[:, i+1], 2)^2
        if i > 1
           g[:, i] .+= 2.0 .* (y[:, i] .- y[:, i-1])
        end
        if i < nframes
            g[:, i] .+= 2.0 .* (y[:, i] .- y[:, i+1])
        end
    end

    g[:, nframes] .+= 2.0 .* (y[:, nframes] .- y[:, nframes-1])
    return Δtf
end

function ∇θ(I, tvinfo; type="backwards")
    diff = zeros(length(I))
    if type=="backward"
        diff .= I .- I[tvinfo[3]]  # x_i - x_i-1
    elseif type=="forward"
        for i=1:length(I)
            eix = tvinfo[1][i][4]
            diff[i] = I[eix] - I[i]  # x_i+1 - x_i
        end
    elseif type=="centered"
        for i=1:length(I)
            wix = tvinfo[1][i][3]
            eix = tvinfo[1][i][4]
            diff[i] = (I[wix] - I[eix]) / 2 
        end
    end
    return diff
end

function F_flow(I, mθ, mϕ, tvinfo)
    ## Computes F_flow ⋅ I
    Im = I .* mθ
    ΔI = -∇θ(Im, tvinfo)
    return I .+ ΔI
end

function evolve_frame(y, mθ, mϕ, δt, δθ, δϕ, tvinfo, star_epoch_geom)
    evolved_frame = evaluate_courant(y, mθ, 0.0.*mϕ, δt, δθ, δϕ, tvinfo, star_epoch_geom)  # x-velocity first
    # evolved_frame .= evaluate_courant(evolved_frame, 0.0.*mθ, mϕ, δt, δθ, δϕ, tvinfo, star_epoch_geom)  # then y-velocity
    return evolved_frame
end

function evaluate_courant(y, mθ, mϕ, δt, δθ, δϕ, tvinfo, star_epoch_geom)
    ## check if Courant condition is satisfied (v ⋅ δt/δx ≤ 1)
    evolved_frame = deepcopy(y)
    if (any(abs.(mθ) .> 1.0) || any(abs.(mϕ) .> 1.0))
        n_artificial_frames = ceil(max(maximum(abs.(mθ)), maximum(abs.(mϕ))))
        # println("nframes(CFL): $(n_artificial_frames)")
        δt_artificial = δt / n_artificial_frames
        mθ_artificial = mθ .* δt_artificial/δt
        mϕ_artificial = mϕ .* δt_artificial/δt
        for k=1:n_artificial_frames
            evolved_frame .= F_flow(evolved_frame, mθ_artificial, mϕ_artificial, tvinfo)
            # mollplot_temperature_longlat(evolved_frame, 40, 80)
            # show()
        end
    else
        evolved_frame .= F_flow(y, mθ, mϕ, tvinfo)  # I_j+1 - F_flow ⋅ I_j
    end

    return evolved_frame
end

function di_crit_fg(x, g, I_j, SNR_j, W_j, star_epoch_geom, stellar_parameters, model_grid, vsini, visible_pixels; regularizers=[], verb = true, wavelength_regions=[])
    ## Compute χ^2 and ∇χ^2

    if wavelength_regions == []  # Specific spectral lines are selected here
        ixs = collect(1:length(model_grid.wavelength))
    else
        ixs = vcat([findall(bounds[1] .< model_grid.wavelength .< bounds[2]) for bounds in wavelength_regions]...)
        sort!(ixs)
    end

    global_line_profiles = calculate_global_profiles(x, star_epoch_geom, stellar_parameters, W_j, model_grid, vsini)
    res = Array{Vector{Float64}}(undef, length(star_epoch_geom))
    for i=1:length(star_epoch_geom)
        # global_line_profiles.intensities[i] .= renorm(global_line_profiles.wavelength, global_line_profiles.intensities[i], F[i], ixs)
        global_line_profiles.intensities[i] ./= (mean(global_line_profiles.intensities[i][ixs]) / mean(I_j[i][ixs]))  # renormalize
        res[i] = I_j[i][ixs] .- global_line_profiles.intensities[i][ixs]
    end
    # plot_profiles(global_line_profiles, didata=didata, line_centers=[6411.6476, 6419.9496, 6421.3495, 6430.8446, 6439.0753], line_labels=["Fe I λ6411", "Fe I λ6419", "Fe I λ6421", "Fe I λ6430", "Ca I λ6439"], vsini=vsini)
    # show()

    chi2_f = calculate_chi2_f(res, SNR_j)  # Calculate χ^2, weighted by SNR
    g .= calculate_chi2_g(x, g, star_epoch_geom, stellar_parameters, model_grid, W_j, vsini, res, SNR_j, wavelength_regions=wavelength_regions)  # calculate ∇χ^2
    reg_f = spheroid_regularization(x, g, regularizers=regularizers, verb = false);  #Apply regularizers; note: this adds to g
    if verb == true
        y = x[visible_pixels]
        low_x = minimum(y);
        high_x = maximum(y);
        mean_x = mean(y);
        printstring = @sprintf("ftot: %8.2f\tχ^2: %8.2f\tχ^2r: %4.2f\tμ*freg:%8.2f\tminT: %5.1f\tmaxT: %5.1f\tmeanT: %5.1f\n", chi2_f+reg_f, chi2_f, chi2_f/length(ixs)/length(res), reg_f, low_x, high_x, mean_x);
    end

    return chi2_f + reg_f, printstring
end
=#