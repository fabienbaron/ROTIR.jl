# =======================================================================================
# Binary movie rendering (matplotlib)
# =======================================================================================
# Split out of src/animation.jl because it is the ONLY part of the animation path that draws:
# `binary_frame_maps`, `_binary_value_range`, `_binary_axis_max` and `frames_to_movie` are all
# pure Julia (or a shell-out to ffmpeg) and stay in the core module. This file is included by
# ROTIRPythonPlotExt, so `using ROTIR` alone no longer pulls matplotlib in on its account.

"""
    binary_movie(bparams, tessels1, params1, tessels2, params2;
                 tstart=bparams.T0, tstop=bparams.T0 + bparams.P, nframes=200,
                 panels=:single, outdir, prefix="frame", fps=20, encode=true, ...)
        -> (outdir, moviepath_or_nothing)

Render a binary through its orbit and (optionally) encode a movie.

`panels`:
* `:single`  — one panel, the irradiated system.
* `:compare` — three panels: irradiation on, irradiation off, and the ΔT difference. This
  is the one that actually shows the effect; on its own the heated star just looks like a
  slightly different star.

Physics keywords (`reflection`, `albedo1/2`, `ldbol1/2`, `method`, `volume_conserving`)
are passed through to [`binary_frame_maps`](@ref); rendering keywords (`intensity`,
`intensity_model`, `band`, `colormap`, `graticules`, …) to `plot2d_binary`.

The colour scale and axis limits are fixed for the whole sequence: a coarse pre-scan over
`nscan` frames finds the global ranges, which are then held constant. Without that a
per-frame autoscale would rescale away exactly the variation being looked for.

`volume_conserving=true` is expensive per frame; the orbit-wide Ω(D) tables from
[`roche_omega_table`](@ref) are built once here and interpolated, so the cost is paid
`2 × npoints` times rather than `2 × nframes`.

Returns the frame directory and the movie path (`nothing` if `encode=false` or `ffmpeg`
is unavailable — the PNGs are kept either way and the command is printed).
"""
function binary_movie(bparams, tessels1::tessellation, params1,
                      tessels2::tessellation, params2;
                      tstart = bparams.T0, tstop = bparams.T0 + bparams.P,
                      nframes::Int = 200, nscan::Int = 24,
                      panels::Symbol = :single,
                      outdir::AbstractString = "movie_frames", prefix = "frame",
                      fps::Int = 20, encode::Bool = true, dpi::Int = 110,
                      reflection::Bool = true, albedo1 = 0.6, albedo2 = 0.6,
                      ldbol1 = (ldtype = 0,), ldbol2 = (ldtype = 0,),
                      method::Symbol = :horvat, volume_conserving::Bool = false,
                      T = promote_type(eltype(tessels1), eltype(tessels2)),
                      kernel_eltype::Type = Float64,
                      intensity::Bool = true, intensity_model::Symbol = :planck,
                      band = 1.65e-6, colormap = "gist_heat",
                      diff_colormap = "inferno", figsize = nothing, pad = 0.25,
                      title_fn = nothing, verbose::Bool = true,
                      plot_kwargs = (;))
    panels in (:single, :compare) ||
        error("binary_movie: panels must be :single or :compare (got $(panels))")
    mkpath(outdir)
    times = collect(range(Float64(tstart), Float64(tstop), length=nframes))

    # Volume-conserving Ω is a root solve over a Romberg quadrature; tabulate it once.
    Ω1tab = Ω2tab = nothing
    if volume_conserving
        verbose && @info "binary_movie: tabulating volume-conserving Ω(D)…"
        Ω1tab = roche_omega_table(tessels1, params1, bparams; secondary=false)
        Ω2tab = roche_omega_table(tessels2, params2, bparams; secondary=true)
    end
    _omegas(t) = begin
        Ω1tab === nothing && return (nothing, nothing)
        D = compute_separation(bparams, t)
        (Ω1tab(D), Ω2tab(D))
    end

    frame(t) = begin
        o1, o2 = _omegas(t)
        binary_frame_maps(tessels1, params1, tessels2, params2, bparams, t;
                          reflection=reflection, albedo1=albedo1, albedo2=albedo2,
                          ldbol1=ldbol1, ldbol2=ldbol2, method=method,
                          volume_conserving=false, omega1=o1, omega2=o2,
                          T=T, kernel_eltype=kernel_eltype)
    end

    # ---- pre-scan: global colour range and axis limits -------------------------------
    verbose && @info "binary_movie: pre-scan over $(min(nscan, nframes)) frames for fixed scales…"
    scan_idx = unique(round.(Int, range(1, nframes, length=min(nscan, nframes))))
    vlo = Inf; vhi = -Inf; dlo = Inf; dhi = -Inf; amax = 0.0
    for k in scan_idx
        s1, s2, m1, m2, i1, i2 = frame(times[k])
        lo, hi = _binary_value_range(s1, m1, s2, m2; intensity=intensity,
                                     intensity_model=intensity_model, band=band)
        vlo = min(vlo, lo); vhi = max(vhi, hi)
        if panels === :compare
            d1 = m1 .- i1; d2 = m2 .- i2
            dlo = min(dlo, minimum(d1[s1.index_quads_visible]),
                           minimum(d2[s2.index_quads_visible]))
            dhi = max(dhi, maximum(d1[s1.index_quads_visible]),
                           maximum(d2[s2.index_quads_visible]))
        end
        amax = max(amax, _binary_axis_max(s1, s2, bparams, times[k]; pad=pad))
    end
    if panels === :compare && !(dhi > dlo)
        dlo, dhi = 0.0, max(dhi, 1.0)     # reflection off, or a flat difference
    end
    verbose && @info "binary_movie: colour range [$(round(vlo, sigdigits=5)), $(round(vhi, sigdigits=5))], axis ±$(round(amax, digits=3)) mas"

    # ---- render ----------------------------------------------------------------------
    P = bparams.P
    fsize = figsize === nothing ? (panels === :single ? (8, 8) : (19, 6.5)) : figsize
    fig = figure("ROTIR binary movie", figsize=fsize, facecolor="White")
    axs = panels === :single ? [fig.add_subplot(1, 1, 1)] :
                               [fig.add_subplot(1, 3, k) for k in 1:3]
    files = String[]
    for (k, t) in enumerate(times)
        s1, s2, m1, m2, i1, i2 = frame(t)
        phase = mod((t - bparams.T0) / P, 1.0)
        sep = projected_separation(bparams, t)
        head = title_fn === nothing ?
               @sprintf("phase %.3f   ρ = %.3f mas", phase, sep) :
               title_fn(t, phase, sep)
        for a in axs; a.clear(); end
        if panels === :single
            plot2d_binary(m1, m2, s1, s2, bparams, t; ax=axs[1], colorbar_on=false,
                          intensity=intensity, intensity_model=intensity_model, band=band,
                          vmin=vlo, vmax=vhi, axis_max=amax, colormap=colormap,
                          figtitle=head, plot_kwargs...)
        else
            plot2d_binary(m1, m2, s1, s2, bparams, t; ax=axs[1], colorbar_on=false,
                          intensity=intensity, intensity_model=intensity_model, band=band,
                          vmin=vlo, vmax=vhi, axis_max=amax, colormap=colormap,
                          figtitle="irradiated — $head", plot_kwargs...)
            plot2d_binary(i1, i2, s1, s2, bparams, t; ax=axs[2], colorbar_on=false,
                          intensity=intensity, intensity_model=intensity_model, band=band,
                          vmin=vlo, vmax=vhi, axis_max=amax, colormap=colormap,
                          figtitle="intrinsic only", plot_kwargs...)
            plot2d_binary(m1 .- i1, m2 .- i2, s1, s2, bparams, t; ax=axs[3],
                          colorbar_on=false, intensity=false,
                          vmin=dlo, vmax=dhi, axis_max=amax, colormap=diff_colormap,
                          figtitle=@sprintf("ΔT (K), %.0f – %.0f", dlo, dhi),
                          plot_kwargs...)
        end
        f = joinpath(outdir, @sprintf("%s_%04d.png", prefix, k))
        fig.savefig(f, dpi=dpi, bbox_inches="tight", facecolor="white")
        push!(files, f)
        verbose && (k % 25 == 0 || k == nframes) && @info "binary_movie: frame $k/$nframes"
    end
    pyplot.close(fig)

    movie = encode ? frames_to_movie(outdir, prefix; fps=fps, verbose=verbose) : nothing
    return outdir, movie
end
