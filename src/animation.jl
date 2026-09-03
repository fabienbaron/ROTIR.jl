# animation.jl
# ---------------------------------------------------------------------------
# Frame sequences and movies of a binary over its orbit.
#
# Everything a frame needs already exists: `create_binary_geometry` gives both meshes in
# one frame at a real epoch, `temperature_map_vonZeipel_roche_single` (via
# `parametric_temperature_map`) gives the gravity-darkened maps, `handle_reflection` heats
# them, and `plot2d_binary` renders them. This file only supplies the loop, plus the two
# things a *sequence* needs that a single plot does not:
#
#   * a colour scale and axis limits pinned across all frames — otherwise the per-frame
#     autoscale makes the star appear to change when only the scale did, which would
#     completely mask a few-percent irradiation signature;
#   * an encode step.
# ---------------------------------------------------------------------------

"""
    binary_frame_maps(tessels1, params1, tessels2, params2, bparams, tepoch_jd;
                      reflection=true, albedo1=0.6, albedo2=0.6,
                      ldbol1=(ldtype=0,), ldbol2=(ldtype=0,), method=:horvat,
                      volume_conserving=false, omega1=nothing, omega2=nothing,
                      T=promote_type(eltype(tessels1), eltype(tessels2)),
                      kernel_eltype=Float64)
        -> (star1, star2, tmap1, tmap2, tmap1_intrinsic, tmap2_intrinsic)

One frame's worth of physics: both meshes in the common orbital frame, their intrinsic
(gravity-darkened) temperature maps, and — if `reflection` — the irradiation-heated maps.

Both the heated and the intrinsic maps are returned so a caller can render the difference
without recomputing anything. With `reflection = false` the two pairs are identical.
"""
function binary_frame_maps(tessels1::tessellation, params1,
                           tessels2::tessellation, params2,
                           bparams, tepoch_jd;
                           reflection::Bool = true, albedo1 = 0.6, albedo2 = 0.6,
                           ldbol1 = (ldtype = 0,), ldbol2 = (ldtype = 0,),
                           method::Symbol = :horvat,
                           volume_conserving::Bool = false,
                           omega1 = nothing, omega2 = nothing,
                           T = promote_type(eltype(tessels1), eltype(tessels2)),
                           kernel_eltype::Type = Float64,
                           kernels = nothing)
    star1, star2 = create_binary_geometry(tessels1, params1, tessels2, params2,
                                          bparams, tepoch_jd;
                                          T=T, volume_conserving=volume_conserving,
                                          omega1=omega1, omega2=omega2)
    t0 = parametric_temperature_map(params1, star1)
    t2 = parametric_temperature_map(params2, star2; secondary=true)
    if reflection
        h1, h2 = handle_reflection(star1, t0, star2, t2;
                                   albedo1=albedo1, albedo2=albedo2,
                                   ldbol1=ldbol1, ldbol2=ldbol2, method=method,
                                   kernel_eltype=kernel_eltype)
        return star1, star2, h1, h2, t0, t2
    else
        return star1, star2, t0, t2, t0, t2
    end
end

# Shared value range across both components, matching what plot2d_binary will colour by.
function _binary_value_range(star1, tmap1, star2, tmap2;
                             intensity::Bool, intensity_model::Symbol, band)
    if intensity
        I1 = intensity_model === :linear ? tmap1 : ROTIR.intensity(tmap1, intensity_model, band)
        I2 = intensity_model === :linear ? tmap2 : ROTIR.intensity(tmap2, intensity_model, band)
        p1 = I1[star1.index_quads_visible] .* star1.ldmap[star1.index_quads_visible]
        p2 = I2[star2.index_quads_visible] .* star2.ldmap[star2.index_quads_visible]
    else
        p1 = tmap1[star1.index_quads_visible]
        p2 = tmap2[star2.index_quads_visible]
    end
    return min(minimum(p1), minimum(p2)), max(maximum(p1), maximum(p2))
end

# Same formula plot2d_binary uses, so a fixed axis_max cannot clip a frame.
function _binary_axis_max(star1, star2, bparams, tepoch; pad = 1.0)
    ow, on = orbit_to_rotir_offset(bparams, Float64(tepoch))
    r1 = maximum(sqrt.(star1.vertices_xyz[:,:,1].^2 .+ star1.vertices_xyz[:,:,2].^2))
    r2 = maximum(sqrt.(star2.vertices_xyz[:,:,1].^2 .+ star2.vertices_xyz[:,:,2].^2))
    return max(r1, r2 + abs(ow), r2 + abs(on)) + pad
end

"""
    frames_to_movie(dir, prefix; fps=20, out=nothing, verbose=true) -> path or nothing

Encode `dir/prefix_%04d.png` into an mp4 with `ffmpeg`. Returns the movie path, or
`nothing` (after printing the command) if `ffmpeg` is not on `PATH` or the encode fails.
The frames are never deleted.
"""
function frames_to_movie(dir::AbstractString, prefix::AbstractString = "frame";
                         fps::Int = 20, out = nothing, verbose::Bool = true)
    outfile = out === nothing ? joinpath(dir, "$(prefix).mp4") : out
    pattern = joinpath(dir, "$(prefix)_%04d.png")
    cmd = `ffmpeg -y -loglevel error -framerate $fps -i $pattern
           -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -c:v libx264 -pix_fmt yuv420p $outfile`
    if Sys.which("ffmpeg") === nothing
        @warn "frames_to_movie: ffmpeg not found on PATH; frames kept in $dir. Encode with:\n  $cmd"
        return nothing
    end
    try
        run(cmd)
        verbose && @info "frames_to_movie: wrote $outfile"
        return outfile
    catch err
        @warn "frames_to_movie: ffmpeg failed ($err); frames kept in $dir. Try:\n  $cmd"
        return nothing
    end
end
