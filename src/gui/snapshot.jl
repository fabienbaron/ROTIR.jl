# ── "Save view": what a plot panel writes to a file ──────────────────────────────────────
#
# NOTHING HERE READS THE FRAMEBUFFER ON SCREEN, and it cannot. QMLMakie swaps the screen's last
# postprocessor for one that renders into Qt's FBO, so `Makie.colorbuffer` on the live screen
# hands back the raw buffer instead of the composited frame, and `Makie.save` on a live figure
# throws before it gets that far — OITOOLS measured both and its src/gui/snapshot.jl carries
# the finding.
#
# A figure that was never DISPLAYED has neither problem: `Makie.save` opens its own hidden
# GLMakie screen for it, which works while Qt holds the context. So a snapshot REBUILDS the
# panel offscreen through the very builders the window uses, and saves that. What lands in the
# file is then what the equivalent script would draw — the same property the command log has,
# and the reason this is not a screen grab.

"""
    _snapshot_figure(sh, which, size) -> Figure or String

Rebuild one panel offscreen. Returns the Figure, or a message beginning with `!`.

`which` names a PLOT AREA, not a tab: a tab showing three views would otherwise have to guess
which was meant, and the answer differs from what is on screen the moment the user switches.
"""
function _snapshot_figure(sh::ShellState, which::String, sz::Tuple{Int,Int})
    fig = Makie.Figure(size = sz)
    got = build_epoch_star(sh)
    d = current_dataset(sh.session)
    m = current_model(sh.session)

    if which in ("sky", "mollweide", "star3d")
        got === nothing && return "! no model to draw — add one on the Model tab"
        star, tmap = got[1], got[2]
        # BOTH COMPONENTS OF A BINARY, exactly as the live canvases take it. This destructured
        # `star, tmap = got` and always called the single-star path, so a binary saved as one
        # star — framed on the primary, and on the primary's colour range. `build_epoch_star`
        # returns five values for a pair, and the extra three are the whole companion.
        bin2 = length(got) >= 5 ? (got[3], got[4], got[5]) : nothing
        allv = surface_values(sh, tmap, star; visible_only = false)
        if which == "sky"
            c = build_sky_canvas(fig)
            c.cbarlabel[] = sh.intensity[] ? "I (arb.)" : "T (K)"
            ttl = d === nothing ? "" :
                  Printf.@sprintf("epoch %d — MJD %.4f", sh.session.current_epoch,
                                  d.mjd[clamp(sh.session.current_epoch, 1, length(d.mjd))])
            if bin2 === nothing
                show_map!(c, star, visible_values(sh, allv, star);
                          title = ttl, star_params = star_params(m), _decor(sh)...)
            else
                star2, tmap2, off = bin2
                show_binary_map!(c, star, visible_values(sh, allv, star), star2,
                                 surface_values(sh, tmap2, star2; visible_only = true), off;
                                 title = ttl, star_params = star_params(m),
                                 star_params2 = star_params(m.companion), _decor(sh)...)
            end
        elseif which == "mollweide"
            # The PRIMARY, deliberately: a Mollweide is a map of one surface, and two of them
            # in one frame is two pictures. The live view makes the same choice.
            c = build_moll_canvas(fig)
            c.cbarlabel[] = whole_surface_label(sh)
            show_mollweide!(c, allv, star)
        else
            c = build_star_canvas(fig)
            c.cbarlabel[] = whole_surface_label(sh)
            if bin2 === nothing
                show_star3d!(c, star, allv)
            else
                star2, tmap2, off = bin2
                cc = m.companion
                # The line-of-sight term the sky view drops, which is the reason to look at a
                # pair in three dimensions at all.
                zoff = cc === nothing ? 0.0 :
                       (cc.place === :offset ? Float64(cc.offset[3]) : 0.0)
                track = if cc !== nothing && cc.place === :orbit
                    try
                        relative_orbit_track(orbit_bparams(sh.orbit; binary = m))
                    catch
                        Makie.Point3f[]
                    end
                else
                    Makie.Point3f[]
                end
                show_binary3d!(c, star, allv, star2,
                               surface_values(sh, tmap2, star2; visible_only = false),
                               (Float64(off[1]), Float64(off[2]), zoff); track = track)
            end
        end
        return fig
    elseif which == "chi2"
        b = epoch_chi2(sh)
        b === nothing && return "! no χ² yet — a dataset and a model are both needed"
        show_chi2!(build_chi2_canvas(fig), b)
        return fig
    elseif which == "obs"
        d === nothing && return "! no dataset loaded"
        O = _oitools_gui()
        O === nothing && return "! the observable plots need OITOOLSGUIExt"
        c, pts = build_obs_canvas(fig)
        c === nothing && return "! the observable canvas could not be built"
        i = clamp(sh.session.current_epoch, 1, length(d.data))
        O.update_canvas!(c, d.data[i], sh.obskind[]; color = sh.obscolor[],
                         logscale = _obs_logscale(sh))
        # The model prediction too, when the panel is showing it: a data plot saved without
        # the overlay is a different figure from the one on screen.
        sh.obsoverlay[] &&
            Makie.scatter!(c.axis, _model_points(sh, d.data[i], i);
                           color = (:black, 0.0), strokecolor = :black, strokewidth = 1.1,
                           markersize = 9, marker = :circle)
        return fig
    elseif which == "orbit"
        c = build_orbit_canvas(fig)
        te = d === nothing ? Float64[] : d.tepochs
        show_orbit!(c, sh.orbit, te; nside_exp = sh.nside_exp[], T = sh.precision[],
                    binary = _orbit_binary(sh))
        return fig
    elseif which == "posterior"
        f = current_fit(sh.session)
        f === nothing && return "! no fit yet"
        isempty(f.samples) && return "! $(f.method) is an optimiser and has no posterior"
        c = build_post_canvas(fig)
        n = size(f.samples, 2)
        i = clamp(sh.postx[], 1, n); j = clamp(sh.posty[], 1, n)
        show_posterior!(c, f.samples, f.names, i, j, f.best[i], f.q16[i], f.q84[i])
        return fig
    elseif which in ("imaging", "imaging_moll", "imaging_3d")
        im = current_image(sh.session)
        im === nothing && return "! no reconstruction yet"
        got === nothing && return "! the map needs the model's geometry, which does not build"
        star, _ = got
        if which == "imaging"
            c = build_sky_canvas(fig)
            c.cbarlabel[] = "T (K)"
            show_map!(c, star, Float64.(im.x[star.index_quads_visible]);
                      title = im.name, star_params = star_params(m), _decor(sh)...)
        elseif which == "imaging_moll"
            show_mollweide!(build_moll_canvas(fig), Float64.(im.x), star)
        else
            show_star3d!(build_star_canvas(fig), star, Float64.(im.x))
        end
        return fig
    end
    return "! no such plot: " * which
end

"""
    shell_save_figure(which, path, width, height) -> String

Write one plot area to a file. Returns `""` on success, or a message beginning with `!`.

`width`/`height` come from the panel on screen, so the file is framed the way the window frames
it rather than at some size chosen here.

The extension decides the format and Makie decides what it can write: `.png` always works
through GLMakie, and a vector format needs CairoMakie loaded, which ROTIR does not depend on —
so an unwritable extension is reported rather than silently turned into a PNG with the wrong
name.
"""
function shell_save_figure(which, path, width, height)
    sh = _sh()
    p = strip(String(path))
    isempty(p) && return "! no file name"
    occursin('.', basename(p)) || (p = p * ".png")
    # Numbered rather than overwritten: every Save view writes the same name, so saving the
    # orthographic view twice to compare two models used to leave only the second.
    p = unique_path(p)
    sz = (max(320, round(Int, Float64(width))), max(240, round(Int, Float64(height))))

    fig = try
        _snapshot_figure(sh, String(which), sz)
    catch err
        msg = "! could not build the figure: $(sprint(showerror, err))"
        console!(sh, msg); sh.status = msg; return msg
    end
    fig isa AbstractString && (console!(sh, fig); sh.status = fig; return fig)

    try
        Makie.save(p, fig)
    catch err
        msg = "! could not write $(p): $(sprint(showerror, err))"
        console!(sh, msg); sh.status = msg; return msg
    finally
        # CLOSE THE HIDDEN SCREEN. `Makie.save` renders the offscreen figure by attaching a
        # GLMakie `Screen{GLFW.Window}` to its scene and leaving it there, registered in
        # `GLMakie.ALL_SCREENS`. That open GLFW window then keeps the PROCESS alive after the
        # Qt window has closed and `QML.exec()` has returned: measured in the click test as a
        # session that drove every panel correctly, saved a PNG, and then never exited — the
        # driver finished and `wait` hung for the rest of its timeout.
        try
            scr = Makie.getscreen(fig.scene)
            scr === nothing || close(scr)
        catch
        end
    end
    log!(sh.session, "save($(repr(p)), figure)"; note = "save $(which)")
    sh.status = "wrote $(p), $(sz[1])×$(sz[2])"
    console!(sh, sh.status)
    return ""
end
