# UI scaling policy.
#
# Adapted from OITOOLS' src/gui/scaling.jl, and for the same reason: the window mixes two
# sizing systems that scale by different rules.
#
#   font.pointSize   points      Qt converts using the font DPI, so text already grows on a
#                                HiDPI screen without anyone asking it to
#   width/spacing    logical px  scaled by the device pixel ratio only, so containers do NOT
#                                grow with the text they hold
#
# Left unscaled, text outgrows the things that contain it. So the pixel quantities are scaled
# by the same factor Qt applies to the text.
#
# DETECTION IS NOT HERE. `Screen.logicalPixelDensity` is a live QML property — correct before
# any window exists, and updated when the window is dragged to a monitor with a different
# scale, neither of which a value computed once in Julia could manage. QML works the factor
# out; this file supplies only the override, which keeps the policy (parsing, range, what
# counts as "auto") in plain Julia where it is testable without a display.
#
# Plot typography on the STATIC path is NOT scaled: `plot2d_makie` and friends keep the sizes
# in `src/oiplot_spheroid_makie.jl` so a figure saved from a script is the figure in the paper.
# The live canvas follows the UI scale, because it is read on screen rather than on paper.

"Environment variable that overrides the detected UI scale."
const UI_SCALE_VAR = "ROTIRGUI_SCALE"

"Sentinel handed to QML meaning \"work it out from the screen\"."
const UI_SCALE_AUTO = 0.0

# Below 0.5 the window is unreadable; above 4.0 it cannot fit a 4K panel. A value outside the
# range is far more likely to be a typo (a stray percentage) than an intention, so it is
# clamped and reported rather than obeyed.
const UI_SCALE_MIN = 0.5
const UI_SCALE_MAX = 4.0

"""
    ui_scale_override(; verbose = true) -> (; scale, auto, reason)

Read `ROTIRGUI_SCALE` and turn it into the value QML expects.

`scale` is `UI_SCALE_AUTO` (0.0) when QML should detect the scale itself, which is the default
and the case for every input that does not name a usable number. `reason` is the one-line
explanation shown at startup.

    ROTIRGUI_SCALE=1.5     # everything 50% larger than the screen asks for
    ROTIRGUI_SCALE=auto    # or unset, or empty: let QML decide

Out-of-range values are clamped with a warning; anything unparseable warns and falls back to
automatic. Neither is an error — a bad scale must not stop the GUI opening, because the GUI is
how you would fix it.
"""
function ui_scale_override(; verbose::Bool = true)
    raw = strip(get(ENV, UI_SCALE_VAR, ""))
    auto(reason) = (scale = UI_SCALE_AUTO, auto = true, reason = reason)

    isempty(raw)                    && return auto("unset; QML scales from the screen")
    lowercase(raw) in ("auto", "0") && return auto("$UI_SCALE_VAR=$raw; QML scales from the screen")

    parsed = tryparse(Float64, raw)
    if parsed === nothing || !isfinite(parsed)
        verbose && @warn "$UI_SCALE_VAR is not a number, ignoring it" value = raw
        return auto("$UI_SCALE_VAR=$raw is not a number; falling back to automatic")
    end
    if parsed < UI_SCALE_MIN || parsed > UI_SCALE_MAX
        clamped = clamp(parsed, UI_SCALE_MIN, UI_SCALE_MAX)
        verbose && @warn "$UI_SCALE_VAR is outside the usable range, clamping" requested =
            parsed used = clamped range = (UI_SCALE_MIN, UI_SCALE_MAX)
        return (scale = clamped, auto = false,
                reason = "$UI_SCALE_VAR=$parsed clamped to $clamped")
    end
    return (scale = parsed, auto = false, reason = "$UI_SCALE_VAR=$parsed")
end

"""
The UI scale QML settled on, pushed back from the window.

QML owns the detection, so the value has to travel back for anything drawn in Julia to match
it. Until the window reports in this holds a desktop-ish default: the canvas is built before
any window exists, and every later refresh restyles from the live value anyway.
"""
const UI_SCALE_LIVE = Ref(1.0)

"Physical screen density in dpi as QML read it; 0 until the window reports in."
const SCREEN_DPI_LIVE = Ref(0.0)

# Plot typography is sized from the DENSITY rather than from the UI scale, so the two can be
# retuned separately. Tying them together means every adjustment to one moves the other.
const PLOT_REF_DPI          = 92.6
const PLOT_SCALE_AT_REF_DPI = 1.19
const PLOT_SCALE_AT_REFERENCE = 1.19

"""
    set_ui_scale!(x, dpi = 0.0) -> Float64

Record what QML measured. Called once, from the window's `Component.onCompleted`.

`dpi` may be 0 where the screen does not report one — a virtual display, say — in which case
plot sizing falls back to the UI scale.
"""
function set_ui_scale!(x::Real, dpi::Real = 0.0)
    isfinite(dpi) && dpi > 0 && (SCREEN_DPI_LIVE[] = Float64(dpi))
    (isfinite(x) && x > 0) || return UI_SCALE_LIVE[]
    return UI_SCALE_LIVE[] = clamp(Float64(x), UI_SCALE_MIN, UI_SCALE_MAX)
end

"Live override from the settings panel. Zero means \"no override, compute it\"."
const PLOT_SCALE_USER = Ref(0.0)

"Set the plot scale by hand. Zero restores the value computed from the screen."
function set_plot_scale!(x::Real)
    PLOT_SCALE_USER[] = (isfinite(x) && x > 0) ? clamp(Float64(x), 0.2, 5.0) : 0.0
    return PLOT_SCALE_USER[]
end

"""
    live_plot_scale() -> Float64

The factor the live canvas multiplies its font sizes by.

Follows the square root of screen density, anchored so a `PLOT_REF_DPI` desktop gets
`PLOT_SCALE_AT_REF_DPI`. Falls back to the UI scale when the screen reports no density, so
plot text still tracks the chrome rather than staying fixed beside it.
"""
function live_plot_scale()
    PLOT_SCALE_USER[] > 0 && return PLOT_SCALE_USER[]
    dpi = SCREEN_DPI_LIVE[]
    dpi > 0 && return PLOT_SCALE_AT_REF_DPI * sqrt(dpi / PLOT_REF_DPI)
    return PLOT_SCALE_AT_REFERENCE * UI_SCALE_LIVE[]
end
