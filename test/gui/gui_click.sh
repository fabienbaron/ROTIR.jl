#!/usr/bin/env bash
# Click-driven GUI test: drives the real widgets with xdotool on a virtual screen.
#
#     test/gui/gui_click.sh              # headless, via Xvfb
#     test/gui/gui_click.sh --display    # against $DISPLAY, e.g. your own machine
#
# WHY THIS EXISTS SEPARATELY FROM test/gui/runtests.jl.
#
# That file calls the shell callbacks directly. It covers everything reachable from Julia and it
# PASSES on a machine where the GUI is visibly broken — because it never opens a POPUP. Under
# QMLMakie a popup gets its own GL context, and from then on inserting a plot allocates buffers
# with none bound:
#
#     gl_renderobject: glGenBuffers returned invalid id. OpenGL Context active?
#
# ROTIR is built to survive that (every plot is created before the window and only Observables
# are assigned afterwards — see src/gui/livecanvas.jl), and this script is what proves it:
# it opens the file picker, the surface-type combo and the fit-method combo, and then checks
# that everything still draws. It reproduces on llvmpipe, so no GPU is needed.
#
# Clicks are placed as FRACTIONS of the window, not fixed pixels: a fixed pixel silently stops
# landing on its control the moment the default window size or the UI scale changes, and a
# missed click shows up as a missing console line rather than as an obvious failure. The scale
# is pinned to 1.0 for the same reason.
set -uo pipefail
cd "$(dirname "$0")/../.."          # repo root: --project=bin needs it

DUMP="$(mktemp)"; JLOG="$(mktemp)"
SHOTDIR="${ROTIRGUI_SHOTS:-$(mktemp -d)}"; mkdir -p "$SHOTDIR"
# Two λ And nights, symlinked into a directory of their own, so the picker listing is
# deterministic: the driver selects rows by ORDER, and demos/data has fifteen entries whose
# order would decide which file the test loads.
PICKDIR="$(mktemp -d)"
ln -sf "$PWD/demos/data/2011Sep02.lam_And_prepped.oifits" "$PICKDIR/a_first.oifits"
ln -sf "$PWD/demos/data/2011Sep06.lam_And_prepped.oifits" "$PICKDIR/b_second.oifits"
FAIL=0; GUI_PID=""; XVFB_PID=""

if [ "${1:-}" = "--display" ]; then
    DISP="${DISPLAY:?--display needs DISPLAY set}"
else
    DISP=":99"
fi

cleanup() {
    [ -n "$GUI_PID" ]  && kill "$GUI_PID"  2>/dev/null
    [ -n "$XVFB_PID" ] && kill "$XVFB_PID" 2>/dev/null
    rm -rf "$PICKDIR"
    return 0
}
trap cleanup EXIT

say() { printf '  [%s] %s\n' "$(date +%H:%M:%S)" "$*"; }

for tool in xdotool import xwininfo; do
    command -v "$tool" >/dev/null || { say "SKIP  $tool is not installed"; exit 77; }
done

say "display      $DISP"
say "julia log    $JLOG"
say "console dump $DUMP"
say "screenshots  $SHOTDIR"
say "picker dir   $PICKDIR"
if [ "$DISP" != ":99" ]; then
    say "NOTE: on a real display this moves your pointer and clicks. Keep hands off for ~2 min."
fi

if [ "$DISP" = ":99" ]; then
    # A display already on :99 is FATAL, not something to attach to. `Xvfb` fails silently when
    # the display is taken, the driver then talks to whatever was already there, and every
    # coordinate below is a fraction of a window that is now the wrong size — which surfaces as
    # six unrelated checks failing with "clicks missed?" and sends you looking for a regression
    # that is not there. Cost an hour once; the guard is one line.
    if [ -e "/tmp/.X11-unix/X99" ] || xdpyinfo -display :99 >/dev/null 2>&1; then
        say "FAIL  :99 is already in use — kill the stray X server first:"
        pgrep -a Xvfb | sed 's/^/        /'
        exit 1
    fi
    # +extension GLX and +render: llvmpipe needs both, and without them the window never maps.
    Xvfb "$DISP" -screen 0 1600x1000x24 +extension GLX +render >/dev/null 2>&1 &
    XVFB_PID=$!
    sleep 3
fi

say "precompiling (the slow part; silent otherwise) ..."
if ! DISPLAY="$DISP" julia --project=bin -e 'using Pkg; Pkg.precompile()' > "$JLOG" 2>&1; then
    say "FAIL  precompilation failed:"; tail -25 "$JLOG"; exit 1
fi
say "precompiled; starting the GUI"

ROTIRGUI_DATA_DIR="$PICKDIR" ROTIRGUI_SCALE=1.0 ROTIRGUI_CONSOLE_DUMP="$DUMP" \
DISPLAY="$DISP" LIBGL_ALWAYS_SOFTWARE=1 QT_QPA_PLATFORM="${QT_QPA_PLATFORM:-xcb}" \
  julia --project=bin -e '
      # The order is forced: Mesa reads its driver variables and GLFW its platform hint when
      # the first GL context is created, so both must happen ABOVE `using GLMakie`.
      using ROTIR
      using OITOOLS: configure_graphics!, configure_qt_platform!, prefer_native_wayland!
      configure_graphics!()
      using GLFW_jll
      wl = prefer_native_wayland!()
      configure_qt_platform!(; match_x11 = !wl.applied)
      using GLMakie, QMLMakie, QML
      @info "packages loaded; opening the window"
      gui(; autoquit_ms = 420000)' > "$JLOG" 2>&1 &
GUI_PID=$!

WID=""
for i in $(seq 1 300); do
    WID=$(DISPLAY="$DISP" xdotool search --name "^ROTIR$" 2>/dev/null | head -1)
    [ -n "$WID" ] && break
    if ! kill -0 "$GUI_PID" 2>/dev/null; then
        say "FAIL  the GUI exited before a window appeared:"; tail -30 "$JLOG"; exit 1
    fi
    [ $((i % 20)) -eq 0 ] && say "waiting for the window (${i}s) ... $(tail -1 "$JLOG" | cut -c1-80)"
    sleep 1
done
if [ -z "$WID" ]; then
    say "FAIL  no window titled ROTIR after 300s."
    if [ -n "${WAYLAND_DISPLAY:-}" ]; then
        say "      WAYLAND_DISPLAY is set, so Qt may have opened a WAYLAND surface, which"
        say "      xdotool cannot see or click. Re-run with QT_QPA_PLATFORM=xcb."
    fi
    say "      xdotool sees $(DISPLAY="$DISP" xdotool search --name . 2>/dev/null | wc -l) X windows"
    tail -30 "$JLOG"; exit 1
fi
say "window $WID found; driving it"
sleep 4

eval "$(DISPLAY="$DISP" xdotool getwindowgeometry --shell "$WID")"
say "geometry ${WIDTH}x${HEIGHT}"
# The fractions below were read off a 1360x920 window. They survive a modest change, but a
# window half that size puts the picker's own buttons off screen, so say so here rather than
# reporting it as six missed clicks.
if [ "$WIDTH" -lt 1100 ] || [ "$HEIGHT" -lt 700 ]; then
    say "FAIL  the window is ${WIDTH}x${HEIGHT}; the click fractions assume at least 1100x700."
    exit 1
fi

# Fractions of the window, resolved once. Read off a screenshot at ROTIRGUI_SCALE=1.0; if the
# layout changes they must be re-read, and a missed click shows up as a failed check below.
FX() { echo $((WIDTH  * $1 / 1000)); }
FY() { echo $((HEIGHT * $1 / 1000)); }
CL()   { DISPLAY="$DISP" xdotool mousemove --window "$WID" "$1" "$2" click 1; sleep "${3:-1.5}"; }
# `--window`, not a bare `key`: Xvfb runs with no window manager, so there is no input
# focus and a bare key goes to whatever is under the pointer.
KEY()  { DISPLAY="$DISP" xdotool key --window "$WID" "$1"; sleep "${2:-1}"; }
SHOT() { DISPLAY="$DISP" import -window "$WID" "$SHOTDIR/$1.png" 2>/dev/null; }

GEAR_X=$(FX 33);     OPEN_X=$(FX 102);   TOP_Y=$(FY 22)
# Inside the picker. Measured off 02_loaded.png at 1360x920, ROTIRGUI_SCALE=1.0.
ROW1_X=$(FX 515);    ROW1_Y=$(FY 276)
ROW2_X=$(FX 515);    ROW2_Y=$(FY 325)
ACTION_Y=$(FY 784)
OPENBTN_X=$(FX 812)      # rightmost button when there is no dataset yet
ADDBTN_X=$(FX 748)       # "Add epoch", once one exists (the row is one button wider)
# Four tabs share the bar, so each centre is an eighth, three eighths, ... of the width.
TAB_DATA=$(FX 128);  TAB_MODEL=$(FX 376); TAB_IMAGING=$(FX 623); TAB_ORBIT=$(FX 871)
TAB_Y=$(FY 55)
sleep 8       # the first frame. xdotool finds the window before Qt has drawn into it,
              # and a shot taken then is a blank 290-byte PNG.
SHOT 00_start

# ── the file picker: opening it is what takes the GL context ─────────────────
CL "$OPEN_X" "$TOP_Y" 3               # "Open OIFITS…"
SHOT 01_picker
# CLICKED, not typed. The picker's keyboard navigation is a Shortcut, and under an Xvfb with
# no window manager there is no input focus for `xdotool key` to target reliably — a synthetic
# key event there is not a fair test of the shortcut, and a driver that cannot tell a broken
# shortcut from a WM-less display is worse than one that clicks.
CL "$ROW1_X" "$ROW1_Y" 2              # select the first file
CL "$OPENBTN_X" "$ACTION_Y" 8         # "Open"
SHOT 02_loaded

# ── Ctrl-click: the modifier path, which was dead until the Qt linter found it ──
# `TapHandler`'s signal argument is a `QEventPoint` and has no `modifiers`, so the picker read
# `undefined` and every click took the plain-click branch. Nothing on screen said so — the row
# highlighted, the selection just did not extend. Two files ctrl-clicked must open as two.
CL "$OPEN_X" "$TOP_Y" 3               # "Open OIFITS…"
CL "$ROW1_X" "$ROW1_Y" 1              # plain click selects the first
DISPLAY="$DISP" xdotool keydown ctrl
CL "$ROW2_X" "$ROW2_Y" 1              # ctrl-click ADDS the second
DISPLAY="$DISP" xdotool keyup ctrl
SHOT 03a_ctrl_click
CL "$(FX 748)" "$ACTION_Y" 8          # the button now reads "Open 2"
SHOT 03b_ctrl_opened

# ── a second night, added to the same dataset ────────────────────────────────
# One Open button now; WHAT happens to the file is chosen in the picker, so this opens it
# again, moves to the second file and presses the picker's own "Add epoch".
CL "$OPEN_X" "$TOP_Y" 3               # "Open OIFITS…"
CL "$ROW2_X" "$ROW2_Y" 2              # the SECOND file
SHOT 02b_second_row
# "Add epoch" appears only once a dataset exists, so the action row is one button wider here
# than it was above and every button in it has shifted left.
CL "$ADDBTN_X" "$ACTION_Y" 8          # "Add epoch"
SHOT 03_two_epochs

# ── the Model tab: a combo POPUP, then a model, then every view ──────────────
CL "$TAB_MODEL" "$TAB_Y" 2
SHOT 04_model_tab
CL "$(FX 132)" "$(FY 113)" 2          # surface-type combo -> popup (second GL surface)
SHOT 05_type_popup
KEY Down 1
KEY Return 2
CL "$(FX 294)" "$(FY 113)" 5          # "+ model"
SHOT 06_model_added

CL "$(FX 491)" "$(FY 89)" 3           # Mollweide
SHOT 07_view_mollweide
CL "$(FX 554)" "$(FY 89)" 4           # Full 3D
SHOT 08_view_3d
CL "$(FX 618)" "$(FY 89)" 3           # posterior — a fourth MakieArea, and a fourth
SHOT 08b_view_posterior             #   GL surface after three popups
CL "$(FX 425)" "$(FY 89)" 3           # back to orthographic
SHOT 09_view_ortho
# Saving REBUILDS the panel offscreen and writes a PNG — it never touches the live framebuffer,
# which under QMLMakie hands back noise (see src/gui/snapshot.jl). Done HERE, on a view that
# has a map on it: from the posterior view with no fit yet the panel correctly refuses, which
# is right behaviour and a useless assertion.
CL "$(FX 966)" "$(FY 88)" 5           # Save view…
SHOT 09a_saved
# Zoom out with the wheel, then reset with a RIGHT-CLICK on the plot — the gesture replaced a
# button, so this is the only way to exercise it.
DISPLAY="$DISP" xdotool mousemove --window "$WID" "$(FX 700)" "$(FY 480)" click 5 click 5
sleep 2; SHOT 09b_zoomed
DISPLAY="$DISP" xdotool mousemove --window "$WID" "$(FX 700)" "$(FY 480)" click 3
sleep 2; SHOT 09c_reset
CL "$(FX 402)" "$(FY 153)" 3          # the intensity tick
SHOT 10a_intensity
CL "$(FX 557)" "$(FY 153)" 3          # the graticules tick
SHOT 10b_graticules
CL "$(FX 557)" "$(FY 120)" 3          # the epoch ▶ arrow
SHOT 10c_next_epoch
CL "$(FX 586)" "$(FY 185)" 3          # a colormap button
SHOT 10_colormap

# "− model" clears it: the surface views must go back to their idle state rather than keeping
# the last picture, which would be a model that is no longer there.
CL "$(FX 357)" "$(FY 113)" 3          # "− model"
SHOT 10d_model_cleared
# ...and back, because the reconstruction below needs a geometry to sit on.
CL "$(FX 294)" "$(FY 113)" 5          # "+ model"

# ── the Imaging tab: tick a regulariser and run one ──────────────────────────
CL "$TAB_IMAGING" "$TAB_Y" 2
SHOT 11_imaging_tab
# The context header — which dataset and which model the run is about — sits above the
# regulariser list, so every row here is ~50px lower than it was before it existed.
# Re-read after the kernel selector was added below the mesh row, which pushed everything
# under it down by one row.
CL "$(FX 23)" "$(FY 226)" 1           # tick sobel2
CL "$(FX 23)" "$(FY 414)" 1           # tick radflat  (Monnier's)
SHOT 12_regularisers
CL "$(FX 38)" "$(FY 677)" 30          # Reconstruct, then let it run
SHOT 13_reconstructed

# ── the settings panel: one more popup, one more GL surface ──────────────────
CL "$GEAR_X" "$TOP_Y" 3               # the gear
SHOT 14_settings
KEY Escape 2
SHOT 15_settings_closed

# ── the Orbit tab: the elements, and the star-model selector below them ──────
CL "$TAB_ORBIT" "$TAB_Y" 3
SHOT 18a_orbit_tab
# Free a fixed element: its bounds appear in the trailing slot, and the columns to its left
# must not move. The screenshot is the check — the alignment is not observable from Julia.
CL "$(FX 188)" "$(FY 265)" 2          # the "e" row's state combo
KEY Down 1
KEY Return 2
SHOT 18b_orbit_free
# The star model: analytic 2-D profiles or ROTIR's tessellated 3-D components. Another popup,
# and the one place the two frames of this tab are chosen between.
CL "$(FX 230)" "$(FY 475)" 2
SHOT 18c_star_model_popup
KEY Down 1
KEY Return 3
SHOT 18d_tessellated

# ── back to Data: everything must still draw after all those popups ──────────
CL "$TAB_DATA" "$TAB_Y" 3
SHOT 16_back_to_data
# The observable views come first now, so the model views are further along the row.
CL "$(FX 330)" "$(FY 89)" 3           # first observable panel
SHOT 17_observable
# The epoch is ONE thing shared by every tab, and the Model tab's arrow moved it to 2 above.
# The table has to be on 2 as well — it used to restore its own remembered row, so it read
# "1 / 2" over a plot of night two. The console line is what says which one Julia is on.
CL "$(FX 405)" "$(FY 156)" 3          # the Data tab's own epoch arrow, back to 1
SHOT 17b_epoch_back
CL "$(FX 400)" "$(FY 89)" 3           # the next one
SHOT 18_observable2

wait "$GUI_PID" 2>/dev/null

echo "-- console transcript --"
sed 's/^/  /' "$DUMP"

echo "-- checks --"
ck() {  # ck <description> <0|1 ok>
    if [ "$2" = "1" ]; then echo "  PASS  $1"; else echo "  FAIL  $1"; FAIL=1; fi
}

grep -qF "ROTIRGUI ready" "$DUMP" && ck "the GUI started" 1 || ck "the GUI started" 0
grep -qF "load:" "$DUMP"          && ck "the picker loaded a file" 1 \
                                  || ck "the picker loaded a file (clicks missed?)" 0
grep -qF "add:" "$DUMP"           && ck "Add epoch appended a second night" 1 \
                                  || ck "Add epoch appended a second night" 0
grep -qE "model .*surface_type" "$DUMP" && ck "a model was created" 1 \
                                        || ck "a model was created" 0
grep -qF "cleared model" "$DUMP" && ck "− model cleared it" 1 || ck "− model cleared it" 0
grep -qF "image_reconstruct_oi" "$DUMP" && ck "a reconstruction was started" 1 \
                                        || ck "a reconstruction was started" 0
grep -qE "χ²ᵣ = " "$DUMP"         && ck "the reconstruction finished and reported χ²" 1 \
                                  || ck "the reconstruction finished and reported χ²" 0
grep -qF "star model: tessellated" "$DUMP" && ck "the 3-D star model was selected" 1 \
                                           || ck "the 3-D star model was selected" 0
grep -qE "wrote .*\.png" "$DUMP" && ck "Save view wrote a PNG" 1 \
                                 || ck "Save view wrote a PNG" 0
grep -qE "now has 2 epoch\(s\)" "$DUMP" && ck "Ctrl-click selected two files at once" 1 \
                                        || ck "Ctrl-click selected two files at once" 0

# THE regression check, and the reason the script exists: anything drawn after a popup must
# still draw.
if grep -qE "glGenBuffers|gl_renderobject" "$DUMP" "$JLOG"; then
    ck "no GL context loss after the popups (the known QMLMakie bug)" 0
else
    ck "no GL context loss after the popups" 1
fi
if grep -qF "exception in render" "$JLOG"; then
    ck "no swallowed render exception" 0
else
    ck "no swallowed render exception" 1
fi
if grep -qE "Qt Warning.*(unavailable|ReferenceError|Duplicate signal)" "$JLOG"; then
    echo "  FAIL  QML errors:"; grep -E "Qt Warning" "$JLOG" | sort -u | head -5; FAIL=1
else
    ck "no QML load errors" 1
fi
grep -qE "failed|could not" "$DUMP" && { ck "no failure lines in the transcript" 0; \
    grep -E "failed|could not" "$DUMP" | head -3; } || ck "no failure lines in the transcript" 1

# A screenshot of a window that never mapped is a few hundred bytes; a real one is ~100 kB.
SMALL=$(find "$SHOTDIR" -name '*.png' -size -20k | wc -l)
[ "$SMALL" -eq 0 ] && ck "every screenshot captured a mapped window" 1 \
                   || ck "$SMALL screenshot(s) captured a blank window" 0

echo "  screenshots: $SHOTDIR"
exit $FAIL
