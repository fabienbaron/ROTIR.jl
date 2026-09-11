// ROTIR — one window, three perspectives over one session.
//
// The context bar, task tray and console sit OUTSIDE the tab stack on purpose: the dataset is
// what the perspectives share, background work must survive a tab switch, and the log records
// across modes. Putting any of them inside a tab would make this three programs sharing a
// title bar.
//
// QML holds no state. Everything is read from, or pushed through, Julia.

import QtQuick
import QtQuick.Window
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import Makie

ApplicationWindow {
    id: win

    // ── a fixed light palette ─────────────────────────────────────────────────
    //
    // Every panel, table and label below paints itself with a fixed light colour, and the
    // Makie canvases are light too. Qt's controls, left alone, follow the SYSTEM theme — so on
    // a dark desktop the window came out mixed: dark chrome around white tables and white
    // plots, with grey-on-grey text that was hard to read.
    //
    // `apply_controls_style!` pins the STYLE for the same reason; this pins the colours the
    // style derives from, which is the half a style choice does not cover. The values are the
    // ones already used throughout this file, so nothing shifts on a light desktop — they only
    // stop a dark one leaking in.
    //
    // If a dark theme is ever wanted, this block plus `style_axis!` for the figures are the two
    // places that would have to change, and the colour literals scattered through the tabs
    // would have to become palette roles first. That is why the Theme box below is one entry
    // and disabled rather than a dropdown that half works.
    palette {
        window:          "#f4f4f4"
        windowText:      "#222222"
        base:            "#ffffff"
        alternateBase:   "#fbfbfb"
        text:            "#222222"
        button:          "#efefef"
        buttonText:      "#222222"
        placeholderText: "#888888"
        light:           "#ffffff"
        midlight:        "#eeeeee"
        mid:             "#dddddd"
        dark:            "#999999"
        shadow:          "#666666"
        highlight:       "#3874d8"
        highlightedText: "#ffffff"
        toolTipBase:     "#ffffe1"
        toolTipText:     "#222222"
        brightText:      "#ffffff"
        link:            "#0645ad"
    }

    // ── one scale factor for all the chrome ───────────────────────────────────
    //
    // Qt already scales `font.pointSize` by the font DPI, so text grows on a HiDPI screen by
    // itself. Pixel quantities do not, so without `dp()` the text outgrows the containers
    // holding it. `qtTextScale` is exactly the ratio Qt is applying to text by itself, so
    // `fontScale` carries only the DEPARTURE from it and the pointSize never gets the screen's
    // DPI twice.
    readonly property real qtTextScale:
        Math.max(0.5, Math.min(4.0, Screen.logicalPixelDensity * 25.4 / 96.0))

    // What the screen actually is, from PHYSICAL pixel density — which comes from the EDID and
    // genuinely differs between machines, where logicalPixelDensity does not (it reads 96 dpi
    // on essentially every Linux screen). Scaling off the logical figure is why one hardcoded
    // factor can fit one machine and be far too large on another.
    readonly property real physicalDpi: Screen.pixelDensity * 25.4

    // Anchored on points judged by eye, not derived: 0.73 on a 92.6 dpi 1920x1080 desktop and,
    // earlier, 1.25 on a laptop panel of roughly 189 dpi. `refDpi`/`refScale`/`dpiExponent` ARE
    // that fit — move them if a screen reads wrong, rather than trying to re-derive a rule.
    //
    // The anchor was 0.875 and is now 0.73, matching OITOOLS: the same judgement, on the same
    // desk, against the same Fusion style, and the two windows are used side by side. It moves
    // the whole curve down by 17%, so the laptop point it predicts is 1.04 rather than the 1.25
    // once judged there; if 1.25 still reads right on that panel, the two anchors want an
    // exponent of 0.75 rather than a different refScale.
    //
    // The exponent falls out as 0.5. Matching physical SIZE would be 1.0 and give 0.61 here,
    // too small: a desktop monitor sits further away than a laptop screen and wants larger
    // elements to subtend the same angle. Viewing distance is not knowable, so the root splits
    // the difference, and it happens to pass through both judgements.
    readonly property real refDpi:      92.6
    readonly property real refScale:    0.73
    readonly property real dpiExponent: 0.5
    readonly property real autoScale:
        Math.max(0.5, Math.min(4.0,
            physicalDpi > 0 ? refScale * Math.pow(physicalDpi / refDpi, dpiExponent)
                            : refScale))

    // ── the window's own font ─────────────────────────────────────────────────
    //
    // Qt resolves a family name through the SYSTEM font database, so naming one only works
    // where the machine has it — stock Windows has no Noto. These two faces come out of the
    // MakieAssets artifact, which the bundle already carries for the plots, and a FontLoader
    // puts them in Qt's database at startup: the same UI font on every platform, nothing extra
    // shipped. Empty when the assets cannot be found, and then the platform's font is used.
    //
    // This is what makes the `dp()` metrics mean one thing: the layout was measured against
    // Noto Sans, and a different family at the same point size is a different width — which on
    // macOS, whose theme font is wider, is what pushed text into the controls beside it.
    readonly property var shippedFontFiles: {
        var s = Julia.shell_ui_font_files()
        return s.length > 0 ? s.split("\n") : []
    }
    FontLoader { id: shippedFontRegular
                 source: win.shippedFontFiles.length > 0 ? win.shippedFontFiles[0] : "" }
    FontLoader { id: shippedFontBold
                 source: win.shippedFontFiles.length > 1 ? win.shippedFontFiles[1] : "" }
    // Bold is loaded into the same family, so `font.bold` picks the real face rather than
    // having Qt smear the regular one.
    readonly property string shippedFontFamily:
        shippedFontRegular.status === FontLoader.Ready ? shippedFontRegular.name : ""

    // The out-of-the-box appearance, in ONE place. Both the initial values below and the
    // settings panel's "Reset to defaults" read it, so the button cannot drift from the
    // declarations — a reset that restored the wrong number would be worse than none.
    //
    // Zero is not a size anywhere here: it is "unset", and the shipped value behind it lives in
    // Julia (`PLOT_SCALE_AT_REF_DPI`, `ZOOM_PER_DETENT`, each plot's own marker size) or is
    // computed from the screen. So the reset sends zero rather than a number, and the boxes are
    // filled from what Julia then reports — repeating 1.19 here would pin today's constant as a
    // user override and survive a change to it.
    readonly property var appearanceDefaults: ({
        uiScaleUser: 0, uiFontFamily: "", baseFontPt: 11, plotScaleUser: 0, markerSizeUser: 0,
        zoomStepUser: 0,
        // The shipped style is Julia's to name (`DEFAULT_CONTROLS_STYLE`), so it is asked for
        // rather than repeated here: a reset has to restore what the window would ship with,
        // not a literal that was true when this line was written.
        controlsStyle: Julia.shell_default_controls_style()
    })

    // Turned by hand in the settings panel; 0 means "no override". It wins over both the
    // startup variable and the screen, because it is the most recent thing the user said.
    property real uiScaleUser: appearanceDefaults.uiScaleUser
    // Empty means the window's own default: the shipped Noto face, or the platform's font if
    // the assets could not be found.
    property string uiFontFamily: appearanceDefaults.uiFontFamily
    // Qt Quick Controls style. Read by `apply_controls_style!` from the settings file BEFORE
    // the QML is loaded, so a change here takes effect at the next launch, not this one.
    property string controlsStyle: appearanceDefaults.controlsStyle
    // What the plot layer was last TOLD, as opposed to what it draws: zero means "computed
    // from the screen". "Save defaults" stores THESE rather than the values in force, so a
    // scale worked out from this monitor's DPI is not pinned onto the next one.
    property real plotScaleUser: appearanceDefaults.plotScaleUser
    // The data-point size in the Data tab's plot; 0 is "each plot's own". Lives in OITOOLS,
    // which owns that canvas, so it is filled from `shell_plot_scale` rather than assumed.
    property real markerSizeUser: appearanceDefaults.markerSizeUser
    // How far one wheel detent zooms a plot. A setting because it depends on the pointing
    // device as much as on taste: a detented wheel, a free-spinning one and a touchpad all
    // deliver different amounts of scroll for the same gesture. Filled from Julia when the
    // panel opens, since it is Julia that zooms.
    property real zoomStepUser: appearanceDefaults.zoomStepUser
    // The model mesh level, owned here so the settings panel and the Model tab agree on it.
    property int nsideExp: 3

    // CLAMPED, like `autoScale` above. Without this a scale of 0.05 — one step up from
    // "auto" in a spin box counting in hundredths — collapsed every widget in the window to a
    // few pixels, with no way back except restarting.
    readonly property real uiScale:
        Math.max(0.5, Math.min(4.0, uiScaleUser > 0     ? uiScaleUser
                                  : uiScaleOverride > 0 ? uiScaleOverride
                                                        : autoScale))
    readonly property real fontScale: uiScale / qtTextScale

    function dp(px)     { return Math.round(px * uiScale) }
    function pt(points) { return points * fontScale }

    property real baseFontPt: appearanceDefaults.baseFontPt
    font.pointSize: pt(baseFontPt)
    // Empty means the window's own default: the shipped Noto face when the assets loaded, and
    // the platform's font when they did not. Naming "Noto Sans" outright was the old default
    // and only worked where the machine happened to have it installed.
    //
    // The last fallback is the family NAME, not "". A FontLoader is ASYNCHRONOUS: until it
    // reports Ready, `shippedFontFamily` is empty, and every control built in that window —
    // which is most of them — resolves "" to Qt's default and keeps it. Measured: the whole
    // window came out in a monospace face while `win.font.family` read "Noto Sans", because
    // the binding settled after the children had already taken their font.
    //
    // Naming it resolves immediately from the system database where Noto is installed, and is
    // the same family the loader then supplies, so nothing shifts when it arrives. Where Noto
    // is absent Qt falls back silently to the platform font, which is the old behaviour.
    //
    // Not `Qt.application.font.family`: under jlqml `Qt.application` is a QQmlApplication with
    // no `font` member, so that expression is a ReferenceError on every re-evaluation.
    font.family: uiFontFamily.length > 0      ? uiFontFamily
               : shippedFontFamily.length > 0 ? shippedFontFamily
                                              : "Noto Sans"

    // Sized from the SCREEN, not through `dp()`. Passing it through the UI scale ties the
    // window to the widgets, and the two want opposite things: asking for smaller chrome
    // should buy MORE room for the plot, not less.
    width:  Math.max(1000, Math.round(Screen.desktopAvailableWidth  * 0.85))
    height: Math.max(640,  Math.round(Screen.desktopAvailableHeight * 0.92))
    visible: true
    title: "ROTIR"

    property string status: initialStatus
    property bool jobRunning: false
    property real jobElapsed: 0
    // What the running engine last reported: evaluations and criterion. Empty when nothing is
    // running, or when the job is one that cannot report (a sampler has no partial answer).
    property string jobProgress: ""
    // Whether anything is loaded, so the picker can offer "Add epoch" only when there is
    // something to add to.
    property bool hasDataset: false

    // ── the polling timers ────────────────────────────────────────────────────
    //
    // Two, at different rates. The job poll is what lets a running engine's output reach the
    // console and, when it finishes, what takes the result up ON THE GUI THREAD — the only
    // place a GL call may happen. 200 ms is fast enough to read as live and slow enough that
    // reading a growing file costs nothing.
    Timer {
        interval: 200
        running: true
        repeat: true
        onTriggered: {
            var r = Julia.shell_job_poll().split("\t")
            var wasRunning = win.jobRunning
            win.jobRunning = r[0] === "1"
            win.jobElapsed = parseFloat(r[1])
            // The fourth field is the running engine's own progress, and reading it is what
            // DRAWS the partial map: `shell_job_poll` puts it on the canvases before
            // returning, on this thread, which is the only place a GL call may happen.
            win.jobProgress = r.length > 3 ? r[3] : ""
            if (win.jobRunning || wasRunning) {
                consolePane.text = Julia.shell_console() +
                                   (r.length > 2 && r[2].length > 0 ? "\n" + r[2] : "")
                if (win.jobProgress.length > 0) imageTab.redraw()
            }
            // The transition from running to not running is where a finished job has just been
            // taken up in Julia, so this is when the canvases have new data in them.
            if (wasRunning && !win.jobRunning) {
                win.status = Julia.shell_status()
                win.refreshAll()
            }
        }
    }
    // The console alone, more slowly: it also changes for things that are not jobs (a load, a
    // validation message), and re-reading a 4000-line pane five times a second is waste.
    Timer {
        interval: 700
        running: true
        repeat: true
        onTriggered: if (!win.jobRunning) consolePane.text = Julia.shell_console()
    }

    // Test hook: an automated run has nobody to close the window.
    Timer {
        interval: autoQuitMs > 0 ? autoQuitMs : 1
        running: autoQuitMs > 0
        repeat: false
        onTriggered: Qt.quit()
    }

    // Two refreshes, and the difference is what makes tab switching instant.
    //
    // `refreshAll()` recomputes: it rebuilds the geometry, the temperature map, the Mollweide
    // resampling and the per-epoch χ². Call it when the DATA changed — a file opened, a
    // parameter edited, a job finished.
    //
    // `repaintAll()` only asks Qt to repaint what the canvases already hold, and re-reads the
    // (cheap, string) tables. Call it when only the VIEW changed. Tab switching is the whole
    // reason it exists: it was calling refreshAll(), which meant a click on a tab header
    // rebuilt every epoch's Fourier setup.
    function refreshAll() {
        win.status = Julia.shell_refresh()
        repaintAll()
    }

    function repaintAll() {
        dataTab.refresh();  dataTab.redraw()
        modelTab.refresh(); modelTab.redraw()
        imageTab.refresh(); imageTab.redraw()
        orbitTab.refresh(); orbitTab.redraw()
    }

    // One tab, for a TAB SWITCH. `repaintAll()` was doing this four times over, and the three
    // invisible ones are not free: `refresh()` rebuilds that tab's ListModels from the shell's
    // tables — the epoch table has one entry per epoch and the parameter form one per field —
    // and a ListModel rebuild is QML-engine work whether or not anyone can see the result.
    //
    // Nothing needs RECOMPUTING either way. The Julia side of all four refreshes measures
    // 0.03 ms together: each tab owns its own Figures and they keep their contents, so a
    // switch has always been a repaint rather than a rebuild. This is about not doing the
    // repaint three extra times.
    //
    // The tabs left behind stay correct because every path that changes the DATA goes through
    // `refreshAll()`, which still refreshes all four.
    // "Save plot…" from any tab: remember WHAT to save, then ask where. The picker returns a
    // path under mode "savefig" and `onAccepted` finishes the job — the write cannot happen
    // before the user has named a file, and the size has to be captured now because the area
    // that knows it is on a tab the picker is about to cover.
    property string saveWhich: ""
    property int saveW: 0
    property int saveH: 0
    function askSaveFigure(which, w, h) {
        win.saveWhich = which; win.saveW = w; win.saveH = h
        picker.purpose = "image"
        picker.canAdd = false
        picker.saveMode = true
        picker.suggestedName = "rotir_" + which + ".png"
        picker.openAt(Julia.image_dir())
    }

    function repaintTab(i) {
        if      (i === 0) { dataTab.refresh();  dataTab.redraw()  }
        else if (i === 1) { modelTab.refresh(); modelTab.redraw() }
        else if (i === 2) { imageTab.refresh(); imageTab.redraw() }
        else if (i === 3) { orbitTab.refresh(); orbitTab.redraw() }
    }

    // Saved appearance defaults, applied BEFORE the first plot is drawn: the plot layer sizes
    // its type from the scale, and applying it afterwards would draw once at the wrong size.
    function applySavedSettings() {
        var txt = Julia.shell_load_settings()
        if (txt.length === 0) return
        var lines = txt.split("\n")
        for (var i = 0; i < lines.length; ++i) {
            var f = lines[i].split("\t")
            if (f.length < 2) continue
            if      (f[0] === "ui_scale")    win.uiScaleUser   = parseFloat(f[1])
            else if (f[0] === "ui_font")     win.uiFontFamily  = f[1]
            else if (f[0] === "controls_style") win.controlsStyle = f[1]
            else if (f[0] === "ui_font_pt")  win.baseFontPt    = parseFloat(f[1])
            else if (f[0] === "marker_size") win.markerSizeUser = parseFloat(f[1])
            else if (f[0] === "plot_scale")  win.plotScaleUser = parseFloat(f[1])
            else if (f[0] === "zoom_step")   { var z = parseFloat(f[1])
                                               if (z > 1) { win.zoomStepUser = z
                                                            zoomSpin.value = Math.round(z * 100)
                                                            Julia.shell_set_zoom_step(z) } }
            else if (f[0] === "precision")   precisionBox.currentIndex =
                                                 (f[1] === "Float64" ? 1 : 0)
        }
        // DEFERRED, every one of them. These setters redraw, and this function runs inside
        // `Component.onCompleted` — before the first frame, with the render thread waiting on
        // the Julia lock. Called directly they deadlock the window black. `Qt.callLater` puts
        // them after the event loop is turning, which is also when a redraw can succeed.
        Qt.callLater(function () {
            if (win.plotScaleUser > 0) Julia.shell_set_plot_scale(win.plotScaleUser)
            if (win.markerSizeUser > 0) Julia.shell_set_marker_size(win.markerSizeUser)
        })
    }

    // Plot scale, marker size and zoom step live in Julia, so ask it rather than trusting the
    // numbers the spin boxes were built with. The first field is the scale in FORCE, which is
    // what the box shows: `plotScaleUser` may be 0 for "computed from the screen", and a box
    // reading 0 would say nothing about what the plots are doing.
    function readPlotSettings() {
        var f = Julia.shell_plot_scale().split("\t")
        if (f.length !== 4) return
        plotScaleSpin.value  = Math.round(parseFloat(f[0]) * 100)
        win.plotScaleUser    = parseFloat(f[1])
        win.markerSizeUser   = parseFloat(f[2])
        markerSpin.value     = Math.round(win.markerSizeUser)
        win.zoomStepUser     = parseFloat(f[3])
        zoomSpin.value       = Math.round(win.zoomStepUser * 100)
    }

    Component.onCompleted: {
        applySavedSettings()
        // QML owns the scale detection because `Screen` is live; Julia needs the value for
        // anything it draws. See src/gui/scaling.jl.
        Julia.shell_ui_scale(uiScale, physicalDpi)
        refreshAll()
    }

    // ── the settings panel ───────────────────────────────────────────────────
    //
    // Appearance only, and saved per USER rather than per project: a scale that suits this
    // monitor suits it for every dataset. Centred and dimmed, because the middle of the window
    // belongs to nothing else and a panel there reads as a dialog rather than as another
    // settings column.
    Popup {
        id: settingsPanel
        // Capped, and clamped, because it can outgrow the window. With no height set a Popup
        // takes its content's, and the content is font-sized: at a large enough UI font the
        // panel ran off the bottom and took the buttons with it — the settings panel being
        // exactly where one goes to undo a setting like that. The content scrolls instead.
        //
        // Against `parent`, NOT `win`: a Popup is positioned inside the window's CONTENT item,
        // which is the window less its header. Measured against `win.height` the panel sits a
        // toolbar's height too low and runs off the bottom even when it would otherwise fit.
        //
        // The width grows with the FONT, not with the content: 430 was chosen at 11 pt, and the
        // rows are text, so a larger font needs proportionally more room or the spin boxes get
        // squeezed. Deriving it from the content's implicit width instead is what Qt reports as
        // a binding loop — the content's width comes back from the panel's, so asking the
        // content how wide the panel should be closes the circle.
        x: Math.max(dp(12), Math.round((parent.width  - width)  / 2))
        y: Math.max(dp(12), Math.round((parent.height - height) / 2))
        width:  Math.min(dp(430) * Math.max(1, baseFontPt / 11), parent.width - dp(24))
        // Height from the CONTENT COLUMN, not from `implicitHeight`. The panel's only child is
        // a ScrollView filling it, and an anchored child contributes nothing to a Popup's
        // implicit size — so `implicitHeight` was padding alone, the panel collapsed to a
        // sliver, and the column drew at its natural size outside it. Asking the column is safe
        // and closes no loop: its implicitHeight depends on its WIDTH, and the width above
        // depends only on the font and the window.
        height: Math.min(settingsColumn.implicitHeight + topPadding + bottomPadding,
                         parent.height - dp(24))
        modal: true
        dim: true
        focus: true
        closePolicy: Popup.CloseOnEscape | Popup.CloseOnPressOutsideParent
        padding: dp(12)

        // The panel is capped at the window's height, so its content has to be able to scroll —
        // and to clip, since a control drawn outside the panel is worse than one below the fold.
        ScrollView {
            id: settingsScroll
            anchors.fill: parent
            clip: true
            rightPadding: dp(12)                   // room for the vertical scrollbar
            // The content is exactly as wide as the view: the rows are a label/control grid
            // that stretches, so there is nothing to scroll to sideways, and any binding that
            // lets the content decide the width closes a loop back through the panel's width.
            contentWidth: availableWidth
            ScrollBar.horizontal.policy: ScrollBar.AlwaysOff

        ColumnLayout {
            id: settingsColumn
            width: settingsScroll.availableWidth
            spacing: dp(10)

            RowLayout {
                Layout.fillWidth: true
                Label { text: "Appearance"; font.bold: true }
                Item { Layout.fillWidth: true }
                Label {
                    text: "plot font and controls style need a restart"
                    color: "#888"
                    font.pointSize: pt(baseFontPt - 2)
                }
            }

            GridLayout {
                columns: 3
                Layout.fillWidth: true
                columnSpacing: dp(8)
                rowSpacing: dp(6)

                Label { text: "UI scale"; color: "#666" }
                SpinBox {
                    id: uiScaleSpin
                    // 0 is "auto"; anything else starts at 50, because below that the window
                    // is unreadable and there is no way back to this panel from inside it.
                    from: 0; to: 400; stepSize: 5
                    // A DOUBLE validator, because the box DISPLAYS "1.10" and an IntValidator
                    // then refuses the "." the user has to type to edit it — the field looked
                    // broken for every value that was not a whole number.
                    validator: DoubleValidator {
                        bottom: 0.0; top: 4.0; decimals: 2; notation: DoubleValidator.StandardNotation
                    }
                    // The value IN FORCE, not the override. `uiScaleUser` is 0 for "work it
                    // out from the screen", and a box reading 0 — or "auto" — says nothing
                    // about how large the window actually is. Going back to automatic is the
                    // button beside it, which is also the only thing that needs to express 0.
                    value: Math.round(win.uiScale * 100)
                    editable: true
                    textFromValue: function (v) { return (v / 100).toFixed(2) }
                    valueFromText: function (s) { return Math.round(parseFloat(s) * 100) }
                    Layout.fillWidth: true
                    // Everything sized through dp()/pt() rebinds, so the window resizes as the
                    // number changes rather than on close.
                    onValueModified: win.uiScaleUser = value / 100
                }
                Button {
                    text: "auto"
                    enabled: win.uiScaleUser > 0
                    ToolTip.visible: hovered
                    ToolTip.text: uiScaleOverride > 0
                        ? "back to ROTIRGUI_SCALE=" + uiScaleOverride.toFixed(2)
                        : "back to the value computed from " + win.physicalDpi.toFixed(0) + " dpi"
                    onClicked: win.uiScaleUser = 0
                }

                Label { text: "Controls style"; color: "#666" }
                ComboBox {
                    id: controlsStyleBox
                    Layout.fillWidth: true
                    // The list comes from Julia rather than being repeated here: a second copy
                    // would be free to drift from CONTROLS_STYLES, and offering a style the
                    // bundled Qt does not carry would fail at the next launch rather than at
                    // the click. Line 1 is the style in force; the rest are the choices.
                    property var styleInfo: Julia.shell_controls_styles().split("\n")
                    model: styleInfo.slice(1)
                    // Line 1 is the style the window is actually RUNNING, which is what the box
                    // has to show: after a change it stays the previous choice until the next
                    // launch, and a box claiming otherwise would misreport the window.
                    currentIndex: Math.max(0, model.indexOf(styleInfo[0]))
                    onActivated: win.controlsStyle = currentText
                    ToolTip.visible: hovered
                    ToolTip.text: "Fusion is what this window was laid out against; Basic is " +
                                  "Qt's own and draws its controls larger. Applies at the next launch."
                }
                Label {
                    text: win.controlsStyle === controlsStyleBox.styleInfo[0] ? "" : "on restart"
                    color: "#888"; font.pointSize: pt(baseFontPt - 2)
                }

                Label { text: "UI font"; color: "#666" }
                ComboBox {
                    id: fontBox
                    Layout.fillWidth: true
                    // Index 0 is "whatever the window would use if you had never set this",
                    // and it NAMES that font rather than saying "default": which one it is
                    // depends on whether the shipped face loaded. When it did, a separate
                    // "Noto Sans" entry would be the same font twice, so the list carries the
                    // platform's own font at the end instead.
                    model: win.shippedFontFamily.length > 0
                         ? ["Noto Sans", "DejaVu Sans", "Liberation Sans", "JuliaMono",
                            "(system font)"]
                         : ["(system default)", "DejaVu Sans", "Noto Sans",
                            "Liberation Sans", "JuliaMono"]
                    currentIndex: Math.max(0, model.indexOf(win.uiFontFamily))
                    // "" is the only way back to the window's own default, and naming a family
                    // that is not installed falls back silently rather than reporting anything
                    // — so a user who has no Noto needs an explicit escape.
                    onActivated: win.uiFontFamily =
                        currentIndex === 0              ? ""
                      : currentText === "(system font)" ? ""
                                                        : currentText
                }
                Item {}

                // Its own row, so the box lands in the SAME column as every other control.
                // Sharing the font row put it in column 3, out of line with the spin boxes and
                // combos above and below it.
                Label { text: "UI font size"; color: "#666" }
                SpinBox {
                    id: uiFontSizeSpin
                    from: 6; to: 24; stepSize: 1
                    value: Math.round(win.baseFontPt)
                    editable: true
                    Layout.fillWidth: true
                    onValueModified: win.baseFontPt = value
                }
                Label { text: "pt"; color: "#888"; font.pointSize: pt(baseFontPt - 2) }

                Label { text: "Plot font"; color: "#666" }
                SpinBox {
                    id: plotScaleSpin
                    from: 0; to: 400; stepSize: 5
                    value: Math.round(win.plotScaleUser * 100)
                    editable: true
                    textFromValue: function (v) { return v === 0 ? "auto" : (v / 100).toFixed(2) }
                    valueFromText: function (s) {
                        return s === "auto" ? 0 : Math.round(parseFloat(s) * 100)
                    }
                    onValueModified: {
                        win.plotScaleUser = value / 100
                        win.status = Julia.shell_set_plot_scale(win.plotScaleUser)
                        win.refreshAll()
                    }
                }
                Label {
                    text: "0 = from the screen"
                    color: "#888"; font.pointSize: pt(baseFontPt - 2)
                }

                // ── plot symbols ──────────────────────────────────────────────
                // The Data tab's plot is OITOOLS' canvas, so this is forwarded there rather
                // than held here — see `marker_size_user` in shell.jl.
                Label { text: "Plot symbols"; color: "#666" }
                SpinBox {
                    id: markerSpin
                    from: 0; to: 30; stepSize: 1
                    value: 0
                    editable: true
                    Layout.fillWidth: true
                    // 0 means "whatever the plot chooses", which differs per view — uv coverage
                    // draws smaller points than an observable plot.
                    textFromValue: function (v) { return v === 0 ? "auto" : String(v) }
                    valueFromText: function (s) { return s === "auto" ? 0 : parseInt(s) }
                    onValueModified: {
                        win.markerSizeUser = value
                        win.status = Julia.shell_set_marker_size(value)
                    }
                }
                Label { text: "px"; color: "#888"; font.pointSize: pt(baseFontPt - 2) }

                // ── wheel zoom ────────────────────────────────────────────────
                // How far one notch goes — a matter of hardware as much as taste, since a
                // detented wheel, a free-spinning one and a touchpad all report different
                // amounts of scroll for the same intent.
                Label { text: "Wheel zoom"; color: "#666" }
                SpinBox {
                    id: zoomSpin
                    from: 102; to: 300; stepSize: 5
                    value: Math.round(win.zoomStepUser * 100)
                    editable: true
                    Layout.fillWidth: true
                    textFromValue: function (v) { return (v / 100).toFixed(2) + "\u00d7" }
                    valueFromText: function (s) { return Math.round(parseFloat(s) * 100) }
                    onValueModified: {
                        win.zoomStepUser = value / 100
                        win.status = Julia.shell_set_zoom_step(win.zoomStepUser)
                    }
                }
                Label {
                    text: "per detent"
                    color: "#888"; font.pointSize: pt(baseFontPt - 2)
                }

                Label { text: "Precision"; color: "#666" }
                ComboBox {
                    id: precisionBox
                    model: ["Float32", "Float64"]
                    font.pointSize: pt(baseFontPt - 1)
                    ToolTip.text: "the float type the mesh and the polygon FT are built in. " +
                                  "Float32 halves the memory and is what the tessellation " +
                                  "defaults to; the analytic shape fit uses Float64 regardless."
                    ToolTip.visible: hovered
                    // The LEVEL comes from the Model tab, which owns the mesh; this panel
                    // owns only the float type it is built in.
                    onActivated: {
                        win.status = Julia.shell_set_tessellation(
                            "healpix", win.nsideExp, currentText)
                        win.refreshAll()
                    }
                }
                Item {}

                // One entry, and it is honest: the window PINS a light palette at the top of
                // this file, and every panel, table and Makie canvas below is drawn light. A
                // dark mode is real work — the colour literals scattered through the tabs
                // would have to become palette roles first — not a dropdown. Listed rather
                // than omitted so that "why is there no dark mode" has an answer on screen.
                Label { text: "Theme"; color: "#666" }
                ComboBox {
                    model: ["Light"]
                    enabled: false
                    ToolTip.text: "the window pins a light palette and the panels are drawn " +
                                  "against it; Dark is not styled yet"
                    ToolTip.visible: hovered
                }
                Item {}
            }

            RowLayout {
                Layout.fillWidth: true
                // The buttons get the row to themselves. What happened is reported on the line
                // BELOW: sharing the row, the message was elided to a middle-truncated path
                // sitting beside the buttons, which reads as a stray fragment rather than as
                // the answer to what the button just did — and it squeezed the buttons besides.
                Button {
                    text: "Save config"
                    Layout.fillWidth: true
                    ToolTip.visible: hovered
                    ToolTip.text: "write these settings to the per-user config file, which " +
                                  "the window reads at every launch"
                    onClicked: {
                        // The OVERRIDES, not the values in force: a scale of 0 means "work it
                        // out from the screen", which is the right thing to carry to a machine
                        // with a different one.
                        var path = Julia.shell_save_settings(
                            [ "ui_scale\t"       + win.uiScaleUser,
                              "ui_font\t"        + win.uiFontFamily,
                              "controls_style\t" + win.controlsStyle,
                              "ui_font_pt\t"     + win.baseFontPt,
                              "plot_scale\t"     + win.plotScaleUser,
                              "marker_size\t"    + win.markerSizeUser,
                              "zoom_step\t"      + win.zoomStepUser,
                              "precision\t"  + precisionBox.currentText ].join("\n"))
                        savedLabel.text = path.length > 0 ? "saved to " + path
                                                          : "could not save — see the console"
                    }
                }
                Button {
                    text: "Reset to defaults"
                    Layout.fillWidth: true
                    ToolTip.visible: hovered
                    ToolTip.text: "back to the out-of-the-box appearance, and delete the saved config"
                    onClicked: {
                        // Both halves, and the file matters more than the window: settings are
                        // applied at startup, so restoring the look while leaving the file in
                        // place would come back tweaked at the next launch.
                        var d = win.appearanceDefaults
                        win.uiScaleUser   = d.uiScaleUser
                        win.uiFontFamily  = d.uiFontFamily
                        win.baseFontPt    = d.baseFontPt
                        win.controlsStyle = d.controlsStyle
                        uiScaleSpin.value = 0
                        fontBox.currentIndex = 0
                        uiFontSizeSpin.value = d.baseFontPt
                        controlsStyleBox.currentIndex =
                            Math.max(0, controlsStyleBox.model.indexOf(d.controlsStyle))
                        // The plot side is Julia's, and zero means "work it out from the
                        // screen" rather than "zero" — so the setters go first and the boxes
                        // are read back from what Julia then computed, not from d.
                        Julia.shell_set_plot_scale(d.plotScaleUser)
                        Julia.shell_set_marker_size(d.markerSizeUser)
                        Julia.shell_set_zoom_step(d.zoomStepUser)
                        readPlotSettings()
                        win.refreshAll()
                        var removed = Julia.shell_reset_settings()
                        savedLabel.text = removed.length > 0 ? "reset · removed " + removed
                                                             : "reset · nothing had been saved"
                    }
                }
                Button { text: "Close"; Layout.fillWidth: true
                         onClicked: settingsPanel.close() }
            }

            // Always exactly one line. Empty until a button was pressed, it made the panel grow
            // the first time anything was saved; wrapped, it grew again on a long path. And it
            // says what it IS, because a truncated path on its own reads as a stray fragment
            // rather than as the answer to what the button just did.
            Label {
                id: savedLabel
                Layout.fillWidth: true
                elide: Text.ElideMiddle
                maximumLineCount: 1
                text: "config: " + Julia.shell_settings_path()
                color: "#888"
                font.pointSize: pt(baseFontPt - 2)
                HoverHandler { id: savedHover }
                ToolTip.visible: savedHover.hovered && savedLabel.truncated
                ToolTip.text: savedLabel.text
            }

            // Which code is actually running. The first thing to establish about any bug
            // report, and the settings panel is where a user already comes to look at what the
            // window is doing rather than at their data.
            Label {
                Layout.fillWidth: true
                horizontalAlignment: Text.AlignRight
                text: Julia.shell_version()
                color: "#888"
                font.pointSize: pt(baseFontPt - 2)
            }
        }
        }
    }

    FilePicker {
        id: picker
        width: Math.min(win.width * 0.8, dp(940))
        height: Math.min(win.height * 0.8, dp(560))
        anchors.centerIn: Overlay.overlay
        fontPt: pt(10)
        onAccepted: function (paths, mode) {
            if (mode === "orbit") {
                win.status = Julia.shell_load_orbit(paths.split("\n")[0])
                win.refreshAll()
                return
            }
            if (mode === "savefig") {
                var path = paths.split("\n")[0]
                var msg = Julia.shell_save_figure(win.saveWhich, path, win.saveW, win.saveH)
                win.status = msg.length > 0 ? msg : "saved " + path
                return
            }
            if (mode === "map") {
                win.status = Julia.shell_load_map(paths.split("\n")[0])
                win.refreshAll()
                return
            }
            // `shell_open_many` handles one file too, and reads a whole set together —
            // `readoifits_multiepochs` in one call rather than six, which is what makes the
            // epoch origin come out right and the log one line instead of six.
            win.status = Julia.shell_open_many(paths, mode === "add" ? "1" : "0",
                                               mode === "single" ? "0" : "1")
            win.hasDataset = true
            win.refreshAll()
        }
    }

    ColumnLayout {
        anchors.fill: parent
        anchors.margins: dp(6)
        spacing: dp(6)

        // ── context bar: shared by every perspective ─────────────────────────
        RowLayout {
            Layout.fillWidth: true
            spacing: dp(6)
            // Settings first: it governs the whole window rather than the dataset.
            Button {
                text: "\u2699"                    // GEAR
                // Sized well above the bar's text, as in the OITOOLS window. This is the only
                // control in the bar that is a SYMBOL rather than a word, and at text size it
                // read as punctuation beside "Open OIFITS…" instead of as a button.
                font.pointSize: pt(baseFontPt + 11)
                // The glyph is taller than the default content box, so give it room rather
                // than letting the bar clip its top and bottom.
                implicitWidth: dp(38)
                implicitHeight: dp(38)
                ToolTip.text: "appearance settings"
                ToolTip.visible: hovered
                onClicked: {
                    if (settingsPanel.opened) { settingsPanel.close(); return }
                    // Read the live values back, so the panel opens showing what is in force
                    // rather than what it last displayed.
                    uiScaleSpin.value = Math.round(win.uiScale * 100)
                    fontBox.currentIndex = Math.max(0, fontBox.model.indexOf(win.uiFontFamily))
                    uiFontSizeSpin.value = Math.round(win.baseFontPt)
                    readPlotSettings()
                    settingsPanel.open()
                }
            }
            // ONE open button. What happens to the file — new dataset, new dataset taken
            // whole, or another epoch of the current one — is chosen in the picker, next to
            // the file it applies to. Two top-bar buttons named the same operation twice and
            // neither name said what the difference was.
            Button {
                text: "Open OIFITS…"
                font.pointSize: pt(10)
                enabled: !win.jobRunning
                onClicked: {
                    picker.purpose = "data"
                    picker.canAdd = win.hasDataset
                    picker.openAt(initialFolder)
                }
            }
            Button {
                text: "Export script…"
                font.pointSize: pt(10)
                enabled: !win.jobRunning
                // Writes next to the working directory rather than opening a save dialog: the
                // picker is a file CHOOSER, and a chooser that has to invent a name for a file
                // that does not exist yet is a different widget. One less thing to build for a
                // one-line action.
                onClicked: win.status = Julia.shell_export("rotirgui_session.jl")
            }
            // The map itself, beside the script that would rebuild it. The script says how the
            // session got here; the FITS says what it arrived at, and carries the HEALPix level
            // and every parameter, so the χ² can be recomputed from the file alone.
            Button {
                text: "Save model map"
                font.pointSize: pt(10)
                enabled: !win.jobRunning
                ToolTip.text: "the reconstructed map if there is one, otherwise the model's " +
                              "own, with the tessellation level and every parameter"
                ToolTip.visible: hovered
                onClicked: win.status = Julia.shell_save_map("rotir_map.fits")
            }
            Button {
                text: "Load model map"
                font.pointSize: pt(10)
                enabled: !win.jobRunning
                onClicked: {
                    picker.purpose = "map"
                    picker.canAdd = false
                    picker.saveMode = false
                    picker.openAt(initialFolder)
                }
            }
            Item { Layout.fillWidth: true }
            BusyIndicator {
                running: win.jobRunning
                visible: win.jobRunning
                implicitWidth: dp(20); implicitHeight: dp(20)
            }
            Label {
                // The engine's own count while it runs, the outcome once it has stopped. A
                // spinner alone says a job exists; this says whether it is getting anywhere.
                text: win.jobRunning && win.jobProgress.length > 0
                      ? win.jobProgress + "   (" + win.jobElapsed.toFixed(0) + " s)"
                      : win.status
                elide: Text.ElideMiddle
                Layout.maximumWidth: win.width * 0.5
                font.pointSize: pt(10)
            }
        }

        TabBar {
            id: tabBar
            Layout.fillWidth: true
            currentIndex: initialTab
            TabButton { text: "Data";    font.pointSize: pt(10) }
            TabButton { text: "Model";   font.pointSize: pt(10) }
            TabButton { text: "Imaging"; font.pointSize: pt(10) }
            TabButton { text: "Orbit";   font.pointSize: pt(10) }
        }

        StackLayout {
            Layout.fillWidth: true
            Layout.fillHeight: true
            currentIndex: tabBar.currentIndex
            // Repaint only. Makie draws on demand, so a MakieArea that was not visible while
            // its Observables were reassigned still holds the previous frame and needs an
            // update() — but nothing needs recomputing, because nothing changed.
            onCurrentIndexChanged: win.repaintTab(currentIndex)

            DataTab  { id: dataTab; fontFamily: win.uiFontFamily;  uiScale: win.uiScale; fontPt: pt(10)
                       onStatusChanged: function (s) { if (s.length > 0) win.status = s }
                       onSaveRequested: function (which, w, h) { win.askSaveFigure(which, w, h) }
                       onRefreshAllRequested: win.refreshAll() }
            ModelTab { id: modelTab; fontFamily: win.uiFontFamily;
                       precision: precisionBox.currentText; uiScale: win.uiScale; fontPt: pt(10)
                       jobRunning: win.jobRunning; jobElapsed: win.jobElapsed
                       jobProgress: win.jobProgress
                       onStatusChanged: function (s) { if (s.length > 0) win.status = s }
                       onSaveRequested: function (which, w, h) { win.askSaveFigure(which, w, h) } }
            ImageTab { id: imageTab; fontFamily: win.uiFontFamily; uiScale: win.uiScale; fontPt: pt(10)
                       jobRunning: win.jobRunning; jobElapsed: win.jobElapsed
                       onStatusChanged: function (s) { if (s.length > 0) win.status = s }
                       onSaveRequested: function (which, w, h) { win.askSaveFigure(which, w, h) } }
            // Orbit last: it is the one perspective that is about the SYSTEM rather than
            // about one star's surface, and it does not depend on the other three.
            OrbitTab { id: orbitTab; fontFamily: win.uiFontFamily
                       uiScale: win.uiScale; fontPt: pt(10)
                       jobRunning: win.jobRunning
                       onStatusChanged: function (s) { if (s.length > 0) win.status = s }
                       onSaveRequested: function (which, w, h) { win.askSaveFigure(which, w, h) }
                       onPickFile: function (mode) {
                           picker.purpose = mode
                           picker.canAdd = false
                    picker.saveMode = false
                           // Orbits have a folder of their own, seeded with the ones that
                           // ship — so "Load orbit…" opens on β Lyr and Spica rather than on
                           // whatever directory Julia was started from.
                           picker.openAt(mode === "orbit" ? Julia.orbit_dir() : initialFolder)
                       } }
        }

        OutputConsole {
            id: consolePane
            Layout.fillWidth: true
            // Collapsed, it keeps only its own header row; the tabs above take the rest.
            Layout.preferredHeight: expanded ? dp(150) : dp(30)
            fontPt: pt(9)
        }
    }
}
