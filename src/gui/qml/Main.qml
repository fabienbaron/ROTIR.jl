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

    // Anchored on judgement rather than derived: 0.875 on a 92.6 dpi desktop. The exponent is
    // 0.5 — matching physical SIZE would be 1.0 and gives elements that are too small, because
    // a desktop monitor sits further away than a laptop panel and wants larger elements to
    // subtend the same angle.
    readonly property real refDpi:      92.6
    readonly property real refScale:    0.875
    readonly property real dpiExponent: 0.5
    readonly property real autoScale:
        Math.max(0.5, Math.min(4.0,
            physicalDpi > 0 ? refScale * Math.pow(physicalDpi / refDpi, dpiExponent)
                            : refScale))

    // Turned by hand in the settings panel; 0 means "no override". It wins over both the
    // startup variable and the screen, because it is the most recent thing the user said.
    property real uiScaleUser: 0
    property string uiFontFamily: ""
    // What the plot layer was last TOLD, as opposed to what it draws: zero means "computed
    // from the screen". "Save defaults" stores THESE rather than the values in force, so a
    // scale worked out from this monitor's DPI is not pinned onto the next one.
    property real plotScaleUser: 0
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

    property real baseFontPt: 11
    font.pointSize: pt(baseFontPt)
    // Empty means whatever the platform theme chose, which is the sane default.
    font.family: uiFontFamily.length > 0 ? uiFontFamily : ""

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
            else if (f[0] === "ui_font_pt")  win.baseFontPt    = parseFloat(f[1])
            else if (f[0] === "plot_scale")  win.plotScaleUser = parseFloat(f[1])
            else if (f[0] === "precision")   precisionBox.currentIndex =
                                                 (f[1] === "Float64" ? 1 : 0)
        }
        if (win.plotScaleUser > 0) Julia.shell_set_plot_scale(win.plotScaleUser)
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
        x: Math.round((win.width  - width)  / 2)
        y: Math.round((win.height - height) / 2)
        width: dp(430)
        modal: true
        dim: true
        focus: true
        closePolicy: Popup.CloseOnEscape | Popup.CloseOnPressOutsideParent
        padding: dp(12)

        ColumnLayout {
            width: parent.width
            spacing: dp(10)

            RowLayout {
                Layout.fillWidth: true
                Label { text: "Appearance"; font.bold: true }
                Item { Layout.fillWidth: true }
                Label {
                    text: "the plot font needs a restart"
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
                    value: Math.round(win.uiScaleUser * 100)
                    editable: true
                    textFromValue: function (v) { return v === 0 ? "auto" : (v / 100).toFixed(2) }
                    valueFromText: function (s) {
                        return s === "auto" ? 0 : Math.round(parseFloat(s) * 100)
                    }
                    onValueModified: {
                        // Step over the unusable range in one go rather than walking through
                        // it: 5 % is not a scale anyone wants, and passing through it makes
                        // the panel unclickable on the way.
                        if (value > 0 && value < 50) value = (value < uiScaleSpin.stepSize * 2) ? 0 : 50
                        win.uiScaleUser = value / 100
                    }
                }
                Label {
                    text: "0 = auto (" + win.autoScale.toFixed(2) + " here)"
                    color: "#888"; font.pointSize: pt(baseFontPt - 2)
                }

                Label { text: "UI font"; color: "#666" }
                ComboBox {
                    id: fontBox
                    Layout.fillWidth: true
                    // "" first: the platform theme's own font is the right default, and
                    // naming a family that is not installed silently falls back to something
                    // else rather than reporting anything.
                    model: ["", "Sans Serif", "DejaVu Sans", "Noto Sans", "Liberation Sans",
                            "Cantarell", "Ubuntu", "monospace"]
                    currentIndex: Math.max(0, model.indexOf(win.uiFontFamily))
                    displayText: currentText.length === 0 ? "system default" : currentText
                    onActivated: win.uiFontFamily = currentText
                }
                SpinBox {
                    from: 6; to: 24; value: Math.round(win.baseFontPt)
                    onValueModified: win.baseFontPt = value
                }

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

                Label { text: "Theme"; color: "#666" }
                ComboBox {
                    model: ["Light"]
                    // One entry, and it is honest: the panels hardcode light colours, so a
                    // Dark option would leave half the window unreadable. Listed rather than
                    // omitted so that "why is there no dark mode" has an answer on screen.
                    enabled: false
                    ToolTip.text: "the panels hardcode light colours; Dark is not styled yet"
                    ToolTip.visible: hovered
                }
                Item {}
            }

            RowLayout {
                Layout.fillWidth: true
                Label {
                    id: savedLabel
                    Layout.fillWidth: true
                    elide: Text.ElideMiddle
                    color: "#888"
                    font.pointSize: pt(baseFontPt - 2)
                }
                Button {
                    text: "Save defaults"
                    ToolTip.visible: hovered
                    ToolTip.text: "write these as the startup defaults for every project"
                    onClicked: {
                        // The OVERRIDES, not the values in force: a scale of 0 means "work it
                        // out from the screen", which is the right thing to carry to a machine
                        // with a different one.
                        var path = Julia.shell_save_settings(
                            [ "ui_scale\t"   + win.uiScaleUser,
                              "ui_font\t"    + win.uiFontFamily,
                              "ui_font_pt\t" + win.baseFontPt,
                              "plot_scale\t" + win.plotScaleUser,
                              "precision\t"  + precisionBox.currentText ].join("\n"))
                        savedLabel.text = path.length > 0 ? "saved to " + path
                                                          : "could not save — see the console"
                    }
                }
                Button { text: "Close"; onClicked: settingsPanel.close() }
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
                font.pointSize: pt(11)
                ToolTip.text: "appearance settings"
                ToolTip.visible: hovered
                onClicked: settingsPanel.opened ? settingsPanel.close() : settingsPanel.open()
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
                       onRefreshAllRequested: win.refreshAll() }
            ModelTab { id: modelTab; fontFamily: win.uiFontFamily;
                       precision: precisionBox.currentText; uiScale: win.uiScale; fontPt: pt(10)
                       jobRunning: win.jobRunning; jobElapsed: win.jobElapsed
                       jobProgress: win.jobProgress
                       onStatusChanged: function (s) { if (s.length > 0) win.status = s } }
            ImageTab { id: imageTab; fontFamily: win.uiFontFamily; uiScale: win.uiScale; fontPt: pt(10)
                       jobRunning: win.jobRunning; jobElapsed: win.jobElapsed
                       onStatusChanged: function (s) { if (s.length > 0) win.status = s } }
            // Orbit last: it is the one perspective that is about the SYSTEM rather than
            // about one star's surface, and it does not depend on the other three.
            OrbitTab { id: orbitTab; fontFamily: win.uiFontFamily
                       uiScale: win.uiScale; fontPt: pt(10)
                       jobRunning: win.jobRunning
                       onStatusChanged: function (s) { if (s.length > 0) win.status = s }
                       onPickFile: function (mode) {
                           picker.purpose = mode
                           picker.canAdd = false
                           picker.openAt(initialFolder)
                       } }
        }

        OutputConsole {
            id: consolePane
            Layout.fillWidth: true
            Layout.preferredHeight: dp(150)
            fontPt: pt(9)
        }
    }
}
