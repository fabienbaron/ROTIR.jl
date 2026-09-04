// Data: what was observed, and how well the current model reproduces it.
//
// Deliberately NOT a generic uv / V² / T3 browser — OITOOLS has those and does them better.
// What ROTIR adds is the EPOCH structure: a surface map is constrained by how the star looks
// at several rotation phases, so the questions this tab answers are "which phase is this",
// "what does the model look like there", and "which epoch does the model fit badly".

import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml            // supplies `Julia`
import Makie            // MakieArea

// A Pane, not an Item: Qt Quick Controls inherit their font from the nearest CONTROL
// ancestor, and a plain Item has no font — so a font chosen in the settings panel reached the
// panel itself (a direct child of the window) and nothing else. `background: null` keeps it
// invisible; it exists only to carry the font down.
Pane {
    id: root
    padding: 0
    background: null
    // "" means INHERIT, and it is what happened anyway: this read
    // `Qt.application.font.family`, but QML.jl's application object is a QQmlApplication,
    // which has no `font` — so the binding raised a TypeError and left the family unset.
    font.family: root.fontFamily.length > 0 ? root.fontFamily : ""
    font.pointSize: root.fontPt
    property string fontFamily: ""
    property real uiScale: 1
    property real fontPt: 10
    function dp(px) { return Math.round(px * uiScale) }

    signal statusChanged(string s)
    // Removing an epoch or a dataset changes what every tab shows, so the window refreshes
    // all of them rather than this tab redrawing itself in isolation.
    signal refreshAllRequested()

    // Whether the three observable columns hold reduced χ² or point counts. Julia decides —
    // there is a χ² only when a model is set — and the header says which, rather than putting
    // "χ²ᵣ" above a column of counts.
    property bool showingChi2: false
    // The observable views this file actually has, then the two model views. Built from
    // OITOOLS' own availability check, so the tab cannot offer a panel for a table that is
    // not in the file.
    property var obsKinds: []
    property var views: []
    property int nObs: 0
    property int viewIndex: 0

    // What the observable canvas should show. Pushed to Julia, which owns the drawing.
    function applyView() {
        if (root.viewIndex < root.nObs) {
            root.statusChanged(Julia.shell_set_obs_view(
                root.obsKinds[root.viewIndex], colorBox.currentText,
                overlayBox.checked ? "1" : "0",
                logBox.enabled && logBox.checked ? "1" : "0",
                maxBaseBox.enabled && maxBaseBox.checked ? "1" : "0"))
        }
        root.redraw()
    }

    // Labels for the observable keys. Kept here rather than in Julia because they are purely
    // presentational — the KEY is what crosses the boundary.
    function obsLabel(k) {
        return k === "uv"     ? "uv coverage"
             : k === "v2"     ? "V²"
             : k === "t3amp"  ? "T3 amp"
             : k === "t3phi"  ? "T3 phase"
             : k === "visamp" ? "vis amp"
             : k === "visphi" ? "vis phase"
             : k === "flux"   ? "flux"
             : k === "diffphi" ? "diff. phase"
             : k === "diffvisamp" ? "diff. vis amp"
             : k
    }

    ListModel { id: epochModel }
    ListModel { id: datasetModel }

    function refresh() {
        // Which observable panels this dataset supports, from OITOOLS' availability flags.
        if (root.obsKinds.length === 0) {
            var flags = Julia.shell_obs_kinds()
            var ks = []
            if (flags.length > 0) {
                var fl = flags.split("\n")
                for (var j = 0; j < fl.length; ++j) {
                    var g = fl[j].split("\t")
                    if (g.length >= 2 && g[1] === "1") ks.push(g[0])
                }
            }
            if (ks.length > 0) {
                root.obsKinds = ks
                root.nObs = ks.length
                var vs = []
                for (var q = 0; q < ks.length; ++q) vs.push(root.obsLabel(ks[q]))
                vs.push("sky map"); vs.push("per-epoch χ²")
                root.views = vs
                root.applyView()
            }
        }
        var drows = Julia.shell_datasets()
        var dkeep = dsBox.currentIndex
        datasetModel.clear()
        if (drows.length > 0) {
            var dlines = drows.split("\n")
            for (var k = 0; k < dlines.length; ++k) {
                var g = dlines[k].split("\t")
                if (g.length < 4) continue
                datasetModel.append({ dname: g[0], nep: g[1], nv2t: g[2], nt3t: g[3] })
            }
        }
        // Follow Julia rather than QML's own memory: `shell_open` can create OR extend a
        // dataset, and which of the two happened is Julia's decision.
        var cur = parseInt(Julia.shell_current_dataset())
        if (cur >= 1 && cur <= datasetModel.count) dsBox.currentIndex = cur - 1
        else if (dkeep >= 0 && dkeep < datasetModel.count) dsBox.currentIndex = dkeep

        var rows = Julia.shell_epochs()
        var keep = epochList.currentIndex
        epochModel.clear()
        if (rows.length > 0) {
            var lines = rows.split("\n")
            for (var i = 0; i < lines.length; ++i) {
                var f = lines[i].split("\t")
                if (f.length < 7) continue
                epochModel.append({ idx: f[0], mjd: f[1], tday: f[2], v2: f[3],
                                    t3amp: f[4], t3phi: f[5], haveChi2: f[6] === "1" })
                root.showingChi2 = (f[6] === "1")
            }
        }
        // FOLLOW JULIA, exactly as the dataset box above does. "Add epoch" moves the session
        // to the epoch it appended, and a table that restored its own previous row instead
        // highlighted epoch 1 while the plot beside it drew epoch 2. Only when Julia has no
        // opinion does the remembered row win — this runs after every parameter edit, and
        // jumping back to epoch 1 on each keystroke makes the tab unusable.
        var ce = parseInt(Julia.shell_current_epoch())
        if (ce >= 1 && ce <= epochModel.count) epochList.currentIndex = ce - 1
        else if (keep >= 0 && keep < epochModel.count) epochList.currentIndex = keep
    }

    RowLayout {
        anchors.fill: parent
        spacing: dp(8)

        // ── left: the epoch table ────────────────────────────────────────────
        ColumnLayout {
            // fillWidth is FALSE explicitly: a nested layout defaults it to true,
            // and this column would then take the whole row and squeeze the canvas
            // beside it to zero width — which reads as a plot that failed to draw.
            Layout.fillWidth: false
            Layout.preferredWidth: dp(450)
            Layout.minimumWidth: dp(450)
            Layout.maximumWidth: dp(450)
            Layout.fillHeight: true
            spacing: dp(4)

            Label { text: "Dataset"; font.bold: true; font.pointSize: root.fontPt }
            ComboBox {
                id: dsBox
                Layout.fillWidth: true
                model: datasetModel
                font.pointSize: root.fontPt
                displayText: currentIndex < 0 ? "no dataset" :
                    (datasetModel.get(currentIndex).dname + "  —  " +
                     datasetModel.get(currentIndex).nep + " epoch(s), " +
                     datasetModel.get(currentIndex).nv2t + " V², " +
                     datasetModel.get(currentIndex).nt3t + " T3")
                delegate: ItemDelegate {
                    width: dsBox.width
                    text: dname + "  (" + nep + " ep)"
                    font.pointSize: root.fontPt
                }
                onActivated: {
                    root.statusChanged(Julia.shell_select_dataset(currentIndex + 1))
                    root.refresh()
                }
            }
            RowLayout {
                Layout.fillWidth: true
                spacing: dp(6)
                Button {
                    text: "Remove epoch"
                    enabled: epochList.currentIndex >= 0 && epochModel.count > 0
                    font.pointSize: root.fontPt - 1
                    ToolTip.text: "drop the selected epoch; removing the last one closes the dataset"
                    ToolTip.visible: hovered
                    onClicked: {
                        root.statusChanged(Julia.shell_remove_epoch(epochList.currentIndex + 1))
                        root.refreshAllRequested()
                    }
                }
                Button {
                    text: "Close dataset"
                    enabled: datasetModel.count > 0
                    font.pointSize: root.fontPt - 1
                    onClicked: {
                        root.statusChanged(Julia.shell_close_dataset())
                        root.refreshAllRequested()
                    }
                }
                Item { Layout.fillWidth: true }
            }

            Label { text: "Epochs"; font.bold: true; font.pointSize: root.fontPt }

            // A header row of its own rather than a real table header: this is six fixed
            // columns, and a TableView with a HorizontalHeaderView costs far more QML than it
            // buys at that size.
            //
            // INSIDE the Frame, and above the ListView, because that is the only way it takes
            // the same inset the rows take. Outside it, the header started at the panel edge
            // while every row started at the Frame's padding plus the ItemDelegate's — a
            // couple of dozen pixels of style-dependent offset that no width in this file
            // could have compensated for, since neither number is fixed here.
            Frame {
                Layout.fillWidth: true
                Layout.fillHeight: true
                ColumnLayout {
                    anchors.fill: parent
                    spacing: dp(2)
                    RowLayout {
                        Layout.fillWidth: true
                        // The rows lose their right edge to the vertical scrollbar and the
                        // header does not, so the last column drifts by its width as soon as
                        // the list overflows — which for epochs is the normal case.
                        Layout.rightMargin: epochScroll.visible ? epochScroll.width : 0
                        spacing: dp(6)
                        Label { text: "#";     Layout.preferredWidth: dp(24); font.pointSize: root.fontPt - 1; color: "#7f8c98" }
                        Label { text: "MJD";   Layout.fillWidth: true;        font.pointSize: root.fontPt - 1; color: "#7f8c98" }
                        Label { text: "t (d)"; Layout.preferredWidth: dp(52); font.pointSize: root.fontPt - 1; color: "#7f8c98"; horizontalAlignment: Text.AlignRight }
                        Label { text: "V²";    Layout.preferredWidth: dp(54); font.pointSize: root.fontPt - 1; color: "#7f8c98"; horizontalAlignment: Text.AlignRight }
                        Label { text: "T3amp"; Layout.preferredWidth: dp(54); font.pointSize: root.fontPt - 1; color: "#7f8c98"; horizontalAlignment: Text.AlignRight }
                        Label { text: "T3phi"; Layout.preferredWidth: dp(54); font.pointSize: root.fontPt - 1; color: "#7f8c98"; horizontalAlignment: Text.AlignRight }
                    }
                ListView {
                    id: epochList
                    Layout.fillWidth: true
                    Layout.fillHeight: true
                    clip: true
                    model: epochModel
                    ScrollBar.vertical: ScrollBar { id: epochScroll }
                    onCurrentIndexChanged: {
                        // The REPAINT is the second half of this and it is not optional:
                        // `shell_select_epoch` refills the Observables, but Makie draws on
                        // demand and Qt repaints only when asked, so without this the V²
                        // curve stays on the previous epoch while every number beside it
                        // has already moved.
                        if (currentIndex >= 0) {
                            root.statusChanged(Julia.shell_select_epoch(currentIndex + 1))
                            root.redraw()
                        }
                    }
                    delegate: ItemDelegate {
                        width: epochList.width
                        // Zero horizontal padding: the header above has none either, and the
                        // default inset is what put the columns out of line with it.
                        leftPadding: 0
                        rightPadding: 0
                        highlighted: ListView.isCurrentItem
                        onClicked: epochList.currentIndex = index
                        contentItem: RowLayout {
                            spacing: dp(6)
                            Label { text: idx;  Layout.preferredWidth: dp(24); font.pointSize: root.fontPt }
                            Label { text: mjd;  Layout.fillWidth: true;        font.pointSize: root.fontPt }
                            Label { text: tday; Layout.preferredWidth: dp(52); font.pointSize: root.fontPt; horizontalAlignment: Text.AlignRight; color: "#5a6772" }
                            // One column per observable. A model can fit V² well and closure
                            // phase badly, and a single combined number hides exactly that.
                            Repeater {
                                model: [v2, t3amp, t3phi]
                                Label {
                                    text: modelData
                                    Layout.preferredWidth: dp(54)
                                    font.pointSize: root.fontPt
                                    horizontalAlignment: Text.AlignRight
                                    // Red past 3 is a judgement, not a threshold in the code:
                                    // it marks the cells worth looking at first. Counts stay
                                    // neutral — a large one is not a bad one.
                                    color: !haveChi2 ? "#5a6772"
                                         : modelData === "—" ? "#7f8c98"
                                         : parseFloat(modelData) > 3 ? "#c0392b" : "#1e824c"
                                }
                            }
                        }
                    }
                }
                }
            }

            Label {
                Layout.fillWidth: true
                font.pointSize: root.fontPt - 1
                color: "#7f8c98"
                visible: epochModel.count > 0
                text: root.showingChi2 ? "V² / T3amp / T3phi columns: reduced χ² per epoch"
                                       : "V² / T3amp / T3phi columns: point counts (no model yet)"
            }

            Label {
                Layout.fillWidth: true
                wrapMode: Text.WordWrap
                font.pointSize: root.fontPt - 1
                color: "#7f8c98"
                text: epochModel.count === 0
                      ? "No dataset. Open an OIFITS — one file is split on gaps in its V² " +
                        "timestamps; several files become several epochs of one dataset."
                      : "t is days since the first epoch, which is what the rotation phase uses."
            }
        }

        // ── right: one view at a time ────────────────────────────────────────
        //
        // The DATA is what this tab is for, so the observable views come first: V² against
        // baseline, closure phase, uv coverage — the OITOOLS GUI's own plots, drawn by its own
        // canvas (see `build_obs_canvas`). The model's prediction goes OVER them, which is the
        // thing ROTIR adds and the reason the tab is not simply OITOOLS' Explore.
        //
        // The sky projection and the per-epoch χ² are the other two views. They are about the
        // model rather than the data, which is why they are not the default here.
        ColumnLayout {
            Layout.fillWidth: true
            Layout.fillHeight: true
            spacing: dp(4)

            RowLayout {
                Layout.fillWidth: true
                spacing: dp(6)
                Repeater {
                    model: root.views
                    Button {
                        text: modelData
                        checkable: true
                        checked: root.viewIndex === index
                        font.pointSize: root.fontPt
                        onClicked: { root.viewIndex = index; root.applyView() }
                    }
                }
                Button {
                    text: "Save view…"
                    font.pointSize: root.fontPt - 1
                    onClicked: {
                        // Rebuilt offscreen, not grabbed: see src/gui/snapshot.jl.
                        var w = root.viewIndex < root.nObs ? "obs"
                              : root.viewIndex === root.nObs ? "sky" : "chi2"
                        var m = Julia.shell_save_figure(w, "rotir_" + w + ".png",
                                                        obsArea.width, obsArea.height)
                        root.statusChanged(m.length > 0 ? m : "saved rotir_" + w + ".png")
                    }
                }
                Item { Layout.fillWidth: true }
                Label {
                    text: root.viewIndex < root.nObs
                          ? "the selected epoch · model prediction over the data"
                          : root.viewIndex === root.nObs
                            ? "the model at the selected epoch · scroll to zoom · right-click resets"
                            : "reduced χ² per epoch: V² / T3amp / T3phi"
                    color: "#7f8c98"
                    font.pointSize: root.fontPt - 1
                }
            }

            RowLayout {
                Layout.fillWidth: true
                spacing: dp(6)
                visible: root.viewIndex < root.nObs
                Label { text: "colour by"; color: "#7f8c98"; font.pointSize: root.fontPt - 1 }
                ComboBox {
                    id: colorBox
                    Layout.preferredWidth: dp(130)
                    model: ["baseline", "wav", "mjd", "station"]
                    font.pointSize: root.fontPt - 1
                    onActivated: root.applyView()
                }
                CheckBox {
                    id: overlayBox
                    text: "model prediction"
                    checked: true
                    font.pointSize: root.fontPt - 1
                    ToolTip.text: "black rings: the current model at this epoch"
                    ToolTip.visible: hovered
                    onToggled: root.applyView()
                }
                CheckBox {
                    id: logBox
                    text: "log y"
                    // Only where it means something. V² and T3amp are positive quantities;
                    // a closure phase takes both signs, and uv coverage is a geometry.
                    enabled: ["v2", "t3amp"].indexOf(root.obsKinds[root.viewIndex]) >= 0
                    font.pointSize: root.fontPt - 1
                    ToolTip.text: "points at or below zero are dropped, not piled on the floor"
                    ToolTip.visible: hovered
                    onToggled: root.applyView()
                }
                CheckBox {
                    id: maxBaseBox
                    text: "max baseline"
                    // Only the closure quantities have a choice to make: a triangle has three
                    // legs, so it has both a mean baseline and a longest one. V² has a single
                    // baseline and uv coverage is a geometry.
                    enabled: ["t3amp", "t3phi"].indexOf(root.obsKinds[root.viewIndex]) >= 0
                    font.pointSize: root.fontPt - 1
                    ToolTip.text: "plot against the LONGEST leg of each triangle instead of " +
                                  "the geometric mean — the resolution the closure actually " +
                                  "probes"
                    ToolTip.visible: hovered
                    onToggled: root.applyView()
                }
                Item { Layout.fillWidth: true }
            }

            // The same intensity tick as the Model tab, for the sky-map view.
            RowLayout {
                Layout.fillWidth: true
                spacing: dp(6)
                visible: root.viewIndex === root.nObs
                CheckBox {
                    id: dataIntensity
                    // Ticked by default, matching the Model and Imaging tabs and Julia's own
                    // default: the map is a temperature, but the intensity is what is measured.
                    checked: true
                    text: "intensity"
                    font.pointSize: root.fontPt - 1
                    ToolTip.text: "limb darkening multiplied in; unticked shows the temperature"
                    ToolTip.visible: hovered
                    onToggled: {
                        root.statusChanged(Julia.shell_set_surface_field(
                            checked ? "1" : "0", dataIntensityModel.currentText, "0"))
                        root.redraw()
                    }
                }
                ComboBox {
                    id: dataIntensityModel
                    visible: dataIntensity.checked
                    Layout.preferredWidth: dp(110)
                    model: ["linear", "planck"]
                    font.pointSize: root.fontPt - 1
                    onActivated: {
                        root.statusChanged(Julia.shell_set_surface_field(
                            dataIntensity.checked ? "1" : "0", currentText, "0"))
                        root.redraw()
                    }
                }
                Item { Layout.fillWidth: true }
            }

            // Stepping through the epochs. ARROWS, not a slider: the action is "the next
            // night", and a slider makes that a drag at a target a few pixels wide. The table
            // on the left is for reading numbers off a particular night; this is for sweeping
            // the rotation and watching the map turn.
            RowLayout {
                Layout.fillWidth: true
                spacing: dp(4)
                visible: root.viewIndex !== root.nObs + 1 && epochModel.count > 1
                Label { text: "epoch"; color: "#7f8c98"; font.pointSize: root.fontPt - 1 }
                Button {
                    text: "◀"
                    enabled: epochList.currentIndex > 0
                    font.pointSize: root.fontPt - 1
                    onClicked: epochList.currentIndex = epochList.currentIndex - 1
                }
                Label {
                    text: (epochList.currentIndex + 1) + " / " + epochModel.count
                    color: "#555"
                    font.pointSize: root.fontPt - 1
                    horizontalAlignment: Text.AlignHCenter
                    Layout.preferredWidth: dp(56)
                }
                Button {
                    text: "▶"
                    enabled: epochList.currentIndex < epochModel.count - 1
                    font.pointSize: root.fontPt - 1
                    onClicked: epochList.currentIndex = epochList.currentIndex + 1
                }
                Item { Layout.fillWidth: true }
            }

            StackLayout {
                Layout.fillWidth: true
                Layout.fillHeight: true
                // Every observable view shares ONE canvas — they differ in what is pushed into
                // it, not in which axis they are — so the stack has three entries however many
                // observables the file has.
                currentIndex: root.viewIndex < root.nObs ? 0 : root.viewIndex - root.nObs + 1
                MakieArea { id: obsArea;  scene: obsPlot }
                MakieArea { id: skyArea;  scene: skyPlot }
                MakieArea { id: chi2Area; scene: chi2Plot }
            }
        }
    }

    function redraw() {
        if (root.viewIndex < root.nObs) obsArea.update()
        else if (root.viewIndex === root.nObs) skyArea.update()
        else chi2Area.update()
    }
}
