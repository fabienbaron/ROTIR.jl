// Explore: what was observed, and how well the current model reproduces it.
//
// The observable plots, their options and the filtering are OITOOLS' own — the same canvas,
// the same plot kinds, per-baseline panels and residuals — so a dataset reads the same in both
// GUIs. What ROTIR adds is the EPOCH structure: a surface map is constrained by how the star
// looks at several rotation phases, so the questions this tab answers are "which phase is
// this", "what does the model look like there", and "which epoch does the model fit badly".

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
    // "Save plot…" goes through the window's file picker rather than writing a fixed name:
    // the window owns the picker, and a tab has no way to reach it. `which` names the plot,
    // and the size is the area's, so the PNG matches what is on screen.
    signal saveRequested(string which, int w, int h)
    // Removing an epoch or a dataset changes what every tab shows, so the window refreshes
    // all of them rather than this tab redrawing itself in isolation.
    signal refreshAllRequested()

    // Whether the three observable columns hold reduced χ² or point counts. Julia decides —
    // there is a χ² only when a model is set — and the header says which, rather than putting
    // "χ²ᵣ" above a column of counts.
    property bool showingChi2: false
    property string filterKey: ""
    property string kindsKey: ""
    // What one panel of the per-group view holds — baseline, triplet or station — from OITOOLS.
    property string groupingNoun: "baseline"
    // "data", "model", "image", "residmodel" or "residimage": one radio group, so at most one
    // is ever on — a residual has the model subtracted already, and overplotting it there
    // would draw the prediction against its own difference.
    property string compareMode: "model"
    property bool haveModelObs: false
    property bool haveImageObs: false
    property string pickText: ""

    // The point a click identified and the selection a drag made are Julia's, from Makie's
    // event stream, so they are polled — as OITOOLS' Explore does; 150 ms is imperceptible.
    Timer {
        interval: 150; running: root.visible; repeat: true
        onTriggered: {
            var t = Julia.shell_pick_text()
            if (t !== root.pickText) root.pickText = t
            if (filterPanel.selectTool !== "none" || filterPanel.selectionText.length > 0) {
                var st = Julia.shell_selection_text()
                if (st !== filterPanel.selectionText) filterPanel.selectionText = st
            }
        }
    }

    function currentKind() {
        return kindBox.currentIndex >= 0 ? kindModel.get(kindBox.currentIndex).key : "uv"
    }
    function predictedKind() { return ["v2", "t3amp", "t3phi"].indexOf(currentKind()) >= 0 }

    // What the observable canvas should show. Pushed to Julia, which owns the drawing.
    function applyView() {
        root.groupingNoun = Julia.shell_obs_grouping_noun(root.currentKind())
        root.statusChanged(Julia.shell_set_obs_view(
            root.currentKind(), colorBox.currentText, root.compareMode,
            logBox.enabled && logBox.checked ? "1" : "0",
            maxBaseBox.enabled && maxBaseBox.checked ? "1" : "0",
            panelsBox.enabled && panelsBox.checked ? "1" : "0"))
        root.redraw()
    }

    // Every OITOOLS plot kind, with whether this dataset can supply it. Shown but disabled
    // when it cannot, as OITOOLS does: which observables a file carries is worth seeing.
    ListModel {
        id: kindModel
        ListElement { key: "uv";         label: "uv coverage";   available: true }
        ListElement { key: "v2";         label: "V²";            available: true }
        ListElement { key: "t3phi";      label: "T3 phase";      available: true }
        ListElement { key: "t3amp";      label: "T3 amp";        available: true }
        ListElement { key: "visamp";     label: "vis amp";       available: true }
        ListElement { key: "visphi";     label: "vis phase";     available: true }
        ListElement { key: "diffphi";    label: "diff. phase";   available: true }
        ListElement { key: "diffvisamp"; label: "diff. vis amp"; available: true }
        ListElement { key: "flux";       label: "flux";          available: true }
    }

    function refreshKinds() {
        var flags = Julia.shell_obs_kinds()
        if (flags.length === 0) return
        var have = {}
        var fl = flags.split("\n")
        for (var j = 0; j < fl.length; ++j) {
            var g = fl[j].split("\t")
            if (g.length >= 2) have[g[0]] = (g[1] === "1")
        }
        for (var r = 0; r < kindModel.count; ++r)
            kindModel.setProperty(r, "available", have[kindModel.get(r).key] === true)
        // A kind the new dataset cannot draw falls back to V², or to uv coverage.
        if (!kindModel.get(kindBox.currentIndex).available)
            kindBox.currentIndex = kindModel.get(1).available ? 1 : 0
    }

    function refreshCompare() {
        var parts = Julia.shell_obs_compare_available().split(",")
        for (var i = 0; i < parts.length; ++i) {
            var kv = parts[i].split("=")
            if (kv[0] === "model") root.haveModelObs = (kv[1] === "1")
            if (kv[0] === "image") root.haveImageObs = (kv[1] === "1")
        }
    }

    ListModel { id: epochModel }
    ListModel { id: datasetModel }

    function refresh() {
        root.refreshCompare()
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
        // The Filter panel offers what THIS dataset holds, and shows a filtered dataset's own
        // settings again — but only when the dataset changes: this runs on every parameter
        // edit, and the choices come from reading the files' headers.
        var dkey = dsBox.currentIndex >= 0 && dsBox.currentIndex < datasetModel.count
                   ? datasetModel.count + ":" + datasetModel.get(dsBox.currentIndex).dname : ""
        if (dkey !== root.filterKey) {
            root.filterKey = dkey
            if (dkey.length > 0) filterPanel.refresh()
        }
        // Which plot kinds this dataset can draw, again only when the dataset changes.
        if (dkey !== root.kindsKey) {
            root.kindsKey = dkey
            if (dkey.length > 0) { root.refreshKinds(); root.applyView() }
        }

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
            //
            // AS TALL AS ITS ROWS, up to the room there is: one epoch gets one row, not half
            // the column. `fillHeight` with the frame's own implicit height as the maximum
            // grows it with the epoch count and lets the Filter panel below share the rest;
            // past that the list scrolls.
            Frame {
                id: epochFrame
                Layout.fillWidth: true
                Layout.fillHeight: true
                Layout.maximumHeight: implicitHeight
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
                    // Every row, or at least three of them once that no longer fits.
                    Layout.preferredHeight: contentHeight
                    Layout.minimumHeight: Math.min(contentHeight, dp(90))
                    clip: true
                    model: epochModel
                    ScrollBar.vertical: ScrollBar { id: epochScroll }
                    // UP AND DOWN STEP THROUGH THE EPOCHS once the table has been clicked:
                    // a click gives the list the focus (the rows take none of their own), and
                    // moving the current row is all a step needs — `onCurrentIndexChanged`
                    // below selects the epoch and redraws.
                    focus: true
                    keyNavigationEnabled: true
                    Keys.onUpPressed: if (currentIndex > 0) currentIndex = currentIndex - 1
                    Keys.onDownPressed: if (currentIndex < count - 1) currentIndex = currentIndex + 1
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
                        focusPolicy: Qt.NoFocus
                        onClicked: {
                            epochList.currentIndex = index
                            epochList.forceActiveFocus()
                        }
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

            // OITOOLS' filtering, epoch by epoch: Apply makes a new dataset and selects it,
            // and the unfiltered one stays in the Dataset box as the undo.
            Label {
                text: "Filter"; font.bold: true; font.pointSize: root.fontPt
                visible: datasetModel.count > 0
            }
            ScrollView {
                id: filterScroll
                visible: datasetModel.count > 0
                Layout.fillWidth: true
                // Takes what the epoch table leaves, up to its own height; scrolls below that.
                Layout.fillHeight: true
                Layout.preferredHeight: filterPanel.implicitHeight
                Layout.maximumHeight: filterPanel.implicitHeight
                Layout.minimumHeight: Math.min(filterPanel.implicitHeight, dp(120))
                clip: true
                contentWidth: availableWidth
                FilterPanel {
                    id: filterPanel
                    width: filterScroll.availableWidth
                    uiScale: root.uiScale
                    baseFontPt: root.fontPt
                    onStatusReported: function (t) { root.statusChanged(t) }
                    onApplied: root.refreshAllRequested()
                    onPreviewChanged: root.redraw()
                    onSelectionCleared: root.redraw()
                }
            }

            // The dataset on screen, every epoch, filtered as it is — OITOOLS' Export.
            RowLayout {
                visible: datasetModel.count > 0
                Layout.fillWidth: true
                spacing: dp(6)
                Button {
                    text: "Export OIFITS…"
                    font.pointSize: root.fontPt - 1
                    ToolTip.text: "write this dataset, as filtered, to a new OIFITS file"
                    ToolTip.visible: hovered
                    onClicked: root.saveRequested("oifits", 0, 0)
                }
                Item { Layout.fillWidth: true }
            }

            // Whatever both leave stays BELOW them, rather than being spread between the rows.
            Item { Layout.fillHeight: true; Layout.preferredHeight: 0 }
        }

        // ── right: the observable plot ───────────────────────────────────────
        //
        // V² against baseline, closure phase, uv coverage — the OITOOLS GUI's own plots, drawn
        // by its own canvas (see `build_obs_canvas`), for the selected epoch. The model's or a
        // reconstruction's prediction goes over them, or their residuals replace them. The
        // model's own pictures are on the Model tab.
        ColumnLayout {
            Layout.fillWidth: true
            Layout.fillHeight: true
            spacing: dp(4)

            RowLayout {
                Layout.fillWidth: true
                spacing: dp(6)
                Button {
                    text: "Save plot…"
                    font.pointSize: root.fontPt - 1
                    // Rebuilt offscreen, not grabbed: see src/gui/snapshot.jl.
                    onClicked: root.saveRequested("obs", obsArea.width, obsArea.height)
                }
                Item { Layout.fillWidth: true }
                Label {
                    text: "the selected epoch · scroll to zoom · right-click resets"
                    color: "#7f8c98"
                    font.pointSize: root.fontPt - 1
                }
            }

            // OITOOLS' Plot options, in the order its Explore lists them.
            RowLayout {
                Layout.fillWidth: true
                spacing: dp(6)
                Label { text: "plot"; color: "#7f8c98"; font.pointSize: root.fontPt - 1 }
                ComboBox {
                    id: kindBox
                    Layout.preferredWidth: dp(130)
                    model: kindModel
                    textRole: "label"
                    currentIndex: 1
                    font.pointSize: root.fontPt - 1
                    delegate: ItemDelegate {
                        id: kindItem
                        required property var model
                        required property int index
                        width: kindBox.width
                        enabled: model.available
                        highlighted: kindBox.highlightedIndex === index
                        // Dimmed by hand: Fusion draws a disabled ItemDelegate like an enabled one.
                        contentItem: Label {
                            text: kindItem.model.label
                            opacity: kindItem.enabled ? 1.0 : 0.4
                            font.pointSize: root.fontPt - 1
                        }
                        ToolTip.visible: hovered && !model.available
                        ToolTip.text: "this dataset has no " + model.label + " data"
                    }
                    onActivated: root.applyView()
                }
                Label { text: "colour by"; color: "#7f8c98"; font.pointSize: root.fontPt - 1 }
                ComboBox {
                    id: colorBox
                    Layout.preferredWidth: dp(110)
                    model: ["baseline", "wav", "mjd", "none"]
                    font.pointSize: root.fontPt - 1
                    onActivated: root.applyView()
                }
                CheckBox {
                    id: logBox
                    text: "log y"
                    // Only where it means something: the positive quantities. A phase takes
                    // both signs, uv coverage is a geometry, and a residual is centred on zero.
                    enabled: ["v2", "visamp", "t3amp", "flux", "diffvisamp"]
                             .indexOf(root.currentKind()) >= 0 &&
                             root.compareMode.indexOf("resid") !== 0
                    font.pointSize: root.fontPt - 1
                    ToolTip.text: "points at or below zero are dropped, not piled on the floor"
                    ToolTip.visible: hovered
                    onToggled: root.applyView()
                }
                CheckBox {
                    id: panelsBox
                    // Labelled from OITOOLS: closure quantities group by triplet and flux by
                    // station, and "per baseline" over a grid of triangles would be wrong.
                    text: "per " + root.groupingNoun
                    enabled: root.currentKind() !== "uv"
                    font.pointSize: root.fontPt - 1
                    ToolTip.text: "one panel per " + root.groupingNoun + ", against wavelength"
                    ToolTip.visible: hovered
                    onToggled: root.applyView()
                }
                CheckBox {
                    id: maxBaseBox
                    text: "max baseline"
                    // Only the closure quantities have a choice to make: a triangle has three
                    // legs, so it has both a mean baseline and a longest one.
                    enabled: ["t3amp", "t3phi"].indexOf(root.currentKind()) >= 0 &&
                             !(panelsBox.enabled && panelsBox.checked)
                    font.pointSize: root.fontPt - 1
                    ToolTip.text: "plot against the LONGEST leg of each triangle instead of " +
                                  "the geometric mean — the resolution the closure actually " +
                                  "probes"
                    ToolTip.visible: hovered
                    onToggled: root.applyView()
                }
                Item { Layout.fillWidth: true }
            }

            // OITOOLS' Compare choices. ROTIR predicts V², T3amp and T3φ, so they apply there;
            // the image is the current reconstruction, on the current model's geometry.
            RowLayout {
                Layout.fillWidth: true
                spacing: dp(6)
                Label { text: "compare"; color: "#7f8c98"; font.pointSize: root.fontPt - 1 }
                ButtonGroup { id: compareGroup }
                Repeater {
                    model: [
                        { key: "data",       label: "data only",         needModel: false },
                        { key: "model",      label: "model",             needModel: true  },
                        { key: "image",      label: "image",             needModel: false },
                        { key: "residmodel", label: "residuals (model)", needModel: true  },
                        { key: "residimage", label: "residuals (image)", needModel: false }
                    ]
                    RadioButton {
                        required property var modelData
                        text: modelData.label
                        font.pointSize: root.fontPt - 1
                        ButtonGroup.group: compareGroup
                        checked: root.compareMode === modelData.key
                        enabled: modelData.key === "data" ||
                                 (root.predictedKind() && !(panelsBox.enabled && panelsBox.checked) &&
                                  (modelData.needModel ? root.haveModelObs : root.haveImageObs))
                        ToolTip.visible: hovered && !enabled
                        ToolTip.text: !root.predictedKind()
                            ? "ROTIR predicts V², T3 amp and T3 phase only"
                            : (panelsBox.enabled && panelsBox.checked)
                              ? "the per-" + root.groupingNoun + " panels show the data alone"
                              : modelData.needModel ? "no model, or its parameters do not build"
                                                    : "no reconstruction yet, or the model is a binary"
                        onClicked: { root.compareMode = modelData.key; root.applyView() }
                    }
                }
                Item { Layout.fillWidth: true }
            }

            // Stepping through the epochs. ARROWS, not a slider: the action is "the next
            // night", and a slider makes that a drag at a target a few pixels wide. The Up and
            // Down keys do the same in the table on the left.
            RowLayout {
                Layout.fillWidth: true
                spacing: dp(4)
                visible: epochModel.count > 1
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

            // Every plot kind shares ONE canvas: they differ in what is pushed into it.
            MakieArea {
                id: obsArea
                Layout.fillWidth: true
                Layout.fillHeight: true
                scene: obsPlot
            }

            // The point a click identified, as OITOOLS shows it under its plot.
            Label {
                Layout.fillWidth: true
                text: root.pickText.length > 0 ? root.pickText
                    : panelsBox.enabled && panelsBox.checked ? ""
                    : filterPanel.selectTool !== "none" ? "drag on the plot to select points"
                    : "click a point to identify it"
                color: root.pickText.length > 0 ? "#000" : "#888"
                font.family: "monospace"
                font.pointSize: root.fontPt - 1
                elide: Text.ElideRight
            }
        }
    }

    function redraw() {
        obsArea.update()
    }
}
