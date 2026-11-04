// Explore: the loaded data, plotted, and compared with the model, the reconstruction or the
// orbit.
//
// The interface is OITOOLS' Explore (`OITOOLS/src/gui/qml/ExploreTab.qml`): a frame of
// controls beside the plot — Dataset, Filter, Plot, Compare, Export — that collapses to give
// the plot the width, and the same canvas, plot kinds, per-group panels and compare views.
// ROTIR adds the EPOCHS, its own section: a surface map is constrained by how the star looks
// at several rotation phases, so the table says which phase is plotted and how well the
// prediction fits each one. And the Orbit tab's orbit is a third thing to compare with.
//
// The dataset box and Open OIFITS are in the window's context bar, as in OITOOLS.

import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml            // supplies `Julia`
import Makie            // MakieArea

// A Pane, not an Item: Qt Quick Controls inherit their font from the nearest CONTROL
// ancestor, and a plain Item has no font. `background: null` keeps it invisible.
Pane {
    id: root
    padding: 0
    background: null
    font.family: root.fontFamily.length > 0 ? root.fontFamily : ""
    font.pointSize: root.fontPt
    property string fontFamily: ""
    property real uiScale: 1
    property real fontPt: 10
    function dp(px) { return Math.round(px * uiScale) }

    signal statusChanged(string s)
    // "Save PNG" goes through the window's file picker; `which` names the plot and the size
    // is the area's, so the PNG matches what is on screen.
    signal saveRequested(string which, int w, int h)
    // Removing an epoch or a dataset changes what every tab shows, so the window refreshes
    // all of them.
    signal refreshAllRequested()

    // The dataset as the window's Dataset box shows it, and what is in it.
    property string datasetName: ""
    property var summaryRows: []
    // Whether the frame of controls is open; the plot keeps its width when it is not.
    property bool frameOpen: true

    // Whether the epoch table's three columns hold reduced χ² or point counts.
    property bool showingChi2: false
    property string filterKey: ""
    property string kindsKey: ""
    property string groupingNoun: ""
    // "data", "model", "image", "orbit", "residmodel", "residimage" or "residorbit": one radio
    // group, so at most one is on — a residual has the prediction subtracted already.
    property string compareMode: "data"
    property bool haveModelObs: false
    property bool haveImageObs: false
    property bool haveOrbitObs: false
    property string pickText: ""

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

    function refreshSummary() {
        var t = Julia.shell_dataset_summary()
        var out = []
        if (t.length > 0) {
            var lines = t.split("\n")
            for (var i = 0; i < lines.length; ++i) {
                var f = lines[i].split("\t")
                if (f.length === 2) out.push({ label: f[0], value: f[1] })
            }
        }
        root.summaryRows = out
    }

    // The kind actually plotted: `t3phi`/`t3amp` have `_max` twins plotted against the
    // longest leg of each triangle; Julia applies the tick.
    function currentKind() { return kindBox.currentText }
    function comparableKind() {
        return ["v2", "t3amp", "t3phi", "visamp", "visphi"].indexOf(kindBox.currentText) >= 0
    }

    // Why a compare choice is disabled. Every reason found, not only the first.
    function compareWhy(need) {
        var why = []
        if (!comparableKind())
            why.push("'" + kindBox.currentText + "' is not a comparable observable; " +
                     "plot v2, t3amp, t3phi, visamp or visphi")
        if (need === "model" && !haveModelObs) why.push("no model defined, or it does not build")
        if (need === "image" && !haveImageObs) why.push("no image reconstructed, or the model is a binary")
        if (need === "orbit" && !haveOrbitObs) why.push("no dataset, or the Orbit tab's star model is 3-D")
        return why.join(" — ")
    }

    // What the observable canvas should show. Pushed to Julia, which owns the drawing.
    function applyView() {
        root.groupingNoun = Julia.shell_obs_grouping_noun(
            (maxBlBox.enabled && maxBlBox.checked) ? kindBox.currentText + "_max"
                                                   : kindBox.currentText)
        root.statusChanged(Julia.shell_set_obs_view(
            kindBox.currentText, colorBox.currentText, root.compareMode,
            logyBox.enabled && logyBox.checked ? "1" : "0",
            maxBlBox.enabled && maxBlBox.checked ? "1" : "0",
            panelsBox.enabled && panelsBox.checked ? "1" : "0"))
        // The epoch table's χ² follows the compare choice (the orbit's, or the model's).
        root.refreshEpochs()
        root.redraw()
    }

    // Plot kinds, with whether the loaded file can supply each; Julia decides availability.
    ListModel {
        id: kindModel
        ListElement { name: "uv";         available: true }
        ListElement { name: "v2";         available: true }
        ListElement { name: "t3phi";      available: true }
        ListElement { name: "t3amp";      available: true }
        ListElement { name: "visamp";     available: true }
        ListElement { name: "visphi";     available: true }
        ListElement { name: "diffphi";    available: true }
        ListElement { name: "diffvisamp"; available: true }
        ListElement { name: "flux";       available: true }
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
            kindModel.setProperty(r, "available", have[kindModel.get(r).name] === true)
        // A kind the new dataset cannot draw falls back to uv, as in OITOOLS.
        if (!kindModel.get(kindBox.currentIndex).available) kindBox.currentIndex = 0
    }

    function refreshCompare() {
        var parts = Julia.shell_obs_compare_available().split(",")
        for (var i = 0; i < parts.length; ++i) {
            var kv = parts[i].split("=")
            if (kv[0] === "model") root.haveModelObs = (kv[1] === "1")
            if (kv[0] === "image") root.haveImageObs = (kv[1] === "1")
            if (kv[0] === "orbit") root.haveOrbitObs = (kv[1] === "1")
        }
        // A choice that has just lost its source stops claiming to show it.
        var m = root.compareMode
        if (((m === "model" || m === "residmodel") && !root.haveModelObs) ||
            ((m === "image" || m === "residimage") && !root.haveImageObs) ||
            ((m === "orbit" || m === "residorbit") && !root.haveOrbitObs)) {
            root.compareMode = "data"
            root.applyView()
        }
    }

    ListModel { id: epochModel }

    function refreshEpochs() {
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
        // Follow Julia: "Add epoch" moves the session to the epoch it appended.
        var ce = parseInt(Julia.shell_current_epoch())
        if (ce >= 1 && ce <= epochModel.count) epochList.currentIndex = ce - 1
        else if (keep >= 0 && keep < epochModel.count) epochList.currentIndex = keep
    }

    function refresh() {
        root.refreshCompare()
        root.refreshSummary()
        // The Filter panel and the kind list follow the DATASET: refreshed only when it
        // changes, since this runs on every parameter edit.
        var dkey = root.datasetName
        if (dkey !== root.filterKey) {
            root.filterKey = dkey
            if (dkey.length > 0) filterPanel.refresh()
        }
        if (dkey !== root.kindsKey) {
            root.kindsKey = dkey
            if (dkey.length > 0) { root.refreshKinds(); root.applyView() }
        }
        root.refreshEpochs()
    }

    function redraw() { obsArea.update() }

    SplitView {
        anchors.fill: parent
        orientation: Qt.Horizontal

        // ── the frame ────────────────────────────────────────────────────────────
        ScrollView {
            id: frame
            visible: root.frameOpen
            SplitView.preferredWidth: dp(460)
            SplitView.minimumWidth: dp(260)
            contentWidth: availableWidth
            clip: true

            ColumnLayout {
                width: frame.availableWidth
                spacing: dp(8)

                GroupBox {
                    title: "Dataset"
                    Layout.fillWidth: true
                    Layout.leftMargin: dp(6); Layout.rightMargin: dp(6); Layout.topMargin: dp(6)
                    GridLayout {
                        anchors.fill: parent
                        columns: 2
                        columnSpacing: dp(8)
                        rowSpacing: dp(2)
                        Label {
                            Layout.columnSpan: 2
                            Layout.fillWidth: true
                            text: root.datasetName.length > 0 ? root.datasetName : "no dataset loaded"
                            color: root.datasetName.length > 0 ? "#000" : "#888"
                            font.bold: root.datasetName.length > 0
                            elide: Text.ElideMiddle
                        }
                        Repeater {
                            model: root.summaryRows
                            delegate: Label {
                                required property var modelData
                                required property int index
                                Layout.row: 1 + index
                                Layout.column: 0
                                Layout.alignment: Qt.AlignTop
                                text: modelData.label
                                color: "#666"
                                font.pointSize: root.fontPt - 1
                            }
                        }
                        Repeater {
                            model: root.summaryRows
                            delegate: Label {
                                required property var modelData
                                required property int index
                                Layout.row: 1 + index
                                Layout.column: 1
                                Layout.fillWidth: true
                                Layout.preferredWidth: 1
                                text: modelData.value
                                wrapMode: Text.Wrap
                                font.pointSize: root.fontPt - 1
                            }
                        }
                    }
                }

                // ROTIR's own section: the epochs, each with its fit.
                GroupBox {
                    title: "Epochs"
                    visible: root.datasetName.length > 0
                    Layout.fillWidth: true
                    Layout.leftMargin: dp(6); Layout.rightMargin: dp(6)
                    ColumnLayout {
                        anchors.fill: parent
                        spacing: dp(2)
                        RowLayout {
                            Layout.fillWidth: true
                            Layout.rightMargin: epochScroll.visible ? epochScroll.width : 0
                            spacing: dp(6)
                            Label { text: "#";     Layout.preferredWidth: dp(24); font.pointSize: root.fontPt - 1; color: "#7f8c98" }
                            Label { text: "MJD";   Layout.fillWidth: true;        font.pointSize: root.fontPt - 1; color: "#7f8c98" }
                            Label { text: "t (d)"; Layout.preferredWidth: dp(52); font.pointSize: root.fontPt - 1; color: "#7f8c98"; horizontalAlignment: Text.AlignRight }
                            Label { text: "V²";    Layout.preferredWidth: dp(50); font.pointSize: root.fontPt - 1; color: "#7f8c98"; horizontalAlignment: Text.AlignRight }
                            Label { text: "T3amp"; Layout.preferredWidth: dp(50); font.pointSize: root.fontPt - 1; color: "#7f8c98"; horizontalAlignment: Text.AlignRight }
                            Label { text: "T3phi"; Layout.preferredWidth: dp(50); font.pointSize: root.fontPt - 1; color: "#7f8c98"; horizontalAlignment: Text.AlignRight }
                        }
                        ListView {
                            id: epochList
                            Layout.fillWidth: true
                            Layout.preferredHeight: Math.min(contentHeight, dp(200))
                            clip: true
                            model: epochModel
                            ScrollBar.vertical: ScrollBar { id: epochScroll }
                            // Up and Down step through the epochs once the table has the focus.
                            focus: true
                            keyNavigationEnabled: true
                            Keys.onUpPressed: if (currentIndex > 0) currentIndex = currentIndex - 1
                            Keys.onDownPressed: if (currentIndex < count - 1) currentIndex = currentIndex + 1
                            onCurrentIndexChanged: {
                                // The repaint is not optional: Makie draws on demand.
                                if (currentIndex >= 0) {
                                    root.statusChanged(Julia.shell_select_epoch(currentIndex + 1))
                                    root.redraw()
                                }
                            }
                            delegate: ItemDelegate {
                                width: epochList.width
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
                                    // One column per observable: a prediction can fit V² well
                                    // and closure phase badly.
                                    Repeater {
                                        model: [v2, t3amp, t3phi]
                                        Label {
                                            text: modelData
                                            Layout.preferredWidth: dp(50)
                                            font.pointSize: root.fontPt
                                            horizontalAlignment: Text.AlignRight
                                            color: !haveChi2 ? "#5a6772"
                                                 : modelData === "—" ? "#7f8c98"
                                                 : parseFloat(modelData) > 3 ? "#c0392b" : "#1e824c"
                                        }
                                    }
                                }
                            }
                        }
                        Label {
                            Layout.fillWidth: true
                            wrapMode: Text.Wrap
                            font.pointSize: root.fontPt - 2
                            color: "#7f8c98"
                            text: !root.showingChi2
                                  ? "V² / T3amp / T3phi: point counts (nothing to compare with)"
                                  : (root.compareMode === "orbit" || root.compareMode === "residorbit")
                                    ? "V² / T3amp / T3phi: reduced χ² of the Orbit tab's orbit"
                                    : "V² / T3amp / T3phi: reduced χ² of the Model tab's model"
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
                                enabled: root.datasetName.length > 0
                                font.pointSize: root.fontPt - 1
                                onClicked: {
                                    root.statusChanged(Julia.shell_close_dataset())
                                    root.refreshAllRequested()
                                }
                            }
                            Item { Layout.fillWidth: true }
                        }
                    }
                }

                // OITOOLS' filtering, epoch by epoch: Apply makes a new dataset and selects
                // it, and the unfiltered one stays in the Dataset box as the undo.
                GroupBox {
                    title: "Filter"
                    visible: root.datasetName.length > 0
                    Layout.fillWidth: true
                    Layout.leftMargin: dp(6); Layout.rightMargin: dp(6)
                    FilterPanel {
                        id: filterPanel
                        anchors.fill: parent
                        uiScale: root.uiScale
                        baseFontPt: root.fontPt
                        onStatusReported: function (t) { root.statusChanged(t) }
                        onApplied: root.refreshAllRequested()
                        onPreviewChanged: root.redraw()
                        onSelectionCleared: root.redraw()
                    }
                }

                GroupBox {
                    title: "Plot"
                    Layout.fillWidth: true
                    Layout.leftMargin: dp(6); Layout.rightMargin: dp(6)
                    GridLayout {
                        anchors.fill: parent
                        columns: 2
                        columnSpacing: dp(8)
                        rowSpacing: dp(4)

                        Label { text: "Plot" }
                        ComboBox {
                            id: kindBox
                            objectName: "kindBox"
                            Layout.fillWidth: true
                            textRole: "name"
                            model: kindModel
                            currentIndex: 1
                            delegate: ItemDelegate {
                                id: kindItem
                                required property var model
                                required property int index
                                width: kindBox.width
                                enabled: model.available
                                highlighted: kindBox.highlightedIndex === index
                                // Dimmed by hand: this style draws a disabled ItemDelegate
                                // like an enabled one.
                                contentItem: Label {
                                    text: kindItem.model.name
                                    opacity: kindItem.enabled ? 1.0 : 0.4
                                }
                                ToolTip.visible: hovered && !model.available
                                ToolTip.text: "this dataset has no " + model.name + " data"
                            }
                            onActivated: root.applyView()
                        }

                        Label { text: "Colour by" }
                        ComboBox {
                            id: colorBox
                            Layout.fillWidth: true
                            model: ["baseline", "wav", "mjd", "none"]
                            onActivated: root.applyView()
                        }

                        CheckBox {
                            id: logyBox
                            Layout.columnSpan: 2
                            text: "log y"
                            // Only the positive quantities: a log axis on a phase drops half
                            // the points rather than rescaling them.
                            enabled: ["v2", "visamp", "t3amp", "flux", "diffvisamp"]
                                     .indexOf(kindBox.currentText) >= 0
                            onToggled: root.applyView()
                        }
                        CheckBox {
                            id: panelsBox
                            Layout.columnSpan: 2
                            // Labelled from OITOOLS: closure quantities group by triplet and
                            // flux by station.
                            text: "per " + (root.groupingNoun.length > 0 ? root.groupingNoun : "group")
                            enabled: root.groupingNoun.length > 0
                            onToggled: root.applyView()
                        }
                        CheckBox {
                            id: maxBlBox
                            Layout.columnSpan: 2
                            text: "max baseline"
                            enabled: ["t3phi", "t3amp"].indexOf(kindBox.currentText) >= 0
                            ToolTip.visible: hovered
                            ToolTip.text: "closure quantities against their longest baseline, " +
                                          "not the geometric mean of the three"
                            onToggled: root.applyView()
                        }
                    }
                }

                // The model's, the reconstruction's or the orbit's observables over the data,
                // or the residuals in their place. One choice at a time; each follows the plot
                // kind, on the same axes and in the same groups.
                GroupBox {
                    title: "Compare"
                    Layout.fillWidth: true
                    Layout.leftMargin: dp(6); Layout.rightMargin: dp(6)
                    ButtonGroup { id: compareGroup }
                    ColumnLayout {
                        anchors.fill: parent
                        spacing: 0
                        Repeater {
                            model: [
                                { key: "data",       label: "data only",          need: ""      },
                                { key: "model",      label: "model observables",  need: "model" },
                                { key: "image",      label: "image observables",  need: "image" },
                                { key: "orbit",      label: "orbit observables",  need: "orbit" },
                                { key: "residmodel", label: "residuals (model)",  need: "model" },
                                { key: "residimage", label: "residuals (image)",  need: "image" },
                                { key: "residorbit", label: "residuals (orbit)",  need: "orbit" }
                            ]
                            RadioButton {
                                required property var modelData
                                objectName: "compare:" + modelData.key
                                text: modelData.label
                                ButtonGroup.group: compareGroup
                                checked: root.compareMode === modelData.key
                                enabled: modelData.key === "data" ||
                                         (root.comparableKind() &&
                                          (modelData.need === "model" ? root.haveModelObs :
                                           modelData.need === "image" ? root.haveImageObs :
                                                                        root.haveOrbitObs))
                                ToolTip.visible: hovered && !enabled
                                ToolTip.text: root.compareWhy(modelData.need)
                                onClicked: { root.compareMode = modelData.key; root.applyView() }
                            }
                        }
                    }
                }

                GroupBox {
                    title: "Export"
                    Layout.fillWidth: true
                    Layout.leftMargin: dp(6); Layout.rightMargin: dp(6)
                    RowLayout {
                        anchors.fill: parent
                        Button {
                            text: "OIFITS…"
                            enabled: root.datasetName.length > 0
                            ToolTip.visible: hovered
                            ToolTip.text: "write this dataset, every epoch, as filtered, to a new OIFITS file"
                            onClicked: root.saveRequested("oifits", 0, 0)
                        }
                    }
                }

                Item { Layout.fillHeight: true }
            }
        }

        // ── the plot ─────────────────────────────────────────────────────────────
        ColumnLayout {
            SplitView.fillWidth: true
            spacing: 0

            Item {
                Layout.fillWidth: true
                Layout.fillHeight: true

                MakieArea {
                    id: obsArea
                    anchors.fill: parent
                    scene: obsPlot
                }

                // Over the plot, as OITOOLS places them: the frame's toggle in one corner,
                // Save PNG in the other.
                Button {
                    anchors.left: parent.left
                    anchors.top: parent.top
                    anchors.margins: dp(6)
                    text: root.frameOpen ? "◀" : "▶"
                    opacity: hovered ? 1.0 : 0.5
                    ToolTip.visible: hovered
                    ToolTip.text: root.frameOpen ? "hide the controls" : "show the controls"
                    onClicked: root.frameOpen = !root.frameOpen
                }
                Button {
                    anchors.right: parent.right
                    anchors.top: parent.top
                    anchors.margins: dp(6)
                    text: "Save PNG"
                    enabled: root.datasetName.length > 0
                    opacity: hovered ? 1.0 : 0.5
                    ToolTip.visible: hovered
                    ToolTip.text: "write this plot to a PNG file"
                    onClicked: root.saveRequested("obs", obsArea.width, obsArea.height)
                }
            }

            // Under the plot: ROTIR's epoch stepper, then the identified point as OITOOLS
            // shows it. Not over the plot — the canvas draws its title there.
            RowLayout {
                Layout.fillWidth: true
                Layout.margins: dp(6)
                spacing: dp(4)
                Button {
                    visible: epochModel.count > 1
                    text: "◀"
                    enabled: epochList.currentIndex > 0
                    ToolTip.visible: hovered; ToolTip.text: "previous epoch"
                    onClicked: epochList.currentIndex = epochList.currentIndex - 1
                }
                Label {
                    visible: epochModel.count > 1
                    text: "epoch " + (epochList.currentIndex + 1) + " / " + epochModel.count
                    color: "#555"
                }
                Button {
                    visible: epochModel.count > 1
                    text: "▶"
                    enabled: epochList.currentIndex < epochModel.count - 1
                    ToolTip.visible: hovered; ToolTip.text: "next epoch"
                    onClicked: epochList.currentIndex = epochList.currentIndex + 1
                }
            Label {
                Layout.fillWidth: true
                Layout.leftMargin: dp(8)
                text: root.pickText.length > 0 ? root.pickText
                    : filterPanel.selectTool !== "none" ? "drag on the plot to select points"
                    : "click a point to identify it"
                color: root.pickText.length > 0 ? "#000" : "#888"
                font.family: "monospace"
                elide: Text.ElideRight
            }
            }
        }
    }
}
