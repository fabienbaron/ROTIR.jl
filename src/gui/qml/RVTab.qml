// RV: a table of radial velocities, drawn against the orbit of the Orbit tab.
//
// The velocity parameters — parallax, mass ratio, gamma, one offset per extra instrument, and the two
// jitters — belong to the ORBIT (they are stored with its elements), so editing them here
// edits the Orbit tab's values and the other way round. Fitting happens in the Orbit tab.
//
// The parameter rows use the Orbit tab's columns and callbacks (`shell_set_orbit_*`), so a
// velocity parameter is fixed, free or tied exactly as an element is.

import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import Makie

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
    signal saveRequested(string which, int w, int h)
    signal pickFile(string mode)          // the window owns the picker

    ListModel { id: paramModel }
    ListModel { id: jitterModel }
    ListModel { id: rowModel }

    property bool hasTable: false

    readonly property int wName:  dp(112)
    readonly property int wValue: dp(96)
    readonly property int wState: dp(76)
    readonly property int wBound: dp(66)

    function fillParams(model, rows) {
        model.clear()
        if (rows.length === 0) return
        var ls = rows.split("\n")
        for (var i = 0; i < ls.length; ++i) {
            var f = ls[i].split("\t")
            if (f.length < 8) continue
            model.append({ pname: f[0], punit: f[1], pvalue: f[2], pstate: f[3],
                           plo: f[4], phi: f[5], ptie: f[6], pdoc: f[7] })
        }
    }

    function refresh() {
        var s = Julia.shell_rv_summary()
        root.hasTable = s.length > 0
        summary.text = root.hasTable ? s : "no RV table — CSV with columns MJD,RV,sigma,star"
        root.fillParams(paramModel, Julia.shell_rv_params())
        root.fillParams(jitterModel, Julia.shell_rv_jitter())
        stats.text = Julia.shell_rv_stats()
        rowModel.clear()
        var rows = Julia.shell_rv_rows()
        if (rows.length > 0) {
            var ls = rows.split("\n")
            for (var i = 0; i < ls.length; ++i) {
                var f = ls[i].split("\t")
                if (f.length < 8) continue
                rowModel.append({ rmjd: f[0], rrv: f[1], rsig: f[2], rstar: f[3],
                                  rinst: f[4], rmodel: f[5], roc: f[6], rnorm: f[7] })
            }
        }
        var v = Julia.shell_rv_view().split("\t")
        if (v.length >= 2) {
            xBox.currentIndex = v[0] === "time" ? 1 : 0
            normBox.checked = v[1] === "1"
        }
    }

    function redraw() { rvArea.update() }

    function pushView() {
        root.statusChanged(Julia.shell_rv_set_view(xBox.currentText, normBox.checked ? "1" : "0"))
        root.redraw()
    }

    // The Orbit tab's parameter row, with the same columns and callbacks.
    Component {
        id: rvParamRow
        RowLayout {
            width: ListView.view.width
            spacing: dp(4)
            Label {
                Layout.preferredWidth: root.wName
                Layout.alignment: Qt.AlignVCenter
                text: pname + (punit.length > 0 ? " (" + punit + ")" : "")
                elide: Text.ElideRight
                font.pointSize: root.fontPt
                ToolTip.text: pdoc
                ToolTip.visible: rma.containsMouse && pdoc.length > 0
                MouseArea { id: rma; anchors.fill: parent; hoverEnabled: true }
            }
            TextField {
                Layout.preferredWidth: root.wValue
                text: pvalue
                font.pointSize: root.fontPt
                selectByMouse: true
                enabled: pstate !== "tied"
                onEditingFinished: {
                    var m = Julia.shell_set_orbit_param(pname, text)
                    if (m.length > 0) root.statusChanged(m)
                    root.refresh(); root.redraw()
                }
            }
            ComboBox {
                Layout.preferredWidth: root.wState
                model: ["fixed", "free", "tied"]
                font.pointSize: root.fontPt - 1
                currentIndex: pstate === "free" ? 1 : pstate === "tied" ? 2 : 0
                onActivated: {
                    root.statusChanged(Julia.shell_set_orbit_state(pname, model[currentIndex]))
                    root.refresh()
                }
            }
            Item {
                Layout.fillWidth: true
                Layout.preferredHeight: rLo.implicitHeight
                TextField {
                    id: rLo
                    visible: pstate === "free"
                    anchors.left: parent.left
                    width: root.wBound
                    text: plo
                    font.pointSize: root.fontPt - 2
                    selectByMouse: true
                    ToolTip.text: "lower bound"; ToolTip.visible: hovered
                    onEditingFinished: root.statusChanged(
                        Julia.shell_set_orbit_bound(pname, text, rHi.text))
                }
                TextField {
                    id: rHi
                    visible: pstate === "free"
                    anchors.left: parent.left
                    anchors.leftMargin: root.wBound + dp(4)
                    width: root.wBound
                    text: phi
                    font.pointSize: root.fontPt - 2
                    selectByMouse: true
                    ToolTip.text: "upper bound"; ToolTip.visible: hovered
                    onEditingFinished: root.statusChanged(
                        Julia.shell_set_orbit_bound(pname, rLo.text, text))
                }
                TextField {
                    visible: pstate === "tied"
                    anchors.left: parent.left
                    anchors.right: parent.right
                    text: ptie
                    placeholderText: "expression, e.g. 1000 / 77"
                    font.pointSize: root.fontPt - 1
                    selectByMouse: true
                    onEditingFinished: {
                        root.statusChanged(Julia.shell_set_orbit_tie(pname, text))
                        root.refresh(); root.redraw()
                    }
                }
            }
        }
    }

    contentItem: RowLayout {
        spacing: dp(8)

        // ── left: the table and the velocity parameters ──────────────────────
        ColumnLayout {
            Layout.fillWidth: false
            Layout.preferredWidth: dp(540)
            Layout.minimumWidth: dp(540)
            Layout.maximumWidth: dp(540)
            Layout.fillHeight: true
            spacing: dp(4)

            RowLayout {
                Layout.fillWidth: true
                spacing: dp(6)
                Button {
                    objectName: "rvOpen"
                    text: "Open RV table…"
                    font.pointSize: root.fontPt - 1
                    ToolTip.text: "CSV with columns MJD (or JD), RV, sigma, star (a/b), and optionally instrument"
                    ToolTip.visible: hovered
                    onClicked: root.pickFile("rv")
                }
                Button {
                    text: "Close"
                    enabled: root.hasTable
                    font.pointSize: root.fontPt - 1
                    onClicked: {
                        root.statusChanged(Julia.shell_rv_close())
                        root.refresh(); root.redraw()
                    }
                }
                Item { Layout.fillWidth: true }
            }
            Label {
                id: summary
                Layout.fillWidth: true
                wrapMode: Text.Wrap
                color: root.hasTable ? palette.text : "#7f8c98"
                font.pointSize: root.fontPt - 1
            }

            Label { text: "Velocity parameters"; font.bold: true; font.pointSize: root.fontPt
                    visible: root.hasTable }
            Label {
                visible: root.hasTable
                Layout.fillWidth: true
                wrapMode: Text.Wrap
                text: "Stored with the orbit, and fitted in the Orbit tab. The amplitudes " +
                      "follow from the parallax and the mass ratio at the orbit's a, i, P and e; " +
                      "they are shown below with the distance and the dynamical masses."
                color: "#7f8c98"
                font.pointSize: root.fontPt - 2
            }
            Frame {
                visible: root.hasTable
                Layout.fillWidth: true
                Layout.preferredHeight: Math.max(1, paramModel.count) * dp(34) + dp(12)
                ListView {
                    anchors.fill: parent
                    clip: true
                    interactive: false
                    model: paramModel
                    spacing: dp(2)
                    delegate: rvParamRow
                }
            }

            Label { text: "Jitter"; font.bold: true; font.pointSize: root.fontPt
                    visible: root.hasTable }
            Label {
                visible: root.hasTable
                Layout.fillWidth: true
                wrapMode: Text.Wrap
                text: "Added in quadrature to each star's errors, √(σ² + s²). A setting: the " +
                      "fit uses these values and does not fit them."
                color: "#7f8c98"
                font.pointSize: root.fontPt - 2
            }
            Frame {
                visible: root.hasTable
                Layout.fillWidth: true
                Layout.preferredHeight: dp(2 * 34 + 12)
                ListView {
                    anchors.fill: parent
                    clip: true
                    interactive: false
                    model: jitterModel
                    spacing: dp(2)
                    // Value only: a jitter is set here, never fixed/free/tied.
                    delegate: RowLayout {
                        width: ListView.view.width
                        spacing: dp(4)
                        Label {
                            Layout.preferredWidth: root.wName
                            text: pname + (punit.length > 0 ? " (" + punit + ")" : "")
                            font.pointSize: root.fontPt
                            ToolTip.text: pdoc
                            ToolTip.visible: jma.containsMouse && pdoc.length > 0
                            MouseArea { id: jma; anchors.fill: parent; hoverEnabled: true }
                        }
                        TextField {
                            Layout.preferredWidth: root.wValue
                            text: pvalue
                            font.pointSize: root.fontPt
                            selectByMouse: true
                            onEditingFinished: {
                                var m = Julia.shell_set_orbit_param(pname, text)
                                if (m.length > 0) root.statusChanged(m)
                                root.refresh(); root.redraw()
                            }
                        }
                        Item { Layout.fillWidth: true }
                    }
                }
            }

            Frame {
                visible: root.hasTable
                Layout.fillWidth: true
                Text {
                    id: stats
                    width: parent.width
                    font.family: "monospace"
                    font.pointSize: root.fontPt - 2
                    color: palette.text
                }
            }

            RowLayout {
                visible: root.hasTable
                Layout.fillWidth: true
                spacing: dp(6)
                Label { text: "x axis"; color: "#7f8c98"; font.pointSize: root.fontPt - 1 }
                ComboBox {
                    id: xBox
                    model: ["phase", "time"]
                    font.pointSize: root.fontPt - 1
                    onActivated: root.pushView()
                }
                CheckBox {
                    id: normBox
                    text: "residuals in σ"
                    font.pointSize: root.fontPt - 1
                    onToggled: root.pushView()
                }
                Item { Layout.fillWidth: true }
            }

            // The table: one row per velocity.
            RowLayout {
                visible: root.hasTable
                Layout.fillWidth: true
                Layout.leftMargin: dp(6)
                spacing: dp(4)
                Repeater {
                    model: [["MJD", 92], ["RV", 62], ["σ", 48], ["", 18], ["instr.", 58],
                            ["model", 62], ["O−C", 58], ["(O−C)/σ'", 64]]
                    Label { Layout.preferredWidth: dp(modelData[1]); text: modelData[0]
                            color: "#7f8c98"; font.pointSize: root.fontPt - 2 }
                }
            }
            Frame {
                visible: root.hasTable
                Layout.fillWidth: true
                Layout.fillHeight: true
                ListView {
                    id: rowList
                    anchors.fill: parent
                    clip: true
                    model: rowModel
                    ScrollBar.vertical: ScrollBar {}
                    delegate: RowLayout {
                        width: ListView.view.width
                        spacing: dp(4)
                        Label { Layout.preferredWidth: dp(92); text: rmjd; font.family: "monospace"
                                font.pointSize: root.fontPt - 2 }
                        Label { Layout.preferredWidth: dp(62); text: rrv; font.family: "monospace"
                                font.pointSize: root.fontPt - 2; horizontalAlignment: Text.AlignRight }
                        Label { Layout.preferredWidth: dp(48); text: rsig; font.family: "monospace"
                                font.pointSize: root.fontPt - 2; horizontalAlignment: Text.AlignRight }
                        Label { Layout.preferredWidth: dp(18); text: rstar
                                color: rstar === "a" ? "#0072b2" : "#e69f00"; font.bold: true
                                font.pointSize: root.fontPt - 2; horizontalAlignment: Text.AlignHCenter }
                        Label { Layout.preferredWidth: dp(58); text: rinst; elide: Text.ElideRight
                                font.pointSize: root.fontPt - 2 }
                        Label { Layout.preferredWidth: dp(62); text: rmodel; font.family: "monospace"
                                font.pointSize: root.fontPt - 2; horizontalAlignment: Text.AlignRight }
                        Label { Layout.preferredWidth: dp(58); text: roc; font.family: "monospace"
                                font.pointSize: root.fontPt - 2; horizontalAlignment: Text.AlignRight }
                        Label { Layout.preferredWidth: dp(64); text: rnorm; font.family: "monospace"
                                font.pointSize: root.fontPt - 2; horizontalAlignment: Text.AlignRight
                                color: Math.abs(parseFloat(rnorm)) > 3 ? "#c0392b" : palette.text }
                    }
                }
            }
            Item { visible: !root.hasTable; Layout.fillHeight: true }
        }

        // ── right: the velocities and the residuals ──────────────────────────
        ColumnLayout {
            Layout.fillWidth: true
            Layout.fillHeight: true
            spacing: dp(4)
            RowLayout {
                Layout.fillWidth: true
                spacing: dp(6)
                Label {
                    text: "a (primary) blue · b (secondary) orange · model from the Orbit tab · " +
                          "marker shape = instrument"
                    color: "#7f8c98"
                    font.pointSize: root.fontPt - 1
                }
                Item { Layout.fillWidth: true }
                Button {
                    text: "Save plot…"
                    enabled: root.hasTable
                    font.pointSize: root.fontPt - 1
                    onClicked: root.saveRequested("rv", rvArea.width, rvArea.height)
                }
            }
            MakieArea {
                id: rvArea
                Layout.fillWidth: true
                Layout.fillHeight: true
                scene: rvPlot
            }
        }
    }
}
