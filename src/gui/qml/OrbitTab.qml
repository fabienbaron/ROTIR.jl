// Orbit: the elements of a binary, a fit, and a picture of the orbit on the sky.
//
// TWO BOXES, NOT ONE. The elements say WHERE the secondary is; the star model says WHAT sits
// at the two positions — analytic 2-D profiles, or ROTIR's tessellated 3-D components. Mixing
// them is how `q` and `rpole` end up looking like orbital elements, which they are not: the
// relative astrometric orbit does not constrain a mass ratio.
//
// The element rows are laid out on FIXED columns rather than a plain row of items, so name,
// value, state and bounds line up down the list whatever state each row is in. A row whose
// bounds appear only when it is free shifts every field on that line otherwise, and a column
// of parameters that jitters as you flip a combo is unreadable.

import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import Makie

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
    property bool jobRunning: false
    function dp(px) { return Math.round(px * uiScale) }

    signal statusChanged(string s)
    signal pickFile(string mode)          // the window owns the picker

    ListModel { id: paramModel }
    ListModel { id: kindModel }
    ListModel { id: optionModel }
    ListModel { id: starModelModel }

    property string starModel: "analytic"

    // The element table's columns. One place, so the header and every row agree.
    readonly property int wName:  dp(112)   // "domega (deg/d)" fits without eliding
    readonly property int wValue: dp(96)
    readonly property int wState: dp(76)
    readonly property int wBound: dp(66)

    function fillStarModels() {
        if (starModelModel.count > 0) return
        var txt = Julia.shell_orbit_star_models()
        if (txt.length === 0) return
        var ls = txt.split("\n")
        for (var i = 0; i < ls.length; ++i) {
            var f = ls[i].split("\t")
            if (f.length < 3) continue
            starModelModel.append({ key: f[0], label: f[1], doc: f[2] })
        }
    }

    function fillKinds() {
        if (kindModel.count > 0) return
        var txt = Julia.shell_orbit_component_kinds()
        if (txt.length === 0) return
        var ls = txt.split("\n")
        for (var i = 0; i < ls.length; ++i) {
            var f = ls[i].split("\t")
            if (f.length < 3) continue
            kindModel.append({ key: f[0], label: f[1] })
        }
    }

    function refresh() {
        fillKinds()
        fillStarModels()
        var rows = Julia.shell_orbit_params()
        paramModel.clear()
        if (rows.length > 0) {
            var ls = rows.split("\n")
            for (var i = 0; i < ls.length; ++i) {
                var f = ls[i].split("\t")
                if (f.length < 8) continue
                paramModel.append({ pname: f[0], punit: f[1], pvalue: f[2], pstate: f[3],
                                    plo: f[4], phi: f[5], ptie: f[6], pdoc: f[7] })
            }
        }
        var opts = Julia.shell_orbit_options()
        optionModel.clear()
        if (opts.length > 0) {
            var ol = opts.split("\n")
            for (var j = 0; j < ol.length; ++j) {
                var g = ol[j].split("\t")
                if (g.length < 4) continue
                optionModel.append({ oname: g[0], olabel: g[1], oon: g[2] === "1", odoc: g[3] })
            }
        }
        root.starModel = Julia.shell_orbit_star_model()
        for (var s = 0; s < modelBox.count; ++s)
            if (modelBox.model.get(s).key === root.starModel) modelBox.currentIndex = s
        var comps = Julia.shell_orbit_components().split("\t")
        if (comps.length >= 2) {
            c1Box.currentIndex = root.kindIndex(comps[0])
            c2Box.currentIndex = root.kindIndex(comps[1])
        }
        var rp = Julia.shell_orbit_render_params().split("\t")
        if (rp.length >= 5) {
            rpole1.text = rp[0]; rpole2.text = rp[1]
            tpole1.text = rp[2]; tpole2.text = rp[3]; qField.text = rp[4]
        }
        resultTable.text = Julia.shell_orbit_result()
    }

    function kindIndex(k) {
        for (var i = 0; i < kindModel.count; ++i)
            if (kindModel.get(i).key === k) return i
        return 0
    }

    function redraw() { orbitArea.update() }

    contentItem: RowLayout {
        spacing: dp(8)

        // ── left: the elements, the components, the switches ─────────────────
        ColumnLayout {
            Layout.fillWidth: false
            Layout.preferredWidth: dp(540)
            Layout.minimumWidth: dp(540)
            Layout.maximumWidth: dp(540)
            Layout.fillHeight: true
            spacing: dp(4)

            Label { text: "Orbital elements"; font.bold: true; font.pointSize: root.fontPt }

            // The column header, which is also what fixes the widths every row aligns to.
            RowLayout {
                Layout.fillWidth: true
                Layout.leftMargin: dp(6)
                Layout.rightMargin: dp(6)
                spacing: dp(4)
                Label { Layout.preferredWidth: root.wName; text: "element"
                        color: "#7f8c98"; font.pointSize: root.fontPt - 2 }
                Label { Layout.preferredWidth: root.wValue; text: "value"
                        color: "#7f8c98"; font.pointSize: root.fontPt - 2 }
                Label { Layout.preferredWidth: root.wState; text: "state"
                        color: "#7f8c98"; font.pointSize: root.fontPt - 2 }
                Label { Layout.fillWidth: true; text: "bounds, or the tie expression"
                        color: "#7f8c98"; font.pointSize: root.fontPt - 2 }
            }

            Frame {
                Layout.fillWidth: true
                Layout.fillHeight: true
                ListView {
                    id: paramList
                    anchors.fill: parent
                    clip: true
                    model: paramModel
                    spacing: dp(2)
                    ScrollBar.vertical: ScrollBar {}
                    // Fixed columns: the trailing slot always occupies the same space, and
                    // swaps its CONTENTS between bounds and a tie expression, so nothing to
                    // its left ever moves.
                    delegate: RowLayout {
                        width: paramList.width
                        spacing: dp(4)
                        Label {
                            Layout.preferredWidth: root.wName
                            Layout.alignment: Qt.AlignVCenter
                            text: pname + (punit.length > 0 ? " (" + punit + ")" : "")
                            elide: Text.ElideRight
                            font.pointSize: root.fontPt
                            ToolTip.text: pdoc
                            ToolTip.visible: oma.containsMouse && pdoc.length > 0
                            MouseArea { id: oma; anchors.fill: parent; hoverEnabled: true }
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
                                root.statusChanged(
                                    Julia.shell_set_orbit_state(pname, model[currentIndex]))
                                root.refresh()
                            }
                        }
                        // ANCHORED, not a nested RowLayout: a layout with only some of its
                        // items visible redistributes the slack, so the upper bound drifted
                        // right by fifty pixels and no two rows agreed on where it sat.
                        Item {
                            Layout.fillWidth: true
                            Layout.preferredHeight: oLo.implicitHeight
                            TextField {
                                id: oLo
                                visible: pstate === "free"
                                anchors.left: parent.left
                                width: root.wBound
                                text: plo
                                font.pointSize: root.fontPt - 2
                                selectByMouse: true
                                ToolTip.text: "lower bound"; ToolTip.visible: hovered
                                onEditingFinished: root.statusChanged(
                                    Julia.shell_set_orbit_bound(pname, text, oHi.text))
                            }
                            TextField {
                                id: oHi
                                visible: pstate === "free"
                                anchors.left: parent.left
                                anchors.leftMargin: root.wBound + dp(4)
                                width: root.wBound
                                text: phi
                                font.pointSize: root.fontPt - 2
                                selectByMouse: true
                                ToolTip.text: "upper bound"; ToolTip.visible: hovered
                                onEditingFinished: root.statusChanged(
                                    Julia.shell_set_orbit_bound(pname, oLo.text, text))
                            }
                            TextField {
                                visible: pstate === "tied"
                                anchors.left: parent.left
                                anchors.right: parent.right
                                text: ptie
                                placeholderText: "expression, e.g. -Omega"
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
            }

            // ── the star model, which is a DIFFERENT frame from the elements ─
            GroupBox {
                Layout.fillWidth: true
                font.pointSize: root.fontPt
                label: Label {
                    text: "Star model"
                    font.bold: true
                    font.pointSize: root.fontPt
                }

                ColumnLayout {
                    anchors.fill: parent
                    spacing: dp(4)

                    RowLayout {
                        Layout.fillWidth: true
                        spacing: dp(6)
                        Label { text: "components"; color: "#7f8c98"
                                font.pointSize: root.fontPt - 1 }
                        ComboBox {
                            id: modelBox
                            Layout.fillWidth: true
                            model: starModelModel
                            textRole: "label"
                            font.pointSize: root.fontPt - 1
                            ToolTip.text: currentIndex >= 0 && starModelModel.count > 0
                                ? starModelModel.get(currentIndex).doc : ""
                            ToolTip.visible: hovered
                            onActivated: {
                                root.statusChanged(Julia.shell_set_orbit_star_model(
                                    starModelModel.get(currentIndex).key))
                                root.refresh(); root.redraw()
                            }
                        }
                    }

                    // 2-D: which analytic profile each component is.
                    RowLayout {
                        Layout.fillWidth: true
                        spacing: dp(6)
                        enabled: root.starModel === "analytic"
                        opacity: enabled ? 1 : 0.45
                        Label { text: "primary"; color: "#7f8c98"
                                font.pointSize: root.fontPt - 1 }
                        ComboBox {
                            id: c1Box
                            Layout.fillWidth: true
                            model: kindModel
                            textRole: "label"
                            font.pointSize: root.fontPt - 1
                            onActivated: {
                                root.statusChanged(Julia.shell_set_orbit_component(
                                    "1", kindModel.get(currentIndex).key))
                                root.refresh(); root.redraw()
                            }
                        }
                        Label { text: "secondary"; color: "#7f8c98"
                                font.pointSize: root.fontPt - 1 }
                        ComboBox {
                            id: c2Box
                            Layout.fillWidth: true
                            model: kindModel
                            textRole: "label"
                            font.pointSize: root.fontPt - 1
                            onActivated: {
                                root.statusChanged(Julia.shell_set_orbit_component(
                                    "2", kindModel.get(currentIndex).key))
                                root.refresh(); root.redraw()
                            }
                        }
                    }

                    // 3-D: the surface parameters, which no analytic profile has. These are
                    // star-model quantities — `q` above all, which the relative orbit does
                    // not constrain and which therefore cannot be an element.
                    GridLayout {
                        columns: 6
                        Layout.fillWidth: true
                        columnSpacing: dp(5)
                        rowSpacing: dp(3)
                        // Greyed under the analytic model for the same reason the profile
                        // combos are greyed under the tessellated one: a radius and a pole
                        // temperature are things only a surface has, and a field that reads
                        // as editable but changes nothing is the frame-mixing this tab is
                        // laid out to avoid.
                        enabled: root.starModel === "tessellated"
                        opacity: enabled ? 1 : 0.45
                        Label { text: "R₁ (mas)"; color: "#7f8c98"; font.pointSize: root.fontPt - 1 }
                        TextField { id: rpole1; Layout.preferredWidth: dp(64); font.pointSize: root.fontPt - 1
                                    selectByMouse: true
                                    onEditingFinished: { Julia.shell_set_orbit_render_param("rpole1", text); root.redraw() } }
                        Label { text: "R₂ (mas)"; color: "#7f8c98"; font.pointSize: root.fontPt - 1 }
                        TextField { id: rpole2; Layout.preferredWidth: dp(64); font.pointSize: root.fontPt - 1
                                    selectByMouse: true
                                    onEditingFinished: { Julia.shell_set_orbit_render_param("rpole2", text); root.redraw() } }
                        Label { text: "q"; color: "#7f8c98"; font.pointSize: root.fontPt - 1
                                ToolTip.text: "mass ratio, an ASSUMPTION — the relative orbit does not constrain it"
                                ToolTip.visible: qma.containsMouse
                                MouseArea { id: qma; anchors.fill: parent; hoverEnabled: true } }
                        TextField { id: qField; Layout.preferredWidth: dp(64); font.pointSize: root.fontPt - 1
                                    selectByMouse: true
                                    onEditingFinished: { Julia.shell_set_orbit_render_param("q", text); root.redraw() } }
                        Label { text: "T₁ (K)"; color: "#7f8c98"; font.pointSize: root.fontPt - 1 }
                        TextField { id: tpole1; Layout.preferredWidth: dp(64); font.pointSize: root.fontPt - 1
                                    selectByMouse: true
                                    onEditingFinished: { Julia.shell_set_orbit_render_param("tpole1", text); root.redraw() } }
                        Label { text: "T₂ (K)"; color: "#7f8c98"; font.pointSize: root.fontPt - 1 }
                        TextField { id: tpole2; Layout.preferredWidth: dp(64); font.pointSize: root.fontPt - 1
                                    selectByMouse: true
                                    onEditingFinished: { Julia.shell_set_orbit_render_param("tpole2", text); root.redraw() } }
                        Item {}
                        Item {}
                    }

                    Flow {
                        Layout.fillWidth: true
                        spacing: dp(8)
                        Repeater {
                            model: optionModel
                            CheckBox {
                                text: olabel
                                checked: oon
                                font.pointSize: root.fontPt - 1
                                ToolTip.text: odoc
                                ToolTip.visible: hovered
                                onToggled: {
                                    root.statusChanged(
                                        Julia.shell_set_orbit_option(oname, checked ? "1" : "0"))
                                    root.redraw()
                                }
                            }
                        }
                    }
                }
            }

            // ── load / save / fit ────────────────────────────────────────────
            RowLayout {
                Layout.fillWidth: true
                spacing: dp(6)
                Button {
                    text: "Load orbit…"
                    font.pointSize: root.fontPt - 1
                    ToolTip.text: "a TOML written by Save: elements, bounds, ties, components"
                    ToolTip.visible: hovered
                    onClicked: root.pickFile("orbit")
                }
                Button {
                    text: "Save orbit"
                    font.pointSize: root.fontPt - 1
                    onClicked: root.statusChanged(Julia.shell_save_orbit("rotir_orbit.toml"))
                }
                Item { Layout.fillWidth: true }
            }

            ComboBox {
                id: orbitMethod
                Layout.fillWidth: true
                model: ["neldermead", "nautilus"]
                font.pointSize: root.fontPt - 1
            }
            RowLayout {
                Layout.fillWidth: true
                spacing: dp(6)
                Label { text: "evaluations"; color: "#7f8c98"; font.pointSize: root.fontPt - 1 }
                SpinBox {
                    id: orbitEval
                    Layout.preferredWidth: dp(120)
                    from: 100; to: 500000; stepSize: 1000; value: 20000
                    editable: true
                    font.pointSize: root.fontPt - 1
                }
                Item { Layout.fillWidth: true }
                Button {
                    text: "Fit orbit"
                    enabled: !root.jobRunning
                    font.pointSize: root.fontPt
                    onClicked: root.statusChanged(
                        Julia.shell_fit_orbit(orbitMethod.currentText, orbitEval.value))
                }
                Button {
                    text: "Stop"
                    enabled: root.jobRunning
                    font.pointSize: root.fontPt
                    onClicked: root.statusChanged(Julia.shell_job_stop())
                }
            }

            Frame {
                Layout.fillWidth: true
                Layout.preferredHeight: dp(90)
                Text {
                    id: resultTable
                    anchors.fill: parent
                    font.family: "monospace"
                    font.pointSize: root.fontPt - 1
                    text: ""
                }
            }
        }

        // ── right: the orbit ─────────────────────────────────────────────────
        ColumnLayout {
            Layout.fillWidth: true
            Layout.fillHeight: true
            spacing: dp(4)
            RowLayout {
                Layout.fillWidth: true
                spacing: dp(6)
                Label {
                    text: "the relative orbit · × is the primary · numbered marks are the " +
                          "observed epochs"
                    color: "#7f8c98"
                    font.pointSize: root.fontPt - 1
                }
                Item { Layout.fillWidth: true }
                Button {
                    text: "Save view…"
                    font.pointSize: root.fontPt - 1
                    onClicked: {
                        var m = Julia.shell_save_figure("orbit", "rotir_orbit.png",
                                                        orbitArea.width, orbitArea.height)
                        root.statusChanged(m.length > 0 ? m : "saved rotir_orbit.png")
                    }
                }
            }
            MakieArea {
                id: orbitArea
                Layout.fillWidth: true
                Layout.fillHeight: true
                scene: orbitPlot
            }
        }
    }
}
