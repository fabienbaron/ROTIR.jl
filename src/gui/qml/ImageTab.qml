// Imaging: a regularised surface reconstruction on the sphere.
//
// NOT OITOOLS' 2-D imaging. The image here is one value per HEALPix tessel on a 3-D surface,
// the Fourier transform is the polygon FT of that surface at each epoch, and the regularisers
// are the ones that mean something on a sphere.
//
// The engine reports by printing, so its trace streams into the console while it runs (see
// `GuiJob` in src/gui/shell.jl) and the map is drawn ONCE, on completion, on the GUI thread.

import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import Makie

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
    property bool jobRunning: false
    property real jobElapsed: 0
    // 0 = the map on the sky, 1 = the whole surface unwrapped.
    property int viewIndex: 0
    // What a run would use, read from Julia — the dataset and the model are chosen elsewhere.
    property string ctxDataset: "—"
    property int    ctxEpochs: 0
    property string ctxModel: "—"
    property string ctxType: "—"
    property string ctxMesh: "—"
    function dp(px) { return Math.round(px * uiScale) }

    signal statusChanged(string s)

    ListModel { id: backendModel }

    function fillBackends() {
        if (backendModel.count > 0) return
        var txt = Julia.shell_polyft_backends()
        if (txt.length === 0) return
        var ls = txt.split("\n")
        for (var i = 0; i < ls.length; ++i) {
            var f = ls[i].split("\t")
            if (f.length < 3) continue
            backendModel.append({ bkey: f[0], blabel: f[1], bdoc: f[2] })
        }
        var cur = Julia.shell_polyft_backend()
        for (var j = 0; j < backendModel.count; ++j)
            if (backendModel.get(j).bkey === cur) backendBox.currentIndex = j
    }
    ListModel { id: regModel }      // every available regulariser, with its default weight
    ListModel { id: imageModel }

    Component.onCompleted: {
        var rows = Julia.shell_regularizer_kinds()
        if (rows.length === 0) return
        var lines = rows.split("\n")
        for (var i = 0; i < lines.length; ++i) {
            var f = lines[i].split("\t")
            if (f.length < 5) continue
            regModel.append({ rname: f[0], weight: f[1], extraLabel: f[2],
                              extra: f[3], rdoc: f[4], on: false })
        }
    }

    function refresh() {
        root.fillBackends()
        var ctx = Julia.shell_imaging_context().split("\t")
        if (ctx.length >= 5) {
            root.ctxDataset = ctx[0]; root.ctxEpochs = parseInt(ctx[1])
            root.ctxModel = ctx[2];   root.ctxType = ctx[3]; root.ctxMesh = ctx[4]
        }
        var rows = Julia.shell_images()
        imageModel.clear()
        if (rows.length === 0) return
        var lines = rows.split("\n")
        for (var i = 0; i < lines.length; ++i) {
            var f = lines[i].split("\t")
            if (f.length < 5) continue
            imageModel.append({ iname: f[0], nexp: f[1], chi2r: f[2], nd: f[3], src: f[4] })
        }
    }

    // The regulariser string the shell parses: `name:weight:a:b` joined by `;`. Built here
    // from the ticked rows and parsed in one place on the Julia side, so the list of what
    // exists has a single authority.
    function regSpec() {
        var parts = []
        for (var i = 0; i < regModel.count; ++i) {
            var r = regModel.get(i)
            if (!r.on) continue
            // name:weight:extra — three fields, because the panel can only choose three.
            // The operators and index sets go with them are built in Julia; see
            // `build_regularizers`.
            parts.push(r.rname + ":" + r.weight + ":" + r.extra)
        }
        return parts.join(";")
    }

    RowLayout {
        anchors.fill: parent
        spacing: dp(8)

        ColumnLayout {
            // fillWidth is FALSE explicitly: a nested layout defaults it to true,
            // and this column would then take the whole row and squeeze the canvas
            // beside it to zero width — which reads as a plot that failed to draw.
            Layout.fillWidth: false
            Layout.preferredWidth: dp(520)
            Layout.minimumWidth: dp(520)
            Layout.maximumWidth: dp(520)
            Layout.fillHeight: true
            spacing: dp(4)

            // What this run would actually use. Both choices are made on OTHER tabs, and
            // without them here a reconstruction started after switching datasets looks
            // exactly like one started before.
            Frame {
                Layout.fillWidth: true
                padding: dp(6)
                ColumnLayout {
                    anchors.fill: parent
                    spacing: dp(2)
                    Label {
                        Layout.fillWidth: true
                        elide: Text.ElideMiddle
                        font.pointSize: root.fontPt
                        color: root.ctxDataset === "—" ? "#c0392b" : "#333"
                        text: "data:  " + root.ctxDataset +
                              (root.ctxEpochs > 0 ? "  (" + root.ctxEpochs + " epoch" +
                                                    (root.ctxEpochs === 1 ? ")" : "s)") : "")
                    }
                    Label {
                        Layout.fillWidth: true
                        elide: Text.ElideMiddle
                        font.pointSize: root.fontPt
                        color: root.ctxModel === "—" ? "#c0392b" : "#333"
                        text: "model: " + root.ctxModel +
                              (root.ctxType === "—" ? "" : "  (" + root.ctxType + ")")
                    }
                    Label {
                        Layout.fillWidth: true
                        elide: Text.ElideRight
                        font.pointSize: root.fontPt - 1
                        color: "#7f8c98"
                        text: "mesh:  " + root.ctxMesh
                    }
                }
            }

            Label { text: "Regularisers"; font.bold: true; font.pointSize: root.fontPt }

            Frame {
                Layout.fillWidth: true
                Layout.preferredHeight: dp(340)
                ListView {
                    id: regList
                    anchors.fill: parent
                    clip: true
                    model: regModel
                    ScrollBar.vertical: ScrollBar {}
                    delegate: RowLayout {
                        width: regList.width
                        spacing: dp(6)
                        CheckBox {
                            checked: on
                            font.pointSize: root.fontPt
                            onToggled: regModel.setProperty(index, "on", checked)
                        }
                        Label {
                            Layout.preferredWidth: dp(80)
                            text: rname
                            font.pointSize: root.fontPt
                            ToolTip.text: rdoc
                            ToolTip.visible: rma.containsMouse
                            MouseArea { id: rma; anchors.fill: parent; hoverEnabled: true }
                        }
                        TextField {
                            Layout.preferredWidth: dp(80)
                            text: weight
                            font.pointSize: root.fontPt - 1
                            selectByMouse: true
                            ToolTip.text: "weight"
                            ToolTip.visible: hovered
                            onEditingFinished: regModel.setProperty(index, "weight", text)
                        }
                        // The one extra scalar, where there is one: B for `bias`, nbins for
                        // the radial pair. Hidden rather than disabled on the other eight
                        // rows, which would otherwise carry a permanently greyed box.
                        Label {
                            visible: extraLabel.length > 0
                            text: extraLabel
                            color: "#7f8c98"
                            font.pointSize: root.fontPt - 1
                        }
                        TextField {
                            visible: extraLabel.length > 0
                            Layout.preferredWidth: dp(52)
                            text: extra
                            font.pointSize: root.fontPt - 1
                            selectByMouse: true
                            onEditingFinished: regModel.setProperty(index, "extra", text)
                        }
                        Label {
                            Layout.fillWidth: true
                            text: rdoc
                            elide: Text.ElideRight
                            color: "#7f8c98"
                            font.pointSize: root.fontPt - 1
                        }
                    }
                }
            }

            RowLayout {
                Layout.fillWidth: true
                spacing: dp(6)
                Label { text: "HEALPix level"; font.pointSize: root.fontPt }
                SpinBox {
                    id: nsideBox
                    // 3 is the practical minimum, the maximum is 6, and 3 is also the
                    // DEFAULT — the same level the model tab starts at, so a reconstruction
                    // sits on the mesh the preview was drawn on unless it is told otherwise.
                    from: 2; to: 6; value: 3
                    font.pointSize: root.fontPt - 1
                    ToolTip.text: "nside = 2^level; npix = 12·nside²"
                    ToolTip.visible: hovered
                }
                Label {
                    text: "= " + (12 * Math.pow(Math.pow(2, nsideBox.value), 2)) + " tessels"
                    color: "#7f8c98"
                    font.pointSize: root.fontPt - 1
                }
                Item { Layout.fillWidth: true }
            }

            // Which forward kernel evaluates the polygon Fourier transform. A COMPUTE choice,
            // not a display one, so it sits with the mesh: those two together are what a χ²
            // costs. Both backends give the exact transform and are asserted against each
            // other in the suite; the vectorised one is 17-19x faster and is the default.
            RowLayout {
                Layout.fillWidth: true
                spacing: dp(6)
                Label { text: "kernel"; color: "#7f8c98"; font.pointSize: root.fontPt - 1 }
                ComboBox {
                    id: backendBox
                    Layout.fillWidth: true
                    model: backendModel
                    textRole: "blabel"
                    font.pointSize: root.fontPt - 1
                    ToolTip.text: currentIndex >= 0 && backendModel.count > 0
                                  ? backendModel.get(currentIndex).bdoc : ""
                    ToolTip.visible: hovered
                    onActivated: root.statusChanged(Julia.shell_set_polyft_backend(
                        backendModel.get(currentIndex).bkey))
                }
            }

            RowLayout {
                Layout.fillWidth: true
                spacing: dp(6)
                Label { text: "optimizer"; font.pointSize: root.fontPt }
                ComboBox {
                    id: optimizerBox
                    Layout.preferredWidth: dp(220)
                    // One entry, and it is the truth: `image_reconstruct_oi` drives
                    // OptimPackNextGen's VMLMB and nothing else. Listed rather than omitted so
                    // that "which optimizer is this" has an answer on screen, and so adding a
                    // second one later is a line here rather than a new control.
                    model: ["VMLMB (OptimPackNextGen)"]
                    font.pointSize: root.fontPt - 1
                    ToolTip.text: "bound-constrained limited-memory quasi-Newton; the only " +
                                  "one image_reconstruct_oi implements"
                    ToolTip.visible: hovered
                }
                Item { Layout.fillWidth: true }
            }

            RowLayout {
                Layout.fillWidth: true
                spacing: dp(6)
                Label { text: "iterations"; font.pointSize: root.fontPt }
                SpinBox {
                    id: iterBox
                    from: 1; to: 100000; stepSize: 50; value: 200
                    editable: true
                    font.pointSize: root.fontPt - 1
                    ToolTip.text: "VMLMB iterations; the trace streams to the console"
                    ToolTip.visible: hovered
                }
                Item { Layout.fillWidth: true }
            }

            RowLayout {
                Layout.fillWidth: true
                spacing: dp(6)
                Button {
                    text: "Reconstruct"
                    enabled: !root.jobRunning
                    font.pointSize: root.fontPt
                    onClicked: root.statusChanged(
                        Julia.shell_reconstruct(nsideBox.value, root.regSpec(), iterBox.value))
                }
                Button {
                    text: "Stop"
                    enabled: root.jobRunning
                    font.pointSize: root.fontPt
                    onClicked: root.statusChanged(Julia.shell_job_stop())
                }
                Item { Layout.fillWidth: true }
                Label {
                    visible: root.jobRunning
                    text: "running " + root.jobElapsed.toFixed(1) + " s"
                    color: "#7f8c98"
                    font.pointSize: root.fontPt - 1
                }
            }

            Label { text: "Maps"; font.bold: true; font.pointSize: root.fontPt }
            Frame {
                Layout.fillWidth: true
                Layout.fillHeight: true
                ListView {
                    id: imageList
                    anchors.fill: parent
                    clip: true
                    model: imageModel
                    ScrollBar.vertical: ScrollBar {}
                    delegate: RowLayout {
                        width: imageList.width
                        spacing: dp(6)
                        Label { Layout.preferredWidth: dp(66); text: iname;  font.pointSize: root.fontPt }
                        Label { Layout.preferredWidth: dp(44); text: "2^" + nexp; font.pointSize: root.fontPt - 1; color: "#7f8c98" }
                        Label { Layout.preferredWidth: dp(80); text: "χ²ᵣ " + chi2r; font.pointSize: root.fontPt }
                        Label { Layout.fillWidth: true; text: src; elide: Text.ElideMiddle; color: "#7f8c98"; font.pointSize: root.fontPt - 1 }
                    }
                }
            }

            Label {
                Layout.fillWidth: true
                wrapMode: Text.WordWrap
                color: "#7f8c98"
                font.pointSize: root.fontPt - 1
                text: "The reconstruction sits on the current model's geometry — set the shape " +
                      "on the Model tab first. The engine's trace appears in the console as it runs."
            }
        }

        ColumnLayout {
            Layout.fillWidth: true
            Layout.fillHeight: true
            spacing: dp(4)

            // THE SAME BAR THE MODEL TAB HAS, minus the posterior — a reconstruction is a
            // map on a star, and there is no reason to look at it differently from the map
            // the model generates, which is the comparison being made most of the time. The
            // decoration ticks and the colormap buttons act on every panel at once, so they
            // need no wiring of their own here.
            RowLayout {
                Layout.fillWidth: true
                spacing: dp(6)
                ViewBar {
                    id: viewBar
                    Layout.fillWidth: true
                    fontPt: root.fontPt
                    uiScale: root.uiScale
                    showPosterior: false
                    // A reconstruction is fitted to every epoch at once, so there is no epoch
                    // to step through here.
                    showEpoch: false
                    onViewChanged: { root.viewIndex = viewBar.viewIndex; root.redraw() }
                    onStatusChanged: function (s) { root.statusChanged(s) }
                }
                Button {
                    text: "Save view…"
                    Layout.alignment: Qt.AlignTop
                    font.pointSize: root.fontPt - 1
                    onClicked: {
                        var w = ["imaging", "imaging_moll", "imaging_3d"][root.viewIndex]
                        var m = Julia.shell_save_figure(w, "rotir_" + w + ".png",
                                                        imSkyArea.width, imSkyArea.height)
                        root.statusChanged(m.length > 0 ? m : "saved rotir_" + w + ".png")
                    }
                }
            }

            StackLayout {
                Layout.fillWidth: true
                Layout.fillHeight: true
                currentIndex: root.viewIndex
                MakieArea { id: imSkyArea;  scene: imSkyPlot }
                MakieArea { id: imMollArea; scene: imMollPlot }
                MakieArea { id: imStarArea; scene: imStarPlot }
            }
        }
    }

    function redraw() {
        if (root.viewIndex === 0)      imSkyArea.update()
        else if (root.viewIndex === 1) imMollArea.update()
        else                           imStarArea.update()
    }
}
