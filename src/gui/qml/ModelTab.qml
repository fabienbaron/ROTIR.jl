// Model: a 3-D surface, its parameters, and a fit.
//
// The form is GENERATED from the schema (src/surface_schema.jl) rather than written out here.
// That is the point of the schema existing: the fields a surface type needs, their units,
// their plausible ranges and which limb-darkening coefficients a given law actually reads are
// all one declaration, so the form cannot offer a field `compute_radii` does not read, and
// adding a surface type does not mean editing QML.
//
// Each parameter is free / fixed / tied. Three-way, not a tick box: a tick cannot express
// "expression", and a tie is exactly an expression in the other parameters.

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
    property string jobProgress: ""
    // How many epochs the current dataset has, for the slider in the view bar.
    property int epochCount: 0
    // Set from the settings panel; the mesh is rebuilt in this float type.
    property string precision: "Float32"
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
    ListModel { id: paramModel }
    ListModel { id: modelListModel }
    ListModel { id: typeModel }
    ListModel { id: methodModel }
    ListModel { id: fitModel }              // every completed fit, for comparison
    ListModel { id: fitParamModel }         // the selected fit's free parameters

    // Filled on the first `refresh()` rather than in `Component.onCompleted`. A tab's
    // onCompleted runs while the QML tree is still being built, and the combos came up empty;
    // populating them from the refresh Main.qml triggers once the window is up is one less
    // ordering assumption to get wrong.
    function fillChoices() {
        if (typeModel.count > 0) return
        var rows = Julia.shell_surface_types()
        if (rows.length > 0) {
            var lines = rows.split("\n")
            for (var i = 0; i < lines.length; ++i) {
                var f = lines[i].split("\t")
                if (f.length < 4) continue
                typeModel.append({ code: f[0], name: f[1], label: f[2], doc: f[3] })
            }
        }
        var mrows = Julia.shell_fit_methods()
        if (mrows.length > 0) {
            var mlines = mrows.split("\n")
            for (var j = 0; j < mlines.length; ++j) {
                var g = mlines[j].split("\t")
                if (g.length < 2) continue
                methodModel.append({ key: g[0], label: g[1] })
            }
        }
        // A ComboBox whose model was EMPTY at construction keeps currentIndex = -1 and renders
        // blank when rows arrive later; it only defaults to 0 for a model that already had
        // rows. Setting it explicitly is the whole fix.
        if (typeBox.currentIndex < 0 && typeModel.count > 0) typeBox.currentIndex = 0
        if (methodBox.currentIndex < 0 && methodModel.count > 0) methodBox.currentIndex = 0
        viewBar.loadColormaps()
        viewBar.loadDecorations()
        viewBar.loadField()
    }

    function pushTessellation() {
        root.statusChanged(Julia.shell_set_tessellation(
            "healpix", nsideBox.value, root.precision))
        root.redraw()
    }

    function refresh() {
        root.fillBackends()
        fillChoices()
        root.refreshFits()
        // The model list. A session can hold several — a binary is two components — and
        // without a selector the second one could be created but never reached again.
        var mrows = Julia.shell_models()
        modelListModel.clear()
        if (mrows.length > 0) {
            var mlines = mrows.split("\n")
            for (var k = 0; k < mlines.length; ++k) {
                var h = mlines[k].split("\t")
                if (h.length < 3) continue
                modelListModel.append({ mname: h[0], mtype: h[1], mfree: h[2] })
            }
        }
        var rows = Julia.shell_params()
        paramModel.clear()
        if (rows.length > 0) {
            var lines = rows.split("\n")
            for (var i = 0; i < lines.length; ++i) {
                var f = lines[i].split("\t")
                if (f.length < 12) continue
                paramModel.append({ pname: f[0], plabel: f[1], punit: f[2], pvalue: f[3],
                                    pstate: f[4], plo: f[5], phi: f[6], ptie: f[7],
                                    pgroup: f[8], pkind: f[9], pchoices: f[10], pdoc: f[11] })
            }
        }
        warnLabel.text = Julia.shell_validate_model()
        fitTable.text = Julia.shell_last_fit()
        var ep = Julia.shell_epochs()
        root.epochCount = ep.length > 0 ? ep.split("\n").length : 0
    }

    RowLayout {
        anchors.fill: parent
        spacing: dp(8)

        // ── left: the generated parameter form ───────────────────────────────
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

            // ONE model at a time, so this is a label rather than a selector. A session with
            // several was offering a choice nobody was making, while the unselected one still
            // decided what the χ² column and the reconstruction were about.
            Label {
                Layout.fillWidth: true
                elide: Text.ElideRight
                font.pointSize: root.fontPt
                color: modelListModel.count > 0 ? "#333" : "#7f8c98"
                text: modelListModel.count === 0 ? "no model — pick a surface type below"
                    : "model: " + modelListModel.get(0).mname + "  (type " +
                      modelListModel.get(0).mtype + ", " + modelListModel.get(0).mfree + " free)"
            }

            RowLayout {
                Layout.fillWidth: true
                spacing: dp(6)
                ComboBox {
                    id: typeBox
                    Layout.fillWidth: true
                    model: typeModel
                    textRole: "label"
                    font.pointSize: root.fontPt
                    ToolTip.text: currentIndex >= 0 ? typeModel.get(currentIndex).doc : ""
                    ToolTip.visible: hovered && ToolTip.text.length > 0
                }
                Button {
                    text: "+ model"
                    font.pointSize: root.fontPt
                    ToolTip.text: "replace the current model with a new one of this type"
                    ToolTip.visible: hovered
                    onClicked: {
                        if (typeBox.currentIndex < 0) return
                        root.statusChanged(
                            Julia.shell_add_model(parseInt(typeModel.get(typeBox.currentIndex).code)))
                        root.refresh()
                        root.redraw()
                    }
                }
                Button {
                    text: "− model"
                    enabled: modelListModel.count > 0
                    font.pointSize: root.fontPt
                    ToolTip.text: "clear it: the surface views go idle and the epoch table " +
                                  "falls back to point counts"
                    ToolTip.visible: hovered
                    onClicked: {
                        root.statusChanged(Julia.shell_clear_model())
                        root.refresh()
                        root.redraw()
                    }
                }
            }

            // The mesh the model is evaluated on. It belongs here rather than in the settings
            // panel because it decides what is being FITTED — a level-3 sphere and a level-5
            // one are different models — not how anything is drawn.
            RowLayout {
                Layout.fillWidth: true
                spacing: dp(6)
                Label { text: "mesh"; color: "#7f8c98"; font.pointSize: root.fontPt - 1 }
                ComboBox {
                    id: tesselBox
                    Layout.preferredWidth: dp(150)
                    model: ["HEALPix", "long-lat"]
                    font.pointSize: root.fontPt - 1
                    // long-lat is offered so the choice is visible, and disabled because it is
                    // not wired: `create_star` builds it, but every regulariser and every
                    // shape gradient is written against the HEALPix neighbour structure.
                    delegate: ItemDelegate {
                        width: tesselBox.width
                        text: modelData
                        enabled: index === 0
                        font.pointSize: root.fontPt - 1
                    }
                    onActivated: {
                        if (currentIndex !== 0) { currentIndex = 0; return }
                        root.pushTessellation()
                    }
                }
                Label { text: "n"; color: "#7f8c98"; font.pointSize: root.fontPt - 1 }
                SpinBox {
                    id: nsideBox
                    // 3 is the default and the practical minimum; 6 the maximum, past which
                    // the polygon FT dominates everything the GUI does.
                    from: 2; to: 6; value: 3
                    font.pointSize: root.fontPt - 1
                    ToolTip.text: "nside = 2^n; npix = 12·nside²"
                    ToolTip.visible: hovered
                    onValueModified: root.pushTessellation()
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
                    // Grouped by section header. `pgroup` comes from the schema, so the order
                    // of the sections is the order the schema lists the fields in — geometry,
                    // then thermal, then limb darkening, then orientation, then orbit.
                    section.property: "pgroup"
                    section.delegate: Label {
                        width: paramList.width
                        text: section
                        font.bold: true
                        font.pointSize: root.fontPt - 1
                        color: "#7f8c98"
                        topPadding: dp(6)
                    }
                    delegate: RowLayout {
                        width: paramList.width
                        spacing: dp(4)

                        Label {
                            Layout.preferredWidth: dp(118)
                            text: plabel + (punit.length > 0 ? " (" + punit + ")" : "")
                            elide: Text.ElideRight
                            font.pointSize: root.fontPt
                            ToolTip.text: pdoc
                            ToolTip.visible: ma.containsMouse && pdoc.length > 0
                            MouseArea { id: ma; anchors.fill: parent; hoverEnabled: true }
                        }

                        // A choice field (`ldtype`) renders as a combo, everything else as a
                        // number. The schema says which; QML does not decide.
                        Loader {
                            Layout.preferredWidth: dp(120)
                            sourceComponent: pkind === "choice" ? choiceField : numberField
                        }

                        Component {
                            id: numberField
                            TextField {
                                text: pvalue
                                font.pointSize: root.fontPt
                                selectByMouse: true
                                enabled: pstate !== "tied"
                                onEditingFinished: {
                                    var msg = Julia.shell_set_param(pname, text)
                                    if (msg.length > 0) root.statusChanged(msg)
                                    root.refresh(); root.redraw()
                                }
                            }
                        }
                        Component {
                            id: choiceField
                            ComboBox {
                                model: pchoices.length > 0 ? pchoices.split("|") : []
                                font.pointSize: root.fontPt
                                currentIndex: {
                                    var opts = pchoices.split("|")
                                    for (var i = 0; i < opts.length; ++i)
                                        if (opts[i].split("=")[0] === String(Math.round(parseFloat(pvalue))))
                                            return i
                                    return -1
                                }
                                displayText: currentIndex >= 0
                                    ? model[currentIndex].split("=").slice(1).join("=") : ""
                                delegate: ItemDelegate {
                                    width: parent.width
                                    text: modelData.split("=").slice(1).join("=")
                                    font.pointSize: root.fontPt
                                }
                                onActivated: {
                                    Julia.shell_set_param(pname, model[currentIndex].split("=")[0])
                                    root.refresh(); root.redraw()
                                }
                            }
                        }

                        ComboBox {
                            Layout.preferredWidth: dp(74)
                            model: ["fixed", "free", "tied"]
                            font.pointSize: root.fontPt - 1
                            currentIndex: pstate === "free" ? 1 : pstate === "tied" ? 2 : 0
                            onActivated: {
                                var msg = Julia.shell_set_param_state(pname, model[currentIndex])
                                if (msg.length > 0) root.statusChanged(msg)
                                root.refresh()
                            }
                        }

                        // The tie expression, shown only when the parameter is TIED, and the
                        // bounds, shown only when it is FREE. Hidden rather than disabled: a
                        // greyed-out box on every row triples the width of a form that is
                        // already twenty rows long, and only one of the three states ever
                        // needs an extra control.
                        TextField {
                            Layout.fillWidth: true
                            visible: pstate === "tied"
                            text: ptie
                            placeholderText: "expression, e.g. P or 180 - i"
                            font.pointSize: root.fontPt - 1
                            selectByMouse: true
                            onEditingFinished: {
                                root.statusChanged(Julia.shell_set_tie(pname, text))
                                root.refresh(); root.redraw()
                            }
                        }
                        TextField {
                            id: loField
                            visible: pstate === "free"
                            Layout.preferredWidth: dp(58)
                            text: plo
                            font.pointSize: root.fontPt - 2
                            selectByMouse: true
                            ToolTip.text: "lower bound"
                            ToolTip.visible: hovered
                            onEditingFinished:
                                root.statusChanged(Julia.shell_set_bound(pname, text, hiField.text))
                        }
                        TextField {
                            id: hiField
                            visible: pstate === "free"
                            Layout.preferredWidth: dp(58)
                            text: phi
                            font.pointSize: root.fontPt - 2
                            selectByMouse: true
                            ToolTip.text: "upper bound"
                            ToolTip.visible: hovered
                            onEditingFinished:
                                root.statusChanged(Julia.shell_set_bound(pname, loField.text, text))
                        }
                        // An AT-BOUND indicator: a free parameter sitting on its bound is not
                        // a fit result, it is a fit that was stopped, and the number alone
                        // gives no sign of that.
                        Label {
                            visible: pstate === "free" &&
                                     (Math.abs(parseFloat(pvalue) - parseFloat(plo)) <
                                      1e-6 * Math.max(1, Math.abs(parseFloat(plo))) ||
                                      Math.abs(parseFloat(pvalue) - parseFloat(phi)) <
                                      1e-6 * Math.max(1, Math.abs(parseFloat(phi))))
                            text: "⚠"
                            color: "#c0392b"
                            font.pointSize: root.fontPt
                            ToolTip.text: "at a bound — widen it or the fit is being held here"
                            // A HoverHandler, not `hovered`: that property belongs to Control,
                            // and a plain Label referencing it raises a ReferenceError on
                            // every re-evaluation of the binding.
                            ToolTip.visible: boundHover.hovered
                            HoverHandler { id: boundHover }
                        }
                        Item { Layout.fillWidth: true; visible: pstate !== "tied" }
                    }
                }
            }

            Label {
                id: warnLabel
                Layout.fillWidth: true
                wrapMode: Text.WordWrap
                color: "#c0392b"
                font.pointSize: root.fontPt - 1
                visible: text.length > 0
            }

            // ── the fit launcher ─────────────────────────────────────────────
            // Two rows, not one: the method labels are sentences ("VMLMB + Zygote gradient
            // (fast; rapid rotator only)"), and four controls across this column squeezed the
            // combo to a width that showed none of it.
            ComboBox {
                id: methodBox
                Layout.fillWidth: true
                model: methodModel
                textRole: "label"
                font.pointSize: root.fontPt
            }
            RowLayout {
                Layout.fillWidth: true
                spacing: dp(6)
                Label { text: "evaluations"; font.pointSize: root.fontPt - 1; color: "#7f8c98" }
                SpinBox {
                    id: evalBox
                    Layout.preferredWidth: dp(120)
                    from: 10; to: 200000; stepSize: 100; value: 2000
                    editable: true
                    font.pointSize: root.fontPt - 1
                    ToolTip.text: "objective evaluations (gradient path: /50 -> iterations)"
                    ToolTip.visible: hovered
                }
                Item { Layout.fillWidth: true }
                Button {
                    text: "Fit"
                    enabled: !root.jobRunning
                    font.pointSize: root.fontPt
                    onClicked: {
                        if (methodBox.currentIndex < 0) return
                        root.statusChanged(Julia.shell_fit(
                            methodModel.get(methodBox.currentIndex).key, evalBox.value))
                    }
                }
                Button {
                    text: "Stop"
                    enabled: root.jobRunning
                    font.pointSize: root.fontPt
                    onClicked: root.statusChanged(Julia.shell_job_stop())
                }
            }

            // The same running readout the Imaging tab has. A sampler can run for minutes and
            // a spinner alone does not say whether it is getting anywhere; the elapsed time is
            // what tells you a fit is worth waiting for or worth stopping.
            Label {
                Layout.fillWidth: true
                visible: root.jobRunning
                text: root.jobProgress.length > 0
                      ? "running — " + root.jobProgress + "   (" +
                        root.jobElapsed.toFixed(0) + " s)"
                      : "running — " + root.jobElapsed.toFixed(0) + " s"
                color: "#7f8c98"
                font.pointSize: root.fontPt - 1
                elide: Text.ElideRight
            }

            Frame {
                Layout.fillWidth: true
                Layout.preferredHeight: dp(110)
                Text {
                    id: fitTable
                    anchors.fill: parent
                    font.family: "monospace"
                    font.pointSize: root.fontPt - 1
                    text: ""
                }
            }

            // ── every fit so far, and the evidence that compares them ────────
            //
            // Δχ² CANNOT compare a sphere with a rapid rotator: more parameters always fit
            // better, and χ²ᵣ ≫ 1 here anyway. log(Z) can, because the prior volume a model
            // spends is already in it — and every nested sampler and Pigeons was computing it
            // all along while the panel showed two numbers per parameter and dropped the rest.
            //
            // Fits ACCUMULATE where models do not: fit a sphere, switch surface type, fit
            // again, and the two rows are the comparison.
            Label { text: "Fits"; font.bold: true; font.pointSize: root.fontPt }
            // The header sits INSIDE the frame, above the list, and the rows carry no padding
            // of their own. With the header outside it, the two were offset by the Frame's
            // padding plus the ItemDelegate's — a couple of dozen pixels of style-dependent
            // inset that no width here could have compensated for, since neither number is
            // fixed by this file.
            Frame {
                Layout.fillWidth: true
                Layout.preferredHeight: dp(118)
                ColumnLayout {
                    anchors.fill: parent
                    spacing: dp(2)
                    RowLayout {
                        Layout.fillWidth: true
                        spacing: dp(4)
                        Label { Layout.preferredWidth: dp(96); text: "method"
                                color: "#7f8c98"; font.pointSize: root.fontPt - 2 }
                        Label { Layout.preferredWidth: dp(34); text: "type"
                                color: "#7f8c98"; font.pointSize: root.fontPt - 2 }
                        Label { Layout.preferredWidth: dp(62); text: "χ²ᵣ"
                                color: "#7f8c98"; font.pointSize: root.fontPt - 2
                                horizontalAlignment: Text.AlignRight }
                        Label { Layout.fillWidth: true; text: "log(Z)"
                                color: "#7f8c98"; font.pointSize: root.fontPt - 2
                                horizontalAlignment: Text.AlignRight }
                        Label { Layout.preferredWidth: dp(56); text: "draws"
                                color: "#7f8c98"; font.pointSize: root.fontPt - 2
                                horizontalAlignment: Text.AlignRight }
                    }
                    ListView {
                        id: fitList
                        Layout.fillWidth: true
                        Layout.fillHeight: true
                        clip: true
                        model: fitModel
                        ScrollBar.vertical: ScrollBar {}
                    onCurrentIndexChanged: {
                        if (currentIndex >= 0) {
                            root.statusChanged(Julia.shell_select_fit(currentIndex + 1))
                            root.refreshFits()
                            postArea.update()
                        }
                    }
                        delegate: ItemDelegate {
                        width: fitList.width
                        // Zero horizontal padding: the header above has none either, and the
                        // default inset is what put the columns out of line with it.
                        leftPadding: 0
                        rightPadding: 0
                        highlighted: ListView.isCurrentItem
                        onClicked: fitList.currentIndex = index
                        ToolTip.text: fdiag.length > 0 ? fname + " — " + fdiag : fname
                        ToolTip.visible: hovered && fname.length > 0
                        contentItem: RowLayout {
                            spacing: dp(4)
                            Label { Layout.preferredWidth: dp(96); text: fmethod
                                    font.pointSize: root.fontPt - 1; elide: Text.ElideRight }
                            Label { Layout.preferredWidth: dp(34); text: ftype
                                    font.pointSize: root.fontPt - 1 }
                            Label { Layout.preferredWidth: dp(62); text: fchi2r
                                    font.pointSize: root.fontPt - 1
                                    horizontalAlignment: Text.AlignRight }
                            Label { Layout.fillWidth: true
                                    // The error beside it, because a log(Z) difference smaller
                                    // than the error is not a preference for either model.
                                    text: flogz === "—" ? "—" : flogz + " ± " + flogzerr
                                    font.pointSize: root.fontPt - 1
                                    horizontalAlignment: Text.AlignRight }
                            Label { Layout.preferredWidth: dp(56); text: fns
                                    color: fns === "0" ? "#a0a6ac" : "#333"
                                    font.pointSize: root.fontPt - 1
                                    horizontalAlignment: Text.AlignRight }
                        }
                    }
                    }
                }
            }
        }

        // ── right: one view at a time ────────────────────────────────────────
        //
        // The ORTHOGRAPHIC projection is the default, and deliberately so: it is what the
        // interferometer measures, at the epoch the slider selects, so it is the view a χ²
        // can be reasoned about from. The 3-D scene shows the far side, which no single epoch
        // constrains, and the Mollweide shows the whole surface at once — both are a step
        // away from the data and are reached on purpose.
        ColumnLayout {
            Layout.fillWidth: true
            Layout.fillHeight: true
            spacing: dp(4)

            RowLayout {
                Layout.fillWidth: true
                spacing: dp(6)
                ViewBar {
                    id: viewBar
                    Layout.fillWidth: true
                fontPt: root.fontPt
                uiScale: root.uiScale
                epochCount: root.epochCount
                onViewChanged: root.redraw()
                onEpochRequested: function (i) {
                    root.statusChanged(Julia.shell_select_epoch(i))
                    root.redraw()
                }
                    onStatusChanged: function (s) { root.statusChanged(s) }
                }
                // The view is REBUILT offscreen to be saved, never grabbed off the screen:
                // QMLMakie renders into Qt's FBO, so reading the live framebuffer hands back
                // noise. See src/gui/snapshot.jl.
                Button {
                    text: "Save view…"
                    Layout.alignment: Qt.AlignTop
                    font.pointSize: root.fontPt - 1
                    ToolTip.text: "write this plot to a PNG beside the working directory"
                    ToolTip.visible: hovered
                    onClicked: root.saveView()
                }
            }

            // The posterior's own controls: which parameter the marginal is for, and which
            // pair the scatter shows. Shown only on that view — they mean nothing on a map.
            RowLayout {
                Layout.fillWidth: true
                spacing: dp(6)
                visible: viewBar.viewIndex === 3
                Label { text: "marginal"; color: "#7f8c98"; font.pointSize: root.fontPt - 1 }
                ComboBox {
                    id: postX
                    Layout.preferredWidth: dp(150)
                    model: fitParamModel
                    textRole: "pname"
                    font.pointSize: root.fontPt - 1
                    onActivated: { root.applyPair(); }
                }
                Label { text: "against"; color: "#7f8c98"; font.pointSize: root.fontPt - 1 }
                ComboBox {
                    id: postY
                    Layout.preferredWidth: dp(150)
                    model: fitParamModel
                    textRole: "pname"
                    font.pointSize: root.fontPt - 1
                    onActivated: { root.applyPair(); }
                }
                Item { Layout.fillWidth: true }
                Label {
                    text: "red: the median · band: 16-84%"
                    color: "#7f8c98"
                    font.pointSize: root.fontPt - 1
                }
            }

            // One MakieArea per view, only one of them ever laid out. A StackLayout rather
            // than toggling `visible`: an area with zero size still holds its scene, so
            // switching back is a repaint rather than a rebuild.
            StackLayout {
                id: viewStack
                Layout.fillWidth: true
                Layout.fillHeight: true
                // The stack order follows the button order: orthographic, Mollweide, 3-D,
                // posterior.
                currentIndex: viewBar.viewIndex
                MakieArea { id: mSkyArea;  scene: modelSkyPlot }
                MakieArea { id: mollArea;  scene: mollPlot }
                MakieArea { id: starArea;  scene: starPlot }
                MakieArea { id: postArea;  scene: postPlot }
            }
        }
    }

    // Only the visible one: an update() on a zero-sized area costs a render pass for nothing,
    // and this is called on every keystroke in the form.
    // The fits table and the parameter selectors that read off it. Rebuilt whole rather than
    // patched: a fit is appended, never edited, so there is nothing to preserve.
    function refreshFits() {
        var keepX = postX.currentIndex, keepY = postY.currentIndex
        fitModel.clear()
        var rows = Julia.shell_fits()
        if (rows.length > 0) {
            var ls = rows.split("\n")
            for (var i = 0; i < ls.length; ++i) {
                var f = ls[i].split("\t")
                if (f.length < 9) continue
                fitModel.append({ fname: f[0], fmethod: f[1], fmodel: f[2], ftype: f[3],
                                  fchi2r: f[4], flogz: f[5], flogzerr: f[6], fns: f[7],
                                  fdiag: f[8] })
            }
        }
        var cur = parseInt(Julia.shell_current_fit())
        if (cur >= 1 && cur <= fitModel.count) fitList.currentIndex = cur - 1

        fitParamModel.clear()
        var prows = Julia.shell_fit_params()
        if (prows.length > 0) {
            var pl = prows.split("\n")
            for (var j = 0; j < pl.length; ++j) {
                var g = pl[j].split("\t")
                if (g.length < 5) continue
                fitParamModel.append({ pname: g[0], pval: g[1], perr: g[2],
                                       pq16: g[3], pq84: g[4] })
            }
        }
        var pair = Julia.shell_posterior_pair().split("\t")
        if (pair.length >= 2) {
            var xi = parseInt(pair[0]) - 1, yi = parseInt(pair[1]) - 1
            postX.currentIndex = (xi >= 0 && xi < fitParamModel.count) ? xi
                               : (keepX >= 0 && keepX < fitParamModel.count ? keepX : 0)
            postY.currentIndex = (yi >= 0 && yi < fitParamModel.count) ? yi
                               : (keepY >= 0 && keepY < fitParamModel.count ? keepY : 0)
        }
    }

    function saveView() {
        var w = root.currentPlotName()
        var m = Julia.shell_save_figure(w, "rotir_" + w + ".png",
                                        viewStack.width, viewStack.height)
        root.statusChanged(m.length > 0 ? m : "saved rotir_" + w + ".png")
    }

    // Which plot area `shell_save_figure` should rebuild. The name follows the VIEW, not the
    // tab: a tab showing four views would otherwise have to guess which was meant.
    function currentPlotName() {
        return ["sky", "mollweide", "star3d", "posterior"][viewBar.viewIndex]
    }

    function applyPair() {
        root.statusChanged(Julia.shell_set_posterior_pair(postX.currentIndex + 1,
                                                          postY.currentIndex + 1))
        postArea.update()
    }

    function redraw() {
        if (viewBar.viewIndex === 3) { postArea.update(); return }
        if (viewBar.viewIndex === 0) mSkyArea.update()
        else if (viewBar.viewIndex === 1) mollArea.update()
        else starArea.update()
    }
}
