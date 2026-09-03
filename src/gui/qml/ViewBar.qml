// The view controls, shared by every tab that shows a surface.
//
// One selector — orthographic / Mollweide / 3-D — and never more than one view on screen.
// Showing all three at once was the first arrangement and it was wrong twice over: each panel
// got a third of the height, and the two the user was not looking at were still being drawn.
//
// The EPOCH SLIDER belongs to the two views that have an epoch. A Mollweide is the whole
// surface at once and has none, so the slider is hidden there rather than left inert — an
// enabled control that does nothing is a worse answer than an absent one.

import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml

// The same reasoning as the tabs: this is a Control so that the font set in the settings
// panel reaches the ticks and buttons inside it.
Pane {
    id: root
    padding: 0
    background: null
    font.pointSize: root.fontPt
    property real fontPt: 10
    property real uiScale: 1
    // 0 = orthographic (the sky projection), 1 = Mollweide, 2 = Full 3-D, 3 = posterior.
    //
    // The first three are a progression away from the data: the orthographic view is what the
    // interferometer measures at one epoch, the Mollweide adds the far side as a flat map, and
    // the 3-D scene is the whole body with the mouse in it. Putting 3-D in the middle made the
    // two flat projections non-adjacent, which is the comparison that gets made most.
    //
    // The posterior is last because it is about the FIT rather than the star: it appears only
    // once a sampler has run, and it has no epoch, no colormap and no decorations.
    property int viewIndex: 0
    property var colormaps: []
    property string colormap: "gist_heat"
    property int epochCount: 0
    property int epoch: 1
    property bool showEpoch: true
    property bool showColormap: true
    // Temperature or emergent intensity. Two different physical quantities: the map IS a
    // temperature, but what the interferometer measures is the intensity, which carries the
    // limb darkening and — through Planck — a strongly non-linear response to temperature.
    property bool showField: true
    // The posterior belongs to a FIT, so it appears where fits are run and nowhere else. The
    // Imaging tab shares every other view — a reconstruction is a map like any other and is
    // looked at the same three ways — but it has no posterior to show.
    property bool showPosterior: true
    property bool intensity: true
    property string intensityModel: "linear"
    property string band: "1.6500"
    property var decorations: []          // {name, label, on} per tick
    property bool graticulesOn: false
    property var graticuleColors: []
    signal viewChanged()
    signal epochRequested(int i)
    signal statusChanged(string s)

    function dp(px) { return Math.round(px * uiScale) }
    function pushField() {
        statusChanged(Julia.shell_set_surface_field(
            intensity ? "1" : "0", modelBox.currentText, bandField.text))
        viewChanged()
    }

    function loadField() {
        var f = Julia.shell_surface_field()
        if (f.length === 0) return
        var g = f.split("\t")
        if (g.length < 3) return
        intensity = g[0] === "1"
        intensityModel = g[1]
        band = g[2]
    }

    function pushGraticule() {
        statusChanged(Julia.shell_set_graticule(gratLat.value, gratLon.value,
                                                gratColor.currentText))
        viewChanged()
    }

    function loadDecorations() {
        if (decorations.length > 0) return
        var txt = Julia.shell_decorations()
        if (txt.length === 0) return
        var out = []
        var ls = txt.split("\n")
        for (var i = 0; i < ls.length; ++i) {
            var f = ls[i].split("\t")
            if (f.length < 3) continue
            out.push({ name: f[0], label: f[1], on: f[2] === "1" })
        }
        decorations = out
        for (var j = 0; j < out.length; ++j)
            if (out[j].name === "graticules") graticulesOn = out[j].on
        if (graticuleColors.length === 0) {
            var cs = Julia.shell_graticule_colors()
            graticuleColors = cs.length > 0 ? cs.split("\n") : ["black"]
        }
        var g = Julia.shell_graticule().split("\t")
        if (g.length >= 3) {
            gratLat.value = Math.round(parseFloat(g[0]))
            gratLon.value = Math.round(parseFloat(g[1]))
            var ci = graticuleColors.indexOf(g[2])
            if (ci >= 0) gratColor.currentIndex = ci
        }
    }

    function loadColormaps() {
        if (colormaps.length > 0) return
        var cm = Julia.shell_colormaps()
        colormaps = cm.length > 0 ? cm.split("\n") : []
    }

    contentItem: ColumnLayout {
    spacing: root.dp(3)

    RowLayout {
        Layout.fillWidth: true
        spacing: root.dp(6)
        Repeater {
            model: root.showPosterior ? ["orthographic", "Mollweide", "Full 3D", "posterior"]
                                      : ["orthographic", "Mollweide", "Full 3D"]
            Button {
                text: modelData
                checkable: true
                checked: root.viewIndex === index
                font.pointSize: root.fontPt
                onClicked: { root.viewIndex = index; root.viewChanged() }
            }
        }
        Item { Layout.fillWidth: true }
        Label {
            text: root.viewIndex === 2 ? "drag to rotate · scroll to zoom · axes in mas"
                : root.viewIndex === 0 ? "scroll to zoom · right-click resets"
                                       : "the whole surface, equal-area"
            color: "#7f8c98"
            font.pointSize: root.fontPt - 1
        }
    }

    // The epoch, where an epoch means something. ARROWS, not a slider: stepping to the next
    // night is the action, and a slider makes that a drag to a target you have to aim at —
    // with six epochs, four pixels wide each. The count beside them says where you are, which
    // is the only thing the slider was contributing.
    RowLayout {
        Layout.fillWidth: true
        spacing: dp(4)
        // Hidden on the Mollweide, which is the whole surface at once and has no epoch. The
        // orthographic and 3-D views both do.
        visible: root.showEpoch && root.viewIndex !== 1 && root.viewIndex !== 3
                 && root.epochCount > 1
        Label { text: "epoch"; color: "#7f8c98"; font.pointSize: root.fontPt - 1 }
        Button {
            text: "◀"
            enabled: root.epoch > 1
            font.pointSize: root.fontPt - 1
            onClicked: { root.epoch = root.epoch - 1; root.epochRequested(root.epoch) }
        }
        Label {
            text: root.epoch + " / " + root.epochCount
            color: "#555"
            font.pointSize: root.fontPt - 1
            horizontalAlignment: Text.AlignHCenter
            Layout.preferredWidth: root.dp(56)
        }
        Button {
            text: "▶"
            enabled: root.epoch < root.epochCount
            font.pointSize: root.fontPt - 1
            onClicked: { root.epoch = root.epoch + 1; root.epochRequested(root.epoch) }
        }
        Item { Layout.fillWidth: true }
    }

    RowLayout {
        Layout.fillWidth: true
        spacing: dp(6)
        visible: root.showField && root.viewIndex !== 3
        CheckBox {
            id: intensityBox
            text: "intensity"
            checked: root.intensity
            font.pointSize: root.fontPt - 1
            ToolTip.text: "limb darkening multiplied in; unticked shows the temperature map"
            ToolTip.visible: hovered
            onToggled: { root.intensity = checked; root.pushField() }
        }
        ComboBox {
            id: modelBox
            visible: root.intensity
            Layout.preferredWidth: root.dp(110)
            model: ["linear", "planck"]
            font.pointSize: root.fontPt - 1
            ToolTip.text: "linear: I proportional to T. planck: a real surface brightness at " +
                          "the observing wavelength"
            ToolTip.visible: hovered
            onActivated: { root.intensityModel = currentText; root.pushField() }
        }
        Label {
            visible: root.intensity && modelBox.currentText === "planck"
            text: "λ (µm)"
            color: "#7f8c98"
            font.pointSize: root.fontPt - 1
        }
        TextField {
            id: bandField
            visible: root.intensity && modelBox.currentText === "planck"
            Layout.preferredWidth: root.dp(66)
            // Filled from Julia with the wavelength that WOULD be used — the data's own mean
            // — rather than left at a sentinel the reader has to know the meaning of.
            text: root.band
            font.pointSize: root.fontPt - 1
            selectByMouse: true
            ToolTip.text: "defaults to the mean wavelength of the loaded observations"
            ToolTip.visible: hovered
            onEditingFinished: root.pushField()
        }
        // Decorations, to the RIGHT of the intensity controls: they annotate the same
        // picture, they are all one-click toggles, and putting them on a row of their own
        // would push the plot down for six tick boxes.
        //
        // `graticules` is LAST in the list Julia returns, because it is the only one with
        // settings of its own — ticking it opens a spacing and a colour to its right, and a
        // control that grows the row should grow it at the end.
        Repeater {
            id: decorRepeater
            model: root.decorations
            CheckBox {
                text: modelData.label
                checked: modelData.on
                font.pointSize: root.fontPt - 1
                onToggled: {
                    if (modelData.name === "graticules") root.graticulesOn = checked
                    root.statusChanged(
                        Julia.shell_set_decoration(modelData.name, checked ? "1" : "0"))
                    root.viewChanged()
                }
            }
        }
        // Spacing in DEGREES, which is how a map is read, and a colour that has to work
        // against both ends of the colormap — hence dark and light both on offer.
        Label {
            visible: root.graticulesOn
            text: "Δlat"; color: "#7f8c98"; font.pointSize: root.fontPt - 1
        }
        SpinBox {
            id: gratLat
            visible: root.graticulesOn
            Layout.preferredWidth: root.dp(78)
            from: 5; to: 90; stepSize: 5; value: 30
            font.pointSize: root.fontPt - 2
            ToolTip.text: "degrees between parallels"
            ToolTip.visible: hovered
            onValueModified: root.pushGraticule()
        }
        Label {
            visible: root.graticulesOn
            text: "Δlon"; color: "#7f8c98"; font.pointSize: root.fontPt - 1
        }
        SpinBox {
            id: gratLon
            visible: root.graticulesOn
            Layout.preferredWidth: root.dp(78)
            from: 5; to: 90; stepSize: 5; value: 30
            font.pointSize: root.fontPt - 2
            ToolTip.text: "degrees between meridians"
            ToolTip.visible: hovered
            onValueModified: root.pushGraticule()
        }
        ComboBox {
            id: gratColor
            visible: root.graticulesOn
            Layout.preferredWidth: root.dp(96)
            model: root.graticuleColors
            font.pointSize: root.fontPt - 2
            ToolTip.text: "a light graticule vanishes on a pale map and a dark one on a dark map"
            ToolTip.visible: hovered
            onActivated: root.pushGraticule()
        }
        Item { Layout.fillWidth: true }
    }

    RowLayout {
        Layout.fillWidth: true
        spacing: dp(3)
        visible: root.showColormap && root.viewIndex !== 3
        Label { text: "colormap"; color: "#7f8c98"; font.pointSize: root.fontPt - 1 }
        Repeater {
            model: root.colormaps
            Button {
                text: modelData
                checkable: true
                checked: root.colormap === modelData
                font.pointSize: root.fontPt - 2
                onClicked: {
                    root.colormap = modelData
                    root.statusChanged(Julia.shell_set_colormap(modelData))
                    root.viewChanged()
                }
            }
        }
        Item { Layout.fillWidth: true }
    }
    }
}
