// Which observables a fit or a reconstruction uses: V², T3amp and T3φ, each in or out.
//
// The ticks are weights of 1 or 0 on that observable's χ², as in OITOOLS. A common sequence
// is to fit the shape on V² and T3amp, whose meaning does not depend on the surface, and add
// the closure phases once the map is being imaged. `tab` picks which of the two sets Julia
// keeps — "model" for the Model tab's fits, "image" for the Imaging tab.

import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml

RowLayout {
    id: root
    property string tab: "model"
    property real fontPt: 10
    signal statusReported(string msg)
    spacing: 6

    function load() {
        var w = Julia.shell_obs_weights(root.tab).split("\t")
        v2Box.checked = w[0] === "1"
        t3ampBox.checked = w[1] === "1"
        t3phiBox.checked = w[2] === "1"
    }
    function push() {
        root.statusReported(Julia.shell_set_obs_weights(root.tab, v2Box.checked,
                                                        t3ampBox.checked, t3phiBox.checked))
    }
    Component.onCompleted: load()

    Label { text: "use"; font.pointSize: root.fontPt }
    CheckBox {
        id: v2Box
        text: "V²"
        checked: true
        font.pointSize: root.fontPt - 1
        onToggled: root.push()
    }
    CheckBox {
        id: t3ampBox
        text: "T3amp"
        checked: true
        font.pointSize: root.fontPt - 1
        onToggled: root.push()
    }
    CheckBox {
        id: t3phiBox
        text: "T3φ"
        checked: true
        font.pointSize: root.fontPt - 1
        ToolTip.text: "closure phases; untick to constrain only the symmetric part of the " +
                      "surface"
        ToolTip.visible: hovered
        onToggled: root.push()
    }
    Item { Layout.fillWidth: true }
}
