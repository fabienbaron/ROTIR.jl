// The console pane: what the session has done, and what a running engine is saying.
//
// Read-only and append-only, and LIGHT — the same treatment as the OITOOLS GUI's console. A
// dark pane under light panels was the one piece of this window that did not match anything
// else in it.
//
// Lines beginning "» " are replayable calls — the same text the command log exports — so the
// pane doubles as a preview of the script the session would produce, which is what keeps the
// two honest about each other.

import QtQuick
import QtQuick.Controls
import QtQuick.Layouts

Rectangle {
    id: root
    color: "#fcfcfc"
    border.color: "#c8ccd0"
    border.width: 1
    radius: 3

    property string text: ""
    property real fontPt: 9
    // Follow the tail unless the user has scrolled up. Snapping back to the bottom while
    // somebody is reading an error further up is the single most annoying thing a log pane can
    // do, and a long optimiser trace makes it happen five times a second.
    property bool follow: true

    onTextChanged: if (follow) Qt.callLater(function () { view.positionViewAtEnd() })

    ColumnLayout {
        anchors.fill: parent
        anchors.margins: 4
        spacing: 2

        RowLayout {
            Layout.fillWidth: true
            spacing: 6
            Label {
                text: "console"
                color: "#666"
                font.pointSize: root.fontPt
                Layout.fillWidth: true
            }
            CheckBox {
                text: "follow"
                checked: root.follow
                font.pointSize: root.fontPt
                onToggled: { root.follow = checked; if (checked) view.positionViewAtEnd() }
            }
        }

        ListView {
            id: view
            Layout.fillWidth: true
            Layout.fillHeight: true
            clip: true
            model: root.text.length > 0 ? root.text.split("\n") : []
            // A ListView rather than a TextArea: a TextArea re-lays out the whole document on
            // every append, and the trace this pane carries is tens of thousands of lines.
            delegate: Text {
                width: view.width
                text: modelData
                wrapMode: Text.NoWrap
                font.family: "monospace"
                font.pointSize: root.fontPt
                // Replayable calls in near-black, narration in grey — the same dark-on-light
                // pairing the OITOOLS console uses. The point of the colour is to make the
                // script stand out from the commentary around it.
                color: modelData.indexOf("» ") === 0 ? "#111"
                     : modelData.indexOf("could not") >= 0 ||
                       modelData.indexOf("failed") >= 0 ? "#c62828"
                     : "#555"
            }
            ScrollBar.vertical: ScrollBar {}
            ScrollBar.horizontal: ScrollBar {}
            // Turning follow back on when the user returns to the bottom by hand, so it is not
            // a switch they have to remember to flip back.
            onMovementEnded: if (atYEnd) root.follow = true
        }
    }
}
