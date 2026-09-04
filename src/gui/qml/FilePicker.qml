// The in-window file picker.
//
// Not `QtQuick.Dialogs.FileDialog` — see the note at the top of src/gui/filepicker.jl. Julia
// supplies the listing, the parent directory and the shortcut places; this file only draws
// them and reports what was chosen.
//
// A Popup inside our own window, so it shares the GL context rather than getting one of its
// own. That is the whole reason it exists.

import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml

Popup {
    id: root
    modal: true
    focus: true
    closePolicy: Popup.CloseOnEscape | Popup.CloseOnPressOutside

    property string folder: ""
    property string selected: ""
    // Every ticked file, newline-separated — the shape `shell_open_many` takes. Multi-select
    // is the normal case here: λ And is six nights in six files and they are one dataset.
    property var checked: []
    property bool showAll: false
    // Which button opened the picker; read back by Main.qml's onAccepted.
    // Whether "Add epoch" is meaningful: there has to be a dataset to add to.
    property bool canAdd: false
    property real fontPt: 10
    // `mode` says WHAT to do with the chosen files, which is the thing the two top-bar
    // buttons were trying to express with their names: "open" starts a dataset and cuts a
    // single file on timestamp gaps, "single" starts one and takes the file whole, "add"
    // appends to the dataset already loaded. Deciding it here, where the files are, is where
    // the decision belongs.
    signal accepted(string paths, string mode)
    // Set by the caller: an orbit file is chosen through the same picker but goes somewhere
    // else entirely, so it bypasses the three dataset modes.
    property string purpose: "data"

    function openAt(startFolder) {
        folder = Julia.picker_start(startFolder)
        selected = ""
        checked = []
        anchorRow = -1
        refresh()
        open()
        moveTo(0)
    }

    // Ticked files if any, otherwise the singly-selected one. That way a single click still
    // works exactly as before and the checkboxes are an addition rather than a new ritual.
    function chosen() {
        return checked.length > 0 ? checked.join("\n") : selected
    }

    // What Return does on a row, shared by the key handler and the double-click.
    function activate(i) {
        if (i < 0 || i >= entries.count) return
        var e = entries.get(i)
        if (e.kind === "dir") {
            folder = Julia.picker_join(folder, e.name)
            refresh()
            fileView.currentIndex = 0
        } else {
            selected = Julia.picker_join(folder, e.name)
            // BOTH arguments: `accepted(string paths, string mode)`. This passed only the
            // paths, so a double-click or Return on a file reached the handler with `mode`
            // undefined — "Insufficient arguments" from the QML engine, and then a silent
            // fall through every branch of the mode dispatch, so the file simply did not open.
            //
            // The mode is the PRIMARY button's: "open" for a dataset, and the purpose itself
            // for an orbit or a map, which is what those pickers' single Open button sends.
            accepted(chosen(), root.purpose === "data" ? "open" : root.purpose)
            close()
        }
    }

    // Anchor for a Shift-range, in ROW terms — a range is between two positions in the
    // listing, not between two paths, and the listing can change under it.
    property int anchorRow: -1

    // The classical selection gestures, on top of the tick boxes rather than instead of them:
    //   plain click  -> select just this one
    //   Ctrl-click   -> add or remove this one
    //   Shift-click  -> everything between the anchor and this one
    // The ticks stay because they survive a change of directory and a Shift-click does not,
    // which is how a selection is built across folders.
    function selectRow(i, mods) {
        if (i < 0 || i >= entries.count) return
        fileView.currentIndex = i
        var e = entries.get(i)
        if (e.kind !== "file") return
        var pth = Julia.picker_join(folder, e.name)
        if (mods & Qt.ShiftModifier && anchorRow >= 0) {
            var lo = Math.min(anchorRow, i), hi = Math.max(anchorRow, i)
            var c = []
            for (var k = lo; k <= hi; ++k) {
                var ek = entries.get(k)
                if (ek.kind === "file") c.push(Julia.picker_join(folder, ek.name))
            }
            checked = c
            selected = pth
        } else if (mods & Qt.ControlModifier) {
            toggle(pth, checked.indexOf(pth) < 0)
            selected = pth
            anchorRow = i
        } else {
            checked = []
            selected = pth
            anchorRow = i
        }
    }

    // Move the cursor and select what it lands on, without disturbing a tick set.
    function moveTo(i) {
        if (entries.count === 0) return
        var j = Math.max(0, Math.min(i, entries.count - 1))
        fileView.currentIndex = j
        var e = entries.get(j)
        selected = e.kind === "file" ? Julia.picker_join(folder, e.name) : ""
        anchorRow = j
    }

    function toggle(path, on) {
        var c = checked.slice()
        var i = c.indexOf(path)
        if (on && i < 0) c.push(path)
        else if (!on && i >= 0) c.splice(i, 1)
        checked = c
    }

    function refresh() {
        var rows = Julia.picker_list(folder, showAll ? "1" : "0", root.purpose)
        entries.clear()
        if (rows.length === 0) return
        var lines = rows.split("\n")
        for (var i = 0; i < lines.length; ++i) {
            var f = lines[i].split("\t")
            if (f.length < 3) continue
            entries.append({ kind: f[0], name: f[1], size: f[2] })
        }
    }

    ListModel { id: entries }
    ListModel { id: places }


    Component.onCompleted: {
        var rows = Julia.picker_places()
        if (rows.length === 0) return
        var lines = rows.split("\n")
        for (var i = 0; i < lines.length; ++i) {
            var f = lines[i].split("\t")
            if (f.length < 2) continue
            places.append({ label: f[0], path: f[1] })
        }
    }

    // Keyboard navigation, as SHORTCUTS rather than `Keys` handlers.
    //
    // `Keys` needs the item to have active focus, and in this popup it never does: the path
    // TextField takes focus when the popup opens and swallows Return through its own
    // `onAccepted`. Chasing the focus chain was tried twice — handlers on the ListView, then
    // on a FocusScope — and both left the row highlighted with Return doing nothing, which
    // reads as a picker that has failed. A Shortcut is window-level and fires whatever has
    // focus, which is what "press Return to open the selected row" should mean.
    //
    // `enabled: root.opened` scopes them: they must not fire while the picker is closed and
    // the main window is being typed into.
    Shortcut {
        sequences: ["Return", "Enter"]
        enabled: root.opened
        onActivated: root.activate(fileView.currentIndex)
    }
    Shortcut {
        sequence: "Down"
        enabled: root.opened
        onActivated: root.moveTo(fileView.currentIndex + 1)
    }
    Shortcut {
        sequence: "Up"
        enabled: root.opened
        onActivated: root.moveTo(fileView.currentIndex - 1)
    }
    Shortcut {
        sequence: "Space"
        enabled: root.opened
        onActivated: {
            var e = entries.get(fileView.currentIndex)
            if (e && e.kind === "file") {
                var pth = Julia.picker_join(root.folder, e.name)
                root.toggle(pth, root.checked.indexOf(pth) < 0)
            }
        }
    }

    ColumnLayout {
        anchors.fill: parent
        spacing: 6

        RowLayout {
            Layout.fillWidth: true
            spacing: 6
            Button {
                text: "↑"
                font.pointSize: root.fontPt
                ToolTip.text: "parent directory"
                ToolTip.visible: hovered
                onClicked: { root.folder = Julia.picker_parent(root.folder); root.refresh() }
            }
            // Editable, because typing or pasting a path is faster than navigating to it and
            // is the only way to reach somewhere the shortcuts do not cover.
            TextField {
                id: pathField
                Layout.fillWidth: true
                text: root.folder
                font.pointSize: root.fontPt
                selectByMouse: true
                onAccepted: { root.folder = text; root.refresh() }
            }
            Label {
                text: "click · ctrl-click · shift-click · or tick"
                color: "#7f8c98"
                font.pointSize: root.fontPt - 1
            }
            CheckBox {
                text: "all files"
                checked: root.showAll
                font.pointSize: root.fontPt
                onToggled: { root.showAll = checked; root.refresh() }
            }
        }

        RowLayout {
            Layout.fillWidth: true
            Layout.fillHeight: true
            spacing: 6

            // Shortcuts on the left, listing on the right.
            Frame {
                Layout.preferredWidth: 170
                Layout.fillHeight: true
                ListView {
                    id: placeView
                    anchors.fill: parent
                    clip: true
                    model: places
                    delegate: ItemDelegate {
                        width: placeView.width
                        text: label
                        font.pointSize: root.fontPt
                        onClicked: { root.folder = path; root.refresh() }
                    }
                }
            }

            Frame {
                Layout.fillWidth: true
                Layout.fillHeight: true
                ListView {
                    id: fileView
                    anchors.fill: parent
                    clip: true
                    model: entries
                    currentIndex: 0
                    highlightMoveDuration: 0
                    ScrollBar.vertical: ScrollBar {}
                    delegate: ItemDelegate {
                        width: fileView.width
                        // Highlight the row the keyboard is on OR the one clicked, so the
                        // current position is visible however it was reached.
                        // Highlight anything in the SELECTION as well as the keyboard row, so
                        // a Shift-range reads as one block rather than as one highlighted row
                        // with some ticks scattered around it.
                        highlighted: ListView.isCurrentItem ||
                                     (kind === "file" &&
                                      (root.selected === Julia.picker_join(root.folder, name) ||
                                       root.checked.indexOf(
                                           Julia.picker_join(root.folder, name)) >= 0))
                        contentItem: RowLayout {
                            spacing: 8
                            CheckBox {
                                visible: kind === "file"
                                checked: root.checked.indexOf(
                                    Julia.picker_join(root.folder, name)) >= 0
                                onToggled: root.toggle(Julia.picker_join(root.folder, name),
                                                       checked)
                            }
                            Label {
                                text: kind === "dir" ? "🗀" : "🗎"
                                font.pointSize: root.fontPt
                            }
                            Label {
                                Layout.fillWidth: true
                                text: name
                                elide: Text.ElideMiddle
                                font.pointSize: root.fontPt
                            }
                            Label {
                                text: size.length > 0 ? (Math.round(size / 1024) + " kB") : ""
                                color: "#7f8c98"
                                font.pointSize: root.fontPt - 1
                            }
                        }
                        // A TapHandler, because an ItemDelegate's `onClicked` does not carry
                        // the keyboard modifiers, and Shift/Ctrl-click is the whole point.
                        //
                        // `point.modifiers`, NOT the signal argument's. The `eventPoint` a
                        // TapHandler hands its signal is a `QEventPoint`, which has no
                        // `modifiers` at all — so `evt.modifiers` was `undefined`, and
                        // `undefined & Qt.ShiftModifier` is 0 in JavaScript, which meant every
                        // click took the plain-click branch and Shift/Ctrl-click silently did
                        // nothing from the day it was written. The handler's OWN `point` is a
                        // `QQuickHandlerPoint`, and that is the type carrying `modifiers`.
                        // Found by the Qt linter; nothing else would have, since the wrong
                        // expression is valid JavaScript that quietly evaluates to 0.
                        TapHandler {
                            id: rowTap
                            acceptedButtons: Qt.LeftButton
                            onSingleTapped: {
                                if (kind === "dir") {
                                    root.folder = Julia.picker_join(root.folder, name)
                                    root.refresh()
                                    root.anchorRow = -1
                                } else {
                                    root.selectRow(index, rowTap.point.modifiers)
                                }
                            }
                        }
                        onDoubleClicked: root.activate(index)
                    }
                }
            }
        }

        RowLayout {
            Layout.fillWidth: true
            spacing: 6
            Item { Layout.fillWidth: true }
            Button {
                text: "Tick all"
                font.pointSize: root.fontPt - 1
                onClicked: {
                    var c = []
                    for (var i = 0; i < entries.count; ++i) {
                        var e = entries.get(i)
                        if (e.kind === "file") c.push(Julia.picker_join(root.folder, e.name))
                    }
                    root.checked = c
                }
            }
            Button {
                text: "Clear"
                font.pointSize: root.fontPt - 1
                enabled: root.checked.length > 0
                onClicked: root.checked = []
            }
            Label {
                Layout.fillWidth: true
                text: root.checked.length > 1 ? root.checked.length + " files selected"
                                              : root.selected
                elide: Text.ElideMiddle
                color: "#7f8c98"
                font.pointSize: root.fontPt
            }
            Button {
                text: "Cancel"
                font.pointSize: root.fontPt
                onClicked: root.close()
            }
            // THREE openings, not one button plus a mode chosen somewhere else. A file is
            // either a new dataset, a new dataset that must not be cut into epochs, or another
            // epoch of the one already loaded — and which of those it is depends on the file,
            // so the choice belongs next to the file.
            Button {
                visible: root.purpose !== "data"
                text: "Open"
                enabled: root.selected.length > 0
                font.pointSize: root.fontPt
                onClicked: { root.accepted(root.chosen(), root.purpose); root.close() }
            }
            Button {
                visible: root.purpose === "data"
                text: "Open as single epoch"
                enabled: root.selected.length > 0 || root.checked.length > 0
                font.pointSize: root.fontPt - 1
                ToolTip.text: "take the file whole instead of cutting it on gaps in its V² " +
                              "timestamps"
                ToolTip.visible: hovered
                onClicked: { root.accepted(root.chosen(), "single"); root.close() }
            }
            Button {
                visible: root.purpose === "data" && root.canAdd
                text: root.checked.length > 1 ? "Add " + root.checked.length + " epochs"
                                              : "Add epoch"
                enabled: root.selected.length > 0 || root.checked.length > 0
                font.pointSize: root.fontPt - 1
                ToolTip.text: "append to the dataset already loaded"
                ToolTip.visible: hovered
                onClicked: { root.accepted(root.chosen(), "add"); root.close() }
            }
            Button {
                visible: root.purpose === "data"
                text: root.checked.length > 1 ? "Open " + root.checked.length : "Open"
                enabled: root.selected.length > 0 || root.checked.length > 0
                font.pointSize: root.fontPt
                onClicked: { root.accepted(root.chosen(), "open"); root.close() }
            }
        }
    }
}
