import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml

// Export the model, or a reconstruction, as SKY IMAGES for OITOOLS: an (nx, nx, channels,
// epochs) FITS cube in OITOOLS's own layout. Everything about what that cube is lives in
// src/sky_export.jl; this only gathers the choices into `shell_sky_export`'s key=value lines.
Popup {
    id: root
    modal: true
    focus: true
    closePolicy: Popup.CloseOnEscape
    property real fontPt: 10
    property real uiScale: 1.0
    signal statusChanged(string s)
    function dp(px) { return Math.round(px * root.uiScale) }

    property var sources: []
    property var epochs: []
    property var channels: []
    property string reason: ""
    property string intensity: ""
    property string error: ""
    property string autoPix: ""
    property string autoNx: ""

    // Fill from Julia each time it opens: the dataset, model or reconstructions may have
    // changed since the last time.
    function openDialog() {
        var t = Julia.shell_sky_export_defaults().split("\n")
        var kv = {}
        for (var i = 0; i < t.length; ++i) {
            var j = t[i].indexOf("\t")
            if (j > 0) kv[t[i].substring(0, j)] = t[i].substring(j + 1)
        }
        root.sources = (kv["sources"] || "model").split("|")
        root.epochs = kv["epochs"] ? kv["epochs"].split("|") : []
        root.channels = kv["channels"] ? kv["channels"].split("|") : []
        root.reason = kv["reason"] || ""
        root.intensity = kv["intensity"] || ""
        root.error = kv["error"] || ""
        root.autoPix = kv["pixsize"] || ""
        root.autoNx = kv["nx"] || ""
        pixField.text = root.autoPix
        nxField.text = root.autoNx
        epochMode.currentIndex = root.epochs.length > 0 ? 0 : 2
        chanMode.currentIndex = root.channels.length > 0 ? 0 : 2
        open()
    }

    function spec() {
        var em = ["data", "data:" + epochText.text, "mjd:" + epochText.text,
                  "phases:" + epochText.text][epochMode.currentIndex]
        var cm = ["data", "mean", "um:" + chanText.text][chanMode.currentIndex]
        return ["path=" + fileField.text,
                "source=" + sourceBox.currentText,
                "epochs=" + em,
                "channels=" + cm,
                "pixsize=" + pixField.text,
                "nx=" + nxField.text,
                "rendering=" + (renderBox.currentIndex === 0 ? "bandlimited" : "area"),
                "normalize=" + (normBox.currentIndex === 0 ? "plane" : "relative"),
                "check=" + (checkBox.checked ? "1" : "0")].join("\n")
    }

    contentItem: ColumnLayout {
        spacing: dp(6)
        Label {
            text: "Export sky images (FITS) for OITOOLS"
            font.bold: true
            font.pointSize: root.fontPt + 1
        }
        Label {
            visible: root.error.length > 0
            Layout.fillWidth: true
            wrapMode: Text.WordWrap
            color: "#c0392b"
            text: root.error
            font.pointSize: root.fontPt
        }
        GridLayout {
            columns: 3
            columnSpacing: dp(8)
            rowSpacing: dp(4)

            Label { text: "Source"; font.pointSize: root.fontPt }
            ComboBox {
                id: sourceBox
                model: root.sources
                Layout.preferredWidth: dp(200)
                font.pointSize: root.fontPt
            }
            Label {
                text: root.intensity
                color: "#7f8c98"
                font.pointSize: root.fontPt - 1
                ToolTip.text: "the model's \"Intensity from\", evaluated at each channel's wavelength"
                ToolTip.visible: maI.containsMouse
                MouseArea { id: maI; anchors.fill: parent; hoverEnabled: true }
            }

            Label { text: "Epochs"; font.pointSize: root.fontPt }
            ComboBox {
                id: epochMode
                model: ["all data epochs", "data epochs (list)", "MJD list", "rotation phases"]
                Layout.preferredWidth: dp(200)
                font.pointSize: root.fontPt
            }
            TextField {
                id: epochText
                enabled: epochMode.currentIndex > 0
                Layout.preferredWidth: dp(220)
                font.pointSize: root.fontPt
                placeholderText: ["", "1,3,5 (1-based)", "59000.1, 59001.2", "N phases"][epochMode.currentIndex]
            }

            Label { text: "Channels"; font.pointSize: root.fontPt }
            ComboBox {
                id: chanMode
                model: ["data channels (" + root.channels.length + ")", "mean wavelength", "custom (µm)"]
                Layout.preferredWidth: dp(200)
                font.pointSize: root.fontPt
            }
            TextField {
                id: chanText
                enabled: chanMode.currentIndex === 2
                Layout.preferredWidth: dp(220)
                font.pointSize: root.fontPt
                placeholderText: "1.55, 1.60, 1.65"
            }

            Label { text: "Pixel (mas)"; font.pointSize: root.fontPt }
            TextField {
                id: pixField
                Layout.preferredWidth: dp(200)
                font.pointSize: root.fontPt
                selectByMouse: true
            }
            Button {
                text: "auto"
                font.pointSize: root.fontPt - 1
                onClicked: { pixField.text = root.autoPix; nxField.text = root.autoNx }
                ToolTip.text: root.reason
                ToolTip.visible: hovered
            }

            Label { text: "Image size (px)"; font.pointSize: root.fontPt }
            TextField {
                id: nxField
                Layout.preferredWidth: dp(200)
                font.pointSize: root.fontPt
                selectByMouse: true
            }
            Label {
                text: root.reason
                color: "#7f8c98"
                elide: Text.ElideRight
                Layout.preferredWidth: dp(220)
                font.pointSize: root.fontPt - 2
            }

            Label { text: "Rendering"; font.pointSize: root.fontPt }
            ComboBox {
                id: renderBox
                model: ["band-limited", "area"]
                Layout.preferredWidth: dp(200)
                font.pointSize: root.fontPt
                ToolTip.text: "band-limited: OITOOLS reproduces ROTIR's visibilities (~1e-5), " +
                              "with ringing at the limb. area: positive, flux-conserving, for " +
                              "display; its visibilities are box-filtered"
                ToolTip.visible: hovered
            }
            Item { width: 1 }

            Label { text: "Normalization"; font.pointSize: root.fontPt }
            ComboBox {
                id: normBox
                model: ["each plane = 1", "relative flux"]
                Layout.preferredWidth: dp(200)
                font.pointSize: root.fontPt
                ToolTip.text: "each plane sums to 1 (OITOOLS imaging), or physical relative " +
                              "flux across channels and epochs (keeps the spectrum and light curve)"
                ToolTip.visible: hovered
            }
            Item { width: 1 }

            Label { text: "File"; font.pointSize: root.fontPt }
            TextField {
                id: fileField
                text: "rotir_sky.fits"
                Layout.preferredWidth: dp(200)
                font.pointSize: root.fontPt
                selectByMouse: true
            }
            CheckBox {
                id: checkBox
                checked: true
                text: "check against ROTIR"
                font.pointSize: root.fontPt
                ToolTip.text: "compute OITOOLS's visibilities from the cube on the data's uv " +
                              "and compare with ROTIR's own (data epochs only)"
                ToolTip.visible: hovered
            }
        }
        RowLayout {
            Layout.alignment: Qt.AlignRight
            Button {
                text: "Export"
                enabled: root.error.length === 0
                font.pointSize: root.fontPt
                onClicked: {
                    root.statusChanged(Julia.shell_sky_export(root.spec()))
                    root.close()
                }
            }
            Button {
                text: "Cancel"
                font.pointSize: root.fontPt
                onClicked: root.close()
            }
        }
    }
}
