#!/bin/sh
#
# Register this bundle with the desktop: menu entry, icon, and .oifits association.
#
# Run once after unpacking, from anywhere:   ./install-desktop.sh
#
# The paths are computed from where THIS script sits, so the bundle can be unpacked anywhere
# and moved afterwards by re-running it. Everything is written under ~/.local/share, so no
# administrator rights are needed and an uninstall is three `rm`s (printed at the end).
set -eu
here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd -P)
apps="$HOME/.local/share/applications"
icons="$HOME/.local/share/icons/hicolor/256x256/apps"
mime="$HOME/.local/share/mime/packages"
mkdir -p "$apps" "$icons" "$mime"

cp "$here/share/rotir/app/assets/rotir.png" "$icons/rotir.png"

# A MIME type for OIFITS: it is FITS underneath, so without this the desktop hands .oifits to
# whatever claims application/fits and "Open with" never offers this application. Declaring the
# same type as OITOOLS is deliberate — both read the same files, and the desktop then offers
# both, which is the behaviour a user wants when they have both installed.
cat > "$mime/oifits.xml" <<XML
<?xml version="1.0" encoding="UTF-8"?>
<mime-info xmlns="http://www.freedesktop.org/standards/shared-mime-info">
  <mime-type type="application/x-oifits">
    <comment>Optical interferometry data (OIFITS)</comment>
    <sub-class-of type="application/fits"/>
    <glob pattern="*.oifits"/>
  </mime-type>
</mime-info>
XML

# Terminal=false matters: with it true the desktop spawns a terminal that serves no purpose,
# since stdout goes nowhere useful and failures are written to the crash log instead.
cat > "$apps/rotir.desktop" <<DESKTOP
[Desktop Entry]
Type=Application
Name=ROTIR
GenericName=Stellar Surface Imaging
Comment=Model and reconstruct stellar surfaces from optical interferometry
Exec=$here/ROTIR %f
Icon=rotir
Terminal=false
Categories=Science;Astronomy;Education;
MimeType=application/x-oifits;
StartupWMClass=ROTIR
DESKTOP

command -v update-desktop-database >/dev/null 2>&1 && update-desktop-database "$apps" || true
command -v update-mime-database    >/dev/null 2>&1 && update-mime-database "$HOME/.local/share/mime" || true
command -v gtk-update-icon-cache   >/dev/null 2>&1 && gtk-update-icon-cache -f -t "$HOME/.local/share/icons/hicolor" 2>/dev/null || true

echo "Installed. ROTIR is in the application menu, and .oifits files can be opened with it."
echo "To undo:  rm $apps/rotir.desktop $icons/rotir.png $mime/oifits.xml"
