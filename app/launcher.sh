#!/bin/sh
#
# Start the application, having settled the windowing system first.
#
# THIS IS NOT OPTIONAL ON LINUX. A bundle is a sysimage, and every package `__init__` in a
# sysimage runs before any user code: GLFW has called `glfwInit` and QML has constructed its
# `QGuiApplication` before `julia_main` gets control. Neither choice can be revised from
# inside the process, so the only moment left is here, before it starts.
#
# Left alone on a Wayland session, GLFW takes its hardcoded X11 and Qt follows the session to
# Wayland. That is two windowing systems and two EGL display connections in one process, and
# it is the configuration OITOOLS' `prefer_native_wayland!` exists to avoid — the call ROTIR
# makes in bin/rotirgui.jl, and which is too late to make from inside a bundle.
#
# Anything set by hand is honoured -- this only fills in a blank, so
#
#     QT_QPA_PLATFORM=xcb ./ROTIR            # both halves on X11 instead
#
# still does what it says.
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd -P)

if [ -n "${WAYLAND_DISPLAY:-}" ] && [ -z "${JULIA_GLFW_PLATFORM:-}" ] && [ -z "${QT_QPA_PLATFORM:-}" ]; then
    JULIA_GLFW_PLATFORM=wayland
    export JULIA_GLFW_PLATFORM
fi

# THE X11 LOCALE DIRECTORY, for the bundled libxkbcommon.
#
# `xkbcommon_jll` was compiled by BinaryBuilder, so its default locale directory is the
# sandbox path it was built in -- `/workspace/destdir/share/X11/locale`, which exists on no
# real machine. It therefore never finds the Compose file for the user's locale and every
# launch prints
#
#     xkbcommon: ERROR: [XKB-679] No Compose file for locale "en_US.UTF-8"
#     GLFW.GLFWError(GLFW.PLATFORM_ERROR, "Wayland: Failed to create XKB compose table")
#
# and dead keys stop working. Installing the host's X11 data does NOT help, because the
# baked-in path is not where it looks. Same shape as the QML import path QMLMakie bakes at
# precompile time -- see register_qml_modules! in app/src/ROTIRApp.jl.
#
# Only filled in when the host actually has the data, and never overriding a user's own value.
if [ -z "${XLOCALEDIR:-}" ]; then
    # The bundle's own copy FIRST -- that is what makes this self-contained. The host paths
    # remain as a fallback for a bundle built where there was nothing to stage.
    for d in "$here/share/X11/locale" /usr/share/X11/locale /usr/local/share/X11/locale; do
        if [ -f "$d/compose.dir" ]; then
            XLOCALEDIR=$d
            export XLOCALEDIR
            break
        fi
    done
fi

exec "$here/bin/ROTIRApp" "$@"
