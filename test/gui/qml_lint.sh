#!/usr/bin/env bash
# Static check of every QML file with Qt's own linter.
#
#     test/gui/qml_lint.sh
#
# WHICH BINARY, because this is a trap. `/usr/bin/qmllint` on this machine is a Qt5-era
# wrapper that exits 0 on ANYTHING — it passes a file containing a bare syntax error
# (`var x = ;`) without a word. Only the versioned Qt6 binary actually parses. So this script
# refuses to run rather than report a clean bill from a linter that checks nothing: a silent
# pass is worse than no check, and that is exactly what it looked like the first time.
#
# WHAT IS FILTERED, and why each is not a defect:
#
#   * `jlqml` and `Makie` fail to import. Both modules are registered by QML.jl at RUN time
#     through `@qmlfunction`/`qmlRegisterType`; neither has a qmldir on disk for a static tool
#     to find. `MakieArea was not found` follows from the second.
#   * `Unqualified access`. QML resolves ids and root-scope functions lexically, and this file
#     set uses that deliberately: `dp(...)` is a function on each tab's root, `modelData`,
#     `idx`, `mjd` and friends are delegate model roles, and `Julia.foo` is the registered
#     singleton. The linter cannot see any of them and flags all of them.
#
# Anything else is a real finding and fails the run.
set -uo pipefail
cd "$(dirname "$0")/../.."
QMLDIR="src/gui/qml"

LINT=""
for cand in /usr/lib/qt6/bin/qmllint /usr/lib/qt6/qmllint "$(command -v qmllint6 2>/dev/null)"; do
    [ -n "$cand" ] && [ -x "$cand" ] || continue
    # The version string is the test: 6.x parses, the Qt5 wrapper reports "qmllint 1.0".
    case "$("$cand" --version 2>&1)" in
        *"qmllint 6"*) LINT="$cand"; break ;;
    esac
done
if [ -z "$LINT" ]; then
    echo "SKIP  no Qt6 qmllint found (/usr/bin/qmllint is a stub that passes everything)"
    exit 77
fi
echo "linter: $LINT ($("$LINT" --version 2>&1))"

FAIL=0
for f in "$QMLDIR"/*.qml; do
    out=$("$LINT" -I "$QMLDIR" "$f" 2>&1 |
          grep -E "^(Warning|Error):" |
          grep -vE 'module "(jlqml|Makie)"' |
          grep -vE 'Failed to import (jlqml|Makie)' |
          grep -vE 'MakieArea was not found' |
          grep -vE 'Unqualified access')
    if [ -n "$out" ]; then
        echo "FAIL  $(basename "$f")"
        echo "$out" | sed 's/^/        /'
        FAIL=1
    else
        echo "  ok  $(basename "$f")"
    fi
done

# SIGNAL ARITY, which qmllint does not check at all. A signal declared with two parameters and
# emitted with one passes every check above while the QML engine reports "Insufficient
# arguments" at run time and the handler sees `undefined` — see the script's own header.
if python3 "$(dirname "$0")/qml_signal_arity.py" "$QMLDIR"; then :; else FAIL=1; fi

# The filter must not be able to hide a real problem: a file with a deliberate syntax error has
# to fail. Without this the whole script degrades to "prints ok" the day the linter changes its
# message format — which is the failure mode it exists to prevent.
PROBE="$(mktemp -d)/Probe.qml"
printf 'import QtQuick\nItem { Component.onCompleted: { var x = ; } }\n' > "$PROBE"
# Captured first, not piped: the linter EXITS NON-ZERO when it finds something, and under
# `pipefail` that exit status is the pipeline's — so `linter | grep -q` was false whenever the
# grep would have matched, and the self-check reported the opposite of the truth.
probe_out="$("$LINT" "$PROBE" 2>&1)"
if echo "$probe_out" | grep -qE "^(Warning|Error):"; then
    echo "  ok  self-check: the linter still rejects a syntax error"
else
    echo "FAIL  self-check: the linter accepted a syntax error — it is not analysing"
    FAIL=1
fi
rm -rf "$(dirname "$PROBE")"

exit $FAIL
