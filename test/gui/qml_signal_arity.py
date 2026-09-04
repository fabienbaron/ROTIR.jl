#!/usr/bin/env python3
"""Check that every QML signal is emitted with the arity it was declared with.

qmllint does NOT do this. `signal accepted(string paths, string mode)` called as
`accepted(chosen())` passed the linter on every run, while the QML engine reported
"Insufficient arguments" at runtime and the handler received `mode` as undefined — so a
double-click on a file fell through every branch of the mode dispatch and silently failed to
open it. A warning in a log nobody reads is how that survived.

    qml_signal_arity.py <qml-dir>

Exits non-zero, listing file:line, when a call's argument count differs from the declaration.
"""
import re
import sys
import pathlib


def top_level_commas(args: str) -> int:
    """Commas at depth zero: nested calls and string literals do not separate arguments."""
    flat = args
    while True:
        stripped = re.sub(r"\([^()]*\)", "", flat)
        if stripped == flat:
            break
        flat = stripped
    flat = re.sub(r"\[[^\[\]]*\]", "", flat)
    flat = re.sub(r'"[^"]*"|\'[^\']*\'', "", flat)
    return flat.count(",")


def main(qmldir: pathlib.Path) -> int:
    files = sorted(qmldir.glob("*.qml"))
    if not files:
        print(f"FAIL  no .qml files under {qmldir}")
        return 1

    declared = {}
    for f in files:
        for m in re.finditer(r"^\s*signal\s+(\w+)\s*\(([^)]*)\)", f.read_text(), re.M):
            params = [a for a in m.group(2).split(",") if a.strip()]
            declared[m.group(1)] = len(params)
    if not declared:
        print("FAIL  no signal declarations found — the check would pass vacuously")
        return 1

    bad = 0
    for f in files:
        src = f.read_text()
        for name, want in declared.items():
            for m in re.finditer(r"(?<![\w.])(?:root\.)?" + name + r"\s*\(", src):
                # The declaration itself, and property/function definitions, are not calls.
                if re.search(r"\b(signal|function)\s+\w*$", src[: m.start()].rstrip("(")):
                    continue
                if re.search(r"\bsignal\s+$", src[: m.start()]):
                    continue
                i, depth, args = m.end(), 1, ""
                while i < len(src) and depth:
                    depth += (src[i] == "(") - (src[i] == ")")
                    if depth:
                        args += src[i]
                    i += 1
                got = 0 if not args.strip() else top_level_commas(args) + 1
                if got != want:
                    line = src[: m.start()].count("\n") + 1
                    print(f"FAIL  {f.name}:{line}  {name}() declared with {want} "
                          f"argument(s), called with {got}")
                    bad += 1
    if bad == 0:
        print(f"  ok  signal calls match their declarations "
              f"({len(declared)} signals: {', '.join(sorted(declared))})")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "src/gui/qml")))
