#!/usr/bin/env python3
"""No docstring may be separated from what it documents.

Inserting a documented function immediately above another function's `function` line puts the
new block BETWEEN that function's docstring and its definition. Julia then tries to document a
string literal and the whole file fails to load:

    ERROR: LoadError: cannot document the following expression

It is a load-time failure, so nothing subtle — but it costs a precompile cycle to find, and it
happened three times in one sitting, always the same way: anchoring an insertion on
`function foo(` rather than on `\"\"\"\\n    foo(`.

    julia_docstrings.py <src-dir>
"""
import re
import sys
import pathlib


def main(root: pathlib.Path) -> int:
    bad = 0
    checked = 0
    for f in sorted(root.rglob("*.jl")):
        lines = f.read_text().split("\n")
        i = 0
        while i < len(lines):
            # A docstring block opening at column 0.
            if lines[i].strip() == '"""' and not lines[i].startswith(" "):
                j = i + 1
                while j < len(lines) and lines[j].strip() != '"""':
                    j += 1
                if j >= len(lines):
                    break
                checked += 1
                # What follows must be a definition, not another docstring or a bare string.
                k = j + 1
                while k < len(lines) and (not lines[k].strip() or
                                          lines[k].lstrip().startswith("#")):
                    k += 1
                if k < len(lines) and lines[k].strip() == '"""':
                    print(f"FAIL  {f.relative_to(root)}:{i + 1}  docstring is followed by "
                          f"another docstring, not a definition")
                    bad += 1
                i = j + 1
            else:
                i += 1
    if bad == 0:
        print(f"  ok  {checked} docstrings each attach to a definition")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "src")))
