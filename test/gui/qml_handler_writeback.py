#!/usr/bin/env python3
"""Find property-change handlers that write back the property they react to.

THE BUG THIS EXISTS FOR. The fits list had

    onCurrentIndexChanged: { ...; root.refreshFits() }        // reacts to currentIndex

and `refreshFits()` ended with

    fitList.currentIndex = cur - 1                            // writes currentIndex

Its `fitModel.clear()` resets the index to -1 on every pass, so the restore is always a real
change, the handler fires again, and the two recurse until the engine gives up with

    RangeError: Maximum call stack size exceeded. (ModelTab.qml:1455181072)

— at a line number that does not exist, because the frame is synthetic. The aborted pass also
leaves the ListModel half-filled, which draws as a stray blank row. Nothing in the Julia-side
test suite can see this: it never evaluates the QML.

WHAT IS REPORTED. For each `on<Prop>Changed` handler, the functions it can reach inside the
same file are followed, and an assignment to `<Prop>` on the same object is a cycle. This does
not try to prove the cycle is unsafe: a re-entrancy guard makes it fine, and that is the usual
fix. So a cycle must be ANNOTATED to pass, with a comment anywhere in the handler body:

    // SAFE-WRITEBACK: guarded by root.syncingFits

which is both the opt-out and the record of why. An unannotated cycle fails.

    qml_handler_writeback.py <qml-dir>
"""
import re
import sys
import pathlib

MARKER = "SAFE-WRITEBACK"


def block_at(text: str, open_brace: int) -> str:
    """The source from `open_brace` to its matching close, braces balanced."""
    depth = 0
    for i in range(open_brace, len(text)):
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
            if depth == 0:
                return text[open_brace:i + 1]
    return text[open_brace:]


def functions(text: str) -> dict:
    """Every `function name(...) { ... }` in the file, by name."""
    out = {}
    for m in re.finditer(r"\bfunction\s+(\w+)\s*\([^)]*\)\s*\{", text):
        out[m.group(1)] = block_at(text, m.end() - 1)
    return out


def enclosing_id(text: str, pos: int) -> str:
    """The `id:` of the object the handler at `pos` belongs to.

    Scanned BACKWARDS for the nearest `id:` that is not closed off before `pos`, which is what
    an indentation-free brace walk would need a parser for. The nearest preceding `id:` at a
    shallower or equal brace depth is right in every shape this codebase uses: a handler sits
    inside the object whose id was declared a few lines above it.
    """
    ids = [(m.start(), m.group(1)) for m in re.finditer(r"^\s*id:\s*(\w+)", text[:pos], re.M)]
    return ids[-1][1] if ids else ""


def writes(body: str, prop: str, obj_id: str) -> bool:
    """Does `body` assign `prop`, bare or through `obj_id`?"""
    # `x.prop = ...` for the same object, or a bare `prop = ...` (same object's scope).
    pats = [rf"\b{re.escape(obj_id)}\.{re.escape(prop)}\s*=(?!=)"] if obj_id else []
    pats.append(rf"(?<![.\w]){re.escape(prop)}\s*=(?!=)")
    return any(re.search(p, body) for p in pats)


def reachable(body: str, funcs: dict, seen=None) -> str:
    """`body` plus the bodies of every local function it can call, transitively."""
    seen = seen if seen is not None else set()
    out = [body]
    for m in re.finditer(r"(?:\broot\.|\bwin\.|\b)(\w+)\s*\(", body):
        name = m.group(1)
        if name in funcs and name not in seen:
            seen.add(name)
            out.append(reachable(funcs[name], funcs, seen))
    return "\n".join(out)


def main(qmldir: pathlib.Path) -> int:
    files = sorted(qmldir.glob("*.qml"))
    if not files:
        print(f"FAIL  no .qml files under {qmldir}")
        return 1

    bad, checked, annotated = [], 0, 0
    for f in files:
        text = f.read_text()
        funcs = functions(text)
        for m in re.finditer(r"^([ \t]*)on([A-Z]\w*)Changed\s*:", text, re.M):
            prop = m.group(2)[0].lower() + m.group(2)[1:]
            brace = text.find("{", m.end())
            # A one-line handler with no block: take the rest of the line.
            if brace < 0 or "\n" in text[m.end():brace]:
                body = text[m.end():text.find("\n", m.end())]
            else:
                body = block_at(text, brace)
            checked += 1
            obj_id = enclosing_id(text, m.start())
            whole = reachable(body, funcs)
            if not writes(whole, prop, obj_id):
                continue
            if MARKER in whole:
                annotated += 1
                continue
            line = text[:m.start()].count("\n") + 1
            bad.append(f"{f.name}:{line}: on{m.group(2)}Changed can reach an assignment to "
                       f"`{prop}`" + (f" on `{obj_id}`" if obj_id else "") +
                       f" — a write-back cycle. Guard it and annotate with "
                       f"`// {MARKER}: <why>`, or move the write out of the handler's reach.")

    if bad:
        print(f"FAIL  {len(bad)} handler(s) can write back the property they react to:")
        for b in bad:
            print("      " + b)
        return 1
    print(f"  ok  {checked} property-change handlers, none writes back unguarded"
          f" ({annotated} annotated)")
    return 0


if __name__ == "__main__":
    sys.exit(main(pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else ".")))
