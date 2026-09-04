#!/usr/bin/env python3
"""Every model role a delegate reads must be one its ListModel actually appends.

A delegate that names a role its model never supplies does not warn: the identifier is simply
unresolved, and the handler containing it throws a ReferenceError at the moment it runs. If
that handler is a signal handler — `onActivated`, `onEditingFinished` — the failure is silent
and total: the control looks like it worked, and nothing behind it moved.

That happened for real. `paramRow` is shared by three ListViews and reads `pcomp` to know which
component it is editing; only one of the three models appended it, so every edit in the other
two forms threw before reaching Julia.
"""
import re, sys, pathlib

def roles_of(src, model):
    """Role names appended to `model`, as a set per append call."""
    out = []
    for m in re.finditer(re.escape(model) + r"\.append\(\{", src):
        i, depth, body = m.end(), 1, ""
        while i < len(src) and depth:
            depth += (src[i] == "{") - (src[i] == "}")
            if depth:
                body += src[i]
            i += 1
        # Comments first: a `//` line between two roles breaks the "comma then name" shape
        # the scan below looks for, and would report a role that is plainly there.
        flat = re.sub(r"//[^\n]*", "", body)
        while True:
            new = re.sub(r"\([^()]*\)|\{[^{}]*\}|\[[^\[\]]*\]", "", flat)
            if new == flat:
                break
            flat = new
        out.append({k.strip() for k in re.findall(r"(?:^|,)\s*(\w+)\s*:", flat)})
    return out

def block_after(src, idx):
    """The braced block starting at or after `idx`."""
    i = src.index("{", idx); depth, j = 1, i + 1
    while j < len(src) and depth:
        depth += (src[j] == "{") - (src[j] == "}")
        j += 1
    return src[i:j]

bad = 0
for f in sorted(pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "src/gui/qml").glob("*.qml")):
    src = f.read_text()
    models = re.findall(r"ListModel\s*\{\s*id:\s*(\w+)", src)
    appended = {m: roles_of(src, m) for m in models}
    known = {r for sets in appended.values() for s in sets for r in s}
    # Named Components, so a shared delegate can be resolved to its body.
    comps = {}
    for m in re.finditer(r"Component\s*\{\s*id:\s*(\w+)", src):
        comps[m.group(1)] = block_after(src, m.end())
    # Every view that binds a model and a delegate.
    for m in re.finditer(r"(ListView|Repeater)\s*\{", src):
        body = block_after(src, m.start())
        mm = re.search(r"\bmodel:\s*(\w+)\s*$", body, re.M)
        dm = re.search(r"\bdelegate:\s*(\w+)\s*$", body, re.M)
        if not mm or mm.group(1) not in appended:
            continue
        model = mm.group(1)
        text = comps.get(dm.group(1), body) if dm else body
        used = set(re.findall(r"(?<![\w.])(\w+)(?![\w(])", text)) & known
        for s in appended[model]:
            missing = sorted(used - s)
            if missing:
                line = src[:m.start()].count("\n") + 1
                print(f"FAIL  {f.name}:{line}  delegate on `{model}` reads "
                      f"{', '.join(missing)} — never appended")
                bad += 1
print("  ok  every delegate role is supplied by its model" if not bad else "")
sys.exit(1 if bad else 0)
