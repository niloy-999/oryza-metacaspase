#!/usr/bin/env python3
"""Label the clade of all tips whose name starts with PREFIX as {Test} (all branches inside the clade) and everything else {Reference},
for HyPhy RELAX. The tree is re-rooted on the branch that separates the clade, so the clade must be monophyletic in the tree.
Usage: label_tree.py in.treefile PREFIX out.nwk"""
import sys, re
def parse(s):
    s = s.strip().rstrip(";"); pos = [0]
    def node():
        if s[pos[0]] == "(":
            pos[0] += 1; ch = []
            while True:
                ch.append(node())
                if s[pos[0]] == ",": pos[0] += 1; continue
                break
            pos[0] += 1; lab = ""
            while pos[0] < len(s) and s[pos[0]] not in ":,);": lab += s[pos[0]]; pos[0] += 1
            n = dict(ch=ch, name="", bl=0.0)
        else:
            st = pos[0]
            while s[pos[0]] not in ":,);": pos[0] += 1
            n = dict(ch=[], name=s[st:pos[0]], bl=0.0)
        if pos[0] < len(s) and s[pos[0]] == ":":
            pos[0] += 1; st = pos[0]
            while pos[0] < len(s) and s[pos[0]] not in ",);": pos[0] += 1
            n["bl"] = float(s[st:pos[0]])
        return n
    return node()
def leaves(n): return [n["name"]] if not n["ch"] else [x for c in n["ch"] for x in leaves(c)]
def to_nw(n, label, top=False):
    if not n["ch"]: return f"{n['name']}{{{label}}}:{n['bl']:.8f}"
    return "(" + ",".join(to_nw(c, label) for c in n["ch"]) + f"){{{label}}}:{n['bl']:.8f}"
def main():
    tree = parse(open(sys.argv[1]).read()); pre = sys.argv[2]
    tips = [t for t in leaves(tree) if t.startswith(pre)]; tipset = set(tips)
    # unrooted adjacency
    adj = {}; names = {}; cnt = [0]
    def build(n, parent):
        i = cnt[0]; cnt[0] += 1; adj[i] = []; names[i] = n["name"]
        if parent is not None: adj[i].append((parent[0], parent[1])); adj[parent[0]].append((i, parent[1]))
        for c in n["ch"]: build(c, (i, c["bl"]))
        return i
    build(tree, None)
    def side(u, v):                     # leaves reachable from v without crossing the edge (u,v)
        stack = [(v, u)]; out = set()
        while stack:
            x, p = stack.pop()
            nb = [(y, w) for y, w in adj[x] if y != p]
            if not nb: out.add(names[x])
            for y, w in nb: stack.append((y, x))
        return out
    best = None
    for u in adj:
        for v, w in adj[u]:
            if side(u, v) == tipset: best = (u, v, w); break
        if best: break
    if not best: sys.exit(f"The {len(tips)} {pre} tips do not form a clade in this tree; cannot label for RELAX.")
    u, v, w = best
    def emit(x, p, label):
        nb = [(y, bl) for y, bl in adj[x] if y != p]
        if not nb: return f"{names[x]}{{{label}}}"
        return "(" + ",".join(emit(y, x, label) + ":" + f"{bl:.8f}" for y, bl in nb) + "){" + label + "}"
    test_side = emit(v, u, "Test"); ref_side = emit(u, v, "Reference")
    open(sys.argv[3], "w").write(f"({test_side}:{w/2:.8f},{ref_side}:{w/2:.8f});\n")
    print(f"{len(tips)} tips labelled Test; {len(leaves(tree))-len(tips)} Reference; root placed on the clade stem branch.")
main()
