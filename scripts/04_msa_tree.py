"""Step 4 - reference-anchored multiple alignment of the catalytic core (all sequences aligned to AtMC9; insertions relative to AtMC9 dropped),
Poisson-corrected distances, neighbour-joining tree with 1000 bootstrap replicates (column resampling), Newick + PNG.
Tree is rooted at the Type I / Type II split (midpoint) - documented limitation; a formal analysis should use MAFFT + IQ-TREE 2 (pipeline_oryza_wide/)."""
import sys, os, random, math
sys.path.insert(0, os.path.dirname(__file__))
from common import *; from align import *; from mcdata import *
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

rice = read_fasta(RICE_PROT); at = read_fasta(AT_PROT, keyfn=lambda h: h.split()[0].split("|")[0])
seqs = OrderedDict()
for k, (loc, t) in OSMC.items(): seqs[f"{k}"] = rice[loc + ".1"].rstrip("*")
for k, loc in AT_MC.items(): seqs[k] = at[loc + ".1"].rstrip("*")
_ex = f"{OUT}/results/sequences/OsMC9_Os10g0565100.fa"
if os.path.exists(_ex): seqs["OsMC9"] = list(read_fasta(_ex).values())[0]
NEG_CTRL = rice["LOC_Os09g36640.1"].rstrip("*")
ref = seqs["AtMC9"]; L = len(ref)
aln = OrderedDict()
for k, s in seqs.items():
    sc, ra, rb, pairs = pair_align(ref, s, "local")
    col = ["-"] * L
    for i, j in pairs:
        if i is not None and j is not None: col[i] = s[j]
    aln[k] = "".join(col)
# columns where >=70% of sequences have a residue (gap-tolerant trimming), on the p20/p10 core
keep = [c for c in range(L) if sum(aln[k][c] != "-" for k in aln) >= 0.7 * len(aln)]
trim = OrderedDict((k, "".join(v[c] for c in keep)) for k, v in aln.items())
write_fasta(aln, f"{OUT}/results/alignments/core_alignment_AtMC9_anchored.fa"); write_fasta(trim, f"{OUT}/results/alignments/core_alignment_trimmed.fa")
print("aligned columns:", L, "kept:", len(keep))
names = list(trim)

def dist_matrix(cols):
    n = len(names); D = np.zeros((n, n))
    for a in range(n):
        for b in range(a+1, n):
            m = d = 0
            for c in cols:
                x, y = trim[names[a]][c], trim[names[b]][c]
                if x != "-" and y != "-": m += 1; d += (x != y)
            p = d / m if m else 0.75
            D[a, b] = D[b, a] = 5.0 if p >= 0.93 else -math.log(1 - p - 0.2 * p * p)   # Kimura 1983 protein distance
    return D

def nj(D, names):
    n = len(names); nodes = [(nm, 0.0) for nm in names]; D = D.copy(); active = list(range(n)); trees = {i: names[i] for i in range(n)}
    Dm = {i: {j: D[i, j] for j in range(n)} for i in range(n)}; nxt = n; bl = {}
    while len(active) > 2:
        r = {i: sum(Dm[i][j] for j in active if j != i) for i in active}; m = len(active)
        best = None
        for ii, i in enumerate(active):
            for j in active[ii+1:]:
                q = (m - 2) * Dm[i][j] - r[i] - r[j]
                if best is None or q < best[0]: best = (q, i, j)
        _, i, j = best; dij = Dm[i][j]
        li = 0.5 * dij + (r[i] - r[j]) / (2 * (m - 2)); lj = dij - li
        li, lj = max(li, 0), max(lj, 0)
        trees[nxt] = f"({trees[i]}:{li:.5f},{trees[j]}:{lj:.5f})"
        Dm[nxt] = {}
        for k in active:
            if k in (i, j): continue
            d = 0.5 * (Dm[i][k] + Dm[j][k] - dij); Dm[nxt][k] = d; Dm[k][nxt] = d
        active = [k for k in active if k not in (i, j)] + [nxt]; nxt += 1
    i, j = active; d = Dm[i][j]
    return f"({trees[i]}:{d/2:.5f},{trees[j]}:{d/2:.5f});"

# --- tree parsing / plotting ---
def parse_newick(s):
    pos = [0]
    def node():
        if s[pos[0]] == "(":
            pos[0] += 1; ch = []
            while True:
                ch.append(node())
                if s[pos[0]] == ",": pos[0] += 1; continue
                break
            assert s[pos[0]] == ")"; pos[0] += 1
            lab = ""; 
            while s[pos[0]] not in ":,);": lab += s[pos[0]]; pos[0] += 1
            n = dict(children=ch, name=lab, bl=0.0)
        else:
            st = pos[0]
            while s[pos[0]] not in ":,);": pos[0] += 1
            n = dict(children=[], name=s[st:pos[0]], bl=0.0)
        if s[pos[0]] == ":":
            pos[0] += 1; st = pos[0]
            while s[pos[0]] not in ",);": pos[0] += 1
            n["bl"] = float(s[st:pos[0]])
        return n
    return node()
def leaves(n): return [n["name"]] if not n["children"] else [x for c in n["children"] for x in leaves(c)]
def splits(n, allset, out):
    if n["children"]:
        for c in n["children"]: splits(c, allset, out)
        ls = frozenset(leaves(n))
        if 1 < len(ls) < len(allset) - 1:
            out.add(min(ls, allset - ls, key=lambda x: (len(x), sorted(x))))
    return out

D0 = dist_matrix(range(len(keep))); nw = nj(D0, names); tree = parse_newick(nw); allset = frozenset(names)
ref_splits = splits(tree, allset, set()); support = {s: 0 for s in ref_splits}
random.seed(20250928); NB = 1000
for b in range(NB):
    cols = [random.randrange(len(keep)) for _ in keep]
    tb = parse_newick(nj(dist_matrix(cols), names)); sp = splits(tb, allset, set())
    for s in ref_splits:
        if s in sp: support[s] += 1
# annotate support onto nodes
def annotate(n):
    for c in n["children"]: annotate(c)
    if n["children"]:
        ls = frozenset(leaves(n)); key = min(ls, allset - ls, key=lambda x: (len(x), sorted(x)))
        n["support"] = round(100 * support.get(key, 0) / NB) if key in support else None
annotate(tree)
def to_nw(n):
    if not n["children"]: return f"{n['name']}:{n['bl']:.5f}"
    lab = "" if n.get("support") is None else str(n["support"])
    return "(" + ",".join(to_nw(c) for c in n["children"]) + f"){lab}:{n['bl']:.5f}"
open(f"{OUT}/results/trees/core_NJ_bootstrap1000.nwk", "w").write(to_nw(tree) + ";\n")

# layout
y = [0]; pos = {}
def layout(n, x):
    n["x"] = x + n["bl"]
    if not n["children"]: n["y"] = y[0]; y[0] += 1
    else:
        for c in n["children"]: layout(c, n["x"])
        n["y"] = sum(c["y"] for c in n["children"]) / len(n["children"])
layout(tree, 0.0)
typ = {**{k: v[1] for k, v in OSMC.items()}, **AT_TYPE, "OsMC9": "I"}
fig, ax = plt.subplots(figsize=(8.5, 6.2))
def draw(n):
    for c in n["children"]:
        ax.plot([n["x"], n["x"]], [n["y"], c["y"]], "k", lw=1); ax.plot([n["x"], c["x"]], [c["y"], c["y"]], "k", lw=1); draw(c)
    if n["children"] and n.get("support") is not None and n["support"] >= 50:
        ax.text(n["x"] - 0.004, n["y"] + 0.12, str(n["support"]), fontsize=7, ha="right", color="#a00")
    if not n["children"]:
        col = {"I": "#c0392b", "II": "#1f5fa8"}[typ[n["name"]]]; sp = "#1b7f3b" if n["name"].startswith("Os") else "#555"
        ax.text(n["x"] + 0.01, n["y"], f"{n['name']} (Type {typ[n['name']]})", va="center", fontsize=9, color=col, fontweight="bold" if n["name"].startswith("Os") else "normal")
draw(tree); ax.set_yticks([]); ax.invert_yaxis(); ax.set_xlabel("Kimura-corrected distance (substitutions/site)")
for sp in ("top", "right", "left"): ax.spines[sp].set_visible(False)
ax.set_xlim(0, max(l["x"] for l in [tree]) + 1.4 if False else None)
xs = [];  
def collect(n): xs.append(n["x"]); [collect(c) for c in n["children"]]
collect(tree); ax.set_xlim(0, max(xs) * 1.55)
ax.set_title("Catalytic-core NJ tree of rice (Os) and Arabidopsis (At) metacaspases\nred = Type I, blue = Type II; numbers = bootstrap % (1000 replicates, shown if ≥50)", fontsize=10)
plt.tight_layout(); plt.savefig(f"{OUT}/results/figures/Fig2_phylogeny_NJ.png", dpi=300); plt.close()
print(to_nw(tree))
print("\nsupport (>=50):", sorted([(round(100*v/NB), sorted(k)[:3]) for k, v in support.items() if v/NB >= .5], reverse=True))
# pairwise identity table
with open(f"{OUT}/results/tables/S3_pairwise_identity_core.tsv", "w") as f:
    f.write("\t" + "\t".join(names) + "\n")
    for a in names:
        row = []
        for b in names:
            m = d = 0
            for x, y2 in zip(trim[a], trim[b]):
                if x != "-" and y2 != "-": m += 1; d += (x == y2)
            row.append(f"{100*d/m:.1f}" if m else "NA")
        f.write(a + "\t" + "\t".join(row) + "\n")
