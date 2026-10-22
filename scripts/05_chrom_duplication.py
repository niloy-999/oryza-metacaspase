"""Step 5 - chromosomal localization + tandem/proximal duplication (within O. sativa only - the only genome present in the
uploaded project). Tandem duplicate = paralog pair on the same chromosome, <=10 intervening genes apart (Holub 2001 criterion),
confirmed by >=60% protein identity over ≥70% aligned length (Gu et al. 2002)."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from common import *; from align import *; from mcdata import *
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

genes, order = parse_gff(RICE_GFF)
CHRLEN = {}
for line in open(RICE_GEN):
    pass
import re
cur = None
with open(RICE_GEN) as f:
    ln = 0
    for line in f:
        if line[0] == ">":
            if cur: CHRLEN[cur] = ln
            cur = line[1:].split()[0]; ln = 0
        else: ln += len(line.strip())
    if cur: CHRLEN[cur] = ln

rice = read_fasta(RICE_PROT)
mc = {k: (loc, rice[loc + ".1"].rstrip("*")) for k, (loc, t) in OSMC.items()}
# gene rank on each chromosome (index among ALL annotated genes, for the "<=10 genes apart" tandem criterion)
rank = {}
for chrom, lst in order.items():
    for i, (s, e, gid) in enumerate(lst): rank[gid] = i

pairs = []
names = list(mc)
for i in range(len(names)):
    for j in range(i+1, len(names)):
        a, b = names[i], names[j]; la, sa = mc[a]; lb, sb = mc[b]
        if genes[la]["chr"] != genes[lb]["chr"]: continue
        sc, ra, rb, al = pair_align(sa, sb, "local"); ident = identity(ra, rb)
        cov = sum(x != "-" and y != "-" for x, y in zip(ra, rb)) / min(len(sa), len(sb))
        gdist = abs(rank[la] - rank[lb]) - 1
        kind = "tandem" if (gdist <= 10 and ident >= 0.60 and cov >= 0.70) else ("proximal" if gdist <= 20 else "dispersed (same chr)")
        pairs.append(dict(pair=f"{a}-{b}", chrom=genes[la]["chr"], intervening_genes=gdist, identity_pct=round(100*ident,1),
                          aligned_coverage_pct=round(100*cov,1), bp_distance=abs(genes[la]["start"]-genes[lb]["start"]), duplication_type=kind))
with open(f"{OUT}/results/tables/Table6_duplication_pairs_within_Osativa.tsv", "w") as f:
    ks = list(pairs[0].keys()); f.write("\t".join(ks)+"\n")
    for p in pairs: f.write("\t".join(str(p[k]) for k in ks)+"\n")
for p in pairs: print(p)

# chromosome map figure
fig, ax = plt.subplots(figsize=(4.6, 7.5))
chrom_order = [f"Chr{i}" for i in range(1,13)]
for i, c in enumerate(chrom_order):
    L = CHRLEN.get(c, 0)
    ax.add_patch(plt.Rectangle((i-0.15, 0), 0.3, L/1e6, color="#d8d8d8", zorder=1))
    ax.text(i, -1.5, c.replace("Chr",""), ha="center", fontsize=9)
# stack labels that would otherwise collide (same chromosome, close together)
by_chr = {}
for k, (loc, seq) in mc.items():
    g = genes[loc]; by_chr.setdefault(g["chr"], []).append((g["start"]/1e6, k))
label_y = {}
for c, items in by_chr.items():
    items.sort()
    ys = [y for y, k in items]
    placed = []
    for y, k in items:
        py = y
        while any(abs(py - p) < 1.6 for p in placed): py += 1.6
        placed.append(py); label_y[k] = py
for k, (loc, seq) in mc.items():
    g = genes[loc]; ci = chrom_order.index(g["chr"]); y = g["start"]/1e6
    col = "#c0392b" if OSMC[k][1] == "I" else "#1f5fa8"
    ax.plot([ci-0.15, ci+0.15], [y, y], color=col, lw=2, zorder=3)
    ly = label_y[k]
    ax.plot([ci+0.16, ci+0.34], [y, ly], color=col, lw=0.5, zorder=3) if abs(ly-y)>0.05 else None
    ax.text(ci+0.38, ly, k, fontsize=7.5, va="center", color=col)
for p in pairs:
    if p["duplication_type"] == "tandem":
        a, b = p["pair"].split("-"); ga, gb = genes[mc[a][0]], genes[mc[b][0]]
        ci = chrom_order.index(ga["chr"]); y1, y2 = ga["start"]/1e6, gb["start"]/1e6
        ax.plot([ci-0.35, ci-0.35], [y1, y2], color="#666", lw=1, ls="--", zorder=2)
ax.set_xlim(-0.6, len(chrom_order)-0.3); ax.set_ylim(max(CHRLEN.values())/1e6+3, -3)
ax.set_ylabel("Position (Mb)"); ax.set_xticks([]); 
for sp in ("top","right"): ax.spines[sp].set_visible(False)
ax.set_title("Chromosomal distribution of rice metacaspase genes\nred=Type I, blue=Type II, dashed=tandem duplicate pair", fontsize=9)
plt.tight_layout(); plt.savefig(f"{OUT}/results/figures/Fig4_chromosome_distribution.png", dpi=300); plt.close()
print("\nchr lens (Mb):", {k: round(v/1e6,1) for k,v in CHRLEN.items() if k.startswith("Chr") and len(k)<=5})
