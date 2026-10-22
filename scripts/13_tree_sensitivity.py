"""Step 13 - sensitivity of the key topological claims to alignment trimming and distance model (NJ, 200 bootstraps per variant).
Claims tested: (A) Type I clade {AtMC1-3, OsMC1-3} monophyletic; (B) OsMC8 joins that clade (A+OsMC8); (C) {OsMC2,OsMC3,OsMC8} grouped;
(D) {OsMC2,OsMC3} grouped; (E) ratio of OsMC1 terminal branch to mean OsMC2/OsMC3 terminal branch."""
import sys, os, math, random
sys.path.insert(0, os.path.dirname(__file__))
from common import *
aln = read_fasta(f"{OUT}/results/alignments/core_alignment_AtMC9_anchored.fa"); names = list(aln); L = len(next(iter(aln.values())))
def dist(cols, model):
    n = len(names); D = np.zeros((n, n))
    for a in range(n):
        for b in range(a+1, n):
            m = d = 0
            for c in cols:
                x, y = aln[names[a]][c], aln[names[b]][c]
                if x != "-" and y != "-": m += 1; d += (x != y)
            p = d / m if m else 0.9
            if model == "p": v = p
            elif model == "poisson": v = 6.0 if p >= 0.95 else -math.log(1 - p)
            else: v = 5.0 if p >= 0.93 else -math.log(1 - p - 0.2 * p * p)
            D[a, b] = D[b, a] = v
    return D
def nj_tree(D):
    n = len(names); act = list(range(n)); Dm = {i: {j: D[i, j] for j in range(n)} for i in range(n)}; nxt = n
    leaves = {i: frozenset([names[i]]) for i in range(n)}; blen = {}; splits = set()
    while len(act) > 3:
        r = {i: sum(Dm[i][j] for j in act if j != i) for i in act}; m = len(act); best = None
        for ii, i in enumerate(act):
            for j in act[ii+1:]:
                q = (m - 2) * Dm[i][j] - r[i] - r[j]
                if best is None or q < best[0]: best = (q, i, j)
        _, i, j = best; dij = Dm[i][j]
        li = max(0, 0.5 * dij + (r[i] - r[j]) / (2 * (m - 2))); lj = max(0, dij - li)
        if len(leaves[i]) == 1: blen[next(iter(leaves[i]))] = li
        if len(leaves[j]) == 1: blen[next(iter(leaves[j]))] = lj
        leaves[nxt] = leaves[i] | leaves[j]; splits.add(leaves[nxt]); Dm[nxt] = {}
        for k in act:
            if k in (i, j): continue
            d = 0.5 * (Dm[i][k] + Dm[j][k] - dij); Dm[nxt][k] = d; Dm[k][nxt] = d
        act = [k for k in act if k not in (i, j)] + [nxt]; nxt += 1
    # last three nodes: terminal branch lengths for leaves among them
    i, j, k = act; d_ij, d_ik, d_jk = Dm[i][j], Dm[i][k], Dm[j][k]
    for node, dd in ((i, (d_ij + d_ik - d_jk) / 2), (j, (d_ij + d_jk - d_ik) / 2), (k, (d_ik + d_jk - d_ij) / 2)):
        if len(leaves[node]) == 1: blen[next(iter(leaves[node]))] = max(0, dd)
    return splits, blen
ALL = frozenset(names)
def has(splits, S):
    S = frozenset(S); return (S in splits) or ((ALL - S) in splits)
TI = ["AtMC1", "AtMC2", "AtMC3", "OsMC1", "OsMC2", "OsMC3"]
claims = {"A_TypeI_clade": TI, "B_TypeI+OsMC8": TI + ["OsMC8"], "C_OsMC2+3+8": ["OsMC2", "OsMC3", "OsMC8"], "D_OsMC2+3": ["OsMC2", "OsMC3"], "E_TypeI+OsMC9": TI + ["OsMC9"], "F_TypeI+OsMC8+OsMC9": TI + ["OsMC8", "OsMC9"], "G_OsMC9+AtMC1+2": ["OsMC9", "AtMC1", "AtMC2"]}
random.seed(7); NB = 200; rows = []
print(f"{'trim':>5s} {'dist':>8s} cols  " + "  ".join(f"{c:>14s}" for c in claims) + "  OsMC1/OsMC2-3_branch")
for occ in (0.5, 0.7, 0.9):
    keep = [c for c in range(L) if sum(aln[n][c] != "-" for n in names) >= occ * len(names)]
    for model in ("p", "poisson", "kimura"):
        base_sp, base_bl = nj_tree(dist(keep, model))
        support = {c: 0 for c in claims}
        for b in range(NB):
            cols = [random.choice(keep) for _ in keep]
            sp, _ = nj_tree(dist(cols, model))
            for c, S in claims.items(): support[c] += has(sp, S)
        ratio = base_bl["OsMC1"] / ((base_bl["OsMC2"] + base_bl["OsMC3"]) / 2)
        obs = {c: ("yes" if has(base_sp, S) else "no") for c, S in claims.items()}
        line = "  ".join(f"{obs[c]:>3s}/{100*support[c]/NB:5.0f}%   " for c in claims)
        print(f"{occ:5.1f} {model:>8s} {len(keep):4d}  {line} {ratio:5.2f}")
        rows.append((occ, model, len(keep), ratio, {c: (obs[c], support[c] / NB) for c in claims}))
with open(f"{OUT}/results/tables/S12_tree_sensitivity.tsv", "w") as f:
    f.write("occupancy\tdistance\tcolumns\tOsMC1_branch_ratio\t" + "\t".join(f"{c}_observed\t{c}_bootstrap" for c in claims) + "\n")
    for occ, model, nc, ratio, res in rows:
        f.write(f"{occ}\t{model}\t{nc}\t{ratio:.3f}\t" + "\t".join(f"{res[c][0]}\t{res[c][1]:.3f}" for c in claims) + "\n")
