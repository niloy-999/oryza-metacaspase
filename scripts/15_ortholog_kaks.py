"""Step 15 - constraint on the OsMC1 lineage: Ka/Ks between Nipponbare and each other genome's ortholog, compared with intact control groups
(OsMC5, OsMC6) computed for exactly the same genome pairs. Nei-Gojobori + Jukes-Cantor, codon alignment from protein alignment."""
import sys, os, math, csv, itertools, statistics as st, collections
sys.path.insert(0, os.path.dirname(__file__))
from common import *; from align import *
NT = "TCAG"
def sites(codon):
    S = 0.0
    for p in range(3):
        syn = 0
        for n in "ACGT":
            if n == codon[p]: continue
            c2 = codon[:p] + n + codon[p+1:]
            if CODON[c2] == "*": continue
            if CODON[c2] == CODON[codon]: syn += 1
        S += syn / 3.0
    return S
def diffs(c1, c2):
    pos = [i for i in range(3) if c1[i] != c2[i]]
    if not pos: return 0.0, 0.0
    Sd = Nd = 0.0; npath = 0
    for perm in itertools.permutations(pos):
        cur = c1; s = n = 0; ok = True
        for p in perm:
            nxt = cur[:p] + c2[p] + cur[p+1:]
            if CODON[nxt] == "*": ok = False; break
            if CODON[nxt] == CODON[cur]: s += 1
            else: n += 1
            cur = nxt
        if ok: Sd += s; Nd += n; npath += 1
    return (Sd / npath, Nd / npath) if npath else None
def jc(p): return None if p >= 0.75 else -0.75 * math.log(1 - 4 * p / 3)
def kaks(cds_a, cds_b):
    pa, pb = translate(cds_a).rstrip("*"), translate(cds_b).rstrip("*")
    sc, ra, rb, pairs = pair_align(pa, pb, "local")
    S = N = Sd = Nd = 0.0; n = 0
    for i, j in pairs:
        if i is None or j is None: continue
        ca, cb = cds_a[3*i:3*i+3], cds_b[3*j:3*j+3]
        if len(ca) < 3 or len(cb) < 3 or "*" in (CODON.get(ca), CODON.get(cb)): continue
        d = diffs(ca, cb)
        if d is None: continue
        sa, sb = sites(ca), sites(cb); S += (sa + sb) / 2; N += 3 - (sa + sb) / 2; Sd += d[0]; Nd += d[1]; n += 1
    if n < 200: return None
    Ks, Ka = jc(Sd / S), jc(Nd / N)
    return n, Ka, Ks, (Ka / Ks if Ks and Ks > 0 and Ka is not None else None)
D = os.environ.get("MC_P2DATA", os.path.join(OUT, "phase2", "phase2_out")); R = os.environ.get("MC_P2RES", os.path.join(OUT, "phase2", "results_stageB"))
rows = list(csv.DictReader(open(f"{R}/all_candidates_assigned.tsv"), delimiter="\t"))
cds = {}
for sp in os.listdir(D):
    p = f"{D}/{sp}/candidates.cds.fna"
    if os.path.exists(p):
        k = None
        for l in open(p):
            if l[0] == ">": k = l[1:].strip().split("|"); k = (k[0], k[1]); cds[k] = ""
            else: cds[k] += l.strip().upper()
def grp(g): return {r["species"]: r for r in rows if r["best_OsMC"] == g and r["RBH"] == "True"}
out = []; summary = {}
for g in ("OsMC1", "OsMC5", "OsMC6"):
    G = grp(g); ref = G["oryza_sativa"]; refc = cds[("oryza_sativa", ref["gene"])]
    vals = []
    for sp, r in G.items():
        if sp == "oryza_sativa": continue
        c = cds.get((sp, r["gene"]))
        if not c or len(c) % 3: continue
        res = kaks(refc, c)
        if res is None: continue
        n, Ka, Ks, ratio = res
        out.append((g, sp, n, Ka, Ks, ratio))
        if Ks is not None and 0.0005 < Ks < 1.0 and ratio is not None: vals.append((sp, Ks, ratio))
    summary[g] = vals
with open(f"{OUT}/results/tables/TableS15_ortholog_KaKs_vs_Nipponbare.tsv", "w") as f:
    f.write("group\tgenome\taligned_codons\tKa\tKs\tKa_Ks\n")
    for g, sp, n, Ka, Ks, ratio in out: f.write(f"{g}\t{sp}\t{n}\t{'NA' if Ka is None else round(Ka,4)}\t{'NA' if Ks is None else round(Ks,4)}\t{'NA' if ratio is None else round(ratio,3)}\n")
print("Ka/Ks of each genome's ortholog vs the Nipponbare ortholog (pairs with 0.0005 < Ks < 1 only)")
for g, v in summary.items():
    r = [x[2] for x in v]
    print(f"  {g}: n={len(r):2d} median Ka/Ks={st.median(r):.3f} mean={st.mean(r):.3f} range={min(r):.3f}-{max(r):.3f} | median Ks={st.median([x[1] for x in v]):.3f}")
# matched comparison: same genomes
common = set.intersection(*[{x[0] for x in summary[g]} for g in ("OsMC1", "OsMC5", "OsMC6")])
print("\nMatched genomes (present in OsMC1, OsMC5 and OsMC6 sets):", len(common))
for g in ("OsMC1", "OsMC5", "OsMC6"):
    r = [x[2] for x in summary[g] if x[0] in common]; print(f"  {g}: median Ka/Ks {st.median(r):.3f}")
# distant pairs
print("\nOsMC1 ortholog pairs with larger divergence:")
for sp, Ks, ratio in sorted(summary["OsMC1"], key=lambda x: -x[1])[:5]: print(f"  {sp:24s} Ks={Ks:.3f} Ka/Ks={ratio:.3f}")
