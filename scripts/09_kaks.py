"""Step 9 - Ka/Ks (Nei-Gojobori 1986, Jukes-Cantor correction) for all 28 OsMC1-8 pairs, from codon alignments derived from the
protein alignment (PAL2NAL logic) of CDS pulled directly from genome+GFF3. Pathway averaging for multi-hit codons; stop-involving pathways excluded."""
import sys, os, math, itertools
sys.path.insert(0, os.path.dirname(__file__))
from common import *; from align import *; from mcdata import *

genes, order = parse_gff(RICE_GFF)
chrc = {}
def cds_of(locus):
    g = genes[locus]; mid = [m for m, v in g["mrna"].items() if v["longest"]][0]
    if g["chr"] not in chrc: chrc[g["chr"]] = load_chrom(RICE_GEN, g["chr"])
    return gene_cds(genes, locus, mid, chrc[g["chr"]])
NT = "TCAG"
def is_syn(c1, c2): return CODON[c1] == CODON[c2]
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
    if npath == 0: return None
    return Sd / npath, Nd / npath
def jc(p): return None if p >= 0.75 else -0.75 * math.log(1 - 4 * p / 3)

cds = {k: cds_of(v[0]) for k, v in OSMC.items()}
prot = {k: translate(c).rstrip("*") for k, c in cds.items()}
rows = []
names = list(OSMC)
for a, b in itertools.combinations(names, 2):
    sc, ra, rb, pairs = pair_align(prot[a], prot[b], "local")
    S = N = Sd = Nd = 0.0; ncod = 0
    for i, j in pairs:
        if i is None or j is None: continue
        ca, cb = cds[a][3*i:3*i+3], cds[b][3*j:3*j+3]
        if "*" in (CODON.get(ca), CODON.get(cb)) or len(ca) < 3 or len(cb) < 3: continue
        d = diffs(ca, cb)
        if d is None: continue
        s_a = sites(ca); s_b = sites(cb)
        S += (s_a + s_b) / 2; N += 3 - (s_a + s_b) / 2; Sd += d[0]; Nd += d[1]; ncod += 1
    pS, pN = Sd / S, Nd / N
    Ks, Ka = jc(pS), jc(pN)
    ratio = (Ka / Ks) if (Ks and Ka is not None and Ks > 0) else None
    sat = (Ks is None) or (Ks > 2.0)
    rows.append(dict(pair=f"{a}-{b}", type_pair=f"{OSMC[a][1]}-{OSMC[b][1]}", codons=ncod, Ka=Ka, Ks=Ks, KaKs=ratio, saturated=sat,
                     same_chr=(genes[OSMC[a][0]]["chr"] == genes[OSMC[b][0]]["chr"])))
def f(x): return "NA" if x is None else f"{x:.3f}"
with open(f"{OUT}/results/tables/Table7_KaKs_all_pairs.tsv", "w") as fh:
    fh.write("pair\ttype_pair\taligned_codons\tKa\tKs\tKa_Ks\tKs_saturated(>2 or undefined)\tsame_chromosome\n")
    for r in rows: fh.write(f"{r['pair']}\t{r['type_pair']}\t{r['codons']}\t{f(r['Ka'])}\t{f(r['Ks'])}\t{f(r['KaKs'])}\t{r['saturated']}\t{r['same_chr']}\n")
for r in rows: print(f"{r['pair']:14s} {r['type_pair']:5s} n={r['codons']:3d} Ka={f(r['Ka'])} Ks={f(r['Ks'])} Ka/Ks={f(r['KaKs'])} sat={r['saturated']} sameChr={r['same_chr']}")
