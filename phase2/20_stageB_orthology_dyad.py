#!/usr/bin/env python3
"""
PHASE 2, STAGE B - run in the analysis sandbox on the unzipped phase2_upload.zip.
For every candidate metacaspase in every genome:
  (1) assign it to an OsMC group by RECIPROCAL BEST HIT against the 8 rice OsMC proteins (Smith-Waterman scores),
      with the score margin to the runner-up (so OsMC1 is never confused with OsMC2/OsMC3/OsMC8);
  (2) test CONSERVED GENOMIC CONTEXT: how many of the ~66 Nipponbare chromosome-3 window genes have their best hit within
      +-15 genes of the candidate;
  (3) read the residues at the two catalytic positions by direct motif search (alignment-independent);
  (4) write the OsMC1-ortholog table, the per-group dyad tables and FASTA sets for the phylogeny / selection stage.
"""
import sys, os, re, gzip, subprocess, collections, csv, tempfile
D = sys.argv[1] if len(sys.argv) > 1 else "phase2_out"
OUT = sys.argv[2] if len(sys.argv) > 2 else "phase2_results"
SW = os.environ.get("SW_BIN", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts", "sw"))
SEEDS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "seeds_all.fa")
os.makedirs(OUT, exist_ok=True)

def read_fa(p):
    d = collections.OrderedDict(); k = None
    for l in open(p):
        if l[0] == ">": k = l[1:].strip(); d[k] = []
        elif k: d[k].append(l.strip())
    return collections.OrderedDict((a, "".join(b)) for a, b in d.items())
seeds = read_fa(SEEDS); osmc = [k for k in seeds if k.startswith("OsMC")]
species = sorted(x for x in os.listdir(D) if os.path.exists(os.path.join(D, x, "candidates.faa")))
print(len(species), "genomes with candidate files")

def sw_scores(qfa, dfa, minscore=60):
    out = subprocess.run([SW, qfa, dfa, str(minscore)], capture_output=True, text=True, check=True).stdout
    sc = collections.defaultdict(dict)
    for l in out.splitlines():
        q, t, s, *_ = l.split("\t"); sc[t][q] = int(s)
    return sc
def motif(seq):
    """residues at the catalytic positions found by motif, requiring His-motif ... Cys-motif spacing of 45-70 aa"""
    best = None
    for m in re.finditer(r"[FYWLMI]S[GA](.)G", seq):
        hpos = m.start(1)
        for c in re.finditer(r"D[AST](.)[HN]S", seq[hpos+40:hpos+80]):
            cpos = hpos + 40 + c.start(1)
            best = (seq[hpos], hpos + 1, seq[cpos], cpos + 1, m.group(0), seq[cpos-2:cpos+3]); break
        if best: break
    return best

rows = []; allfa = collections.OrderedDict(); window = {}
tmp = tempfile.mkdtemp()
for sp in species:
    cand = read_fa(os.path.join(D, sp, "candidates.faa"))
    if not cand: continue
    qf = os.path.join(tmp, "q.fa"); df = os.path.join(tmp, "d.fa")
    with open(qf, "w") as f:
        for k in osmc: f.write(f">{k.split('|')[0]}\n{seeds[k]}\n")
    with open(df, "w") as f:
        for k, s in cand.items(): f.write(f">{k}\n{s}\n")
    sc = sw_scores(qf, df)
    # reciprocal best: for each OsMC seed, best candidate in this genome
    best_cand = {}
    for q in (k.split("|")[0] for k in osmc):
        ranked = sorted(((sc[t].get(q, 0), t) for t in cand), reverse=True)
        best_cand[q] = ranked[0][1] if ranked and ranked[0][0] > 0 else None
    # synteny window
    wh = {}
    wp = os.path.join(D, sp, "window_hits.tsv"); gop = os.path.join(D, sp, "gene_order.tsv.gz")
    idx = {}
    if os.path.exists(gop):
        for l in gzip.open(gop, "rt"):
            f = l.rstrip("\n").split("\t")
            if f[0] != "idx": idx[f[1]] = int(f[0])
    if os.path.exists(wp):
        for l in open(wp):
            f = l.rstrip("\n").split("\t")
            if f[0] != "nipponbare_gene" and f[4] in idx: wh[f[0]] = idx[f[4]]
    for t, s in cand.items():
        scores = sorted(((sc[t].get(q.split("|")[0], 0), q.split("|")[0]) for q in osmc), reverse=True)
        top, second = scores[0], scores[1]
        rbh = best_cand.get(top[1]) == t
        gene = t.split("|")[1]; gi = idx.get(gene)
        near = sum(1 for g, i in wh.items() if gi is not None and abs(i - gi) <= 15) if gi is not None else 0
        m = motif(s)
        ident = None
        rows.append(dict(species=sp, gene=gene, protein=t.split("|")[2], locus=t.split("|")[3], aa=len(s), best_OsMC=top[1], score=top[0],
                         runner_up=second[1], runner_score=second[0], margin=top[0] - second[0], RBH=rbh, window_genes_within15=near,
                         His_res=(m[0] if m else "?"), His_pos=(m[1] if m else ""), Cys_res=(m[2] if m else "?"), Cys_pos=(m[3] if m else ""),
                         His_motif=(m[4] if m else ""), Cys_motif=(m[5] if m else "")))
        allfa[f"{sp}|{gene}|{top[1]}"] = s
fields = list(rows[0].keys())
with open(os.path.join(OUT, "all_candidates_assigned.tsv"), "w") as f:
    f.write("\t".join(fields) + "\n")
    for r in rows: f.write("\t".join(str(r[k]) for k in fields) + "\n")
with open(os.path.join(OUT, "all_candidates.faa"), "w") as f:
    for k, s in allfa.items(): f.write(f">{k}\n{s}\n")

print("\n=== OsMC1 orthologs (reciprocal best hit with OsMC1, margin over runner-up) ===")
print(f"{'species':28s} {'gene':20s} {'aa':>4s} {'margin':>6s} {'ctx':>3s}  His  Cys   motifs")
n1 = both = hisonly = cysonly = intact = 0
for r in sorted(rows, key=lambda r: r["species"]):
    if r["best_OsMC"] == "OsMC1" and r["RBH"]:
        n1 += 1
        h, c = r["His_res"], r["Cys_res"]
        both += (h != "H" and c not in ("C",) and h != "?"); intact += (h == "H" and c == "C")
        print(f"{r['species']:28s} {r['gene']:20s} {r['aa']:4d} {r['margin']:6d} {r['window_genes_within15']:3d}   {h}    {c}    {r['His_motif']} / {r['Cys_motif']}")
print(f"\nOsMC1 orthologs found: {n1}; both catalytic residues substituted: {both}; both intact: {intact}")
print("\n=== per-group summary: species with RBH ortholog, and how many have intact H and C ===")
for g in ["OsMC%d" % i for i in range(1, 10)]:
    rr = [r for r in rows if r["best_OsMC"] == g and r["RBH"]]
    ok = sum(1 for r in rr if r["His_res"] == "H" and r["Cys_res"] == "C")
    print(f"{g}: {len(rr):2d} orthologs, {ok:2d} with intact His+Cys")
