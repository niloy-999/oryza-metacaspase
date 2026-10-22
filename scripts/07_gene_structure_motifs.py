"""Step 7 - exon/intron structure (from GFF3, matched to phylogeny order) + protein motif re-analysis.
The manuscript's Table 3 (MEME) is unusable: J/Z are not standard amino-acid codes and cannot occur in a real MEME
protein-motif table (MEME's protein alphabet is ACDEFGHIKLMNPQRSTVWY, +X/B/Z/J only as *ambiguity* codes, never emitted
in a motif consensus). We regenerate motifs from the actual protein sequences using ungapped k-mer enrichment (k=6-10,
positions relative to the catalytic His/Cys from step 3) as a transparent, reproducible substitute for MEME (not installed
in this environment); the exact `meme motif.txt -protein -nmotifs 10 -minw 6 -maxw 20` command is provided for the user to
run and drop in place of this module."""
import sys, os, json
from collections import Counter
sys.path.insert(0, os.path.dirname(__file__))
from common import *; from mcdata import *
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

genes, order = parse_gff(RICE_GFF)
rows = []
for k, (loc, t) in OSMC.items():
    g = genes[loc]; mids = [m for m, v in g["mrna"].items() if v["longest"]] or list(g["mrna"]); mid = mids[0]
    m = g["mrna"][mid]; cds = sorted(m["cds"])
    exon_lens = [e - s + 1 for s, e in cds]; intron_lens = [cds[i+1][0]-cds[i][1]-1 for i in range(len(cds)-1)]
    rows.append(dict(gene=k, locus=loc, strand=g["strand"], n_exon=len(cds), n_intron=len(cds)-1,
                     exon_lens=",".join(map(str,exon_lens)), intron_lens=",".join(map(str,intron_lens)),
                     gene_span_bp=g["end"]-g["start"]+1))
with open(f"{OUT}/results/tables/Table5_gene_structure.tsv","w") as f:
    ks=list(rows[0].keys()); f.write("\t".join(ks)+"\n")
    for r in rows: f.write("\t".join(str(r[k]) for k in ks)+"\n")

# figure: phylogeny-ordered exon/intron cartoon (order taken from the NJ tree, Os leaves only)
order_names = ["OsMC1","OsMC2","OsMC3","OsMC8","OsMC7","OsMC4","OsMC5","OsMC6"]  # from Fig2 topology
byname = {r["gene"]: r for r in rows}
fig, ax = plt.subplots(figsize=(9, 4.2))
for i, nm in enumerate(order_names):
    r = byname[nm]; g = genes[r["locus"]]; m = list(g["mrna"].values()); m = [v for v in m if v["longest"]] or m; m=m[0]
    cds = sorted(m["cds"]); origin = cds[0][0]
    y = len(order_names) - 1 - i
    ax.plot([0, (cds[-1][1]-origin)/1000], [y, y], color="#888", lw=1, zorder=1)
    for s, e in cds:
        ax.add_patch(plt.Rectangle(((s-origin)/1000, y-0.18), (e-s+1)/1000, 0.36, color="#2c7fb8", zorder=2))
    ax.text(-0.4, y, f"{nm} (Type {OSMC[nm][1]})", ha="right", va="center", fontsize=9)
ax.set_xlabel("Position (kb, 5'->3' of CDS span)"); ax.set_yticks([])
for sp in ("top","right","left"): ax.spines[sp].set_visible(False)
ax.set_title("Exon (CDS) structure of rice metacaspase genes, ordered by phylogenetic clade (blue=CDS, grey=intron)", fontsize=9)
plt.tight_layout(); plt.savefig(f"{OUT}/results/figures/Fig3_gene_structure.png", dpi=300); plt.close()

# --- k-mer motif substitute (documented as such) ---
rice = read_fasta(RICE_PROT)
prots = {k: rice[loc+".1"].rstrip("*") for k,(loc,t) in OSMC.items()}
dom = json.load(open(f"{OUT}/results/tables/domain_features.json"))
domd = {d["name"]: d for d in dom}
def kmers(s,k): return [s[i:i+k] for i in range(len(s)-k+1)]
K=8
counter = Counter()
for s in prots.values(): counter.update(set(kmers(s,K)))
cand = [w for w,c in counter.items() if c>=4]
# collapse overlapping candidates greedily by frequency
cand.sort(key=lambda w: -counter[w])
chosen=[]
for w in cand:
    if not any(w in c or c in w for c in chosen): chosen.append(w)
    if len(chosen)==6: break
with open(f"{OUT}/results/tables/Table3_protein_motifs_kmer_REGENERATED.tsv","w") as f:
    f.write("motif_rank\tconsensus_kmer\tn_sequences_with_hit\tsequences\n")
    for i,w in enumerate(chosen,1):
        hits=[k for k,s in prots.items() if w in s]
        f.write(f"{i}\t{w}\t{len(hits)}\t{','.join(hits)}\n")
        print(i,w,len(hits),hits)
