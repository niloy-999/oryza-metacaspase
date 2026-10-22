"""Step 6 - independently re-extract the 2 kb upstream (ATG-relative) promoter regions from the O. sativa genome+GFF3
(reproducible from raw data, not copied from the uploaded PlantCARE table) and scan them with a curated library of
literature cis-element consensus motifs (PlantCARE/PLACE core motifs). Separate from the MEME PROTEIN-motif analysis (04b)."""
import sys, os, re
sys.path.insert(0, os.path.dirname(__file__))
from common import *; from mcdata import *
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import numpy as np

# consensus IUPAC motifs, literature-curated (PlantCARE/PLACE), category-tagged
IUPAC = {"R":"[AG]","Y":"[CT]","W":"[AT]","S":"[GC]","K":"[GT]","M":"[AC]","N":"[ACGT]"}
def iupac_re(m):
    return re.compile("".join(IUPAC.get(c, c) for c in m))
ELEMENTS = {
 "TATA-box": ("core promoter", "TATAWAWR"), "CAAT-box": ("core promoter", "CCAAT"),
 "ABRE": ("hormone: ABA", "ACGTG"), "G-box": ("light/hormone", "CACGTG"),
 "MYB": ("stress/dev TF site", "WAACCA"), "MYB-core": ("stress/dev TF site", "CNGTTR"),
 "MYC": ("dehydration/JA", "CATGTG"), "W-box": ("defense (WRKY)", "TTGACY"),
 "STRE": ("general stress", "AGGGG"), "ARE": ("anaerobic induction", "TGGTTT"),
 "LTR": ("low-temperature", "CCGAAA"), "TC-rich repeats": ("defense/stress", "GTTTTCTTAC"),
 "GARE": ("gibberellin", "TAACAAR"), "TGA-element": ("auxin", "AACGAC"),
 "CGTCA-motif": ("MeJA", "CGTCA"), "TCA-element": ("salicylic acid", "CCATCTTTTT"),
 "ERE": ("ethylene", "ATTTCAAA"), "AuxRR-core": ("auxin", "GGTCCAT"),
 "GT1-motif": ("light", "GGTTAA"), "Box 4": ("light", "ATTAAT"),
 "circadian": ("circadian control", "CAANNNNATC"), "WUN-motif": ("wounding", "AAATTTCTT"),
 "O2-site": ("zein metabolism reg.", "GATGATGTGG"), "P-box": ("gibberellin", "CCTTTTG"),
}
CAT_COLORS = {"hormone: ABA":"#8e44ad","MeJA":"#8e44ad","salicylic acid":"#8e44ad","ethylene":"#8e44ad","gibberellin":"#8e44ad","auxin":"#8e44ad",
              "light/hormone":"#f39c12","light":"#f39c12","general stress":"#c0392b","anaerobic induction":"#c0392b","low-temperature":"#c0392b",
              "defense/stress":"#c0392b","defense (WRKY)":"#c0392b","dehydration/JA":"#c0392b","wounding":"#c0392b","stress/dev TF site":"#16a085",
              "core promoter":"#7f8c8d","circadian control":"#2980b9","zein metabolism reg.":"#7f8c8d"}
CRE = {k: iupac_re(v[1]) for k, v in ELEMENTS.items()}

genes, order = parse_gff(RICE_GFF)
chrcache = {}
UP = 2000
proms = {}
for k, (loc, t) in OSMC.items():
    g = genes[loc]
    if g["chr"] not in chrcache: chrcache[g["chr"]] = load_chrom(RICE_GEN, g["chr"])
    cs = chrcache[g["chr"]]
    if g["strand"] == "+":
        s = max(1, g["start"] - UP); seq = cs[s-1:g["start"]-1]
    else:
        e = min(len(cs), g["end"] + UP); seq = revcomp(cs[g["end"]:e])
    proms[k] = seq
write_fasta(proms, f"{OUT}/results/sequences/promoters_2kb_upstream.fa")

counts = {k: {e: 0 for e in ELEMENTS} for k in proms}
for k, seq in proms.items():
    for e, rx in CRE.items():
        counts[k][e] = len(rx.findall(seq)) + len(rx.findall(revcomp(seq)))

with open(f"{OUT}/results/tables/Table4_promoter_cis_elements.tsv", "w") as f:
    f.write("gene\t" + "\t".join(ELEMENTS) + "\n")
    for k in proms: f.write(k + "\t" + "\t".join(str(counts[k][e]) for e in ELEMENTS) + "\n")

mat = np.array([[counts[k][e] for e in ELEMENTS] for k in proms])
fig, ax = plt.subplots(figsize=(10.5, 4.2))
im = ax.imshow(mat, cmap="YlOrRd", aspect="auto")
ax.set_xticks(range(len(ELEMENTS))); ax.set_xticklabels(list(ELEMENTS), rotation=60, ha="right", fontsize=7.5)
ax.set_yticks(range(len(proms))); ax.set_yticklabels(list(proms), fontsize=8)
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        if mat[i,j]>0: ax.text(j, i, int(mat[i,j]), ha="center", va="center", fontsize=6.5, color="black" if mat[i,j]<mat.max()*0.6 else "white")
for j,e in enumerate(ELEMENTS):
    ax.get_xticklabels()[j].set_color(CAT_COLORS.get(ELEMENTS[e][0], "black"))
plt.colorbar(im, ax=ax, label="element count / 2 kb promoter", shrink=0.8)
ax.set_title("Cis-regulatory element counts, 2 kb promoters (re-extracted from genome; PlantCARE/PLACE consensus motifs)", fontsize=9)
plt.tight_layout(); plt.savefig(f"{OUT}/results/figures/Fig7_promoter_cis_elements_heatmap.png", dpi=300); plt.close()
for k in proms: print(k, sum(counts[k].values()), {e:c for e,c in counts[k].items() if c})
