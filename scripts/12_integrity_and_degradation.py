"""Step 12 - (a) gene-model integrity: ORF, start/stop, canonical GT-AG splice sites for all OsMC genes (OsMC1 vs the rest);
(b) is OsMC1 degraded only at the dyad, or across otherwise invariant sites?  Leave-one-out invariant-site analysis on the
AtMC9-anchored catalytic-core alignment (rice + Arabidopsis, 17 proteins)."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from common import *; from mcdata import *
genes, order = parse_gff(RICE_GFF); chrc = {}
print("(a) GENE-MODEL INTEGRITY")
rows = []
for k, (loc, t) in OSMC.items():
    g = genes[loc]; mid = [m for m, v in g["mrna"].items() if v["longest"]][0]; m = g["mrna"][mid]
    if g["chr"] not in chrc: chrc[g["chr"]] = load_chrom(RICE_GEN, g["chr"])
    cs = chrc[g["chr"]]; cds = sorted(m["cds"])
    sp = []
    for i in range(len(cds) - 1):
        s_int, e_int = cds[i][1] + 1, cds[i+1][0] - 1
        seq = cs[s_int-1:e_int]
        if g["strand"] == "-": seq = revcomp(seq)
        sp.append(seq[:2] + ".." + seq[-2:])
    full = gene_cds(genes, loc, mid, cs); pep = translate(full)
    internal_stops = pep[:-1].count("*")
    rows.append((k, len(cds), full[:3], pep[-1], internal_stops, ",".join(sp) if sp else "intronless"))
    print(f"{k:6s} exons={len(cds)} start={full[:3]} stop={'*' if pep[-1]=='*' else 'NONE'} internal_stops={internal_stops} introns(donor..acceptor)={rows[-1][5]} {'ALL_CANONICAL' if all(x in ('GT..AG','GC..AG') for x in sp) else ''}")
open(f"{OUT}/results/tables/S10_gene_model_integrity.tsv","w").write("gene\texons\tstart\tstop\tinternal_stops\tintron_dinucleotides\n"+"\n".join("\t".join(map(str,r)) for r in rows)+"\n")

print("\n(b) LEAVE-ONE-OUT INVARIANT-SITE DEGRADATION")
aln = read_fasta(f"{OUT}/results/alignments/core_alignment_AtMC9_anchored.fa")
names = list(aln); L = len(next(iter(aln.values())))
def analyse(subset_names, label):
    out = []
    for target in subset_names:
        others = [n for n in subset_names if n != target]
        inv = 0; dev = 0
        for c in range(L):
            col = [aln[n][c] for n in others]
            if "-" in col or aln[target][c] == "-": continue
            if len(set(col)) == 1:
                inv += 1
                if aln[target][c] != col[0]: dev += 1
        out.append((target, inv, dev, dev / inv if inv else float("nan")))
    print(f"-- {label}: columns invariant among the other proteins, and the fraction where the left-out protein deviates")
    for t, inv, dev, f in sorted(out, key=lambda x: -x[3]): print(f"   {t:7s} invariant_cols={inv:3d} deviations={dev:3d} fraction={f:.3f}")
    return out
allr = analyse(names, "all 17 proteins")
typeI = [n for n in names if (OSMC.get(n, (0, ''))[1] == "I") or AT_TYPE.get(n) == "I"]
t1 = analyse(typeI, "Type I only (OsMC1-3, AtMC1-3)")
with open(f"{OUT}/results/tables/S11_invariant_site_deviation.tsv","w") as f:
    f.write("set\tprotein\tinvariant_columns\tdeviations\tfraction\n")
    for lab, res in (("all17", allr), ("typeI", t1)):
        for t, inv, dev, fr in res: f.write(f"{lab}\t{t}\t{inv}\t{dev}\t{fr:.4f}\n")
