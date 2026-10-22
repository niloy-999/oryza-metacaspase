"""Step 2 - verify every number in the manuscript's Table 2 against the annotation (GFF3), the genome and the proteome,
and recompute physicochemical properties.  Writes Table 2 (corrected) and a discrepancy log."""
import sys, os, csv, math
sys.path.insert(0, os.path.dirname(__file__))
from common import *
from mcdata import *

MASS = dict(A=71.0788,R=156.1875,N=114.1038,D=115.0886,C=103.1388,E=129.1155,Q=128.1307,G=57.0519,H=137.1411,I=113.1594,
            L=113.1594,K=128.1741,M=131.1926,F=147.1766,P=97.1167,S=87.0782,T=101.1051,W=186.2132,Y=163.1760,V=99.1326)
KD = dict(A=1.8,R=-4.5,N=-3.5,D=-3.5,C=2.5,Q=-3.5,E=-3.5,G=-0.4,H=-3.2,I=4.5,L=3.8,K=-3.9,M=1.9,F=2.8,P=-1.6,S=-0.8,T=-0.7,W=-0.9,Y=-1.3,V=4.2)
def mw(p): return (sum(MASS[a] for a in p if a in MASS) + 18.01524) / 1000
def gravy(p): return sum(KD[a] for a in p if a in KD) / len(p)
def aliphatic(p):
    n = len(p); f = lambda a: 100 * p.count(a) / n
    return f("A") + 2.9 * f("V") + 3.9 * (f("I") + f("L"))
def pI(p):
    pos = {"K": 10.0, "R": 12.0, "H": 5.98}; neg = {"D": 4.05, "E": 4.45, "C": 9.0, "Y": 10.0}
    nterm = {"A":7.59,"M":7.0,"S":6.93,"P":8.36,"T":6.82,"V":7.44,"E":7.7}.get(p[0], 9.0)
    cterm = {"D":4.55,"E":4.75}.get(p[-1], 2.0)
    cnt = {a: p.count(a) for a in "KRHDECY"}
    def charge(ph):
        c = 1 / (1 + 10 ** (ph - nterm)) - 1 / (1 + 10 ** (cterm - ph))
        for a, pk in pos.items(): c += cnt[a] / (1 + 10 ** (ph - pk))
        for a, pk in neg.items(): c -= cnt[a] / (1 + 10 ** (pk - ph))
        return c
    lo, hi = 0.0, 14.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if charge(mid) > 0: lo = mid
        else: hi = mid
    return (lo + hi) / 2

if __name__ == "__main__":
    genes, order = parse_gff(RICE_GFF)
    prot = read_fasta(RICE_PROT)
    user_cds = read_fasta(f"{PROJ}/CDS_seq.txt", keyfn=lambda h: h.split()[0])
    chrcache = {}
    rows, log = [], []
    loci = list(MS_NAME.items())
    for ms, loc in loci:
        g = genes[loc]; mids = [m for m, v in g["mrna"].items() if v["longest"]] or list(g["mrna"])
        mid = mids[0]; m = g["mrna"][mid]
        if g["chr"] not in chrcache: chrcache[g["chr"]] = load_chrom(RICE_GEN, g["chr"])
        cds = gene_cds(genes, loc, mid, chrcache[g["chr"]])
        pep = translate(cds)
        pep_nostop = pep.rstrip("*")
        prot_seq = prot[mid].rstrip("*")
        ok_prot = pep_nostop == prot_seq
        ok_start = cds.startswith("ATG"); ok_stop = pep.endswith("*")
        uc = user_cds.get(loc, "")
        ok_user = (uc == cds) if uc else None
        ncds = len(m["cds"]); nex = len(m["exon"]); nint = ncds - 1
        u5 = sum(e - s + 1 for s, e in m["utr5"]); u3 = sum(e - s + 1 for s, e in m["utr3"])
        exp = MS_TABLE2[ms]
        name = LOCUS2OSMC.get(loc, "NOT_A_METACASPASE")
        row = dict(manuscript_name=ms, OsMC_name=name, locus=loc, transcript=mid, chromosome=g["chr"], start_1based=g["start"], end=g["end"],
                   strand=g["strand"], gene_len_bp=g["end"] - g["start"] + 1, CDS_nt=len(cds), protein_aa=len(prot_seq), exons=ncds, introns=nint,
                   utr5_nt=u5, utr3_nt=u3, MW_kDa=round(mw(prot_seq), 2), pI=round(pI(prot_seq), 2), GRAVY=round(gravy(prot_seq), 3),
                   aliphatic_index=round(aliphatic(prot_seq), 1), CDS_translates_to_proteome=ok_prot, has_ATG=ok_start, has_stop=ok_stop,
                   matches_user_CDS_txt=ok_user)
        rows.append(row)
        # discrepancy log vs manuscript
        d = []
        if exp[0] != int(g["chr"].replace("Chr", "")): d.append(f"chromosome {exp[0]} vs {g['chr']}")
        if exp[1] != g["start"]: d.append(f"start {exp[1]} vs {g['start']} (GFF3 is 1-based; off by {g['start']-exp[1]})")
        if exp[2] != g["end"]: d.append(f"end {exp[2]} vs {g['end']}")
        if exp[3] != g["strand"]: d.append(f"strand {exp[3]} vs {g['strand']}")
        if exp[4] != len(cds): d.append(f"CDS {exp[4]} vs {len(cds)}")
        if exp[5] != len(prot_seq): d.append(f"'protein length' {exp[5]} counts the stop codon; actual protein = {len(prot_seq)} aa")
        if abs(exp[6] - mw(prot_seq)) > 0.05: d.append(f"MW {exp[6]} vs recomputed {mw(prot_seq):.2f}")
        if abs(exp[7] - pI(prot_seq)) > 0.15: d.append(f"pI {exp[7]} vs recomputed {pI(prot_seq):.2f}")
        if abs(exp[8] - gravy(prot_seq)) > 0.01: d.append(f"GRAVY {exp[8]} vs recomputed {gravy(prot_seq):.3f}")
        if exp[9] != f"{nint}:{nex:02d}" and exp[9] != f"{nint}:{ncds:02d}": d.append(f"intron:exon '{exp[9]}' vs {nint} introns / {ncds} exons")
        log.append((ms, name, loc, d))
    keys = list(rows[0].keys())
    with open(f"{OUT}/results/tables/Table2_gene_models_physicochemical_VERIFIED.tsv", "w") as f:
        f.write("\t".join(keys) + "\n")
        for r in rows: f.write("\t".join(str(r[k]) for k in keys) + "\n")
    with open(f"{OUT}/results/tables/S2_table2_discrepancies_vs_manuscript.txt", "w") as f:
        for ms, name, loc, d in log:
            f.write(f"{ms} ({name}, {loc}): " + ("; ".join(d) if d else "all values reproduce") + "\n")
    for r in rows: print({k: r[k] for k in ("manuscript_name","OsMC_name","locus","chromosome","start_1based","end","strand","CDS_nt","protein_aa","exons","MW_kDa","pI","GRAVY","CDS_translates_to_proteome","matches_user_CDS_txt")})
    print(); print(open(f"{OUT}/results/tables/S2_table2_discrepancies_vs_manuscript.txt").read())
