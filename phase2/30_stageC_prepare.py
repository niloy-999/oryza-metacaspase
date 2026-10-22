#!/usr/bin/env python3
"""STAGE C preparation (runs in the sandbox). Builds the sequence sets for (1) a family-wide ML tree and (2) a Type I codon alignment for a
selection test, from the Stage B ortholog assignments. Only reciprocal-best-hit orthologs with plausible model length are used; the Nipponbare
sequences of OsMC1-8 come from the MSU models (the RAP-DB/Ensembl models of OsMC2/OsMC3/OsMC8 are fragmentary)."""
import sys, os, csv, collections
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from common import *; from mcdata import *
R = os.environ.get("MC_P2RES", os.path.join(ROOT, "phase2", "results_stageB")); D = os.environ.get("MC_P2DATA", os.path.join(ROOT, "phase2", "phase2_out")); O = os.path.join(ROOT, "phase2", "stageC")
rows = list(csv.DictReader(open(f"{R}/all_candidates_assigned.tsv"), delimiter="\t"))
prot = {}; k = None
for l in open(f"{R}/all_candidates.faa"):
    if l[0] == ">": k = l[1:].strip(); prot[k] = ""
    else: prot[k] += l.strip()
cds = {}
for sp in os.listdir(D):
    p = f"{D}/{sp}/candidates.cds.fna"
    if os.path.exists(p):
        k = None
        for l in open(p):
            if l[0] == ">": t = l[1:].strip().split("|"); k = (t[0], t[1]); cds[k] = ""
            else: cds[k] += l.strip().upper()
LEN = {"OsMC1": (250, 520), "OsMC2": (280, 520), "OsMC3": (280, 520), "OsMC4": (300, 520), "OsMC5": (300, 520), "OsMC6": (300, 520), "OsMC7": (260, 450), "OsMC8": (240, 420), "OsMC9": (250, 600)}
fam = collections.OrderedDict(); t1p = collections.OrderedDict(); t1c = collections.OrderedDict(); drop = collections.Counter(); kept = collections.Counter()
genes, order = parse_gff(RICE_GFF); chrc = {}
def msu_cds(loc):
    g = genes[loc]; mid = [m for m, v in g["mrna"].items() if v["longest"]][0]
    if g["chr"] not in chrc: chrc[g["chr"]] = load_chrom(RICE_GEN, g["chr"])
    return gene_cds(genes, loc, mid, chrc[g["chr"]])
# Nipponbare from MSU for OsMC1-8
rice = read_fasta(RICE_PROT)
for name, (loc, typ) in OSMC.items():
    fam[f"{name}__oryza_sativa_MSU__{loc}"] = rice[loc + ".1"].rstrip("*"); kept[name] += 1
    if name in ("OsMC1", "OsMC2", "OsMC3", "OsMC8"):
        c = msu_cds(loc); t1p[f"{name}__oryza_sativa_MSU__{loc}"] = rice[loc + ".1"].rstrip("*"); t1c[f"{name}__oryza_sativa_MSU__{loc}"] = c[:-3] if translate(c).endswith("*") else c
# Arabidopsis references
at = read_fasta(os.path.join(os.path.dirname(O), "seeds_all_17.fa"))
for k, v in at.items():
    if k.startswith("AtMC"): fam[k.replace("|", "__At__")] = v
for r in rows:
    g = r["best_OsMC"]; sp = r["species"]
    if r["RBH"] != "True" or g not in LEN: continue
    if sp == "oryza_sativa" and g != "OsMC9": continue          # MSU model used instead
    aa = int(r["aa"]); lo, hi = LEN[g]
    if not (lo <= aa <= hi): drop[g] += 1; continue
    key = f"{sp}|{r['gene']}|{g}"; seq = prot[key]
    fam[f"{g}__{sp}__{r['gene']}"] = seq; kept[g] += 1
    if g in ("OsMC1", "OsMC2", "OsMC3", "OsMC8", "OsMC9") and sp != "leersia_perrieri_x":
        c = cds.get((sp, r["gene"]))
        if c and len(c) % 3 == 0:
            tp = translate(c).rstrip("*")
            if tp == seq: t1p[f"{g}__{sp}__{r['gene']}"] = seq; t1c[f"{g}__{sp}__{r['gene']}"] = c[:3*len(seq)]
            else: drop[g + "_cds_mismatch"] += 1
        else: drop[g + "_no_cds"] += 1
os.makedirs(O, exist_ok=True)
write_fasta(fam, f"{O}/family_proteins.faa"); write_fasta(t1p, f"{O}/typeI_proteins.faa"); write_fasta(t1c, f"{O}/typeI_cds.fna")
print("family-wide set:", len(fam), "sequences | kept per group:", dict(kept)); print("dropped (length/cds):", dict(drop))
print("Type I codon set:", len(t1p), dict(collections.Counter(k.split('__')[0] for k in t1p)))
