"""Shared helpers for the Oryza/Arabidopsis metacaspase re-analysis (pure Python + numpy; no Biopython needed)."""
import re, os, math
from collections import defaultdict, OrderedDict
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROJ = os.environ.get("MC_PROJ", os.path.join(ROOT, "data", "raw"))   # Phytozome files, see README
OUT  = os.environ.get("MC_OUT", ROOT)
RICE_PROT = f"{PROJ}/O.Sativa/Osativa_323_v7.0.protein.fa"
RICE_GFF  = f"{PROJ}/O.Sativa/Osativa_323_v7.0.gene.gff3"
RICE_GEN  = f"{PROJ}/O.Sativa/Osativa_323_v7.0.fa"
AT_PROT   = f"{PROJ}/Arabidopsis thaliana/Athaliana_167_protein.fa"
AT_GFF    = f"{PROJ}/Arabidopsis thaliana/Athaliana_167_gene.gff3"

AA = "ARNDCQEGHILKMFPSTWYV"
_B62 = """
 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3
-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2
-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4
"""
B62 = np.array([[int(x) for x in l.split()] for l in _B62.strip().splitlines()], dtype=int)
assert B62.shape == (20, 20) and (B62 == B62.T).all(), "BLOSUM62 typo"
AAIDX = {a: i for i, a in enumerate(AA)}

def read_fasta(path, keep=None, keyfn=lambda h: h.split()[0]):
    seqs, name, buf = OrderedDict(), None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None and (keep is None or name in keep): seqs[name] = "".join(buf)
                name = keyfn(line[1:].strip()); buf = []
            else:
                buf.append(line.strip())
        if name is not None and (keep is None or name in keep): seqs[name] = "".join(buf)
    return seqs

def write_fasta(seqs, path, width=60):
    with open(path, "w") as fh:
        for k, v in seqs.items():
            fh.write(f">{k}\n")
            for i in range(0, len(v), width): fh.write(v[i:i+width] + "\n")

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
def revcomp(s): return s.translate(COMP)[::-1]

CODON = {}
_b = "TCAG"; _aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_i = 0
for a in _b:
    for b in _b:
        for c in _b:
            CODON[a+b+c] = _aa[_i]; _i += 1
def translate(cds):
    return "".join(CODON.get(cds[i:i+3], "X") for i in range(0, len(cds) - len(cds) % 3, 3))

def parse_gff(path):
    """gene_id -> dict(chr,start,end,strand,mrna={id:{exon,cds,utr5,utr3,longest}}); plus per-chrom ordered gene list."""
    genes, order, raw2gene, mrna2gene = {}, defaultdict(list), {}, {}
    with open(path) as fh:
        for line in fh:
            if line[0] == "#": continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9: continue
            chrom, src, typ, s, e, _, strand, _, attr = f
            if typ not in ("gene", "mRNA", "exon", "CDS", "five_prime_UTR", "three_prime_UTR"): continue
            a = dict(kv.split("=", 1) for kv in attr.split(";") if "=" in kv)
            if typ == "gene":
                gid = a.get("Name") or a["ID"]
                genes[gid] = dict(chr=chrom, start=int(s), end=int(e), strand=strand, mrna=OrderedDict())
                raw2gene[a["ID"]] = gid
                order[chrom].append((int(s), int(e), gid))
            elif typ == "mRNA":
                pg = raw2gene.get(a["Parent"])
                if pg is None: continue
                mid = a.get("Name") or a["ID"]
                genes[pg]["mrna"][mid] = dict(exon=[], cds=[], utr5=[], utr3=[], longest=a.get("longest", "0") == "1", pacid=a.get("pacid"))
                mrna2gene[a["ID"]] = (pg, mid)
            else:
                par = a.get("Parent")
                if par in mrna2gene:
                    pg, mid = mrna2gene[par]
                    key = dict(exon="exon", CDS="cds", five_prime_UTR="utr5", three_prime_UTR="utr3")[typ]
                    genes[pg]["mrna"][mid][key].append((int(s), int(e)))
    for c in order: order[c].sort()
    return genes, order

def load_chrom(path, name):
    out, on = [], False
    with open(path) as fh:
        for line in fh:
            if line[0] == ">":
                if on: break
                on = (line[1:].split()[0] == name)
            elif on: out.append(line.strip())
    return "".join(out).upper()

def gene_cds(genes, gid, mid, chrseq):
    m = genes[gid]["mrna"][mid]; strand = genes[gid]["strand"]
    segs = sorted(m["cds"], key=lambda x: x[0])          # ascending genomic order
    seq = "".join(chrseq[s-1:e] for s, e in segs)
    return revcomp(seq) if strand == "-" else seq          # revcomp of the whole spliced CDS for minus-strand genes
