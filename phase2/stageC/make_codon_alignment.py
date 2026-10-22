#!/usr/bin/env python3
"""Protein alignment (aligned FASTA) + unaligned CDS (same names, stop codon removed) -> codon alignment (PAL2NAL logic).
Checks that every CDS translates exactly to its aligned protein. Usage: make_codon_alignment.py prot.aln.faa cds.fna out.codon.fna"""
import sys
CODON = {}
b = "TCAG"; aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"; i = 0
for x in b:
    for y in b:
        for z in b: CODON[x+y+z] = aa[i]; i += 1
def read(p):
    d = {}; k = None
    for l in open(p):
        if l[0] == ">": k = l[1:].strip().split()[0]; d[k] = []
        elif k: d[k].append(l.strip())
    return {a: "".join(b) for a, b in d.items()}
prot, cds = read(sys.argv[1]), read(sys.argv[2]); out = {}; bad = 0
for k, p in prot.items():
    c = cds.get(k)
    if c is None: print("no CDS for", k, file=sys.stderr); bad += 1; continue
    res = []; j = 0
    for ch in p:
        if ch == "-": res.append("---")
        else:
            cod = c[3*j:3*j+3]
            if CODON.get(cod, "X") != ch.upper(): print(f"translation mismatch in {k} at residue {j+1}: {ch} vs {CODON.get(cod,'X')}", file=sys.stderr); bad += 1; break
            res.append(cod); j += 1
    else:
        if j * 3 != len(c): print(f"{k}: CDS longer than protein ({len(c)} vs {3*j})", file=sys.stderr)
        out[k] = "".join(res)
with open(sys.argv[3], "w") as f:
    for k, s in out.items(): f.write(f">{k}\n{s}\n")
print(f"codon alignment: {len(out)} sequences, {len(next(iter(out.values())))//3 if out else 0} codons; problems: {bad}")
