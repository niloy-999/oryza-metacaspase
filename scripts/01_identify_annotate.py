"""Step 1 - genome-wide re-identification of metacaspases (rice v7.0 and Arabidopsis TAIR10/Phytozome 167) with a
Smith-Waterman search (BLOSUM62, gap 11/1) seeded by the nine experimentally characterised Arabidopsis metacaspases,
followed by catalytic-dyad / domain / type annotation.  Positive control: the same search must return exactly AtMC1-9 in Arabidopsis."""
import sys, math, re, csv, os
sys.path.insert(0, os.path.dirname(__file__))
from common import *
from align import pair_align, map_positions, identity

LAMBDA, K = 0.267, 0.041   # Karlin-Altschul parameters, BLOSUM62 gap 11/1 (BLAST defaults)
def bits(s): return (LAMBDA * s - math.log(K)) / math.log(2)
def evalue(s, m, n): return m * n * 2 ** (-bits(s))

def load_hits(path, dbsize):
    best = {}
    for l in open(path):
        q, t, s, m, n, qe, te = l.rstrip().split("\t"); s = int(s); m = int(m)
        locus = t.rsplit(".", 1)[0]
        e = evalue(s, m, dbsize)
        if locus not in best or e < best[locus]["E"]:
            best[locus] = dict(query=q.split("|")[0], transcript=t, raw=s, bits=bits(s), E=e, tlen=int(n))
    return best

rice = read_fasta(RICE_PROT); at = read_fasta(AT_PROT, keyfn=lambda h: h.split()[0].split("|")[0])
rice_n = sum(len(v.rstrip("*")) for v in rice.values()); at_n = sum(len(v.rstrip("*")) for v in at.values())
rh = load_hits(os.environ.get("MC_SWDIR","/tmp") + "/sw_rice.tsv", rice_n); ah = load_hits(os.environ.get("MC_SWDIR","/tmp") + "/sw_at.tsv", at_n)

def table(hits, label, thr=1e-3):
    rows = sorted(hits.items(), key=lambda kv: kv[1]["E"])
    out = [(k, v) for k, v in rows if v["E"] < thr]
    print(f"\n{label}: loci with E<{thr}: {len(out)}")
    for k, v in out[:25]: print(f"  {k:22s} {v['query']:6s} bits={v['bits']:7.1f} E={v['E']:.1e} len={v['tlen']}")
    return out
ro = table(rh, "RICE"); ao = table(ah, "ARABIDOPSIS (positive control)")

# significance gap: report the E-value ladder around the cut-off
def ladder(out): return [f"{v['E']:.0e}" for k, v in out[:14]]
print("rice E ladder:", ladder(ro)); print("at   E ladder:", ladder(ao))
with open(f"{OUT}/results/tables/S1_search_hits_rice.tsv", "w") as f:
    f.write("locus\tbest_query\ttranscript\traw_score\tbits\tE_value\ttarget_len\n")
    for k, v in sorted(rh.items(), key=lambda kv: kv[1]["E"])[:60]:
        f.write(f"{k}\t{v['query']}\t{v['transcript']}\t{v['raw']}\t{v['bits']:.1f}\t{v['E']:.2e}\t{v['tlen']}\n")
with open(f"{OUT}/results/tables/S1_search_hits_arabidopsis_positive_control.tsv", "w") as f:
    f.write("locus\tbest_query\ttranscript\traw_score\tbits\tE_value\ttarget_len\n")
    for k, v in sorted(ah.items(), key=lambda kv: kv[1]["E"])[:40]:
        f.write(f"{k}\t{v['query']}\t{v['transcript']}\t{v['raw']}\t{v['bits']:.1f}\t{v['E']:.2e}\t{v['tlen']}\n")
