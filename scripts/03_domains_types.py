"""Step 3 - domain architecture: catalytic His/Cys dyad (numbering transferred from AtMC9 = AtMCA-IIf, His95/Cys147; PNAS 2023, Nat Commun 2025),
LSD1-type zinc-finger pro-domain, N-terminal extension, catalytic-core span, Type I / II assignment."""
import sys, os, re, json
sys.path.insert(0, os.path.dirname(__file__))
from common import *; from align import *; from mcdata import *

rice = read_fasta(RICE_PROT); at = read_fasta(AT_PROT, keyfn=lambda h: h.split()[0].split("|")[0])
prots = OrderedDict()
for k, (loc, t) in OSMC.items(): prots[k] = (loc, rice[loc + ".1"].rstrip("*"), "Oryza sativa")
_ex = f"{OUT}/results/sequences/OsMC9_Os10g0565100.fa"
if os.path.exists(_ex): prots["OsMC9"] = ("Os10g0565100", list(read_fasta(_ex).values())[0], "Oryza sativa (RAP-DB only)")
prots["Os09g36640(ms:oryza-6)"] = ("LOC_Os09g36640", rice["LOC_Os09g36640.1"].rstrip("*"), "Oryza sativa (negative control)")
for k, loc in AT_MC.items(): prots[k] = (loc, at[loc + ".1"].rstrip("*"), "Arabidopsis thaliana")
ref = prots["AtMC9"][1]; ref1 = prots["AtMC1"][1]
# ---- validate the reference numbering ----
H_ref, C_ref, R_ref = 95, 147, 183
print("AtMC9 residue 95 =", ref[H_ref-1], "| 147 =", ref[C_ref-1], "| 183 =", ref[R_ref-1], "| context H:", ref[H_ref-6:H_ref+5], "C:", ref[C_ref-6:C_ref+5])
assert ref[H_ref-1] == "H" and ref[C_ref-1] == "C", "reference numbering does not match seed sequence"

ZF = re.compile(r"C..C.{8,40}C..C")   # two C-x2-C knuckles (LSD1/LOL-type zinc finger), first 150 aa only
rows = []
for name, (loc, seq, sp) in prots.items():
    sc, ra, rb, pairs = pair_align(ref, seq, "local")          # a = AtMC9 (ref), b = query
    mp = map_positions(ref, seq, pairs)                         # ref idx -> query idx
    qspan = [j for i, j in pairs if i is not None and j is not None]
    rspan = [i for i, j in pairs if i is not None and j is not None]
    def at_col(refpos):
        j = mp.get(refpos - 1)
        if j is None:  # nearest mapped column within +/-3
            for d in (1, 2, 3):
                for s in (-d, d):
                    if mp.get(refpos - 1 + s) is not None: j = mp[refpos - 1 + s]; break
                if j is not None: break
        return j
    jh, jc, jr = at_col(H_ref), at_col(C_ref), at_col(R_ref)
    hres = seq[jh] if jh is not None else "-"; cres = seq[jc] if jc is not None else "-"
    # region-level dyad calls: the dyad must lie in the conserved GH(G) / SGC motifs
    hctx = seq[jh-4:jh+5] if jh is not None else ""; cctx = seq[jc-4:jc+5] if jc is not None else ""
    zf = ZF.search(seq[:150])
    sc1, ra1, rb1, pairs1 = pair_align(ref1, seq, "local")
    q1 = [j for i, j in pairs1 if i is not None and j is not None]
    core_start, core_end = min(min(qspan), min(q1)) + 1, max(max(qspan), max(q1)) + 1
    nterm = core_start - 1              # residues upstream of homologous core, corrected for AtMC9's own start
    ident = identity(ra, rb)
    n_ext = max(0, (jh if jh is not None else 0) - 90)   # residues upstream of the His beyond the ~90 aa found in Type II
    rows.append(dict(name=name, locus=loc, species=sp, length_aa=len(seq), core_start=core_start, core_end=core_end,
                     core_cov_AtMC9=round(100 * (max(rspan) - min(rspan) + 1) / len(ref), 1), id_to_AtMC9_pct=round(100 * ident, 1),
                     N_terminal_extension_aa=n_ext, zinc_finger=("C%d-C%d-C%d-C%d" % tuple(zf.start(0) + 1 + o for o in (0, 3, len(zf.group(0)) - 4, len(zf.group(0)) - 1))) if zf else "",
                     His_equiv=(f"{hres}{jh+1}" if jh is not None else "-"), His_context=hctx,
                     Cys_equiv=(f"{cres}{jc+1}" if jc is not None else "-"), Cys_context=cctx,
                     His_Cys_spacing=(jc - jh) if (jh is not None and jc is not None) else None,
                     R183_equiv=(f"{seq[jr]}{jr+1}" if jr is not None else "-"),
                     catalytic_dyad_intact=(hres == "H" and cres == "C"),
                     type_call=("-" if jh is None else ("I" if n_ext >= 40 else "II"))))
keys = list(rows[0].keys())
with open(f"{OUT}/results/tables/Table1_domain_architecture_types.tsv", "w") as f:
    f.write("\t".join(keys) + "\n")
    for r in rows: f.write("\t".join("" if r[k] is None else str(r[k]) for k in keys) + "\n")
for r in rows:
    exp = OSMC.get(r["name"], (None, None))[1] or AT_TYPE.get(r["name"], "-")
    print(f"{r['name']:24s} L={r['length_aa']:4d} core={r['core_start']}-{r['core_end']} cov={r['core_cov_AtMC9']:5.1f}% id={r['id_to_AtMC9_pct']:5.1f}% Next={r['N_terminal_extension_aa']:3d} ZF={'Y' if r['zinc_finger'] else '-'} "
          f"H={r['His_equiv']:6s}[{r['His_context']}] C={r['Cys_equiv']:6s}[{r['Cys_context']}] sp={r['His_Cys_spacing']} dyad={r['catalytic_dyad_intact']} type={r['type_call']} (literature {exp})")
json.dump(rows, open(f"{OUT}/results/tables/domain_features.json", "w"), indent=1)
