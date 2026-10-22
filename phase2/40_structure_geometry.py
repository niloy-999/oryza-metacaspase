#!/usr/bin/env python3
"""
PHASE 2, structure comparison. Input: predicted (AlphaFold DB .pdb, AlphaFold Server .cif) or experimental models.
For each model you give the file, a label, and the residue numbers of the catalytic His-equivalent and Cys-equivalent.
Reports, per model: residue identities, confidence (pLDDT from the B-factor column of AlphaFold PDBs) at the dyad, the His-Cys
CA-CA and side-chain heavy-atom distances, and the local pocket superposition RMSD of every model against the first (reference).
USAGE:  python3 40_structure_geometry.py REF.pdb:OsMC2:173:229  OsMC1.pdb:OsMC1:173:229  OsMC3.pdb:OsMC3:195:251 [AtMC9.pdb:AtMC9:95:147]
NOTE: a predicted model cannot prove inactivity; distances in low-confidence regions (pLDDT < 70) must not be over-interpreted.
"""
import sys, math, re
import numpy as np
AA3 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
def load(path):
    res = {}
    if path.lower().endswith((".cif", ".mmcif")):
        cols = []; inloop = False
        for l in open(path):
            if l.startswith("_atom_site."): cols.append(l.strip().split(".")[1]); inloop = True; continue
            if inloop and (l.startswith("#") or l.startswith("loop_") or l.startswith("_")) and cols and not l.startswith("_atom_site."): break
            if inloop and l.startswith(("ATOM", "HETATM")):
                f = l.split(); d = dict(zip(cols, f))
                if d.get("group_PDB") != "ATOM": continue
                key = (d.get("auth_asym_id", d.get("label_asym_id")), int(d.get("auth_seq_id", d.get("label_seq_id"))))
                r = res.setdefault(key, dict(name=d["label_comp_id"], atoms={}, b=None))
                r["atoms"][d["label_atom_id"]] = np.array([float(d["Cartn_x"]), float(d["Cartn_y"]), float(d["Cartn_z"])])
                if d["label_atom_id"] == "CA": r["b"] = float(d.get("B_iso_or_equiv", "nan"))
    else:
        for l in open(path):
            if l.startswith("ATOM"):
                key = (l[21], int(l[22:26])); r = res.setdefault(key, dict(name=l[17:20].strip(), atoms={}, b=None))
                nm = l[12:16].strip(); r["atoms"][nm] = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
                if nm == "CA": r["b"] = float(l[60:66])
    return res
def kabsch(P, Q):
    Pc, Qc = P - P.mean(0), Q - Q.mean(0); H = Pc.T @ Qc; U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T)); Dm = np.diag([1, 1, d]); R = Vt.T @ Dm @ U.T
    return math.sqrt(((Pc @ R.T - Qc) ** 2).sum() / len(P))
def seq_of(res, chain): return "".join(AA3.get(res[k]["name"], "X") for k in sorted(k for k in res if k[0] == chain))
def main():
    models = []
    for a in sys.argv[1:]:
        p, lab, h, c = a.split(":"); res = load(p); chain = sorted({k[0] for k in res})[0]
        models.append(dict(path=p, label=lab, res=res, chain=chain, h=int(h), c=int(c)))
    ref = models[0]
    print(f"{'model':8s} {'His-pos':>8s} {'Cys-pos':>8s} {'pLDDT_H':>8s} {'pLDDT_C':>8s} {'CA-CA(H,C)':>10s} {'minSCdist':>9s} {'pocketRMSD':>10s} {'n_pocket':>8s}")
    for m in models:
        r, ch = m["res"], m["chain"]; H, C = r[(ch, m["h"])], r[(ch, m["c"])]
        dca = float(np.linalg.norm(H["atoms"]["CA"] - C["atoms"]["CA"]))
        sc = lambda R: [v for k, v in R["atoms"].items() if k not in ("N", "C", "O", "CA")]
        mn = min(float(np.linalg.norm(a - b)) for a in sc(H) for b in sc(C)) if sc(H) and sc(C) else float("nan")
        # pocket = residues of the reference with any atom within 8 A of the reference His/Cys side chains
        rr, rch = ref["res"], ref["chain"]
        cen = [v for k in (ref["h"], ref["c"]) for v in rr[(rch, k)]["atoms"].values()]
        pocket = [k for k in sorted(kk for kk in rr if kk[0] == rch) if any(np.linalg.norm(a - b) <= 8.0 for a in rr[k]["atoms"].values() for b in cen)]
        # map reference residue numbers onto this model by sequence alignment (simple global NW with identity scoring)
        sa, sb = seq_of(rr, rch), seq_of(r, ch)
        ka = sorted(k for k in rr if k[0] == rch); kb = sorted(k for k in r if k[0] == ch)
        n, mm = len(sa), len(sb); S = np.zeros((n+1, mm+1)); T = np.zeros((n+1, mm+1), int)
        for i in range(1, n+1): S[i, 0] = -2*i; T[i, 0] = 1
        for j in range(1, mm+1): S[0, j] = -2*j; T[0, j] = 2
        for i in range(1, n+1):
            for j in range(1, mm+1):
                opts = (S[i-1, j-1] + (2 if sa[i-1] == sb[j-1] else -1), S[i-1, j] - 2, S[i, j-1] - 2); t = int(np.argmax(opts)); S[i, j] = opts[t]; T[i, j] = t
        i, j = n, mm; mp = {}
        while i > 0 or j > 0:
            t = T[i, j] if i > 0 and j > 0 else (1 if i > 0 else 2)
            if t == 0: mp[ka[i-1]] = kb[j-1]; i -= 1; j -= 1
            elif t == 1: i -= 1
            else: j -= 1
        pairs = [(k, mp[k]) for k in pocket if k in mp]
        if len(pairs) >= 4:
            P = np.array([rr[a]["atoms"]["CA"] for a, b in pairs]); Q = np.array([r[b]["atoms"]["CA"] for a, b in pairs]); rmsd = kabsch(P, Q)
        else: rmsd = float("nan")
        print(f"{m['label']:8s} {AA3[H['name']]+str(m['h']):>8s} {AA3[C['name']]+str(m['c']):>8s} {H['b'] or float('nan'):8.1f} {C['b'] or float('nan'):8.1f} {dca:10.2f} {mn:9.2f} {rmsd:10.2f} {len(pairs):8d}")
if __name__ == "__main__": main()
