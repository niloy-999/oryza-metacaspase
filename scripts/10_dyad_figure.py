"""Step 10 - Fig. 5: catalytic-dyad motif alignment (His-window and Cys-window) across OsMC1-8 and AtMC1-9, taken directly from sequence."""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
from common import *; from mcdata import *
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
rice = read_fasta(RICE_PROT); at = read_fasta(AT_PROT, keyfn=lambda h: h.split()[0].split("|")[0])
dom = {d["name"]: d for d in json.load(open(f"{OUT}/results/tables/domain_features.json"))}
_ex9 = list(read_fasta(f"{OUT}/results/sequences/OsMC9_Os10g0565100.fa").values())[0]
order = [(k, rice[v[0]+".1"].rstrip("*"), v[1], "Os") for k, v in OSMC.items() if v[1]=="I"] + [("OsMC9", _ex9, "I", "Os")] + \
        [(k, at[v+".1"].rstrip("*"), AT_TYPE[k], "At") for k, v in AT_MC.items() if AT_TYPE[k]=="I"] + \
        [(k, rice[v[0]+".1"].rstrip("*"), v[1], "Os") for k, v in OSMC.items() if v[1]=="II"] + \
        [(k, at[v+".1"].rstrip("*"), AT_TYPE[k], "At") for k, v in AT_MC.items() if AT_TYPE[k]=="II"]
W = 8
def window(name, seq, key):
    pos = int(dom[name][key][1:]) - 1
    return seq[pos-W:pos+W+1], pos+1
fig, ax = plt.subplots(figsize=(10.5, 6.6)); ax.axis("off")
y0 = len(order)
x_h, x_c = 0.0, 22.0
ax.text(x_h+W, y0+1.2, "Catalytic His window", ha="center", fontsize=10, fontweight="bold")
ax.text(x_c+W, y0+1.2, "Catalytic Cys window", ha="center", fontsize=10, fontweight="bold")
prev_type = None; off = 0.0
for i, (name, seq, typ, sp) in enumerate(order):
    if typ != prev_type and prev_type is not None: off += 0.7
    y = y0 - i - off
    prev_type = typ
    col_name = "#c0392b" if typ == "I" else "#1f5fa8"
    ax.text(-2.0, y, f"{name} (T{typ})", ha="right", va="center", fontsize=9, color=col_name, fontweight="bold" if sp=="Os" else "normal")
    for x0, key, want in ((x_h, "His_equiv", "H"), (x_c, "Cys_equiv", "C")):
        w, p = window(name, seq, key)
        for j, aa in enumerate(w):
            cat = (j == W)
            bad = cat and aa != want
            ax.text(x0+j, y, aa, ha="center", va="center", family="monospace", fontsize=9.5,
                    color=("white" if cat else "black"), fontweight="bold" if cat else "normal",
                    bbox=(dict(boxstyle="square,pad=0.15", fc=("#d62728" if bad else "#2e7d32"), ec="none") if cat else None))
ax.set_xlim(-12, 41); ax.set_ylim(-1.2, y0+2)
ax.text(20, -1.0, "Green = catalytic residue present (His / Cys); red = substituted. Windows are ±8 residues centred on the position homologous to AtMC9 His95 / Cys147.", ha="center", fontsize=8)
plt.tight_layout(); plt.savefig(f"{OUT}/results/figures/Fig5_catalytic_dyad_alignment.png", dpi=300); plt.close()
print("ok")
