"""Step 8 - publication-style domain-architecture figure (Type I vs Type II) built from the verified positions in Table1."""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
from common import *
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
dom = json.load(open(f"{OUT}/results/tables/domain_features.json")); d = {x["name"]: x for x in dom}
t1, t2 = d["OsMC3"], d["OsMC4"]   # representative Type I / Type II (both catalytic-dyad intact, well-covered)
fig, ax = plt.subplots(figsize=(8, 2.6))
def draw(rec, y, label):
    L = rec["length_aa"]; ax.add_patch(plt.Rectangle((0, y-0.15), L, 0.3, color="#eee", ec="k", lw=0.8, zorder=1))
    if rec["N_terminal_extension_aa"] > 10:
        ax.add_patch(plt.Rectangle((0, y-0.15), rec["N_terminal_extension_aa"], 0.3, color="#f1c40f", ec="k", lw=0.6, zorder=2))
        ax.text(rec["N_terminal_extension_aa"]/2, y+0.32, "N-terminal\npro-domain", ha="center", fontsize=7)
    p20s, p20e = rec["core_start"], rec["core_start"] + int(0.55*(rec["core_end"]-rec["core_start"]))
    p10s, p10e = p20e+8, rec["core_end"]
    ax.add_patch(plt.Rectangle((p20s, y-0.15), p20e-p20s, 0.3, color="#2c7fb8", ec="k", lw=0.6, zorder=2))
    ax.add_patch(plt.Rectangle((p10s, y-0.15), p10e-p10s, 0.3, color="#41ab5d", ec="k", lw=0.6, zorder=2))
    ax.text((p20s+p20e)/2, y, "p20", ha="center", va="center", fontsize=7, color="white", fontweight="bold")
    ax.text((p10s+p10e)/2, y, "p10", ha="center", va="center", fontsize=7, color="white", fontweight="bold")
    hj = int(rec["His_equiv"][1:]); cj = int(rec["Cys_equiv"][1:])
    ax.plot([hj], [y+0.2], "v", color="#c0392b", ms=6, zorder=3); ax.text(hj, y+0.34, "His", ha="center", fontsize=6.5, color="#c0392b")
    ax.plot([cj], [y+0.2], "v", color="#7d3c98", ms=6, zorder=3); ax.text(cj, y+0.34, "Cys", ha="center", fontsize=6.5, color="#7d3c98")
    ax.text(-10, y, f"{label}\n({L} aa)", ha="right", va="center", fontsize=8.5, fontweight="bold")
draw(t1, 1.0, "OsMC3\nType I"); draw(t2, 0.0, "OsMC4\nType II")
ax.set_xlim(-90, 430); ax.set_ylim(-0.6, 1.6); ax.axis("off")
ax.set_title("Domain architecture of representative Type I and Type II rice metacaspases", fontsize=10)
plt.tight_layout(); plt.savefig(f"{OUT}/results/figures/Fig1_domain_architecture.png", dpi=300); plt.close()
print("done")
