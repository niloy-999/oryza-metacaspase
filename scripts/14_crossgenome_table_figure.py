"""Step 14 - Table S14 (OsMC1 orthologs in 27 genomes) and Fig. 7 (dyad status by group + OsMC1 residues per genome)."""
import csv, re, collections, os, sys
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT = ROOT
R = os.environ.get("MC_P2RES", os.path.join(ROOT, "phase2", "results_stageB")); D = os.environ.get("MC_P2DATA", os.path.join(ROOT, "phase2", "phase2_out"))
rows = list(csv.DictReader(open(f"{R}/all_candidates_assigned.tsv"), delimiter="\t"))
man = {}
for l in open(f"{D}/manifest.tsv"):
    f = l.rstrip("\n").split("\t")
    if len(f) >= 2:
        m = re.search(r"/pep/[^.]+\.(.+)\.pep\.all", f[1]); man[f[0]] = m.group(1) if m else ""
order = ["oryza_sativa"] + sorted(s for s in {r["species"] for r in rows} if s.startswith("oryza_sativa_")) + ["oryza_indica"] + \
        sorted(s for s in {r["species"] for r in rows} if not s.startswith("oryza_sativa") and s not in ("oryza_indica", "leersia_perrieri")) + ["leersia_perrieri"]
o1 = {r["species"]: r for r in rows if r["best_OsMC"] == "OsMC1" and r["RBH"] == "True"}
with open(f"{OUT}/results/tables/TableS14_OsMC1_orthologs_27_genomes.tsv", "w") as f:
    f.write("genome\tassembly\tOsMC1_ortholog_gene\tprotein_aa\tHis_position_residue\tCys_position_residue\tHis_motif\tCys_motif\tscore_margin_over_runner_up\tcluster_window_genes_within_15\tlocus\n")
    for sp in order:
        r = o1[sp]; f.write(f"{sp}\t{man.get(sp,'')}\t{r['gene']}\t{r['aa']}\t{r['His_res']}{r['His_pos']}\t{r['Cys_res']}{r['Cys_pos']}\t{r['His_motif']}\t{r['Cys_motif']}\t{r['margin']}\t{r['window_genes_within_15'] if 'window_genes_within_15' in r else r['window_genes_within15']}\t{r['locus']}\n")
# status per group
groups = ["OsMC%d" % i for i in range(1, 10)]; st = {g: [0, 0, 0] for g in groups}
for r in rows:
    if r["RBH"] == "True" and r["best_OsMC"] in st:
        if r["His_res"] == "?": st[r["best_OsMC"]][2] += 1
        elif r["His_res"] == "H" and r["Cys_res"] == "C": st[r["best_OsMC"]][0] += 1
        else: st[r["best_OsMC"]][1] += 1
fig, (a, b) = plt.subplots(1, 2, figsize=(11.5, 6.6), gridspec_kw={"width_ratios": [1, 1.25]})
x = range(len(groups))
a.bar(x, [st[g][0] for g in groups], color="#2e7d32", label="both His and Cys present")
a.bar(x, [st[g][1] for g in groups], bottom=[st[g][0] for g in groups], color="#d62728", label="substituted")
a.bar(x, [st[g][2] for g in groups], bottom=[st[g][0] + st[g][1] for g in groups], color="#bdbdbd", label="not assessable (fragmentary gene model)")
a.set_xticks(list(x)); a.set_xticklabels(groups, rotation=45, ha="right"); a.set_ylabel("orthologs (of up to 27 genomes)")
a.set_title("A  Catalytic-dyad status by ortholog group", loc="left", fontsize=10, fontweight="bold"); a.set_ylim(0, 36); a.legend(fontsize=7.5, loc="upper right", frameon=False)
b.axis("off"); b.set_title("B  OsMC1 orthologs: residues at the catalytic positions", loc="left", fontsize=10, fontweight="bold")
n = len(order)
for i, sp in enumerate(order):
    y = n - i; r = o1[sp]
    lab = sp.replace("oryza_", "O. ").replace("sativa_", "sativa ").replace("leersia_perrieri", "Leersia perrieri").replace("O. sativa", "O. sativa").replace("_", " ")
    b.text(0.0, y, lab, ha="left", va="center", fontsize=7.5, style="italic" if sp.startswith(("oryza_", "leersia")) and "sativa " not in lab else "normal")
    for xx, res, ref in ((0.62, r["His_res"], "H"), (0.80, r["Cys_res"], "C")):
        b.text(xx, y, res, ha="center", va="center", fontsize=8, color="white", fontweight="bold", bbox=dict(boxstyle="square,pad=0.25", fc="#d62728" if res != ref else "#2e7d32", ec="none"))
    b.text(0.92, y, f"{r['aa']} aa", ha="left", va="center", fontsize=6.5, color="#555")
b.text(0.62, n + 1.3, "His", ha="center", fontsize=8, fontweight="bold"); b.text(0.80, n + 1.3, "Cys", ha="center", fontsize=8, fontweight="bold")
b.set_xlim(0, 1.1); b.set_ylim(0, n + 2)
plt.tight_layout(); plt.savefig(f"{OUT}/results/figures/Fig8_crossgenome_dyad_status.png", dpi=300); plt.close()
print("ok", len(order), st)
