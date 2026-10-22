#!/usr/bin/env python3
"""
PHASE 2, STAGE A  -  run this on YOUR computer (needs internet + NCBI BLAST+ or DIAMOND).

What it does, per Oryza/Leersia genome found on the Gramene / Ensembl Plants FTP:
  1. downloads the protein (pep.all) and CDS (cds.all) FASTA files,
  2. keeps the longest protein per gene and writes a gene-order table (for synteny),
  3. searches the proteome with the 17 seed metacaspases (AtMC1-9, OsMC1-8) and keeps every candidate (E <= 1e-20),
  4. (optional) runs HMMER against Pfam PF00656 if hmmsearch is installed,
  5. BLASTs the ~60 genes around the chromosome-3 OsMC cluster of Nipponbare against each genome (synteny evidence),
  6. writes small per-species result files and a manifest, then zips them (phase2_upload.zip, usually < 60 MB).

USAGE
  conda install -c bioconda -c conda-forge blast diamond hmmer        # blast OR diamond is enough; hmmer optional
  python3 01_fetch_and_extract.py                  # all Oryza/Leersia genomes found
  python3 01_fetch_and_extract.py --max-species 4  # quick test first
  python3 01_fetch_and_extract.py --local-dir DIR  # DIR/<species>/pep.fa.gz and DIR/<species>/cds.fa.gz placed by hand
Re-running skips species that are already finished.
"""
import argparse, gzip, io, os, re, shutil, ssl, subprocess, sys, tempfile, time, zipfile, urllib.request, urllib.error
from collections import OrderedDict, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
BASES = [  # tried in this order; each is a directory that contains one sub-folder per species
    "https://ftp.gramene.org/release-69/fasta/",
    "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-69/fasta/",
    "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-63/fasta/",
    "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta/",
]
WANT = re.compile(r"^(oryza|leersia)", re.I)   # overridden by --want
E_CUTOFF = 1e-20

try:
    import certifi
    CTX = ssl.create_default_context(cafile=certifi.where())
except Exception:
    CTX = ssl.create_default_context()

def http_get(url, binary=False, tries=3):
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"}), timeout=120, context=CTX) as r:
                d = r.read(); return d if binary else d.decode("utf-8", "replace")
        except Exception as e:
            err = e; time.sleep(2 * (i + 1))
    raise RuntimeError(f"{url}: {err}")

def list_dir(url):
    html = http_get(url)
    return sorted(set(re.findall(r'href="([^"?/][^"?]*)"', html)))

def download(url, dest, tries=6):
    if os.path.exists(dest) and os.path.getsize(dest) > 0: return
    err = None
    for i in range(tries):
        try:
            print("   downloading", url, ("(retry %d)" % i if i else ""))
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=180, context=CTX) as r, open(dest + ".part", "wb") as f: shutil.copyfileobj(r, f, 1 << 20)
            os.replace(dest + ".part", dest); return
        except Exception as e:
            err = e; print("   network problem:", str(e)[:80], "- waiting", 15 * (i + 1), "s"); time.sleep(15 * (i + 1))
            if os.path.exists(dest + ".part"): os.remove(dest + ".part")
    raise RuntimeError(f"download failed after {tries} tries: {err}")

def discover():
    """returns {species: (pep_url, cds_url, base)}"""
    found = OrderedDict()
    for base in BASES:
        try: items = list_dir(base)
        except Exception as e: print("  (cannot read", base, "->", str(e)[:80] + ")"); continue
        sp_dirs = [x.strip("/") for x in items if x.endswith("/") and WANT.match(x.strip("/"))]
        print(f"  {base}: {len(sp_dirs)} Oryza/Leersia folders")
        for sp in sp_dirs:
            if sp in found: continue
            try:
                pep = [x for x in list_dir(f"{base}{sp}/pep/") if x.endswith("pep.all.fa.gz")]
                cds = [x for x in list_dir(f"{base}{sp}/cds/") if x.endswith("cds.all.fa.gz")]
            except Exception: continue
            if pep and cds: found[sp] = (f"{base}{sp}/pep/{pep[0]}", f"{base}{sp}/cds/{cds[0]}", base)
    return found

HDR = re.compile(r"(?:pep|cds)\s+(\w+):([^:\s]+):([^:\s]+):(\d+):(\d+):(-?1)")
def parse_fasta(path):
    op = gzip.open if path.endswith(".gz") else open
    name = hdr = None; buf = []
    with op(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name: yield name, hdr, "".join(buf)
                hdr = line[1:].strip(); name = hdr.split()[0]; buf = []
            else: buf.append(line.strip())
    if name: yield name, hdr, "".join(buf)

def load_proteome(path):
    genes = {}
    for pid, hdr, seq in parse_fasta(path):
        seq = seq.rstrip("*")
        m = HDR.search(hdr); g = re.search(r"gene:(\S+)", hdr); t = re.search(r"transcript:(\S+)", hdr)
        gene = g.group(1) if g else pid; tid = t.group(1) if t else pid
        rec = dict(gene=gene, pid=pid, tid=tid, seq=seq, chr=(m.group(3) if m else "NA"), start=int(m.group(4)) if m else 0,
                   end=int(m.group(5)) if m else 0, strand=(m.group(6) if m else "0"), hdr=hdr)
        if gene not in genes or len(seq) > len(genes[gene]["seq"]): genes[gene] = rec
    return genes

def write_fa(path, recs):
    with open(path, "w") as f:
        for k, s in recs:
            f.write(f">{k}\n")
            for i in range(0, len(s), 60): f.write(s[i:i+60] + "\n")

def which(x): return shutil.which(x)

def search(engine, query_fa, subject_fa, out_tsv, workdir, evalue=1e-5):
    """Runs BLAST+/DIAMOND inside a fresh temp directory (BLAST+ breaks on folder names that contain spaces)."""
    tmp = tempfile.mkdtemp(prefix="oryza_blast_")
    try:
        shutil.copy(query_fa, os.path.join(tmp, "q.fa")); shutil.copy(subject_fa, os.path.join(tmp, "s.fa"))
        env = dict(os.environ); env["BLASTDB"] = tmp
        if engine == "blast":
            subprocess.run(["makeblastdb", "-in", "s.fa", "-dbtype", "prot", "-out", "db"], check=True, cwd=tmp, env=env, stdout=subprocess.DEVNULL)
            subprocess.run(["blastp", "-query", "q.fa", "-db", "db", "-evalue", str(evalue), "-outfmt", "6 qseqid sseqid evalue bitscore pident length qlen slen",
                            "-max_target_seqs", "500", "-out", "out.tsv", "-num_threads", "4"], check=True, cwd=tmp, env=env)
        else:
            subprocess.run(["diamond", "makedb", "--in", "s.fa", "-d", "db", "--quiet"], check=True, cwd=tmp, env=env)
            subprocess.run(["diamond", "blastp", "-q", "q.fa", "-d", "db", "-e", str(evalue), "--sensitive", "-k", "500", "--quiet",
                            "--outfmt", "6", "qseqid", "sseqid", "evalue", "bitscore", "pident", "length", "qlen", "slen", "-o", "out.tsv", "-p", "4"], check=True, cwd=tmp, env=env)
        shutil.copy(os.path.join(tmp, "out.tsv"), out_tsv)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="phase2_out"); ap.add_argument("--max-species", type=int, default=0)
    ap.add_argument("--local-dir", default=None); ap.add_argument("--keep-downloads", action="store_true")
    ap.add_argument("--species", default="", help="comma list to restrict")
    ap.add_argument("--want", default="", help="regex of Ensembl folder names to fetch instead of Oryza/Leersia, e.g. \"brachypodium_distachyon|sorghum_bicolor|zea_mays|hordeum_vulgare|setaria_italica|oropetium_thomaeum\" (oryza_sativa is always added as the synteny reference)")
    a = ap.parse_args()
    engine = "blast" if which("blastp") and which("makeblastdb") else ("diamond" if which("diamond") else None)
    if not engine: sys.exit("Install BLAST+ or DIAMOND first:  conda install -c bioconda -c conda-forge blast diamond")
    print("search engine:", engine, "| hmmsearch:", bool(which("hmmsearch")))
    os.makedirs(a.out, exist_ok=True); dl = os.path.join(a.out, "_downloads"); os.makedirs(dl, exist_ok=True)
    seeds = os.path.join(HERE, "seeds_all.fa")
    plan = OrderedDict()
    if a.local_dir:
        for sp in sorted(os.listdir(a.local_dir)):
            p, c = os.path.join(a.local_dir, sp, "pep.fa.gz"), os.path.join(a.local_dir, sp, "cds.fa.gz")
            if os.path.exists(p): plan[sp] = (p, c if os.path.exists(c) else None, "local")
    else:
        global WANT
        if a.want: WANT = re.compile(r"^(" + a.want + r")$|^oryza_sativa$", re.I)
        print("Discovering genomes ..."); found = discover()
        if not found: sys.exit("No Oryza folders could be read from the FTP servers. Use --local-dir with manually downloaded files (see README).")
        plan = found
    if a.species: plan = OrderedDict((k, v) for k, v in plan.items() if k in a.species.split(","))
    if "oryza_sativa" in plan: plan = OrderedDict([("oryza_sativa", plan["oryza_sativa"])] + [(k, v) for k, v in plan.items() if k != "oryza_sativa"])
    if a.max_species: plan = OrderedDict(list(plan.items())[:a.max_species])
    if a.max_species and "oryza_sativa" in plan and list(plan)[0] != "oryza_sativa": pass
    print(len(plan), "genomes:", ", ".join(plan))
    hmm = os.path.join(a.out, "PF00656.hmm")
    if which("hmmsearch") and not os.path.exists(hmm):
        try:
            raw = http_get("https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/PF00656?annotation=hmm", binary=True)
            open(hmm, "wb").write(gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw)
        except Exception as e: print("  could not fetch Pfam HMM (skipping HMMER):", str(e)[:80])
    rice_window = None
    manifest = []
    for sp, (pep_u, cds_u, base) in plan.items():
        sdir = os.path.join(a.out, sp); done = os.path.join(sdir, "DONE")
        if os.path.exists(done): print("skip (done):", sp); continue
        os.makedirs(sdir, exist_ok=True); print("==", sp)
        try:
            if a.local_dir: pep, cds = pep_u, cds_u
            else:
                pep = os.path.join(dl, sp + ".pep.fa.gz"); cds = os.path.join(dl, sp + ".cds.fa.gz")
                download(pep_u, pep); download(cds_u, cds)
            genes = load_proteome(pep)
            order = sorted(genes.values(), key=lambda r: (r["chr"], r["start"]))
            with gzip.open(os.path.join(sdir, "gene_order.tsv.gz"), "wt") as f:
                f.write("idx\tgene\tprotein\ttranscript\tchr\tstart\tend\tstrand\taa\n")
                for i, r in enumerate(order): f.write(f"{i}\t{r['gene']}\t{r['pid']}\t{r['tid']}\t{r['chr']}\t{r['start']}\t{r['end']}\t{r['strand']}\t{len(r['seq'])}\n")
            subj = os.path.join(sdir, "_proteome.faa"); write_fa(subj, [(r["pid"], r["seq"]) for r in genes.values()])
            hits = os.path.join(sdir, "seed_hits.tsv"); search(engine, seeds, subj, hits, sdir)
            cand = defaultdict(lambda: [1.0, 0.0, ""])
            for l in open(hits):
                q, s, e, bs, *_ = l.split("\t"); e = float(e); bs = float(bs)
                if e <= E_CUTOFF and e < cand[s][0]: cand[s] = [e, bs, q]
            tag = {s: "blast" for s in cand}
            if os.path.exists(hmm) and which("hmmsearch"):
                dt = os.path.join(sdir, "pfam.domtbl")
                subprocess.run(["hmmsearch", "--cut_ga", "--domtblout", dt, "-o", os.devnull, hmm, subj], check=False)
                if os.path.exists(dt):
                    for l in open(dt):
                        if l[0] == "#": continue
                        f = l.split(); pid = f[0]
                        if pid not in cand: cand[pid] = [float(f[6]), float(f[7]), "PF00656"]; tag[pid] = "hmmer_only"
                        else: tag[pid] = "blast+hmmer"
            bypid = {r["pid"]: r for r in genes.values()}
            recs = [(f"{sp}|{bypid[p]['gene']}|{p}|{bypid[p]['chr']}:{bypid[p]['start']}-{bypid[p]['end']}({bypid[p]['strand']})|{tag.get(p,'')}", bypid[p]["seq"]) for p in cand if p in bypid]
            write_fa(os.path.join(sdir, "candidates.faa"), recs)
            with open(os.path.join(sdir, "candidates.tsv"), "w") as f:
                f.write("protein\tgene\ttranscript\tchr\tstart\tend\tstrand\taa\tbest_E\tbest_bits\tbest_seed\tevidence\n")
                for p, (e, bs, q) in cand.items():
                    if p in bypid:
                        r = bypid[p]; f.write(f"{p}\t{r['gene']}\t{r['tid']}\t{r['chr']}\t{r['start']}\t{r['end']}\t{r['strand']}\t{len(r['seq'])}\t{e:.2e}\t{bs:.1f}\t{q}\t{tag.get(p,'')}\n")
            if cds:
                want = {bypid[p]["tid"]: p for p in cand if p in bypid}; out = []
                for nm, hdr, seq in parse_fasta(cds):
                    if nm in want: out.append((f"{sp}|{bypid[want[nm]]['gene']}|{want[nm]}", seq))
                write_fa(os.path.join(sdir, "candidates.cds.fna"), out)
            # synteny window: Nipponbare chr3 cluster genes (+-30 genes) searched against this genome
            if rice_window is None and sp == "oryza_sativa":
                ids = [i for i, r in enumerate(order) if r["gene"] in ("Os03g0388900", "Os03g0389000", "Os03g0389100", "Os03g0389400")]
                if ids:
                    lo, hi = max(0, min(ids) - 30), min(len(order) - 1, max(ids) + 30)
                    rice_window = [(f"{order[i]['gene']}", order[i]["seq"]) for i in range(lo, hi + 1)]
                    write_fa(os.path.join(a.out, "nipponbare_chr3_window.faa"), rice_window); print("   rice window:", len(rice_window), "genes")
            if rice_window:
                win = os.path.join(a.out, "nipponbare_chr3_window.faa"); wh = os.path.join(sdir, "window_hits_raw.tsv")
                search(engine, win, subj, wh, sdir, evalue=1e-10)
                best = {}
                for l in open(wh):
                    q, s, e, bs, pi, ln, ql, sl = l.rstrip().split("\t")
                    if q not in best or float(bs) > best[q][1]:
                        r = bypid.get(s)
                        if r: best[q] = (s, float(bs), pi, r["gene"], r["chr"], r["start"])
                with open(os.path.join(sdir, "window_hits.tsv"), "w") as f:
                    f.write("nipponbare_gene\thit_protein\tbits\tidentity\thit_gene\tchr\tstart\n")
                    for q, (s, bs, pi, g, c, st) in best.items(): f.write(f"{q}\t{s}\t{bs}\t{pi}\t{g}\t{c}\t{st}\n")
            for tmp in ("_proteome.faa", "seed_hits.tsv", "window_hits_raw.tsv", "pfam.domtbl"):
                p = os.path.join(sdir, tmp)
                if os.path.exists(p): os.remove(p)
            for p in os.listdir(sdir):
                if p.startswith("db"): os.remove(os.path.join(sdir, p))
            open(done, "w").write("ok")
            if os.path.exists(os.path.join(sdir, "FAILED.txt")): os.remove(os.path.join(sdir, "FAILED.txt"))
            manifest.append((sp, pep_u, cds_u, len(genes), len(cand)))
            print(f"   genes={len(genes)} candidates={len(cand)}")
            if not a.keep_downloads and not a.local_dir:
                for p in (pep, cds):
                    if os.path.exists(p): os.remove(p)
        except Exception as e:
            print("   FAILED:", e); open(os.path.join(sdir, "FAILED.txt"), "w").write(str(e))
        finally:
            for tmpn in ("_proteome.faa", "seed_hits.tsv", "window_hits_raw.tsv", "pfam.domtbl"):
                pp = os.path.join(sdir, tmpn)
                if os.path.exists(pp): os.remove(pp)
    with open(os.path.join(a.out, "manifest.tsv"), "a") as f:
        for m in manifest: f.write("\t".join(map(str, m)) + "\n")
    zp = "phase2_upload.zip"
    with zipfile.ZipFile(zp, "w", zipfile.ZIP_DEFLATED) as z:
        for root, _, files in os.walk(a.out):
            if "_downloads" in root: continue
            for fn in files:
                if fn.endswith((".hmm", ".faa")) and fn != "candidates.faa": continue
                if fn.startswith("_"): continue
                z.write(os.path.join(root, fn), os.path.relpath(os.path.join(root, fn), a.out))
    failed = [sp for sp in plan if os.path.exists(os.path.join(a.out, sp, "FAILED.txt")) and not os.path.exists(os.path.join(a.out, sp, "DONE"))]
    done_n = sum(1 for sp in plan if os.path.exists(os.path.join(a.out, sp, "DONE")))
    print(f"\nFinished genomes: {done_n} of {len(plan)}")
    if failed: print("NOT finished (run the same command again to retry only these):", ", ".join(failed))
    print("\nDone. Upload this file to the chat:", os.path.abspath(zp), f"({os.path.getsize(zp)/1e6:.1f} MB)")

if __name__ == "__main__": main()
