#!/bin/bash
# Step 0: compile the Smith-Waterman search tool and run the seeded searches used by 01_identify_annotate.py.
# Requires the Phytozome files described in the README under data/raw (or set MC_PROJ).
set -e
cd "$(dirname "$0")"
TMP="${MC_SWDIR:-/tmp}"; mkdir -p "$TMP"
gcc -O3 -o sw sw.c
python3 - <<PY
import os, sys
sys.path.insert(0, ".")
from common import *
os.makedirs(os.environ.get("MC_SWDIR", "/tmp"), exist_ok=True)
t = os.environ.get("MC_SWDIR", "/tmp")
write_fasta({k: v.rstrip("*") for k, v in read_fasta(RICE_PROT).items()}, f"{t}/rice_db.fa")
write_fasta({k: v.rstrip("*") for k, v in read_fasta(AT_PROT, keyfn=lambda h: h.split()[0].split("|")[0]).items()}, f"{t}/at_db.fa")
PY
SEEDS=../results/sequences/seeds_Arabidopsis_AtMC1-9.fa
./sw "$SEEDS" "$TMP/rice_db.fa" 50 > "$TMP/sw_rice.tsv"
./sw "$SEEDS" "$TMP/at_db.fa" 50 > "$TMP/sw_at.tsv"
echo "search done; now run: python3 01_identify_annotate.py"
