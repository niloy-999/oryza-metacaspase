#!/bin/bash
# STAGE C - run on YOUR computer, inside this folder.   Usage:  bash run_stageC.sh
# Needs: mafft, iqtree2 (or iqtree), clipkit (or trimal), python3; hyphy is optional (for the selection test).
set -u
THREADS=${THREADS:-4}
mkdir -p out
IQ=$(command -v iqtree2 || command -v iqtree || true)
[ -z "$IQ" ] && { echo "iqtree2 not found. Install:  conda install -c bioconda -c conda-forge mafft iqtree clipkit hyphy"; exit 1; }
command -v mafft >/dev/null || { echo "mafft not found."; exit 1; }
trim() {   # $1 in  $2 out
  if command -v clipkit >/dev/null; then clipkit "$1" -m kpic-gappy -o "$2" >/dev/null
  elif command -v trimal >/dev/null; then trimal -in "$1" -out "$2" -automated1
  else echo "WARNING: no trimming tool found; using the untrimmed alignment"; cp "$1" "$2"; fi
}
echo "== 1. family-wide protein tree (ModelFinder + 1000 ultrafast bootstraps + SH-aLRT)"
mafft --auto --thread "$THREADS" family_proteins.faa > out/family.aln.faa
trim out/family.aln.faa out/family.trim.faa
"$IQ" -s out/family.trim.faa -m MFP -B 1000 --alrt 1000 -T "$THREADS" --prefix out/family -redo
echo "== 1b. sensitivity: fixed LG+G4 model and an alignment trimmed less strictly"
"$IQ" -s out/family.trim.faa -m LG+G4 -B 1000 --alrt 1000 -T "$THREADS" --prefix out/family_LGG4 -redo
"$IQ" -s out/family.aln.faa -m MFP -B 1000 -T "$THREADS" --prefix out/family_untrimmed -redo

echo "== 2. Type I codon alignment and selection test"
mafft --auto --thread "$THREADS" typeI_proteins.faa > out/typeI.aln.faa
python3 make_codon_alignment.py out/typeI.aln.faa typeI_cds.fna out/typeI.codon.fna
"$IQ" -s out/typeI.codon.fna -m MFP -B 1000 -T "$THREADS" --prefix out/typeI_codon -redo
if python3 label_tree.py out/typeI_codon.treefile OsMC1__ out/typeI.labeled.nwk; then
  if command -v hyphy >/dev/null; then
    hyphy relax CPU="$THREADS" --alignment out/typeI.codon.fna --tree out/typeI.labeled.nwk --test Test --reference Reference --output out/relax.json > out/relax.log 2>&1 \
      || echo "RELAX failed - see out/relax.log (send me the last 20 lines)"
  else echo "hyphy not installed - skipping RELAX (install with: conda install -c bioconda -c conda-forge hyphy)"; fi
else echo "OsMC1 sequences are not a clade in the codon tree - RELAX not run (this is itself a result; send the files)"; fi

echo "== 3. packaging"
rm -f stageC_upload.zip
zip -q -r stageC_upload.zip out -x "out/*.ckp.gz" "out/*.uniqueseq.phy" "out/*.splits.nex" "out/*.mldist" "out/*.bionj" "out/*.model.gz"
echo "Done. Upload:  $(pwd)/stageC_upload.zip"
