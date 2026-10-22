# Stage C - model-based phylogeny and selection test (run on your Mac)

What this adds to the paper: (1) a proper maximum-likelihood tree (ModelFinder, 1000 ultrafast bootstraps, SH-aLRT) with sensitivity runs,
replacing the neighbour-joining prototype; (2) a branch test of whether the OsMC1 lineage evolves under relaxed selection, using HyPhy RELAX
on a Type I codon alignment (OsMC1 branches = Test, other Type I branches = Reference).

## One-time install (same style as before)
    conda install -c bioconda -c conda-forge mafft iqtree clipkit hyphy
If hyphy will not install, install the others and run anyway (RELAX is skipped and I will use another approach). Do NOT use `conda create` with
`CONDA_SUBDIR=osx-64` - that failed for you before; installing into the base environment worked for BLAST.

## Run
    cd "/Volumes/movies and others/Oryza_MC_v2/phase2/stageC"
    caffeinate -i bash run_stageC.sh
Takes roughly 10-60 minutes. It prints a progress line for each of three steps and finally the path of `stageC_upload.zip`. Upload that file.

## Files
* `family_proteins.faa` - 160 proteins: reciprocal-best-hit orthologs of OsMC1-9 from 27 genomes (full-length-ish models only), the MSU Nipponbare models, and AtMC1-9.
* `typeI_proteins.faa`, `typeI_cds.fna` - 86 Type I orthologs (OsMC1, 2, 3, 8, 9) for the selection test.
* `make_codon_alignment.py`, `label_tree.py` - tested here on real and synthetic data; the external tools (MAFFT, IQ-TREE, HyPhy) were NOT run in my sandbox.

## Caveats I will report honestly
* If the OsMC1 sequences are not a monophyletic clade in the codon tree, RELAX is not run; that is itself informative.
* RELAX results depend on sampling (OsMC8 has only 6 usable sequences); I will not over-interpret a borderline result.
* Gene models that were fragmentary were left out, which reduces the OsMC2, OsMC3, OsMC8 and OsMC9 sample.
