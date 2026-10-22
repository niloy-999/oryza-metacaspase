# Phase 2 - make the OsMC1 finding robust (cross-accession, cross-species, tree, structure, selection)

Why this exists: the one weakness a reviewer will press is that the OsMC1 substitutions (His->Ile, Cys->Glu) are confirmed only in Nipponbare.
The sandbox I work in has no internet and cannot run IQ-TREE / MAFFT / codeml, so the work is split: **you run downloads and heavy tools, I analyse and write.**

## Stage A (you, ~1-2 h, mostly waiting for downloads)
> Run it as `caffeinate -i python3 01_fetch_and_extract.py` so the Mac does not sleep and drop the connection. It resumes safely if interrupted; re-run the same command to retry any genome that failed. (v2 of the script works in folders whose names contain spaces.)
1. Open Terminal in this `phase2` folder.
2. `conda install -c bioconda -c conda-forge blast hmmer`   (diamond works instead of blast; hmmer is optional but recommended)
3. Quick test on two genomes:  `python3 01_fetch_and_extract.py --max-species 2`
4. If it prints genomes and candidates, run everything:  `python3 01_fetch_and_extract.py`
   * It reads the GrameneOryza / Ensembl Plants FTP (28 Oryza assemblies incl. MH63, IR64, ZS97, Aus, Azucena, ARC, LiuXu and wild species, plus *Leersia* as outgroup, if present).
   * Needs about 1-2 GB of temporary downloads (deleted automatically). It resumes if interrupted.
5. Upload the file it prints at the end: `phase2_upload.zip` (small, tens of MB).
If it cannot reach the FTP servers: download `pep.all.fa.gz` and `cds.all.fa.gz` for each genome by hand from the species pages on plants.ensembl.org or ensembl.gramene.org ("Download genes, cDNAs, ncRNA, proteins - FASTA"), put them in `DIR/<species>/pep.fa.gz` and `DIR/<species>/cds.fa.gz`, and run with `--local-dir DIR`. Tell me the error text if anything fails.

## Stage B (me) - done on your upload
Reciprocal-best-hit orthology of every candidate to OsMC1-8 with the score margin to the runner-up (so OsMC1 cannot be confused with OsMC2/3/8), conserved-genomic-context test (Nipponbare chromosome-3 window genes within +-15 genes), and a direct read of the catalytic residues in every accession and species. Script: `20_stageB_orthology_dyad.py` (already tested on mock data).

## Stage C (you, after Stage B) - model-based tree and selection test
I will send a ready-to-run script. Install once: `conda install -c bioconda mafft iqtree clipkit paml hyphy`.
It will run MAFFT -> ClipKIT -> IQ-TREE 2 (ModelFinder, 1000 ultrafast bootstraps, SH-aLRT) with an alternative model as a sensitivity check, and, only if there are enough OsMC1-lineage orthologs and Ks is not saturated, a branch-model test (codeml) and/or HyPhy RELAX. Upload the .treefile/.iqtree/.log and the codeml/RELAX outputs.

## Structure (you, 20 minutes)
1. Go to alphafold.ebi.ac.uk, search the UniProt entries for rice OsMC1 (Os03g0388900 / LOC_Os03g27120), OsMC2 (Os03g0389400) and OsMC3 (Os03g0389000), download each model as PDB. If a gene is not found, run its sequence from `seeds_all.fa` in AlphaFold Server and download the mmCIF.
2. Download an experimental plant metacaspase structure from rcsb.org (search "metacaspase 9 Arabidopsis") if one exists, as a reference.
3. Upload the files. I measure the His/Cys equivalents (OsMC1 173/229, OsMC2 173/229, OsMC3 195/251, AtMC9 95/147), per-residue confidence (pLDDT), and local pocket RMSD with `40_structure_geometry.py` (tested on a toy model). A prediction cannot prove inactivity; it can show whether the pocket is disrupted.

## Expression / transcript support (you, 10 minutes)
Open the RAP-DB locus page for Os03g0388900 and note (a) whether full-length cDNA / transcript models support the gene, and (b) the RiceXPro or Expression Atlas expression values. Screenshot or paste them. Published qRT-PCR for OsMC1 (Huang et al. 2015) already shows it is expressed; public RNA-seq would add independent support. Do not claim more than the data shows.

## How I will interpret the Stage B result
| Outcome | Meaning | Paper |
|---|---|---|
| Both substitutions in the OsMC1 ortholog of most/all sativa accessions and in wild species | Ancient, conserved loss of catalysis | Strongest story: conserved pseudo-metacaspase in *Oryza* |
| Present in all japonica-type accessions, absent in indica / wild | Lineage- or subspecies-specific loss | Still publishable; reframe as a recent, lineage-restricted event |
| Only in Nipponbare | Cultivar-specific variant | Weaker; report honestly as a variant, drop the "pseudoenzyme" headline |
| Ortholog not annotated in a genome | Inconclusive for that genome | I will search the genome directly (needs the genome FASTA) before drawing conclusions |
Selection tests are run only if taxon sampling and synonymous-site saturation allow; otherwise they are reported as not performed.


---
## STATUS UPDATE (Stage A and B completed)
Stage A ran on 27 genomes (Ensembl Plants release 63) and Stage B analysed them. Result: the OsMC1 ortholog lacks both catalytic residues in all 27 genomes (His: Ile in 25, Leu in 2; Cys: Glu in 27), with conserved gene order in 24/24 testable genomes and Ka/Ks 0.19-0.46 in the four most divergent pairs. Outputs: `results_stageB/`, `manuscript/supplementary/TableS14*`, `TableS15*`, manuscript Figure 2.

## OPTIONAL NEXT: extend sampling to other grasses (dates the origin of the loss)
```
caffeinate -i python3 01_fetch_and_extract.py --out phase2_grasses --want "brachypodium_distachyon|sorghum_bicolor|zea_mays|hordeum_vulgare|setaria_italica|oropetium_thomaeum"
```
(oryza_sativa is added automatically as the synteny reference.) Upload the new `phase2_upload.zip`. If the OsMC1 orthologs in these grasses keep both residues, the loss arose in the ancestor of Oryza + Leersia; if they lost them too, it is older.

## STILL TO DO (Stage C): model-based tree and branch test
Needs MAFFT + IQ-TREE 2 (+ HyPhy RELAX or codeml) on your computer; I will prepare the ortholog sets and scripts once the grass sampling is in.


## UPDATE: ninth rice metacaspase
The Stage B audit of extra candidates revealed Os10g0565100 (chromosome 10), now included as a ninth seed (`seeds_all.fa`; the old 17-seed file is `seeds_all_17.fa`). `results_stageB/` contains the 9-seed assignments. Orthologs of this gene: 24 of 27 genomes, all with intact His and Cys.
