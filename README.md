# Rice metacaspase family: catalytic-residue audit and cross-genome conservation

Reproducible pipeline for a residue-level audit of the metacaspase (C14 peptidase) family in rice (*Oryza sativa*), and a test of whether the
catalytic-residue substitutions found in **OsMC1** are conserved across 27 *Oryza* and *Leersia* genome annotations.

> Status: analysis code and derived results for a manuscript in preparation. Numbers may change after peer review.

## Main findings (see `results/`)
* The rice family has **nine** metacaspases: OsMC1-8 (MSU v7.0) plus a chromosome-10 gene annotated only in RAP-DB (Os10g0565100, provisionally OsMC9),
  the closest rice relative of *Arabidopsis* AtMC1 (about 77% identity).
* OsMC1 has lost **both** catalytic residues (His -> Ile/Leu, Cys -> Glu). Its reciprocal-best-hit ortholog carries the same substitutions in **all 27 genomes**
  examined (17 *O. sativa* assemblies, 9 other *Oryza* species and *Leersia perrieri*), with conserved gene order in the chromosome-3 cluster;
  none of 133 assessable orthologs of the other eight groups has a substitution.
* Ka/Ks across divergent OsMC1 orthologs is 0.19-0.46 (purifying selection, weaker than on intact paralogs).

## Repository layout
| Path | Content |
|---|---|
| `scripts/` | Steps 00-15 for the *O. sativa* / *Arabidopsis* analysis (pure Python 3 + a small C aligner) |
| `phase2/` | Cross-genome pipeline: downloader/extractor (`01_`), ortholog and catalytic-residue analysis (`20_`), structure geometry (`40_`), model-based tree and selection test package (`stageC/`) |
| `results/` | Derived tables, figures, alignments, trees, sequences |
| `docs/extension_to_all_oryza/` | Recipes for extending the analysis (HMMER, MCScanX, Ka/Ks, expression) |

## Requirements
Python >= 3.10 with `numpy` and `matplotlib`; a C compiler. For the cross-genome and tree steps: BLAST+, HMMER (optional), MAFFT, IQ-TREE 2, ClipKIT, HyPhy (optional).
`conda env create -f environment.yml`

## Input data (not included; third-party)
Download from Phytozome v13 and place under `data/raw/` (or set `MC_PROJ`):
`O.Sativa/Osativa_323_v7.0.{fa,gene.gff3,protein.fa}` and `Arabidopsis thaliana/Athaliana_167_{protein.fa,gene.gff3}`.
The cross-genome step downloads Ensembl Plants release 63 FASTA files itself.

## Run order
```bash
bash scripts/00_search.sh                 # compile aligner, run seeded searches
for s in 01_identify_annotate 02_gene_models_physchem 03_domains_types 04_msa_tree 05_chrom_duplication \
         06_promoters 07_gene_structure_motifs 08_domain_figure 09_kaks 10_dyad_figure 12_integrity_and_degradation 13_tree_sensitivity; do
  python3 scripts/$s.py; done
# cross-genome analysis (needs internet; ~1-2 GB temporary downloads)
python3 phase2/01_fetch_and_extract.py
python3 phase2/20_stageB_orthology_dyad.py phase2/phase2_out phase2/results_stageB
python3 scripts/14_crossgenome_table_figure.py && python3 scripts/15_ortholog_kaks.py
# model-based tree and selection test: see phase2/stageC/README_STAGEC.md
```
All steps are deterministic (fixed random seeds). An independent re-run on a second computer reproduced 14 of 16 result files byte-for-byte.

## Limitations
Computational study of annotations only; no new expression, structural or functional data. Gene models are fragmentary in some genomes.
The initial search is a custom Smith-Waterman implementation (positively controlled with *Arabidopsis*; cross-checked with BLAST+/HMMER in the cross-genome step).

## Citation / license
See `CITATION.cff`. Code: MIT. Derived tables and figures: please cite the associated paper when available.
