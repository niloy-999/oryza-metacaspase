# Step 1 — genus-wide metacaspase identification (all 11 Oryza species)

Prerequisite: genome FASTA + GFF3/GTF annotation + protein FASTA for each species, from
Phytozome, EnsemblPlants, or the Rice Genome Annotation Project (RGAP). Record for each:
species, genome assembly version, annotation version, database accession, download date —
put this in `Oryza_MC/genomes/MANIFEST.tsv` before running anything (this is the "single
source of truth" table the earlier review recommended).

```bash
mkdir -p Oryza_MC/{genomes,proteins,hmmer,blast,domains,identification}
cd Oryza_MC

# 1. HMM search — primary evidence, not BLAST identity
#    PF00656 = Peptidase_C14, the metacaspase/paracaspase catalytic domain
wget https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/PF00656?annotation=hmm -O PF00656.hmm
for sp in Osativa Oglaberrima Orufipogon Onivara Obarthii Oglumaepatula Omeridionalis Obrachyantha Opunctata Olongistaminata OspAA; do
  hmmsearch --domtblout hmmer/${sp}.domtblout --cut_ga PF00656.hmm proteins/${sp}.protein.fa
done

# 2. Reciprocal BLASTP confirmation against the 9 validated AtMC + 8 validated OsMC proteins
#    (seeds already extracted in ../results/sequences/seeds_Arabidopsis_AtMC1-9.fa —
#     add the 8 OsMC1-8 sequences from ../results/tables/Table1_domain_architecture_types.tsv)
cat seeds_AtMC1-9.fa seeds_OsMC1-8.fa > seeds_combined.fa
makeblastdb -in seeds_combined.fa -dbtype prot
for sp in ...; do
  blastp -query proteins/${sp}.protein.fa -db seeds_combined.fa -evalue 1e-10 \
         -outfmt 6 -max_target_seqs 5 -num_threads 8 > blast/${sp}_vs_seeds.tsv
done

# 3. Domain validation — every HMM hit must ALSO have Pfam/CDD/InterPro Peptidase_C14
#    (use InterProScan locally, or NCBI Batch CD-Search as you already did for rice)
interproscan.sh -i proteins/${sp}.protein.fa -appl Pfam -f tsv -o domains/${sp}_pfam.tsv

# 4. Candidate = HMM hit (--cut_ga) AND reciprocal BLAST hit (E<1e-10) AND Pfam Peptidase_C14 present
#    Merge these three files per species with a simple join on protein ID; keep only the
#    intersection. This is the step your original pipeline skipped, and it is exactly what
#    would have caught LOC_Os09g36640 before it ever reached Table 2.

# 5. Redundancy handling — genomic LOCUS is the gene identifier, not transcript/isoform.
#    Use the "longest=1" / canonical-transcript flag in the GFF3 (same field this
#    re-analysis used for O. sativa) to pick one representative protein per locus.
```

Output: `identification/Table_S1_full_MC_inventory.tsv` with columns
`species, locus_id, transcript_id, chromosome, start, end, strand, protein_length,
HMM_bitscore, HMM_evalue, reciprocal_best_hit, Pfam_domain, Pfam_evalue`.

**Sanity check before proceeding to any downstream analysis:** the total gene count you
report (previously "93") should emerge from this table, not be assumed going in. If your
final count differs from 93, that is fine and expected — say so plainly rather than tuning
thresholds to hit a pre-committed number. (The 2020 Oryza MC paper's count of 93 across ten
species is a useful cross-check, not a target to reproduce exactly.)
