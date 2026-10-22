# Step 2 — synteny and duplication (within- and between-species)

## Within-species (tandem/proximal/segmental), per species
```bash
# MCScanX input: <species>.gff (BLAST-derived gene order, MCScanX's own format, NOT GFF3)
#   col1=chr  col2=gene_id  col3=start  col4=end   (tab-separated, no header)
awk -F'\t' '$3=="gene"{split($9,a,";"); gsub("ID=","",a[1]); print $1"\t"a[1]"\t"$4"\t"$5}' \
    annotation/Osativa.gff3 > mcscanx/Osativa.gff

# all-vs-all BLASTP within the species (MCScanX needs its own blast format, -outfmt 6)
blastp -query proteins/Osativa.fa -db proteins/Osativa.fa -evalue 1e-10 \
       -outfmt 6 -num_threads 8 > mcscanx/Osativa.blast

MCScanX mcscanx/Osativa   # reads Osativa.gff + Osativa.blast, writes Osativa.collinearity
# classify each MC-gene pair's duplicate_gene_classifier output as: singleton / dispersed /
# proximal / tandem / WGD-or-segmental  -- this replaces the ad hoc tandem/segmental split
# in the original manuscript with MCScanX's own, citable classification.
duplicate_gene_classifier mcscanx/Osativa
```

## Between-species (cross-Oryza synteny / orthologous blocks)
```bash
# concatenate two species' gff + blast (query vs subject cross-species blastp),
# label chromosomes with a species prefix so MCScanX can tell them apart, e.g. "Os_Chr1" vs
# "Og_Chr1"
cat mcscanx/Osativa.gff mcscanx/Oglaberrima.gff > mcscanx/Os_Og.gff
blastp -query proteins/Osativa.fa -db proteins/Oglaberrima.fa -evalue 1e-10 -outfmt 6 \
       -num_threads 8 > mcscanx/Os_Og.blast
MCScanX mcscanx/Os_Og
# repeat for every species pair of interest (or run once against a concatenated
# all-species proteome for a single genus-wide synteny map — heavier but gives Figure 5
# directly: "cross-Oryza synteny of the MC locus")
```

Then intersect the `.collinearity` blocks with the MC gene inventory from Step 1 to answer:
*which MC loci have conserved syntenic positions across Oryza, and which are lineage-specific
insertions/losses?* That comparison — not a bare "duplication contributed to expansion"
sentence — is what turns this into Figure 5 of the paper.
