# Step 4 — developmental / stress expression (public RNA-seq)

No expression data of any kind (raw, processed, or as an accession number) was present in
the uploaded project, so nothing here can be verified yet. This is a ready-to-run recipe,
not a re-analysis of anything you sent.

## Candidate public datasets (verify current accession status before use — public
## repositories occasionally withdraw or supersede series)
- Rice expression atlas (developmental): search ArrayExpress/EBI Expression Atlas or NCBI
  GEO for "Oryza sativa developmental RNA-seq atlas" (e.g. series built on the RGASP /
  Sato et al. rice tissue panel).
- Abiotic stress (drought/salt/cold): search GEO for rice drought/salinity RNA-seq time
  courses; several multi-tissue, multi-stress series exist — pick ones with clear tissue,
  timepoint, and replicate metadata, and cite the accession explicitly in Methods.

## Processing (fully reproducible, put every parameter in Methods)
```bash
# QC
fastqc raw/*.fastq.gz
# pseudoalignment against the SAME transcriptome version used for the GFF3 in Step 1
salmon index -t Osativa_323_v7.0.transcript.fa -i salmon_index
salmon quant -i salmon_index -l A -1 R1.fastq.gz -2 R2.fastq.gz -p 8 -o quant/sample1
# gene-level counts -> DESeq2
tximport + DESeq2 in R:
#   dds <- DESeqDataSetFromTximport(txi, colData, ~condition)
#   dds <- DESeq(dds); res <- results(dds, contrast=c("condition","stress","control"))
#   report log2FC + BH-adjusted P, with the exact fold-change/P-value thresholds stated
```

Report explicitly, per sample group: tissue, developmental stage or stress type, duration,
number of biological replicates, normalization method, and the differential-expression
threshold used (e.g. |log2FC|>1, padj<0.05). Build the heatmap from normalized (e.g. VST or
TPM) values, not raw counts. Language for the Results/Discussion: "OsMCx showed
stress-associated expression changes (log2FC=…, padj=…) under drought at Xh" — never "OsMCx
regulates drought tolerance" without functional/genetic evidence.
