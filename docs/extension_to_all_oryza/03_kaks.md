# Step 3 — Ka/Ks (selective pressure on duplicated / orthologous MC pairs)

Only run this on pairs already established as real duplicates/orthologs in Step 2 (MCScanX
output) — not on every possible gene pair, which invites nonsense ratios from saturated or
non-homologous comparisons.

```bash
# For each pair (protein_a, protein_b) with CDS_a.fa / CDS_b.fa:
mafft --auto <(cat protein_a.fa protein_b.fa) > pair.aln.fa   # protein alignment
perl pal2nal.pl pair.aln.fa <(cat CDS_a.fa CDS_b.fa) -output paml > pair.codon.aln

# KaKs_Calculator 2.0, model averaging (NG, YN, MLWL, etc.) — report YN00 by default,
# it's the most commonly cited in plant gene-family papers
KaKs_Calculator -i pair.codon.aln -o pair.kaks -m YN

# Batch over every duplicate/ortholog pair, then tabulate:
#   gene_pair  Ka  Ks  Ka/Ks  Ks_synonymous_saturation_flag(Ks>2 -> exclude, too saturated)
```

Interpretation to put in the paper (do not oversimplify a single ratio, per the earlier
review): Ka/Ks « 1 across most pairs = purifying selection maintaining the catalytic
machinery; any individual pair with Ka/Ks > 1 needs its own sentence (which codons, is it
in the pro-domain vs the catalytic core) rather than a blanket "positive selection" claim.
Also report the Ks distribution (histogram) across all pairs as a rough duplication-age
proxy — cite the synonymous-substitution-rate caveats for grasses rather than converting Ks
to an absolute date unless you have a rice-specific clock calibration to cite.
