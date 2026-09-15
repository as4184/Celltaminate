# Score interpretation

Celltaminate reports a 0-100 score for each retained genus- or species-level taxon. Lower scores indicate stronger support for prioritization, while higher scores indicate stronger background-like evidence.

The deployed model integrates retained terms from abundance, unique k-mer support, taxonomic ambiguity, sterile-reference background, decontamination, curated clinical-pathogen and kitome/background membership, and selected interactions between these evidence layers.

The default operating thresholds are:

```text
score <= 1        prioritized / likely true
1 < score < 75    uncertain
score >= 75       likely background
```

The prioritization threshold can be changed with `--fp-true-cutoff`, and the upper background threshold can be changed with `--fp-falsepos-cutoff`.

Example:

```bash
bin/celltaminate run \
  --input-list input_list.txt \
  --outdir results \
  --fp-true-cutoff 1 \
  --fp-falsepos-cutoff 75
```

The score is a logistic mapping of the fitted contamination model and should not be interpreted as a calibrated probability.
