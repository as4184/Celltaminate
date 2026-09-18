# Score interpretation

Celltaminate reports a score from 0 to 100 for each retained genus- or species-level taxon. Lower scores indicate stronger support for prioritization, whereas higher scores indicate stronger background-like evidence.

The deployed model integrates retained terms from abundance, unique k-mer support, taxonomic ambiguity, sterile-reference background, decontamination, curated clinical-pathogen and kitome/background membership, and selected interactions between these evidence layers.

The default operating rule is:

```text
score <= 1    prioritized
score > 1     not prioritized
```

The prioritization threshold can be changed with `--fp-true-cutoff`. Changing the threshold changes the binary call but does not change the underlying Celltaminate score.

Example:

```bash
bin/celltaminate run \
  --input-list input_list.txt \
  --outdir results \
  --fp-true-cutoff 1
```

The Celltaminate score is a logistic mapping of the fitted regularized model and should not be interpreted as a calibrated probability.
