# Score interpretation

Celltaminate reports a score from 0 to 100 for each genus or species-level taxon. Lower scores indicate stronger support for prioritization, whereas higher scores indicate stronger support for contaminant or technical artifact.

Celltaminate integrates features from taxonomy reports (e.g., Kraken2, KrakenUniq) such as abundance, unique k-mer support, taxonomic ambiguity, sterile-reference background, decontamination, curated clinical-pathogen and kitome membership, and selected interactions among these features.

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
