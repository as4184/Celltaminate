# Score interpretation

Celltaminate reports a 0-100 score for each retained genus or species-level taxon.

```text
lower score  = stronger evidence for likely true microbial signal
higher score = stronger evidence for likely false-positive/background signal
```

Default interpretation:

```text
score <= 5   likely true
score >= 75  likely false positive/background
otherwise    uncertain
```

These thresholds are configurable:

```bash
bin/celltaminate run       --input-list input_list.txt       --outdir results       --fp-true-cutoff 5       --fp-falsepos-cutoff 75       --fp-aggressiveness 1.15
```

## Main evidence layers

Celltaminate combines:

- low read-count and low RPMM penalties;
- unique k-mer support when available;
- reference background abundance and prevalence;
- fold enrichment above reference background;
- optional matched-control enrichment and prevalence;
- cohort recurrence and ubiquity;
- microbial biomass and non-microbial fraction relationships;
- taxonomic ambiguity;
- kitome/background blacklist membership;
- clinically important pathogen membership.

The score is intended for prioritization and triage. It should not be interpreted as a calibrated clinical probability without additional validation.
