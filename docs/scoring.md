# Score interpretation

Celltaminate reports a 0-100 score for each retained genus or species-level taxon.

```text
lower score  = stronger evidence for likely true microbial signal
higher score = stronger evidence for likely false-positive/background signal
```

These thresholds are configurable:

```bash
bin/celltaminate run       --input-list input_list.txt       --outdir results       --fp-true-cutoff 5       --fp-falsepos-cutoff 75       --fp-aggressiveness 1.15
```
