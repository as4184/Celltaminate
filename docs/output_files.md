# Output files

A Celltaminate run writes tables, plots, cleaned reports, and a run summary.

The main taxon-level outputs are in `output_dir/tables/per_sample/`. Important columns include taxon name, clade reads, RPM, RPMM, decontaminated RPMM, `Celltaminate_Score`, `Call`, and `Call_reason`.

The binary call is:

```text
score <= 1    Prioritized
score > 1     Not prioritized
```

Cleaned Kraken-style reports retain prioritized taxa by default. The original taxonomic report is not modified.
