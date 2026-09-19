# Output files

For each provided taxonomy report, Celltaminate writes intuitive tables, plots, cleaned reports, and a run summary. Cleaned Kraken-style reports retain prioritized taxa by default, but a user can also choose to retain selected organisms from the shiny app. The original taxonomic report is not modified by Celltaminate.

The main taxon-level outputs are in `output_dir/tables/per_sample/`. Important columns include taxon name, clade reads, RPM, RPMM, decontaminated RPMM, `Celltaminate_Score`, `Call`, and `Call_reason`.

The binary call for a taxon is:

```text
score <= 1    Prioritized
score > 1     Not prioritized
```
