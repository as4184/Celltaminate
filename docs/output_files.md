# Output files

A Celltaminate run writes tables, plots, cleaned reports, and a run summary.

```text
output_dir/
├── tables/
│   ├── all_taxa_table.tsv
│   ├── per_sample/
│   └── summaries/
├── plots/
├── cleaned_reports/
├── cleaned_reanalysis/
├── run_config/
└── run_summary.tsv
```

## Main result table

The most useful tables are usually in:

```text
output_dir/tables/per_sample/
```

Important columns include:

- `Name`: taxon name;
- `Reads`: clade reads;
- `RPM`: reads per million total reads;
- `RPMM`: reads per million microbial reads;
- `Decontaminated_RPMM`: background-subtracted microbial abundance;
- `Celltaminate_Score`: 0-100 contamination-aware score;
- `Call`: likely true, uncertain, or likely false positive/background;
- `Call_reason`: short explanation of the call.

## Cleaned reports

By default, Celltaminate writes cleaned Kraken-style reports after removing taxa called as likely false positive/background or uncertain. The cleaned reports are useful for downstream exploratory analysis, but they should not replace the original taxonomic reports.
