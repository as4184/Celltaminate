# Output files

A Celltaminate run writes tables, plots, optional cleaned reports, and a run summary.

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

## Main result tables

The main taxon-level outputs are in:

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
- `Call`: three-tier software call;
- `Call_reason`: concise explanation of the call.

## Cleaned reports

Celltaminate can write cleaned Kraken-style reports after removing user-selected call categories. The original taxonomic report is not modified.
