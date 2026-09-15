# Input formats

## Kraken reports

Celltaminate accepts:

1. Kraken2 6-column reports;
2. Kraken2 8-column reports with minimizer data;
3. KrakenUniq 8-column reports.

The original indentation of taxon names must be preserved because it is used to reconstruct the taxonomic lineage.

Eight-column reports are preferred because the fitted score uses unique k-mer/minimizer support.

## Input list

`input_list.txt` contains one report path per line:

```text
sample1.kraken.report.txt
sample2.kraken.report.txt
sample3.kraken.report.txt
```

Relative paths are interpreted relative to the input-list location.

## Metadata

Optional metadata must be tab-separated:

```text
sample	group	sample_type	is_control
sample1.kraken.report.txt	Case	Nanopore	FALSE
sample2.kraken.report.txt	Control	Nanopore	TRUE
```

Columns:

- `sample`: report basename;
- `group`: biological or study group;
- `sample_type`: sample or sequencing description;
- `is_control`: `TRUE` or `FALSE`.
