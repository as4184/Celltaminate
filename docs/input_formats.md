# Input formats

## Kraken report files

Celltaminate starts from taxonomic summary reports. It accepts:

1. Kraken2 6-column report files;
2. Kraken2 8-column report files;
3. KrakenUniq 8-column report files.

The workflow uses the indentation structure of the report to reconstruct lineage relationships. For this reason, the original taxon-name indentation should be preserved.

## Input list

`input_list.txt` should contain one report path per line:

```text
sample1.kraken.report.txt
sample2.kraken.report.txt
sample3.kraken.report.txt
```

Relative paths are interpreted relative to the location of the input list by the wrapper.

## Metadata table

`metadata.tsv` must be tab-separated:

```text
sample	group	sample_type	is_control
sample1.kraken.report.txt	Case	Nanopore	FALSE
sample2.kraken.report.txt	Control	Nanopore	TRUE
```

Required columns:

- `sample`: report basename;
- `group`: biological or study group;
- `sample_type`: sequencing or sample description;
- `is_control`: `TRUE` or `FALSE`.
