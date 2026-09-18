# Input formats

## Kraken reports

Celltaminate accepts Kraken2 6-column reports, Kraken2 8-column reports with minimizer data, and KrakenUniq 8-column reports. The original indentation of taxon names must be preserved because it is used to reconstruct the taxonomic lineage. Eight-column reports are preferred because the fitted score uses unique k-mer/minimizer support.

## Input list

`input_list.txt` contains one report path per line:

```text
sample1.kraken.report.txt
sample2.kraken.report.txt
```

Relative paths are interpreted relative to the input-list location.

## Metadata

Metadata are optional and tab-separated:

```text
sample	group	sample_type
sample1.kraken.report.txt	Case	Nanopore
sample2.kraken.report.txt	Case	Nanopore
```

Columns:

- `sample`: report basename;
- `group`: biological or study group;
- `sample_type`: sample or sequencing description.
