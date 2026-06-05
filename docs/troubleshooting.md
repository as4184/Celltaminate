# Troubleshooting

## `Rscript` is not found

Activate the conda environment first:

```bash
conda activate celltaminate
```

Or provide the path explicitly:

```bash
bin/celltaminate run --rscript /path/to/Rscript ...
```

## Metadata rows are ignored

The `sample` column must match the report basename exactly. If your report path is:

```text
/data/project/S6.kraken.report.txt
```

the metadata sample value should be:

```text
S6.kraken.report.txt
```

## The full reference table is missing

Use:

```bash
--ref-bg-tsv data/reference/refined_cell.lines.tsv
```

The bundled `refined_cell.lines.example.tsv` is only for testing. For real analysis, download the full table from Zenodo DOI `10.5281/zenodo.20560460`.

## Scores look too permissive or too strict

Adjust the score thresholds or scoring strictness:

```bash
--fp-true-cutoff 5     --fp-falsepos-cutoff 75     --fp-aggressiveness 1.15
```

Higher aggressiveness makes the score more willing to classify background-like taxa as likely false positive/background.

## Kraken report format error

Confirm that your input is a Kraken report, not a Kraken output file. Celltaminate expects summary report rows with percentage, clade reads, direct reads, rank, taxid, and taxon name.
