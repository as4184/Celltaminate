# Troubleshooting

## Metadata rows are ignored

The `sample` column must match the report basename exactly. If your report path is:

```text
/data/project/S6.kraken.report.txt
```

the metadata sample value should be:

```text
S6.kraken.report.txt
```
## Scores look too permissive or too strict

Adjust the score thresholds or scoring strictness:

```bash
--fp-true-cutoff 5     --fp-falsepos-cutoff 75     --fp-aggressiveness 1.15
```

Higher aggressiveness makes the score more willing to classify background-like taxa as likely false positives/background.

## Kraken report format error

Confirm that your input is a Kraken report, not a Kraken output file. Celltaminate expects summary report rows with percentage, clade reads, direct reads, rank, taxid, and taxon name.
