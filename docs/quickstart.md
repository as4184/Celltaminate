# Quick start

## Run the example

```bash
conda activate celltaminate
bash examples/lung_nanopore/run_example.sh
```

## Reference background table

For real analyses, download the full reference background table from Zenodo before running Celltaminate:

```bash
bin/celltaminate download-reference --out data/reference/refined_cell.lines.tsv
```

Alternatively:

```bash
mkdir -p data/reference
curl -L -o data/reference/refined_cell.lines.tsv \
  "https://zenodo.org/records/20560460/files/refined_cell.lines.tsv?download=1"
```

DOI: `10.5281/zenodo.20560460`

## Run your own Kraken reports

```bash
bin/celltaminate run \
  --reports sample1.kraken.report.txt sample2.kraken.report.txt \
  --metadata metadata.tsv \
  --outdir celltaminate_results \
  --tax-level S \
  --ref-bg-tsv data/reference/refined_cell.lines.tsv
```

## Run without metadata

Metadata is optional. Without metadata, every report is assigned to `Group 1`, sample type `Clinical / sterile`, and `is_control = FALSE`.

```bash
bin/celltaminate run \
  --reports sample1.kraken.report.txt \
  --outdir celltaminate_results \
  --tax-level S
```

This is acceptable for a first technical run, but reference background and metadata are recommended for real interpretation.

## Create metadata templates

```bash
bin/celltaminate init --outdir template_files
```
