# Quick start

## Run the example

```bash
conda activate celltaminate
bash examples/lung_nanopore/run_example.sh
```

## Reference background table

Download the full reference background table from Zenodo before running Celltaminate:

```bash
bin/celltaminate download-reference --out data/reference/refined_cell.lines.tsv
```

If the helper cannot resolve the file URL, download `refined_cell.lines.tsv` manually from the Zenodo DOI landing page and place it at `data/reference/refined_cell.lines.tsv`.

DOI: `10.5281/zenodo.20560460`

## Start from FASTQ files

Celltaminate starts from Kraken2/KrakenUniq reports. To run the upstream FASTQ classification step and then Celltaminate:

```bash
python scripts/celltaminate_full_workflow.py \
  --samplesheet samples.tsv \
  --kraken2-db /path/to/kraken2_database \
  --outdir celltaminate_full_run \
  --threads 16 \
  --report-minimizer-data \
  --ref-bg-tsv data/reference/refined_cell.lines.tsv
```

The same workflow is available through the main wrapper:

```bash
bin/celltaminate from-fastq --samplesheet samples.tsv --kraken2-db /path/to/kraken2_database --outdir celltaminate_full_run
```

See [FASTQ to Celltaminate workflow](full_fastq_workflow.md) for the samplesheet format and standalone Kraken2 commands.

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

## Create metadata templates

```bash
bin/celltaminate init --outdir template_files
```
