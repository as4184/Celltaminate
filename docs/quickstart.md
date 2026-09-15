# Quick start

## 1. Activate the environment

```bash
conda activate celltaminate
```

## 2. Download the reference background table

```bash
bin/celltaminate download-reference
```

## 3. Run the bundled example

```bash
bash examples/lung_nanopore/run_example.sh
```

## 4. Run your own Kraken reports

```bash
bin/celltaminate run \
  --reports sample1.kraken.report.txt sample2.kraken.report.txt \
  --metadata metadata.tsv \
  --outdir celltaminate_results
```

Metadata are optional. Without metadata, reports are assigned to a default group.

## 5. Launch the Shiny application

```bash
bin/celltaminate app
```
