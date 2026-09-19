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
  --outdir celltaminate_results
```

Metadata are optional. Supply `--metadata metadata.tsv` when group or sample-type annotations are useful.

## 5. Launch the Shiny application

```bash
bin/celltaminate app
```
