#!/usr/bin/env bash
set -euo pipefail

# Run this from the repository root:
# bash examples/lung_nanopore/run_example.sh

bin/celltaminate run       --input-list examples/lung_nanopore/input_list.txt       --metadata examples/lung_nanopore/metadata.tsv       --outdir examples/lung_nanopore/results       --tax-level S       --ref-bg-tsv data/reference/refined_cell.lines.tsv
