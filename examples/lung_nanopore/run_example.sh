#!/usr/bin/env bash
set -euo pipefail

# Run this from the repository root:
# bash examples/lung_nanopore/run_example.sh

REF="data/reference/refined_cell.lines.tsv"
REF_ARGS=()

if [[ -f "${REF}" ]]; then
  REF_ARGS=(--ref-bg-tsv "${REF}")
else
  echo "WARNING: ${REF} was not found. Running example without the reference background layer." >&2
  echo "For full scoring, download refined_cell.lines.tsv from Zenodo DOI 10.5281/zenodo.20560460 and place it at ${REF}." >&2
fi

bin/celltaminate run   --input-list examples/lung_nanopore/input_list.txt   --metadata examples/lung_nanopore/metadata.tsv   --outdir examples/lung_nanopore/results   --tax-level S   "${REF_ARGS[@]}"
