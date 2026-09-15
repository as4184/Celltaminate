#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
EXAMPLE_DIR="${ROOT_DIR}/examples/lung_nanopore"
OUTDIR="${EXAMPLE_DIR}/celltaminate_example_output"

if [[ ! -f "${ROOT_DIR}/data/reference/refined_cell.lines.tsv" ]]; then
  echo "Reference background table not found. Downloading it first."
  "${ROOT_DIR}/bin/celltaminate" download-reference
fi

"${ROOT_DIR}/bin/celltaminate" run \
  --input-list "${EXAMPLE_DIR}/input_list.txt" \
  --metadata "${EXAMPLE_DIR}/metadata.tsv" \
  --outdir "${OUTDIR}"

echo "Example completed: ${OUTDIR}"
