# Reference data

Celltaminate can use a reference background table built from sterile human cell-line microbial profiles. The expected columns are:

```text
sample	rank	taxid	name	reads	min	uniq	rpm	rpmm
```

## Full reference table

The full `refined_cell.lines.tsv` reference background table is available from Zenodo:

```text
DOI: 10.5281/zenodo.20560460
```

Download it before real analysis:

```bash
bin/celltaminate download-reference --out data/reference/refined_cell.lines.tsv
```

Alternatively:

```bash
mkdir -p data/reference
curl -L -o data/reference/refined_cell.lines.tsv \
  "https://zenodo.org/records/20560460/files/refined_cell.lines.tsv?download=1"
```

Then pass the file explicitly:

```bash
bin/celltaminate run \
  --input-list input_list.txt \
  --metadata metadata.tsv \
  --outdir results \
  --ref-bg-tsv data/reference/refined_cell.lines.tsv
```

## Example reference table

This repository includes a small example file:

```text
data/reference/refined_cell.lines.example.tsv
```

This file is only for testing the command-line workflow. It is not a replacement for the full reference background table.

## Storage note

The full reference table should usually stay outside normal Git tracking. The repository `.gitignore` excludes `data/reference/refined_cell.lines.tsv` by default so users can download it locally without accidentally committing the large table.
