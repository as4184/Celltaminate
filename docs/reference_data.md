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

Download it before running an analysis:

```bash
bin/celltaminate download-reference --out data/reference/refined_cell.lines.tsv
```

Alternatively:

```bash
mkdir -p data/reference
curl -L -o data/reference/refined_cell.lines.tsv \
  "https://zenodo.org/records/20560460/files/refined_cell.lines.tsv?download=1"
```

Then pass the file like this:

```bash
bin/celltaminate run \
  --input-list input_list.txt \
  --metadata metadata.tsv \
  --outdir results \
  --ref-bg-tsv data/reference/refined_cell.lines.tsv
```
