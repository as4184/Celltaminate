# Reference data

Celltaminate uses a processed reference background derived from sterile human cell-line sequencing experiments. The source resource contains 2,491 sterile experiments, and the processed reference table used by the current implementation contains 1,040 unique sample identifiers.

The expected columns are:

```text
sample	rank	taxid	name	reads	min	uniq	rpm	rpmm
```

The full `refined_cell.lines.tsv` table is distributed through Zenodo:

```text
DOI: 10.5281/zenodo.20560460
```

Download it to the default location with:

```bash
bin/celltaminate download-reference
```

or provide an explicit file at run time:

```bash
bin/celltaminate run \
  --input-list input_list.txt \
  --outdir results \
  --ref-bg-tsv /path/to/refined_cell.lines.tsv
```
