# Reference data

Celltaminate uses a reference background dataset derived from sterile human cell-line sequencing experiments that contains 1,040 unique sample identifiers. This dataset is the same that we used in SAHMI (https://github.com/sjdlabgroup/SAHMI)

Format of the reference dataset file:

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

or provide a custom reference background dataset at run time:

```bash
bin/celltaminate run \
  --input-list input_list.txt \
  --outdir results \
  --ref-bg-tsv /path/to/refined_cell.lines.tsv
```
