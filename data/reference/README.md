# Reference background data

`refined_cell.lines.example.tsv` is a small example subset of the reference background table. It is included only so that the repository contains a lightweight runnable example.

The full `refined_cell.lines.tsv` reference background table is available from Zenodo:

```text
DOI: 10.5281/zenodo.20560460
```

Download it locally before real analysis:

```bash
curl -L -o data/reference/refined_cell.lines.tsv \
  "https://zenodo.org/records/20560460/files/refined_cell.lines.tsv?download=1"
```

The repository `.gitignore` excludes `data/reference/refined_cell.lines.tsv` by default.
