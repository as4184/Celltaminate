# Troubleshooting

## The reference background is not found

Download it to the default location:

```bash
bin/celltaminate download-reference
```

or provide `--ref-bg-tsv /path/to/refined_cell.lines.tsv`.

## Metadata sample names do not match

The `sample` column in `metadata.tsv` must match the basename of the corresponding Kraken report exactly.

## Unique k-mer information is missing

Kraken2 6-column reports can be processed, but 8-column Kraken2 reports with minimizer data or KrakenUniq reports are preferred because unique k-mer/minimizer support is used by the fitted model.

## The Shiny application cannot find the curated panels

Launch the application through `bin/celltaminate app`. The launcher and application resolve the bundled files from the repository automatically.

## R packages are missing

Create the supplied conda environment again:

```bash
conda env create -f environment.yml
conda activate celltaminate
```
