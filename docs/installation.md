# Installation

## Conda

```bash
git clone https://github.com/as4184/Celltaminate.git
cd Celltaminate
conda env create -f environment.yml
conda activate celltaminate
bin/celltaminate --help
```

## Reference background

Download the full reference background table after installation:

```bash
bin/celltaminate download-reference
```

This writes:

```text
data/reference/refined_cell.lines.tsv
```

Reference DOI: `10.5281/zenodo.20560460`.

## Docker

```bash
docker build -t celltaminate:latest .
docker run --rm celltaminate:latest --help
```

For local input files, mount a working directory and the reference table as needed.

## Apptainer

```bash
apptainer build celltaminate.sif apptainer.def
apptainer run celltaminate.sif --help
```

## R packages

The command-line workflow uses `dplyr`, `ggplot2`, `ggrepel`, `fmsb`, `writexl`, `tibble`, and `tidyr`. The Shiny application additionally uses `shiny`, `shinyWidgets`, `DT`, `httr`, `jsonlite`, `zip`, and `xml2`.
