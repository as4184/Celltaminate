# Installation

## Conda installation

```bash
git clone https://github.com/<organization>/Celltaminate.git
cd Celltaminate
conda env create -f environment.yml
conda activate celltaminate
bin/celltaminate --help
```

## Docker

```bash
docker build -t celltaminate:latest .
docker run --rm celltaminate:latest --help
```

To run local data with Docker:

```bash
docker run --rm -v "$PWD":/work -w /work celltaminate:latest run       --input-list examples/lung_nanopore/input_list.txt       --metadata examples/lung_nanopore/metadata.tsv       --outdir examples/lung_nanopore/results       --use-example-reference
```

## Apptainer

```bash
apptainer build celltaminate.sif apptainer.def
apptainer run celltaminate.sif --help
```

## Required R packages

The command-line workflow uses `dplyr`, `ggplot2`, `ggrepel`, `fmsb`, `writexl`, `tibble`, and `tidyr`. The Shiny app additionally uses `shiny`, `shinyWidgets`, `DT`, `httr`, `jsonlite`, `zip`, and `xml2`.

## Reference background table

The full `refined_cell.lines.tsv` reference background table is available from Zenodo, DOI `10.5281/zenodo.20560460`. Download it locally and provide it with `--ref-bg-tsv`.


## Windows executable launcher

A Windows `.exe` launcher can be built with PyInstaller:

```bat
windows\build_windows_exe.bat
```

The generated `dist\CelltaminateRunner.exe` runs the workflow orchestration code. Kraken2, the Kraken2 database, R/Rscript, and the reference background table must still be installed or provided separately. See [Windows executable launcher](../windows/README.md).
