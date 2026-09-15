# Celltaminate

**Celltaminate** (**C**ontamination **E**valuation for **L**ow-Level **T**axon-**A**bundance in **M**icrobial **I**nference with **N**oise-**A**ware **T**riage **E**stimate) is an R-based command-line and Shiny application for contamination-aware prioritization of microbial taxa from host-dominated sequencing data.

Celltaminate starts from Kraken2 or KrakenUniq taxonomic reports and evaluates whether each microbial taxon is supported by quantitative sequence evidence or is more consistent with recurrent background, contamination, or taxonomic ambiguity.

## Scoring framework

The deployed score integrates retained evidence from:

- clade-read and RPMM support;
- unique k-mer support and taxonomic ambiguity;
- sterile-reference prevalence, abundance, and abundance quantiles;
- background-subtracted abundance and decontamination ratio;
- curated clinically important pathogen and kitome/background panels;
- selected interaction terms between these evidence layers.

The score ranges from 0 to 100. Lower scores indicate stronger support for prioritization. The default prioritization threshold is **score <= 1**. The software also retains an intermediate uncertain range and a high-score background category for users who want three-tier calls.

## Installation

Clone the repository and create the conda environment:

```bash
git clone https://github.com/as4184/Celltaminate.git
cd Celltaminate
conda env create -f environment.yml
conda activate celltaminate
bin/celltaminate --help
```

The full sterile-reference background table is distributed separately through Zenodo:

```text
DOI: 10.5281/zenodo.20560460
```

Download it into the expected repository location:

```bash
bin/celltaminate download-reference
```

## Quick start

Run the bundled Kraken-report example:

```bash
bash examples/lung_nanopore/run_example.sh
```

Run Celltaminate on your own reports:

```bash
bin/celltaminate run \
  --reports sample1.kraken.report.txt sample2.kraken.report.txt \
  --metadata metadata.tsv \
  --outdir celltaminate_results
```

Or provide one report path per line:

```bash
bin/celltaminate run \
  --input-list input_list.txt \
  --metadata metadata.tsv \
  --outdir celltaminate_results
```

The bundled kitome/background and clinical-pathogen panels are loaded automatically. If `data/reference/refined_cell.lines.tsv` is present, the reference background table is also loaded automatically.

## Shiny application

Launch the interactive application with:

```bash
bin/celltaminate app
```

The Shiny interface accepts Kraken2 or KrakenUniq reports and provides interactive score summaries, taxon tables, abundance views, and contextual taxon information.

## Input

Celltaminate accepts:

- Kraken2 6-column reports;
- Kraken2 8-column reports generated with minimizer data;
- KrakenUniq 8-column reports.

Eight-column reports are preferred because unique k-mer/minimizer support is part of the fitted scoring model.

Optional metadata must be tab-separated and contain:

```text
sample	group	sample_type	is_control
sample1.kraken.report.txt	Group 1	Clinical	FALSE
```

The `sample` value must match the report basename.

## Output

A standard run creates:

```text
celltaminate_results/
├── tables/
│   ├── all_taxa_table.tsv
│   ├── per_sample/
│   └── summaries/
├── plots/
├── cleaned_reports/
├── cleaned_reanalysis/
├── run_config/
└── run_summary.tsv
```

The per-sample taxa tables contain the main taxon-level evidence, Celltaminate score, and call.

## Containers

A Dockerfile and Apptainer definition are included for containerized execution.

```bash
docker build -t celltaminate:latest .
docker run --rm celltaminate:latest --help
```

## Documentation

- [Installation](docs/installation.md)
- [Quick start](docs/quickstart.md)
- [Input formats](docs/input_formats.md)
- [Output files](docs/output_files.md)
- [Reference data](docs/reference_data.md)
- [Score interpretation](docs/scoring.md)
- [Shiny application](docs/shiny_app.md)
- [Troubleshooting](docs/troubleshooting.md)

## Citation metadata

See [`CITATION.cff`](CITATION.cff).

## License

See [`LICENSE.md`](LICENSE.md).
