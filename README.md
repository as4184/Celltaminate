# Celltaminate

**Celltaminate** is an R-based command-line and Shiny application for contamination-aware prioritization of microbial taxa from host-dominated sequencing data.

Celltaminate starts from Kraken2 or KrakenUniq taxonomic reports and integrates quantitative sequence evidence, taxonomic ambiguity, sterile-reference background, decontamination, curated clinical-pathogen and kitome/background evidence, and selected interaction terms into a continuous score from 0 to 100.

Lower scores indicate stronger support for prioritization. The default operating rule is:

```text
score <= 1    prioritized
score > 1     not prioritized
```

The score is a logistic mapping of the fitted regularized model and is not a calibrated probability.

## Installation

```bash
git clone https://github.com/as4184/Celltaminate.git
cd Celltaminate
conda env create -f environment.yml
conda activate celltaminate
bin/celltaminate --help
```

Download the sterile-reference background table:

```bash
bin/celltaminate download-reference
```

Reference DOI: `10.5281/zenodo.20560460`

## Quick start

Run the bundled example:

```bash
bash examples/lung_nanopore/run_example.sh
```

Run Celltaminate on your own reports:

```bash
bin/celltaminate run \
  --reports sample1.kraken.report.txt sample2.kraken.report.txt \
  --outdir celltaminate_results
```

Optional metadata can be supplied with `--metadata metadata.tsv`. Metadata contain `sample`, `group`, and `sample_type`.

## Shiny application

```bash
bin/celltaminate app
```

## Input

Celltaminate accepts Kraken2 6-column reports, Kraken2 8-column reports with minimizer data, and KrakenUniq 8-column reports. Eight-column reports are preferred because unique k-mer/minimizer support is part of the fitted scoring model.

## Output

The primary per-sample tables report taxon-level evidence, the Celltaminate score, and the binary call `Prioritized` or `Not prioritized`. Cleaned reports retain prioritized taxa by default.

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
