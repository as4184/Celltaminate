
## Unreleased

- Added FASTQ-to-Kraken2-to-Celltaminate workflow launcher.
- Added documentation for upstream Kraken2 report generation from FASTQ files.
- Added Windows PyInstaller build script for a `CelltaminateRunner.exe` launcher.
- Added `bin/celltaminate from-fastq` delegation command.

# Changelog

## 0.1.0

- Added command-line wrapper `bin/celltaminate`.
- Added Shiny app entry point under `app/`.
- Added lung Nanopore example input.
- Added curated kitome/background and clinically important pathogen panels.
- Added documentation, container definitions, and smoke-test workflow.

## Unreleased

- Updated Shiny app source.
- Added Zenodo reference background table DOI: `10.5281/zenodo.20560460`.
- Added `bin/celltaminate download-reference` helper.

- Revised the FASTQ workflow to run Kraken2 with explicit FASTQ input handling, optional compression flags, and minimizer-enhanced reports by default.
