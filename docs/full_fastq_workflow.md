# FASTQ to Celltaminate workflow

Celltaminate itself starts from Kraken2 or KrakenUniq taxonomic report files. If you are starting from raw FASTQ files, first classify reads with Kraken2, then pass the generated report files into Celltaminate.

## Step 1: Prepare a Kraken2 database

Use a Kraken2 database appropriate for the experiment. For host-dominated sequencing, the database should include the host genome and microbial genomes whenever possible, so that host-like reads are not forced into microbial taxa.

Example database path:

```bash
KRAKEN2_DB=/path/to/kraken2_database
```

The database directory should contain Kraken2 database files such as:

```text
hash.k2d
opts.k2d
taxo.k2d
```

## Step 2: Create a FASTQ samplesheet

Create `samples.tsv`:

```text
sample	r1	r2	group	sample_type	is_control
Sample1	/path/to/Sample1_R1.fastq.gz	/path/to/Sample1_R2.fastq.gz	Case	Illumina	FALSE
Sample2	/path/to/Sample2.fastq.gz		Case	Nanopore	FALSE
Control1	/path/to/Control1_R1.fastq.gz	/path/to/Control1_R2.fastq.gz	Control	Illumina	TRUE
```

Required columns:

- `sample`: sample identifier with no `/` or `\` characters;
- `r1`: FASTQ or FASTQ.GZ path.

Optional columns:

- `r2`: second FASTQ for paired-end Illumina data; leave blank for single-end/Nanopore data;
- `group`: case/control/cohort label;
- `sample_type`: sequencing/sample description;
- `is_control`: `TRUE` or `FALSE`;
- `report`: existing Kraken2 report path, used only with `--skip-kraken2`.

## Step 3A: Run the complete workflow with one command

```bash
python scripts/celltaminate_full_workflow.py \
  --samplesheet samples.tsv \
  --kraken2-db "$KRAKEN2_DB" \
  --outdir celltaminate_full_run \
  --threads 16 \
  --input-compression auto \
  --report-minimizer-data \
  --ref-bg-tsv data/reference/refined_cell.lines.tsv
```

This writes:

```text
celltaminate_full_run/
├── kraken2_reports/
├── kraken2_output/
├── workflow_inputs/
└── celltaminate_results/
```

## Step 3B: Run Kraken2 first, then run Celltaminate manually

For paired-end Illumina data:

```bash
mkdir -p kraken2_reports kraken2_output

kraken2 \
  --db "$KRAKEN2_DB" \
  --threads 16 \
  --fastq-input \
  --gzip-compressed \
  --report kraken2_reports/sample.kraken.report.txt \
  --report-minimizer-data \
  --output kraken2_output/sample.kraken.output.txt \
  --paired sample_R1.fastq.gz sample_R2.fastq.gz
```

For single-end or Nanopore data:

```bash
mkdir -p kraken2_reports kraken2_output

kraken2 \
  --db "$KRAKEN2_DB" \
  --threads 16 \
  --fastq-input \
  --gzip-compressed \
  --report kraken2_reports/sample.kraken.report.txt \
  --report-minimizer-data \
  --output kraken2_output/sample.kraken.output.txt \
  sample.fastq.gz
```

Then create `input_list.txt`:

```text
kraken2_reports/sample.kraken.report.txt
```

Create `metadata.tsv`:

```text
sample	group	sample_type	is_control
sample.kraken.report.txt	Case	Illumina	FALSE
```

Run Celltaminate:

```bash
bin/celltaminate run \
  --input-list input_list.txt \
  --metadata metadata.tsv \
  --outdir celltaminate_results \
  --tax-level S \
  --ref-bg-tsv data/reference/refined_cell.lines.tsv
```

## Step 4: Check the main output table

```bash
find celltaminate_full_run/celltaminate_results/tables/per_sample -type f -name '*_taxa_table.tsv'
```

The per-sample taxa tables contain the Celltaminate score, call, reference background values, kitome/clinical panel annotations, and cleaned abundance values.

## Notes

- Celltaminate uses Kraken report files, not the raw Kraken read-classification output. The `--output` file is useful for audit/troubleshooting but is not the main Celltaminate input.
- Use `--report-minimizer-data` when running Kraken2. Celltaminate can parse 6-column Kraken2 reports, but the 8-column minimizer-enhanced report preserves k-mer/minimizer support and is preferred for scoring.
- Use `--fastq-input` explicitly for FASTQ files and add `--gzip-compressed` or `--bzip2-compressed` for compressed input. The full workflow script can infer this with `--input-compression auto`.
- Keep the original Kraken report indentation intact.
- For low-biomass host-dominated sequencing, use a database that includes the host genome when possible.
- Matched negative controls can be listed in the same samplesheet with `is_control = TRUE`.
