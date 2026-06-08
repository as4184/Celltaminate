# Input formats

## Kraken report files

Celltaminate starts from taxonomic summary reports. It accepts:

1. Kraken2 6-column report files;
2. Kraken2 8-column report files;
3. KrakenUniq 8-column report files.

The workflow uses the indentation structure of the report to reconstruct lineage relationships. For this reason, the original taxon-name indentation should be preserved.

## Input list

`input_list.txt` should contain one report path per line:

```text
sample1.kraken.report.txt
sample2.kraken.report.txt
sample3.kraken.report.txt
```

Relative paths are interpreted relative to the location of the input list by the wrapper.

## Metadata table

`metadata.tsv` must be tab-separated:

```text
sample	group	sample_type	is_control
sample1.kraken.report.txt	Case	Nanopore	FALSE
sample2.kraken.report.txt	Control	Nanopore	TRUE
```

Required columns:

- `sample`: report basename;
- `group`: biological or study group;
- `sample_type`: sequencing or sample description;
- `is_control`: `TRUE` or `FALSE`.


## FASTQ samplesheet for the full workflow

If using `scripts/celltaminate_full_workflow.py` or `bin/celltaminate from-fastq`, provide a TSV/CSV samplesheet with:

```text
sample	r1	r2	group	sample_type	is_control
Sample1	/path/to/Sample1_R1.fastq.gz	/path/to/Sample1_R2.fastq.gz	Case	Illumina	FALSE
Sample2	/path/to/Sample2.fastq.gz		Case	Nanopore	FALSE
```

Required columns:

- `sample`: sample identifier with no path separators;
- `r1`: first FASTQ/FASTQ.GZ file, or the only FASTQ file for single-end/Nanopore data.

Optional columns:

- `r2`: second FASTQ/FASTQ.GZ file for paired-end data;
- `group`: biological group;
- `sample_type`: sample or sequencing type;
- `is_control`: `TRUE` or `FALSE`;
- `report`: existing Kraken2 report path when `--skip-kraken2` is used.
