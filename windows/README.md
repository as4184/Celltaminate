# Windows executable launcher

`CelltaminateRunner.exe` is a Windows command-line launcher for the full workflow:

```text
FASTQ / FASTQ.GZ -> Kraken2 report files -> Celltaminate tables and plots
```

The executable bundles the Python orchestration code. It does **not** bundle the large Kraken2 database, Kraken2 itself, R, or the Celltaminate reference background table. Those files/tools remain external because they are large and platform-specific.

## Build the executable

On Windows, from the repository root:

```bat
windows\build_windows_exe.bat
```

This creates:

```text
dist\CelltaminateRunner.exe
```

Test it:

```bat
dist\CelltaminateRunner.exe --help
```

## FASTQ samplesheet

Create `samples.tsv`:

```text
sample	r1	r2	group	sample_type	is_control
Sample1	C:\data\Sample1_R1.fastq.gz	C:\data\Sample1_R2.fastq.gz	Case	Illumina	FALSE
Sample2	C:\data\Sample2.fastq.gz		Case	Nanopore	FALSE
```

For single-end or Nanopore data, leave `r2` blank.

## Run the full workflow

Example using tools available on the Windows PATH:

```bat
dist\CelltaminateRunner.exe ^
  --repo-root C:\path\to\Celltaminate ^
  --samplesheet C:\path\to\samples.tsv ^
  --kraken2-db D:\kraken2_db ^
  --outdir C:\path\to\celltaminate_full_run ^
  --threads 8 ^
  --rscript "C:\Program Files\R\R-4.4.0\bin\Rscript.exe" ^
  --ref-bg-tsv C:\path\to\refined_cell.lines.tsv
```

Outputs:

```text
celltaminate_full_run\
├── kraken2_reports\
├── kraken2_output\
├── workflow_inputs\
└── celltaminate_results\
```

## Practical recommendation

For large FASTQ files and Kraken2 databases, Linux, WSL2, HPC, or Docker is usually more reliable than native Windows. The Windows executable is intended as a launcher/orchestrator, not as a replacement for Kraken2, R, or the Kraken2 database.
