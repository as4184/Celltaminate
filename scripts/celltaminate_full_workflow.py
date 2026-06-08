#!/usr/bin/env python3
"""
Run the full Celltaminate upstream workflow:

FASTQ/FASTQ.GZ files -> Kraken2 taxonomic reports -> Celltaminate R workflow.

The script is intentionally dependency-light. It orchestrates external tools
(Kraken2 and Rscript) rather than reimplementing them. It can also be compiled
into a Windows .exe launcher with PyInstaller.
"""

from __future__ import annotations

import argparse
import csv
import os
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

VERSION = "0.1.0"

REQUIRED_SAMPLE_COLUMNS = {"sample", "r1"}
OPTIONAL_METADATA_COLUMNS = ["group", "sample_type", "is_control"]


def script_path() -> Path:
    return Path(__file__).resolve()


def guess_repo_root(user_value: Optional[str] = None) -> Path:
    if user_value:
        return Path(user_value).expanduser().resolve()

    env_value = os.environ.get("CELLTAMINATE_HOME", "").strip()
    if env_value:
        return Path(env_value).expanduser().resolve()

    here = script_path()
    candidates = [
        Path.cwd(),
        here.parent.parent,
        here.parent,
    ]
    for candidate in candidates:
        candidate = candidate.resolve()
        if (candidate / "scripts" / "Celltaminate.R").exists():
            return candidate

    return Path.cwd().resolve()


def detect_delimiter(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix in {".csv"}:
        return ","
    if suffix in {".tsv", ".txt"}:
        return "\t"

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        sample = handle.read(4096)
    try:
        return csv.Sniffer().sniff(sample, delimiters=["\t", ","]).delimiter
    except csv.Error:
        return "\t"


def norm(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()


def read_samplesheet(path: Path) -> List[Dict[str, str]]:
    if not path.exists() or not path.is_file():
        raise SystemExit(f"ERROR: samplesheet not found: {path}")

    delimiter = detect_delimiter(path)
    rows: List[Dict[str, str]] = []

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None:
            raise SystemExit(f"ERROR: samplesheet has no header: {path}")

        fieldnames = {norm(x) for x in reader.fieldnames}
        missing = REQUIRED_SAMPLE_COLUMNS.difference(fieldnames)
        if missing:
            raise SystemExit(
                "ERROR: samplesheet is missing required columns: "
                + ", ".join(sorted(missing))
                + "\nRequired columns: sample, r1\nOptional columns: r2, group, sample_type, is_control, report"
            )

        for i, raw in enumerate(reader, start=2):
            row = {norm(k): norm(v) for k, v in raw.items() if k is not None}
            if not any(row.values()):
                continue

            sample = row.get("sample", "")
            if not sample:
                raise SystemExit(f"ERROR: blank sample name in samplesheet line {i}")
            if any(x in sample for x in ["/", "\\"]):
                raise SystemExit(f"ERROR: sample name must not contain path separators in line {i}: {sample}")

            r1 = row.get("r1", "")
            if not r1:
                raise SystemExit(f"ERROR: blank r1 path for sample {sample} in line {i}")

            rows.append(row)

    if not rows:
        raise SystemExit(f"ERROR: samplesheet contains no samples: {path}")

    seen = set()
    duplicates = []
    for row in rows:
        sample = row["sample"]
        if sample in seen:
            duplicates.append(sample)
        seen.add(sample)
    if duplicates:
        raise SystemExit("ERROR: duplicate sample names in samplesheet: " + ", ".join(sorted(set(duplicates))))

    return rows


def resolve_path(value: str, base: Path) -> Path:
    p = Path(value).expanduser()
    if not p.is_absolute():
        p = (base / p).resolve()
    else:
        p = p.resolve()
    return p


def validate_fastqs(rows: Iterable[Dict[str, str]], samplesheet_dir: Path, skip_kraken2: bool) -> None:
    if skip_kraken2:
        return

    missing: List[str] = []
    for row in rows:
        for col in ["r1", "r2"]:
            value = row.get(col, "")
            if not value:
                continue
            path = resolve_path(value, samplesheet_dir)
            if not path.exists() or not path.is_file():
                missing.append(f"{row['sample']}:{col} -> {path}")

    if missing:
        raise SystemExit("ERROR: missing FASTQ files:\n" + "\n".join(missing))


def split_command(value: str) -> List[str]:
    value = norm(value)
    if not value:
        raise SystemExit("ERROR: empty command value")
    return shlex.split(value, posix=(os.name != "nt"))


def quote_command(cmd: Sequence[str]) -> str:
    if os.name == "nt":
        return subprocess.list2cmdline([str(x) for x in cmd])
    return " ".join(shlex.quote(str(x)) for x in cmd)


def run_command(cmd: Sequence[str], dry_run: bool = False) -> None:
    print("\nRunning:")
    print(quote_command(cmd))
    if dry_run:
        return
    rc = subprocess.call([str(x) for x in cmd])
    if rc != 0:
        raise SystemExit(rc)


def kraken2_report_path(report_dir: Path, sample: str) -> Path:
    return report_dir / f"{sample}.kraken.report.txt"


def kraken2_output_path(output_dir: Path, sample: str) -> Path:
    return output_dir / f"{sample}.kraken.output.txt"


def detect_compression_for_kraken2(paths: Sequence[Path], mode: str) -> Optional[str]:
    """Return the Kraken2 compression flag for the selected input mode.

    Kraken2 can read compressed inputs, but the official manual documents
    explicit switches for gzip and bzip2 input. The workflow therefore adds
    these flags in auto mode when extensions indicate compressed FASTQ files.
    """
    mode = norm(mode).lower()
    if mode in {"none", "uncompressed"}:
        return None
    if mode in {"gzip", "gz"}:
        return "--gzip-compressed"
    if mode in {"bzip2", "bz2"}:
        return "--bzip2-compressed"
    if mode != "auto":
        raise SystemExit(f"ERROR: unsupported --input-compression value: {mode}")

    suffixes = [str(path).lower() for path in paths]
    has_gz = any(x.endswith(".gz") or x.endswith(".gzip") for x in suffixes)
    has_bz2 = any(x.endswith(".bz2") or x.endswith(".bzip2") for x in suffixes)

    if has_gz and has_bz2:
        raise SystemExit("ERROR: mixed gzip and bzip2 inputs are not supported in one Kraken2 command")
    if has_gz:
        return "--gzip-compressed"
    if has_bz2:
        return "--bzip2-compressed"
    return None


def run_kraken2_for_sample(
    row: Dict[str, str],
    samplesheet_dir: Path,
    report_dir: Path,
    output_dir: Path,
    kraken2_command: str,
    kraken2_db: Path,
    threads: int,
    confidence: Optional[float],
    minimum_hit_groups: Optional[int],
    input_compression: str,
    report_minimizer_data: bool,
    use_names: bool,
    extra_args: str,
    skip_existing: bool,
    dry_run: bool,
) -> Path:
    sample = row["sample"]
    report = kraken2_report_path(report_dir, sample)
    output = kraken2_output_path(output_dir, sample)

    if skip_existing and report.exists() and report.stat().st_size > 0:
        print(f"Skipping Kraken2 for {sample}; existing report found: {report}")
        return report

    r1 = resolve_path(row["r1"], samplesheet_dir)
    r2_value = row.get("r2", "")
    paired = bool(r2_value)
    r2 = resolve_path(r2_value, samplesheet_dir) if paired else None

    input_paths = [r1] + ([r2] if paired and r2 is not None else [])
    compression_flag = detect_compression_for_kraken2(input_paths, input_compression)

    cmd = split_command(kraken2_command)
    cmd += [
        "--db",
        str(kraken2_db),
        "--threads",
        str(threads),
        "--fastq-input",
        "--report",
        str(report),
    ]

    # Celltaminate can read ordinary 6-column Kraken2 reports, but the
    # minimizer-enhanced 8-column report is strongly preferred because it
    # preserves the distinct-minimizer/k-mer support used by the scoring model.
    if report_minimizer_data:
        cmd += ["--report-minimizer-data"]

    cmd += ["--output", str(output)]

    if use_names:
        cmd += ["--use-names"]
    if compression_flag:
        cmd += [compression_flag]
    if confidence is not None:
        cmd += ["--confidence", str(confidence)]
    if minimum_hit_groups is not None:
        cmd += ["--minimum-hit-groups", str(minimum_hit_groups)]
    if extra_args:
        cmd += shlex.split(extra_args, posix=(os.name != "nt"))

    if paired:
        cmd += ["--paired", str(r1), str(r2)]
    else:
        cmd += [str(r1)]

    run_command(cmd, dry_run=dry_run)

    if not dry_run and (not report.exists() or report.stat().st_size == 0):
        raise SystemExit(f"ERROR: Kraken2 did not create a non-empty report for {sample}: {report}")

    return report


def write_input_list(path: Path, reports: Iterable[Path]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for report in reports:
            handle.write(str(report.resolve()) + "\n")


def write_metadata(path: Path, rows: Iterable[Dict[str, str]], report_paths: Dict[str, Path]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample", "group", "sample_type", "is_control"], delimiter="\t")
        writer.writeheader()
        for row in rows:
            sample = row["sample"]
            writer.writerow(
                {
                    "sample": report_paths[sample].name,
                    "group": row.get("group", "Group 1") or "Group 1",
                    "sample_type": row.get("sample_type", "FASTQ") or "FASTQ",
                    "is_control": row.get("is_control", "FALSE") or "FALSE",
                }
            )


def first_existing(paths: Iterable[Path]) -> Optional[Path]:
    for path in paths:
        if path.exists() and path.is_file():
            return path.resolve()
    return None


def bool_to_r(value: bool) -> str:
    return "TRUE" if value else "FALSE"


def build_celltaminate_command(args: argparse.Namespace, repo: Path, input_list: Path, metadata: Path) -> List[str]:
    celltaminate_r = Path(args.celltaminate_r).expanduser().resolve() if args.celltaminate_r else repo / "scripts" / "Celltaminate.R"
    if not celltaminate_r.exists():
        raise SystemExit(f"ERROR: Celltaminate.R not found: {celltaminate_r}")

    cmd = split_command(args.rscript)
    cmd += [
        "--vanilla",
        str(celltaminate_r),
        "--input_list",
        str(input_list),
        "--metadata_tsv",
        str(metadata),
        "--output_dir",
        str(Path(args.outdir).expanduser().resolve() / "celltaminate_results"),
        "--tax_level",
        args.tax_level,
        "--host_species",
        args.host_species,
        "--include_protozoa",
        bool_to_r(args.include_protozoa),
        "--collapse_species",
        bool_to_r(args.collapse_species),
        "--collapse_min_genus_reads",
        str(args.collapse_min_genus_reads),
        "--collapse_top_frac",
        str(args.collapse_top_frac),
        "--fp_falsepos_cutoff",
        str(args.fp_falsepos_cutoff),
        "--fp_true_cutoff",
        str(args.fp_true_cutoff),
        "--fp_aggressiveness",
        str(args.fp_aggressiveness),
        "--quick_remove",
        args.quick_remove,
        "--keep_descendants",
        bool_to_r(args.keep_descendants),
        "--show_fp_breakdown",
        bool_to_r(args.show_fp_breakdown),
        "--save_cleaned_reanalysis",
        bool_to_r(args.save_cleaned_reanalysis),
    ]

    ref_bg = Path(args.ref_bg_tsv).expanduser().resolve() if args.ref_bg_tsv else first_existing([
        repo / "data" / "reference" / "refined_cell.lines.tsv",
        Path.cwd() / "data" / "reference" / "refined_cell.lines.tsv",
    ])
    if ref_bg is not None:
        cmd += ["--ref_bg_tsv", str(ref_bg)]
    else:
        print("WARNING: no reference background TSV was supplied or found. Celltaminate will run without the cell-line reference background layer.")

    kitome = Path(args.kitome_blacklist_tsv).expanduser().resolve() if args.kitome_blacklist_tsv else first_existing([
        repo / "data" / "panels" / "kitome_and_background_blacklist.tsv",
        Path.cwd() / "data" / "panels" / "kitome_and_background_blacklist.tsv",
    ])
    if kitome is not None:
        cmd += ["--kitome_blacklist_tsv", str(kitome)]

    clinical = Path(args.clinical_panel_tsv).expanduser().resolve() if args.clinical_panel_tsv else first_existing([
        repo / "data" / "panels" / "clinically_important_pathogens.tsv",
        Path.cwd() / "data" / "panels" / "clinically_important_pathogens.tsv",
    ])
    if clinical is not None:
        cmd += ["--clinical_panel_tsv", str(clinical)]

    if args.protected_taxa_tsv:
        cmd += ["--protected_taxa_tsv", str(Path(args.protected_taxa_tsv).expanduser().resolve())]
    if args.host_taxids_manual:
        cmd += ["--host_taxids_manual", args.host_taxids_manual]

    return cmd


def cmd_full_workflow(args: argparse.Namespace) -> int:
    repo = guess_repo_root(args.repo_root)
    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    samplesheet = Path(args.samplesheet).expanduser().resolve()
    samplesheet_dir = samplesheet.parent
    rows = read_samplesheet(samplesheet)
    validate_fastqs(rows, samplesheet_dir, skip_kraken2=args.skip_kraken2)

    report_dir = outdir / "kraken2_reports"
    output_dir = outdir / "kraken2_output"
    workflow_dir = outdir / "workflow_inputs"
    report_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    workflow_dir.mkdir(parents=True, exist_ok=True)

    report_paths: Dict[str, Path] = {}

    if args.skip_kraken2:
        for row in rows:
            sample = row["sample"]
            if row.get("report", ""):
                report = resolve_path(row["report"], samplesheet_dir)
            else:
                report = kraken2_report_path(report_dir, sample)
            if not report.exists() or not report.is_file():
                raise SystemExit(f"ERROR: --skip-kraken2 was used but report is missing for {sample}: {report}")
            report_paths[sample] = report.resolve()
    else:
        if not args.kraken2_db:
            raise SystemExit("ERROR: --kraken2-db is required unless --skip-kraken2 is used")
        kraken2_db = Path(args.kraken2_db).expanduser().resolve()
        if not kraken2_db.exists() or not kraken2_db.is_dir():
            raise SystemExit(f"ERROR: Kraken2 database directory not found: {kraken2_db}")

        for row in rows:
            report = run_kraken2_for_sample(
                row=row,
                samplesheet_dir=samplesheet_dir,
                report_dir=report_dir,
                output_dir=output_dir,
                kraken2_command=args.kraken2,
                kraken2_db=kraken2_db,
                threads=args.threads,
                confidence=args.confidence,
                minimum_hit_groups=args.minimum_hit_groups,
                input_compression=args.input_compression,
                report_minimizer_data=args.report_minimizer_data,
                use_names=args.use_names,
                extra_args=args.extra_kraken2_args,
                skip_existing=args.skip_existing_reports,
                dry_run=args.dry_run,
            )
            report_paths[row["sample"]] = report.resolve()

    input_list = workflow_dir / "input_list.txt"
    metadata = workflow_dir / "metadata.tsv"
    write_input_list(input_list, [report_paths[row["sample"]] for row in rows])
    write_metadata(metadata, rows, report_paths)

    print(f"\nWrote Celltaminate input list: {input_list}")
    print(f"Wrote Celltaminate metadata: {metadata}")

    if args.skip_celltaminate:
        print("Skipping Celltaminate because --skip-celltaminate was used.")
        return 0

    cmd = build_celltaminate_command(args, repo=repo, input_list=input_list, metadata=metadata)
    run_command(cmd, dry_run=args.dry_run)
    return 0


def add_bool_pair(parser: argparse.ArgumentParser, name: str, default: bool, help_text: str) -> None:
    dest = name.lstrip("-").replace("-", "_")
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument(name, dest=dest, action="store_true", help=help_text)
    group.add_argument("--no-" + name.lstrip("-"), dest=dest, action="store_false")
    parser.set_defaults(**{dest: default})


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="celltaminate_full_workflow",
        description="Run FASTQ -> Kraken2 reports -> Celltaminate.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"Celltaminate full workflow {VERSION}")

    parser.add_argument("--samplesheet", required=True, help="TSV/CSV with columns: sample, r1, optional r2, group, sample_type, is_control, report")
    parser.add_argument("--outdir", required=True, help="Output directory for Kraken2 reports and Celltaminate results")
    parser.add_argument("--repo-root", help="Celltaminate repository root. Defaults to $CELLTAMINATE_HOME or current/source location")

    parser.add_argument("--kraken2", default="kraken2", help="Kraken2 command or command prefix, e.g. 'kraken2' or 'wsl kraken2'")
    parser.add_argument("--kraken2-db", help="Kraken2 database directory")
    parser.add_argument("--threads", type=int, default=8, help="Threads passed to Kraken2")
    parser.add_argument("--confidence", type=float, help="Optional Kraken2 --confidence value")
    parser.add_argument("--minimum-hit-groups", type=int, help="Optional Kraken2 --minimum-hit-groups value")
    parser.add_argument("--input-compression", choices=["auto", "none", "gzip", "bzip2"], default="auto", help="Input compression mode for Kraken2 FASTQ files")
    add_bool_pair(parser, "--report-minimizer-data", default=True, help_text="Write Kraken2 minimizer data columns into the report; recommended for Celltaminate scoring")
    add_bool_pair(parser, "--use-names", default=False, help_text="Pass --use-names to Kraken2 output; not needed by Celltaminate")
    parser.add_argument("--extra-kraken2-args", default="", help="Additional raw arguments appended to every Kraken2 command")
    parser.add_argument("--skip-existing-reports", action="store_true", help="Reuse non-empty report files already present in outdir/kraken2_reports")
    parser.add_argument("--skip-kraken2", action="store_true", help="Skip Kraken2 and use a report column or existing reports")
    parser.add_argument("--skip-celltaminate", action="store_true", help="Generate Kraken2 reports and metadata only")

    parser.add_argument("--rscript", default=os.environ.get("RSCRIPT", "Rscript"), help="Rscript command or path")
    parser.add_argument("--celltaminate-r", help="Path to Celltaminate.R. Defaults to scripts/Celltaminate.R under repo root")
    parser.add_argument("--tax-level", choices=["G", "S"], default="S", help="Taxonomic level for Celltaminate")
    parser.add_argument("--host-species", default="human", help="Comma-separated host names: human,mouse,drosophila")
    parser.add_argument("--host-taxids-manual", default="", help="Additional host taxids, comma-separated")
    add_bool_pair(parser, "--include-protozoa", default=True, help_text="Include non-host/non-plant eukaryotes as microbial candidates")
    add_bool_pair(parser, "--collapse-species", default=True, help_text="Collapse ambiguous species within a genus")
    parser.add_argument("--collapse-min-genus-reads", type=int, default=30)
    parser.add_argument("--collapse-top-frac", type=float, default=0.85)
    parser.add_argument("--fp-falsepos-cutoff", type=float, default=75)
    parser.add_argument("--fp-true-cutoff", type=float, default=5)
    parser.add_argument("--fp-aggressiveness", type=float, default=1.15)
    parser.add_argument("--quick-remove", default="false,unc")
    add_bool_pair(parser, "--keep-descendants", default=False, help_text="Keep descendants of removed taxa in cleaned reports")
    add_bool_pair(parser, "--show-fp-breakdown", default=False, help_text="Add false-positive component columns to output tables")
    add_bool_pair(parser, "--save-cleaned-reanalysis", default=True, help_text="Save Celltaminate cleaned reanalysis outputs")
    parser.add_argument("--ref-bg-tsv", help="Full refined_cell.lines.tsv background reference")
    parser.add_argument("--kitome-blacklist-tsv", help="Kitome/background blacklist TSV")
    parser.add_argument("--clinical-panel-tsv", help="Clinically important pathogen TSV")
    parser.add_argument("--protected-taxa-tsv", help="Optional TSV with sample and taxon columns")

    parser.add_argument("--dry-run", action="store_true", help="Print commands and write input files, but do not execute Kraken2/R")
    return parser


def main() -> int:
    parser = make_parser()
    args = parser.parse_args()
    return cmd_full_workflow(args)


if __name__ == "__main__":
    raise SystemExit(main())
