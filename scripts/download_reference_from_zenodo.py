#!/usr/bin/env python3
"""Download the full Celltaminate reference background table from Zenodo."""

from __future__ import annotations

import argparse
import sys
import urllib.error
import urllib.request
from pathlib import Path

DEFAULT_URL = "https://zenodo.org/records/20560460/files/refined_cell.lines.tsv?download=1"
DEFAULT_OUT = Path("data/reference/refined_cell.lines.tsv")


def download(url: str, out: Path) -> None:
    out = out.expanduser().resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(out.suffix + ".tmp")

    request = urllib.request.Request(url, headers={"User-Agent": "Celltaminate reference downloader"})

    try:
        with urllib.request.urlopen(request) as response, tmp.open("wb") as handle:
            total_raw = response.headers.get("Content-Length")
            total = int(total_raw) if total_raw and total_raw.isdigit() else None
            copied = 0
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                handle.write(chunk)
                copied += len(chunk)
                if total:
                    pct = 100 * copied / total
                    print(f"Downloaded {copied / (1024 ** 2):.1f} / {total / (1024 ** 2):.1f} MB ({pct:.1f}%)", end="\r")
                else:
                    print(f"Downloaded {copied / (1024 ** 2):.1f} MB", end="\r")
        print()
        tmp.replace(out)
    except (urllib.error.URLError, urllib.error.HTTPError, OSError) as exc:
        if tmp.exists():
            tmp.unlink()
        raise SystemExit(f"ERROR: failed to download reference table: {exc}") from exc


def main() -> int:
    parser = argparse.ArgumentParser(description="Download refined_cell.lines.tsv from Zenodo.")
    parser.add_argument("--url", default=DEFAULT_URL, help="Download URL")
    parser.add_argument("--out", default=str(DEFAULT_OUT), help="Output TSV path")
    args = parser.parse_args()

    download(args.url, Path(args.out))
    print(f"Reference table written to: {Path(args.out).expanduser().resolve()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
