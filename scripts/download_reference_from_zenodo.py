#!/usr/bin/env python3
"""Download the full Celltaminate reference background table from Zenodo."""

from __future__ import annotations

import argparse
import json
import sys
import urllib.error
import urllib.request
from pathlib import Path
from typing import List

REFERENCE_DOI = "10.5281/zenodo.20560460"
DEFAULT_URL = "https://zenodo.org/records/20560460/files/refined_cell.lines.tsv?download=1"
DEFAULT_OUT = Path("data/reference/refined_cell.lines.tsv")


def candidate_urls(primary_url: str = DEFAULT_URL, doi: str = REFERENCE_DOI) -> List[str]:
    urls: List[str] = []
    if primary_url:
        urls.append(primary_url)

    record_id = doi.rstrip("/").split(".")[-1]
    if record_id.isdigit():
        api_url = f"https://zenodo.org/api/records/{record_id}"
        try:
            request = urllib.request.Request(api_url, headers={"User-Agent": "Celltaminate reference downloader"})
            with urllib.request.urlopen(request, timeout=30) as response:
                payload = json.loads(response.read().decode("utf-8"))
            files = payload.get("files", []) or []
            ranked = []
            for f in files:
                key = str(f.get("key") or f.get("filename") or "")
                links = f.get("links", {}) or {}
                link = links.get("self") or links.get("download")
                size = int(f.get("size") or 0)
                if not link:
                    continue
                score = 0
                if key == "refined_cell.lines.tsv":
                    score += 100
                if key.endswith(".tsv"):
                    score += 10
                ranked.append((score, size, link))
            for _, _, link in sorted(ranked, reverse=True):
                if link not in urls:
                    urls.append(link)
        except Exception:
            pass
    return urls


def download_one(url: str, out: Path) -> None:
    tmp = out.with_suffix(out.suffix + ".tmp")
    request = urllib.request.Request(url, headers={"User-Agent": "Celltaminate reference downloader"})
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


def download(url: str, out: Path, force: bool = False) -> None:
    out = out.expanduser().resolve()
    out.parent.mkdir(parents=True, exist_ok=True)

    if out.exists() and not force:
        raise SystemExit(f"ERROR: output already exists: {out}\nUse --force to overwrite it.")

    last_error = None
    for candidate in candidate_urls(url):
        tmp = out.with_suffix(out.suffix + ".tmp")
        try:
            print(f"Trying: {candidate}")
            download_one(candidate, out)
            return
        except (urllib.error.URLError, urllib.error.HTTPError, OSError) as exc:
            last_error = exc
            if tmp.exists():
                tmp.unlink()

    raise SystemExit(
        "ERROR: failed to download reference table. "
        f"Last error: {last_error}\n"
        "Confirm that DOI 10.5281/zenodo.20560460 is public and that the file name is refined_cell.lines.tsv, "
        "or pass an explicit --url from the Zenodo file download link."
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Download refined_cell.lines.tsv from Zenodo.")
    parser.add_argument("--url", default=DEFAULT_URL, help="Download URL. If this fails, the script also tries the Zenodo record API.")
    parser.add_argument("--out", default=str(DEFAULT_OUT), help="Output TSV path")
    parser.add_argument("--force", action="store_true", help="Overwrite an existing output file")
    args = parser.parse_args()

    print(f"Downloading reference background table from DOI {REFERENCE_DOI}")
    print(f"Output: {Path(args.out).expanduser().resolve()}")
    download(args.url, Path(args.out), force=args.force)
    print(f"Reference table written to: {Path(args.out).expanduser().resolve()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
