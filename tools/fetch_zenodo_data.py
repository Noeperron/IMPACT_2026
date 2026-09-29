#!/usr/bin/env python3
"""Download the deposited data from Zenodo into the layout the figure scripts expect.

Purpose:      The Zenodo record stores every file at the top level with its directory encoded
              in the filename (for example elisa_elisa_mmr.csv). The figure scripts read a
              nested tree under data/ (data/elisa/elisa_mmr.csv). This script downloads the
              record and unflattens it, so no file has to be renamed or moved by hand.
Inputs:       The Zenodo record, resolved from the concept DOI so it always fetches the latest
              version. Nothing local is required.
Outputs:      A populated data/ directory (40 files, about 10.4 GB in full, about 190 MB with
              --tabular-only). Existing files whose checksum already matches are left alone,
              so the script is safe to re-run and to resume after an interrupted download.
Dependencies: Python + requests (standard library otherwise).

Usage:
    python3 tools/fetch_zenodo_data.py                  # everything
    python3 tools/fetch_zenodo_data.py --tabular-only   # skip the single-cell objects
    python3 tools/fetch_zenodo_data.py --dry-run        # list what would be fetched
    python3 tools/fetch_zenodo_data.py --token <TOKEN>  # while the record is under embargo

While the record is embargoed the files are restricted and a share token is required; pass it
with --token or set ZENODO_TOKEN. Once the embargo lifts no token is needed.
"""
import argparse
import hashlib
import os
import sys
from pathlib import Path

try:
    import requests
except ImportError:
    sys.exit("This script needs 'requests'. Install it with: pip install requests")

CONCEPT_DOI = "10.5281/zenodo.18989222"   # always resolves to the latest version
CONCEPT_ID = CONCEPT_DOI.rsplit(".", 1)[1]
API = "https://zenodo.org/api"

# The record stores files flat, with the directory encoded as a filename prefix. Longest
# prefix first, so external_zavidij_bm_ is matched before external_.
PREFIX_MAP = [
    ("external_zavidij_bm_", "external/zavidij_bm"),
    ("elisa_", "elisa"),
    ("external_", "external"),
    ("metadata_", "metadata"),
    ("olink_", "olink"),
    ("subclusters_", "subclusters"),
    ("tcr_", "tcr"),
]

# Only these scripts need the single-cell objects; everything else runs from the tabular
# inputs, which is why --tabular-only exists.
LARGE_SUFFIX = ".h5ad"


def local_path(deposit_name: str) -> str:
    """Map a flat deposit filename back to its path under data/."""
    for prefix, directory in PREFIX_MAP:
        if deposit_name.startswith(prefix):
            return f"{directory}/{deposit_name[len(prefix):]}"
    return deposit_name


def md5(path: Path, chunk: int = 8 << 20) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while True:
            block = fh.read(chunk)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def human(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} GB"


def resolve_record(session, record_id, params):
    """Return the record JSON, following the concept id to the latest version if needed."""
    r = session.get(f"{API}/records/{record_id}", params=params, timeout=60)
    if r.status_code == 404:
        sys.exit(f"Record {record_id} not found.")
    if r.status_code in (401, 403):
        sys.exit("Access denied. The record is probably still under embargo; pass --token "
                 "or set ZENODO_TOKEN.")
    r.raise_for_status()
    return r.json()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dest", default=None,
                    help="destination directory (default: the repo's data/ directory)")
    ap.add_argument("--record", default=None,
                    help="fetch a specific record id instead of the latest version")
    ap.add_argument("--token", default=os.environ.get("ZENODO_TOKEN"),
                    help="Zenodo share token, needed while the record is embargoed")
    ap.add_argument("--tabular-only", action="store_true",
                    help="skip the single-cell .h5ad objects (about 10.2 GB of the 10.4 GB)")
    ap.add_argument("--dry-run", action="store_true", help="list files without downloading")
    args = ap.parse_args()

    dest = Path(args.dest) if args.dest else Path(__file__).resolve().parent.parent / "data"
    params = {"token": args.token} if args.token else {}
    session = requests.Session()

    rec = resolve_record(session, args.record or CONCEPT_ID, params)
    print(f"Record {rec['id']}  DOI {rec['doi']}  ({rec['metadata'].get('publication_date','')})")

    files = rec.get("files", [])
    if not files:
        sys.exit("The record lists no files. It is probably embargoed; pass --token or set "
                 "ZENODO_TOKEN.")
    if args.tabular_only:
        files = [f for f in files if not f["key"].endswith(LARGE_SUFFIX)]

    total = sum(f["size"] for f in files)
    print(f"{len(files)} files, {human(total)} -> {dest}\n")

    done = skipped = 0
    for f in sorted(files, key=lambda x: x["key"]):
        rel = local_path(f["key"])
        out = dest / rel
        want = f["checksum"].split(":", 1)[-1]

        # Size is checked first: a mismatch means the file is absent or truncated, and there
        # is no point hashing 9.6 GB to learn that. A dry run stops at the size check, so it
        # stays fast on a fully populated data directory.
        same_size = out.exists() and out.stat().st_size == f["size"]
        if same_size and (args.dry_run or md5(out) == want):
            print(f"  ok      {rel}")
            skipped += 1
            continue
        if args.dry_run:
            print(f"  fetch   {rel}  ({human(f['size'])})")
            continue

        out.parent.mkdir(parents=True, exist_ok=True)
        url = f["links"]["self"]
        print(f"  get     {rel}  ({human(f['size'])}) ... ", end="", flush=True)
        tmp = out.with_suffix(out.suffix + ".part")
        with session.get(url, params=params, stream=True, timeout=300) as resp:
            resp.raise_for_status()
            with open(tmp, "wb") as fh:
                for block in resp.iter_content(chunk_size=8 << 20):
                    fh.write(block)
        got = md5(tmp)
        if got != want:
            tmp.unlink(missing_ok=True)
            sys.exit(f"\nchecksum mismatch for {rel}: expected {want}, got {got}")
        tmp.rename(out)
        print("ok")
        done += 1

    if args.dry_run:
        print(f"\nDry run. {len(files) - skipped} file(s) would be downloaded.")
    else:
        print(f"\nDone. {done} downloaded, {skipped} already present and verified.")
        if args.tabular_only:
            print("Single-cell objects were skipped; re-run without --tabular-only for "
                  "Figure 3B, 3C, 5C, Supplementary Figures 2, 3B and 6.")


if __name__ == "__main__":
    main()
