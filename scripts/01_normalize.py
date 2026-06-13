"""Normalize TPMs and average per age bracket, for every tissue.

Reads:
    data/parquet/{version}/tpm_matrix.parquet
    data/parquet/{version}/metadata.parquet

Writes:
    output/{version}/normalized/{tissue}.parquet  (one per tissue)

Usage
-----
    python scripts/01_normalize.py --version v10
    python scripts/01_normalize.py --version v10 --tissue "Liver"      # one tissue
    python scripts/01_normalize.py --version v10 --epsilon 0.01        # custom filter
    python scripts/01_normalize.py --version v10 --resume              # skip done

Memory notes
------------
Loading the full TPM matrix takes ~2 GB of RAM (float32, 60k x 19k).
For a single tissue, only the relevant sample columns are read from
Parquet (selective column read), so peak RAM stays well below 1 GB.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from stamp.config import EPSILON, GtexVersion, paths_for
from stamp.io import (
    _safe_filename,
    load_metadata,
    load_tpm_matrix,
    save_normalized_tissue,
)
from stamp.normalize import normalize_tissue
from stamp.tissues import tissues_with_complete_age_bins


def get_tissues_to_process(
    metadata: pd.DataFrame,
    version: GtexVersion,
    only: str | None,
    resume: bool,
    complete: bool = False,
) -> list[str]:
    """Return the list of tissues that need normalization in this run."""
    all_tissues = sorted(metadata["tissue"].dropna().unique().tolist())

    # Complete age-bin coverage: keep only tissues with a sample in all 6 brackets.
    if complete:
        complete_set = set(tissues_with_complete_age_bins(metadata))
        excluded = [t for t in all_tissues if t not in complete_set]
        all_tissues = [t for t in all_tissues if t in complete_set]
        print(
            f"  Complete age-bins mode: {len(all_tissues)} complete tissues "
            f"({len(excluded)} excluded: {sorted(excluded)})"
        )

    if only is not None:
        if only not in all_tissues:
            raise ValueError(
                f"Tissue '{only}' not in metadata (or excluded by --complete-age-bins). "
                f"Available tissues: {all_tissues}"
            )
        tissues = [only]
    else:
        tissues = all_tissues

    if resume:
        done_dir = paths_for(version, complete)["normalized"]
        if done_dir.exists():
            done = {p.stem for p in done_dir.glob("*.parquet")}
            already = [t for t in tissues if _safe_filename(t) in done]
            if already:
                print(f"  Resume: skipping {len(already)} already-normalized tissues")
            tissues = [t for t in tissues if _safe_filename(t) not in done]

    return tissues


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", choices=["v8", "v10"], required=True)
    parser.add_argument("--tissue", default=None,
                        help="Process only this tissue (exact match).")
    parser.add_argument("--epsilon", type=float, default=EPSILON,
                        help=f"Filter genes with (max-min) < epsilon. Default: {EPSILON}")
    parser.add_argument("--resume", action="store_true",
                        help="Skip tissues already present in the output dir.")
    parser.add_argument("--complete-age-bins", action="store_true",
                        help="Process only tissues with a sample in all 6 GTEx "
                             "age brackets; write to output/{version}_complete/.")
    args = parser.parse_args()

    mode = " [complete age-bins]" if args.complete_age_bins else ""
    print(f"=== Normalization (GTEx {args.version}){mode} ===")
    print(f"  Loading metadata...")
    metadata = load_metadata(args.version)
    print(f"  Metadata: {len(metadata)} samples")

    tissues = get_tissues_to_process(
        metadata, args.version, args.tissue, args.resume,
        complete=args.complete_age_bins,
    )
    if not tissues:
        print("  Nothing to do.")
        return

    print(f"  Tissues to process: {len(tissues)}")
    print(f"  Loading TPM matrix...")
    t_load = time.time()
    tpm = load_tpm_matrix(args.version)
    print(f"  TPM matrix shape: {tpm.shape} (loaded in {time.time()-t_load:.1f}s)")

    failures: list[tuple[str, str]] = []
    for tissue in tqdm(tissues, desc="Tissues", unit="tissue"):
        try:
            normalized = normalize_tissue(
                tpm, metadata, tissue, epsilon=args.epsilon
            )
            save_normalized_tissue(
                args.version, tissue, normalized, complete=args.complete_age_bins
            )
        except Exception as e:
            failures.append((tissue, str(e)))
            tqdm.write(f"  [FAIL] {tissue}: {e}")

    print(f"\nDone.")
    print(f"  Successful: {len(tissues) - len(failures)} / {len(tissues)}")
    if failures:
        print(f"  Failures:")
        for t, msg in failures:
            print(f"    - {t}: {msg}")
        sys.exit(1)


if __name__ == "__main__":
    main()