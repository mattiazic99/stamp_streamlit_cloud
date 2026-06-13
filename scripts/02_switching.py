"""Identify switching genes per tissue and write 5-line sets text files.

For each tissue with a normalized matrix in output/{version}/normalized,
applies a binarisation threshold and Definition 1 of switching genes
(see stamp.switching) to produce a 5-line plain-text file:

    line 1 -> genes switching at [30-39]
    line 2 -> genes switching at [40-49]
    line 3 -> genes switching at [50-59]
    line 4 -> genes switching at [60-69]
    line 5 -> genes switching at [70-79]

Genes are space-separated; empty brackets produce empty lines.

Usage
-----
    python scripts/02_switching.py --version v10
    python scripts/02_switching.py --version v10 --threshold 0.3
    python scripts/02_switching.py --version v10 --tissue "Liver"

Tissues with missing brackets
-----------------------------
A tissue with fewer than 6 age brackets in its normalized matrix
(e.g. Cervix_Ectocervix has no [70-79]) still produces a 5-line file:
brackets that have no samples produce empty lines, since the corresponding
bracket index never gets reached by any switching event in that tissue.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from stamp.config import DEFAULT_THRESHOLD, GtexVersion, paths_for
from stamp.io import _safe_filename, load_normalized_tissue, save_sets_txt
from stamp.switching import identify_switching_genes


def get_tissues_to_process(
    version: GtexVersion,
    only: list[str] | None,
    complete: bool = False,
) -> list[str]:
    """Return tissues that have a normalized Parquet, in alphabetical order."""
    norm_dir = paths_for(version, complete)["normalized"]
    if not norm_dir.exists():
        raise FileNotFoundError(
            f"No normalized data found at {norm_dir}. "
            f"Run scripts/01_normalize.py --version {version}"
            f"{' --complete-age-bins' if complete else ''} first."
        )

    available = sorted(p.stem for p in norm_dir.glob("*.parquet"))
    if not available:
        raise FileNotFoundError(f"No .parquet files in {norm_dir}.")

    if only is not None:
        targets = [_safe_filename(t) for t in only]
        missing = [t for t in targets if t not in available]
        if missing:
            raise ValueError(
                f"Tissues not found: {missing}. "
                f"Available: {available}"
            )
        return targets

    return available


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", choices=["v8", "v10"], required=True)
    parser.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD,
                        help=f"Binarisation threshold tau in [0, 1]. "
                        f"Default: {DEFAULT_THRESHOLD}")
    parser.add_argument("--tissues", nargs="+", default=None,
                        help="Process only these tissues (display names).")
    parser.add_argument("--complete-age-bins", action="store_true",
                        help="Use the complete-age-bins outputs "
                             "(output/{version}_complete/).")
    args = parser.parse_args()

    if not 0.0 <= args.threshold <= 1.0:
        print(f"Error: threshold must be in [0, 1], got {args.threshold}")
        sys.exit(1)

    mode = " [complete age-bins]" if args.complete_age_bins else ""
    print(f"=== Switching detection (GTEx {args.version}, tau={args.threshold}){mode} ===")

    tissues = get_tissues_to_process(args.version, args.tissues, complete=args.complete_age_bins)
    print(f"  Tissues to process: {len(tissues)}")

    # Stats accumulator: {tissue: total_switching_genes}
    stats: dict[str, int] = {}
    failures: list[tuple[str, str]] = []

    t0 = time.time()
    for tissue_safe in tqdm(tissues, desc="Tissues", unit="tissue"):
        try:
            df = load_normalized_tissue(args.version, tissue_safe, complete=args.complete_age_bins)
            sets = identify_switching_genes(df, threshold=args.threshold)
            save_sets_txt(args.version, tissue_safe, sets, complete=args.complete_age_bins)
            stats[tissue_safe] = sum(len(v) for v in sets.values())
        except Exception as e:
            failures.append((tissue_safe, str(e)))
            tqdm.write(f"  [FAIL] {tissue_safe}: {e}")

    elapsed = time.time() - t0
    print(f"\nDone in {elapsed:.1f}s.")
    print(f"  Successful: {len(tissues) - len(failures)} / {len(tissues)}")

    if stats:
        print(f"\nSwitching gene counts (tau={args.threshold}):")
        for tissue, n in sorted(stats.items(), key=lambda x: -x[1]):
            print(f"  {tissue:<55s} {n:>6d}")

    if failures:
        print(f"\nFailures:")
        for t, msg in failures:
            print(f"  - {t}: {msg}")
        sys.exit(1)


if __name__ == "__main__":
    main()