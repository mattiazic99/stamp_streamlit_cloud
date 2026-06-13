"""Compute Jaccard similarity matrices between tissues.

Reads:
    output/{version}/sets/*_sets.txt

Writes:
    output/{version}/jaccard/jaccard_age.csv
    output/{version}/jaccard/jaccard_life.csv

Usage
-----
    python scripts/03_jaccard.py --version v10
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import pandas as pd

from stamp.config import GtexVersion, paths_for
from stamp.io import load_sets_txt
from stamp.jaccard import tissue_similarity_age, tissue_similarity_life


def load_all_sets(
    version: GtexVersion, complete: bool = False
) -> dict[str, dict[str, list[str]]]:
    """Load every _sets.txt file in output/{version}[_complete]/sets/."""
    sets_dir = paths_for(version, complete)["sets"]
    if not sets_dir.exists():
        raise FileNotFoundError(
            f"No sets directory at {sets_dir}. "
            f"Run scripts/02_switching.py --version {version}"
            f"{' --complete-age-bins' if complete else ''} first."
        )

    files = sorted(sets_dir.glob("*_sets.txt"))
    if not files:
        raise FileNotFoundError(f"No *_sets.txt files in {sets_dir}.")

    out: dict[str, dict[str, list[str]]] = {}
    for f in files:
        tissue = f.stem.replace("_sets", "")
        out[tissue] = load_sets_txt(version, tissue, complete=complete)
    return out


def save_matrix(mat: pd.DataFrame, path: Path) -> None:
    """Save a similarity matrix as CSV with 6-decimal precision."""
    path.parent.mkdir(parents=True, exist_ok=True)
    mat.round(6).to_csv(path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", choices=["v8", "v10"], required=True)
    parser.add_argument("--complete-age-bins", action="store_true",
                        help="Use the complete-age-bins outputs "
                             "(output/{version}_complete/).")
    args = parser.parse_args()

    mode = " [complete age-bins]" if args.complete_age_bins else ""
    print(f"=== Jaccard similarity (GTEx {args.version}){mode} ===")
    t0 = time.time()

    print("  Loading sets files...")
    sets_by_tissue = load_all_sets(args.version, complete=args.complete_age_bins)
    print(f"  Loaded {len(sets_by_tissue)} tissues.")

    print("  Computing J_age matrix...")
    mat_age = tissue_similarity_age(sets_by_tissue)
    out_age = paths_for(args.version, args.complete_age_bins)["jaccard"] / "jaccard_age.csv"
    save_matrix(mat_age, out_age)
    print(f"    saved to {out_age}")

    print("  Computing J_life matrix...")
    mat_life = tissue_similarity_life(sets_by_tissue)
    out_life = paths_for(args.version, args.complete_age_bins)["jaccard"] / "jaccard_life.csv"
    save_matrix(mat_life, out_life)
    print(f"    saved to {out_life}")

    elapsed = time.time() - t0
    print(f"\nDone in {elapsed:.1f}s.")

    # Print a few highlights
    print("\nTop 10 tissue pairs by J_life similarity (excluding diagonal):")
    flat = mat_life.where(~_diag_mask(mat_life)).stack().sort_values(ascending=False)
    seen = set()
    shown = 0
    for (a, b), v in flat.items():
        if (b, a) in seen:
            continue
        seen.add((a, b))
        print(f"  {a:<40s}  <-> {b:<40s}  {v:.4f}")
        shown += 1
        if shown >= 10:
            break


def _diag_mask(mat: pd.DataFrame):
    import numpy as np
    return pd.DataFrame(
        np.eye(len(mat), dtype=bool),
        index=mat.index, columns=mat.columns
    )


if __name__ == "__main__":
    main()