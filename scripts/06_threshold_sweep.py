"""Sweep tau across [0.05, 0.95] and count switching genes per tissue.

Reads:
    output/{version}/normalized/*.parquet

Writes:
    output/{version}/threshold_sweep/sweep_results.csv
    output/{version}/threshold_sweep/stability_scores.csv

Usage
-----
    python scripts/06_threshold_sweep.py --version v10
    python scripts/06_threshold_sweep.py --version v8
    python scripts/06_threshold_sweep.py --version v10 --taus 0.3 0.4 0.5 0.6 0.7
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from stamp.config import GtexVersion, paths_for
from stamp.io import load_normalized_tissue
from stamp.threshold import DEFAULT_TAUS, stability_score, sweep_all_tissues


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", choices=["v8", "v10"], required=True)
    parser.add_argument(
        "--taus", nargs="+", type=float, default=None,
        help="Custom tau values to test (default: 0.05 to 0.95 step 0.05)."
    )
    args = parser.parse_args()

    taus = args.taus or DEFAULT_TAUS
    print(f"=== Threshold sweep (GTEx {args.version}) ===")
    print(f"  Tau grid: {len(taus)} values from {min(taus):.2f} to {max(taus):.2f}")

    norm_dir = paths_for(args.version)["normalized"]
    if not norm_dir.exists():
        print(f"[ERROR] {norm_dir} not found. Run 01_normalize.py first.")
        sys.exit(1)

    tissue_files = sorted(norm_dir.glob("*.parquet"))
    if not tissue_files:
        print(f"[ERROR] No .parquet files in {norm_dir}.")
        sys.exit(1)

    print(f"  Tissues: {len(tissue_files)}")
    print(f"  Total combinations: {len(tissue_files) * len(taus)}")

    # Load all normalized matrices
    print("  Loading normalized matrices...")
    tissues_normalized: dict[str, pd.DataFrame] = {}
    for f in tqdm(tissue_files, desc="Loading", unit="tissue"):
        tissue = f.stem
        df = load_normalized_tissue(args.version, tissue)
        tissues_normalized[tissue] = df

    # Run sweep
    print("  Running sweep...")
    t0 = time.time()
    sweep_df = sweep_all_tissues(tissues_normalized, taus=taus)
    elapsed = time.time() - t0
    print(f"  Done in {elapsed:.1f}s.")

    # Save sweep results
    out_dir = paths_for(args.version)["threshold_sweep"]
    out_dir.mkdir(parents=True, exist_ok=True)

    sweep_path = out_dir / "sweep_results.csv"
    sweep_df.to_csv(sweep_path, index=False)
    print(f"  Saved: {sweep_path}")

    # Compute and save stability scores
    stability = stability_score(sweep_df)
    stab_path = out_dir / "stability_scores.csv"
    stability.to_csv(stab_path, header=True)
    print(f"  Saved: {stab_path}")

    # Print summary
    print(f"\nSweep summary at tau=0.5:")
    at_05 = sweep_df[sweep_df["tau"] == 0.5].set_index("tissue")["n_switching"]
    print(f"  Max: {at_05.idxmax()} ({at_05.max()})")
    print(f"  Min: {at_05.idxmin()} ({at_05.min()})")
    print(f"  Mean across tissues: {at_05.mean():.0f}")

    print(f"\nTop 10 most STABLE tissues (low CV of n_switching, tau 0.3-0.7):")
    for tissue, cv in stability.head(10).items():
        print(f"  {tissue:<50s} CV={cv:.3f}")

    print(f"\nTop 10 most UNSTABLE tissues (high CV):")
    for tissue, cv in stability.tail(10).items():
        print(f"  {tissue:<50s} CV={cv:.3f}")


if __name__ == "__main__":
    main()