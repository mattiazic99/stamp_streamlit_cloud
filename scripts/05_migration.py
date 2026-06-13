"""Compute cross-tissue gene migration paths.

Reads:
    output/{version}/sets/*_sets.txt

Writes:
    output/{version}/migration/migrations_raw.csv
    output/{version}/migration/migration_paths.csv
    output/{version}/migration/migration_summary.json

Usage
-----
    python scripts/05_migration.py --version v10
    python scripts/05_migration.py --version v8
    python scripts/05_migration.py --version v10 --top-n 200
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

from stamp.config import GtexVersion, paths_for
from stamp.io import load_sets_txt
from stamp.migration import find_migrations, migration_paths, migration_summary


def load_all_sets(version: GtexVersion) -> dict[str, dict[str, list[str]]]:
    sets_dir = paths_for(version)["sets"]
    if not sets_dir.exists():
        raise FileNotFoundError(f"No sets directory at {sets_dir}.")
    files = sorted(sets_dir.glob("*_sets.txt"))
    if not files:
        raise FileNotFoundError(f"No *_sets.txt files in {sets_dir}.")
    out = {}
    for f in files:
        tissue = f.stem.replace("_sets", "")
        out[tissue] = load_sets_txt(version, tissue)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", choices=["v8", "v10"], required=True)
    parser.add_argument("--top-n", type=int, default=100,
                        help="Top N migration paths to save. Default: 100.")
    args = parser.parse_args()

    print(f"=== Migration analysis (GTEx {args.version}) ===")
    t0 = time.time()

    print("  Loading sets files...")
    sets_by_tissue = load_all_sets(args.version)
    print(f"  Loaded {len(sets_by_tissue)} tissues.")

    print("  Finding migration events...")
    migrations = find_migrations(sets_by_tissue)
    elapsed = time.time() - t0
    print(f"  Done in {elapsed:.1f}s.")

    summary = migration_summary(migrations)
    print(f"\n  Migration events:   {summary['n_migration_events']:,}")
    print(f"  Migrating genes:    {summary['n_migrating_genes']:,}")
    print(f"  Source tissues:     {summary.get('n_source_tissues', 0)}")
    print(f"  Target tissues:     {summary.get('n_target_tissues', 0)}")
    if summary["n_migration_events"] > 0:
        print(f"  Top source tissue:  {summary.get('top_source_tissue', 'N/A')}")
        print(f"  Top target tissue:  {summary.get('top_target_tissue', 'N/A')}")

    out_dir = paths_for(args.version)["migration"]
    out_dir.mkdir(parents=True, exist_ok=True)

    # Save raw migrations (can be large)
    raw_path = out_dir / "migrations_raw.csv"
    migrations.to_csv(raw_path, index=False)
    size_mb = raw_path.stat().st_size / 1e6
    print(f"\n  Saved raw: {raw_path} ({size_mb:.1f} MB)")

    # Save aggregated paths
    paths_df = migration_paths(migrations, top_n=args.top_n)
    paths_path = out_dir / "migration_paths.csv"
    paths_df.to_csv(paths_path, index=False)
    print(f"  Saved paths: {paths_path}")

    # Save summary as JSON
    summ_path = out_dir / "migration_summary.json"
    summ_path.write_text(json.dumps(summary, indent=2))
    print(f"  Saved summary: {summ_path}")

    # Print top 10 paths
    if not paths_df.empty:
        print(f"\nTop 10 migration paths (tau=0.5):")
        print(f"  {'Source tissue':<40s} {'Bracket':<8s} -> "
              f"{'Target tissue':<40s} {'Bracket':<8s} {'Genes':>6s}")
        print("  " + "-" * 108)
        for _, row in paths_df.head(10).iterrows():
            print(f"  {row['source_tissue']:<40s} {row['source_bracket']:<8s} -> "
                  f"  {row['target_tissue']:<40s} {row['target_bracket']:<8s} "
                  f"{row['n_genes']:>6d}")


if __name__ == "__main__":
    main()