"""Check GTEx metadata and select the 50 tissues with complete age-bin coverage."""

from __future__ import annotations

import sys

import pandas as pd

from stamp.config import AGE_BRACKETS
from stamp.io import load_metadata
from stamp.tissues import incomplete_tissues

from .common import (
    base_parser,
    enforce_final_scope,
    parse_tissue_arguments,
    select_tissues,
    validation_root,
    write_json,
)


def run(args) -> int:
    metadata = load_metadata(args.version)
    requested = parse_tissue_arguments(args.tissue)
    all_complete = select_tissues(metadata)
    selected = select_tissues(metadata, requested=requested)
    enforce_final_scope(selected, all_complete, args.run_label, "Audit")
    out_dir = validation_root(
        args.version, args.output_root, args.run_label
    ) / "01_audit"
    out_dir.mkdir(parents=True, exist_ok=True)

    valid = metadata.loc[metadata["age_bracket"].isin(AGE_BRACKETS)].copy()
    counts = (
        valid.pivot_table(
            index="tissue",
            columns="age_bracket",
            values="sample_id",
            aggfunc="count",
            fill_value=0,
        )
        .reindex(columns=AGE_BRACKETS, fill_value=0)
        .astype(int)
    )
    counts["n_total"] = counts.sum(axis=1)
    counts["n_min"] = counts[list(AGE_BRACKETS)].min(axis=1)
    counts["complete_six_bins"] = (counts[list(AGE_BRACKETS)] > 0).all(axis=1)
    counts["selected_primary"] = counts.index.isin(selected)
    counts = counts.reset_index()
    counts.to_csv(out_dir / "tissue_age_bin_counts.csv", index=False)

    duplicate_table = (
        valid.groupby(["tissue", "subject_id"], observed=True)
        .agg(n_samples=("sample_id", "size"))
        .reset_index()
    )
    duplicate_table = duplicate_table.loc[duplicate_table["n_samples"] > 1]
    duplicate_table.to_csv(
        out_dir / "duplicate_subject_tissue_samples.csv", index=False
    )

    selected_counts = counts.loc[counts["selected_primary"]].copy()
    selected_counts.to_csv(out_dir / "selected_tissues.csv", index=False)

    summary = {
        "version": args.version,
        "n_metadata_rows": int(metadata.shape[0]),
        "n_subjects": int(metadata["subject_id"].nunique()),
        "n_all_tissues": int(metadata["tissue"].nunique()),
        "n_complete_tissues": len(all_complete),
        "n_selected_tissues": len(selected),
        "selected_tissues": selected,
        "incomplete_tissues": incomplete_tissues(metadata),
        "n_duplicate_subject_tissue_groups": int(duplicate_table.shape[0]),
        "age_brackets": list(AGE_BRACKETS),
    }
    write_json(out_dir / "audit_summary.json", summary)

    print(f"GTEx {args.version}: {summary['n_all_tissues']} tissues found")
    print(f"Complete-tissue universe: {len(all_complete)} tissues")
    print(f"Tissues selected for this run: {len(selected)}")
    print(f"Incomplete tissues: {summary['incomplete_tissues']}")
    print(
        "Duplicate donor+tissue groups: "
        f"{summary['n_duplicate_subject_tissue_groups']}"
    )
    print(f"Audit written to: {out_dir}")

    if duplicate_table.shape[0] > 0:
        print(
            "ERROR: multiple samples from the same donor and tissue require "
            "aggregation or blocked resampling before validation.",
            file=sys.stderr,
        )
        return 2
    return 0


def build_parser():
    return base_parser(
        "Audit sample counts, donor uniqueness and complete age-bin coverage."
    )


def main() -> None:
    args = build_parser().parse_args()
    raise SystemExit(run(args))


if __name__ == "__main__":
    main()
