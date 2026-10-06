"""Check calculated switching sets against the published STAMP atlas."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from stamp.config import DEFAULT_THRESHOLD, EPSILON, SWITCHING_BRACKETS
from stamp.io import load_metadata, load_sets_txt

from .common import (
    age_means_from_codes,
    base_parser,
    load_tissue_data,
    parse_tissue_arguments,
    select_tissues,
    switching_calls_frame,
    switching_state,
    validation_root,
)



def _check_tissue(version: str, tissue: str, tau: float, epsilon: float) -> dict:
    data = load_tissue_data(version, tissue, epsilon=epsilon)
    state = switching_state(
        age_means_from_codes(data.normalized_values, data.age_codes),
        tau=tau,
    )
    calls = switching_calls_frame(data.gene_ids, state)
    canonical = load_sets_txt(version, tissue, complete=True)

    observed_by_bracket = {
        bracket: set(
            calls.loc[calls["switch_bracket"] == bracket, "gene_id"].astype(str)
        )
        for bracket in SWITCHING_BRACKETS
    }
    canonical_by_bracket = {
        bracket: set(canonical[bracket]) for bracket in SWITCHING_BRACKETS
    }
    observed_count = int(sum(len(genes) for genes in observed_by_bracket.values()))
    canonical_count = int(sum(len(genes) for genes in canonical_by_bracket.values()))
    exact_match = all(
        observed_by_bracket[bracket] == canonical_by_bracket[bracket]
        for bracket in SWITCHING_BRACKETS
    )
    return {
        "version": version,
        "tissue": tissue,
        "observed_switching": observed_count,
        "canonical_switching": canonical_count,
        "count_match": observed_count == canonical_count,
        "exact_gene_and_bracket_match": exact_match,
    }


def run(args) -> Path:
    if not np.isclose(args.tau, DEFAULT_THRESHOLD):
        raise ValueError(
            "Canonical regression requires tau=0.5; use --skip-regression "
            "for an explicitly non-canonical sensitivity analysis."
        )
    if not np.isclose(args.epsilon, EPSILON):
        raise ValueError(
            "Canonical regression requires epsilon=0.01; use "
            "--skip-regression for an explicitly non-canonical analysis."
        )

    metadata = load_metadata(args.version)
    requested = parse_tissue_arguments(args.tissue)
    tissues = select_tissues(metadata, requested=requested)

    rows = []
    for index, tissue in enumerate(tissues, start=1):
        print(f"[{index}/{len(tissues)}] Canonical regression: {tissue}")
        rows.append(_check_tissue(args.version, tissue, args.tau, args.epsilon))

    frame = pd.DataFrame(rows).sort_values("tissue", kind="stable")
    out_dir = validation_root(
        args.version, args.output_root, args.run_label
    ) / "01_regression"
    out_dir.mkdir(parents=True, exist_ok=True)
    result_path = out_dir / "canonical_regression_results.csv"
    frame.to_csv(result_path, index=False)

    failures = frame.loc[
        ~(frame["count_match"] & frame["exact_gene_and_bracket_match"])
    ]
    if not failures.empty:
        raise RuntimeError(
            "Canonical regression failed for: "
            f"{failures['tissue'].tolist()}. Long analyses must not be started."
        )
    print(f"Canonical regression passed for {frame.shape[0]} tissue(s).")
    print(f"Regression results written to: {result_path}")
    return result_path


def build_parser():
    return base_parser(
        "Compare new switching calls with canonical complete-analysis sets."
    )


def main() -> None:
    run(build_parser().parse_args())


if __name__ == "__main__":
    main()
