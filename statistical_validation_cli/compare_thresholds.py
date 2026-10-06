"""Apply significance cutoffs to existing primary and SMOTE results.

Reads the saved Benjamini-Hochberg q-values and writes comparison tables.
Permutation tests and q-values are not recalculated.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence


DEFAULT_THRESHOLDS = (0.01, 0.05, 0.10)


@dataclass(frozen=True)
class DatasetSpec:
    dataset_id: str
    dataset_label: str
    input_path: Path
    p_column: str
    q_column: str
    eligible_column: str
    observed_column: str
    permutations_column: str


def _parse_bool(value: str) -> bool:
    return value.strip().lower() in {"true", "1", "yes"}


def _parse_optional_float(value: str) -> float | None:
    stripped = value.strip()
    return float(stripped) if stripped else None


def _parse_optional_int(value: str) -> int | None:
    stripped = value.strip()
    return int(float(stripped)) if stripped else None


def _format_threshold(value: float) -> str:
    return f"{value:.2f}".replace(".", "_")


def _threshold_label(value: float) -> str:
    return f"q<={value:.2f}"


def _evidence_tier(q_value: float | None, eligible: bool) -> str:
    if not eligible or q_value is None:
        return "not_eligible"
    if q_value <= 0.01:
        return "q_le_0.01"
    if q_value <= 0.05:
        return "0.01_lt_q_le_0.05"
    if q_value <= 0.10:
        return "0.05_lt_q_le_0.10_exploratory"
    return "q_gt_0.10"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(f"Result file not found: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"Result file is empty: {path}")
    return rows


def _write_csv(
    path: Path,
    rows: Sequence[dict[str, Any]],
    fieldnames: Sequence[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _build_detail_rows(
    spec: DatasetSpec,
    thresholds: Sequence[float],
) -> list[dict[str, Any]]:
    input_rows = _read_rows(spec.input_path)
    output_rows: list[dict[str, Any]] = []

    required = {
        "version",
        "tissue",
        spec.p_column,
        spec.q_column,
        spec.eligible_column,
        spec.observed_column,
        spec.permutations_column,
    }
    missing = required.difference(input_rows[0])
    if missing:
        raise ValueError(
            f"{spec.dataset_id}: missing required column(s): {sorted(missing)}"
        )

    for row in input_rows:
        eligible = _parse_bool(row[spec.eligible_column])
        p_value = _parse_optional_float(row[spec.p_column])
        q_value = _parse_optional_float(row[spec.q_column])
        if eligible and q_value is None:
            raise ValueError(
                f"{spec.dataset_id}/{row['tissue']}: eligible row has no q-value"
            )
        if not eligible and q_value is not None:
            raise ValueError(
                f"{spec.dataset_id}/{row['tissue']}: ineligible row has a q-value"
            )

        result: dict[str, Any] = {
            "dataset_id": spec.dataset_id,
            "dataset_label": spec.dataset_label,
            "version": row["version"],
            "tissue": row["tissue"],
            "bh_eligible": eligible,
            "empirical_p_value": p_value,
            "q_value_bh": q_value,
            "evidence_tier": _evidence_tier(q_value, eligible),
            "observed_switching": _parse_optional_float(row[spec.observed_column]),
            "n_permutations": _parse_optional_int(row[spec.permutations_column]),
        }
        for threshold in thresholds:
            suffix = _format_threshold(threshold)
            result[f"significant_q_le_{suffix}"] = bool(
                eligible and q_value is not None and q_value <= threshold
            )
        output_rows.append(result)

    tissues = [row["tissue"] for row in output_rows]
    if len(tissues) != len(set(tissues)):
        raise ValueError(f"{spec.dataset_id}: duplicated tissue names detected")
    return sorted(output_rows, key=lambda row: str(row["tissue"]).lower())


def _build_summary_rows(
    details_by_dataset: dict[str, list[dict[str, Any]]],
    specs: Sequence[DatasetSpec],
    thresholds: Sequence[float],
) -> list[dict[str, Any]]:
    summary: list[dict[str, Any]] = []
    for spec in specs:
        rows = details_by_dataset[spec.dataset_id]
        eligible_count = sum(bool(row["bh_eligible"]) for row in rows)
        for threshold in thresholds:
            flag = f"significant_q_le_{_format_threshold(threshold)}"
            significant_count = sum(bool(row[flag]) for row in rows)
            summary.append(
                {
                    "dataset_id": spec.dataset_id,
                    "dataset_label": spec.dataset_label,
                    "threshold": threshold,
                    "threshold_label": _threshold_label(threshold),
                    "role": (
                        "strict"
                        if threshold == 0.01
                        else "primary"
                        if threshold == 0.05
                        else "exploratory"
                    ),
                    "n_tissues_total": len(rows),
                    "n_tissues_bh_eligible": eligible_count,
                    "n_tissues_significant": significant_count,
                    "proportion_of_eligible": (
                        significant_count / eligible_count if eligible_count else None
                    ),
                }
            )
    return summary


def _build_long_significant_rows(
    details_by_dataset: dict[str, list[dict[str, Any]]],
    specs: Sequence[DatasetSpec],
    thresholds: Sequence[float],
) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for spec in specs:
        for threshold in thresholds:
            flag = f"significant_q_le_{_format_threshold(threshold)}"
            for row in details_by_dataset[spec.dataset_id]:
                if row[flag]:
                    output.append(
                        {
                            "dataset_id": spec.dataset_id,
                            "dataset_label": spec.dataset_label,
                            "threshold": threshold,
                            "threshold_label": _threshold_label(threshold),
                            "tissue": row["tissue"],
                            "empirical_p_value": row["empirical_p_value"],
                            "q_value_bh": row["q_value_bh"],
                            "evidence_tier": row["evidence_tier"],
                        }
                    )
    return output


def _build_tissue_comparison(
    details_by_dataset: dict[str, list[dict[str, Any]]],
    specs: Sequence[DatasetSpec],
    thresholds: Sequence[float],
) -> list[dict[str, Any]]:
    lookup = {
        dataset_id: {str(row["tissue"]): row for row in rows}
        for dataset_id, rows in details_by_dataset.items()
    }
    tissue_sets = {dataset_id: set(rows) for dataset_id, rows in lookup.items()}
    expected = tissue_sets[specs[0].dataset_id]
    for dataset_id, tissue_set in tissue_sets.items():
        if tissue_set != expected:
            missing = sorted(expected.difference(tissue_set))
            extra = sorted(tissue_set.difference(expected))
            raise ValueError(
                f"{dataset_id}: tissue universe mismatch; missing={missing}, extra={extra}"
            )

    comparison: list[dict[str, Any]] = []
    for tissue in sorted(expected, key=str.lower):
        row: dict[str, Any] = {"tissue": tissue}
        for spec in specs:
            item = lookup[spec.dataset_id][tissue]
            prefix = spec.dataset_id
            row[f"{prefix}_bh_eligible"] = item["bh_eligible"]
            row[f"{prefix}_empirical_p_value"] = item["empirical_p_value"]
            row[f"{prefix}_q_value_bh"] = item["q_value_bh"]
            row[f"{prefix}_evidence_tier"] = item["evidence_tier"]
            row[f"{prefix}_observed_switching"] = item["observed_switching"]
            for threshold in thresholds:
                suffix = _format_threshold(threshold)
                row[f"{prefix}_significant_q_le_{suffix}"] = item[
                    f"significant_q_le_{suffix}"
                ]

        for threshold in thresholds:
            suffix = _format_threshold(threshold)
            primary = bool(row[f"primary_significant_q_le_{suffix}"])
            fixed = bool(row[f"smote_fixed_significant_q_le_{suffix}"])
            adaptive = bool(row[f"smote_adaptive_significant_q_le_{suffix}"])
            n_significant = int(primary) + int(fixed) + int(adaptive)
            row[f"n_datasets_significant_q_le_{suffix}"] = n_significant
            row[f"all_three_significant_q_le_{suffix}"] = n_significant == 3
            row[f"either_smote_significant_q_le_{suffix}"] = fixed or adaptive

        comparison.append(row)
    return comparison


def _fieldnames(rows: Sequence[dict[str, Any]]) -> list[str]:
    if not rows:
        raise ValueError("Cannot determine columns from an empty result table")
    return list(rows[0])


def _validate_thresholds(values: Iterable[float]) -> tuple[float, ...]:
    thresholds = tuple(sorted(set(values)))
    if not thresholds or any(value <= 0 or value >= 1 for value in thresholds):
        raise ValueError("Thresholds must be unique values strictly between 0 and 1")
    return thresholds


def run_comparison(
    specs: Sequence[DatasetSpec],
    output_dir: Path,
    thresholds: Sequence[float] = DEFAULT_THRESHOLDS,
) -> dict[str, Path]:
    thresholds = _validate_thresholds(thresholds)
    output_dir.mkdir(parents=True, exist_ok=True)

    details_by_dataset = {
        spec.dataset_id: _build_detail_rows(spec, thresholds) for spec in specs
    }
    summary_rows = _build_summary_rows(details_by_dataset, specs, thresholds)
    significant_rows = _build_long_significant_rows(
        details_by_dataset, specs, thresholds
    )
    comparison_rows = _build_tissue_comparison(
        details_by_dataset, specs, thresholds
    )

    output_paths: dict[str, Path] = {}
    for spec in specs:
        path = output_dir / f"{spec.dataset_id}_significance_by_threshold.csv"
        rows = details_by_dataset[spec.dataset_id]
        _write_csv(path, rows, _fieldnames(rows))
        output_paths[f"{spec.dataset_id}_detail"] = path

    summary_path = output_dir / "threshold_summary.csv"
    _write_csv(summary_path, summary_rows, _fieldnames(summary_rows))
    output_paths["summary"] = summary_path

    significant_path = output_dir / "significant_tissues_long.csv"
    _write_csv(significant_path, significant_rows, _fieldnames(significant_rows))
    output_paths["significant_long"] = significant_path

    comparison_path = output_dir / "three_dataset_tissue_comparison.csv"
    _write_csv(comparison_path, comparison_rows, _fieldnames(comparison_rows))
    output_paths["comparison"] = comparison_path

    metadata = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "operation": (
            "Reclassification of pre-existing Benjamini-Hochberg q-values; "
            "no permutation test and no BH adjustment rerun."
        ),
        "thresholds": list(thresholds),
        "threshold_roles": {
            "0.01": "strict",
            "0.05": "primary",
            "0.10": "exploratory",
        },
        "datasets": [
            {
                "dataset_id": spec.dataset_id,
                "dataset_label": spec.dataset_label,
                "source_path": str(spec.input_path.resolve()),
                "source_sha256": _sha256(spec.input_path),
                "n_tissues_total": len(details_by_dataset[spec.dataset_id]),
                "n_tissues_bh_eligible": sum(
                    bool(row["bh_eligible"])
                    for row in details_by_dataset[spec.dataset_id]
                ),
            }
            for spec in specs
        ],
        "notes": [
            "The primary dataset has 50 BH-eligible tissues.",
            (
                "Each SMOTE analysis has 39 BH-eligible tissues; 11 tissues in "
                "which SMOTE was not applied remain explicitly marked not_eligible."
            ),
            (
                "q<=0.05 remains the prespecified primary decision threshold; "
                "q<=0.10 is reported only as exploratory sensitivity analysis."
            ),
        ],
    }
    metadata_path = output_dir / "comparison_metadata.json"
    metadata_path.write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    output_paths["metadata"] = metadata_path

    return output_paths


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Reclassify existing primary and SMOTE BH q-values at multiple "
            "thresholds without rerunning permutations."
        )
    )
    parser.add_argument(
        "--primary-results",
        type=Path,
        default=Path(
            "output/v10/statistical_validation_cli/runs/final/"
            "02_primary/primary_permutation_results.csv"
        ),
    )
    parser.add_argument(
        "--fixed-results",
        type=Path,
        default=Path(
            "output/v10/statistical_validation_cli/runs/"
            "smote-permutation-fixed-all-9999/03_smote/"
            "smote_sensitivity_results.csv"
        ),
    )
    parser.add_argument(
        "--adaptive-results",
        type=Path,
        default=Path(
            "output/v10/statistical_validation_cli/runs/"
            "smote-permutation-adaptive-all-9999/03_smote/"
            "smote_sensitivity_results.csv"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(
            "output/v10/statistical_validation_cli/runs/"
            "significance-threshold-comparison"
        ),
    )
    parser.add_argument(
        "--thresholds",
        type=float,
        nargs="+",
        default=list(DEFAULT_THRESHOLDS),
        help="BH q-value decision thresholds (default: 0.01 0.05 0.10).",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    specs = (
        DatasetSpec(
            dataset_id="primary",
            dataset_label="Original real-donor data",
            input_path=args.primary_results,
            p_column="empirical_p_value",
            q_column="q_value_bh",
            eligible_column="bh_applied",
            observed_column="observed_switching",
            permutations_column="n_permutations",
        ),
        DatasetSpec(
            dataset_id="smote_fixed",
            dataset_label="SMOTE fixed floor (minimum 20)",
            input_path=args.fixed_results,
            p_column="smote_empirical_p_value",
            q_column="smote_q_value_bh",
            eligible_column="smote_bh_eligible",
            observed_column="smote_test_observed",
            permutations_column="n_smote_permutations",
        ),
        DatasetSpec(
            dataset_id="smote_adaptive",
            dataset_label="SMOTE adaptive (second-smallest class)",
            input_path=args.adaptive_results,
            p_column="smote_empirical_p_value",
            q_column="smote_q_value_bh",
            eligible_column="smote_bh_eligible",
            observed_column="smote_test_observed",
            permutations_column="n_smote_permutations",
        ),
    )
    output_paths = run_comparison(
        specs=specs,
        output_dir=args.output_dir,
        thresholds=args.thresholds,
    )
    print("Threshold comparison completed without rerunning permutations.")
    for label, path in output_paths.items():
        print(f"{label}: {path.resolve()}")


if __name__ == "__main__":
    main()
