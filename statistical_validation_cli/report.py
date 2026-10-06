"""Collect validation results and plot summary figures."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from stamp.io import load_metadata

from .common import base_parser, select_tissues, validation_root


def _read_required(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(
            f"Required result not found: {path}. Run the preceding CLI step."
        )
    return pd.read_csv(path)


def _merge_outputs(root: Path) -> pd.DataFrame:
    primary = _read_required(
        root / "02_primary" / "primary_permutation_results.csv"
    )
    merged = primary.copy()

    smote_path = root / "03_smote" / "smote_sensitivity_results.csv"
    if smote_path.exists():
        smote = pd.read_csv(smote_path)
        smote_columns = [
            "tissue",
            "smote_applied",
            "smote_target_policy",
            "smote_target_n",
            "smote_targeted_bins",
            "smote_sampling_strategy",
            "singleton_bins_not_oversampled",
            "k_neighbors",
            "n_smote_seeds",
            "smote_mean",
            "smote_std",
            "smote_median",
            "smote_q025",
            "smote_q975",
            "jaccard_mean",
            "jaccard_min",
            "jaccard_max",
            "smote_test_observed",
            "smote_empirical_p_value",
            "smote_q_value_bh",
            "smote_significant_bh_0_05",
            "smote_bh_eligible",
            "smote_bh_family_size",
            "n_smote_permutations",
        ]
        available = [column for column in smote_columns if column in smote]
        merged = merged.merge(smote[available], on="tissue", how="left")

    bootstrap_path = root / "04_bootstrap" / "bootstrap_stability_results.csv"
    if bootstrap_path.exists():
        bootstrap = pd.read_csv(bootstrap_path)
        bootstrap_columns = [
            "tissue",
            "normalization_mode",
            "n_bootstrap",
            "bootstrap_mean",
            "bootstrap_std",
            "bootstrap_median",
            "bootstrap_q025",
            "bootstrap_q975",
            "observed_gene_exact_stability_mean",
            "observed_genes_stable_0_80",
            "observed_genes_stable_0_90",
        ]
        available = [column for column in bootstrap_columns if column in bootstrap]
        merged = merged.merge(bootstrap[available], on="tissue", how="left")

    return merged


def _plot_primary(frame: pd.DataFrame, path: Path) -> None:
    ordered = frame.sort_values(
        "z_score_descriptive", ascending=True, na_position="first"
    ).copy()
    colors = np.where(ordered["significant_bh_0_05"], "#1976D2", "#9E9E9E")
    y = np.arange(ordered.shape[0])
    plt.figure(figsize=(10, max(8, ordered.shape[0] * 0.23)))
    plt.hlines(
        y,
        ordered["null_mean"],
        ordered["observed_switching"],
        color=colors,
        linewidth=1.5,
    )
    plt.scatter(
        ordered["null_mean"],
        y,
        color="#616161",
        s=22,
        label="Null mean",
        zorder=3,
    )
    plt.scatter(
        ordered["observed_switching"],
        y,
        color=colors,
        s=28,
        label="Observed",
        zorder=3,
    )
    plt.yticks(y, ordered["tissue"], fontsize=7)
    plt.xlabel("Number of switching genes")
    plt.title("Real-donor STAMP counts versus permutation null")
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(path, dpi=250, bbox_inches="tight")
    plt.close()


def _plot_bootstrap(frame: pd.DataFrame, path: Path) -> None:
    required = {"bootstrap_q025", "bootstrap_q975", "bootstrap_median"}
    if not required.issubset(frame.columns):
        return
    plot_frame = frame.dropna(subset=list(required)).copy()
    if plot_frame.empty:
        return
    plot_frame = plot_frame.sort_values("bootstrap_median", ascending=True)
    y = np.arange(plot_frame.shape[0])
    lower = plot_frame["bootstrap_median"] - plot_frame["bootstrap_q025"]
    upper = plot_frame["bootstrap_q975"] - plot_frame["bootstrap_median"]
    plt.figure(figsize=(10, max(8, plot_frame.shape[0] * 0.23)))
    plt.errorbar(
        plot_frame["bootstrap_median"],
        y,
        xerr=np.vstack([lower, upper]),
        fmt="o",
        color="#2E7D32",
        ecolor="#81C784",
        capsize=2,
        markersize=3,
    )
    plt.scatter(
        plot_frame["observed_switching"],
        y,
        marker="x",
        color="#212121",
        s=25,
        label="Observed",
        zorder=3,
    )
    plt.yticks(y, plot_frame["tissue"], fontsize=7)
    plt.xlabel("Number of switching genes")
    plt.title("Stratified-bootstrap stability intervals")
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(path, dpi=250, bbox_inches="tight")
    plt.close()


def _plot_smote(frame: pd.DataFrame, path: Path) -> None:
    if "smote_median" not in frame.columns:
        return
    plot_frame = frame.dropna(subset=["smote_median"]).copy()
    if plot_frame.empty:
        return
    maximum = float(
        max(
            plot_frame["observed_switching"].max(),
            plot_frame["smote_median"].max(),
        )
    )
    plt.figure(figsize=(7, 7))
    scatter = plt.scatter(
        plot_frame["observed_switching"],
        plot_frame["smote_median"],
        c=plot_frame.get("jaccard_mean", pd.Series(0.0, index=plot_frame.index)),
        cmap="viridis",
        vmin=0,
        vmax=1,
        s=45,
        alpha=0.85,
    )
    plt.plot([0, maximum], [0, maximum], "--", color="#616161", linewidth=1)
    plt.xlabel("Switching genes: real donors")
    plt.ylabel("Median switching genes after SMOTE")
    plt.title("Exploratory SMOTE sensitivity")
    colorbar = plt.colorbar(scatter)
    colorbar.set_label("Jaccard with real switching set")
    plt.tight_layout()
    plt.savefig(path, dpi=250, bbox_inches="tight")
    plt.close()


def _write_markdown(frame: pd.DataFrame, path: Path, run_label: str) -> None:
    significant = frame.loc[frame["significant_bh_0_05"]].sort_values(
        "q_value_bh"
    )
    lines = [
        "# STAMP statistical validation summary",
        "",
        f"- Run label: **{run_label}**",
        f"- Tissues in the primary table: **{frame.shape[0]}**",
        f"- BH family size: **{int(frame['bh_family_size'].max())}**",
        (
            "- BH-significant tissues at q <= 0.05: "
            f"**{significant.shape[0]}**"
        ),
        (
            "- Primary Monte Carlo permutations per tissue: "
            f"**{int(frame['n_permutations'].iloc[0])}**"
        ),
        "",
        "The primary claims must use the real-donor permutation p-values and "
        "BH q-values. Bootstrap columns quantify stability and do not represent "
        "null-hypothesis p-values. SMOTE columns are exploratory and refer to "
        "synthetic interpolation, not additional independent donors.",
        "",
        "## BH-significant tissues",
        "",
    ]
    if significant.empty:
        lines.append("None.")
    else:
        lines.extend(
            [
                "| Tissue | Observed | Null mean | p | q |",
                "|---|---:|---:|---:|---:|",
            ]
        )
        for row in significant.itertuples(index=False):
            lines.append(
                f"| {row.tissue} | {row.observed_switching} | "
                f"{row.null_mean:.2f} | {row.empirical_p_value:.6g} | "
                f"{row.q_value_bh:.6g} |"
            )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(args) -> Path:
    root = validation_root(args.version, args.output_root, args.run_label)
    frame = _merge_outputs(root)
    expected_tissues = set(select_tissues(load_metadata(args.version)))
    observed_tissues = set(frame["tissue"].astype(str))
    allow_incomplete = bool(getattr(args, "allow_incomplete_report", False))
    if args.run_label == "final" and not allow_incomplete:
        problems = []
        if observed_tissues != expected_tissues:
            problems.append(
                f"expected {len(expected_tissues)} complete tissues, found "
                f"{len(observed_tissues)}"
            )
        if not frame.get("bh_applied", pd.Series(False, index=frame.index)).all():
            problems.append("BH was not applied to the complete tissue family")
        if frame["bh_family_size"].nunique() != 1 or (
            int(frame["bh_family_size"].iloc[0]) != len(expected_tissues)
        ):
            problems.append("BH family size does not match the complete set")
        if int(frame["n_permutations"].min()) < 9_999:
            problems.append("fewer than 9999 primary permutations were used")
        if frame["q_value_bh"].isna().any():
            problems.append("one or more primary BH q-values are missing")
        if "n_bootstrap" in frame.columns:
            bootstrap_rows = frame["n_bootstrap"].notna()
            if bootstrap_rows.any() and not bootstrap_rows.all():
                problems.append("bootstrap output covers only part of the family")
            if bootstrap_rows.any() and (
                int(frame.loc[bootstrap_rows, "n_bootstrap"].min()) < 1_000
            ):
                problems.append("fewer than 1000 bootstrap replicates were used")
        if "n_smote_seeds" in frame.columns:
            smote_rows = frame["n_smote_seeds"].notna()
            if smote_rows.any() and not smote_rows.all():
                problems.append("SMOTE output covers only part of the family")
            if smote_rows.any() and (
                int(frame.loc[smote_rows, "n_smote_seeds"].min()) < 20
            ):
                problems.append("fewer than 20 SMOTE seeds were used")
        if problems:
            raise RuntimeError(
                "Protected final report refused: " + "; ".join(problems) + ". "
                "Use a non-final --run-label for pilots/smoke tests."
            )
    report_dir = root / "05_report"
    report_dir.mkdir(parents=True, exist_ok=True)
    integrated_path = report_dir / "integrated_validation_results.csv"
    frame.to_csv(integrated_path, index=False)

    _plot_primary(frame, report_dir / "primary_observed_vs_null.png")
    _plot_bootstrap(frame, report_dir / "bootstrap_intervals.png")
    _plot_smote(frame, report_dir / "smote_real_comparison.png")
    _write_markdown(
        frame, report_dir / "validation_summary.md", args.run_label
    )

    print(f"Integrated table written to: {integrated_path}")
    print(f"Figures and Markdown summary written to: {report_dir}")
    return integrated_path


def build_parser():
    parser = base_parser("Merge validation outputs and create compact figures.")
    parser.add_argument(
        "--allow-incomplete-report",
        action="store_true",
        help="Override final-report completeness checks (not recommended).",
    )
    return parser


def main() -> None:
    run(build_parser().parse_args())


if __name__ == "__main__":
    main()
