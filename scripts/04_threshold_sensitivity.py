"""Measure switching-call sensitivity to the expression threshold tau.

For each primary tissue, compare calls at each threshold with tau = 0.5 using
Jaccard overlap of gene IDs and of (gene, age-bracket) pairs. The latter requires
the switching bracket to match as well as the gene.

Default mode reads output/{version}_complete/normalized/*.parquet. The manuscript
runner uses --sample-level --primary-summary to rebuild age means from checked
inputs and select primary tissues from the newly calculated complete BH family.

Writes per_tissue.csv (one row per tissue and tau) and summary.csv (median per tau)
under output/{version}_complete/threshold_sensitivity/ by default.

Examples:
    python scripts/04_threshold_sensitivity.py --version v10
    python scripts/04_threshold_sensitivity.py --version v10 --taus 0.4 0.45 0.55 0.6
    python scripts/04_threshold_sensitivity.py --version v10 --all-tissues
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from stamp.config import AGE_BRACKETS, DEFAULT_THRESHOLD, SWITCHING_BRACKETS, paths_for  # noqa: E402
from stamp.io import load_normalized_tissue  # noqa: E402
from stamp.switching import identify_switching_genes  # noqa: E402

# Primary tissues at FDR 5% under the age-label permutation null.
# Use --all-tissues to include the rest of the atlas.
PRIMARY_TISSUES: tuple[str, ...] = (
    "Artery - Tibial",
    "Colon - Sigmoid",
    "Nerve - Tibial",
    "Colon - Transverse",
    "Brain - Cerebellar Hemisphere",
    "Liver",
    "Thyroid",
    "Whole Blood",
    "Lung",
    "Ovary",
    "Artery - Aorta",
    "Muscle - Skeletal",
    "Brain - Cerebellum",
)

DEFAULT_TAUS: tuple[float, ...] = (0.40, 0.45, 0.55, 0.60)


def _jaccard(a: set, b: set) -> float:
    """Jaccard index; 0.0 when both sets are empty."""
    union = a | b
    return len(a & b) / len(union) if union else 0.0


def _gene_set(sets: dict[str, list[str]]) -> set[str]:
    return {gene for bracket in SWITCHING_BRACKETS for gene in sets.get(bracket, [])}


def _gene_bracket_set(sets: dict[str, list[str]]) -> set[tuple[str, str]]:
    return {
        (gene, bracket)
        for bracket in SWITCHING_BRACKETS
        for gene in sets.get(bracket, [])
    }


def sensitivity_for_tissue(
    version: str, tissue: str, taus: tuple[float, ...], reference_tau: float,
    sample_level: bool = False, epsilon: float = 0.01,
) -> list[dict]:
    """Compare the calls at each tau against the reference tau, for one tissue."""
    if sample_level:
        from statistical_validation_cli.common import load_tissue_data, age_means_from_codes
        from threadpoolctl import threadpool_limits
        data = load_tissue_data(version, tissue, epsilon=epsilon)
        with threadpool_limits(limits=1):
            means = age_means_from_codes(data.normalized_values, data.age_codes)
        df = pd.DataFrame(means, index=data.gene_ids, columns=AGE_BRACKETS)
    else:
        df = load_normalized_tissue(version, tissue, complete=True)

    reference = identify_switching_genes(df, reference_tau)
    ref_genes = _gene_set(reference)
    ref_pairs = _gene_bracket_set(reference)

    rows = []
    for tau in taus:
        current = identify_switching_genes(df, tau)
        genes = _gene_set(current)
        pairs = _gene_bracket_set(current)
        rows.append(
            {
                "tissue": tissue,
                "tau": tau,
                "n_switching": len(genes),
                "n_switching_reference": len(ref_genes),
                "jaccard_gene": _jaccard(genes, ref_genes),
                "jaccard_gene_bracket": _jaccard(pairs, ref_pairs),
            }
        )
    return rows


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", choices=["v8", "v10"], default="v10")
    parser.add_argument(
        "--taus",
        type=float,
        nargs="+",
        default=list(DEFAULT_TAUS),
        help="Thresholds to compare against the reference (default: 0.4 0.45 0.55 0.6).",
    )
    parser.add_argument(
        "--reference-tau",
        type=float,
        default=DEFAULT_THRESHOLD,
        help="Threshold the others are compared against (default: 0.5).",
    )
    parser.add_argument(
        "--all-tissues",
        action="store_true",
        help="Run over every age-complete tissue instead of the 13 primary ones.",
    )
    parser.add_argument("--output-dir", type=Path, default=None, help="Separate directory for reproduction outputs.")
    parser.add_argument("--sample-level", action="store_true", help="Recompute averages from sample-level TPM inputs.")
    parser.add_argument("--epsilon", type=float, default=0.01)
    parser.add_argument("--primary-summary", type=Path, help="Select q <= 0.05 tissues from a complete primary run.")
    return parser


def main() -> None:
    args = build_parser().parse_args()

    out_dir = args.output_dir or paths_for(args.version, complete=True)["threshold_sensitivity"]
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.primary_summary:
        primary = pd.read_csv(args.primary_summary)
        if (primary.shape[0] != 50 or primary["tissue"].duplicated().any()
                or primary["q_value_bh"].isna().any()
                or not (primary["bh_family_size"] == 50).all()):
            raise ValueError("Threshold selection requires the complete 50-tissue BH family")
        tissues = sorted(primary.loc[primary.q_value_bh <= 0.05, "tissue"])
        if set(tissues) != set(PRIMARY_TISSUES):
            raise ValueError("Fresh primary discoveries differ from the manuscript's 13 tissues")
    elif args.all_tissues:
        # List age-complete tissues from the atlas, without loading raw inputs.
        normalized_dir = paths_for(args.version, complete=True)["normalized"]
        tissues = sorted(p.stem.replace("_", " ") for p in normalized_dir.glob("*.parquet"))
    else:
        tissues = list(PRIMARY_TISSUES)

    taus = tuple(args.taus)
    print(f"GTEx {args.version} — {len(tissues)} tissues, tau in {list(taus)}, "
          f"reference tau = {args.reference_tau}")

    rows: list[dict] = []
    for i, tissue in enumerate(tissues, 1):
        print(f"  [{i:>2}/{len(tissues)}] {tissue}", flush=True)
        rows.extend(
            sensitivity_for_tissue(args.version, tissue, taus, args.reference_tau,
                                   args.sample_level, args.epsilon)
        )

    per_tissue = pd.DataFrame(rows)
    summary = (
        per_tissue.groupby("tau")
        .agg(
            n_tissues=("tissue", "size"),
            median_jaccard_gene=("jaccard_gene", "median"),
            median_jaccard_gene_bracket=("jaccard_gene_bracket", "median"),
            median_n_switching=("n_switching", "median"),
        )
        .reset_index()
    )

    per_tissue_path = out_dir / "per_tissue.csv"
    summary_path = out_dir / "summary.csv"
    per_tissue.to_csv(per_tissue_path, index=False)
    summary.to_csv(summary_path, index=False)

    print("\nMedian Jaccard against tau = "
          f"{args.reference_tau} across {len(tissues)} tissues:\n")
    print(f"  {'tau':>6}  {'gene identity':>14}  {'gene + bracket':>15}")
    for row in summary.itertuples():
        print(f"  {row.tau:>6.2f}  {row.median_jaccard_gene:>14.2f}  "
              f"{row.median_jaccard_gene_bracket:>15.2f}")

    print(f"\nWritten:\n  {per_tissue_path}\n  {summary_path}")


if __name__ == "__main__":
    main()
