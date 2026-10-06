"""Run the validation stages in sequence using their command-line options."""

from __future__ import annotations

import argparse
from pathlib import Path
from types import SimpleNamespace

from . import audit, bootstrap, primary, regression, report, smote


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the complete STAMP statistical validation workflow."
    )
    parser.add_argument("--version", choices=["v8", "v10"], default="v10")
    parser.add_argument("--output-root", type=Path, default=None)
    parser.add_argument(
        "--run-label",
        default=None,
        help="Output label; defaults to final, or smoke with --quick.",
    )
    parser.add_argument("--tissue", action="append", default=None)
    parser.add_argument("--tau", type=float, default=0.5)
    parser.add_argument("--epsilon", type=float, default=0.01)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--primary-permutations", type=int, default=9_999)
    parser.add_argument("--smote-seeds", type=int, default=20)
    parser.add_argument("--smote-target-n", type=int, default=20)
    parser.add_argument(
        "--smote-target-policy",
        choices=["fixed-floor", "adaptive-second-smallest"],
        default="fixed-floor",
    )
    parser.add_argument(
        "--smote-permutations",
        type=int,
        default=0,
        help="0 keeps SMOTE descriptive; 9999 runs the optional compound test.",
    )
    parser.add_argument("--bootstrap-replicates", type=int, default=1_000)
    parser.add_argument(
        "--bootstrap-normalization",
        choices=["fixed", "refit"],
        default="fixed",
    )
    parser.add_argument("--skip-regression", action="store_true")
    parser.add_argument("--skip-smote", action="store_true")
    parser.add_argument("--skip-bootstrap", action="store_true")
    parser.add_argument(
        "--allow-partial-family",
        action="store_true",
        help="Apply BH to an explicitly selected pilot family.",
    )
    parser.add_argument(
        "--allow-expensive-smote-test",
        action="store_true",
        help="Override the compound SMOTE runtime guard.",
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help=(
            "Smoke-test settings: 9 primary permutations, 2 SMOTE seeds and "
            "5 bootstrap replicates. Tissue selection is unchanged."
        ),
    )
    return parser


def run(args) -> None:
    run_label = args.run_label or ("smoke" if args.quick else "final")
    if args.quick and run_label == "final":
        raise ValueError(
            "--quick cannot write into the protected final run. "
            "Use --run-label smoke or omit --run-label."
        )
    primary_permutations = 9 if args.quick else args.primary_permutations
    smote_seeds = 2 if args.quick else args.smote_seeds
    smote_permutations = 0 if args.quick else args.smote_permutations
    bootstrap_replicates = 5 if args.quick else args.bootstrap_replicates

    common = {
        "version": args.version,
        "output_root": args.output_root,
        "run_label": run_label,
        "tissue": args.tissue,
        "tau": args.tau,
        "epsilon": args.epsilon,
        "seed": args.seed,
        "jobs": args.jobs,
        "force": args.force,
    }

    print("\n=== STEP 01 / AUDIT ===")
    audit_code = audit.run(SimpleNamespace(**common))
    if audit_code != 0:
        raise RuntimeError("Input audit failed; validation stopped.")

    if not args.skip_regression:
        print("\n=== CANONICAL REGRESSION CHECK ===")
        regression.run(SimpleNamespace(**common))

    print("\n=== STEP 02 / PRIMARY PERMUTATION ===")
    primary.run(
        SimpleNamespace(
            **common,
            permutations=primary_permutations,
            allow_partial_family=args.allow_partial_family,
        )
    )

    if not args.skip_smote:
        print("\n=== STEP 03 / SMOTE SENSITIVITY ===")
        smote.run(
            SimpleNamespace(
                **common,
                seeds=smote_seeds,
                target_n=args.smote_target_n,
                target_policy=args.smote_target_policy,
                permutations=smote_permutations,
                allow_expensive_test=args.allow_expensive_smote_test,
            )
        )

    if not args.skip_bootstrap:
        print("\n=== STEP 04 / BOOTSTRAP STABILITY ===")
        bootstrap.run(
            SimpleNamespace(
                **common,
                replicates=bootstrap_replicates,
                normalization=args.bootstrap_normalization,
            )
        )

    print("\n=== STEP 05 / INTEGRATED REPORT ===")
    report.run(SimpleNamespace(**common))
    print("\nValidation workflow completed.")


def main() -> None:
    run(build_parser().parse_args())


if __name__ == "__main__":
    main()
