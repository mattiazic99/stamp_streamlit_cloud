# Statistical validation CLI

A self-contained pipeline that validates STAMP **at the tissue level**. The
primary result answers one question: is the observed number of switching genes
in a tissue greater than expected when the age brackets of the donors are
permuted at random?

Step-by-step commands, expected results and the isolated-environment setup are
in the [main README](../README.md#reproducing-the-manuscript-statistical-validation).
This file documents the behaviour of the modules themselves.

Run everything from the repository root.

## The three blocks are kept separate

| Module | What it produces | What it does **not** produce |
|--------|------------------|------------------------------|
| `primary` | Permutation p-values and BH q-values | Gene-level inference |
| `smote` | Exploratory sensitivity to age-bracket imbalance | Additional significant tissues |
| `bootstrap` | Stability frequencies over donor resampling | p-values or q-values |

Supporting modules: `audit` (bracket coverage, sample counts, donor
uniqueness), `regression` (checks the pipeline still reproduces the committed
atlas), `report` (tables and figures), `run_all` (the whole chain).

## Run labels

Each run writes under its own label, so an experiment can never overwrite the
published numbers:

```
output/{version}/statistical_validation_cli/runs/final/
output/{version}/statistical_validation_cli/runs/smoke/
output/{version}/statistical_validation_cli/runs/pilot/
```

`final` is **protected**. It is rejected unless it covers every age-complete
tissue, applies Benjamini–Hochberg once to the complete family, and uses at
least 9,999 permutations per tissue. `--quick` switches to `smoke`
automatically. For a manual trial on a few tissues, always pass an explicit
non-final label, for example `--run-label pilot`.

## Partial families

In a partial run the per-tissue p-value is computed but BH is **not** applied:
`q_value_bh` is left empty. This is deliberate — it prevents silently obtaining
`q = p` from a family of size one. `--allow-partial-family` exists only for an
explicitly defined pilot family and still requires a non-final run label.

## Checkpoints

Per-tissue checkpoints let an interrupted run resume, provided the
configuration, input hashes, code hashes and numerical environment are identical.
Old checkpoints without this provenance are recomputed. The manuscript wrapper
rejects a non-empty output directory without provenance unless `--force` is used.
Its success record distinguishes recomputed and reused tissues. `--force` ignores them and recomputes. For memory
and I/O, start at `--jobs 2` or `--jobs 4`; higher parallelism struggles on the
large tissues.

## SMOTE policies

- **fixed-floor** (default): brackets with 2–19 donors are raised to 20;
  brackets already at ≥ 20 are untouched; singleton brackets cannot be
  interpolated and are left alone; empty brackets are already excluded.
- **adaptive-second-smallest**: if only the smallest bracket is below
  `target_n`, it is raised to the size of the second smallest; if at least two
  are below, the fixed floor applies.

SMOTE interpolates across roughly 54,000 genes with few samples in the
minority classes, where nearest-neighbour distances are fragile. Synthetic
profiles are interpolations, not independent donors, so this stage is treated
as a secondary sensitivity analysis.

The composed SMOTE + STAMP permutation test recomputes SMOTE **inside** every
permutation and is very expensive. Restrict it to 2–4 representative tissues
and 999 permutations. Runs beyond 1,000 permutations on more than four tissues
are blocked without an explicit override, in the general-purpose CLI. The manuscript reproduction explicitly runs
9,999 compound permutations under both policies on all modified tissues;
use `paper run` for that procedure.

## Bootstrap normalization modes

- `fixed` — conditions on the original scale.
- `refit` — recomputes the ε filter and the min–max rescaling, measuring the
  instability of the entire pipeline including the effect of repeated donors.
  Heavier; run it under a separate label.

## Exact manuscript reproduction

See [the reproduction guide](../reproducibility/README.md) and use
`python -m statistical_validation_cli.paper run --jobs 2`. It includes both
full compound SMOTE policies and verifies every published statistical result.
The ordinary `run_all` defaults to descriptive SMOTE.
