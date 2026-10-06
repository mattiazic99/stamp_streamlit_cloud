# Reproduce the manuscript's statistical validation

This guide targets the **Statistical validation** section of the STAMP manuscript.
The historical source is recorded as `main (24).tex` in the reference manifest.
Repository: https://github.com/mattiazic99/stamp_streamlit_cloud.
The reference manifest records that manuscript's SHA-256 and the exact numerical
parameters. The reproduction command is fixed to **GTEx v10**.

## What is checked

| Experiment | Exact configuration | Manuscript result |
|---|---|---|
| Real-donor permutation test | 50 age-complete tissues; seed 42; tau 0.5; epsilon 0.01; 9,999 permutations per tissue; BH over 50 tissues | 4 tissues at q <= 0.01, 13 at q <= 0.05, 18 at q <= 0.10; all 18 published table rows |
| Fixed-floor SMOTE + STAMP | Minimum 20 profiles; singleton brackets unchanged; 20 seeds; 9,999 permutations, with SMOTE rerun inside each permutation | 39 modified tissues; 11 unchanged; BH over the 39 modified tissues |
| Adaptive SMOTE + STAMP | Raise an isolated undersized bracket to the second-smallest; otherwise use the fixed floor; same seed and repetitions | 3 original tissues confirmed by both policies at q <= 0.05; 7 of the original 18 confirmed by both at q <= 0.10, one more by adaptive only |
| Stratified bootstrap | All 50 tissues; seed 42; 1,000 replicates; original filtering and normalization fixed | 2,782 of 24,649 primary gene-tissue associations retain bracket and direction in >= 80% of replicates |
| Threshold sensitivity | The 13 primary tissues; tau 0.40, 0.45, 0.55, 0.60 against 0.50 | Median gene Jaccard 0.15, 0.33, 0.30, 0.05; gene-and-bracket Jaccard 0.01, 0.09, 0.08, 0.00 |

`statistical_validation_cli.run_all` is a general-purpose workflow whose default
SMOTE analysis is descriptive. **Use `statistical_validation_cli.paper run` to
reproduce the manuscript**, including the two compound permutation tests.
The manuscript command specifies every numerical setting and checks every tissue,
every numerical table column, complete Monte Carlo sequences and bootstrap
gene-level stability frequencies, observed gene/bracket/direction lists and
SMOTE gene-level frequencies. Matching only the number of significant tissues
is insufficient.

## 1. Create a clean numerical environment

Run these commands from the repository root. Use **Python 3.11**, with no inherited
system packages. Install the lock file without mixing the application environment
into this environment. The CLI runs from the checkout; no editable installation
of the application or Cassandra connection is needed.

Windows PowerShell:

```powershell
py -3.11 -m venv .venv-reproduction
$reproPython = (Resolve-Path .venv-reproduction\Scripts\python.exe).Path
& $reproPython -m pip install -r requirements-reproduction.lock.txt
& $reproPython -m statistical_validation_cli.paper check-reference
```

Linux/macOS:

```bash
python3.11 -m venv .venv-reproduction
source .venv-reproduction/bin/activate
python -m pip install -r requirements-reproduction.lock.txt
python -m statistical_validation_cli.paper check-reference
```

`check-reference` verifies the integrity and consistency of the committed reference
tables. It **does not rerun experiments**. The recorded environment was captured
when this reproduction procedure was prepared; the historical runs did not record
a package manifest. Fresh replay and full verification determine compatibility
with those runs rather than assuming it from package versions alone.

## 2. Quick fresh replay using the bundled Ovary data

The repository contains sample-level Ovary data: 193 donors and the 54,897 genes
retained by the original epsilon filter. Filtering is independent of age labels;
removing discarded genes preserves these experiments on Ovary. The demo is read
in place, so it does not replace your full input files.

```powershell
& $reproPython -m statistical_validation_cli.paper replay --demo
```

```bash
python -m statistical_validation_cli.paper replay --demo
```

This recomputes the observed count (**3,300**), up/down counts, the first 25 primary
permutations, the first 25 bootstrap counts, the first two SMOTE seed runs, and the
first three compound permutations under each policy. It compares each value with
the beginning of the historical full runs, without using cached checkpoints.
It verifies executable reproducibility on one tissue. It does **not** reproduce
the 50-tissue BH family or replace the full statistical analysis.

## 3. Obtain the exact complete inputs

The committed atlas contains six age-bracket averages. Permutations, SMOTE and
bootstrap require **sample-level TPM and donor metadata**, which cannot be
reconstructed from those averages.

`inputs.json` identifies the two exact original Parquet files by size and SHA-256.
The full matrix is approximately 4.01 GB. Its gene order, sample IDs, values and
donor-to-age assignments are part of the reference input.

The accompanying input bundle has eight matrix parts of at most 512 MiB and one
metadata part, listed in `input_bundle.json`. These large parts are excluded from
Git. They must be made available as release assets alongside the code release.
Download **all nine `.partNN` files** to one directory, for example `input-parts`.
Restore without overwriting different existing data:

```powershell
& $reproPython -m statistical_validation_cli.input_bundle --bundle-dir input-parts
& $reproPython -m statistical_validation_cli.paper check-inputs
```

```bash
python -m statistical_validation_cli.input_bundle --bundle-dir input-parts
python -m statistical_validation_cli.paper check-inputs
```

Restore checks each part and then checks the complete files. Corrupted, missing
or different inputs cause an error. If you already have the original Parquet
files, place them under `data/parquet/v10/` and run `check-inputs` directly.
Use `--input-dir /absolute/path/to/v10` to read an existing input directory.

**Publication status:** the parts are prepared locally for the forthcoming
release. This guide must be updated with the actual release-download location
when those assets are published. A code-only GitHub push provides the demo,
reference checks and scripts, but cannot provide the complete input bundle.

The general pipeline can build Parquet files from public GTEx downloads.
Until those reconstructed inputs have been compared with the historical input,
that route is not certified to return these exact Monte Carlo results. This
manuscript reproduction uses the checked original inputs.

## 4. Fresh replay on all 50 tissues

```powershell
& $reproPython -m statistical_validation_cli.paper replay --all-tissues
```

```bash
python -m statistical_validation_cli.paper replay --all-tissues
```

This runs the same sequence-prefix checks as the demo across the complete atlas.
Replay is useful before committing to the much more expensive full run. It still
does not recompute all 9,999 permutations. `--tissue "Liver"` selects a replay
tissue; the complete `run` command deliberately has no subset mode.

## 5. Recompute and verify the complete manuscript experiments

```powershell
& $reproPython -m statistical_validation_cli.paper run --jobs 2
```

```bash
python -m statistical_validation_cli.paper run --jobs 2
```

The command checks input hashes, installed package versions and the 50-tissue
universe, runs the audit and canonical regression, then executes:

1. The original-donor test with 9,999 permutations and BH over all 50 tissues.
2. The fixed-normalization bootstrap with 1,000 replicates on all 50 tissues.
3. The fixed-floor and adaptive SMOTE tests, each with 20 seeds and 9,999 compound
   permutations per modified tissue, BH over its 39 eligible tissues.
4. Threshold sensitivity recomputed from the checked sample-level inputs;
   the 13 tissues are selected from the fresh primary BH results.
5. Comparison with the manuscript references and an integrated report.

Both compound SMOTE tests recompute SMOTE inside every permutation and are
substantially more expensive than replay. `--jobs 2` limits simultaneous tissue
jobs; adjust it to available RAM. Numerical library threads are limited to one
per worker. Stable seeds are derived separately per tissue and experiment.

Results go into the separate directory `output/reproduction/`:

```text
paper_provenance.json
v10/statistical_validation_cli/runs/paper/02_primary/
v10/statistical_validation_cli/runs/paper/04_bootstrap/
v10/statistical_validation_cli/runs/paper-smote-fixed/03_smote/
v10/statistical_validation_cli/runs/paper-smote-adaptive/03_smote/
paper-threshold/
paper_verified.json
```

Per-tissue checkpoints allow resuming with identical inputs, code and environment.
Input, source and environment hashes are included in each worker checkpoint.
An output directory without provenance is rejected if it is non-empty.
`paper_verified.json` lists recomputed and reused tissues for each experiment.
The bundled atlas used by regression is also hashed.
If provenance changes, the command refuses reuse; choose a new `--output-root` or
explicitly pass `--force` to recompute. Reference files are never overwritten.
`paper_verified.json` is written only after the full comparisons pass.

## 6. Verify an existing full reproduction

```powershell
& $reproPython -m statistical_validation_cli.paper verify --output-root output/reproduction
```

```bash
python -m statistical_validation_cli.paper verify --output-root output/reproduction
```

Count and sequence checks are exact. Continuous table summaries allow only a
relative tolerance of `1e-12`, with no absolute tolerance; missingness and tissue
identities must match. Bootstrap gene identities and frequencies are checked
with a deterministic digest. Reference SHA-256 values refer to artifacts;
Parquet/ZIP byte layouts are not scientific results, so Monte Carlo arrays and
gene-level frequencies are compared after decoding.

The published table rounds p-values to four decimals and q-values to three.
The machine-readable reference retains full precision for verification. It also
records which three tissues were confirmed by SMOTE, which seven/eighth tissue
were recovered at the exploratory threshold, and the six unchanged tissues.

These are tissue-level permutation tests. Bootstrap measures gene-call
stability, and SMOTE is a secondary sensitivity analysis. Neither adds independent
donors nor enlarges the primary 13-tissue discovery set.

Exact fresh replay has been reported on two different Intel CPU configurations.
[GitHub CI](https://github.com/mattiazic99/stamp_streamlit_cloud/actions) has passed
the reference checks and fresh demo on Ubuntu and Windows, plus the interface
container checks. Complete cross-platform numerical identity and a fresh full
50-tissue reproduction have not yet been established.
