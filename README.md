# STAMP

**Spatio-Temporal Analysis and Mapping of gene-expression Patterns**

STAMP identifies and compares *candidate switching genes* in age-stratified
transcriptomic data. For each gene and tissue it normalizes expression,
summarizes it across the six ordered GTEx age brackets, and keeps the profiles
that cross an activation threshold exactly once — making the timing of a
persistent transcriptional change explicit and directly comparable across
tissues.

This repository contains everything behind the STAMP application note: the
analysis library, the interactive web application, and the statistical
validation pipeline.

> **Live application:** https://stampgui.streamlit.app/

Because GTEx is cross-sectional, STAMP describes population-level age-stratified
patterns. Its calls are **candidate events**, not longitudinal trajectories or
causal transitions.

---


## Reproducing the manuscript statistical validation

The exact procedure for the manuscript, including both 9,999-permutation SMOTE
policies, is documented in [reproducibility/README.md](reproducibility/README.md).
It provides a clean pinned environment, a fresh Ovary replay, complete input
checksums, machine-readable reference results and an automatic numerical verifier.
The complete original input parts are prepared separately for release publication;
they are not included in the code-only checkout.

```bash
python -m pip install -r requirements-reproduction.lock.txt
python -m statistical_validation_cli.paper check-reference
python -m statistical_validation_cli.paper replay --demo
```

Use a dedicated Python 3.11 environment as described in the guide. Use
`python -m statistical_validation_cli.paper run --jobs 2` for the full manuscript
analysis after restoring the exact inputs. The general `run_all` command has
different default SMOTE settings and does not reproduce the compound tests.


## Contents

| Path | What it is |
|------|------------|
| `gui/` | The Streamlit application; one module per page in `gui/analysis_modules/` |
| `stamp/` | The analysis library — normalization, switching detection, Jaccard, I/O |
| `scripts/` | The pipeline, one numbered script per stage |
| `statistical_validation_cli/` | Permutation test, SMOTE sensitivity, stratified bootstrap, report |
| `output/v10_complete/`, `output/v8_complete/` | Pre-computed atlas the application reads |
| `data/demo/` | One-tissue sample-level subset, so the validation runs with no download |
| `tests/` | Unit tests for the analysis library |
| `validation/` | Unfinished R comparison stubs; not an executed validation |

The numerical core — range filtering, per-gene min–max normalization,
age-bracket averaging, binarization and switching detection — lives in `stamp/`
and is shared by the interface, the pipeline and the validation modules, so all
three apply exactly the same rule.

The repository is scoped to what the application note describes. The
Alzheimer's disease case study of Section 5 is **not** included: the functional
enrichment, module assignment and biological interpretation were done with
external tools (g:Profiler) and are not part of this codebase.

### A note on the two distributions

This is the **full** repository: interface *and* the code to reproduce the
analysis. A separate, interface-only Docker distribution is also provided for
readers who just want to run the application; it contains `gui/`, `stamp/`,
`output/` and the Docker files, and none of the pipeline or validation code.
If you are holding that one, the sections
[Quick start](#quick-start), [Docker](#running-with-docker) and
[Using the application](#using-the-application) apply and the rest does not.

---

## Quick start

Use **Python 3.11**, the version tested with the pinned application environment.

```bash
git clone https://github.com/mattiazic99/stamp_streamlit_cloud.git
cd stamp_streamlit_cloud

python3.11 -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate

pip install -r requirements.txt
pip install -e .

streamlit run gui/main.py
```

The application opens at <http://localhost:8501>. Nothing to download: the
pre-computed atlas ships with the repository.

<details>
<summary>Windows (PowerShell)</summary>

```powershell
git clone https://github.com/mattiazic99/stamp_streamlit_cloud.git
Set-Location stamp_streamlit_cloud

py -3.11 -m venv .venv
.\.venv\Scripts\Activate.ps1

pip install -r requirements.txt
pip install -e .

streamlit run gui/main.py
```
</details>

---

## Running with Docker

Docker runs only the updated web interface, with its normalized GTEx v10 atlas
and the GTEx v8 data used by Panel Explorer. It does not execute the manuscript
statistical validation. The Python 3.11 base image is pinned by digest, and the interface dependencies are fixed in
`requirements-interface.lock.txt`.

Install and start Docker Desktop with Linux containers, then run from this directory:

```bash
docker compose up --build --wait    # http://localhost:8501
docker compose logs --tail 50 stamp # inspect startup errors
docker compose down                # stop; session downloads are discarded
```

The atlas is embedded in the image and remains read-only. Generator exports and
ZIP archives are held only in each user's Streamlit session memory: the interface
does not write them to disk, mount a downloads volume or send them to GitHub.
Download the ZIP to save a copy on your own computer. **Start New Generation**
clears the application's current exports and ZIP before the page renders;
generating again replaces the previous run. A download keeps the current exports
available for further selections.

Reloading the browser starts fresh session state, but closing or reloading a tab
does not guarantee immediate reclamation of its old data. Streamlit may retain
disconnected session data and download buffers in server memory until session
cleanup and subsequent server activity. Streamlit also handles reclamation of its
own download buffers after a reset. Restarting the server clears all in-memory
exports. These buffers are never written to application storage or GitHub.

After changing the interface, repeat `docker compose up --build --wait` to rebuild
and replace the container. Older installations may still have an unused
`stamp-interface_generated-files` Docker volume; this configuration no longer
mounts or writes to it. Existing files in that old volume are left untouched.

The Home page offers the statistical reproduction guide and dependency lock as
downloads. To run those experiments, use the full repository and a separate
Python environment following [reproducibility/README.md](reproducibility/README.md).
The statistical CLI, raw data, reference results and local research pages are excluded
from the Docker build context and image.

Without Compose:

```bash
docker build -t stamp-interface .
docker run --rm -p 127.0.0.1:8501:8501 stamp-interface
```

---

### Check the container after building

```powershell
Get-Content -Raw scripts/check_interface_container.py | docker compose exec -T stamp python -
```

```bash
docker compose exec -T stamp python - < scripts/check_interface_container.py
```

This checks every page, generations at tau 0.7 and 0.5 in separate sessions,
unchanged atlas hashes, read-only atlas permissions and excluded components.
GitHub CI includes this container check on Ubuntu.

## Deploying to Streamlit Community Cloud

1. Push this repository to GitHub.
2. Go to <https://share.streamlit.io> → **New app** → select the repository.
3. Set **Main file path** to `streamlit_app.py` (or `gui/main.py`; both work).
4. In **Advanced settings**, select **Python 3.11**.
5. Deploy. `requirements.txt` is installed automatically.

For an existing app, verify the deployed Python version in its build logs before
pushing the pinned requirements. Python 3.13 is not the tested environment for
this application lock. Changing Python for an existing Community Cloud app
requires deleting and redeploying the app; changing repository files alone does
not change its interpreter. See the [official Python-version instructions](https://docs.streamlit.io/deploy/streamlit-community-cloud/manage-your-app/upgrade-python).

No secrets and no external services are required. The app reads only the
pre-computed files under `output/{version}_complete/`, which are committed, so
the deployment works from a fresh clone.

---

## Using the application

The interface is organised around the three-step workflow stated on the landing
page: **generate** the switching-gene files, **load** them into an analysis
page, **explore and download** the results.

The exchange format is a plain-text file of five space-separated gene lists,
one per switching bracket, in chronological order — 30–39, 40–49, 50–59, 60–69,
70–79. This is what STAMP produces and what every analysis page accepts:

```
APOE TP53 BRCA1 EGFR
MYC PTEN RB1 VHL
APC KRAS PIK3CA IDH1
CDKN2A ATM SMAD4
MLH1 MSH2 MSH6 PMS2
```

| Page | Question it answers |
|------|--------------------|
| **STAMP Generator** | Produce the five-line files, for one, several or all tissues, at a chosen threshold |
| **Upload & Single Analysis** | How many genes switch in each bracket, for one tissue |
| **Tissue Comparison** | Which genes two tissues share, and when, bracket by bracket |
| **Multi-Tissue Analysis** | Jaccard similarity matrix and hierarchical clustering over many tissues |
| **Age-Specific Analysis** | Which tissues are active in a given bracket |
| **Gene Sharing Analysis** | Pairwise overlap statistics, most similar pairs, arbitrary intersections |
| **Single Gene Analysis** | Where and when a given gene is called switching |
| **Group Comparison** | Two user-defined groups of tissues, compared bracket by bracket |
| **Panel Explorer** | Map a gene panel onto the atlas, with GTEx v8/v10 concordance |

Every table exports to CSV and every figure to PNG.

### Dataset

The application is pinned to **GTEx v10** and to the **50 age-complete
tissues** — those with at least one donor in each of the six age brackets.
GTEx v10 provides 19,788 samples from 948 donors over 54 tissues; in four of
them (Cervix–Ectocervix, Cervix–Endocervix, Fallopian Tube, Kidney–Medulla) at
least one bracket has no donor at all, so the six-value profile cannot be
completed and those tissues are excluded rather than imputed.

**Panel Explorer is the only page that reads GTEx v8**, and it does so purely as
a release-level reproducibility check, colour-coding each event as conserved
(same bracket in both releases), shifted, or single-release.

---

## The switching-gene model

Let `x_gs` be the TPM value of gene *g* in sample *s* of a tissue. Genes whose
raw TPM range within the tissue does not exceed **ε = 0.01** are discarded,
since constant profiles cannot be rescaled stably. The rest are min–max
normalized across the samples of the tissue:

```
z_gs = (x_gs − min_s x_gs) / (max_s x_gs − min_s x_gs)
```

Normalization is per gene and per tissue, **before** age-stratified
aggregation, so the profile describes the relative position of a gene within
its own observed range rather than an absolute expression level. Normalized
values are then averaged within each bracket, giving a six-value profile.

A threshold **τ ∈ [0,1]** turns the profile into an ordered binary vector,
`b_gi = 1(z̄_gi ≥ τ)`. Gene *g* is a candidate switching gene when that vector
contains **exactly one transition**. Two attributes are recorded: the
*switching bracket* (the first bracket in which the new state appears) and the
*direction* (low→high or high→low). Constant vectors, vectors with two or more
transitions, and genes with a missing bracket are excluded. Since the first
bracket has no predecessor, each tissue yields five ordered gene sets.

τ = 0.5 is the midpoint used for primary reporting, **not** a
threshold-invariant biological cutoff.

---

## Reproducing the results

The exact sample-level historical TPM matrix is approximately **4.01 GB**.
The interface atlas and the Ovary demo are included in the repository. The
complete inputs are prepared as nine separate release assets, whose publication
is still pending. Follow [the manuscript reproduction guide](reproducibility/README.md).

| Level | Required inputs | Scope |
|-------|-----------------|-------|
| **1** | Committed atlas | Interface and analyses of the bundled age averages |
| **2** | Bundled Ovary demo | Fresh observed counts and short Monte Carlo replay |
| **3** | Exact inputs identified by SHA-256 | Full manuscript experiments using `paper run` |

### Level 1 — from the committed atlas

The switching calls for all 50 age-complete tissues are committed, so the
published atlas can be checked immediately:

```bash
pip install -r requirements.txt && pip install -e .
pytest tests -q                                         # public core tests
python scripts/04_threshold_sensitivity.py --version v10
```

The last command recomputes the switching calls in the 13 primary tissues at
τ = 0.40, 0.45, 0.55, 0.60 and compares each against τ = 0.5, on gene identity
and on (gene, bracket) pairs. Expected output — these are the numbers reported
in the application note:

```
   tau   gene identity   gene + bracket
  0.40            0.15             0.01
  0.45            0.33             0.09
  0.55            0.30             0.08
  0.60            0.05             0.00
```

Counts vary markedly across thresholds, which is why τ = 0.5 is described as
the midpoint used for primary reporting rather than a threshold-invariant
cutoff. Together with the bootstrap, this cautions against treating individual
gene–bracket assignments as stable, and anchors the primary inference to the
tissue-level permutation analysis at τ = 0.5.

Because the normalized matrices are committed too, switching can be re-run at
any threshold without the raw data:

```bash
python scripts/02_switching.py --version v10 --complete-age-bins --threshold 0.6
python scripts/03_jaccard.py   --version v10 --complete-age-bins
```

### Level 2 — the validation chain on the demo tissue

The permutation test, SMOTE and the bootstrap all resample **donors**, so they
need sample-level data, which the atlas does not contain. `data/demo/v10/`
holds exactly that for one tissue: **Ovary** — 193 donors, one of the smallest
age-complete tissues, and one of the 13 primary tissues, so its result appears
in the paper.

Use the demo in place; never copy it over the complete input files:

```bash
python -m statistical_validation_cli.paper replay --demo
```

Create the dedicated Python 3.11 environment with
`requirements-reproduction.lock.txt` as described in the reproduction guide.
The replay checks the observed **3,300** calls, primary/bootstrap sequence
prefixes and both SMOTE policies against the historical experiments. It does
not recompute the complete 50-tissue BH family or certify the full analysis.
No demo cleanup or deletion of `data/parquet/v10/` is needed.

### Level 3 — exact historical inputs and complete experiments

Restore the nine original input parts and check their SHA-256, then run:

```bash
python -m statistical_validation_cli.paper run --jobs 2
```

This includes both full compound SMOTE policies. A new complete 50-tissue run
has not yet been performed. The input assets and public download links are
still pending publication.

The following general ingestion pipeline can build inputs from GTEx downloads.
Those reconstructed files have not been certified identical to the historical
inputs, so this route does not guarantee the exact published Monte Carlo results.

#### Get the data

From the GTEx portal Downloads page,
<https://gtexportal.org/home/downloads/adult-gtex>, take three open-access
files per release: the **gene-level TPM matrix** (GCT, gzipped), the **sample
attributes** and the **subject phenotypes**. Put them in
`data/external/{version}/` renamed as:

```
data/external/v10/gene_tpm.gct.gz
data/external/v10/SampleAttributesDS.txt
data/external/v10/SubjectPhenotypesDS.txt
```

Only gene-level TPM and the open-access annotations are used. Individual
genotypes, which require a dbGaP application (accession `phs000424`), are not
needed anywhere in STAMP.

#### Build the atlas

```bash
pip install -r requirements-pipeline.txt

python scripts/00_build_parquet_from_gtex.py --version v10
python scripts/01_normalize.py --version v10 --complete-age-bins
python scripts/02_switching.py --version v10 --complete-age-bins --threshold 0.5
python scripts/03_jaccard.py   --version v10 --complete-age-bins
```

Repeat with `--version v8` if you also want the v8 side of the Panel Explorer
concordance check. Stage 0 streams the GCT so peak memory stays around 100 MB
regardless of release size, but it does write several GB of Parquet.
The exact historical inputs used for the reported validation contain 19,616
samples, 946 donors and 59,033 genes; the advertised release totals in the
manuscript must be distinguished from this analyzed subset. These counts
describe the verified inputs rather than all samples advertised for the release.

#### Isolated environment for the validation

`imbalanced-learn` pins a narrow `scikit-learn` range, so keep it away from the
application environment:

```bash
python -m venv .venv-validation
source .venv-validation/bin/activate
pip install -r requirements.txt -r requirements-validation.txt
python -c "import imblearn, sklearn; print(imblearn.__version__, sklearn.__version__)"
```

Expected: `0.14.2 1.8.0`.

#### Output isolation

Every run writes under its own label, so an experiment can never overwrite the
published numbers:

```
output/v10/statistical_validation_cli/runs/final/
output/v10/statistical_validation_cli/runs/smoke/
output/v10/statistical_validation_cli/runs/pilot/
```

`final` is protected: it is rejected unless it covers all 50 age-complete
tissues, applies BH once to the complete family, and uses at least 9,999
permutations per tissue. `--quick` switches to `smoke` automatically.

#### The three blocks

They are deliberately kept separate, and answer different questions.

**1. Permutation test — the primary inference.** Age-bracket labels are
reassigned at random among the donors of each tissue, preserving the six
observed bracket sizes, and the whole STAMP analysis is repeated B = 9,999
times. The empirical right-tailed p-value uses the Monte Carlo correction
`p = (1 + #{T_b ≥ T_obs}) / (B + 1)`; BH is applied **once** to the family of
50 tissues.

```bash
python -m statistical_validation_cli.primary \
    --version v10 --run-label final --permutations 9999 --jobs 4
```

Expected: **13 tissues significant at q ≤ 0.05**, four of them at q ≤ 0.01;
five further tissues meet only the exploratory q ≤ 0.10 threshold, for 18 in
total. The 13 are the primary result; the other five are exploratory, and the
set of 18 must not be presented as 18 discoveries at FDR 5%. The test and the
q-values are at the **tissue** level — they do not validate any individual gene.

**2. SMOTE — exploratory sensitivity to age-bracket imbalance.**

```bash
python -m statistical_validation_cli.smote \
    --version v10 --run-label final \
    --seeds 20 --target-n 20 --permutations 0 --jobs 2

python -m statistical_validation_cli.smote \
    --version v10 --run-label smote-adaptive \
    --target-policy adaptive-second-smallest \
    --seeds 20 --target-n 20 --permutations 0 --jobs 4
```

Fixed-floor: brackets with 2–19 donors are raised to 20; singletons cannot be
interpolated and are left alone; brackets already at ≥ 20 are untouched.
Overall 39 of the 50 tissues are modified and 11 are unchanged. Synthetic
profiles are interpolations, not independent donors, so SMOTE is a secondary
sensitivity analysis and never enlarges the primary set of significant tissues.

**3. Stratified bootstrap — stability, not significance.**

```bash
python -m statistical_validation_cli.bootstrap \
    --version v10 --run-label final \
    --replicates 1000 --normalization fixed --jobs 4
```

Donors are resampled with replacement within each bracket, preserving the
original bracket sizes. Expected: among the 24,649 observed gene–tissue
associations in the 13 primary tissues, **2,782 (11.3 %)** keep the same
switching bracket and direction in at least 80 % of replicates. This is a
gene-level stability measure and answers a different question from the
tissue-level permutation test: it does not contradict tissue significance, but
it does mean individual gene–bracket–direction assignments should always be
reported with their bootstrap frequency.

`--normalization refit` recomputes the ε filter and min–max on each replicate
and measures the instability of the whole pipeline; run it under a separate
label.

**Report.**

```bash
python -m statistical_validation_cli.report --version v10 --run-label final
```

#### Complete manuscript reproduction

Use the exact checked inputs and pinned reproduction environment:

```bash
python -m statistical_validation_cli.paper run --jobs 2
```

The general `run_all` workflow uses descriptive SMOTE by default and is not
an equivalent command. Worker checkpoints require matching input, source and
environment hashes. The final success record lists reused and recomputed tissues.

See [`statistical_validation_cli/README.md`](statistical_validation_cli/README.md)
for the behaviour of the individual modules.

---

## Pipeline reference

| Stage | Script | Reads | Writes |
|-------|--------|-------|--------|
| 0 | `00_build_parquet_from_gtex.py` | `data/external/{v}/` | `data/parquet/{v}/` |
| 1 | `01_normalize.py` | `data/parquet/{v}/` | `output/{v}_complete/normalized/` |
| 2 | `02_switching.py` | `normalized/` | `output/{v}_complete/sets/` |
| 3 | `03_jaccard.py` | `sets/` | `output/{v}_complete/jaccard/` |
| 4 | `04_threshold_sensitivity.py` | `normalized/` | `output/{v}_complete/threshold_sensitivity/` |

Stages 2–4 need no raw data: they run off the committed `normalized/` matrices.

---

## Tests

```bash
pip install -e ".[dev]"
pytest
```

`validation/` contains unfinished comparison stubs, not an executed R validation.
The implemented checks are the unit tests, canonical regression and manuscript
reference/replay procedures described above.

---

## Citation

Kahveci O., Khatib T., Zicarelli M., Guzzi P.H., Veltri P., Mirabelli F.,
Kahveci T. *STAMP: an interactive web application for exploring spatiotemporal
switching-gene dynamics across human aging.* Bioinformatics Advances, 2026.

Software citation metadata is provided in `CITATION.cff`. The manuscript DOI
will be added when assigned.

## License

> **TO BE COMPLETED:** no license file is included yet. Without one the code is
> "all rights reserved" by default and nobody may legally reuse it, which
> conflicts with the availability statement in the application note. Agree a
> license with the co-authors (MIT and BSD-3-Clause are the usual choices for a
> tool like this) and add it as `LICENSE`.

## Contact

Pietro Hiram Guzzi — hguzzi@unicz.it
