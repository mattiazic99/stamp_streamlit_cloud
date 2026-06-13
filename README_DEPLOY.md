# STAMP — Streamlit Cloud deployment

This folder is a self-contained, deploy-ready copy of the STAMP app.

## What's inside
- `stamp/` — analysis package
- `gui/` — Streamlit app (entry point: `gui/main.py`)
- `scripts/` — CLI to (re)generate results for any threshold
- `output/v8_complete/`, `output/v10_complete/` — pre-computed results the app reads
- `requirements.txt`, `pyproject.toml`

## Deploy on Streamlit Community Cloud
1. Create a new GitHub repo and push this folder.
2. On https://share.streamlit.io → New app → pick the repo.
3. Set **Main file path** to `gui/main.py`.
4. Deploy. `requirements.txt` is installed automatically.

## Notes
- The app runs in **complete age-bins** mode and reads `output/{version}_complete/`.
- Analysis pages are **upload-only** (the user uploads STAMP `.txt` files).
- "Pipeline Analysis" pages are hidden behind a checkbox in the sidebar.
- Raw GTEx TPM (`data/parquet`) is intentionally NOT included (not needed online).

## Regenerate results for a different threshold (local, no raw TPM needed)
The normalized matrices are included, so you can re-run switching at any tau:

    python scripts/02_switching.py --version v10 --complete-age-bins --threshold 0.6
    python scripts/03_jaccard.py   --version v10 --complete-age-bins

(Use `--version v8` for v8. Results are written back into `output/{version}_complete/`.)
