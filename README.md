# STAMP

[Open the application](https://stampgui.streamlit.app/)

STAMP explores gene switching across age groups using GTEx v10. Panel Explorer
also compares events with GTEx v8.

## Docker

Start Docker Desktop with Linux containers. From the repository root, run:

```bash
docker compose up --build --wait
```

Open [localhost:8501](http://localhost:8501). To stop the container:

```bash
docker compose down
```

Docker runs the interface with the bundled atlas. Generated downloads are held
in memory; download the ZIP to keep a local copy. **Start New Generation** clears
the current results. Statistical experiments run separately.

## Local installation

Use **Python 3.11**. From the repository root, run:

```bash
python -m pip install -r requirements-interface.lock.txt
python -m streamlit run gui/main.py --server.port=8502
```

Open [localhost:8502](http://localhost:8502). Press **Ctrl+C** to stop the app.
This installs the interface dependencies into the Python environment in use.

For Streamlit Community Cloud, select `streamlit_app.py` and **Python 3.11** in
**Advanced settings**. An existing app must be redeployed to change its Python
version: [deployment instructions](https://docs.streamlit.io/deploy/streamlit-community-cloud/manage-your-app/upgrade-python).

## Statistical reproducibility

Use a separate **Python 3.11** environment with the pinned numerical dependencies.
Run all commands from the repository root.

```bash
python -m venv .venv-reproduction
```

Activate it with `.\.venv-reproduction\Scripts\Activate.ps1` on Windows PowerShell,
or `source .venv-reproduction/bin/activate` on Linux/macOS. Then run:

```bash
python -m pip install -r requirements-reproduction.lock.txt
python -m statistical_validation_cli.paper check-reference
```

To check a fresh replay using the bundled Ovary data:

```bash
python -m statistical_validation_cli.paper replay --demo
```

This checks short sequence prefixes on one tissue, not the full manuscript analysis.

For the complete analysis, obtain the nine input parts listed in
`reproducibility/input_bundle.json` and place them in `input-parts/`.
**The complete input download is not yet published; it is separate from the code.**

```bash
python -m statistical_validation_cli.input_bundle --bundle-dir input-parts
python -m statistical_validation_cli.paper check-inputs
python -m statistical_validation_cli.paper run --jobs 2
```

The full run recomputes permutation tests, bootstrap stability, both SMOTE policies
and threshold sensitivity, then checks the results against the manuscript references.
Outputs are saved in `output/reproduction/`; `paper_verified.json` is written after
successful verification. A full run requires substantially more time and memory
than the demo.

Input hashes, numerical settings and verification details are documented in
[the reproduction protocol](reproducibility/README.md).
