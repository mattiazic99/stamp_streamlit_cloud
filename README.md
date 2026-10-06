# STAMP

[Open the application](https://stampv2.streamlit.app/)

STAMP explores gene switching across age groups using GTEx v10. Panel Explorer
also compares events with GTEx v8.

## Docker

Docker can run STAMP on **Windows, macOS and Linux**. On Windows and macOS,
install and start [Docker Desktop](https://docs.docker.com/get-started/get-docker/).
On Linux, use Docker Engine with Compose or Docker Desktop.
The image runs Linux internally; Docker Desktop provides that environment on
Windows and macOS. On Windows, select Linux containers. Python and the interface
dependencies are included. Statistical experiments run separately, as described below.

### Ready-to-run image

Download [stamp-interface-linux-amd64.tar](https://github.com/mattiazic99/stamp_streamlit_cloud/releases/download/review-2026-10-06/stamp-interface-linux-amd64.tar)
and its [SHA-256 checksum](https://github.com/mattiazic99/stamp_streamlit_cloud/releases/download/review-2026-10-06/stamp-interface-linux-amd64.tar.sha256).
The image contains the interface and atlas for `linux/amd64`;
Apple Silicon Macs require amd64 emulation.
To verify the download on PowerShell, run
`Get-FileHash stamp-interface-linux-amd64.tar -Algorithm SHA256` and compare the
hash with the first value in the `.sha256` file (ignoring letter case).
Open a terminal in the folder containing the downloaded archive and run:

```bash
docker load --input stamp-interface-linux-amd64.tar
docker run --detach --name stamp-interface --publish 127.0.0.1:8501:8501 stamp-interface:latest
```

After startup, open [localhost:8501](http://localhost:8501).
To stop or restart this container:

```bash
docker stop stamp-interface
docker start stamp-interface
```

Before loading a newer release, stop the old container and remove it with
`docker rm stamp-interface`, then repeat the load and run commands.

### Build from source

Alternatively, download and extract this repository, or clone it. Open a terminal
in the folder containing `Dockerfile` and `docker-compose.yml`, then run:

```bash
docker compose up --build --wait
```

The first build requires an internet connection. When the command returns, open
[localhost:8501](http://localhost:8501). To stop and remove this container, run
`docker compose down` from the same folder. Use one Docker method at a time.

## Local installation

Use **Python 3.11**. Check that `python --version` reports 3.11 before installing.
From the repository root, run:

```bash
python -m pip install -r requirements-interface.lock.txt
python -m streamlit run gui/main.py --server.port=8502
```

Open [localhost:8502](http://localhost:8502). Press **Ctrl+C** to stop the app.
This installs the interface dependencies into the Python environment in use.

For Streamlit Community Cloud, select `streamlit_app.py` and **Python 3.11** in
**Advanced settings**. An existing app must be redeployed to change its Python
version: [deployment instructions](https://docs.streamlit.io/deploy/streamlit-community-cloud/manage-your-app/upgrade-python).

In all interface modes, generated files stay in session memory. Download the ZIP
to save them on your computer; **Start New Generation** clears the current results.

## Statistical reproducibility

Use a separate **Python 3.11** environment with the pinned numerical dependencies.
Run all commands from the repository root. Confirm that `python --version`
reports 3.11 before creating the environment.

```bash
python -m venv .venv-reproduction
```

Use the environment's executable directly; activation is not required.
The commands below use Windows PowerShell. On Linux/macOS, replace
`.\.venv-reproduction\Scripts\python.exe` with `.venv-reproduction/bin/python`.

```powershell
.\.venv-reproduction\Scripts\python.exe -m pip install -r requirements-reproduction.lock.txt
.\.venv-reproduction\Scripts\python.exe -m statistical_validation_cli.paper check-reference
```

`check-reference` checks the committed reference tables; it does not rerun experiments.

To check a fresh replay using the bundled Ovary data:

```powershell
.\.venv-reproduction\Scripts\python.exe -m statistical_validation_cli.paper replay --demo
```

This checks short sequence prefixes on one tissue, not the full manuscript analysis.

For the complete analysis, obtain the nine input parts listed in
`reproducibility/input_bundle.json` and place them in `input-parts/`.
**The complete input download is not yet published; it is separate from the code.**

```powershell
.\.venv-reproduction\Scripts\python.exe -m statistical_validation_cli.input_bundle --bundle-dir input-parts
.\.venv-reproduction\Scripts\python.exe -m statistical_validation_cli.paper check-inputs
.\.venv-reproduction\Scripts\python.exe -m statistical_validation_cli.paper run --jobs 2
```

The full run recomputes permutation tests, bootstrap stability, both SMOTE policies
and threshold sensitivity, then checks the results against the manuscript references.
Outputs are saved in `output/reproduction/`; `paper_verified.json` is written after
successful verification. A full run requires substantially more time and memory
than the demo.

**A fresh complete 50-tissue run with this reproduction command has not yet been
completed.** [GitHub CI](https://github.com/mattiazic99/stamp_streamlit_cloud/actions)
has passed the reference checks, fresh demo and interface container checks.

Input hashes, numerical settings and verification details are documented in
[the reproduction protocol](reproducibility/README.md).
