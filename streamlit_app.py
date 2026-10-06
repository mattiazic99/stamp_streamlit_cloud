"""Launch gui/main.py from the repository root for Streamlit Community Cloud.

Add the project and gui directories to the import path before loading the app.
For local use, streamlit run gui/main.py launches the same interface.
"""
import runpy
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GUI = ROOT / "gui"

for path in (str(ROOT), str(GUI)):
    if path not in sys.path:
        sys.path.insert(0, path)

runpy.run_path(str(GUI / "main.py"), run_name="__main__")
