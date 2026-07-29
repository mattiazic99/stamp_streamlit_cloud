"""Cloud entry adapter for the SMOTE Panel Explorer."""

from __future__ import annotations

import sys
from pathlib import Path


GUI_DIR = Path(__file__).resolve().parent.parent
if str(GUI_DIR) not in sys.path:
    sys.path.insert(0, str(GUI_DIR))

import smote_panel_data_cloud  # noqa: E402

# The shared UI imports ``smote_panel_data``. Route that import to the cloud
# adapter before loading the page, so no UI code is duplicated.
sys.modules["smote_panel_data"] = smote_panel_data_cloud

from analysis_modules.smote_panel_explorer import show  # noqa: E402,F401


__all__ = ["show"]
