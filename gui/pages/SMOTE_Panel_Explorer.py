"""Streamlit Cloud entry point for the GTEx v10 SMOTE Panel Explorer."""

from __future__ import annotations

import sys
from pathlib import Path

import streamlit as st


GUI_DIR = Path(__file__).resolve().parent.parent
if str(GUI_DIR) not in sys.path:
    sys.path.insert(0, str(GUI_DIR))

from analysis_modules.smote_panel_explorer_cloud import show  # noqa: E402


st.set_page_config(
    page_title="SMOTE Panel Explorer — GTEx v10",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.sidebar.markdown("## 🧪 SMOTE Panel Explorer")
st.sidebar.success("Dataset: **GTEx v10 only**")
st.sidebar.caption(
    "Fixed-floor and adaptive-second-smallest mapped switching-gene files."
)

show()
