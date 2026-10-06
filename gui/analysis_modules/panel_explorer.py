"""Map a gene panel onto switching events in GTEx v8 and v10.

Events are classified as conserved, shifted to another age bracket, or present
in only one release. Unresolved symbols and genes absent from the atlas are
reported separately: absence is not evidence that a gene does not change with age.

At the default threshold, events come from the bundled sets files. Other
thresholds recompute events from normalized matrices using the backend switching
rule. This is the only page that reads v8; the other pages use v10.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from data_loader import (  # noqa: E402
    COMPLETE_AGE_BINS,
    get_available_tissues,
    get_sets_for_tissue,
    get_normalized_tissue,
)

try:
    from stamp.config import SWITCHING_BRACKETS, DEFAULT_THRESHOLD  # noqa: E402
except Exception:  # pragma: no cover - defensive fallback
    SWITCHING_BRACKETS = ["30-39", "40-49", "50-59", "60-69", "70-79"]
    DEFAULT_THRESHOLD = 0.5

try:
    from stamp.switching import identify_switching_genes  # noqa: E402
except Exception:  # pragma: no cover - defensive fallback
    identify_switching_genes = None


# Defaults: the Alzheimer's-disease example panel
_DEFAULT_PANEL = [
    "APP", "PSEN1", "MAPT", "APOE", "TREM2", "BIN1",   # neurodegeneration
    "IL6", "TNF", "NFKB1", "STAT3",                     # neuroinflammation
    "ACTB", "GAPDH",                                    # housekeeping controls
]
_DEFAULT_TISSUES = [
    "Brain - Cortex", "Brain - Hippocampus", "Brain - Frontal Cortex (BA9)",
    "Brain - Amygdala", "Brain - Substantia nigra",
    "Whole Blood", "Spleen",
]

# Status codes used by the heatmap
_NONE, _V8_ONLY, _V10_ONLY, _SHIFTED, _CONSERVED = 0, 1, 2, 3, 4
_STATUS_LABEL = {
    _V8_ONLY: "v8 only", _V10_ONLY: "v10 only",
    _SHIFTED: "Shifted", _CONSERVED: "Conserved",
}
_STATUS_COLOR = {
    _NONE: "#eceff1", _V8_ONLY: "#e57373", _V10_ONLY: "#64b5f6",
    _SHIFTED: "#ffb74d", _CONSERVED: "#81c784",
}


def _plotly_cfg(fn="panel_explorer"):
    return {
        "toImageButtonOptions": {"format": "png", "scale": 2, "filename": fn},
        "displayModeBar": True,
    }


# Symbol ↔ Ensembl mapping (built once from gui/data/all_genes.txt)
@st.cache_data(ttl=3600)
def _load_symbol_map() -> tuple[dict[str, set[str]], dict[str, str]]:
    """Return (symbol_upper -> {unversioned ENSG}, unversioned ENSG -> symbol).

    The mapping file stores lines like ``ENSG00000142192.22 (APP)``.
    Versions are stripped so the same symbol matches across GTEx releases
    (v8 and v10 use different GENCODE versions).
    """
    path = Path(__file__).resolve().parent.parent / "data" / "all_genes.txt"
    sym2ens: dict[str, set[str]] = {}
    ens2sym: dict[str, str] = {}
    pat = re.compile(r"(ENSG\d+)\.\d+\s+\(([^)]+)\)")
    try:
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                m = pat.match(line.strip())
                if not m:
                    continue
                ens, sym = m.group(1), m.group(2)
                sym2ens.setdefault(sym.upper(), set()).add(ens)
                ens2sym.setdefault(ens, sym)
    except FileNotFoundError:
        pass
    return sym2ens, ens2sym


def _strip_version(gene_id: str) -> str:
    return gene_id.split(".")[0]


@st.cache_data(ttl=600)
def _events_for_version(
    version: str, tissues: tuple[str, ...], complete: bool, tau: float,
):
    """Return {(gene_ens_unversioned, tissue_display): bracket} for a version.

    At the default threshold the pre-computed ``*_sets.txt`` files are used
    (this exactly reproduces the published atlas). At any other threshold the
    switching events are recomputed on the fly from the normalized expression
    matrices using the same backend rule (`identify_switching_genes`), so the
    page stays consistent with the rest of the app.
    """
    use_precomputed = abs(tau - DEFAULT_THRESHOLD) < 1e-9 or identify_switching_genes is None
    out: dict[tuple[str, str], str] = {}
    for t in tissues:
        try:
            if use_precomputed:
                sets = get_sets_for_tissue(version, t, complete)
            else:
                df = get_normalized_tissue(version, t, complete)
                sets = identify_switching_genes(df, tau)
        except Exception:
            continue
        for bracket in SWITCHING_BRACKETS:
            for g in sets.get(bracket, []):
                out[(_strip_version(g), t)] = bracket
    return out


def _parse_panel(raw: str) -> list[str]:
    """Split a free-text panel (newline / comma / space separated) into symbols."""
    tokens = re.split(r"[\s,;]+", raw.strip())
    seen, panel = set(), []
    for tok in tokens:
        s = tok.strip().upper()
        if s and s not in seen:
            seen.add(s)
            panel.append(s)
    return panel


def show():
    st.header("🧩 Panel & Reproducibility Explorer")
    st.markdown(
        "Query the switching atlas with a **custom gene panel** and see its "
        "selectivity and **v8/v10 reproducibility** in one view. This is the "
        "only page that reads GTEx v8, and it does so purely as a "
        "release-level reproducibility check on the v10 results."
    )

    complete = COMPLETE_AGE_BINS
    if complete:
        st.info(
            "🧪 **Complete age-bins mode**: only tissues with samples in all "
            "six age brackets are available."
        )

    sym2ens, ens2sym = _load_symbol_map()
    if not sym2ens:
        st.error(
            "❌ Could not load the gene symbol map (`gui/data/all_genes.txt`). "
            "Panel genes can only be matched by Ensembl ID without it."
        )

    # Tissues available in either release
    try:
        tissues_v8 = set(get_available_tissues("v8", complete))
        tissues_v10 = set(get_available_tissues("v10", complete))
    except Exception as e:
        st.error(f"❌ Cannot load tissue lists: {e}")
        return
    all_tissues = sorted(tissues_v8 | tissues_v10)

    # Inputs
    c1, c2 = st.columns([1, 1])
    with c1:
        panel_raw = st.text_area(
            "🧬 Gene panel (one symbol per line, or comma/space separated):",
            value="\n".join(_DEFAULT_PANEL),
            height=200,
            help="Default = the Alzheimer's-disease example panel.",
        )
    with c2:
        default_tissues = [t for t in _DEFAULT_TISSUES if t in all_tissues]
        sel_tissues = st.multiselect(
            "🧪 Tissues to query:",
            options=all_tissues,
            default=default_tissues or all_tissues[:7],
            help="Default = AD-relevant CNS regions + peripheral immune tissues.",
        )

    # Threshold (tau)
    tau = st.slider(
        "🎚️ Switching threshold τ (binarisation cut-off)",
        min_value=0.05, max_value=0.95, value=float(DEFAULT_THRESHOLD), step=0.05,
        help="Expression is binarised as 1 if value ≥ τ. τ = 0.5 reproduces the "
             "published atlas (pre-computed sets); other values are recomputed "
             "on the fly from the normalized matrices.",
    )
    if abs(tau - DEFAULT_THRESHOLD) < 1e-9:
        st.caption(f"τ = {tau:.2f} — **default**: reproduces the published switching atlas.")
    elif identify_switching_genes is None:
        st.warning("Backend switching function unavailable — falling back to τ = 0.5 (atlas).")
        tau = float(DEFAULT_THRESHOLD)
    else:
        st.caption(f"τ = {tau:.2f} — events recomputed on the fly (sensitivity analysis).")

    panel = _parse_panel(panel_raw)
    if not panel:
        st.warning("Enter at least one gene symbol.")
        return
    if not sel_tissues:
        st.warning("Select at least one tissue.")
        return

    tissues_tuple = tuple(sel_tissues)
    ev_v8 = _events_for_version("v8", tissues_tuple, complete, tau)
    ev_v10 = _events_for_version("v10", tissues_tuple, complete, tau)

    # Resolve symbols
    unresolved = [g for g in panel if g not in sym2ens]
    resolved = [g for g in panel if g in sym2ens]

    # Build per-(gene, tissue) status + long-form records
    records = []                      # long table rows
    status = {}                       # (gene, tissue) -> code
    celltext = {}                     # (gene, tissue) -> annotation
    genes_with_event = set()

    for g in resolved:
        ens_set = sym2ens[g]
        for t in sel_tissues:
            br8 = next((ev_v8[(e, t)] for e in ens_set if (e, t) in ev_v8), None)
            br10 = next((ev_v10[(e, t)] for e in ens_set if (e, t) in ev_v10), None)
            if br8 is None and br10 is None:
                status[(g, t)] = _NONE
                celltext[(g, t)] = ""
                continue
            genes_with_event.add(g)
            if br8 is not None and br10 is not None:
                if br8 == br10:
                    code, txt = _CONSERVED, br10
                else:
                    code, txt = _SHIFTED, f"v8:{br8}<br>v10:{br10}"
            elif br10 is not None:
                code, txt = _V10_ONLY, f"v10:{br10}"
            else:
                code, txt = _V8_ONLY, f"v8:{br8}"
            status[(g, t)] = code
            celltext[(g, t)] = txt
            records.append({
                "Gene": g, "Tissue": t,
                "Bracket v8": br8 or "—", "Bracket v10": br10 or "—",
                "Status": _STATUS_LABEL[code],
            })

    present = sorted(genes_with_event)
    absent = [g for g in resolved if g not in genes_with_event]

    # Selectivity summary
    st.markdown('<div class="analysis-section"><h2>🎯 Selectivity</h2></div>',
                unsafe_allow_html=True)
    m1, m2, m3, m4 = st.columns(4)
    for col, val, lab in [
        (m1, len(panel), "Genes queried"),
        (m2, len(present), "Present in atlas"),
        (m3, len(absent), "Queried · absent"),
        (m4, len(unresolved), "Symbol not found"),
    ]:
        col.markdown(
            f'<div class="metric-card"><h3>{val}</h3><p>{lab}</p></div>',
            unsafe_allow_html=True,
        )

    if absent:
        st.markdown(
            "**Queried but not represented in the atlas** (under the applied "
            f"criteria — *not* evidence that they are age-invariant): "
            f"{', '.join(absent)}"
        )
    if unresolved:
        st.caption(f"⚠️ Symbol not found in annotation: {', '.join(unresolved)}")

    if not present:
        st.info("None of the resolved panel genes switch in the selected tissues.")
        return

    # Reproducibility heatmap (genes × tissues)
    st.markdown('<div class="analysis-section"><h2>🗺️ Event map & v8/v10 '
                'reproducibility</h2></div>', unsafe_allow_html=True)

    # Keep tissues with at least one event as columns (clearer figure),
    # but always keep the full resolved-gene list as rows so absences show.
    cols = [t for t in sel_tissues
            if any(status.get((g, t), _NONE) != _NONE for g in resolved)]
    if not cols:
        cols = list(sel_tissues)
    rows = resolved  # show every resolved gene, including all-absent ones

    z = [[status.get((g, t), _NONE) for t in cols] for g in rows]
    txt = [[celltext.get((g, t), "") for t in cols] for g in rows]

    colorscale = [
        [0.0, _STATUS_COLOR[_NONE]], [0.2, _STATUS_COLOR[_NONE]],
        [0.2, _STATUS_COLOR[_V8_ONLY]], [0.4, _STATUS_COLOR[_V8_ONLY]],
        [0.4, _STATUS_COLOR[_V10_ONLY]], [0.6, _STATUS_COLOR[_V10_ONLY]],
        [0.6, _STATUS_COLOR[_SHIFTED]], [0.8, _STATUS_COLOR[_SHIFTED]],
        [0.8, _STATUS_COLOR[_CONSERVED]], [1.0, _STATUS_COLOR[_CONSERVED]],
    ]
    heat = go.Figure(go.Heatmap(
        z=z, x=[c.replace("Brain - ", "") for c in cols], y=rows,
        text=txt, texttemplate="%{text}", textfont=dict(size=11),
        colorscale=colorscale, zmin=-0.5, zmax=4.5, showscale=False,
        xgap=3, ygap=3,
        hovertemplate="<b>%{y}</b> × <b>%{x}</b><br>%{text}<extra></extra>",
    ))
    heat.update_layout(
        height=max(320, 46 * len(rows) + 120),
        margin=dict(l=80, r=20, t=30, b=120),
        plot_bgcolor="white",
        xaxis=dict(side="bottom", tickangle=-35),
        yaxis=dict(autorange="reversed"),
    )
    st.plotly_chart(heat, use_container_width=True, key="pe_heat",
                    config=_plotly_cfg())

    # Discrete legend
    legend = " &nbsp; ".join(
        f'<span style="background:{_STATUS_COLOR[c]};padding:2px 10px;'
        f'border-radius:4px;">&nbsp;</span> {lab}'
        for c, lab in [
            (_CONSERVED, "Conserved (same bracket)"),
            (_SHIFTED, "Shifted (both, diff. bracket)"),
            (_V10_ONLY, "v10 only"),
            (_V8_ONLY, "v8 only"),
            (_NONE, "No event"),
        ]
    )
    st.markdown(legend, unsafe_allow_html=True)

    # Per-gene reproducibility class
    st.markdown('<div class="analysis-section"><h2>📋 Per-gene reproducibility'
                '</h2></div>', unsafe_allow_html=True)

    def _gene_class(g):
        codes = [status[(g, t)] for t in sel_tissues
                 if status.get((g, t), _NONE) != _NONE]
        if not codes:
            return "Absent", 0
        n_single = sum(c in (_V8_ONLY, _V10_ONLY) for c in codes)
        n_shift = sum(c == _SHIFTED for c in codes)
        if n_single:
            return "Context-dependent", len(codes)
        if n_shift:
            return "Mostly conserved", len(codes)
        return "Conserved", len(codes)

    cls_rows = []
    for g in resolved:
        cls, n = _gene_class(g)
        cls_rows.append({"Gene": g, "Reproducibility": cls, "n events": n})
    cls_df = pd.DataFrame(cls_rows)
    order = {"Conserved": 0, "Mostly conserved": 1,
             "Context-dependent": 2, "Absent": 3}
    cls_df = cls_df.sort_values(
        by=["Reproducibility", "Gene"],
        key=lambda s: s.map(order) if s.name == "Reproducibility" else s,
    ).reset_index(drop=True)
    st.dataframe(cls_df, use_container_width=True, hide_index=True)

    # Long event table + download
    if records:
        ev_df = pd.DataFrame(records).sort_values(
            ["Gene", "Tissue"]).reset_index(drop=True)
        with st.expander("📄 Full event table"):
            st.dataframe(ev_df, use_container_width=True, hide_index=True)
        st.download_button(
            "⬇️ Download panel events CSV",
            ev_df.to_csv(index=False).encode("utf-8"),
            "panel_switching_events.csv", "text/csv", key="pe_dl",
        )

    # Interpretation caveats
    st.info(
        "ℹ️ **How to read this page.** An *absent* gene is one not represented "
        "in the switching atlas under the applied criteria — this is **not** "
        "evidence that the gene is invariant with age. The age bracket marks "
        "when a persistent change of state is **first** detected in "
        "cross-sectional data; it does not denote disease onset, causality, or "
        "an irreversible transition. v8/v10 concordance is a robustness check, "
        "not an independent replication cohort."
    )
