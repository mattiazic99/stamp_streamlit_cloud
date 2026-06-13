import streamlit as st
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.express as px
import io
import base64
from utils.parsing import parse_stamp_file, extract_tissue_name
from components.downloads import create_download_button
from components.styling import apply_plot_style


# Canonical age-group order (used to keep x-axis consistent)
_AGE_ORDER = ["30–39", "40–49", "50–59", "60–69", "70–79"]


# ── helper: Plotly bar chart (single tissue) ────────────────────────────
def _interactive_bar_single(selected_ages, counts, tissue_name):
    """Interactive bar chart for a single tissue with hover details."""
    # Ensure age groups follow canonical order
    ordered = [(a, c) for a, c in zip(selected_ages, counts)]
    ordered.sort(key=lambda t: _AGE_ORDER.index(t[0]) if t[0] in _AGE_ORDER else 99)
    ages_sorted = [t[0] for t in ordered]
    counts_sorted = [t[1] for t in ordered]

    total = sum(counts_sorted)
    percentages = [(c / total * 100) if total > 0 else 0 for c in counts_sorted]

    n = len(ages_sorted)
    # Fixed bar width: keep bars slim
    bar_width = 0.28 if n >= 4 else (0.22 if n >= 2 else 0.15)

    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=ages_sorted,
        y=counts_sorted,
        width=[bar_width] * n,
        marker=dict(
            color=counts_sorted,
            colorscale="Blues",
            line=dict(color="rgb(8,48,107)", width=1.2),
        ),
        text=[f"<b>{c}</b>" for c in counts_sorted],
        textposition="outside",
        hovertemplate=(
            "<b>%{x}</b><br>"
            "Genes: <b>%{y}</b><br>"
            "Percentage: <b>%{customdata:.1f}%</b>"
            "<extra></extra>"
        ),
        customdata=percentages,
    ))

    y_max = max(counts_sorted) if counts_sorted else 1
    y_upper = y_max * 1.20 + 3

    fig.update_layout(
        title=dict(
            text=f"Distribution of Switching Genes by Age Group<br><sup>{tissue_name}</sup>",
            font=dict(size=15),
            x=0.5,
        ),
        xaxis_title="Age Group",
        yaxis_title="Number of Switching Genes",
        height=380,
        margin=dict(l=50, r=30, t=70, b=50),
        plot_bgcolor="white",
        yaxis=dict(gridcolor="rgba(0,0,0,0.08)", zeroline=False,
                   range=[0, y_upper]),
        xaxis=dict(
            gridcolor="rgba(0,0,0,0.08)",
            type="category",
            categoryorder="array",
            categoryarray=_AGE_ORDER,
        ),
    )
    return fig


# ── helper: Plotly grouped bar chart (two tissues) ──────────────────────
def _interactive_bar_comparison(selected_ages, counts1, counts2,
                                tissue1_name, tissue2_name, show_tissue):
    """Interactive comparison bar chart with optional tissue toggle."""
    # Ensure canonical order
    paired = list(zip(selected_ages, counts1, counts2))
    paired.sort(key=lambda t: _AGE_ORDER.index(t[0]) if t[0] in _AGE_ORDER else 99)
    ages_sorted = [t[0] for t in paired]
    c1_sorted = [t[1] for t in paired]
    c2_sorted = [t[2] for t in paired]

    n = len(ages_sorted)
    # Keep bars narrow: single-tissue width when only one shown, thinner when grouped
    is_both = show_tissue == "both"
    # Keep bars slim
    bar_w = (0.22 if n >= 4 else 0.16) if is_both else (0.26 if n >= 4 else 0.18)

    fig = go.Figure()

    if show_tissue in ("both", tissue1_name):
        total1 = sum(c1_sorted) or 1
        pct1 = [(c / total1 * 100) for c in c1_sorted]
        fig.add_trace(go.Bar(
            name=tissue1_name,
            x=ages_sorted,
            y=c1_sorted,
            width=[bar_w] * n,
            marker=dict(color="rgba(100,170,230,0.85)",
                        line=dict(color="rgb(8,48,107)", width=1)),
            text=[f"<b>{c}</b>" for c in c1_sorted],
            textposition="outside",
            hovertemplate=(
                f"<b>{tissue1_name}</b><br>"
                "Age: <b>%{x}</b><br>"
                "Genes: <b>%{y}</b><br>"
                "Tissue %: <b>%{customdata:.1f}%</b>"
                "<extra></extra>"
            ),
            customdata=pct1,
        ))

    if show_tissue in ("both", tissue2_name):
        total2 = sum(c2_sorted) or 1
        pct2 = [(c / total2 * 100) for c in c2_sorted]
        fig.add_trace(go.Bar(
            name=tissue2_name,
            x=ages_sorted,
            y=c2_sorted,
            width=[bar_w] * n,
            marker=dict(color="rgba(250,128,114,0.85)",
                        line=dict(color="darkred", width=1)),
            text=[f"<b>{c}</b>" for c in c2_sorted],
            textposition="outside",
            hovertemplate=(
                f"<b>{tissue2_name}</b><br>"
                "Age: <b>%{x}</b><br>"
                "Genes: <b>%{y}</b><br>"
                "Tissue %: <b>%{customdata:.1f}%</b>"
                "<extra></extra>"
            ),
            customdata=pct2,
        ))

    all_vals = c1_sorted + c2_sorted
    y_max = max(all_vals) if all_vals else 1
    y_upper = y_max * 1.20 + 3

    fig.update_layout(
        barmode="group",
        title=dict(
            text="Comparison of Switching Genes Between Tissues",
            font=dict(size=15),
            x=0.5,
        ),
        xaxis_title="Age Group",
        yaxis_title="Number of Switching Genes",
        height=380,
        margin=dict(l=50, r=30, t=70, b=50),
        plot_bgcolor="white",
        yaxis=dict(gridcolor="rgba(0,0,0,0.08)", zeroline=False,
                   range=[0, y_upper]),
        xaxis=dict(
            gridcolor="rgba(0,0,0,0.08)",
            type="category",
            categoryorder="array",
            categoryarray=_AGE_ORDER,
        ),
        legend=dict(orientation="h", yanchor="bottom", y=1.02,
                    xanchor="center", x=0.5),
    )
    return fig


# ── helper: Plotly heatmap (intra-tissue overlap) ───────────────────────
def _interactive_heatmap(matrix, age_groups, tissue_name):
    """Interactive heatmap with hover showing shared-gene count."""
    # Build hover text
    hover = []
    for i, ai in enumerate(age_groups):
        row = []
        for j, aj in enumerate(age_groups):
            val = int(matrix[i, j])
            if i == j:
                row.append(f"<b>{ai}</b><br>Total genes: {val}")
            else:
                row.append(f"<b>{ai}</b> ∩ <b>{aj}</b><br>Shared genes: {val}")
        hover.append(row)

    fig = go.Figure(data=go.Heatmap(
        z=matrix,
        x=age_groups,
        y=age_groups,
        colorscale="YlGnBu",
        text=matrix.astype(int).astype(str),
        texttemplate="<b>%{text}</b>",
        textfont=dict(size=13),
        hovertext=hover,
        hovertemplate="%{hovertext}<extra></extra>",
        colorbar=dict(title="Shared<br>Genes", thickness=14, len=0.75),
        xgap=2,
        ygap=2,
    ))

    fig.update_layout(
        title=dict(
            text=f"Gene Overlap Between Age Groups — {tissue_name}",
            font=dict(size=14),
            x=0.5,
        ),
        height=480,
        width=520,
        margin=dict(l=70, r=40, t=60, b=60),
        yaxis=dict(autorange="reversed", scaleanchor="x", scaleratio=1,
                   constrain="domain"),
        xaxis=dict(constrain="domain"),
        plot_bgcolor="white",
    )
    return fig


# ── Plotly chart config ──
def _plotly_cfg(filename="chart"):
    return {
        "toImageButtonOptions": {"format": "png", "scale": 2, "filename": filename.replace(".png", "")},
        "displayModeBar": True,
    }


def _download_plotly_as_png(plotly_fig, filename, button_key=None):
    """Render a working Download PNG button using Plotly.js from CDN."""
    import streamlit.components.v1 as _components
    safe_name = filename.replace(".png", "").replace("'", r"\'")
    fig_json = plotly_fig.to_json().replace("</", r"<\/")
    html = (
        '<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>'
        '<div id="hc" style="position:absolute;left:-9999px;width:1200px;height:600px;"></div>'
        '<div style="background:linear-gradient(90deg,#e8f4fd,#f0f7ff);border:1px solid #b8daff;'
        'border-radius:8px;padding:8px 14px;display:flex;align-items:center;gap:10px;">'
        '<button id="dlbtn" disabled style="background:#0d6efd;color:white;border:none;'
        'border-radius:6px;padding:6px 16px;cursor:pointer;font-size:13px;font-weight:600;'
        'white-space:nowrap;opacity:0.6;">Loading\u2026</button>'
        '<span style="color:#555;font-size:12px;">' + filename + '</span></div>'
        '<script>'
        'var _fig=' + fig_json + ';'
        'function _init(){if(typeof Plotly==="undefined"){setTimeout(_init,100);return;}'
        'var b=document.getElementById("dlbtn");b.textContent="\u2b07\ufe0f Download PNG";'
        'b.disabled=false;b.style.opacity="1";'
        'b.onclick=function(){'
        'b.textContent="Generating\u2026";b.disabled=true;'
        'Plotly.newPlot("hc",_fig.data,_fig.layout).then(function(){'
        'return Plotly.downloadImage("hc",{format:"png",scale:2,filename:"' + safe_name + '"});'
        '}).then(function(){Plotly.purge("hc");b.textContent="\u2b07\ufe0f Download PNG";b.disabled=false;});'
        '};}_init();'
        '</script>'
    )
    _components.html(html, height=50)


# ═══════════════════════════════════════════════════════════════════════
def show():
    """Upload and Single Tissue Analysis Page"""

    st.header("📤 Upload & Single Tissue Analysis")
    st.markdown("Upload STAMP gene switching files for individual tissue analysis.")

    # Global age groups
    age_groups = ["30–39", "40–49", "50–59", "60–69", "70–79"]

    # ── Data source toggle ─────────────────────────────────────────
    import sys as _sys
    from pathlib import Path as _Path
    _sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
    from data_loader import get_available_tissues, get_sets_for_tissue, sets_to_gui_format, complete_mode

    # Pre-loaded (GTEx) mode removed by design: users always upload their files.
    data_source = "📂 Upload files"

    complete = complete_mode()

    file1 = file2 = None  # defaults

    if data_source == "📦 Pre-loaded (GTEx)":
        version = st.session_state.get("gtex_version", "v10")
        tissues = get_available_tissues(version, complete)

        st.markdown(f"""
        <div class="analysis-section">
            <h3>🧪 Select Tissue(s) — GTEx {version}</h3>
        </div>
        """, unsafe_allow_html=True)

        col1, col2 = st.columns(2)
        with col1:
            tissue1 = st.selectbox("🧪 Primary tissue", tissues, key="ua_tissue1")
        with col2:
            tissue2 = st.selectbox(
                "🧪 Comparison tissue (optional)",
                ["— none —"] + tissues,
                key="ua_tissue2",
            )
        tissue2 = None if tissue2 == "— none —" else tissue2

        # Load and convert
        sets1 = get_sets_for_tissue(version, tissue1, complete)
        fasce1, counts1, df1 = sets_to_gui_format(sets1)

        if tissue2:
            sets2 = get_sets_for_tissue(version, tissue2, complete)
            fasce2, counts2, df2 = sets_to_gui_format(sets2)
            # Simulate file objects for tissue_name extraction
            class _FakeFile:
                def __init__(self, name): self.name = name
            file1 = _FakeFile(f"{tissue1}.txt")
            file2 = _FakeFile(f"{tissue2}.txt")
        else:
            fasce2, counts2, df2 = [], [], pd.DataFrame(columns=["Age", "Gene"])
            class _FakeFile:
                def __init__(self, name): self.name = name
            file1 = _FakeFile(f"{tissue1}.txt")

    else:
        # Original upload mode
        st.markdown("""
        <div class="analysis-section">
            <h3>📂 File Upload</h3>
        </div>
        """, unsafe_allow_html=True)

        col1, col2 = st.columns(2)
        with col1:
            file1 = st.file_uploader(
                "📂 Upload the first file",
                type=["txt"],
                key="file1",
                help="Upload a .txt file with gene switching data"
            )
        with col2:
            file2 = st.file_uploader(
                "📂 Upload the second file (optional)",
                type=["txt"],
                key="file2",
                help="Optional: Upload a second file for comparison"
            )

    # Age group selection
    st.markdown("### 🎯 Age Group Selection")
    selected_ages = st.multiselect(
        "Select age groups to display:",
        age_groups,
        default=age_groups,
        help="Choose which age groups to include in the analysis"
    )

    if not selected_ages:
        st.warning("⚠️ Please select at least one age group.")
        return

    # Parse files (only in upload mode — pre-loaded already has data)
    if data_source != "📦 Pre-loaded (GTEx)":
        fasce1, counts1, df1 = parse_stamp_file(file1, age_groups)
        fasce2, counts2, df2 = parse_stamp_file(file2, age_groups)

    # ── Show validation errors if any ──────────────────────────────
    _stop = False
    for _file, _df, _label in [(file1, df1, "File 1"), (file2, df2, "File 2")]:
        if _file is None:
            continue
        _val = _df.attrs.get('_validation', None)
        if _val and not _val.get('is_valid', True):
            st.error(
                f"❌ **{_label}** (`{_file.name}`) is not a valid STAMP file:\n\n"
                + "\n".join(f"- {e}" for e in _val.get('errors', []))
            )
            _stop = True
        elif _val and _val.get('warnings'):
            for w in _val['warnings']:
                st.warning(f"⚠️ {_label} (`{_file.name}`): {w}")
    if _stop:
        st.info(
            "ℹ️ **STAMP format**: each file must have exactly **5 lines** "
            "(one per age group), with **space-separated gene names** on each line."
        )
        return

    # ================================================================
    #  SINGLE TISSUE ANALYSIS
    # ================================================================
    if file1 and not file2:
        st.markdown("""
        <div class="analysis-section">
            <h2>📊 Single Tissue Analysis Results</h2>
        </div>
        """, unsafe_allow_html=True)

        tissue_name = extract_tissue_name(file1.name)

        # Filter data
        df1_filtered = df1[df1["Age"].isin(selected_ages)]
        counts1_filtered = [
            df1_filtered[df1_filtered["Age"] == age].shape[0]
            for age in age_groups if age in selected_ages
        ]

        # Display metrics
        col1, col2, col3 = st.columns(3)
        with col1:
            st.markdown(f"""
            <div class="metric-card">
                <h3>{len(df1)}</h3>
                <p>Total Genes</p>
            </div>
            """, unsafe_allow_html=True)

        with col2:
            st.markdown(f"""
            <div class="metric-card">
                <h3>{len(selected_ages)}</h3>
                <p>Age Groups</p>
            </div>
            """, unsafe_allow_html=True)

        with col3:
            avg_genes = np.mean(counts1_filtered) if counts1_filtered else 0
            st.markdown(f"""
            <div class="metric-card">
                <h3>{avg_genes:.1f}</h3>
                <p>Avg Genes/Group</p>
            </div>
            """, unsafe_allow_html=True)

        # ── Interactive bar chart ───────────────────────────────────
        st.markdown("### 📊 Gene Distribution by Age Group")

        bar_fig = _interactive_bar_single(selected_ages, counts1_filtered, tissue_name)
        st.plotly_chart(bar_fig, use_container_width=True, key="bar_single", config=_plotly_cfg())

        _download_plotly_as_png(
            bar_fig,
            f"gene_distribution_{tissue_name.replace(' ', '_').replace('-', '_')}.png",
            button_key="dl_bar_single"
        )

        # Download CSV
        st.markdown("""
        <div class="download-section">
            <h4>📥 Download Data</h4>
        </div>
        """, unsafe_allow_html=True)

        csv1 = df1.to_csv(index=False).encode('utf-8')
        st.download_button(
            label="⬇️ Download Complete Gene List (CSV)",
            data=csv1,
            file_name=f"switching_genes_{tissue_name.replace(' ', '_').replace('-', '_')}.csv",
            mime='text/csv'
        )

        # Gene lists by age group
        st.markdown("### 📋 Switching Genes by Age Group")

        for age in selected_ages:
            age_genes = df1[df1["Age"] == age]["Gene"].tolist()
            if age_genes:
                with st.expander(f"🔹 Age Group {age} – {len(age_genes)} genes"):
                    st.write(", ".join(age_genes))

                    age_df = pd.DataFrame(age_genes, columns=["Gene"])
                    csv_age = age_df.to_csv(index=False).encode('utf-8')
                    st.download_button(
                        label=f"⬇️ Download {age} genes (CSV)",
                        data=csv_age,
                        file_name=f"genes_{age.replace('–','_')}_{tissue_name.replace(' ', '_').replace('-', '_')}.csv",
                        mime='text/csv',
                        key=f"download_{age}"
                    )

        # ── Interactive heatmap ─────────────────────────────────────
        if not df1.empty:
            st.markdown("""
            <div class="analysis-section">
                <h2>Intra-Tissue Heatmap - Age Group Overlap</h2>
            </div>
            """, unsafe_allow_html=True)

            gene_sets = {
                age: set(df1[df1["Age"] == age]["Gene"])
                for age in age_groups
            }

            matrix = np.zeros((len(age_groups), len(age_groups)), dtype=np.int64)
            for i, age_i in enumerate(age_groups):
                for j, age_j in enumerate(age_groups):
                    matrix[i, j] = int(len(gene_sets[age_i] & gene_sets[age_j]))

            # Interactive Plotly heatmap — centred, compact
            hm_fig = _interactive_heatmap(matrix, age_groups, tissue_name)

            # Centre the heatmap with columns so it doesn't stretch full width
            _left, _centre, _right = st.columns([1, 3, 1])
            with _centre:
                st.plotly_chart(hm_fig, use_container_width=True, key="heatmap_single", config=_plotly_cfg())

            _download_plotly_as_png(
                hm_fig,
                f"age_overlap_heatmap_{tissue_name.replace(' ', '_').replace('-', '_')}.png",
                button_key="dl_heatmap_single"
            )

            # Download overlap matrix CSV
            df_matrix = pd.DataFrame(matrix, index=age_groups, columns=age_groups)
            csv_matrix = df_matrix.to_csv(index=True).encode('utf-8')
            st.download_button(
                label="⬇️ Download Overlap Matrix (CSV)",
                data=csv_matrix,
                file_name=f"overlap_matrix_{tissue_name.replace(' ', '_').replace('-', '_')}.csv",
                mime='text/csv'
            )

    # ================================================================
    #  TWO TISSUE COMPARISON
    # ================================================================
    elif file1 and file2:
        st.markdown("""
        <div class="analysis-section">
            <h2>📊 Comparison Between Two Tissues</h2>
        </div>
        """, unsafe_allow_html=True)

        tissue1_name = extract_tissue_name(file1.name)
        tissue2_name = extract_tissue_name(file2.name)

        # Filter data
        df1_filtered = df1[df1["Age"].isin(selected_ages)]
        df2_filtered = df2[df2["Age"].isin(selected_ages)]

        counts1_filtered = [
            df1_filtered[df1_filtered["Age"] == age].shape[0]
            for age in age_groups if age in selected_ages
        ]
        counts2_filtered = [
            df2_filtered[df2_filtered["Age"] == age].shape[0]
            for age in age_groups if age in selected_ages
        ]

        # Comparison metrics
        col1, col2, col3 = st.columns(3)
        with col1:
            st.markdown(f"""
            <div class="metric-card">
                <h3>{len(df1)}</h3>
                <p>{tissue1_name}</p>
            </div>
            """, unsafe_allow_html=True)

        with col2:
            st.markdown(f"""
            <div class="metric-card">
                <h3>{len(df2)}</h3>
                <p>{tissue2_name}</p>
            </div>
            """, unsafe_allow_html=True)

        with col3:
            total_unique = len(set(df1["Gene"]) | set(df2["Gene"]))
            st.markdown(f"""
            <div class="metric-card">
                <h3>{total_unique}</h3>
                <p>Unique Genes</p>
            </div>
            """, unsafe_allow_html=True)

        # ── Tissue toggle + interactive comparison chart ────────────
        st.markdown("### 📊 Gene Count Comparison")

        show_tissue = st.radio(
            "Show:",
            options=["both", tissue1_name, tissue2_name],
            index=0,
            horizontal=True,
            key="tissue_toggle",
        )

        comp_fig = _interactive_bar_comparison(
            selected_ages, counts1_filtered, counts2_filtered,
            tissue1_name, tissue2_name, show_tissue,
        )
        st.plotly_chart(comp_fig, use_container_width=True, key="bar_comparison", config=_plotly_cfg())

        # ── Download section: only current view ───────────────────────
        t1_safe = tissue1_name.replace(' ', '_').replace('-', '_')
        t2_safe = tissue2_name.replace(' ', '_').replace('-', '_')

        # Map toggle to filename suffix
        if show_tissue == "both":
            dl_suffix = f"both_{t1_safe}_vs_{t2_safe}"
        elif show_tissue == tissue1_name:
            dl_suffix = f"{t1_safe}_only"
        else:
            dl_suffix = f"{t2_safe}_only"

        _download_plotly_as_png(
            comp_fig,
            f"comparison_{dl_suffix}.png",
            button_key="dl_comp_current"
        )

        # Download both datasets
        st.markdown("""
        <div class="download-section">
            <h4>📥 Download Datasets</h4>
        </div>
        """, unsafe_allow_html=True)

        col1, col2 = st.columns(2)
        with col1:
            csv1 = df1.to_csv(index=False).encode('utf-8')
            st.download_button(
                label=f"⬇️ Download {tissue1_name} (CSV)",
                data=csv1,
                file_name=f"dataset1_{t1_safe}.csv",
                mime='text/csv'
            )

        with col2:
            csv2 = df2.to_csv(index=False).encode('utf-8')
            st.download_button(
                label=f"⬇️ Download {tissue2_name} (CSV)",
                data=csv2,
                file_name=f"dataset2_{t2_safe}.csv",
                mime='text/csv'
            )

    elif not file1:
        st.info("👆 Please upload at least one file to start the analysis.")

        st.markdown("### 📋 Expected File Format")
        st.code("""
GENE1 GENE2 GENE3 GENE4
GENE5 GENE6 GENE7
GENE8 GENE9 GENE10 GENE11 GENE12
GENE13 GENE14
GENE15 GENE16 GENE17 GENE18
        """)
        st.caption("Each line represents an age group (30-39, 40-49, 50-59, 60-69, 70-79)")