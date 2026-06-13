"""v8 vs v10 Comparison Page.

Robustness analysis: always loads BOTH GTEx versions regardless of the
global version selector. Shows scatter plots, Jaccard correlation, and
top similar tissue pairs for both versions side by side.
"""
import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from data_loader import (
    get_available_tissues_safe,
    get_jaccard_matrix,
    get_sets_for_tissue,
    safe_to_display,
    complete_mode,
)
from stamp.config import SWITCHING_BRACKETS


def _plotly_cfg(fn="chart"):
    return {"toImageButtonOptions": {"format": "png", "scale": 2, "filename": fn.replace(".png", "")}, "displayModeBar": True}


def show():
    st.header("⚖️ v8 vs v10 Comparison")
    st.markdown("Robustness analysis comparing switching gene results across GTEx releases.")

    complete = complete_mode()
    if complete:
        st.info(
            "🧪 **Complete age-bins mode**: comparing only tissues with samples "
            "in all six age brackets in both releases."
        )

    # Load tissues for both versions
    try:
        tissues_v8 = get_available_tissues_safe("v8", complete)
        tissues_v10 = get_available_tissues_safe("v10", complete)
    except Exception as e:
        st.error(f"❌ Cannot load tissue lists: {e}")
        return

    common_tissues = sorted(set(tissues_v8) & set(tissues_v10))
    if not common_tissues:
        st.error("❌ No common tissues found between v8 and v10.")
        return

    # Summary
    st.markdown('<div class="analysis-section"><h2>📊 Version Overview</h2></div>', unsafe_allow_html=True)
    c1, c2, c3 = st.columns(3)
    with c1:
        st.markdown(f'<div class="metric-card"><h3>{len(tissues_v8)}</h3><p>v8 Tissues</p></div>', unsafe_allow_html=True)
    with c2:
        st.markdown(f'<div class="metric-card"><h3>{len(tissues_v10)}</h3><p>v10 Tissues</p></div>', unsafe_allow_html=True)
    with c3:
        st.markdown(f'<div class="metric-card"><h3>{len(common_tissues)}</h3><p>Common Tissues</p></div>', unsafe_allow_html=True)

    # Compute n_switching per tissue per version
    @st.cache_data(ttl=600)
    def _compute_switching_counts(common_tissues, complete):
        rows = []
        for t in common_tissues:
            try:
                sets_v8 = get_sets_for_tissue.__wrapped__("v8", safe_to_display(t), complete)
                sets_v10 = get_sets_for_tissue.__wrapped__("v10", safe_to_display(t), complete)
                n_v8 = sum(len(v) for v in sets_v8.values())
                n_v10 = sum(len(v) for v in sets_v10.values())
                rows.append({"tissue": t, "display": safe_to_display(t), "n_v8": n_v8, "n_v10": n_v10})
            except Exception:
                continue
        return pd.DataFrame(rows)

    counts_df = _compute_switching_counts(common_tissues, complete)
    if counts_df.empty:
        st.error("❌ Could not load switching gene data.")
        return

    # Scatter: v8 vs v10
    st.markdown('<div class="analysis-section"><h2>📈 Switching Gene Count: v8 vs v10</h2></div>', unsafe_allow_html=True)
    scatter = go.Figure()
    scatter.add_trace(go.Scatter(
        x=counts_df["n_v8"], y=counts_df["n_v10"],
        mode="markers+text", text=counts_df["display"],
        textposition="top center", textfont=dict(size=8),
        marker=dict(size=8, color="#667eea", line=dict(color="white", width=1)),
        hovertemplate="<b>%{text}</b><br>v8: <b>%{x:,}</b><br>v10: <b>%{y:,}</b><extra></extra>",
    ))
    # Add y=x reference line
    max_val = max(counts_df["n_v8"].max(), counts_df["n_v10"].max()) * 1.1
    scatter.add_trace(go.Scatter(
        x=[0, max_val], y=[0, max_val], mode="lines",
        line=dict(color="rgba(0,0,0,0.2)", dash="dash", width=1),
        showlegend=False, hoverinfo="skip",
    ))
    scatter.update_layout(
        title=dict(text="n_switching: v8 vs v10 (one point per tissue)", font=dict(size=14), x=0.5),
        xaxis_title="n_switching (v8)", yaxis_title="n_switching (v10)",
        height=550, margin=dict(l=60, r=30, t=70, b=60),
        plot_bgcolor="white",
        xaxis=dict(gridcolor="rgba(0,0,0,0.08)"), yaxis=dict(gridcolor="rgba(0,0,0,0.08)"),
    )
    st.plotly_chart(scatter, use_container_width=True, key="vc_scatter", config=_plotly_cfg())

    # Pearson correlation of counts
    corr_counts = counts_df[["n_v8", "n_v10"]].corr().iloc[0, 1]
    st.markdown(f"**Pearson correlation (n_switching):** `{corr_counts:.4f}`")

    # Jaccard matrices
    st.markdown('<div class="analysis-section"><h2>🔗 Jaccard Matrix Correlation</h2></div>', unsafe_allow_html=True)

    try:
        jac_v8 = get_jaccard_matrix("v8", "life", complete)
        jac_v10 = get_jaccard_matrix("v10", "life", complete)
    except FileNotFoundError as e:
        st.error(f"❌ Jaccard matrix not found: {e}")
        return

    # Align matrices to common tissues
    common_in_jac = sorted(set(jac_v8.index) & set(jac_v10.index) & set(common_tissues))
    jac_v8_aligned = jac_v8.loc[common_in_jac, common_in_jac]
    jac_v10_aligned = jac_v10.loc[common_in_jac, common_in_jac]

    # Upper triangle correlation
    mask = np.triu_indices(len(common_in_jac), k=1)
    flat_v8 = jac_v8_aligned.values[mask]
    flat_v10 = jac_v10_aligned.values[mask]
    pearson_jac = np.corrcoef(flat_v8, flat_v10)[0, 1]

    c1, c2 = st.columns(2)
    with c1:
        st.markdown(f"""
        <div class="metric-card">
            <h3>{pearson_jac:.4f}</h3>
            <p>Pearson Correlation<br>(J<sub>life</sub> matrices)</p>
        </div>
        """, unsafe_allow_html=True)
    with c2:
        st.markdown(f"""
        <div class="metric-card">
            <h3>{len(common_in_jac)}</h3>
            <p>Tissues in Both<br>Jaccard Matrices</p>
        </div>
        """, unsafe_allow_html=True)

    # Top 10 similar pairs side by side
    st.markdown("### 🏆 Top 10 Most Similar Tissue Pairs")

    def _top_pairs(jac_df, n=10):
        rows = []
        tissues = jac_df.index.tolist()
        for i in range(len(tissues)):
            for j in range(i + 1, len(tissues)):
                rows.append({
                    "Tissue A": safe_to_display(tissues[i]),
                    "Tissue B": safe_to_display(tissues[j]),
                    "J_life": jac_df.iloc[i, j],
                })
        return pd.DataFrame(rows).sort_values("J_life", ascending=False).head(n).reset_index(drop=True)

    top_v8 = _top_pairs(jac_v8_aligned)
    top_v10 = _top_pairs(jac_v10_aligned)

    c1, c2 = st.columns(2)
    with c1:
        st.markdown("#### v8")
        st.dataframe(top_v8, use_container_width=True, hide_index=True,
                     column_config={"J_life": st.column_config.NumberColumn(format="%.4f")})
    with c2:
        st.markdown("#### v10")
        st.dataframe(top_v10, use_container_width=True, hide_index=True,
                     column_config={"J_life": st.column_config.NumberColumn(format="%.4f")})

    # Jaccard scatter
    st.markdown("### 📊 Jaccard Values: v8 vs v10 (pairwise)")
    jac_scatter = go.Figure()
    jac_scatter.add_trace(go.Scatter(
        x=flat_v8, y=flat_v10, mode="markers",
        marker=dict(size=3, color="#764ba2", opacity=0.5),
        hovertemplate="v8: <b>%{x:.3f}</b><br>v10: <b>%{y:.3f}</b><extra></extra>",
    ))
    jac_scatter.add_trace(go.Scatter(
        x=[0, 1], y=[0, 1], mode="lines",
        line=dict(color="rgba(0,0,0,0.2)", dash="dash", width=1),
        showlegend=False, hoverinfo="skip",
    ))
    jac_scatter.update_layout(
        title=dict(text=f"Pairwise J_life: v8 vs v10 (r = {pearson_jac:.4f})", font=dict(size=14), x=0.5),
        xaxis_title="J_life (v8)", yaxis_title="J_life (v10)",
        height=450, margin=dict(l=60, r=30, t=70, b=60),
        plot_bgcolor="white",
        xaxis=dict(gridcolor="rgba(0,0,0,0.08)"), yaxis=dict(gridcolor="rgba(0,0,0,0.08)"),
    )
    st.plotly_chart(jac_scatter, use_container_width=True, key="vc_jac_scatter", config=_plotly_cfg())

    # Download
    st.markdown('<div class="download-section"><h4>📥 Download Comparison Data</h4></div>', unsafe_allow_html=True)
    st.download_button("⬇️ Download Switching Counts CSV", counts_df.to_csv(index=False).encode("utf-8"),
                        "v8_vs_v10_switching_counts.csv", "text/csv", key="dl_vc_counts")
