import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
import plotly.figure_factory as ff
from utils.analysis import compute_jaccard_matrix, compute_linkage, compute_common_genes_matrix, compute_percent_overlap_matrix
from utils.parsing import parse_multiple_stamp_files, extract_tissue_name  # AGGIUNTO IMPORT
from components.downloads import create_csv_download, display_download_section


# ── Plotly chart config ──
def _plotly_cfg(filename="chart"):
    return {
        "toImageButtonOptions": {"format": "png", "scale": 2, "filename": filename.replace(".png", "")},
        "displayModeBar": True,
    }

def _download_plotly_as_png(plotly_fig, filename):
    """Render a working Download PNG button using Plotly.js from CDN."""
    import streamlit.components.v1 as _components
    safe_name = filename.replace(".png", "").replace("'", r"\'")
    fig_json = plotly_fig.to_json().replace("</", r"<\/")   # prevent HTML injection
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


# ── helper: Plotly heatmap (tissue × tissue, optional triangular mask) ──
def _interactive_heatmap_tissue(matrix, labels, title, colorbar_label,
                                 cmap="YlGnBu", mask_upper=True,
                                 fmt="float", key=None, center=None):
    n = matrix.shape[0]
    z_vals = matrix.copy().astype(float)
    text_vals, hover_vals = [], []
    for i in range(n):
        row_t, row_h = [], []
        for j in range(n):
            if mask_upper and j > i:
                z_vals[i, j] = None
                row_t.append("")
                row_h.append("")
            else:
                val = matrix[i, j]
                if fmt == "int":
                    row_t.append(f"{int(val)}")
                    row_h.append(f"<b>{labels[i]}</b> ∩ <b>{labels[j]}</b><br>{colorbar_label}: {int(val)}")
                else:
                    row_t.append(f"{val:.3f}")
                    row_h.append(f"<b>{labels[i]}</b> ∩ <b>{labels[j]}</b><br>{colorbar_label}: {val:.3f}")
        text_vals.append(row_t)
        hover_vals.append(row_h)

    hm_kwargs = dict(
        z=z_vals, x=labels, y=labels,
        colorscale=cmap,
        text=text_vals,
        texttemplate="<b>%{text}</b>",
        textfont=dict(size=11),
        hovertext=hover_vals,
        hovertemplate="%{hovertext}<extra></extra>",
        colorbar=dict(title=colorbar_label, thickness=14, len=0.75),
        xgap=2, ygap=2,
    )
    if center is not None:
        hm_kwargs['zmid'] = center

    fig = go.Figure(data=go.Heatmap(**hm_kwargs))

    side = min(520, max(380, n * 70))
    fig.update_layout(
        title=dict(text=title, font=dict(size=14), x=0.5),
        height=side, width=side + 40,
        margin=dict(l=80, r=40, t=60, b=80),
        yaxis=dict(autorange="reversed", scaleanchor="x", scaleratio=1, constrain="domain"),
        xaxis=dict(constrain="domain", tickangle=45),
        plot_bgcolor="white",
    )

    _left, _centre, _right = st.columns([1, 3, 1])
    with _centre:
        st.plotly_chart(fig, use_container_width=True, key=key, config=_plotly_cfg())
    _download_plotly_as_png(fig, f"{title.replace(' ', '_')[:60]}.png")


def _precomputed_jaccard_submatrix(version, tissues_display, metric, complete=False):
    """Canonical Jaccard sub-matrix from output/{version}[_complete]/jaccard/.

    Returns a numpy matrix for the given tissues (display names) in the given
    order, or None if any tissue is missing from the matrix or on any error
    (caller then falls back to the local computation).
    """
    try:
        from data_loader import get_jaccard_matrix, display_to_safe
        m = get_jaccard_matrix(version, metric, complete)
        safe = [display_to_safe(t) for t in tissues_display]
        if not all(s in m.index and s in m.columns for s in safe):
            return None
        return m.loc[safe, safe].to_numpy(dtype=float)
    except Exception:
        return None


def show():
    """Multi-Tissue Analysis Page"""
    
    st.header("🧬 Multi-Tissue Analysis")
    st.markdown("Comprehensive similarity analysis and hierarchical clustering across multiple tissues.")
    
    age_groups = ["30–39", "40–49", "50–59", "60–69", "70–79"]

    # Track pre-computed (GTEx) mode so we can read the canonical Jaccard
    # matrices instead of recomputing them client-side.
    preloaded = False
    version = None

    # ── Data source toggle ─────────────────────────────────────────
    import sys as _sys
    from pathlib import Path as _Path
    _sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
    from data_loader import get_available_tissues, get_sets_for_tissue, sets_to_gui_format, sets_to_gene_sets, complete_mode

    # Pre-loaded (GTEx) mode removed by design: users always upload their files.
    data_source = "📂 Upload files"

    complete = complete_mode()

    if data_source == "📦 Pre-loaded (GTEx)":
        preloaded = True
        version = st.session_state.get("gtex_version", "v10")
        all_tissues = get_available_tissues(version, complete)

        st.markdown(f"""
        <div class="analysis-section">
            <h3>🧪 Select Tissues — GTEx {version}</h3>
            <p>Select at least 3 tissues for comprehensive multi-tissue analysis</p>
        </div>
        """, unsafe_allow_html=True)

        selected_tissues = st.multiselect(
            "Select tissues:",
            all_tissues,
            default=all_tissues[:5],
            key="mt_tissues",
        )

        if not selected_tissues or len(selected_tissues) < 2:
            st.info("👆 Please select at least 2 tissues for multi-tissue analysis.")
            return

        # Build the same data dict as parse_multiple_stamp_files
        data = {}
        for tissue in selected_tissues:
            sets_dict = get_sets_for_tissue(version, tissue, complete)
            _, counts, df = sets_to_gui_format(sets_dict)
            gene_sets_list = sets_to_gene_sets(sets_dict)
            data[tissue] = {
                'gene_sets': gene_sets_list,
                'dataframe': df,
                'counts': counts,
                'total_genes': sum(len(s) for s in gene_sets_list),
                'original_filename': f"{tissue}.txt",
                'clean_name': tissue,
            }
        summary = {
            'total_files': len(selected_tissues),
            'successful_parses': len(selected_tissues),
            'failed_parses': 0,
            'total_unique_genes': len(set().union(*(set().union(*d['gene_sets']) for d in data.values()))),
            'tissue_names': list(data.keys()),
            'rejected_files': [],
        }
        st.success(f"✅ {len(selected_tissues)} tissues loaded from GTEx {version}!")

    else:
        # Original upload mode
        st.markdown("""
        <div class="analysis-section">
            <h3>📂 Upload Multiple Tissue Files</h3>
            <p>Upload at least 3 tissue files for comprehensive multi-tissue analysis</p>
        </div>
        """, unsafe_allow_html=True)

        uploaded_files = st.file_uploader(
            "📂 Upload STAMP .txt files", 
            type=["txt"],
            accept_multiple_files=True, 
            key="multi_tissue_files",
            help="Upload multiple tissue gene switching files for comparison"
        )

        if not uploaded_files or len(uploaded_files) < 2:
            st.info("👆 Please upload at least 2 tissue files for multi-tissue analysis.")
            
            st.markdown("### 📋 Analysis Features")
            col1, col2 = st.columns(2)
            with col1:
                st.markdown("""
                **🔍 Similarity Analysis:**
                - Jaccard similarity matrices
                - Age-averaged comparisons
                - Lifetime gene overlap
                - Correlation heatmaps
                """)
            with col2:
                st.markdown("""
                **🌳 Clustering Analysis:**
                - Hierarchical clustering
                - Dendrogram visualization
                - Distance matrices
                - Tissue grouping
                """)
            return

        st.success(f"✅ {len(uploaded_files)} tissue files loaded successfully!")

        # Parse all files
        result = parse_multiple_stamp_files(uploaded_files, age_groups)
        # ── Show rejected files (STAMP validation errors) ──────────────
        _rejected = result['summary'].get('rejected_files', [])
        if _rejected:
            for _rej in _rejected:
                st.error(
                    f"❌ **{_rej['filename']}** is not a valid STAMP file:\n\n"
                    + "\n".join(f"- {e}" for e in _rej.get('errors', []))
                )
            st.info(
                "ℹ️ **STAMP format**: each file must have exactly **5 lines** "
                "(one per age group), with **space-separated gene names** on each line."
            )

        data = result['data']
        summary = result['summary']
    
    tissues = list(data.keys())
    
    # Analysis options
    st.markdown("### ⚙️ Analysis Options")
    col1, col2 = st.columns(2)
    
    with col1:
        show_dendrograms = st.checkbox("🌳 Show dendrograms", value=True)
        similarity_metric = st.selectbox("📊 Similarity Metric", 
                                       ["Jaccard Index", "Overlap Coefficient"])
    
    with col2:
        clustering_method = st.selectbox("🔗 Clustering Method", 
                                       ["average", "complete", "single", "ward"])
        color_palette = st.selectbox("🎨 Color Palette", 
                                   ["rdbu", "viridis", "rdylbu", "plasma"])
    
    # === TISSUE OVERVIEW ===
    st.markdown("""
    <div class="analysis-section">
        <h2>📊 Tissue Overview</h2>
    </div>
    """, unsafe_allow_html=True)
    
    tissue_stats = []
    for tissue in tissues:
        gene_sets = data[tissue]['gene_sets']
        total_genes = len(set.union(*gene_sets)) if gene_sets else 0
        avg_genes_per_age = np.mean([len(age_set) for age_set in gene_sets]) if gene_sets else 0
        max_genes_age = max([len(age_set) for age_set in gene_sets]) if gene_sets else 0
        min_genes_age = min([len(age_set) for age_set in gene_sets]) if gene_sets else 0
        
        tissue_stats.append({
            'Tissue': tissue,
            'Total Unique Genes': total_genes,
            'Avg Genes/Age': f"{avg_genes_per_age:.1f}",
            'Max Genes (Age)': max_genes_age,
            'Min Genes (Age)': min_genes_age
        })
    
    df_stats = pd.DataFrame(tissue_stats)
    st.dataframe(df_stats, use_container_width=True)
    
    data_for_analysis = {tissue: data[tissue]['gene_sets'] for tissue in tissues}
    
    # Bar chart of total genes per tissue (Plotly)
    total_genes_list = [int(row['Total Unique Genes']) for row in tissue_stats]
    n_t = len(tissues)
    bar_w = 0.28 if n_t >= 4 else (0.22 if n_t >= 2 else 0.15)
    total_sum = sum(total_genes_list) or 1
    pcts = [(c / total_sum * 100) for c in total_genes_list]
    colors_list = px.colors.qualitative.Set3[:n_t]

    bar_fig = go.Figure()
    bar_fig.add_trace(go.Bar(
        x=tissues, y=total_genes_list,
        width=[bar_w] * n_t,
        marker=dict(color=colors_list, line=dict(color="black", width=1.2)),
        text=[f"<b>{c}</b>" for c in total_genes_list],
        textposition="outside",
        hovertemplate=(
            "<b>%{x}</b><br>"
            "Genes: <b>%{y}</b><br>"
            "Percentage: <b>%{customdata:.1f}%</b>"
            "<extra></extra>"
        ),
        customdata=pcts,
    ))
    y_max = max(total_genes_list) if total_genes_list else 1
    bar_fig.update_layout(
        title=dict(text="Total Unique Genes per Tissue", font=dict(size=15), x=0.5),
        xaxis_title="Tissue", yaxis_title="Number of Unique Genes",
        height=400, margin=dict(l=50, r=30, t=70, b=80), plot_bgcolor="white",
        yaxis=dict(gridcolor="rgba(0,0,0,0.08)", zeroline=False, range=[0, y_max * 1.2 + 3]),
        xaxis=dict(tickangle=45),
    )
    st.plotly_chart(bar_fig, use_container_width=True, key="mt_bar_overview", config=_plotly_cfg())
    _download_plotly_as_png(bar_fig, "tissue_overview.png")
    
    # === SIMILARITY MATRICES ===
    st.markdown("""
    <div class="analysis-section">
        <h2>📊 Similarity Analysis</h2>
    </div>
    """, unsafe_allow_html=True)
    
    # Jaccard similarity matrices: in pre-loaded mode read the canonical
    # pre-computed matrices (output/{version}/jaccard/jaccard_{age,life}.csv),
    # the same ones used by the paper and the v8-vs-v10 page. Fall back to the
    # local computation in upload mode (or if a tissue is missing in the matrix).
    tissues_sorted = list(data_for_analysis.keys())
    matrix_age = matrix_life = None
    jaccard_precomputed = False
    if preloaded:
        _ma = _precomputed_jaccard_submatrix(version, tissues_sorted, "age", complete)
        _ml = _precomputed_jaccard_submatrix(version, tissues_sorted, "life", complete)
        if _ma is not None and _ml is not None:
            matrix_age, matrix_life, jaccard_precomputed = _ma, _ml, True
    if not jaccard_precomputed:
        matrix_age, tissues_sorted = compute_jaccard_matrix(data_for_analysis, mode="age")
        matrix_life, _ = compute_jaccard_matrix(data_for_analysis, mode="life")

    if jaccard_precomputed:
        st.caption(
            f"Jaccard matrices read from pre-computed `output/{version}/jaccard/` "
            "— identical to the v8-vs-v10 page and the paper.",
        )

    # Age-averaged similarity heatmap (Plotly)
    st.markdown("### 📊 Average Similarity Across Age Groups")
    _interactive_heatmap_tissue(
        matrix_age, tissues_sorted,
        "Tissue Similarity (Averaged Over Age Groups)",
        "Jaccard Similarity", cmap=color_palette, mask_upper=True,
        fmt="float", key="mt_hm_age"
    )
    
    # Lifetime similarity heatmap (Plotly)
    st.markdown("### 📊 Lifetime Similarity Matrix")
    _interactive_heatmap_tissue(
        matrix_life, tissues_sorted,
        "Tissue Similarity (Whole Lifetime)",
        "Jaccard Similarity", cmap="YlGnBu", mask_upper=True,
        fmt="float", key="mt_hm_life"
    )
    
    # === HIERARCHICAL CLUSTERING (Plotly dendrograms) ===
    if show_dendrograms:
        st.markdown("""
        <div class="analysis-section">
            <h2>🌳 Hierarchical Clustering</h2>
        </div>
        """, unsafe_allow_html=True)
        
        linkage_age = compute_linkage(matrix_age, method=clustering_method)
        linkage_life = compute_linkage(matrix_life, method=clustering_method)
        
        if linkage_age is not None:
            st.markdown("### 🌿 Dendrogram - Age-Based Similarity")
            dist_age = 1 - matrix_age
            np.fill_diagonal(dist_age, 0)
            dist_age = np.clip(dist_age, 0, None)
            dendro_age = ff.create_dendrogram(
                dist_age, labels=tissues_sorted,
                linkagefun=lambda x: linkage_age,
            )
            dendro_age.update_layout(
                title=dict(text=f"Hierarchical Clustering (Age-Based, {clustering_method.title()} Linkage)",
                           font=dict(size=14), x=0.5),
                xaxis_title='Tissue', yaxis_title='Distance (1 − Similarity)',
                height=380, margin=dict(l=50, r=30, t=60, b=80),
                plot_bgcolor='white',
                yaxis=dict(gridcolor='rgba(0,0,0,0.08)'),
            )
            st.plotly_chart(dendro_age, use_container_width=True, key='mt_dendro_age', config=_plotly_cfg())
            _download_plotly_as_png(dendro_age, 'dendrogram_age_based.png')
        
        if linkage_life is not None:
            st.markdown("### 🌿 Dendrogram - Lifetime-Based Similarity")
            dist_life = 1 - matrix_life
            np.fill_diagonal(dist_life, 0)
            dist_life = np.clip(dist_life, 0, None)
            dendro_life = ff.create_dendrogram(
                dist_life, labels=tissues_sorted,
                linkagefun=lambda x: linkage_life,
            )
            dendro_life.update_layout(
                title=dict(text=f"Hierarchical Clustering (Lifetime-Based, {clustering_method.title()} Linkage)",
                           font=dict(size=14), x=0.5),
                xaxis_title='Tissue', yaxis_title='Distance (1 − Similarity)',
                height=380, margin=dict(l=50, r=30, t=60, b=80),
                plot_bgcolor='white',
                yaxis=dict(gridcolor='rgba(0,0,0,0.08)'),
            )
            st.plotly_chart(dendro_life, use_container_width=True, key='mt_dendro_life', config=_plotly_cfg())
            _download_plotly_as_png(dendro_life, 'dendrogram_lifetime_based.png')
    
    # === SHARED GENES ANALYSIS ===
    st.markdown("""
    <div class="analysis-section">
        <h2>🤝 Shared Genes Analysis</h2>
    </div>
    """, unsafe_allow_html=True)
    
    st.markdown("### 🔗 Absolute Shared Gene Counts")
    matrix_common, _ = compute_common_genes_matrix(data_for_analysis)
    if matrix_common.size > 0:
        _interactive_heatmap_tissue(
            matrix_common.astype(np.int64), tissues_sorted,
            "Shared Genes Between Tissues (Absolute Counts)",
            "Shared Gene Count", cmap="Purples", mask_upper=True,
            fmt="int", key="mt_hm_common"
        )
    
    st.markdown("### 📈 Percentage Overlap Analysis")
    matrix_percent, _ = compute_percent_overlap_matrix(data_for_analysis)
    if matrix_percent.size > 0:
        _interactive_heatmap_tissue(
            matrix_percent, tissues_sorted,
            "Gene Overlap Percentage Between Tissues",
            "Overlap %", cmap="teal", mask_upper=True,
            fmt="float", key="mt_hm_pct"
        )
    
    # === TISSUE RANKING ===
    st.markdown("""
    <div class="analysis-section">
        <h2>🏆 Tissue Similarity Ranking</h2>
    </div>
    """, unsafe_allow_html=True)
    
    if len(tissues) >= 3:
        ref_tissue = st.selectbox("🔍 Choose reference tissue for ranking:", 
                                 tissues_sorted, key="ranking_ref_multi")
        
        st.markdown(f"### 📋 Similarity Ranking relative to **{ref_tissue}**")
        
        ref_idx = tissues_sorted.index(ref_tissue)
        similarities = []
        for i, tissue in enumerate(tissues_sorted):
            if tissue != ref_tissue:
                similarity = matrix_life[ref_idx, i] * 100
                similarities.append((tissue, similarity))
        similarities.sort(key=lambda x: x[1], reverse=True)
        
        ranking_df = pd.DataFrame(similarities, columns=['Tissue', 'Similarity (%)'])
        ranking_df['Rank'] = range(1, len(ranking_df) + 1)
        ranking_df = ranking_df[['Rank', 'Tissue', 'Similarity (%)']]
        st.dataframe(ranking_df, use_container_width=True)
        
        # Ranking bar chart (Plotly horizontal bar)
        tissues_rank = [s[0] for s in similarities]
        sims_rank = [s[1] for s in similarities]

        rank_fig = go.Figure()
        rank_fig.add_trace(go.Bar(
            y=tissues_rank, x=sims_rank, orientation="h",
            marker=dict(color=sims_rank, colorscale="Viridis",
                        line=dict(color="black", width=1)),
            text=[f"<b>{v:.1f}%</b>" for v in sims_rank],
            textposition="outside",
            hovertemplate="<b>%{y}</b><br>Similarity: <b>%{x:.1f}%</b><extra></extra>",
        ))
        rank_fig.update_layout(
            title=dict(text=f"Tissue Similarity Ranking (Reference: {ref_tissue})",
                       font=dict(size=14), x=0.5),
            xaxis_title="Similarity (%)",
            height=max(300, len(similarities) * 40 + 100),
            margin=dict(l=120, r=60, t=60, b=50), plot_bgcolor="white",
            xaxis=dict(gridcolor="rgba(0,0,0,0.08)"),
            yaxis=dict(autorange="reversed"),
        )
        st.plotly_chart(rank_fig, use_container_width=True, key="mt_bar_ranking", config=_plotly_cfg())
        _download_plotly_as_png(rank_fig, f"ranking_{ref_tissue}.png")
        create_csv_download(ranking_df, f"similarity_ranking_{ref_tissue}.csv", 
                           "⬇️ Download Ranking CSV")
    
    # === DOWNLOAD SECTION ===
    display_download_section("📥 Download Analysis Results")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.markdown("**📊 Similarity Matrices**")
        df_age_sim = pd.DataFrame(matrix_age, index=tissues_sorted, columns=tissues_sorted)
        create_csv_download(df_age_sim, "similarity_matrix_age_averaged.csv", 
                           "⬇️ Age-Averaged Matrix CSV")
        df_life_sim = pd.DataFrame(matrix_life, index=tissues_sorted, columns=tissues_sorted)
        create_csv_download(df_life_sim, "similarity_matrix_lifetime.csv", 
                           "⬇️ Lifetime Matrix CSV")
    
    with col2:
        st.markdown("**🤝 Shared Gene Matrices**")
        if matrix_common.size > 0:
            df_common = pd.DataFrame(matrix_common, index=tissues_sorted, columns=tissues_sorted)
            create_csv_download(df_common, "shared_genes_matrix.csv", 
                               "⬇️ Shared Genes Matrix CSV")
        if matrix_percent.size > 0:
            df_percent = pd.DataFrame(matrix_percent, index=tissues_sorted, columns=tissues_sorted)
            create_csv_download(df_percent, "overlap_percentage_matrix.csv", 
                               "⬇️ Overlap Percentage CSV")
    
    with col3:
        st.markdown("**📈 Summary Statistics**")
        create_csv_download(df_stats, "tissue_statistics_summary.csv", 
                           "⬇️ Tissue Statistics CSV")
        
        summary_data = {
            'Analysis Type': ['Multi-Tissue Analysis'],
            'Number of Tissues': [len(tissues)],
            'Tissues': [', '.join(tissues_sorted)],
            'Similarity Metric': [similarity_metric],
            'Clustering Method': [clustering_method],
            'Total Unique Genes': [summary.get('total_unique_genes', 0)],
            'Average Similarity': [f"{np.mean(matrix_life[np.triu_indices_from(matrix_life, k=1)]):.3f}"],
            'Max Similarity': [f"{np.max(matrix_life[np.triu_indices_from(matrix_life, k=1)]):.3f}"],
            'Min Similarity': [f"{np.min(matrix_life[np.triu_indices_from(matrix_life, k=1)]):.3f}"]
        }
        summary_df = pd.DataFrame(summary_data)
        create_csv_download(summary_df, "analysis_summary.csv", 
                           "⬇️ Analysis Summary CSV")
    
    # === ADVANCED ANALYSIS ===
    if st.checkbox("🔬 Show Advanced Analysis", value=False):
        st.markdown("""
        <div class="analysis-section">
            <h2>🔬 Advanced Multi-Tissue Analysis</h2>
        </div>
        """, unsafe_allow_html=True)
        
        # Correlation heatmap (Plotly)
        st.markdown("### 📊 Gene Count Correlations Across Age Groups")
        
        count_matrix = np.zeros((len(tissues), len(age_groups)))
        for i, tissue in enumerate(tissues):
            gene_sets = data[tissue]['gene_sets']
            for j, age_idx in enumerate(range(5)):
                if age_idx < len(gene_sets):
                    count_matrix[i, j] = len(gene_sets[age_idx])
        
        corr_matrix = np.corrcoef(count_matrix)
        
        _interactive_heatmap_tissue(
            corr_matrix, tissues,
            "Gene Count Correlation Between Tissues",
            "Correlation", cmap="RdBu_r", mask_upper=True,
            fmt="float", key="mt_hm_corr", center=0
        )
        
        # Age group diversity analysis (Plotly bar)
        st.markdown("### 📅 Age Group Diversity Analysis")
        
        diversities = []
        for tissue in tissues:
            gene_sets = data[tissue]['gene_sets']
            counts = [len(gene_sets[i]) for i in range(len(gene_sets))]
            total = sum(counts)
            if total > 0:
                proportions = [c/total for c in counts if c > 0]
                shannon = -sum(p * np.log(p) for p in proportions)
                diversities.append(shannon)
            else:
                diversities.append(0)
        
        colors_div = px.colors.sequential.Plasma[:n_t]

        div_fig = go.Figure()
        div_fig.add_trace(go.Bar(
            x=tissues, y=diversities,
            marker=dict(color=diversities, colorscale="Plasma",
                        line=dict(color="black", width=1)),
            text=[f"<b>{d:.3f}</b>" for d in diversities],
            textposition="outside",
            hovertemplate="<b>%{x}</b><br>Shannon Diversity: <b>%{y:.3f}</b><extra></extra>",
        ))
        y_max_d = max(diversities) if diversities else 1
        div_fig.update_layout(
            title=dict(text="Age Group Diversity by Tissue (Shannon Index)",
                       font=dict(size=14), x=0.5),
            xaxis_title="Tissue", yaxis_title="Shannon Diversity Index",
            height=400, margin=dict(l=50, r=30, t=70, b=80), plot_bgcolor="white",
            yaxis=dict(gridcolor="rgba(0,0,0,0.08)", zeroline=False,
                       range=[0, y_max_d * 1.2 + 0.1]),
            xaxis=dict(tickangle=45),
        )
        st.plotly_chart(div_fig, use_container_width=True, key="mt_bar_diversity", config=_plotly_cfg())
        _download_plotly_as_png(div_fig, "age_diversity_analysis.png")
        
        diversity_df = pd.DataFrame({'Tissue': tissues, 'Shannon Diversity': diversities})
        create_csv_download(diversity_df, "shannon_diversity.csv", 
                           "⬇️ Download Diversity Data CSV")