import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from utils.analysis import compute_common_genes_matrix, compute_percent_overlap_matrix
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
                                 fmt="int", key=None):
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
                    row_t.append(f"{val:.1f}")
                    row_h.append(f"<b>{labels[i]}</b> ∩ <b>{labels[j]}</b><br>{colorbar_label}: {val:.1f}")
        text_vals.append(row_t)
        hover_vals.append(row_h)

    fig = go.Figure(data=go.Heatmap(
        z=z_vals, x=labels, y=labels,
        colorscale=cmap,
        text=text_vals,
        texttemplate="<b>%{text}</b>",
        textfont=dict(size=11),
        hovertext=hover_vals,
        hovertemplate="%{hovertext}<extra></extra>",
        colorbar=dict(title=colorbar_label, thickness=14, len=0.75),
        xgap=2, ygap=2,
    ))

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


def show():
    """Gene Sharing Analysis Page"""
    
    st.header("🤝 Gene Sharing Analysis")
    st.markdown("Comprehensive analysis of shared and exclusive genes between tissues and age groups.")
    
    age_groups = ["30–39", "40–49", "50–59", "60–69", "70–79"]

    # Track pre-computed (GTEx) mode so the pairwise Jaccard can be read from
    # the canonical pre-computed matrix instead of recomputed client-side.
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
            <p>Select tissues to analyze gene sharing patterns</p>
        </div>
        """, unsafe_allow_html=True)

        selected_tissues = st.multiselect(
            "Select tissues:",
            all_tissues,
            default=all_tissues[:5],
            key="gs_tissues",
        )

        if not selected_tissues or len(selected_tissues) < 2:
            st.info("👆 Please select at least 2 tissues for gene sharing analysis.")
            return

        data_parsed = {}
        for tissue in selected_tissues:
            sets_dict = get_sets_for_tissue(version, tissue, complete)
            _, counts, df = sets_to_gui_format(sets_dict)
            gene_sets_list = sets_to_gene_sets(sets_dict)
            data_parsed[tissue] = {
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
            'total_unique_genes': len(set().union(*(set().union(*d['gene_sets']) for d in data_parsed.values()))),
            'tissue_names': list(data_parsed.keys()),
            'rejected_files': [],
        }
        st.success(f"✅ {len(selected_tissues)} tissues loaded from GTEx {version}!")

    else:
        st.markdown("""
        <div class="analysis-section">
            <h3>📂 Upload Multiple Tissue Files</h3>
            <p>Upload tissue files to analyze gene sharing patterns</p>
        </div>
        """, unsafe_allow_html=True)

        uploaded_files = st.file_uploader(
            "📂 Upload STAMP .txt files", 
            type=["txt"],
            accept_multiple_files=True, 
            key="gene_sharing_files",
            help="Upload multiple tissue gene switching files"
        )

        if not uploaded_files or len(uploaded_files) < 2:
            st.info("👆 Please upload at least 2 tissue files for gene sharing analysis.")

            st.markdown("### 🔍 Gene Sharing Analysis Features")
            col1, col2 = st.columns(2)
            with col1:
                st.markdown("""
                **📊 Sharing Metrics:**
                - Absolute shared gene counts
                - Percentage overlap analysis
                - Jaccard similarity indices
                - Exclusive gene identification
                """)
            with col2:
                st.markdown("""
                **🔎 Analysis Options:**
                - Whole lifespan sharing
                - Age-specific sharing
                - Pairwise comparisons
                - Multi-tissue intersections
                """)
            return

        st.success(f"✅ {len(uploaded_files)} tissue files loaded successfully!")

        result = parse_multiple_stamp_files(uploaded_files, age_groups)
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

        data_parsed = result['data']
        summary = result['summary']

    # Converti in formato compatibile con le funzioni di analisi
    data = {tissue: data_parsed[tissue]['gene_sets'] for tissue in data_parsed.keys()}
    tissues = list(data.keys())  # Nomi già puliti
    
    # Analysis options
    st.markdown("### ⚙️ Analysis Options")
    col1, col2 = st.columns(2)
    
    with col1:
        analysis_scope = st.selectbox("📊 Analysis Scope", 
                                    ["Whole Lifespan", "Age-Specific", "Both"],
                                    help="Choose the temporal scope for analysis")
        sharing_metric = st.selectbox("📈 Sharing Metric", 
                                    ["Absolute Count", "Percentage", "Both"],
                                    help="Choose how to measure gene sharing")
    
    with col2:
        comparison_type = st.selectbox("🔍 Comparison Type", 
                                     ["All Pairwise", "Selected Pairs", "Multi-way"],
                                     help="Choose comparison methodology")
        if analysis_scope in ["Age-Specific", "Both"]:
            selected_age = st.selectbox("🎯 Age Group", age_groups, index=2)
    
    # === WHOLE LIFESPAN ANALYSIS ===
    if analysis_scope in ["Whole Lifespan", "Both"]:
        st.markdown("""
        <div class="analysis-section">
            <h2>🔗 Whole Lifespan Gene Sharing</h2>
        </div>
        """, unsafe_allow_html=True)
        
        # Calculate lifespan gene sets for each tissue
        lifespan_genes = {}
        for tissue in tissues:
            lifespan_genes[tissue] = set.union(*data[tissue])
        
        # Overall statistics
        col1, col2, col3, col4 = st.columns(4)
        
        total_unique = len(set.union(*lifespan_genes.values()))
        avg_genes_per_tissue = np.mean([len(genes) for genes in lifespan_genes.values()])
        max_genes = max([len(genes) for genes in lifespan_genes.values()])
        min_genes = min([len(genes) for genes in lifespan_genes.values()])
        
        with col1:
            st.markdown(f"""
            <div class="metric-card">
                <h3>{total_unique}</h3>
                <p>Total Unique<br>Genes</p>
            </div>
            """, unsafe_allow_html=True)
        
        with col2:
            st.markdown(f"""
            <div class="metric-card">
                <h3>{avg_genes_per_tissue:.0f}</h3>
                <p>Avg Genes<br>per Tissue</p>
            </div>
            """, unsafe_allow_html=True)
        
        with col3:
            st.markdown(f"""
            <div class="metric-card">
                <h3>{max_genes}</h3>
                <p>Max Genes<br>(One Tissue)</p>
            </div>
            """, unsafe_allow_html=True)
        
        with col4:
            st.markdown(f"""
            <div class="metric-card">
                <h3>{min_genes}</h3>
                <p>Min Genes<br>(One Tissue)</p>
            </div>
            """, unsafe_allow_html=True)
        
        # Absolute shared genes matrix (Plotly)
        if sharing_metric in ["Absolute Count", "Both"]:
            st.markdown("### 🔢 Absolute Shared Gene Counts")
            matrix_common, tissues_sorted = compute_common_genes_matrix(data)
            if matrix_common.size > 0:
                _interactive_heatmap_tissue(
                    matrix_common.astype(np.int64), tissues_sorted,
                    "Shared Genes Between Tissues (Absolute Counts)",
                    "Shared Gene Count", cmap="Purples", mask_upper=True,
                    fmt="int", key="gs_hm_abs_life"
                )
        
        # Percentage overlap matrix (Plotly)
        if sharing_metric in ["Percentage", "Both"]:
            st.markdown("### 📊 Percentage Overlap Analysis")
            matrix_percent, tissues_sorted = compute_percent_overlap_matrix(data)
            if matrix_percent.size > 0:
                _interactive_heatmap_tissue(
                    matrix_percent, tissues_sorted,
                    "Gene Overlap Percentage Between Tissues",
                    "Overlap %", cmap="teal", mask_upper=True,
                    fmt="float", key="gs_hm_pct_life"
                )
        
        # Detailed pairwise sharing table
        st.markdown("### 📋 Detailed Pairwise Sharing Statistics")
        
        # In pre-loaded mode, read the canonical J_life from the pre-computed
        # matrix so the table matches the paper / v8-vs-v10 page exactly.
        _jlife_pre = None
        _d2s = None
        if preloaded:
            try:
                from data_loader import get_jaccard_matrix, display_to_safe as _d2s
                _jlife_pre = get_jaccard_matrix(version, "life", complete)
            except Exception:
                _jlife_pre = None

        pairwise_data = []
        for i, tissue1 in enumerate(tissues_sorted):
            for j, tissue2 in enumerate(tissues_sorted):
                if i < j:  # Only upper triangle
                    genes1 = lifespan_genes[tissue1]
                    genes2 = lifespan_genes[tissue2]

                    shared = len(genes1 & genes2)
                    exclusive1 = len(genes1 - genes2)
                    exclusive2 = len(genes2 - genes1)
                    total_union = len(genes1 | genes2)
                    jaccard = shared / total_union if total_union > 0 else 0
                    # Prefer the pre-computed canonical J_life when available.
                    if _jlife_pre is not None and _d2s is not None:
                        _s1, _s2 = _d2s(tissue1), _d2s(tissue2)
                        if _s1 in _jlife_pre.index and _s2 in _jlife_pre.columns:
                            jaccard = float(_jlife_pre.loc[_s1, _s2])

                    pairwise_data.append({
                        'Tissue 1': tissue1,
                        'Tissue 2': tissue2,
                        'Shared Genes': shared,
                        f'{tissue1} Exclusive': exclusive1,
                        f'{tissue2} Exclusive': exclusive2,
                        'Jaccard Similarity': f"{jaccard:.3f}",
                        'Total Union': total_union
                    })
        
        df_pairwise = pd.DataFrame(pairwise_data)
        st.dataframe(df_pairwise, use_container_width=True)
        if _jlife_pre is not None:
            st.caption(
                f"Jaccard Similarity column = pre-computed J_life from "
                f"`output/{version}/jaccard/jaccard_life.csv` (shared/exclusive "
                "counts computed from the gene sets)."
            )
        
        # Top sharing pairs (Plotly bar)
        st.markdown("### 🏆 Top Gene Sharing Pairs")
        
        df_sorted = df_pairwise.sort_values('Shared Genes', ascending=False)
        top_pairs = df_sorted.head(5)
        
        pair_labels = [f"{row['Tissue 1']}\nvs\n{row['Tissue 2']}" 
                      for _, row in top_pairs.iterrows()]
        shared_counts = top_pairs['Shared Genes'].values
        n_p = len(top_pairs)
        bar_w = 0.28 if n_p >= 4 else (0.22 if n_p >= 2 else 0.15)

        colors_bar = px.colors.sequential.Viridis[:n_p]

        top_fig = go.Figure()
        top_fig.add_trace(go.Bar(
            x=pair_labels, y=shared_counts,
            width=[bar_w] * n_p,
            marker=dict(color=colors_bar, line=dict(color="black", width=1.2)),
            text=[f"<b>{int(c)}</b>" for c in shared_counts],
            textposition="outside",
            hovertemplate="<b>%{x}</b><br>Shared Genes: <b>%{y}</b><extra></extra>",
        ))
        y_max = max(shared_counts) if len(shared_counts) > 0 else 1
        top_fig.update_layout(
            title=dict(text="Top 5 Tissue Pairs by Shared Genes", font=dict(size=14), x=0.5),
            yaxis_title="Number of Shared Genes",
            height=400,
            margin=dict(l=50, r=30, t=70, b=80),
            plot_bgcolor="white",
            yaxis=dict(gridcolor="rgba(0,0,0,0.08)", zeroline=False,
                       range=[0, y_max * 1.20 + 3]),
        )
        st.plotly_chart(top_fig, use_container_width=True, key="gs_bar_top_pairs", config=_plotly_cfg())
        _download_plotly_as_png(top_fig, "top_sharing_pairs_lifespan.png")
    
    # === AGE-SPECIFIC ANALYSIS ===
    if analysis_scope in ["Age-Specific", "Both"]:
        st.markdown(f"""
        <div class="analysis-section">
            <h2>📅 Age-Specific Gene Sharing (Age {selected_age})</h2>
        </div>
        """, unsafe_allow_html=True)
        
        selected_idx = age_groups.index(selected_age)
        
        # Age-specific statistics
        age_genes = {tissue: data[tissue][selected_idx] for tissue in tissues}
        
        col1, col2, col3, col4 = st.columns(4)
        
        total_age_unique = len(set.union(*age_genes.values()))
        avg_age_genes = np.mean([len(genes) for genes in age_genes.values()])
        max_age_genes = max([len(genes) for genes in age_genes.values()])
        min_age_genes = min([len(genes) for genes in age_genes.values()])
        
        with col1:
            st.markdown(f"""
            <div class="metric-card">
                <h3>{total_age_unique}</h3>
                <p>Unique Genes<br>in Age {selected_age}</p>
            </div>
            """, unsafe_allow_html=True)
        
        with col2:
            st.markdown(f"""
            <div class="metric-card">
                <h3>{avg_age_genes:.0f}</h3>
                <p>Avg Genes<br>per Tissue</p>
            </div>
            """, unsafe_allow_html=True)
        
        with col3:
            st.markdown(f"""
            <div class="metric-card">
                <h3>{max_age_genes}</h3>
                <p>Max Genes<br>(One Tissue)</p>
            </div>
            """, unsafe_allow_html=True)
        
        with col4:
            st.markdown(f"""
            <div class="metric-card">
                <h3>{min_age_genes}</h3>
                <p>Min Genes<br>(One Tissue)</p>
            </div>
            """, unsafe_allow_html=True)
        
        # Age-specific shared genes matrix (Plotly)
        if sharing_metric in ["Absolute Count", "Both"]:
            st.markdown(f"### 🔢 Shared Gene Counts (Age {selected_age})")
            matrix_age_common, tissues_sorted = compute_common_genes_matrix(data, idx=selected_idx)
            if matrix_age_common.size > 0:
                _interactive_heatmap_tissue(
                    matrix_age_common.astype(np.int64), tissues_sorted,
                    f"Shared Genes — Age Group {selected_age}",
                    "Shared Gene Count", cmap="Purples", mask_upper=True,
                    fmt="int", key="gs_hm_abs_age"
                )
        
        # Age-specific percentage overlap (Plotly)
        if sharing_metric in ["Percentage", "Both"]:
            st.markdown(f"### 📊 Overlap Percentage (Age {selected_age})")
            matrix_age_percent, tissues_sorted = compute_percent_overlap_matrix(data, idx=selected_idx)
            if matrix_age_percent.size > 0:
                _interactive_heatmap_tissue(
                    matrix_age_percent, tissues_sorted,
                    f"Gene Overlap Percentage — Age {selected_age}",
                    "Overlap %", cmap="RdBu_r", mask_upper=True,
                    fmt="float", key="gs_hm_pct_age"
                )
    
    # === PAIRWISE DETAILED COMPARISON ===
    if comparison_type in ["Selected Pairs", "All Pairwise"]:
        st.markdown("""
        <div class="analysis-section">
            <h2>🔍 Detailed Pairwise Comparison</h2>
        </div>
        """, unsafe_allow_html=True)
        
        if comparison_type == "Selected Pairs":
            col1, col2 = st.columns(2)
            with col1:
                tissue1 = st.selectbox("🧪 Select First Tissue", tissues, key="pair_t1")
            with col2:
                tissue2 = st.selectbox("🧪 Select Second Tissue", 
                                      [t for t in tissues if t != tissue1], key="pair_t2")
        else:
            # Show top 3 most similar pairs
            if 'df_pairwise' in locals():
                top_3_pairs = df_sorted.head(3)
                st.markdown("### 🏆 Top 3 Most Similar Tissue Pairs")
                
                for idx, (_, row) in enumerate(top_3_pairs.iterrows()):
                    tissue1, tissue2 = row['Tissue 1'], row['Tissue 2']
                    
                    with st.expander(f"#{idx+1}: {tissue1} vs {tissue2} ({row['Shared Genes']} shared genes)"):
                        show_pairwise_analysis(tissue1, tissue2, data, age_groups, 
                                             analysis_scope, selected_age if 'selected_age' in locals() else None)
        
        if comparison_type == "Selected Pairs":
            show_pairwise_analysis(tissue1, tissue2, data, age_groups, 
                                 analysis_scope, selected_age if 'selected_age' in locals() else None)
    
    # === MULTI-WAY INTERSECTIONS ===
    if comparison_type == "Multi-way" and len(tissues) >= 3:
        st.markdown("""
        <div class="analysis-section">
            <h2>🔄 Multi-way Gene Intersections</h2>
        </div>
        """, unsafe_allow_html=True)
        
        # Allow user to select tissues for intersection
        selected_tissues = st.multiselect(
            "🎯 Select tissues for intersection analysis:",
            tissues,
            default=tissues[:3] if len(tissues) >= 3 else tissues,
            help="Choose 3+ tissues to analyze their gene intersections"
        )
        
        if len(selected_tissues) >= 3:
            show_multiway_analysis(selected_tissues, data, age_groups, 
                                 analysis_scope, selected_age if 'selected_age' in locals() else None)
    
    # === DOWNLOAD SECTION ===
    display_download_section("📥 Download Gene Sharing Analysis Results")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.markdown("**📊 Sharing Matrices**")
        if 'matrix_common' in locals() and matrix_common.size > 0:
            df_common = pd.DataFrame(matrix_common, index=tissues_sorted, columns=tissues_sorted)
            create_csv_download(df_common, "shared_genes_matrix_lifespan.csv", 
                               "⬇️ Shared Genes Matrix CSV")
        
        if 'matrix_percent' in locals() and matrix_percent.size > 0:
            df_percent = pd.DataFrame(matrix_percent, index=tissues_sorted, columns=tissues_sorted)
            create_csv_download(df_percent, "overlap_percentage_matrix_lifespan.csv", 
                               "⬇️ Overlap Percentage CSV")
    
    with col2:
        st.markdown("**📋 Detailed Tables**")
        if 'df_pairwise' in locals():
            create_csv_download(df_pairwise, "pairwise_sharing_statistics.csv", 
                               "⬇️ Pairwise Statistics CSV")
    
    with col3:
        st.markdown("**🧬 Gene Lists**")
        # Create comprehensive gene sharing summary
        summary_data = []
        for tissue in tissues:
            genes_lifespan = set.union(*data[tissue])
            for gene in genes_lifespan:
                # Count in how many tissues this gene appears
                tissues_with_gene = [t for t in tissues if gene in set.union(*data[t])]
                summary_data.append({
                    'Gene': gene,
                    'Primary_Tissue': tissue,
                    'Total_Tissues': len(tissues_with_gene),
                    'Tissues_List': ', '.join(tissues_with_gene)
                })
        
        if summary_data:
            df_gene_summary = pd.DataFrame(summary_data)
            df_gene_summary = df_gene_summary.drop_duplicates('Gene')
            create_csv_download(df_gene_summary, "gene_sharing_summary.csv", 
                               "⬇️ Gene Sharing Summary CSV")


def show_pairwise_analysis(tissue1, tissue2, data, age_groups, analysis_scope, selected_age):
    """Show detailed pairwise analysis between two tissues"""
    
    # Whole lifespan comparison
    if analysis_scope in ["Whole Lifespan", "Both"]:
        genes1_all = set.union(*data[tissue1])
        genes2_all = set.union(*data[tissue2])
        
        shared_all = genes1_all & genes2_all
        exclusive1_all = genes1_all - genes2_all
        exclusive2_all = genes2_all - genes1_all
        
        st.markdown(f"#### 🔗 Whole Lifespan: {tissue1} vs {tissue2}")
        
        col1, col2, col3 = st.columns(3)
        with col1:
            st.markdown(f"""
            <div class="metric-card">
                <h3>{len(shared_all)}</h3>
                <p>Shared Genes</p>
            </div>
            """, unsafe_allow_html=True)
        
        with col2:
            st.markdown(f"""
            <div class="metric-card">
                <h3>{len(exclusive1_all)}</h3>
                <p>{tissue1} Exclusive</p>
            </div>
            """, unsafe_allow_html=True)
        
        with col3:
            st.markdown(f"""
            <div class="metric-card">
                <h3>{len(exclusive2_all)}</h3>
                <p>{tissue2} Exclusive</p>
            </div>
            """, unsafe_allow_html=True)
        
        # Pie chart (Plotly)
        if len(shared_all) + len(exclusive1_all) + len(exclusive2_all) > 0:
            sizes = [len(shared_all), len(exclusive1_all), len(exclusive2_all)]
            labels_pie = [f'Shared ({len(shared_all)})',
                     f'{tissue1} Exclusive ({len(exclusive1_all)})',
                     f'{tissue2} Exclusive ({len(exclusive2_all)})']
            colors_pie = ['#2ecc71', '#3498db', '#e74c3c']

            pie_fig = go.Figure(data=go.Pie(
                labels=labels_pie, values=sizes,
                marker=dict(colors=colors_pie, line=dict(color='white', width=2)),
                textinfo='percent+label', textposition='inside',
                hovertemplate='<b>%{label}</b><br>Genes: %{value}<br>Percentage: %{percent}<extra></extra>',
            ))
            pie_fig.update_layout(
                title=dict(text=f"Gene Distribution: {tissue1} vs {tissue2}",
                           font=dict(size=13), x=0.5),
                height=350, margin=dict(l=30, r=30, t=60, b=30),
                showlegend=False,
            )
            _left, _c, _r = st.columns([1, 2, 1])
            with _c:
                st.plotly_chart(pie_fig, use_container_width=True, config=_plotly_cfg())
            _download_plotly_as_png(pie_fig, f"pie_chart_{tissue1}_vs_{tissue2}_lifespan.png")
    
    # Age-specific comparison
    if analysis_scope in ["Age-Specific", "Both"] and selected_age:
        selected_idx = age_groups.index(selected_age)
        
        genes1_age = data[tissue1][selected_idx]
        genes2_age = data[tissue2][selected_idx]
        
        shared_age = genes1_age & genes2_age
        exclusive1_age = genes1_age - genes2_age
        exclusive2_age = genes2_age - genes1_age
        
        st.markdown(f"#### 📅 Age {selected_age}: {tissue1} vs {tissue2}")
        
        col1, col2, col3 = st.columns(3)
        with col1:
            st.markdown(f"""
            <div class="metric-card">
                <h3>{len(shared_age)}</h3>
                <p>Shared Genes</p>
            </div>
            """, unsafe_allow_html=True)
        
        with col2:
            st.markdown(f"""
            <div class="metric-card">
                <h3>{len(exclusive1_age)}</h3>
                <p>{tissue1} Exclusive</p>
            </div>
            """, unsafe_allow_html=True)
        
        with col3:
            st.markdown(f"""
            <div class="metric-card">
                <h3>{len(exclusive2_age)}</h3>
                <p>{tissue2} Exclusive</p>
            </div>
            """, unsafe_allow_html=True)


def show_multiway_analysis(selected_tissues, data, age_groups, analysis_scope, selected_age):
    """Show multi-way intersection analysis"""
    
    st.markdown(f"### 🔄 Intersection Analysis: {', '.join(selected_tissues)}")
    
    if analysis_scope in ["Whole Lifespan", "Both"]:
        # Calculate intersections for whole lifespan
        tissue_gene_sets = [set.union(*data[tissue]) for tissue in selected_tissues]
        
        # Core intersection (genes in ALL tissues)
        core_intersection = set.intersection(*tissue_gene_sets)
        
        # Union (genes in ANY tissue)
        total_union = set.union(*tissue_gene_sets)
        
        st.markdown(f"#### 🔗 Whole Lifespan Multi-way Analysis")
        
        col1, col2, col3 = st.columns(3)
        with col1:
            st.markdown(f"""
            <div class="metric-card">
                <h3>{len(core_intersection)}</h3>
                <p>Core Genes<br>(All Tissues)</p>
            </div>
            """, unsafe_allow_html=True)
        
        with col2:
            st.markdown(f"""
            <div class="metric-card">
                <h3>{len(total_union)}</h3>
                <p>Total Unique<br>Genes</p>
            </div>
            """, unsafe_allow_html=True)
        
        with col3:
            overlap_ratio = len(core_intersection) / len(total_union) if len(total_union) > 0 else 0
            st.markdown(f"""
            <div class="metric-card">
                <h3>{overlap_ratio:.2%}</h3>
                <p>Core Overlap<br>Ratio</p>
            </div>
            """, unsafe_allow_html=True)
        
        if len(core_intersection) > 0:
            st.markdown("**🎯 Core Genes (present in all selected tissues):**")
            st.text(", ".join(sorted(core_intersection)[:50]) + ("..." if len(core_intersection) > 50 else ""))
            
            # Download core genes
            core_df = pd.DataFrame(sorted(core_intersection), columns=['Gene'])
            create_csv_download(core_df, f"core_genes_{'_'.join(selected_tissues[:3])}.csv", 
                               "⬇️ Download Core Genes CSV")