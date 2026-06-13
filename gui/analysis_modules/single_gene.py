import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from utils.parsing import parse_multiple_stamp_files  # USO PARSING CENTRALIZZATO
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


def show():
    """Single Gene Analysis Page"""
    
    st.header("🔍 Single Gene Analysis")
    st.markdown("Track individual gene expression patterns across tissues and age groups.")
    
    age_groups = ["30–39", "40–49", "50–59", "60–69", "70–79"]
    
    # ── Data source toggle ─────────────────────────────────────────
    import sys as _sys
    from pathlib import Path as _Path
    _sys.path.insert(0, str(_Path(__file__).resolve().parent.parent))
    from data_loader import get_available_tissues, get_sets_for_tissue, sets_to_gui_format, sets_to_gene_sets, complete_mode

    # Pre-loaded (GTEx) mode removed by design: users always upload their files.
    data_source = "📂 Upload files"

    complete = complete_mode()

    if data_source == "📦 Pre-loaded (GTEx)":
        version = st.session_state.get("gtex_version", "v10")
        all_tissues = get_available_tissues(version, complete)

        st.markdown(f"""
        <div class="analysis-section">
            <h3>🧪 Select Tissues — GTEx {version}</h3>
            <p>Select tissues to track gene expression patterns</p>
        </div>
        """, unsafe_allow_html=True)

        selected_tissues = st.multiselect(
            "Select tissues:",
            all_tissues,
            default=all_tissues[:5],
            key="sg_tissues",
        )

        if not selected_tissues or len(selected_tissues) < 2:
            st.info("👆 Please select at least 2 tissues for single gene analysis.")
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
        st.success(f"✅ {len(selected_tissues)} tissues loaded from GTEx {version}!")

    else:
        st.markdown("""
        <div class="analysis-section">
            <h3>📂 Upload Multiple Tissue Files</h3>
            <p>Upload tissue files to track gene expression patterns</p>
        </div>
        """, unsafe_allow_html=True)

        uploaded_files = st.file_uploader(
            "📂 Upload STAMP .txt files", 
            type=["txt"],
            accept_multiple_files=True, 
            key="single_gene_files",
            help="Upload multiple tissue gene switching files"
        )

        if not uploaded_files or len(uploaded_files) < 2:
            st.info("👆 Please upload at least 2 tissue files for single gene analysis.")

            st.markdown("### 🔍 Single Gene Analysis Features")
            col1, col2 = st.columns(2)
            with col1:
                st.markdown("""
                **🎯 Gene Tracking:**
                - Gene presence across tissues
                - Age-specific expression patterns
                - Tissue-specific appearance
                - Expression timeline analysis
                """)
            with col2:
                st.markdown("""
                **📊 Visualization Options:**
                - Presence heatmaps
                - Timeline plots
                - Tissue distribution charts
                - Age pattern analysis
                """)

            st.markdown("### 💡 Example Gene Searches")
            example_genes = ["APOE", "TP53", "BRCA1", "EGFR", "MYC", "PTEN", "RB1", "VHL"]
            st.info(f"Common genes to search: {', '.join(example_genes)}")
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

    # Estrai tessuti e converti formato per compatibilità
    tissues = list(data_parsed.keys())
    data = {}
    for tissue in tissues:
        data[tissue] = data_parsed[tissue]['gene_sets']
    
    # Create comprehensive gene list
    all_genes = set()
    for tissue in tissues:
        for age_set in data[tissue]:
            all_genes.update(age_set)
    
    all_genes = sorted(list(all_genes))
    
    # Gene search section
    st.markdown("""
    <div class="analysis-section">
        <h2>🔬 Gene Search and Analysis</h2>
    </div>
    """, unsafe_allow_html=True)
    
    col1, col2 = st.columns([2, 1])
    
    with col1:
        # Gene input methods
        search_method = st.radio(
            "🔍 Gene Search Method:",
            ["Manual Input", "Select from List", "Multiple Genes"],
            help="Choose how to specify genes for analysis"
        )
        
        if search_method == "Manual Input":
            gene_input = st.text_input(
                "🧬 Enter gene name (e.g., APOE):",
                placeholder="Type gene name...",
                help="Enter a single gene name to analyze"
            ).strip().upper()
            genes_to_analyze = [gene_input] if gene_input else []
            
        elif search_method == "Select from List":
            gene_input = st.selectbox(
                "🧬 Select gene from available list:",
                [""] + all_genes,
                help="Choose from genes present in your datasets"
            )
            genes_to_analyze = [gene_input] if gene_input else []
            
        else:  # Multiple Genes
            genes_input = st.text_area(
                "🧬 Enter multiple gene names (one per line):",
                placeholder="APOE\nTP53\nBRCA1",
                help="Enter multiple gene names, one per line"
            )
            genes_to_analyze = [g.strip().upper() for g in genes_input.split('\n') if g.strip()]
    
    with col2:
        st.markdown("### 📊 Dataset Overview")
        st.markdown(f"""
        <div class="metric-card">
            <h3>{len(all_genes)}</h3>
            <p>Total Unique Genes</p>
        </div>
        """, unsafe_allow_html=True)
        
        st.markdown(f"""
        <div class="metric-card">
            <h3>{len(tissues)}</h3>
            <p>Tissues Available</p>
        </div>
        """, unsafe_allow_html=True)
        
        st.markdown(f"""
        <div class="metric-card">
            <h3>{len(age_groups)}</h3>
            <p>Age Groups</p>
        </div>
        """, unsafe_allow_html=True)
    
    # Analysis options
    if genes_to_analyze and genes_to_analyze != ['']:
        st.markdown("### ⚙️ Analysis Options")
        
        col1, col2 = st.columns(2)
        with col1:
            show_heatmap = st.checkbox("🔥 Show presence heatmap", value=True)
            show_timeline = st.checkbox("📈 Show timeline analysis", value=True)
        with col2:
            show_statistics = st.checkbox("📊 Show detailed statistics", value=True)
            show_tissue_distribution = st.checkbox("🧪 Show tissue distribution", value=True)
        
        # Filter genes that exist in the dataset
        valid_genes = [gene for gene in genes_to_analyze if gene in all_genes]
        invalid_genes = [gene for gene in genes_to_analyze if gene not in all_genes]
        
        if invalid_genes:
            st.warning(f"⚠️ Genes not found in dataset: {', '.join(invalid_genes)}")
        
        if not valid_genes:
            st.error("❌ No valid genes found in the dataset. Please check gene names.")
            return
        
        st.success(f"✅ Analyzing {len(valid_genes)} gene(s): {', '.join(valid_genes)}")
        
        # === MAIN ANALYSIS ===
        for gene_idx, gene in enumerate(valid_genes):
            if len(valid_genes) > 1:
                st.markdown(f"""
                <div class="analysis-section">
                    <h2>🧬 Analysis for Gene: {gene}</h2>
                </div>
                """, unsafe_allow_html=True)
            
            # Create presence matrix for this gene
            presence_matrix = np.zeros((len(tissues), len(age_groups)), dtype=int)
            presence_text = pd.DataFrame(index=tissues, columns=age_groups)
            
            for i, tissue in enumerate(tissues):
                for j, age_idx in enumerate(range(len(age_groups))):
                    is_present = gene in data[tissue][age_idx]
                    presence_matrix[i, j] = 1 if is_present else 0
                    presence_text.iloc[i, j] = "✅" if is_present else "❌"
            
            # Gene statistics
            total_appearances = np.sum(presence_matrix)
            tissues_with_gene = np.sum(np.any(presence_matrix, axis=1))
            ages_with_gene = np.sum(np.any(presence_matrix, axis=0))
            
            # Display statistics
            col1, col2, col3, col4 = st.columns(4)
            with col1:
                st.markdown(f"""
                <div class="metric-card">
                    <h3>{total_appearances}</h3>
                    <p>Total<br>Appearances</p>
                </div>
                """, unsafe_allow_html=True)
            
            with col2:
                st.markdown(f"""
                <div class="metric-card">
                    <h3>{tissues_with_gene}</h3>
                    <p>Tissues with<br>Gene</p>
                </div>
                """, unsafe_allow_html=True)
            
            with col3:
                st.markdown(f"""
                <div class="metric-card">
                    <h3>{ages_with_gene}</h3>
                    <p>Age Groups<br>with Gene</p>
                </div>
                """, unsafe_allow_html=True)
            
            with col4:
                expression_rate = (total_appearances / (len(tissues) * len(age_groups))) * 100
                st.markdown(f"""
                <div class="metric-card">
                    <h3>{expression_rate:.1f}%</h3>
                    <p>Expression<br>Rate</p>
                </div>
                """, unsafe_allow_html=True)
            
            # Presence heatmap (Plotly)
            if show_heatmap:
                st.markdown(f"### 🔥 Presence Heatmap - {gene}")
                
                # Build hover text & display text
                hover = []
                text_disp = []
                for i, tissue in enumerate(tissues):
                    row_h = []
                    row_t = []
                    for j, age in enumerate(age_groups):
                        val = presence_matrix[i, j]
                        symbol = "✅ Present" if val == 1 else "❌ Absent"
                        row_h.append(f"<b>{tissue}</b> × <b>{age}</b><br>{symbol}")
                        row_t.append("✅" if val == 1 else "❌")
                    hover.append(row_h)
                    text_disp.append(row_t)

                hm_fig = go.Figure(data=go.Heatmap(
                    z=presence_matrix,
                    x=age_groups, y=tissues,
                    colorscale="RdYlGn",
                    text=text_disp,
                    texttemplate="%{text}",
                    textfont=dict(size=14),
                    hovertext=hover,
                    hovertemplate="%{hovertext}<extra></extra>",
                    colorbar=dict(title="Presence", thickness=14, len=0.75,
                                  tickvals=[0, 1], ticktext=["Absent", "Present"]),
                    xgap=2, ygap=2,
                    zmin=0, zmax=1,
                ))
                hm_fig.update_layout(
                    title=dict(text=f"Gene Presence Pattern: {gene}",
                               font=dict(size=14), x=0.5),
                    height=max(350, len(tissues) * 40 + 120),
                    margin=dict(l=100, r=40, t=60, b=60),
                    yaxis=dict(autorange="reversed"),
                    plot_bgcolor="white",
                )
                _left, _centre, _right = st.columns([1, 3, 1])
                with _centre:
                    st.plotly_chart(hm_fig, use_container_width=True,
                                   key=f"sg_hm_{gene}_{gene_idx}", config=_plotly_cfg())
                _download_plotly_as_png(hm_fig, f"heatmap_{gene}_presence.png")
            
            # Timeline analysis
            if show_timeline:
                st.markdown(f"### 📈 Timeline Analysis - {gene}")
                
                # Count appearances per age group
                age_counts = np.sum(presence_matrix, axis=0)
                
                # Line plot (Plotly)
                age_counts_list = age_counts.tolist()

                line_fig = go.Figure()
                line_fig.add_trace(go.Scatter(
                    x=age_groups, y=age_counts_list,
                    mode='lines+markers+text',
                    line=dict(color='#e74c3c', width=3),
                    marker=dict(size=9, color='#e74c3c'),
                    fill='tozeroy', fillcolor='rgba(231,76,60,0.18)',
                    text=[str(int(c)) for c in age_counts_list],
                    textposition='top center',
                    hovertemplate='<b>%{x}</b><br>Tissues: <b>%{y}</b><extra></extra>',
                ))
                line_fig.update_layout(
                    title=dict(text=f"Gene Expression Timeline: {gene}",
                               font=dict(size=14), x=0.5),
                    xaxis_title='Age Group', yaxis_title='Number of Tissues',
                    height=350, margin=dict(l=50, r=30, t=60, b=50),
                    plot_bgcolor='white',
                    yaxis=dict(gridcolor='rgba(0,0,0,0.08)', zeroline=False,
                               range=[0, len(tissues) + 0.5]),
                    xaxis=dict(gridcolor='rgba(0,0,0,0.08)'),
                )
                st.plotly_chart(line_fig, use_container_width=True,
                               key=f'sg_line_{gene}_{gene_idx}', config=_plotly_cfg())
                _download_plotly_as_png(line_fig, f"timeline_{gene}_line.png")
                
                # Bar plot of tissue distribution (Plotly horizontal bar)
                tissue_counts = np.sum(presence_matrix, axis=1)
                
                barh_fig = go.Figure()
                barh_fig.add_trace(go.Bar(
                    y=tissues, x=tissue_counts.tolist(), orientation="h",
                    marker=dict(color=tissue_counts.tolist(), colorscale="Viridis",
                                line=dict(color="black", width=1)),
                    text=[f"<b>{int(c)}</b>" for c in tissue_counts],
                    textposition="outside",
                    hovertemplate="<b>%{y}</b><br>Age Groups: <b>%{x}</b><extra></extra>",
                ))
                barh_fig.update_layout(
                    title=dict(text=f"Tissue Distribution: {gene}",
                               font=dict(size=14), x=0.5),
                    xaxis_title="Number of Age Groups",
                    height=max(300, len(tissues) * 35 + 100),
                    margin=dict(l=120, r=60, t=60, b=50), plot_bgcolor="white",
                    xaxis=dict(gridcolor="rgba(0,0,0,0.08)", range=[0, len(age_groups) + 0.5]),
                    yaxis=dict(autorange="reversed"),
                )
                st.plotly_chart(barh_fig, use_container_width=True,
                               key=f"sg_barh_tissue_{gene}_{gene_idx}", config=_plotly_cfg())
                _download_plotly_as_png(barh_fig, f"tissue_distribution_{gene}.png")
            
            # Detailed statistics table
            if show_statistics:
                st.markdown(f"### 📊 Detailed Statistics - {gene}")
                
                stats_data = []
                for tissue in tissues:
                    tissue_appearances = np.sum(presence_matrix[tissues.index(tissue), :])
                    age_list = [age_groups[j] for j in range(len(age_groups)) 
                              if presence_matrix[tissues.index(tissue), j] == 1]
                    
                    stats_data.append({
                        'Tissue': tissue,
                        'Appearances': tissue_appearances,
                        'Expression Rate': f"{(tissue_appearances/len(age_groups)*100):.1f}%",
                        'Age Groups': ', '.join(age_list) if age_list else 'None'
                    })
                
                df_stats = pd.DataFrame(stats_data)
                st.dataframe(df_stats, use_container_width=True)
                
                # Age group analysis
                st.markdown("#### 📅 Age Group Analysis")
                age_stats_data = []
                for i, age in enumerate(age_groups):
                    age_appearances = np.sum(presence_matrix[:, i])
                    tissue_list = [tissues[j] for j in range(len(tissues)) 
                                 if presence_matrix[j, i] == 1]
                    
                    age_stats_data.append({
                        'Age Group': age,
                        'Tissues with Gene': age_appearances,
                        'Expression Rate': f"{(age_appearances/len(tissues)*100):.1f}%",
                        'Tissues': ', '.join(tissue_list) if tissue_list else 'None'
                    })
                
                df_age_stats = pd.DataFrame(age_stats_data)
                st.dataframe(df_age_stats, use_container_width=True)
                create_csv_download(df_age_stats, f"age_expression_{gene}.csv", "⬇️ Download Age Expression Data (CSV)")
            
            # Tissue distribution pie chart + bar
            if show_tissue_distribution and total_appearances > 0:
                st.markdown(f"### 🧪 Tissue Expression Distribution - {gene}")
                
                expressing_tissues = [(tissues[i], np.sum(presence_matrix[i, :])) 
                                    for i in range(len(tissues)) 
                                    if np.sum(presence_matrix[i, :]) > 0]
                
                if expressing_tissues:
                    # Pie chart (Plotly)
                    tissue_names = [t[0] for t in expressing_tissues]
                    counts = [t[1] for t in expressing_tissues]
                    colors_pie = px.colors.qualitative.Set3[:len(expressing_tissues)]

                    pie_fig = go.Figure(data=go.Pie(
                        labels=tissue_names, values=counts,
                        marker=dict(colors=colors_pie, line=dict(color='white', width=2)),
                        textinfo='percent+label', textposition='inside',
                        hovertemplate='<b>%{label}</b><br>Age Groups: %{value}<br>Percentage: %{percent}<extra></extra>',
                    ))
                    pie_fig.update_layout(
                        title=dict(text=f"Expression Distribution: {gene}",
                                   font=dict(size=13), x=0.5),
                        height=350, margin=dict(l=30, r=30, t=60, b=30),
                        showlegend=False,
                    )
                    _left, _c, _r = st.columns([1, 2, 1])
                    with _c:
                        st.plotly_chart(pie_fig, use_container_width=True,
                                       key=f'sg_pie_{gene}_{gene_idx}', config=_plotly_cfg())
                    _download_plotly_as_png(pie_fig, f"pie_{gene}_tissues.png")

                    # Bar chart (Plotly)
                    n_e = len(expressing_tissues)
                    bar_w = 0.28 if n_e >= 4 else (0.22 if n_e >= 2 else 0.15)
                    colors_bar = px.colors.qualitative.Set3[:n_e]

                    expr_fig = go.Figure()
                    expr_fig.add_trace(go.Bar(
                        x=tissue_names, y=counts,
                        width=[bar_w] * n_e,
                        marker=dict(color=colors_bar, line=dict(color="black", width=1.2)),
                        text=[f"<b>{c}</b>" for c in counts],
                        textposition="outside",
                        hovertemplate="<b>%{x}</b><br>Age Groups: <b>%{y}</b><extra></extra>",
                    ))
                    y_max = max(counts) if counts else 1
                    expr_fig.update_layout(
                        title=dict(text=f"Expression Frequency by Tissue: {gene}",
                                   font=dict(size=14), x=0.5),
                        yaxis_title="Number of Age Groups",
                        xaxis_title="Tissue",
                        height=380,
                        margin=dict(l=50, r=30, t=60, b=80), plot_bgcolor="white",
                        yaxis=dict(gridcolor="rgba(0,0,0,0.08)", zeroline=False,
                                   range=[0, y_max * 1.2 + 0.5]),
                        xaxis=dict(tickangle=45),
                    )
                    st.plotly_chart(expr_fig, use_container_width=True,
                                   key=f"sg_bar_expr_{gene}_{gene_idx}", config=_plotly_cfg())
                    _download_plotly_as_png(expr_fig, f"distribution_{gene}_tissues.png")
            
            # Download section for this gene
            if len(valid_genes) > 1:
                st.markdown(f"#### 📥 Downloads for {gene}")
            
            col1, col2 = st.columns(2)
            with col1:
                create_csv_download(presence_text, f"presence_matrix_{gene}.csv", 
                                   f"⬇️ {gene} Presence Matrix CSV")
            with col2:
                if 'df_stats' in locals():
                    create_csv_download(df_stats, f"statistics_{gene}.csv", 
                                       f"⬇️ {gene} Statistics CSV")
            
            if gene_idx < len(valid_genes) - 1:
                st.markdown("---")
        
        # === MULTI-GENE COMPARISON ===
        if len(valid_genes) > 1:
            st.markdown("""
            <div class="analysis-section">
                <h2>🔄 Multi-Gene Comparison</h2>
            </div>
            """, unsafe_allow_html=True)
            
            comparison_data = []
            for gene in valid_genes:
                gene_presence = np.zeros((len(tissues), len(age_groups)), dtype=int)
                for i, tissue in enumerate(tissues):
                    for j, age_idx in enumerate(range(len(age_groups))):
                        gene_presence[i, j] = 1 if gene in data[tissue][age_idx] else 0
                
                total_appearances = np.sum(gene_presence)
                tissues_with_gene = np.sum(np.any(gene_presence, axis=1))
                ages_with_gene = np.sum(np.any(gene_presence, axis=0))
                expression_rate = (total_appearances / (len(tissues) * len(age_groups))) * 100
                
                comparison_data.append({
                    'Gene': gene,
                    'Total Appearances': total_appearances,
                    'Tissues with Gene': tissues_with_gene,
                    'Age Groups with Gene': ages_with_gene,
                    'Expression Rate (%)': f"{expression_rate:.1f}%"
                })
            
            df_comparison = pd.DataFrame(comparison_data)
            st.dataframe(df_comparison, use_container_width=True)
            
            # Comparison grouped bar chart (Plotly)
            appearances = [int(row['Total Appearances']) for _, row in df_comparison.iterrows()]
            tissues_count = [int(row['Tissues with Gene']) for _, row in df_comparison.iterrows()]
            ages_count = [int(row['Age Groups with Gene']) for _, row in df_comparison.iterrows()]
            
            n_g = len(valid_genes)
            bar_w = 0.22 if n_g >= 4 else (0.16 if n_g >= 2 else 0.12)

            comp_fig = go.Figure()
            comp_fig.add_trace(go.Bar(
                name='Total Appearances', x=valid_genes, y=appearances,
                width=[bar_w] * n_g,
                marker=dict(color="rgba(52,152,219,0.85)", line=dict(color="rgb(8,48,107)", width=1)),
                text=[f"<b>{v}</b>" for v in appearances], textposition="outside",
                hovertemplate="<b>%{x}</b><br>Total Appearances: <b>%{y}</b><extra></extra>",
            ))
            comp_fig.add_trace(go.Bar(
                name='Tissues with Gene', x=valid_genes, y=tissues_count,
                width=[bar_w] * n_g,
                marker=dict(color="rgba(231,76,60,0.85)", line=dict(color="darkred", width=1)),
                text=[f"<b>{v}</b>" for v in tissues_count], textposition="outside",
                hovertemplate="<b>%{x}</b><br>Tissues: <b>%{y}</b><extra></extra>",
            ))
            comp_fig.add_trace(go.Bar(
                name='Age Groups with Gene', x=valid_genes, y=ages_count,
                width=[bar_w] * n_g,
                marker=dict(color="rgba(46,204,113,0.85)", line=dict(color="darkgreen", width=1)),
                text=[f"<b>{v}</b>" for v in ages_count], textposition="outside",
                hovertemplate="<b>%{x}</b><br>Age Groups: <b>%{y}</b><extra></extra>",
            ))
            
            all_vals = appearances + tissues_count + ages_count
            y_max = max(all_vals) if all_vals else 1
            comp_fig.update_layout(
                barmode="group",
                title=dict(text="Multi-Gene Expression Comparison", font=dict(size=14), x=0.5),
                xaxis_title="Genes", yaxis_title="Count",
                height=400, margin=dict(l=50, r=30, t=60, b=50),
                plot_bgcolor="white",
                yaxis=dict(gridcolor="rgba(0,0,0,0.08)", zeroline=False,
                           range=[0, y_max * 1.2 + 1]),
                legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5),
            )
            st.plotly_chart(comp_fig, use_container_width=True, key="sg_bar_multigene", config=_plotly_cfg())
            _download_plotly_as_png(comp_fig, "multi_gene_comparison.png")
            
            create_csv_download(df_comparison, "multi_gene_comparison.csv", 
                               "⬇️ Multi-Gene Comparison CSV")
        
        # === GLOBAL DOWNLOAD SECTION ===
        if len(valid_genes) > 0:
            display_download_section("📥 Download All Results")
            
            summary_data = {
                'Analysis_Type': ['Single Gene Analysis'],
                'Genes_Analyzed': [', '.join(valid_genes)],
                'Number_of_Genes': [len(valid_genes)],
                'Number_of_Tissues': [len(tissues)],
                'Number_of_Age_Groups': [len(age_groups)],
                'Total_Data_Points': [len(tissues) * len(age_groups) * len(valid_genes)]
            }
            summary_df = pd.DataFrame(summary_data)
            create_csv_download(summary_df, "gene_analysis_summary.csv", 
                               "⬇️ Analysis Summary CSV")
    
    else:
        st.info("👆 Please enter at least one gene name to start the analysis.")
        
        if all_genes:
            st.markdown("### 📊 Available Genes Overview")
            
            sample_size = min(20, len(all_genes))
            sample_genes = np.random.choice(all_genes, sample_size, replace=False)
            
            st.markdown(f"**Sample of available genes ({sample_size} of {len(all_genes)}):**")
            st.text(", ".join(sorted(sample_genes)))
            
            gene_frequency = {}
            for gene in all_genes:
                count = 0
                for tissue in tissues:
                    for age_set in data[tissue]:
                        if gene in age_set:
                            count += 1
                gene_frequency[gene] = count
            
            top_genes = sorted(gene_frequency.items(), key=lambda x: x[1], reverse=True)[:10]
            
            if top_genes:
                st.markdown("### 🏆 Most Frequently Expressed Genes")
                freq_df = pd.DataFrame(top_genes, columns=['Gene', 'Frequency'])
                freq_df['Expression Rate (%)'] = (freq_df['Frequency'] / (len(tissues) * len(age_groups)) * 100).round(1)
                st.dataframe(freq_df, use_container_width=True)
                create_csv_download(freq_df, "top_genes_frequency.csv", "⬇️ Download Top Genes Frequency (CSV)")