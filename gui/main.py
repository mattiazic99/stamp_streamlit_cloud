import streamlit as st
import pandas as pd
import numpy as np

# Import page modules
from analysis_modules import (
    welcome,
    upload_analysis,
    tissue_comparison,
    multi_tissue,
    age_specific,
    gene_sharing,
    single_gene,
    group_comparison,
    stamp_generator,
    threshold_sweep,
    migration_page,
    aging_timeline,
    switch_direction,
)

# Configure Streamlit page
st.set_page_config(
    page_title="STAMP - Gene Switching Explorer", 
    layout="wide",
    initial_sidebar_state="expanded"
)

# Initialize session state for navigation if not exists
if 'current_page' not in st.session_state:
    st.session_state.current_page = "welcome"

# Custom CSS for enhanced styling
st.markdown("""
<style>
    /* Sidebar generale */
    .css-1d391kg {
        background: linear-gradient(180deg, #f8f9fa 0%, #e9ecef 100%);
        border-right: 3px solid #667eea;
    }
    
    /* Header principale */
    .main-header {
        background: linear-gradient(90deg, #667eea 0%, #764ba2 100%);
        padding: 1.5rem;
        border-radius: 10px;
        color: white;
        text-align: center;
        margin-bottom: 2rem;
        box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
    }
    .main-header h1 {
        margin: 0;
        font-size: 2.5rem;
        text-shadow: 2px 2px 4px rgba(0,0,0,0.3);
    }
    .main-header p {
        margin: 0.5rem 0 0 0;
        font-size: 1.1rem;
        opacity: 0.9;
    }
    
    /* Analysis sections styling */
    .analysis-section {
        background: #f8f9fa;
        padding: 1.5rem;
        border-radius: 10px;
        margin: 1rem 0;
        border-left: 5px solid #667eea;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    }
    .metric-card {
        background: white;
        padding: 1rem;
        border-radius: 8px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        text-align: center;
        margin: 0.5rem 0;
    }
    .download-section {
        background: #e8f5e8;
        padding: 1rem;
        border-radius: 8px;
        margin: 1rem 0;
        border: 1px solid #c3e6c3;
    }
    .stSelectbox > div > div {
        background-color: white;
    }
    
    /* Navigation Section Styling - MIGLIORATO */
    .nav-section {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        padding: 25px;
        border-radius: 20px;
        margin: 25px 0;
        box-shadow: 0 12px 35px rgba(102, 126, 234, 0.4);
        border: 2px solid rgba(255, 255, 255, 0.1);
        backdrop-filter: blur(10px);
        position: relative;
        overflow: hidden;
    }
    
    .nav-section::before {
        content: '';
        position: absolute;
        top: 0;
        left: 0;
        right: 0;
        height: 3px;
        background: linear-gradient(90deg, rgba(255,255,255,0.8) 0%, rgba(255,255,255,0.3) 100%);
    }
    
    .nav-title {
        color: white;
        font-size: 1.5rem;
        font-weight: bold;
        margin-bottom: 18px;
        text-align: center;
        text-shadow: 2px 2px 6px rgba(0,0,0,0.4);
        display: flex;
        align-items: center;
        justify-content: center;
        gap: 12px;
        letter-spacing: 0.5px;
    }
    
    .nav-subtitle {
        color: rgba(255, 255, 255, 0.9);
        font-size: 1rem;
        text-align: center;
        margin-bottom: 25px;
        font-style: italic;
        font-weight: 300;
        text-shadow: 1px 1px 3px rgba(0,0,0,0.3);
    }
    
    /* Generator Section Styling - MIGLIORATO */
    .generator-section {
        background: linear-gradient(135deg, #ff6b35 0%, #f7931e 100%);
        padding: 25px;
        border-radius: 20px;
        margin: 25px 0;
        box-shadow: 0 12px 35px rgba(255, 107, 53, 0.4);
        border: 2px solid rgba(255, 255, 255, 0.1);
        backdrop-filter: blur(10px);
        position: relative;
        overflow: hidden;
    }
    
    .generator-section::before {
        content: '';
        position: absolute;
        top: 0;
        left: 0;
        right: 0;
        height: 3px;
        background: linear-gradient(90deg, rgba(255,255,255,0.8) 0%, rgba(255,255,255,0.3) 100%);
    }
    
    .generator-title {
        color: white;
        font-size: 1.5rem;
        font-weight: bold;
        margin-bottom: 18px;
        text-align: center;
        text-shadow: 2px 2px 6px rgba(0,0,0,0.4);
        display: flex;
        align-items: center;
        justify-content: center;
        gap: 12px;
        letter-spacing: 0.5px;
    }
    
    .generator-subtitle {
        color: rgba(255, 255, 255, 0.9);
        font-size: 1rem;
        text-align: center;
        margin-bottom: 25px;
        font-style: italic;
        font-weight: 300;
        text-shadow: 1px 1px 3px rgba(0,0,0,0.3);
    }
    
    /* Bottoni di navigazione - ELEGANTI */
    .stButton > button {
        width: 100% !important;
        margin: 8px 0 !important;
        padding: 16px 20px !important;
        border-radius: 15px !important;
        font-weight: 600 !important;
        font-size: 1rem !important;
        transition: all 0.3s ease !important;
        border: 2px solid transparent !important;
        box-shadow: 0 4px 12px rgba(0, 0, 0, 0.1) !important;
        position: relative !important;
        overflow: hidden !important;
    }
    
    .stButton > button:hover {
        transform: translateY(-2px) !important;
        box-shadow: 0 8px 25px rgba(0, 0, 0, 0.15) !important;
    }
    
    .stButton > button::before {
        content: '';
        position: absolute;
        top: 0;
        left: -100%;
        width: 100%;
        height: 100%;
        background: linear-gradient(90deg, transparent, rgba(255,255,255,0.2), transparent);
        transition: left 0.5s;
    }
    
    .stButton > button:hover::before {
        left: 100%;
    }
    
    /* Bottoni secondary (non attivi) */
    .stButton > button[kind="secondary"] {
        background: rgba(255, 255, 255, 0.9) !important;
        color: #667eea !important;
        border: 2px solid rgba(102, 126, 234, 0.2) !important;
    }
    
    .stButton > button[kind="secondary"]:hover {
        background: rgba(102, 126, 234, 0.1) !important;
        border: 2px solid #667eea !important;
    }
    
    /* Bottoni primary (attivi) */
    .stButton > button[kind="primary"] {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%) !important;
        color: white !important;
        border: 2px solid rgba(255, 255, 255, 0.3) !important;
        box-shadow: 0 6px 20px rgba(102, 126, 234, 0.4) !important;
    }
    
    .stButton > button[kind="primary"]:hover {
        background: linear-gradient(135deg, #5a6fd8 0%, #6a4190 100%) !important;
        box-shadow: 0 8px 30px rgba(102, 126, 234, 0.5) !important;
    }
    
    /* Info sections - MIGLIORATE */
    .info-section {
        background: linear-gradient(135deg, #ffffff 0%, #f8f9fa 100%);
        padding: 25px;
        border-radius: 20px;
        margin: 25px 0;
        box-shadow: 0 8px 25px rgba(0, 0, 0, 0.12);
        border-left: 6px solid #667eea;
        border: 2px solid rgba(102, 126, 234, 0.1);
        position: relative;
        overflow: hidden;
    }
    
    .info-section::before {
        content: '';
        position: absolute;
        top: 0;
        left: 0;
        width: 6px;
        height: 100%;
        background: linear-gradient(180deg, #667eea 0%, #764ba2 100%);
    }
    
    .info-title {
        color: #2c3e50;
        font-size: 1.3rem;
        font-weight: bold;
        margin-bottom: 18px;
        display: flex;
        align-items: center;
        gap: 10px;
        text-shadow: 0 2px 4px rgba(0,0,0,0.1);
    }
    
    .info-content {
        color: #34495e;
        line-height: 1.7;
        font-size: 1rem;
    }
    
    .feature-list {
        list-style: none;
        padding: 0;
        margin: 15px 0;
    }
    
    .feature-item {
        padding: 12px 0;
        border-bottom: 1px solid rgba(102, 126, 234, 0.15);
        display: flex;
        align-items: center;
        gap: 12px;
        color: #2c3e50;
        transition: all 0.4s ease;
        font-weight: 500;
        border-radius: 8px;
        margin: 2px 0;
        padding-left: 8px;
    }
    
    .feature-item:hover {
        color: #667eea;
        transform: translateX(8px);
        background: rgba(102, 126, 234, 0.05);
        padding-left: 16px;
    }
    
    .feature-item:last-child {
        border-bottom: none;
    }
    
    /* File format section - MIGLIORATA */
    .format-section {
        background: linear-gradient(135deg, #e8f5e8 0%, #d4edda 100%);
        padding: 25px;
        border-radius: 20px;
        margin: 25px 0;
        box-shadow: 0 8px 25px rgba(0, 0, 0, 0.12);
        border-left: 6px solid #28a745;
        border: 2px solid rgba(40, 167, 69, 0.1);
        position: relative;
        overflow: hidden;
    }
    
    .format-section::before {
        content: '';
        position: absolute;
        top: 0;
        left: 0;
        width: 6px;
        height: 100%;
        background: linear-gradient(180deg, #28a745 0%, #20c997 100%);
    }
    
    .format-title {
        color: #155724;
        font-size: 1.3rem;
        font-weight: bold;
        margin-bottom: 18px;
        display: flex;
        align-items: center;
        gap: 10px;
        text-shadow: 0 2px 4px rgba(0,0,0,0.1);
    }
    
    .format-code {
        background: #ffffff;
        padding: 15px;
        border-radius: 8px;
        font-family: 'Monaco', 'Consolas', monospace;
        font-size: 0.85rem;
        border: 1px solid #c3e6cb;
        color: #155724;
    }
    
    .format-description {
        color: #155724;
        font-size: 0.9rem;
        margin-top: 10px;
        font-style: italic;
    }
    
    /* Responsive improvements */
    @media (max-width: 768px) {
        .nav-section, .info-section, .format-section, .generator-section {
            margin: 15px 0;
            padding: 15px;
        }
        
        .nav-title, .info-title, .format-title, .generator-title {
            font-size: 1.1rem;
        }
    }
    
    /* Scrollbar styling */
    .css-1d391kg::-webkit-scrollbar {
        width: 8px;
    }
    
    .css-1d391kg::-webkit-scrollbar-track {
        background: #f1f1f1;
        border-radius: 10px;
    }
    
    .css-1d391kg::-webkit-scrollbar-thumb {
        background: linear-gradient(180deg, #667eea 0%, #764ba2 100%);
        border-radius: 10px;
    }
    
    .css-1d391kg::-webkit-scrollbar-thumb:hover {
        background: linear-gradient(180deg, #5a6fd8 0%, #6a4190 100%);
    }

    /* Responsive figure containers */
    [data-testid="stImage"], .stPlotlyChart {
        max-width: 100%;
        overflow-x: auto;
    }
    
    /* Responsive metric cards */
    .metric-card h3 {
        font-size: clamp(1.2rem, 2.5vw, 2rem);
    }
    .metric-card p {
        font-size: clamp(0.75rem, 1.5vw, 0.95rem);
    }
    
    /* Responsive dataframes */
    [data-testid="stDataFrame"] {
        max-width: 100%;
        overflow-x: auto;
    }
    
    /* Ensure pyplot figures fill container width */
    [data-testid="stImage"] img {
        width: 100%;
        height: auto;
    }

    /* ── Plotly modebar: always visible ── */
    .modebar-container {
        opacity: 1 !important;
    }
    .modebar-group {
        display: flex !important;
    }
</style>
""", unsafe_allow_html=True)

# Main header
st.markdown("""
<div class="main-header">
    <h1>🧬 STAMP - Gene Switching Explorer</h1>
    <p>Spatio-Temporal Migration of gene Patterns in chronic Pathologies</p>
</div>
""", unsafe_allow_html=True)

# ── GTEx version selector (persists across all pages) ──────────
# The key="gtex_version" automatically syncs to st.session_state
st.sidebar.selectbox(
    "🧬 GTEx version",
    ["v10", "v8"],
    key="gtex_version",
)

# ── Tissue set: complete age-bin coverage only ──────────────────
# Only tissues with a sample in all six GTEx age brackets are used
# (the "all tissues" mode has been removed by design). Reads from
# output/{version}_complete/.
st.session_state["tissue_mode"] = "complete"
st.sidebar.caption(
    "🧪 Tissue set: **complete age bins only** — tissues with samples in all "
    "six GTEx age brackets (20-29 … 70-79)."
)

st.sidebar.markdown("---")

# Bottone Home in cima alla sidebar
home_active = st.session_state.current_page == "welcome"
if st.sidebar.button(
    "🏠 Home",
    key="nav_welcome",
    use_container_width=True,
    type="primary" if home_active else "secondary"
):
    if st.session_state.current_page != "welcome":
        st.session_state.current_page = "welcome"
        st.rerun()

# SEZIONE GENERATOR - PRIMA (ordine temporale)
st.sidebar.markdown("""
<div class="generator-section">
    <div class="generator-title">
        🛠️ Data Generation
    </div>
    <div class="generator-subtitle">
        Create STAMP datasets from raw data
    </div>
</div>
""", unsafe_allow_html=True)

# Generator navigation
generator_active = st.session_state.current_page == "stamp_generator"
if st.sidebar.button(
    "🛠️ STAMP Generator", 
    use_container_width=True,
    type="primary" if generator_active else "secondary"
):
    st.session_state.current_page = "stamp_generator"
    st.rerun()

# SEZIONE ANALISI - DOPO
st.sidebar.markdown("""
<div class="nav-section">
    <div class="nav-title">
        📊 Data Analysis
    </div>
    <div class="nav-subtitle">
        Analyze existing STAMP datasets
    </div>
</div>
""", unsafe_allow_html=True)

# Analysis pages
analysis_pages = [
    ("📤 Upload & Single Analysis", "upload_analysis"),
    ("🔄 Tissue Comparison", "tissue_comparison"), 
    ("🧬 Multi-Tissue Analysis", "multi_tissue"),
    ("📅 Age-Specific Analysis", "age_specific"),
    ("🤝 Gene Sharing Analysis", "gene_sharing"),
    ("🔍 Single Gene Analysis", "single_gene"),
    ("👥 Group Comparison", "group_comparison")
]

# Pipeline Analysis pages (new v8/v10 pages)
pipeline_pages = [
    ("🗺️ Aging Timeline", "aging_timeline"),
    ("🔼🔽 Switch Direction", "switch_direction"),
    ("📊 Threshold Sweep", "threshold_sweep"),
    ("🔀 Migration Paths", "migration_page"),
]

# Crea la navigazione con bottoni eleganti
for name, key in analysis_pages:
    is_active = st.session_state.current_page == key
    
    # Estrai emoji e nome pulito
    emoji = name.split()[0]
    clean_name = name.replace(emoji + " ", "")
    
    # Crea bottone con styling condizionale
    if st.sidebar.button(
        f"{emoji} {clean_name}",
        key=f"nav_{key}",
        use_container_width=True,
        type="primary" if is_active else "secondary"
    ):
        if st.session_state.current_page != key:
            st.session_state.current_page = key
            st.rerun()

# ── PIPELINE ANALYSIS SECTION (hidden behind a toggle) ──────────
# These extra, pre-computed pages are not shown by default: tick the
# checkbox to reveal them and choose one.
show_pipeline = st.sidebar.checkbox(
    "📈 Show Pipeline Analysis pages",
    value=False,
    key="show_pipeline_pages",
    help="Pre-computed pages: Aging Timeline, Switch Direction, Threshold "
         "Sweep, Migration Paths, v8 vs v10.",
)
if show_pipeline:
    st.sidebar.markdown("""
    <div class="nav-section" style="background: linear-gradient(135deg, #00b09b 0%, #96c93d 100%); box-shadow: 0 12px 35px rgba(0, 176, 155, 0.4);">
        <div class="nav-title">
            📈 Pipeline Analysis
        </div>
        <div class="nav-subtitle">
            Pre-computed v8/v10 results
        </div>
    </div>
    """, unsafe_allow_html=True)

    for name, key in pipeline_pages:
        is_active = st.session_state.current_page == key
        emoji = name.split()[0]
        clean_name = name.replace(emoji + " ", "")
        if st.sidebar.button(
            f"{emoji} {clean_name}",
            key=f"nav_{key}",
            use_container_width=True,
            type="primary" if is_active else "secondary"
        ):
            if st.session_state.current_page != key:
                st.session_state.current_page = key
                st.rerun()

# About STAMP Section
st.sidebar.markdown("""
<div class="info-section">
    <div class="info-title">
        ℹ️ About STAMP
    </div>
    <div class="info-content">
        <strong>Gene Switching Analysis Tool</strong>
        <br><br>
        Analyze gene expression patterns across:
        <ul class="feature-list">
            <li class="feature-item">🎯 Age groups (30-79 years)</li>
            <li class="feature-item">🧬 Different tissues</li>
            <li class="feature-item">📊 Statistical comparisons</li>
            <li class="feature-item">📈 Hierarchical clustering</li>
        </ul>
    </div>
</div>
""", unsafe_allow_html=True)

# File Format Section (solo per le pagine di analisi)
if st.session_state.current_page != "stamp_generator":
    st.sidebar.markdown("""
    <div class="format-section">
        <div class="format-title">
            📋 File Format
        </div>
        <div class="info-content">
            Upload <code>.txt</code> files with:
            <ul class="feature-list">
                <li class="feature-item">📝 5 lines (one per age group)</li>
                <li class="feature-item">🔤 Space-separated gene names</li>
                <li class="feature-item">📊 Format: gene1 gene2 gene3...</li>
            </ul>
        </div>
    </div>
    """, unsafe_allow_html=True)
    
    # Example code block separato per evitare problemi di rendering
    st.sidebar.code("""APOE TP53 BRCA1 EGFR
MYC PTEN RB1 VHL
APC KRAS PIK3CA IDH1
CDKN2A ATM SMAD4
MLH1 MSH2 MSH6 PMS2""", language="text")
    
    st.sidebar.caption("📝 Example: 5 age groups with gene names")
else:
    # Info section specifica per il generator
    st.sidebar.markdown("""
    <div class="format-section">
        <div class="format-title">
            🛠️ Generator Info
        </div>
        <div class="info-content">
            The STAMP Generator creates analysis-ready files from:
            <ul class="feature-list">
                <li class="feature-item">📊 Raw TPM expression data</li>
                <li class="feature-item">🧬 Gene ID mapping files</li>
                <li class="feature-item">⚙️ Configurable thresholds</li>
                <li class="feature-item">📁 Batch processing</li>
            </ul>
        </div>
    </div>
    """, unsafe_allow_html=True)
    
    st.sidebar.markdown("**Output:** Ready-to-use STAMP files for analysis", 
                       help="Generated files are compatible with all analysis modules")

# Current Page Indicator
if st.session_state.current_page == "welcome":
    current_emoji = "🏠"
    current_display_name = "Home"
elif st.session_state.current_page == "stamp_generator":
    current_emoji = "🛠️"
    current_display_name = "STAMP Generator"
else:
    # Search both analysis and pipeline pages
    all_nav_pages = analysis_pages + pipeline_pages
    current_emoji = next((name.split()[0] for name, key in all_nav_pages if key == st.session_state.current_page), "📤")
    current_display_name = next((name for name, key in all_nav_pages if key == st.session_state.current_page), "Upload Analysis")
    # Pulisci il nome per il display (remove leading emoji)
    import re as _re
    current_display_name = _re.sub(r'^\S+\s+', '', current_display_name)

st.sidebar.markdown(f"""
<div style="text-align: center; padding: 15px; background: linear-gradient(90deg, #667eea, #764ba2); 
     color: white; border-radius: 10px; margin: 20px 0; box-shadow: 0 4px 15px rgba(102, 126, 234, 0.3);">
    <div style="font-size: 1.2rem; margin-bottom: 5px;">{current_emoji}</div>
    <div style="font-weight: bold;">Current: {current_display_name}</div>
</div>
""", unsafe_allow_html=True)

# Main content area
st.markdown("---")

# Route to selected page based on session state
page_key = st.session_state.current_page

if page_key == "welcome":
    welcome.show()
elif page_key == "stamp_generator":
    stamp_generator.show()
elif page_key == "upload_analysis":
    upload_analysis.show()
elif page_key == "tissue_comparison":
    tissue_comparison.show()
elif page_key == "multi_tissue":
    multi_tissue.show()
elif page_key == "age_specific":
    age_specific.show()
elif page_key == "gene_sharing":
    gene_sharing.show()
elif page_key == "single_gene":
    single_gene.show()
elif page_key == "group_comparison":
    group_comparison.show()
elif page_key == "aging_timeline":
    aging_timeline.show()
elif page_key == "switch_direction":
    switch_direction.show()
elif page_key == "threshold_sweep":
    threshold_sweep.show()
elif page_key == "migration_page":
    migration_page.show()
else:
    # Fallback alla pagina upload se non riconosce la pagina
    upload_analysis.show()

# Footer
st.markdown("---")
st.markdown("""
<div style="text-align: center; color: #666; padding: 1rem;">
    <p>🧬 STAMP - Gene Switching Explorer | Built with Streamlit</p>
</div>
""", unsafe_allow_html=True)