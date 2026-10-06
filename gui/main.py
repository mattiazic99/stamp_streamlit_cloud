"""Streamlit entry point: page setup, sidebar and navigation.

Each page in analysis_modules exposes a show() function. Analysis pages accept
five-line STAMP files; Panel Explorer also reads the bundled switching atlas.

Run with: streamlit run gui/main.py
"""
import streamlit as st

from analysis_modules import (
    age_specific,
    gene_sharing,
    group_comparison,
    multi_tissue,
    panel_explorer,
    single_gene,
    stamp_generator,
    tissue_comparison,
    upload_analysis,
    welcome,
)
from components import theme
from data_loader import GTEX_VERSION

# Only supported public pages are registered here.
# Pages
# (session key, sidebar label, module), in sidebar order.
HOME_PAGE = ("welcome", "🏠 Home", welcome)
GENERATOR_PAGE = ("stamp_generator", "🛠️ STAMP Generator", stamp_generator)
ANALYSIS_PAGES = [
    ("upload_analysis", "📤 Upload & Single Analysis", upload_analysis),
    ("tissue_comparison", "🔄 Tissue Comparison", tissue_comparison),
    ("multi_tissue", "🧬 Multi-Tissue Analysis", multi_tissue),
    ("age_specific", "📅 Age-Specific Analysis", age_specific),
    ("gene_sharing", "🤝 Gene Sharing Analysis", gene_sharing),
    ("single_gene", "🔍 Single Gene Analysis", single_gene),
    ("group_comparison", "👥 Group Comparison", group_comparison),
    ("panel_explorer", "🧩 Panel Explorer", panel_explorer),
]

ALL_PAGES = [HOME_PAGE, GENERATOR_PAGE] + ANALYSIS_PAGES
PAGE_MODULES = {key: module for key, _, module in ALL_PAGES}
PAGE_LABELS = {key: label for key, label, _ in ALL_PAGES}

DEFAULT_PAGE = HOME_PAGE[0]

# Example of the exchange format, shown in the sidebar of every analysis page:
# five space-separated gene lists, one per switching bracket, in chronological
# order (30-39, 40-49, 50-59, 60-69, 70-79).
EXAMPLE_STAMP_FILE = """APOE TP53 BRCA1 EGFR
MYC PTEN RB1 VHL
APC KRAS PIK3CA IDH1
CDKN2A ATM SMAD4
MLH1 MSH2 MSH6 PMS2"""


def navigate_to(page_key: str) -> None:
    """Switch to ``page_key`` and rerun, unless it is already displayed."""
    if st.session_state.current_page != page_key:
        st.session_state.current_page = page_key
        st.rerun()


def sidebar_nav_button(page_key: str, label: str) -> None:
    """Render one sidebar entry, highlighted when its page is displayed."""
    is_active = st.session_state.current_page == page_key
    if st.sidebar.button(
        label,
        key=f"nav_{page_key}",
        use_container_width=True,
        type="primary" if is_active else "secondary",
    ):
        navigate_to(page_key)


def render_sidebar() -> None:
    """Draw the persistent sidebar: dataset banner, navigation, help."""
    st.sidebar.caption(
        f"🧬 Dataset: **GTEx {GTEX_VERSION}**, age-complete tissues only — "
        "those with samples in all six GTEx age brackets (20-29 … 70-79)."
    )
    st.sidebar.markdown("---")

    sidebar_nav_button(*HOME_PAGE[:2])

    # Data generation comes first: the analysis pages need its output.
    st.sidebar.markdown(
        '<div class="generator-section">'
        '<div class="generator-title">🛠️ Data Generation</div>'
        '<div class="generator-subtitle">Create STAMP datasets from raw data</div>'
        '</div>',
        unsafe_allow_html=True,
    )
    sidebar_nav_button(*GENERATOR_PAGE[:2])

    st.sidebar.markdown(
        '<div class="nav-section">'
        '<div class="nav-title">📊 Data Analysis</div>'
        '<div class="nav-subtitle">Analyze existing STAMP datasets</div>'
        '</div>',
        unsafe_allow_html=True,
    )
    for page_key, label, _ in ANALYSIS_PAGES:
        sidebar_nav_button(page_key, label)

    st.sidebar.markdown(
        '<div class="info-section">'
        '<div class="info-title">ℹ️ About STAMP</div>'
        '<div class="info-content">'
        '<strong>Gene Switching Analysis Tool</strong><br><br>'
        'Analyze gene expression patterns across:'
        '<ul class="feature-list">'
        '<li class="feature-item">🎯 Age groups (30-79 years)</li>'
        '<li class="feature-item">🧬 Different tissues</li>'
        '<li class="feature-item">📊 Statistical comparisons</li>'
        '<li class="feature-item">📈 Hierarchical clustering</li>'
        '</ul>'
        '</div></div>',
        unsafe_allow_html=True,
    )

    if st.session_state.current_page == GENERATOR_PAGE[0]:
        render_generator_help()
    else:
        render_file_format_help()

    render_current_page_indicator()


def render_file_format_help() -> None:
    """Describe the STAMP exchange format accepted by the analysis pages."""
    st.sidebar.markdown(
        '<div class="format-section">'
        '<div class="format-title">📋 File Format</div>'
        '<div class="info-content">'
        'Upload <code>.txt</code> files with:'
        '<ul class="feature-list">'
        '<li class="feature-item">📝 5 lines (one per age group)</li>'
        '<li class="feature-item">🔤 Space-separated gene names</li>'
        '<li class="feature-item">📊 Format: gene1 gene2 gene3...</li>'
        '</ul>'
        '</div></div>',
        unsafe_allow_html=True,
    )
    # Kept out of the block above: st.code renders the example reliably,
    # while a <pre> inside the same markdown call does not.
    st.sidebar.code(EXAMPLE_STAMP_FILE, language="text")
    st.sidebar.caption("📝 Example: 5 age groups with gene names")


def render_generator_help() -> None:
    """Summarise what the generator consumes and produces."""
    st.sidebar.markdown(
        '<div class="format-section">'
        '<div class="format-title">🛠️ Generator Info</div>'
        '<div class="info-content">'
        'The STAMP Generator creates analysis-ready files from:'
        '<ul class="feature-list">'
        '<li class="feature-item">📊 Raw TPM expression data</li>'
        '<li class="feature-item">🧬 Gene ID mapping files</li>'
        '<li class="feature-item">⚙️ Configurable thresholds</li>'
        '<li class="feature-item">📁 Batch processing</li>'
        '</ul>'
        '</div></div>',
        unsafe_allow_html=True,
    )
    st.sidebar.markdown(
        "**Output:** Ready-to-use STAMP files for analysis",
        help="Generated files are compatible with all analysis modules",
    )


def render_current_page_indicator() -> None:
    """Show which page is displayed, at the foot of the sidebar."""
    label = PAGE_LABELS.get(st.session_state.current_page, PAGE_LABELS[DEFAULT_PAGE])
    emoji, _, name = label.partition(" ")
    st.sidebar.markdown(
        '<div style="text-align: center; padding: 15px; background: '
        'linear-gradient(90deg, #667eea, #764ba2); color: white; '
        'border-radius: 10px; margin: 20px 0; '
        'box-shadow: 0 4px 15px rgba(102, 126, 234, 0.3);">'
        f'<div style="font-size: 1.2rem; margin-bottom: 5px;">{emoji}</div>'
        f'<div style="font-weight: bold;">Current: {name}</div>'
        '</div>',
        unsafe_allow_html=True,
    )


def main() -> None:
    st.set_page_config(
        page_title="STAMP - Gene Switching Explorer",
        layout="wide",
        initial_sidebar_state="expanded",
    )
    st.session_state.setdefault("current_page", DEFAULT_PAGE)
    if st.session_state.current_page not in PAGE_MODULES:
        st.session_state.current_page = DEFAULT_PAGE

    theme.inject()
    st.markdown(
        '<div class="main-header">'
        '<h1>🧬 STAMP - Gene Switching Explorer</h1>'
        '<p>Spatio-Temporal Analysis and Mapping of gene-expression Patterns</p>'
        '</div>',
        unsafe_allow_html=True,
    )

    render_sidebar()

    st.markdown("---")
    page = PAGE_MODULES.get(st.session_state.current_page, PAGE_MODULES[DEFAULT_PAGE])
    page.show()

    st.markdown("---")
    st.markdown(
        '<div style="text-align: center; color: #666; padding: 1rem;">'
        '<p>🧬 STAMP - Gene Switching Explorer | Built with Streamlit</p>'
        '</div>',
        unsafe_allow_html=True,
    )


main()
