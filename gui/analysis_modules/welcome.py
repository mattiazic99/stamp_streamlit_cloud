"""Home page: generation, upload, analysis and the five-line STAMP file format."""
from pathlib import Path

import streamlit as st

# Palette, shared with the sidebar banners in components/theme.py
_INDIGO = "#667eea"
_CARD_BG = "rgba(102, 126, 234, 0.1)"
_CARD_BORDER = "rgba(102, 126, 234, 0.2)"

# The three steps, in the order they are performed. Each is (number, gradient,
# shadow colour, icon + title, body).
_STEPS = [
    (
        "1&#xFE0F;&#x20E3;",
        "#ff6b35 0%, #f7931e 100%",
        "rgba(255, 107, 53, 0.35)",
        "&#x1F6E0;&#xFE0F; Generate Files",
        "Use the <strong>STAMP Generator</strong> to create analysis files "
        "from raw expression data (TPM).<br><br>"
        "&#x26A0;&#xFE0F; <em>This step must be done <strong>at least once"
        "</strong> before you can analyze any data.</em>",
    ),
    (
        "2&#xFE0F;&#x20E3;",
        "#667eea 0%, #764ba2 100%",
        "rgba(102, 126, 234, 0.35)",
        "&#x1F4CA; Analyze Results",
        "Load the generated files into the <strong>analysis pages</strong> "
        "to explore gene switching patterns.<br><br>"
        "&#x1F504; <em>You can reuse the generated files <strong>as many "
        "times as you want</strong>.</em>",
    ),
    (
        "3&#xFE0F;&#x20E3;",
        "#28a745 0%, #20c997 100%",
        "rgba(40, 167, 69, 0.35)",
        "&#x1F4E5; Explore &amp; Download",
        "View <strong>heatmaps, charts, and tables</strong>. Then download "
        "your results as images or CSV files.<br><br>"
        "&#x1F4C8; <em>Each page offers <strong>interactive visualizations"
        "</strong> and downloads.</em>",
    ),
]

# One entry per page in the sidebar, in the same order.
_PAGES = [
    ("&#x1F6E0;&#xFE0F; STAMP Generator",
     "Generate STAMP files from raw TPM expression data"),
    ("&#x1F4E4; Upload &amp; Single Analysis",
     "Upload a file and analyze a single tissue"),
    ("&#x1F504; Tissue Comparison",
     "Compare gene switching across two tissues"),
    ("&#x1F9EC; Multi-Tissue Analysis",
     "Simultaneous analysis across multiple tissues"),
    ("&#x1F4C5; Age-Specific Analysis",
     "Analysis for specific age groups"),
    ("&#x1F91D; Gene Sharing Analysis",
     "Discover shared genes across different tissues"),
    ("&#x1F50D; Single Gene Analysis",
     "Detailed analysis of a single gene"),
    ("&#x1F465; Group Comparison",
     "Compare gene groups across tissues"),
    ("&#x1F9E9; Panel Explorer",
     "Map a gene panel onto the atlas, with v8/v10 concordance"),
]

# Five space-separated gene lists, one per switching bracket, in chronological
# order — what the generator writes and what the analysis pages read.
_EXAMPLE_FILE = (
    "APOE TP53 BRCA1 EGFR<br>"
    "MYC PTEN RB1 VHL<br>"
    "APC KRAS PIK3CA IDH1<br>"
    "CDKN2A ATM SMAD4<br>"
    "MLH1 MSH2 MSH6 PMS2"
)

_PANEL = (
    "background: linear-gradient(135deg, #ffffff 0%, #f8f9fa 100%);"
    "padding: 30px; border-radius: 20px; margin-bottom: 25px;"
    "box-shadow: 0 8px 25px rgba(0, 0, 0, 0.1);"
    "border: 2px solid rgba(102, 126, 234, 0.15);"
)


def _step_card(number: str, gradient: str, shadow: str, title: str, body: str) -> str:
    return (
        f'<div style="background: linear-gradient(135deg, {gradient});'
        'padding: 25px; border-radius: 18px; text-align: center; color: white;'
        f'min-height: 280px; box-shadow: 0 10px 30px {shadow};'
        'border: 2px solid rgba(255, 255, 255, 0.15);">'
        f'<div style="font-size: 3rem; margin-bottom: 10px;">{number}</div>'
        '<h3 style="margin: 0 0 12px 0; font-size: 1.3rem;'
        f' text-shadow: 1px 1px 3px rgba(0,0,0,0.2);">{title}</h3>'
        '<div style="font-size: 0.95rem; line-height: 1.6; opacity: 0.95;">'
        f'{body}</div></div>'
    )


def _page_card(title: str, description: str) -> str:
    return (
        f'<div style="background: {_CARD_BG}; padding: 18px;'
        f' border-radius: 14px; border: 1px solid {_CARD_BORDER};">'
        f'<strong style="color: {_INDIGO}; font-size: 1.1rem;">{title}</strong><br>'
        f'<span style="font-size: 0.9rem; color: #555;">{description}</span>'
        '</div>'
    )


def show():
    """Welcome page — explains how to use the STAMP application."""
    st.markdown(
        '<div style="background: linear-gradient(135deg, #11998e 0%, #38ef7d 100%);'
        'padding: 40px; border-radius: 20px; text-align: center; color: white;'
        'margin-bottom: 30px; box-shadow: 0 12px 35px rgba(17, 153, 142, 0.4);">'
        '<h1 style="margin: 0; font-size: 2.5rem;'
        ' text-shadow: 2px 2px 6px rgba(0,0,0,0.3);">'
        '&#x1F44B; Welcome to STAMP</h1>'
        '<div style="font-size: 1.2rem; margin-top: 10px; opacity: 0.95;">'
        'Gene Switching Explorer &#8212; Analyze gene expression patterns '
        'in chronic pathologies</div></div>',
        unsafe_allow_html=True,
    )

    st.markdown(
        f'<div style="{_PANEL}">'
        '<h2 style="color: #2c3e50; text-align: center; margin-bottom: 15px;'
        ' font-size: 1.8rem;">&#x1F680; How does it work?</h2>'
        '<div style="color: #555; text-align: center; font-size: 1.05rem;'
        ' margin-bottom: 10px; line-height: 1.6;">'
        'STAMP works in <strong>3 simple steps</strong>. First, generate your '
        'analysis files, then you can use them as many times as you want '
        'across the different analysis pages.</div></div>',
        unsafe_allow_html=True,
    )

    for column, step in zip(st.columns(len(_STEPS)), _STEPS):
        with column:
            st.markdown(_step_card(*step), unsafe_allow_html=True)

    st.markdown("<br>", unsafe_allow_html=True)

    st.markdown(
        '<div style="background: linear-gradient(135deg, #fff3cd 0%, #ffeeba 100%);'
        'padding: 30px; border-radius: 20px; margin-bottom: 25px;'
        'box-shadow: 0 8px 25px rgba(0, 0, 0, 0.1);'
        'border: 2px solid rgba(255, 193, 7, 0.3);">'
        '<h2 style="color: #856404; text-align: center; margin-bottom: 20px;'
        ' font-size: 1.8rem;">&#x1F4C4; What are STAMP files?</h2>'
        '<div style="color: #856404; font-size: 1.05rem; line-height: 1.7;'
        ' margin-bottom: 20px; text-align: center;">'
        'The analysis pages <strong>only accept files in STAMP format</strong>. '
        'These are <code>.txt</code> files with a specific structure, generated '
        'automatically by the <strong>STAMP Generator</strong>.</div>'
        '<div style="display: grid; grid-template-columns: 1fr 1fr; gap: 20px;">'
        '<div style="background: white; padding: 20px; border-radius: 14px;'
        ' border: 1px solid rgba(133, 100, 4, 0.15);">'
        '<div style="color: #856404; font-weight: bold; font-size: 1.1rem;'
        ' margin-bottom: 12px;">&#x1F4CB; File Structure</div>'
        '<div style="color: #555; font-size: 0.95rem; line-height: 1.7;">'
        '&#x2022; Each file is a <code>.txt</code> file representing '
        '<strong>one tissue</strong><br>'
        '&#x2022; Contains <strong>5 lines</strong>, one for each age group '
        '(30-39, 40-49, 50-59, 60-69, 70-79)<br>'
        '&#x2022; Each line contains <strong>space-separated gene names</strong> '
        'that switch in that age group<br>'
        '&#x2022; These files are the output of the STAMP Generator</div></div>'
        '<div style="background: white; padding: 20px; border-radius: 14px;'
        ' border: 1px solid rgba(133, 100, 4, 0.15);">'
        '<div style="color: #856404; font-weight: bold; font-size: 1.1rem;'
        ' margin-bottom: 12px;">&#x1F4A1; Example (tissue_brain.txt)</div>'
        '<div style="background: #f8f9fa; padding: 12px; border-radius: 8px;'
        ' font-size: 0.85rem; color: #155724; border: 1px solid #c3e6cb;'
        ' font-family: Monaco, Consolas, monospace; line-height: 1.8;">'
        f'{_EXAMPLE_FILE}</div>'
        '<div style="color: #888; font-size: 0.8rem; margin-top: 8px;'
        ' font-style: italic;">'
        'Line 1 = Age 30-39, Line 2 = Age 40-49, ... Line 5 = Age 70-79'
        '</div></div></div></div>',
        unsafe_allow_html=True,
    )

    st.markdown(
        f'<div style="{_PANEL}">'
        '<h2 style="color: #2c3e50; text-align: center; margin-bottom: 25px;'
        ' font-size: 1.8rem;">&#x1F4CB; Available Pages</h2>'
        '<div style="display: grid; grid-template-columns: 1fr 1fr; gap: 15px;">'
        + "".join(_page_card(title, desc) for title, desc in _PAGES)
        + '</div></div>',
        unsafe_allow_html=True,
    )

    st.markdown(
        '<div style="background: linear-gradient(135deg, #e8f5e8 0%, #d4edda 100%);'
        'padding: 25px 30px; border-radius: 18px; margin: 25px 0;'
        'border: 2px solid rgba(40, 167, 69, 0.3);'
        'box-shadow: 0 6px 20px rgba(40, 167, 69, 0.15); text-align: center;">'
        '<div style="font-size: 1.5rem; margin-bottom: 8px;">&#x1F4A1;</div>'
        '<div style="color: #155724; font-size: 1.15rem; font-weight: bold;'
        ' margin-bottom: 8px;">Already generated your STAMP files?</div>'
        '<div style="color: #155724; font-size: 1.05rem; line-height: 1.6;">'
        'If you have already run the STAMP Generator at least once, your files '
        'are ready! You can skip this step and go directly to any '
        '<strong>analysis page</strong> using the navigation menu on the left.'
        '</div></div>',
        unsafe_allow_html=True,
    )

    st.markdown("<br>", unsafe_allow_html=True)

    _left, centre, _right = st.columns([1, 2, 1])
    with centre:
        if st.button(
            "\U0001F680 Get Started — Go to STAMP Generator",
            use_container_width=True,
            type="primary",
        ):
            st.session_state.current_page = "stamp_generator"
            st.rerun()

    with st.expander("Reproduce the statistical validation"):
        st.write(
            "Reproduce the statistical results reported in the paper using "
            "the command-line tools, reference results and checked GTEx v10 inputs. "
            "The guide covers a quick Ovary replay and the complete experiments. "
            "Run these separately from the interface, without Docker, "
            "using the full repository and its dedicated Python environment."
        )
        st.markdown("[Full repository and command-line tools](https://github.com/mattiazic99/stamp_streamlit_cloud)")
        project_root = Path(__file__).resolve().parent.parent.parent
        guide = project_root / "reproducibility" / "README.md"
        lock = project_root / "requirements-reproduction.lock.txt"
        if guide.is_file():
            st.download_button(
                "Download reproduction guide",
                data=guide.read_bytes(),
                file_name="STAMP_statistical_reproduction.md",
                mime="text/markdown",
                key="download_reproduction_guide",
            )
        if lock.is_file():
            st.download_button(
                "Download dependency versions",
                data=lock.read_bytes(),
                file_name="requirements-reproduction.lock.txt",
                mime="text/plain",
                key="download_reproduction_lock",
            )

