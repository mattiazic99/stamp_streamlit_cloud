"""Shared stylesheet, injected by gui.main on each rerun.

Analysis pages use the analysis-section, metric-card and download-section classes
when rendering HTML through st.markdown.
"""
import streamlit as st

# Brand palette. The two gradients carry the whole interface: indigo for
# analysis, amber for data generation.
_INDIGO, _VIOLET = "#667eea", "#764ba2"
_AMBER_DARK, _AMBER = "#ff6b35", "#f7931e"
_GREEN, _TEAL = "#28a745", "#20c997"

_CSS = f"""
<style>
    /* ── Sidebar shell ─────────────────────────────────────────────── */
    .css-1d391kg {{
        background: linear-gradient(180deg, #f8f9fa 0%, #e9ecef 100%);
        border-right: 3px solid {_INDIGO};
    }}
    .css-1d391kg::-webkit-scrollbar {{
        width: 8px;
    }}
    .css-1d391kg::-webkit-scrollbar-track {{
        background: #f1f1f1;
        border-radius: 10px;
    }}
    .css-1d391kg::-webkit-scrollbar-thumb {{
        background: linear-gradient(180deg, {_INDIGO} 0%, {_VIOLET} 100%);
        border-radius: 10px;
    }}
    .css-1d391kg::-webkit-scrollbar-thumb:hover {{
        background: linear-gradient(180deg, #5a6fd8 0%, #6a4190 100%);
    }}

    /* ── Page header ───────────────────────────────────────────────── */
    .main-header {{
        background: linear-gradient(90deg, {_INDIGO} 0%, {_VIOLET} 100%);
        padding: 1.5rem;
        border-radius: 10px;
        color: white;
        text-align: center;
        margin-bottom: 2rem;
        box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
    }}
    .main-header h1 {{
        margin: 0;
        font-size: 2.5rem;
        text-shadow: 2px 2px 4px rgba(0, 0, 0, 0.3);
    }}
    .main-header p {{
        margin: 0.5rem 0 0 0;
        font-size: 1.1rem;
        opacity: 0.9;
    }}

    /* ── Blocks emitted by the analysis pages ──────────────────────── */
    .analysis-section {{
        background: #f8f9fa;
        padding: 1.5rem;
        border-radius: 10px;
        margin: 1rem 0;
        border-left: 5px solid {_INDIGO};
        box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
    }}
    .metric-card {{
        background: white;
        padding: 1rem;
        border-radius: 8px;
        box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
        text-align: center;
        margin: 0.5rem 0;
    }}
    .metric-card h3 {{
        font-size: clamp(1.2rem, 2.5vw, 2rem);
    }}
    .metric-card p {{
        font-size: clamp(0.75rem, 1.5vw, 0.95rem);
    }}
    .download-section {{
        background: #e8f5e8;
        padding: 1rem;
        border-radius: 8px;
        margin: 1rem 0;
        border: 1px solid #c3e6c3;
    }}
    .stSelectbox > div > div {{
        background-color: white;
    }}

    /* ── Sidebar section banners ───────────────────────────────────── */
    .nav-section, .generator-section {{
        padding: 25px;
        border-radius: 20px;
        margin: 25px 0;
        border: 2px solid rgba(255, 255, 255, 0.1);
        backdrop-filter: blur(10px);
        position: relative;
        overflow: hidden;
    }}
    .nav-section {{
        background: linear-gradient(135deg, {_INDIGO} 0%, {_VIOLET} 100%);
        box-shadow: 0 12px 35px rgba(102, 126, 234, 0.4);
    }}
    .generator-section {{
        background: linear-gradient(135deg, {_AMBER_DARK} 0%, {_AMBER} 100%);
        box-shadow: 0 12px 35px rgba(255, 107, 53, 0.4);
    }}
    .nav-section::before, .generator-section::before {{
        content: '';
        position: absolute;
        top: 0;
        left: 0;
        right: 0;
        height: 3px;
        background: linear-gradient(90deg,
            rgba(255, 255, 255, 0.8) 0%, rgba(255, 255, 255, 0.3) 100%);
    }}
    .nav-title, .generator-title {{
        color: white;
        font-size: 1.5rem;
        font-weight: bold;
        margin-bottom: 18px;
        text-align: center;
        text-shadow: 2px 2px 6px rgba(0, 0, 0, 0.4);
        display: flex;
        align-items: center;
        justify-content: center;
        gap: 12px;
        letter-spacing: 0.5px;
    }}
    .nav-subtitle, .generator-subtitle {{
        color: rgba(255, 255, 255, 0.9);
        font-size: 1rem;
        text-align: center;
        margin-bottom: 25px;
        font-style: italic;
        font-weight: 300;
        text-shadow: 1px 1px 3px rgba(0, 0, 0, 0.3);
    }}

    /* ── Navigation buttons ────────────────────────────────────────── */
    .stButton > button {{
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
    }}
    .stButton > button:hover {{
        transform: translateY(-2px) !important;
        box-shadow: 0 8px 25px rgba(0, 0, 0, 0.15) !important;
    }}
    /* Sheen that sweeps across the button on hover. */
    .stButton > button::before {{
        content: '';
        position: absolute;
        top: 0;
        left: -100%;
        width: 100%;
        height: 100%;
        background: linear-gradient(90deg,
            transparent, rgba(255, 255, 255, 0.2), transparent);
        transition: left 0.5s;
    }}
    .stButton > button:hover::before {{
        left: 100%;
    }}
    /* secondary = inactive page, primary = the page being displayed */
    .stButton > button[kind="secondary"] {{
        background: rgba(255, 255, 255, 0.9) !important;
        color: {_INDIGO} !important;
        border: 2px solid rgba(102, 126, 234, 0.2) !important;
    }}
    .stButton > button[kind="secondary"]:hover {{
        background: rgba(102, 126, 234, 0.1) !important;
        border: 2px solid {_INDIGO} !important;
    }}
    .stButton > button[kind="primary"] {{
        background: linear-gradient(135deg, {_INDIGO} 0%, {_VIOLET} 100%) !important;
        color: white !important;
        border: 2px solid rgba(255, 255, 255, 0.3) !important;
        box-shadow: 0 6px 20px rgba(102, 126, 234, 0.4) !important;
    }}
    .stButton > button[kind="primary"]:hover {{
        background: linear-gradient(135deg, #5a6fd8 0%, #6a4190 100%) !important;
        box-shadow: 0 8px 30px rgba(102, 126, 234, 0.5) !important;
    }}

    /* ── Sidebar information cards ─────────────────────────────────── */
    .info-section, .format-section {{
        padding: 25px;
        border-radius: 20px;
        margin: 25px 0;
        box-shadow: 0 8px 25px rgba(0, 0, 0, 0.12);
        position: relative;
        overflow: hidden;
    }}
    .info-section {{
        background: linear-gradient(135deg, #ffffff 0%, #f8f9fa 100%);
        border: 2px solid rgba(102, 126, 234, 0.1);
        border-left: 6px solid {_INDIGO};
    }}
    .format-section {{
        background: linear-gradient(135deg, #e8f5e8 0%, #d4edda 100%);
        border: 2px solid rgba(40, 167, 69, 0.1);
        border-left: 6px solid {_GREEN};
    }}
    .info-section::before, .format-section::before {{
        content: '';
        position: absolute;
        top: 0;
        left: 0;
        width: 6px;
        height: 100%;
    }}
    .info-section::before {{
        background: linear-gradient(180deg, {_INDIGO} 0%, {_VIOLET} 100%);
    }}
    .format-section::before {{
        background: linear-gradient(180deg, {_GREEN} 0%, {_TEAL} 100%);
    }}
    .info-title, .format-title {{
        font-size: 1.3rem;
        font-weight: bold;
        margin-bottom: 18px;
        display: flex;
        align-items: center;
        gap: 10px;
        text-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
    }}
    .info-title {{
        color: #2c3e50;
    }}
    .format-title {{
        color: #155724;
    }}
    .info-content {{
        color: #34495e;
        line-height: 1.7;
        font-size: 1rem;
    }}
    .feature-list {{
        list-style: none;
        padding: 0;
        margin: 15px 0;
    }}
    .feature-item {{
        padding: 12px 0 12px 8px;
        border-bottom: 1px solid rgba(102, 126, 234, 0.15);
        display: flex;
        align-items: center;
        gap: 12px;
        color: #2c3e50;
        transition: all 0.4s ease;
        font-weight: 500;
        border-radius: 8px;
        margin: 2px 0;
    }}
    .feature-item:hover {{
        color: {_INDIGO};
        transform: translateX(8px);
        background: rgba(102, 126, 234, 0.05);
        padding-left: 16px;
    }}
    .feature-item:last-child {{
        border-bottom: none;
    }}

    /* ── Figures and tables fill the container and scroll if narrower ── */
    [data-testid="stImage"], .stPlotlyChart, [data-testid="stDataFrame"] {{
        max-width: 100%;
        overflow-x: auto;
    }}
    [data-testid="stImage"] img {{
        width: 100%;
        height: auto;
    }}
    /* Plotly toolbar: visible without hovering, so export is discoverable. */
    .modebar-container {{
        opacity: 1 !important;
    }}
    .modebar-group {{
        display: flex !important;
    }}

    @media (max-width: 768px) {{
        .nav-section, .generator-section, .info-section, .format-section {{
            margin: 15px 0;
            padding: 15px;
        }}
        .nav-title, .generator-title, .info-title, .format-title {{
            font-size: 1.1rem;
        }}
    }}
</style>
"""


def inject() -> None:
    """Write the stylesheet into the page."""
    st.markdown(_CSS, unsafe_allow_html=True)
