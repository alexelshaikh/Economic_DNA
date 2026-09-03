from __future__ import annotations

import pandas as pd
import streamlit as st

from economic_dna import (
    Scenario,
    find_crossover_years,
    simulate_dna_unit_costs,
    simulate_scenario,
    simulate_start_years,
)
from economic_dna.assumptions import load_assumptions
from economic_dna.visualization import (
    breakdown_exports,
    breakdown_chart,
    dna_unit_cost_chart,
    dna_unit_cost_exports,
    lifecycle_chart,
    lifecycle_exports,
    projection_chart,
    projection_exports,
)


st.set_page_config(
    page_title="DNA Storage Cost Explorer",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown(
    """
    <style>
    :root {
        --canvas: #f4f6f5;
        --surface: #ffffff;
        --surface-subtle: #eef2f0;
        --ink: #18211f;
        --muted: #64706c;
        --line: #d9dfdc;
        --line-strong: #c4cdc9;
        --accent: #bd4b38;
        --accent-strong: #963829;
        --accent-soft: #faece8;
        --teal: #1f716a;
        --gold: #b78322;
    }

    html, body {
        font-family: "Aptos", "Segoe UI", Arial, sans-serif;
        letter-spacing: 0;
    }
    .stApp {
        background: var(--canvas);
        color: var(--ink);
        font-family: "Aptos", "Segoe UI", Arial, sans-serif;
    }
    [data-testid="stHeader"] {
        background: rgba(244, 246, 245, 0.96);
        border-bottom: 1px solid rgba(217, 223, 220, 0.8);
    }
    .block-container {
        max-width: 1440px;
        padding: 4.75rem 2.2rem 4rem;
    }

    h1, h2, h3, p { letter-spacing: 0; }
    h1 {
        color: var(--ink);
        font-size: 2.15rem;
        line-height: 1.14;
        margin: 0.15rem 0 0.35rem;
    }
    h2 { color: var(--ink); font-size: 1.25rem; }
    h3 { color: var(--ink); font-size: 1.05rem; }
    [data-testid="stCaptionContainer"] {
        color: var(--muted);
        font-size: 0.86rem;
        line-height: 1.55;
    }

    section[data-testid="stSidebar"] {
        width: 460px !important;
        background: #fbfcfb;
        border-right: 1px solid var(--line);
    }
    section[data-testid="stSidebar"] > div { width: 460px !important; }
    section[data-testid="stSidebar"] [data-testid="stSidebarContent"] {
        overflow-x: hidden;
        overflow-y: hidden;
        padding: 1.375rem 0.7rem 0 1.2rem;
    }
    section[data-testid="stSidebar"] [data-testid="stSidebarHeader"] {
        margin-bottom: 0;
    }
    section[data-testid="stSidebar"] [data-testid="stSidebarUserContent"] {
        padding-bottom: 0;
    }
    .sidebar-kicker, .page-kicker, .workspace-kicker, .export-label {
        color: var(--accent);
        font-size: 0.72rem;
        font-weight: 750;
        line-height: 1.2;
        letter-spacing: 0;
        text-transform: uppercase;
    }
    .sidebar-title {
        color: var(--ink);
        font-size: 1.35rem;
        font-weight: 720;
        line-height: 1.25;
        margin: 0.22rem 0 0.35rem;
    }
    .sidebar-copy {
        color: var(--muted);
        font-size: 0.86rem;
        line-height: 1.5;
        margin: 0 0 0.9rem;
    }
    .sidebar-section {
        border-top: 1px solid var(--line);
        color: var(--ink);
        font-size: 0.77rem;
        font-weight: 720;
        letter-spacing: 0;
        margin: 1.1rem 0 0.65rem;
        padding-top: 0.9rem;
        text-transform: uppercase;
    }
    section[data-testid="stSidebar"] label p {
        color: #34413d;
        font-size: 0.82rem;
        font-weight: 600;
    }
    section[data-testid="stSidebar"] [data-baseweb="input"],
    section[data-testid="stSidebar"] [data-baseweb="select"] > div {
        background: var(--surface);
        border-color: var(--line-strong);
        border-radius: 6px;
    }
    section[data-testid="stSidebar"] [data-baseweb="input"]:focus-within,
    section[data-testid="stSidebar"] [data-baseweb="select"] > div:focus-within {
        border-color: var(--accent);
        box-shadow: 0 0 0 1px var(--accent);
    }
    section[data-testid="stSidebar"] [data-testid="stExpander"] {
        background: var(--surface);
        border: 1px solid var(--line);
        border-radius: 6px;
        margin-bottom: 0.55rem;
    }
    section[data-testid="stSidebar"] [data-testid="stExpander"] details summary {
        min-height: 2.75rem;
    }
    section[data-testid="stSidebar"] [data-testid="stVerticalBlockBorderWrapper"] {
        background: var(--surface);
        border-color: var(--line);
        border-radius: 6px;
    }

    .stButton button, .stDownloadButton button {
        border-radius: 6px;
        font-weight: 650;
        min-height: 2.35rem;
        transition: border-color 120ms ease, background-color 120ms ease, color 120ms ease;
    }
    .stButton button:hover, .stDownloadButton button:hover {
        border-color: var(--accent);
        color: var(--accent-strong);
    }
    [data-testid="stFormSubmitButton"] button {
        min-height: 2.75rem;
        background: var(--accent);
        border-color: var(--accent);
        color: #ffffff;
        font-weight: 700;
    }
    [data-testid="stFormSubmitButton"] button:hover {
        background: var(--accent-strong);
        border-color: var(--accent-strong);
        color: #ffffff;
    }
    /* Two-column sidebar: the left column holds all inputs and scrolls on its
       own; the right column holds only the Calculate button. The sidebar's
       own scroll container is disabled (position: sticky is unreliable inside
       it in this DOM: the absolutely-positioned flex app container breaks it),
       so instead the input column becomes the scroller and the button rail
       sits beside it as a fixed sibling that never moves. The in-flow sidebar
       header occupies the top 82px (22px content padding + 60px header in
       Streamlit 1.63), which is reserved so nothing overlaps the collapse
       control. */
    section[data-testid="stSidebar"] [data-testid="stHorizontalBlock"]:has(.st-key-scenario-action-rail) {
        align-items: stretch;
        gap: 0.7rem !important;
    }
    /* Note the child combinator: a descendant selector would also match the
       first column of every inner st.columns() row inside the input column,
       stretching each widget row to the full visible height. */
    section[data-testid="stSidebar"] [data-testid="stHorizontalBlock"]:has(.st-key-scenario-action-rail)
    > [data-testid="stColumn"]:first-child {
        height: calc(100vh - 82px);
        overflow-y: auto;
        overscroll-behavior-y: contain;
        padding: 0 0 1.5rem;
        scrollbar-gutter: stable;
    }
    section[data-testid="stSidebar"] [data-testid="stHorizontalBlock"]:has(.st-key-scenario-action-rail)
    > [data-testid="stColumn"]:first-child::-webkit-scrollbar {
        width: 8px;
    }
    section[data-testid="stSidebar"] [data-testid="stHorizontalBlock"]:has(.st-key-scenario-action-rail)
    > [data-testid="stColumn"]:first-child::-webkit-scrollbar-thumb {
        background: var(--line-strong);
        border-radius: 4px;
    }
    section[data-testid="stSidebar"] .st-key-scenario-action-rail {
        flex: 0 0 auto;
        height: calc(100vh - 82px);
        min-height: calc(100vh - 82px);
        overflow: visible;
        position: relative;
        z-index: 1;
    }
    /* Percentage heights cannot resolve here (the rail sits in a flex-content-
       sized column), so the button fills the rail with absolute positioning
       instead of a height: 100% chain. */
    section[data-testid="stSidebar"] .st-key-scenario-action-rail
    [data-testid="stElementContainer"] {
        position: absolute;
        inset: 0;
    }
    section[data-testid="stSidebar"] .st-key-scenario-action-rail .stButton button {
        align-items: center;
        background: var(--accent);
        border-color: var(--accent);
        border-radius: 6px;
        box-shadow: 0 3px 10px rgba(24, 33, 31, 0.14);
        color: #ffffff;
        display: flex;
        flex-direction: column;
        inset: 0;
        justify-content: space-between;
        padding: 1rem 0.45rem;
        position: absolute;
    }
    section[data-testid="stSidebar"] .st-key-scenario-action-rail .stButton button:hover {
        background: var(--accent-strong);
        border-color: var(--accent-strong);
        color: #ffffff;
    }
    section[data-testid="stSidebar"] [data-testid="stSidebarCollapseButton"],
    [data-testid="stExpandSidebarButton"],
    [data-testid="stSidebarCollapsedControl"] {
        pointer-events: auto !important;
        z-index: 5 !important;
    }
    /* When the sidebar is collapsed its content slides off-canvas but still
       overlaps the main area (the sidebar section keeps a high z-index), so
       the pinned button would paint over the main content and the header
       would intercept clicks aimed at the expand control. Hide the rail and
       let pointer events pass through to the main area. */
    section[data-testid="stSidebar"][aria-expanded="false"] .st-key-scenario-action-rail {
        visibility: hidden;
    }
    section[data-testid="stSidebar"][aria-expanded="false"] [data-testid="stSidebarContent"] {
        pointer-events: none;
    }
    section[data-testid="stSidebar"] .st-key-scenario-action-rail
    .stButton button::before,
    section[data-testid="stSidebar"] .st-key-scenario-action-rail
    .stButton button::after {
        content: "\\2192";
        display: block;
        flex: 0 0 auto;
        font-family: "Segoe UI Symbol", sans-serif;
        font-size: 1.55rem;
        line-height: 1;
    }
    section[data-testid="stSidebar"] .st-key-scenario-action-rail .stButton button p {
        font-size: 0.9rem;
        line-height: 1.25;
        white-space: nowrap;
        writing-mode: vertical-rl;
        transform: rotate(180deg);
    }
    .stDownloadButton button {
        background: var(--surface);
        border-color: var(--line);
        color: #40504b;
    }

    .page-kicker { margin-bottom: 0.35rem; }
    .page-deck {
        color: var(--muted);
        font-size: 1rem;
        line-height: 1.55;
        margin: 0 0 1.15rem;
        max-width: 760px;
    }
    .model-strip {
        align-items: center;
        background: var(--surface);
        border: 1px solid var(--line);
        border-left: 3px solid var(--accent);
        border-radius: 6px;
        color: var(--muted);
        display: flex;
        flex-wrap: wrap;
        font-size: 0.8rem;
        gap: 0.35rem 0;
        margin: 0 0 1.25rem;
        min-height: 42px;
        padding: 0.65rem 0.85rem;
    }
    .model-strip strong { color: var(--ink); font-weight: 700; }
    .model-item {
        border-right: 1px solid var(--line);
        margin-right: 0.75rem;
        padding-right: 0.75rem;
    }
    .model-item:last-child { border-right: 0; margin-right: 0; padding-right: 0; }

    .scenario-bar {
        align-items: center;
        background: var(--surface-subtle);
        border: 1px solid var(--line);
        border-radius: 6px;
        color: #475650;
        display: flex;
        flex-wrap: wrap;
        font-size: 0.82rem;
        gap: 0.45rem 1.15rem;
        margin: 0 0 1rem;
        min-height: 44px;
        padding: 0.65rem 0.85rem;
    }
    .scenario-bar strong { color: var(--ink); }
    .scenario-label { color: var(--teal); font-weight: 750; text-transform: uppercase; }
    .scenario-pending {
        background: #f6ecd8;
        border: 1px solid #e3cf9e;
        border-radius: 4px;
        color: #8a6414;
        font-weight: 650;
        padding: 0.1rem 0.55rem;
    }

    [data-testid="stMetric"] {
        background: var(--surface);
        border: 1px solid var(--line);
        border-radius: 8px;
        box-shadow: 0 1px 2px rgba(24, 33, 31, 0.04);
        min-height: 114px;
        overflow: hidden;
        padding: 0.95rem 1rem;
        position: relative;
    }
    [data-testid="stMetric"]::before {
        background: var(--accent);
        content: "";
        height: 3px;
        left: 0;
        position: absolute;
        right: 0;
        top: 0;
    }
    [data-testid="stMetricLabel"] p {
        color: var(--muted);
        font-size: 0.79rem;
        font-weight: 650;
        line-height: 1.35;
        white-space: normal;
    }
    [data-testid="stMetricValue"] {
        color: var(--ink);
        font-size: 1.48rem;
        line-height: 1.2;
    }
    [data-testid="stMetricDelta"] { color: var(--muted); }

    .workspace-kicker { margin: 1.55rem 0 0.35rem; }
    div[data-testid="stTabs"] div[role="tablist"] {
        background: transparent;
        border-bottom: 1px solid var(--line-strong);
        border-radius: 0;
        gap: 1.6rem;
        margin: 0 0 1.2rem;
        overflow-x: auto;
        padding: 0;
    }
    div[data-testid="stTabs"] div[data-testid="stTab"][role="tab"] {
        background: transparent;
        border: 0;
        border-bottom: 2px solid transparent;
        border-radius: 0;
        color: var(--muted);
        flex: 0 0 auto;
        font-weight: 650;
        min-height: 44px;
        padding: 0.6rem 0.1rem 0.7rem;
    }
    div[data-testid="stTabs"] div[data-testid="stTab"][role="tab"]:hover {
        background: transparent;
        color: var(--accent-strong);
    }
    div[data-testid="stTabs"] div[data-testid="stTab"][role="tab"][aria-selected="true"] {
        background: transparent;
        border-bottom-color: var(--accent);
        color: var(--accent-strong);
    }
    div[data-testid="stTabs"] div[data-testid="stTab"][role="tab"][aria-selected="true"] p {
        color: var(--accent-strong);
    }
    div[data-testid="stTabs"] .react-aria-SelectionIndicator,
    div[data-testid="stTabs"] div[data-baseweb="tab-highlight"],
    div[data-testid="stTabs"] div[data-baseweb="tab-border"] { display: none; }

    .tab-intro { margin: 0.15rem 0 1rem; max-width: 900px; }
    .tab-intro h2 {
        font-size: 1.18rem;
        line-height: 1.3;
        margin: 0 0 0.25rem;
    }
    .tab-intro p {
        color: var(--muted);
        font-size: 0.88rem;
        line-height: 1.55;
        margin: 0;
    }
    .chart-divider {
        border-top: 1px solid var(--line);
        margin: 1.65rem 0 1rem;
        padding-top: 1.15rem;
    }
    .chart-divider h3 { font-size: 1rem; margin: 0 0 0.2rem; }
    .chart-divider p { color: var(--muted); font-size: 0.84rem; margin: 0; }
    div[data-testid="stPlotlyChart"] {
        background: var(--surface);
        border: 1px solid var(--line);
        border-radius: 8px;
        box-shadow: 0 1px 2px rgba(24, 33, 31, 0.03);
        overflow: hidden;
        padding: 0.2rem;
    }
    [data-testid="stDataFrame"] {
        border: 1px solid var(--line);
        border-radius: 6px;
        overflow: hidden;
    }
    .export-label { color: var(--muted); margin: 0.65rem 0 0.35rem; }

    @media (max-width: 900px) {
        section[data-testid="stSidebar"], section[data-testid="stSidebar"] > div {
            width: min(460px, 94vw) !important;
        }
        section[data-testid="stSidebar"] [data-testid="stHorizontalBlock"]:has(.st-key-scenario-action-rail) {
            gap: 0.5rem !important;
        }
        .block-container { padding: 4.5rem 1rem 3rem; }
        .model-item { border-right: 0; }
    }
    @media (max-width: 640px) {
        .block-container { padding: 4.25rem 0.75rem 2.5rem; }
        h1 { font-size: 1.8rem; }
        .page-deck { font-size: 0.92rem; }
        [data-testid="stMetric"] { min-height: 98px; }
        div[data-testid="stTabs"] div[role="tablist"] { gap: 1.1rem; }
        .scenario-bar { align-items: flex-start; flex-direction: column; gap: 0.25rem; }
    }
    </style>
    """,
    unsafe_allow_html=True,
)


WIDGET_KEYS = [
    "archive_value",
    "archive_unit",
    "asset_value",
    "asset_unit",
    "retrieval",
    "start_year_widget",
    "horizon",
    "tech_dna",
    "tech_amazon",
    "tech_azure",
    "tech_tape",
    "tech_custom",
    "discount",
    "dna_cost_base_year",
    "dna_synthesis_cost",
    "dna_sequencing_cost",
    "synthesis_decline",
    "sequencing_decline",
    "amazon_base_year",
    "amazon_put_per_1000",
    "amazon_restore_per_1000",
    "amazon_retrieval_per_tb",
    "amazon_storage_per_tb_month",
    "amazon_decline",
    "azure_base_year",
    "azure_write_per_1000",
    "azure_read_per_1000",
    "azure_retrieval_per_tb",
    "azure_storage_per_tb_month",
    "azure_decline",
    "tape_base_year",
    "tape_media_per_tb",
    "tape_hardware_per_tb",
    "tape_energy_per_tb_year",
    "tape_media_decline",
    "tape_hardware_decline",
    "tape_energy_decline",
    "dna_durability",
    "tape_durability",
    "custom_name",
    "custom_base_year",
    "custom_write_tb",
    "custom_write_asset",
    "custom_storage_tb_year",
    "custom_retrieval_tb",
    "custom_retrieval_asset",
    "custom_decline",
    "custom_replacement",
    "projection_end",
    "log_scale",
]


def _default_projection_end(scenario: Scenario, query: dict[str, str]) -> int:
    return min(2500, max(scenario.start_year + 25, int(query.get("projection_end", 2350))))


def _query_mapping() -> dict[str, str]:
    return {key: str(value) for key, value in st.query_params.items()}


def _initial_scenario() -> Scenario:
    try:
        return Scenario.from_mapping(_query_mapping())
    except (TypeError, ValueError):
        st.warning("The shared URL contained invalid values. The paper baseline has been loaded instead.")
        return Scenario()


def _display_units(scenario: Scenario) -> tuple[float, str, float, str]:
    if scenario.archive_size_tb >= 1_000_000:
        archive_value, archive_unit = scenario.archive_size_tb / 1_000_000, "EB"
    elif scenario.archive_size_tb >= 1_000:
        archive_value, archive_unit = scenario.archive_size_tb / 1_000, "PB"
    else:
        archive_value, archive_unit = scenario.archive_size_tb, "TB"
    if scenario.average_asset_size_mb >= 1000:
        asset_value, asset_unit = scenario.average_asset_size_mb / 1000, "GB"
    else:
        asset_value, asset_unit = scenario.average_asset_size_mb, "MB"
    return archive_value, archive_unit, asset_value, asset_unit


def _widget_state_from_scenario(
    scenario: Scenario, query: dict[str, str] | None = None
) -> dict[str, bool | float | int | str]:
    archive_value, archive_unit, asset_value, asset_unit = _display_units(scenario)
    query = query or {}
    return {
        "archive_value": archive_value,
        "archive_unit": archive_unit,
        "asset_value": asset_value,
        "asset_unit": asset_unit,
        "retrieval": scenario.annual_retrieval_percent,
        "start_year_widget": scenario.start_year,
        "horizon": scenario.horizon_years,
        "tech_dna": "DNA" in scenario.technologies,
        "tech_amazon": "Amazon Deep Archive" in scenario.technologies,
        "tech_azure": "Azure Blob Archive" in scenario.technologies,
        "tech_tape": "Tape On-premise" in scenario.technologies,
        "tech_custom": "Custom storage" in scenario.technologies,
        "discount": scenario.discount_rate_percent,
        "dna_cost_base_year": scenario.dna_cost_base_year,
        "dna_synthesis_cost": scenario.dna_synthesis_cost_per_mb,
        "dna_sequencing_cost": scenario.dna_sequencing_cost_per_mb,
        "synthesis_decline": scenario.synthesis_decline_percent,
        "sequencing_decline": scenario.sequencing_decline_percent,
        "amazon_base_year": scenario.amazon_price_base_year,
        "amazon_put_per_1000": scenario.amazon_put_usd_per_request * 1_000,
        "amazon_restore_per_1000": scenario.amazon_bulk_restore_usd_per_request * 1_000,
        "amazon_retrieval_per_tb": scenario.amazon_bulk_retrieval_usd_per_mb * 1_000_000,
        "amazon_storage_per_tb_month": scenario.amazon_storage_usd_per_mb_month * 1_000_000,
        "amazon_decline": scenario.amazon_decline_percent,
        "azure_base_year": scenario.azure_price_base_year,
        "azure_write_per_1000": scenario.azure_write_usd_per_request * 1_000,
        "azure_read_per_1000": scenario.azure_read_usd_per_request * 1_000,
        "azure_retrieval_per_tb": scenario.azure_retrieval_usd_per_mb * 1_000_000,
        "azure_storage_per_tb_month": scenario.azure_storage_usd_per_mb_month * 1_000_000,
        "azure_decline": scenario.azure_decline_percent,
        "tape_base_year": scenario.tape_price_base_year,
        "tape_media_per_tb": scenario.tape_media_usd_per_tb,
        "tape_hardware_per_tb": scenario.tape_hardware_usd_per_tb,
        "tape_energy_per_tb_year": scenario.tape_energy_usd_per_tb_year,
        "tape_media_decline": scenario.tape_media_decline_percent,
        "tape_hardware_decline": scenario.tape_hardware_decline_percent,
        "tape_energy_decline": scenario.tape_energy_decline_percent,
        "dna_durability": scenario.dna_durability_years,
        "tape_durability": scenario.tape_durability_years,
        "custom_name": scenario.custom_storage_name,
        "custom_base_year": scenario.custom_cost_base_year,
        "custom_write_tb": scenario.custom_write_cost_per_tb,
        "custom_write_asset": scenario.custom_write_cost_per_asset,
        "custom_storage_tb_year": scenario.custom_storage_cost_per_tb_year,
        "custom_retrieval_tb": scenario.custom_retrieval_cost_per_tb,
        "custom_retrieval_asset": scenario.custom_retrieval_cost_per_asset,
        "custom_decline": scenario.custom_decline_percent,
        "custom_replacement": scenario.custom_replacement_years,
        "projection_end": _default_projection_end(scenario, query),
        "log_scale": query.get("log_scale", "True").lower() == "true",
    }


def _reset_to_paper_baseline() -> None:
    baseline_state = _widget_state_from_scenario(Scenario())
    st.session_state.update({key: baseline_state[key] for key in WIDGET_KEYS})
    st.session_state["committed_widgets"] = {key: st.session_state[key] for key in WIDGET_KEYS}
    st.session_state["form_generation"] = st.session_state.get("form_generation", 0) + 1
    st.query_params.clear()


def _snapshot_widgets() -> dict[str, bool | float | int | str]:
    return {key: st.session_state[key] for key in WIDGET_KEYS}


def _scenario_from_widgets(widgets: dict[str, bool | float | int | str]) -> Scenario:
    archive_multiplier = {"TB": 1, "PB": 1000, "EB": 1_000_000}[widgets["archive_unit"]]
    asset_multiplier = {"MB": 1, "GB": 1000}[widgets["asset_unit"]]
    technologies = tuple(
        technology
        for selected, technology in (
            (widgets["tech_dna"], "DNA"),
            (widgets["tech_amazon"], "Amazon Deep Archive"),
            (widgets["tech_azure"], "Azure Blob Archive"),
            (widgets["tech_tape"], "Tape On-premise"),
            (widgets["tech_custom"], "Custom storage"),
        )
        if selected
    )
    return Scenario(
        archive_size_tb=widgets["archive_value"] * archive_multiplier,
        average_asset_size_mb=widgets["asset_value"] * asset_multiplier,
        annual_retrieval_percent=widgets["retrieval"],
        start_year=int(widgets["start_year_widget"]),
        horizon_years=int(widgets["horizon"]),
        discount_rate_percent=widgets["discount"],
        dna_cost_base_year=int(widgets["dna_cost_base_year"]),
        dna_synthesis_cost_per_mb=widgets["dna_synthesis_cost"],
        dna_sequencing_cost_per_mb=widgets["dna_sequencing_cost"],
        synthesis_decline_percent=widgets["synthesis_decline"],
        sequencing_decline_percent=widgets["sequencing_decline"],
        amazon_price_base_year=int(widgets["amazon_base_year"]),
        amazon_put_usd_per_request=widgets["amazon_put_per_1000"] / 1_000,
        amazon_bulk_restore_usd_per_request=widgets["amazon_restore_per_1000"] / 1_000,
        amazon_bulk_retrieval_usd_per_mb=widgets["amazon_retrieval_per_tb"] / 1_000_000,
        amazon_storage_usd_per_mb_month=widgets["amazon_storage_per_tb_month"] / 1_000_000,
        amazon_decline_percent=widgets["amazon_decline"],
        azure_price_base_year=int(widgets["azure_base_year"]),
        azure_write_usd_per_request=widgets["azure_write_per_1000"] / 1_000,
        azure_read_usd_per_request=widgets["azure_read_per_1000"] / 1_000,
        azure_retrieval_usd_per_mb=widgets["azure_retrieval_per_tb"] / 1_000_000,
        azure_storage_usd_per_mb_month=widgets["azure_storage_per_tb_month"] / 1_000_000,
        azure_decline_percent=widgets["azure_decline"],
        tape_price_base_year=int(widgets["tape_base_year"]),
        tape_media_usd_per_tb=widgets["tape_media_per_tb"],
        tape_hardware_usd_per_tb=widgets["tape_hardware_per_tb"],
        tape_energy_usd_per_tb_year=widgets["tape_energy_per_tb_year"],
        tape_media_decline_percent=widgets["tape_media_decline"],
        tape_hardware_decline_percent=widgets["tape_hardware_decline"],
        tape_energy_decline_percent=widgets["tape_energy_decline"],
        dna_durability_years=int(widgets["dna_durability"]),
        tape_durability_years=int(widgets["tape_durability"]),
        custom_storage_name=widgets["custom_name"],
        custom_cost_base_year=int(widgets["custom_base_year"]),
        custom_write_cost_per_tb=widgets["custom_write_tb"],
        custom_write_cost_per_asset=widgets["custom_write_asset"],
        custom_storage_cost_per_tb_year=widgets["custom_storage_tb_year"],
        custom_retrieval_cost_per_tb=widgets["custom_retrieval_tb"],
        custom_retrieval_cost_per_asset=widgets["custom_retrieval_asset"],
        custom_decline_percent=widgets["custom_decline"],
        custom_replacement_years=int(widgets["custom_replacement"]),
        technologies=technologies,
    )


def _money(value: float) -> str:
    absolute = abs(value)
    for threshold, suffix in ((1e12, "T"), (1e9, "B"), (1e6, "M"), (1e3, "K")):
        if absolute >= threshold:
            return f"${value / threshold:,.2f}{suffix}"
    if absolute >= 1:
        return f"${value:,.2f}"
    return f"${value:,.3g}"


def _quantity(value: float) -> str:
    absolute = abs(value)
    for threshold, suffix in ((1e15, "Q"), (1e12, "T"), (1e9, "B"), (1e6, "M"), (1e3, "K")):
        if absolute >= threshold:
            return f"{value / threshold:,.2f}{suffix}"
    return f"{value:,.0f}"


def _section_intro(title: str, description: str) -> None:
    st.markdown(
        f'<div class="tab-intro"><h2>{title}</h2><p>{description}</p></div>',
        unsafe_allow_html=True,
    )


def _chart_downloads(
    key: str,
    csv_data: bytes,
    filename_base: str,
    image_exports: dict[str, bytes],
) -> None:
    st.markdown('<div class="export-label">Export chart</div>', unsafe_allow_html=True)
    columns = st.columns(3)
    columns[0].download_button(
        "CSV",
        csv_data,
        f"{filename_base}.csv",
        "text/csv",
        key=f"download_{key}_csv",
        icon=":material/download:",
        width="stretch",
        help="Download the data shown in this graph.",
    )
    columns[1].download_button(
        "PNG",
        image_exports["png"],
        f"{filename_base}.png",
        "image/png",
        key=f"download_{key}_png",
        icon=":material/download:",
        width="stretch",
        help="Download this graph as a high-resolution PNG image.",
    )
    columns[2].download_button(
        "SVG",
        image_exports["svg"],
        f"{filename_base}.svg",
        "image/svg+xml",
        key=f"download_{key}_svg",
        icon=":material/download:",
        width="stretch",
        help="Download this graph as an editable vector image.",
    )


def _plot_config(filename: str) -> dict:
    return {
        "displaylogo": False,
        "toImageButtonOptions": {"format": "png", "filename": filename, "scale": 2},
    }


@st.cache_data(show_spinner=False)
def _cached_simulation(scenario: Scenario):
    return simulate_scenario(scenario)


@st.cache_data(show_spinner=False)
def _cached_projection(scenario: Scenario, final_year: int) -> pd.DataFrame:
    return simulate_start_years(scenario, final_year)


@st.cache_data(show_spinner=False)
def _cached_dna_costs(scenario: Scenario, final_year: int) -> pd.DataFrame:
    return simulate_dna_unit_costs(scenario, final_year)


initial = _initial_scenario()
query = _query_mapping()
initial_widgets = _widget_state_from_scenario(initial, query)
for widget_key, default_value in initial_widgets.items():
    st.session_state.setdefault(widget_key, default_value)
# The graphs only follow the last calculated inputs: the initial load counts
# as the first calculation, and the Calculate button is the only other way
# to commit widget changes.
st.session_state.setdefault("committed_widgets", dict(initial_widgets))
committed_widgets = st.session_state["committed_widgets"]
pending = _snapshot_widgets() != committed_widgets

with st.sidebar:
    input_column, action_column = st.columns([6, 1], gap="small")
    with input_column:
        st.markdown(
            """
            <div class="sidebar-kicker">Scenario builder</div>
            <div class="sidebar-title">Model inputs</div>
            <p class="sidebar-copy">Archive workload, time horizon, and technology assumptions.</p>
            """,
            unsafe_allow_html=True,
        )
        st.button(
            "Reset to paper baseline",
            key="reset_baseline",
            on_click=_reset_to_paper_baseline,
            icon=":material/refresh:",
            width="stretch",
        )
        st.markdown('<div class="sidebar-section">Workload and time</div>', unsafe_allow_html=True)
        col_a, col_b = st.columns([2, 1])
        with col_a:
            archive_input = st.number_input(
                "Archive size",
                min_value=0.001,
                key="archive_value",
                help="Total logical data stored in the collection, before DNA coding, redundancy, or provider replication.",
            )
        with col_b:
            archive_unit_input = st.selectbox(
                "Unit", ["TB", "PB", "EB"], key="archive_unit",
                help="Decimal capacity unit: 1 PB = 1,000 TB and 1 EB = 1,000,000 TB.",
            )

        col_a, col_b = st.columns([2, 1])
        with col_a:
            asset_input = st.number_input(
                "Average asset size",
                min_value=0.001,
                key="asset_value",
                help="Average data object size. It determines object count.",
            )
        with col_b:
            asset_unit_input = st.selectbox(
                "Unit ", ["MB", "GB"], key="asset_unit",
                help="Unit used for the average size of one asset.",
            )

        time_col_a, time_col_b = st.columns(2)
        with time_col_a:
            start_year = st.number_input(
                "Start year", min_value=2025, max_value=2500,
                key="start_year_widget", help="Calendar year in which the archive is first written.",
            )
        with time_col_b:
            horizon = st.number_input(
                "Retention (years)", min_value=1, max_value=10_000,
                key="horizon",
                help="Number of charged storage years, including the start year.",
            )

        finance_col_a, finance_col_b = st.columns(2)
        with finance_col_a:
            retrieval = st.number_input(
                "Annual retrieval (%)", min_value=0.0, max_value=10_000.0,
                step=0.25, key="retrieval",
                help="Expected share of the logical archive retrieved each year. 1% means reading 10 TB per year from a 1 PB archive.",
            )
        with finance_col_b:
            discount = st.number_input(
                "Discount rate (%)", min_value=0.0, max_value=99.0,
                step=0.25, key="discount",
                help="Real rate used to discount future payments to the storage start year. Use 0 for undiscounted costs.",
            )

        st.markdown('<div class="sidebar-section">Display</div>', unsafe_allow_html=True)
        chart_col_a, chart_col_b = st.columns([1.4, 1])
        with chart_col_a:
            projection_end = st.number_input(
                "Outlook end year", min_value=int(start_year), max_value=2500,
                key="projection_end",
                help="Final archive start year included in the start-year outlook chart.",
            )
        with chart_col_b:
            log_scale = st.toggle(
                "Log", key="log_scale",
                help="Recommended when technologies differ by several orders of magnitude.",
            )

        st.markdown('<div class="sidebar-section">Technologies</div>', unsafe_allow_html=True)
        with st.container(border=True):
            tech_col_a, tech_col_b = st.columns(2)
            with tech_col_a:
                tech_dna = st.checkbox(
                    "DNA", key="tech_dna", help="Archival storage using DNA synthesis and sequencing.",
                )
                tech_amazon = st.checkbox(
                    "Amazon S3", key="tech_amazon", help="Amazon S3 Glacier Deep Archive.",
                )
                tech_tape = st.checkbox(
                    "Tape", key="tech_tape", help="On-premise tape storage with periodic media replacement.",
                )
            with tech_col_b:
                tech_azure = st.checkbox(
                    "Azure Blob", key="tech_azure", help="Microsoft Azure Blob Storage Archive tier.",
                )
                tech_custom = st.checkbox(
                    "Custom", key="tech_custom", help="User-defined storage cost model.",
                )

        st.markdown('<div class="sidebar-section">Cost assumptions</div>', unsafe_allow_html=True)
        with st.expander("DNA cost assumptions"):
            dna_cost_base_year = st.number_input(
                "DNA cost base year", min_value=2000, max_value=2500,
                key="dna_cost_base_year",
                help="Year to which the editable synthesis and sequencing unit costs apply.",
            )
            dna_synthesis_cost = st.number_input(
                "Synthesis cost (USD/MB)", min_value=0.0,
                format="%.6f", key="dna_synthesis_cost",
                help="Cost in the DNA cost base year to synthesize enough bases for 1 MB of logical data, before redundancy and indexing overhead.",
            )
            dna_sequencing_cost = st.number_input(
                "Sequencing cost (USD/MB)", min_value=0.0,
                format="%.8f", key="dna_sequencing_cost",
                help="Cost in the DNA cost base year to sequence 1 MB of retrieved logical data.",
            )
            synthesis_decline = st.number_input(
                "Synthesis annual decline (%)", min_value=0.0, max_value=99.99,
                key="synthesis_decline",
                help="Percentage by which synthesis cost is assumed to fall each calendar year.",
            )
            sequencing_decline = st.number_input(
                "Sequencing annual decline (%)", min_value=0.0, max_value=99.99,
                key="sequencing_decline",
                help="Percentage by which sequencing cost is assumed to fall each calendar year.",
            )
            dna_durability = st.number_input(
                "DNA durability (years)", min_value=1, max_value=10_000,
                key="dna_durability",
                help="Years before the archive must be synthesized again. No replacement occurs at the exact end of the horizon.",
            )

        with st.expander("Amazon Deep Archive assumptions"):
            st.caption("Price reference")
            amazon_base_year = st.number_input(
                "Amazon price base year", min_value=2000, max_value=2500,
                key="amazon_base_year",
                help="Calendar year to which all Amazon prices below apply.",
            )
            amazon_decline = st.number_input(
                "Amazon annual price decline (%)", min_value=0.0, max_value=99.99,
                key="amazon_decline",
                help="Annual reduction applied to Amazon request, retrieval, and storage prices.",
            )
            st.caption("Base-year prices")
            amazon_put_per_1000 = st.number_input(
                "Write requests (USD/1,000)", min_value=0.0,
                format="%.6f",
                key="amazon_put_per_1000",
                help="Charge for 1,000 requests when the archive is initially written.",
            )
            amazon_restore_per_1000 = st.number_input(
                "Bulk restore requests (USD/1,000)", min_value=0.0,
                format="%.6f",
                key="amazon_restore_per_1000",
                help="Charge for 1,000 bulk restore-job requests. Asset size determines the request count.",
            )
            amazon_retrieval_per_tb = st.number_input(
                "Bulk data retrieval (USD/TB)", min_value=0.0,
                format="%.6f",
                key="amazon_retrieval_per_tb",
                help="Capacity charge for retrieving one TB of archived data.",
            )
            amazon_storage_per_tb_month = st.number_input(
                "Storage (USD/TB/month)", min_value=0.0,
                format="%.6f",
                key="amazon_storage_per_tb_month",
                help="Recurring monthly charge to retain one TB in Deep Archive.",
            )

        with st.expander("Azure Blob Archive assumptions"):
            st.caption("Price reference")
            azure_base_year = st.number_input(
                "Azure price base year", min_value=2000, max_value=2500,
                key="azure_base_year",
                help="Calendar year to which all Azure prices below apply.",
            )
            azure_decline = st.number_input(
                "Azure annual price decline (%)", min_value=0.0, max_value=99.99,
                key="azure_decline",
                help="Annual reduction applied to Azure request, retrieval, and storage prices.",
            )
            st.caption("Base-year prices")
            azure_write_per_1000 = st.number_input(
                "Write requests (USD/1,000)", min_value=0.0,
                format="%.6f",
                key="azure_write_per_1000",
                help="Charge for 1,000 requests when the archive is initially written.",
            )
            azure_read_per_1000 = st.number_input(
                "Read requests (USD/1,000)", min_value=0.0,
                format="%.6f",
                key="azure_read_per_1000",
                help="Charge for 1,000 retrieval requests. Asset size determines the request count.",
            )
            azure_retrieval_per_tb = st.number_input(
                "Data retrieval (USD/TB)", min_value=0.0,
                format="%.6f",
                key="azure_retrieval_per_tb",
                help="Capacity charge for retrieving one TB from the Archive tier.",
            )
            azure_storage_per_tb_month = st.number_input(
                "Storage (USD/TB/month)", min_value=0.0,
                format="%.6f",
                key="azure_storage_per_tb_month",
                help="Recurring monthly charge to retain one TB in the Archive tier.",
            )

        with st.expander("Tape on-premise assumptions"):
            st.caption("Price reference")
            tape_base_year = st.number_input(
                "Tape price base year", min_value=2000, max_value=2500,
                key="tape_base_year",
                help="Calendar year to which the tape media, hardware, and energy prices apply.",
            )
            tape_durability = st.number_input(
                "Tape durability (years)", min_value=1, max_value=1_000,
                key="tape_durability",
                help="Years between complete tape media replacement writes.",
            )
            st.caption("Base-year prices")
            tape_media_per_tb = st.number_input(
                "Tape media (USD/TB)", min_value=0.0,
                format="%.6f",
                key="tape_media_per_tb",
                help="Media purchase cost for one TB, charged on the initial write and every replacement.",
            )
            tape_hardware_per_tb = st.number_input(
                "Hardware (USD/TB)", min_value=0.0,
                format="%.6f",
                key="tape_hardware_per_tb",
                help="Hardware cost allocated per TB and amortized over the selected durability period.",
            )
            tape_energy_per_tb_year = st.number_input(
                "Energy (USD/TB/year)", min_value=0.0,
                format="%.6f",
                key="tape_energy_per_tb_year",
                help="Annual energy cost to retain one TB in the tape system.",
            )
            st.caption("Annual price declines")
            tape_media_decline = st.number_input(
                "Tape media decline (%)", min_value=0.0, max_value=99.99,
                key="tape_media_decline",
                help="Annual reduction applied to tape media purchase prices.",
            )
            tape_hardware_decline = st.number_input(
                "Tape hardware decline (%)", min_value=0.0, max_value=99.99,
                key="tape_hardware_decline",
                help="Annual reduction applied to amortized tape hardware costs.",
            )
            tape_energy_decline = st.number_input(
                "Tape energy decline (%)", min_value=0.0, max_value=99.99,
                key="tape_energy_decline",
                help="Annual reduction applied to tape energy costs.",
            )

        with st.expander("Custom storage assumptions"):
            custom_name = st.text_input(
                "Display name", key="custom_name",
                help="Name used for the custom technology in charts, tables, and downloads.",
            )
            custom_base_year = st.number_input(
                "Price base year", min_value=2000, max_value=2500,
                key="custom_base_year", help="Year to which all custom prices apply.",
            )
            custom_write_tb = st.number_input(
                "Initial write cost (USD/TB)", min_value=0.0,
                key="custom_write_tb", help="Capacity-based charge to write or replace one TB.",
            )
            custom_write_asset = st.number_input(
                "Write request cost (USD/asset)", min_value=0.0,
                key="custom_write_asset",
                help="Per-file or per-object charge applied when the archive is written or replaced.",
            )
            custom_storage_tb_year = st.number_input(
                "Annual storage cost (USD/TB)", min_value=0.0,
                key="custom_storage_tb_year",
                help="Recurring cost to retain one TB for one year.",
            )
            custom_retrieval_tb = st.number_input(
                "Retrieval cost (USD/TB)", min_value=0.0,
                key="custom_retrieval_tb",
                help="Capacity-based charge for each TB retrieved.",
            )
            custom_retrieval_asset = st.number_input(
                "Retrieval request cost (USD/asset)", min_value=0.0,
                key="custom_retrieval_asset",
                help="Per-file or per-object charge for the expected assets retrieved each year.",
            )
            custom_decline = st.number_input(
                "Annual price decline (%)", min_value=0.0, max_value=99.99,
                key="custom_decline",
                help="Annual percentage reduction applied to every custom price.",
            )
            custom_replacement = st.number_input(
                "Replacement interval (years)", min_value=0, max_value=10_000,
                key="custom_replacement",
                help="Years between complete rewrites. Use 0 for a service with no replacement writes.",
            )

    with action_column:
        with st.container(key="scenario-action-rail"):
            submitted = st.button(
                "Calculate scenario",
                type="primary",
                key="calculate_scenario",
                help=(
                    "Graphs are out of date — apply the current model inputs."
                    if pending
                    else "Apply the current model inputs and update all graphs."
                ),
                width="stretch",
            )

if submitted:
    technologies = tuple(
        technology
        for selected, technology in (
            (tech_dna, "DNA"),
            (tech_amazon, "Amazon Deep Archive"),
            (tech_azure, "Azure Blob Archive"),
            (tech_tape, "Tape On-premise"),
            (tech_custom, "Custom storage"),
        )
        if selected
    )
    archive_multiplier = {"TB": 1, "PB": 1000, "EB": 1_000_000}[archive_unit_input]
    asset_multiplier = {"MB": 1, "GB": 1000}[asset_unit_input]
    try:
        candidate = Scenario(
            archive_size_tb=archive_input * archive_multiplier,
            average_asset_size_mb=asset_input * asset_multiplier,
            annual_retrieval_percent=retrieval,
            start_year=int(start_year),
            horizon_years=int(horizon),
            discount_rate_percent=discount,
            dna_cost_base_year=int(dna_cost_base_year),
            dna_synthesis_cost_per_mb=dna_synthesis_cost,
            dna_sequencing_cost_per_mb=dna_sequencing_cost,
            synthesis_decline_percent=synthesis_decline,
            sequencing_decline_percent=sequencing_decline,
            amazon_price_base_year=int(amazon_base_year),
            amazon_put_usd_per_request=amazon_put_per_1000 / 1_000,
            amazon_bulk_restore_usd_per_request=amazon_restore_per_1000 / 1_000,
            amazon_bulk_retrieval_usd_per_mb=amazon_retrieval_per_tb / 1_000_000,
            amazon_storage_usd_per_mb_month=amazon_storage_per_tb_month / 1_000_000,
            amazon_decline_percent=amazon_decline,
            azure_price_base_year=int(azure_base_year),
            azure_write_usd_per_request=azure_write_per_1000 / 1_000,
            azure_read_usd_per_request=azure_read_per_1000 / 1_000,
            azure_retrieval_usd_per_mb=azure_retrieval_per_tb / 1_000_000,
            azure_storage_usd_per_mb_month=azure_storage_per_tb_month / 1_000_000,
            azure_decline_percent=azure_decline,
            tape_price_base_year=int(tape_base_year),
            tape_media_usd_per_tb=tape_media_per_tb,
            tape_hardware_usd_per_tb=tape_hardware_per_tb,
            tape_energy_usd_per_tb_year=tape_energy_per_tb_year,
            tape_media_decline_percent=tape_media_decline,
            tape_hardware_decline_percent=tape_hardware_decline,
            tape_energy_decline_percent=tape_energy_decline,
            dna_durability_years=int(dna_durability),
            tape_durability_years=int(tape_durability),
            custom_storage_name=custom_name,
            custom_cost_base_year=int(custom_base_year),
            custom_write_cost_per_tb=custom_write_tb,
            custom_write_cost_per_asset=custom_write_asset,
            custom_storage_cost_per_tb_year=custom_storage_tb_year,
            custom_retrieval_cost_per_tb=custom_retrieval_tb,
            custom_retrieval_cost_per_asset=custom_retrieval_asset,
            custom_decline_percent=custom_decline,
            custom_replacement_years=int(custom_replacement),
            technologies=technologies,
        )
    except ValueError as error:
        st.error(str(error))
        st.stop()
    if not candidate.technologies:
        st.error("Select at least one storage technology.")
        st.stop()
    st.session_state["committed_widgets"] = _snapshot_widgets()
    params = candidate.to_query_params()
    params.update({"projection_end": str(int(projection_end)), "log_scale": str(log_scale)})
    st.query_params.from_dict(params)

committed_widgets = st.session_state["committed_widgets"]
pending = _snapshot_widgets() != committed_widgets

# All graphs render from the last calculated inputs; live widget edits do not
# touch them until the Calculate button commits a new snapshot.
try:
    scenario = _scenario_from_widgets(committed_widgets)
except ValueError as error:
    st.error(str(error))
    st.stop()
if not scenario.technologies:
    st.error("Select at least one storage technology.")
    st.stop()
projection_end = int(committed_widgets["projection_end"])
log_scale = bool(committed_widgets["log_scale"])

result = _cached_simulation(scenario)
projection = _cached_projection(scenario, int(projection_end))
dna_curve_end = max(
    scenario.dna_cost_base_year,
    min(2500, scenario.start_year + scenario.horizon_years - 1),
)
dna_costs = _cached_dna_costs(scenario, dna_curve_end)
use_present_value = scenario.discount_rate_percent > 0
value_column = "present_value_usd" if use_present_value else "total_cost_usd"

st.markdown('<div class="page-kicker">Archival storage economics</div>', unsafe_allow_html=True)
st.title("DNA Storage Cost Explorer")
st.markdown(
    '<p class="page-deck">Compare the long-run cost of DNA, cloud archive, tape, and custom storage '
    'under one consistent workload.</p>',
    unsafe_allow_html=True,
)
st.markdown(
    f"""
    <div class="model-strip">
        <span class="model-item"><strong>Model v{result.metadata['model_version']}</strong></span>
        <span class="model-item">{result.metadata['currency']}</span>
        <span class="model-item">Reviewed {result.metadata['last_reviewed']}</span>
        <span class="model-item">{result.metadata['disclaimer']}</span>
    </div>
    """,
    unsafe_allow_html=True,
)

totals = result.totals.sort_values(value_column)
cheapest = totals.iloc[0]
cheapest_name = str(cheapest["technology"])
cheapest_label = cheapest_name if len(cheapest_name) <= 26 else f"{cheapest_name[:25]}..."
dna_rows = totals[totals["technology"] == "DNA"]
dna_total = float(dna_rows.iloc[0][value_column]) if not dna_rows.empty else None
period_end = scenario.start_year + scenario.horizon_years - 1
cost_basis = "Present value" if use_present_value else "Undiscounted"
scale_label = "Log scale" if log_scale else "Linear scale"

pending_chip = (
    '<span class="scenario-pending">Changes pending — press Calculate</span>' if pending else ""
)
st.markdown(
    f"""
    <div class="scenario-bar">
        <span class="scenario-label">Active scenario</span>
        <span><strong>{scenario.start_year}-{period_end}</strong> archive period</span>
        <span><strong>{len(scenario.technologies)}</strong> technologies</span>
        <span>{cost_basis}</span>
        <span>{scale_label}</span>
        {pending_chip}
    </div>
    """,
    unsafe_allow_html=True,
)

metric_columns = st.columns(4)
metric_columns[0].metric("Archive", f"{scenario.archive_size_tb:,.3g} TB")
metric_columns[1].metric(
    "Data objects",
    _quantity(scenario.number_of_assets),
    help=f"{scenario.number_of_assets:,.0f} total data objects",
)
metric_columns[2].metric(
    f"Lowest: {cheapest_label}",
    _money(float(cheapest[value_column])),
    help=f"Lowest lifecycle cost: {cheapest_name}",
)
if dna_total is None:
    metric_columns[3].metric("Annual retrieval", f"{scenario.annual_retrieval_percent:,.2f}%")
else:
    lowest_cost = float(cheapest[value_column])
    if lowest_cost > 0:
        ratio = dna_total / lowest_cost
        metric_columns[3].metric(
            "DNA lifecycle cost", _money(dna_total), f"{ratio:,.2g}x lowest", delta_color="off"
        )
    else:
        metric_columns[3].metric("DNA lifecycle cost", _money(dna_total))

st.markdown('<div class="workspace-kicker">Analysis workspace</div>', unsafe_allow_html=True)
overview_tab, outlook_tab, dna_cost_tab, assumptions_tab = st.tabs(
    ["Lifecycle", "Start-year outlook", "DNA unit costs", "Assumptions"]
)

with overview_tab:
    _section_intro(
        "Lifecycle comparison",
        "Cumulative lifecycle cost for one archive opened in the selected start year. "
        "Each line includes initial and replacement writes, storage or operation, and expected retrieval.",
    )
    lifecycle_figure = lifecycle_chart(result, use_present_value, log_scale)
    st.plotly_chart(
        lifecycle_figure,
        width="stretch",
        config=_plot_config("dna-storage-lifecycle"),
    )
    lifecycle_files = lifecycle_exports(result, use_present_value, log_scale)
    _chart_downloads(
        "lifecycle",
        result.yearly.to_csv(index=False).encode("utf-8"),
        "dna-storage-lifecycle",
        lifecycle_files,
    )
    st.markdown(
        """
        <div class="chart-divider">
            <h3>Cost composition</h3>
            <p>Undiscounted write and replacement, retrieval, and recurring storage or operating costs.</p>
        </div>
        """,
        unsafe_allow_html=True,
    )
    breakdown_figure = breakdown_chart(result, log_scale)
    st.plotly_chart(
        breakdown_figure,
        width="stretch",
        config=_plot_config("dna-storage-cost-components"),
    )
    breakdown_columns = ["technology", *[component for component in (
        "write_cost_usd", "read_cost_usd", "maintenance_cost_usd", "total_cost_usd"
    )]]
    _chart_downloads(
        "breakdown",
        result.totals[breakdown_columns].to_csv(index=False).encode("utf-8"),
        "dna-storage-cost-components",
        breakdown_exports(result, log_scale),
    )

with outlook_tab:
    _section_intro(
        "Start-year sensitivity",
        "The same archive workload and retention horizon are recalculated for every possible storage "
        "start year. A crossover is the first start year for which DNA's lifecycle cost is no greater "
        "than the comparison technology.",
    )
    projection_figure = projection_chart(projection, use_present_value, log_scale)
    st.plotly_chart(
        projection_figure,
        width="stretch",
        config=_plot_config("dna-storage-start-year-outlook"),
    )
    _chart_downloads(
        "projection",
        projection.to_csv(index=False).encode("utf-8"),
        "dna-storage-start-year-outlook",
        projection_exports(projection, use_present_value, log_scale, result.metadata),
    )
    st.markdown(
        """
        <div class="chart-divider">
            <h3>Crossover years</h3>
            <p>First modeled start year in which DNA reaches or undercuts each comparison technology.</p>
        </div>
        """,
        unsafe_allow_html=True,
    )
    crossovers = find_crossover_years(projection)
    if crossovers:
        crossover_rows = [
            {
                "Comparison": f"DNA <= {technology}",
                "First start year": str(year) if year is not None else f"Not by {projection_end}",
            }
            for technology, year in crossovers.items()
        ]
        st.dataframe(pd.DataFrame(crossover_rows), hide_index=True, width="stretch")
    else:
        st.info("Include DNA and at least one comparison technology to calculate crossover years.")

with dna_cost_tab:
    _section_intro(
        "DNA unit economics",
        f"Both curves begin with the editable {scenario.dna_cost_base_year} unit costs and apply the "
        f"selected annual decline rates through {dna_curve_end}. They are unit-cost assumptions, "
        "not lifecycle totals.",
    )
    chart_columns = st.columns(2)
    with chart_columns[0]:
        synthesis_title = "DNA synthesis cost trajectory"
        st.plotly_chart(
            dna_unit_cost_chart(
                dna_costs,
                "synthesis_cost_usd_per_mb",
                synthesis_title,
                "#bd4b38",
                log_scale,
            ),
            width="stretch",
            config=_plot_config("dna-synthesis-cost-trajectory"),
        )
        st.caption(
            "Modeled cost to synthesize enough DNA bases for 1 MB of logical data. "
            "Redundancy, indexing, and coding overhead are not added here."
        )
        synthesis_data = dna_costs[["year", "synthesis_cost_usd_per_mb"]]
        _chart_downloads(
            "dna_synthesis",
            synthesis_data.to_csv(index=False).encode("utf-8"),
            "dna-synthesis-cost-trajectory",
            dna_unit_cost_exports(
                dna_costs,
                "synthesis_cost_usd_per_mb",
                synthesis_title,
                "#bd4b38",
                log_scale,
                result.metadata,
            ),
        )
    with chart_columns[1]:
        sequencing_title = "DNA sequencing cost trajectory"
        st.plotly_chart(
            dna_unit_cost_chart(
                dna_costs,
                "sequencing_cost_usd_per_mb",
                sequencing_title,
                "#176b87",
                log_scale,
            ),
            width="stretch",
            config=_plot_config("dna-sequencing-cost-trajectory"),
        )
        st.caption(
            "Modeled cost to sequence and retrieve 1 MB of logical data. "
            "It is applied to the share of the archive retrieved each year."
        )
        sequencing_data = dna_costs[["year", "sequencing_cost_usd_per_mb"]]
        _chart_downloads(
            "dna_sequencing",
            sequencing_data.to_csv(index=False).encode("utf-8"),
            "dna-sequencing-cost-trajectory",
            dna_unit_cost_exports(
                dna_costs,
                "sequencing_cost_usd_per_mb",
                sequencing_title,
                "#176b87",
                log_scale,
                result.metadata,
            ),
        )

with assumptions_tab:
    assumptions = load_assumptions()
    _section_intro(
        "Assumptions and scope",
        "A concise record of the active scenario, included cost categories, and source references.",
    )
    contract = pd.DataFrame(
        [
            ("Archive period", f"{scenario.start_year}-{scenario.start_year + scenario.horizon_years - 1}"),
            ("Annual retrieval", f"{scenario.annual_assets_retrieved:,.2f} assets ({scenario.annual_retrieval_percent:,.2f}%)"),
            ("Discounting", f"{scenario.discount_rate_percent:,.2f}% real; present value at storage start"),
            ("Currency", result.metadata["currency"]),
            (
                "DNA synthesis baseline",
                f"${scenario.dna_synthesis_cost_per_mb:,.6g}/MB in {scenario.dna_cost_base_year}",
            ),
            (
                "DNA sequencing baseline",
                f"${scenario.dna_sequencing_cost_per_mb:,.6g}/MB in {scenario.dna_cost_base_year}",
            ),
            (
                "Amazon baseline",
                f"${scenario.amazon_storage_usd_per_mb_month * 1_000_000:,.6g}/TB/month storage; "
                f"${scenario.amazon_bulk_retrieval_usd_per_mb * 1_000_000:,.6g}/TB retrieval "
                f"in {scenario.amazon_price_base_year}",
            ),
            (
                "Azure baseline",
                f"${scenario.azure_storage_usd_per_mb_month * 1_000_000:,.6g}/TB/month storage; "
                f"${scenario.azure_retrieval_usd_per_mb * 1_000_000:,.6g}/TB retrieval "
                f"in {scenario.azure_price_base_year}",
            ),
            (
                "Tape baseline",
                f"${scenario.tape_media_usd_per_tb:,.6g}/TB media; "
                f"${scenario.tape_hardware_usd_per_tb:,.6g}/TB hardware; "
                f"${scenario.tape_energy_usd_per_tb_year:,.6g}/TB/year energy "
                f"in {scenario.tape_price_base_year}",
            ),
            ("Horizon convention", result.metadata["horizon_convention"]),
        ],
        columns=["Item", "Value"],
    )
    contract_column, context_column = st.columns([1.45, 1], gap="large")
    with contract_column:
        st.subheader("Scenario contract")
        st.dataframe(contract, hide_index=True, width="stretch")
    with context_column:
        st.subheader("Cost scope")
        st.write(
            "Included: DNA synthesis and sequencing, cloud write/retrieval/storage charges, "
            "tape media/hardware/energy assumptions, and the selected custom capacity/request charges. "
            "Excluded: labor, cloud egress, taxes, retrieval latency, minimum-storage penalties, "
            "facilities, and unmodeled migrations."
        )
        st.subheader("Sources")
        for source in assumptions["sources"].values():
            st.markdown(f"- [{source['label']}]({source['url']})")
