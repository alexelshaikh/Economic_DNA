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
    :root { --ink: #202326; --muted: #62686d; --line: #dfe2e4; --accent: #b44737; }
    .stApp { background: #f7f8f8; color: var(--ink); }
    [data-testid="stSidebar"] { background: #ffffff; border-right: 1px solid var(--line); }
    [data-testid="stMetric"] {
        background: #ffffff; border: 1px solid var(--line); border-radius: 6px;
        padding: 14px 16px; min-height: 112px;
    }
    [data-testid="stMetricLabel"] { color: var(--muted); }
    [data-testid="stMetricValue"] { font-size: 1.55rem; }
    div[data-testid="stPlotlyChart"] { background: #ffffff; border: 1px solid var(--line); border-radius: 6px; }
    .model-strip {
        background: #ffffff; border-left: 4px solid var(--accent); border-top: 1px solid var(--line);
        border-right: 1px solid var(--line); border-bottom: 1px solid var(--line);
        padding: 10px 14px; margin: 4px 0 18px 0; color: var(--muted); font-size: 0.9rem;
    }
    .block-container { max-width: 1480px; padding-top: 1.8rem; }
    h1, h2, h3 { letter-spacing: 0; }
    h1 { font-size: 2rem; margin-bottom: 0.15rem; }
    .stButton > button, .stDownloadButton > button { border-radius: 5px; }
    div[data-testid="stTabs"] div[role="tablist"] {
        gap: 8px; padding: 6px; margin: 4px 0 14px 0; overflow-x: auto;
        background: #eef0f1; border: 1px solid #d7dade; border-radius: 6px;
    }
    div[data-testid="stTabs"] div[data-testid="stTab"][role="tab"] {
        flex: 0 0 auto; min-height: 42px; padding: 8px 16px;
        background: #ffffff; border: 1px solid #c8cdd1; border-radius: 5px;
        color: #34393d; font-weight: 600; cursor: pointer;
    }
    div[data-testid="stTabs"] div[data-testid="stTab"][role="tab"]:hover {
        background: #f8ece9; border-color: var(--accent); color: #923829;
    }
    div[data-testid="stTabs"] div[data-testid="stTab"][role="tab"][aria-selected="true"] {
        background: var(--accent); border-color: var(--accent); color: #ffffff;
    }
    div[data-testid="stTabs"] div[data-testid="stTab"][role="tab"][aria-selected="true"] p {
        color: #ffffff;
    }
    div[data-testid="stTabs"] .react-aria-SelectionIndicator { display: none; }
    div[data-testid="stTabs"] div[data-baseweb="tab-highlight"],
    div[data-testid="stTabs"] div[data-baseweb="tab-border"] { display: none; }
    @media (max-width: 700px) {
        .block-container { padding: 1rem 0.75rem; }
        h1 { font-size: 1.65rem; }
        [data-testid="stMetric"] { min-height: 94px; }
        div[data-testid="stTabs"] div[data-testid="stTab"][role="tab"] {
            min-height: 40px; padding: 7px 12px;
        }
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


def _money(value: float) -> str:
    absolute = abs(value)
    for threshold, suffix in ((1e12, "T"), (1e9, "B"), (1e6, "M"), (1e3, "K")):
        if absolute >= threshold:
            return f"${value / threshold:,.2f}{suffix}"
    if absolute >= 1:
        return f"${value:,.2f}"
    return f"${value:,.3g}"


def _chart_downloads(
    key: str,
    csv_data: bytes,
    filename_base: str,
    image_exports: dict[str, bytes],
) -> None:
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
archive_value, archive_unit, asset_value, asset_unit = _display_units(initial)
query = _query_mapping()
default_projection_end = min(2500, max(initial.start_year + 25, int(query.get("projection_end", 2350))))
default_log_scale = query.get("log_scale", "True").lower() == "true"

with st.sidebar:
    st.header("Scenario")
    if st.button("Reset to paper baseline", width="stretch"):
        for key in WIDGET_KEYS:
            st.session_state.pop(key, None)
        st.query_params.clear()
        st.rerun()

    with st.form("scenario_form"):
        st.caption("Workload and time")
        col_a, col_b = st.columns([2, 1])
        with col_a:
            archive_input = st.number_input(
                "Archive size",
                min_value=0.001,
                value=float(archive_value),
                key="archive_value",
                help="Total logical data stored in the collection, before DNA coding, redundancy, or provider replication.",
            )
        with col_b:
            archive_unit_input = st.selectbox(
                "Unit", ["TB", "PB", "EB"], index=["TB", "PB", "EB"].index(archive_unit), key="archive_unit",
                help="Decimal capacity unit: 1 PB = 1,000 TB and 1 EB = 1,000,000 TB.",
            )

        col_a, col_b = st.columns([2, 1])
        with col_a:
            asset_input = st.number_input(
                "Average asset size",
                min_value=0.001,
                value=float(asset_value),
                key="asset_value",
                help="Average size of one file or object, such as one video master. It determines asset count and per-request charges.",
            )
        with col_b:
            asset_unit_input = st.selectbox(
                "Unit ", ["MB", "GB"], index=["MB", "GB"].index(asset_unit), key="asset_unit",
                help="Unit used for the average size of one asset.",
            )

        time_col_a, time_col_b = st.columns(2)
        with time_col_a:
            start_year = st.number_input(
                "Start year", min_value=2025, max_value=2500, value=initial.start_year,
                key="start_year_widget", help="Calendar year in which the archive is first written.",
            )
        with time_col_b:
            horizon = st.number_input(
                "Retention (years)", min_value=1, max_value=10_000,
                value=initial.horizon_years, key="horizon",
                help="Number of charged storage years, including the start year.",
            )

        finance_col_a, finance_col_b = st.columns(2)
        with finance_col_a:
            retrieval = st.number_input(
                "Annual retrieval (%)", min_value=0.0, max_value=10_000.0,
                value=float(initial.annual_retrieval_percent), step=0.25, key="retrieval",
                help="Expected share of the logical archive retrieved each year. 1% means reading 10 TB per year from a 1 PB archive.",
            )
        with finance_col_b:
            discount = st.number_input(
                "Discount rate (%)", min_value=0.0, max_value=99.0,
                value=float(initial.discount_rate_percent), step=0.25, key="discount",
                help="Real rate used to discount future payments to the storage start year. Use 0 for undiscounted costs.",
            )

        st.caption("Chart display")
        chart_col_a, chart_col_b = st.columns([1.4, 1])
        with chart_col_a:
            projection_end = st.number_input(
                "Outlook end year", min_value=int(start_year), max_value=2500,
                value=default_projection_end, key="projection_end",
                help="Last hypothetical archive start year shown in the start-year outlook.",
            )
        with chart_col_b:
            log_scale = st.toggle(
                "Log", value=default_log_scale, key="log_scale",
                help="Recommended when technologies differ by several orders of magnitude.",
            )

        st.caption("Storage technologies")
        with st.container(border=True):
            tech_col_a, tech_col_b = st.columns(2)
            with tech_col_a:
                tech_dna = st.checkbox(
                    "DNA", value="DNA" in initial.technologies, key="tech_dna",
                )
                tech_amazon = st.checkbox(
                    "Amazon", value="Amazon Deep Archive" in initial.technologies, key="tech_amazon",
                )
                tech_tape = st.checkbox(
                    "Tape", value="Tape On-premise" in initial.technologies, key="tech_tape",
                )
            with tech_col_b:
                tech_azure = st.checkbox(
                    "Azure", value="Azure Blob Archive" in initial.technologies, key="tech_azure",
                )
                tech_custom = st.checkbox(
                    "Custom", value="Custom storage" in initial.technologies, key="tech_custom",
                )

        with st.expander("DNA cost assumptions"):
            dna_cost_base_year = st.number_input(
                "DNA cost base year", min_value=2000, max_value=2500, value=initial.dna_cost_base_year,
                key="dna_cost_base_year",
                help="Year to which the editable synthesis and sequencing unit costs apply.",
            )
            dna_synthesis_cost = st.number_input(
                "Synthesis cost (USD/MB)", min_value=0.0, value=float(initial.dna_synthesis_cost_per_mb),
                format="%.6f", key="dna_synthesis_cost",
                help="Cost in the DNA cost base year to synthesize enough bases for 1 MB of logical data, before redundancy and indexing overhead.",
            )
            dna_sequencing_cost = st.number_input(
                "Sequencing cost (USD/MB)", min_value=0.0, value=float(initial.dna_sequencing_cost_per_mb),
                format="%.8f", key="dna_sequencing_cost",
                help="Cost in the DNA cost base year to sequence 1 MB of retrieved logical data.",
            )
            synthesis_decline = st.number_input(
                "Synthesis annual decline (%)", min_value=0.0, max_value=99.99,
                value=float(initial.synthesis_decline_percent), key="synthesis_decline",
                help="Percentage by which synthesis cost is assumed to fall each calendar year.",
            )
            sequencing_decline = st.number_input(
                "Sequencing annual decline (%)", min_value=0.0, max_value=99.99,
                value=float(initial.sequencing_decline_percent), key="sequencing_decline",
                help="Percentage by which sequencing cost is assumed to fall each calendar year.",
            )
            dna_durability = st.number_input(
                "DNA durability (years)", min_value=1, max_value=10_000,
                value=initial.dna_durability_years, key="dna_durability",
                help="Years before the archive must be synthesized again. No replacement occurs at the exact end of the horizon.",
            )

        with st.expander("Amazon Deep Archive assumptions"):
            st.caption("Price reference")
            amazon_base_year = st.number_input(
                "Amazon price base year", min_value=2000, max_value=2500,
                value=initial.amazon_price_base_year, key="amazon_base_year",
                help="Calendar year to which all Amazon prices below apply.",
            )
            amazon_decline = st.number_input(
                "Amazon annual price decline (%)", min_value=0.0, max_value=99.99,
                value=float(initial.amazon_decline_percent), key="amazon_decline",
                help="Annual reduction applied to Amazon request, retrieval, and storage prices.",
            )
            st.caption("Base-year prices")
            amazon_put_per_1000 = st.number_input(
                "Write requests (USD/1,000)", min_value=0.0,
                value=float(initial.amazon_put_usd_per_request * 1_000), format="%.6f",
                key="amazon_put_per_1000",
                help="Charge for 1,000 requests when the archive is initially written.",
            )
            amazon_restore_per_1000 = st.number_input(
                "Bulk restore requests (USD/1,000)", min_value=0.0,
                value=float(initial.amazon_bulk_restore_usd_per_request * 1_000), format="%.6f",
                key="amazon_restore_per_1000",
                help="Charge for 1,000 bulk restore-job requests. Asset size determines the request count.",
            )
            amazon_retrieval_per_tb = st.number_input(
                "Bulk data retrieval (USD/TB)", min_value=0.0,
                value=float(initial.amazon_bulk_retrieval_usd_per_mb * 1_000_000), format="%.6f",
                key="amazon_retrieval_per_tb",
                help="Capacity charge for retrieving one TB of archived data.",
            )
            amazon_storage_per_tb_month = st.number_input(
                "Storage (USD/TB/month)", min_value=0.0,
                value=float(initial.amazon_storage_usd_per_mb_month * 1_000_000), format="%.6f",
                key="amazon_storage_per_tb_month",
                help="Recurring monthly charge to retain one TB in Deep Archive.",
            )

        with st.expander("Azure Blob Archive assumptions"):
            st.caption("Price reference")
            azure_base_year = st.number_input(
                "Azure price base year", min_value=2000, max_value=2500,
                value=initial.azure_price_base_year, key="azure_base_year",
                help="Calendar year to which all Azure prices below apply.",
            )
            azure_decline = st.number_input(
                "Azure annual price decline (%)", min_value=0.0, max_value=99.99,
                value=float(initial.azure_decline_percent), key="azure_decline",
                help="Annual reduction applied to Azure request, retrieval, and storage prices.",
            )
            st.caption("Base-year prices")
            azure_write_per_1000 = st.number_input(
                "Write requests (USD/1,000)", min_value=0.0,
                value=float(initial.azure_write_usd_per_request * 1_000), format="%.6f",
                key="azure_write_per_1000",
                help="Charge for 1,000 requests when the archive is initially written.",
            )
            azure_read_per_1000 = st.number_input(
                "Read requests (USD/1,000)", min_value=0.0,
                value=float(initial.azure_read_usd_per_request * 1_000), format="%.6f",
                key="azure_read_per_1000",
                help="Charge for 1,000 retrieval requests. Asset size determines the request count.",
            )
            azure_retrieval_per_tb = st.number_input(
                "Data retrieval (USD/TB)", min_value=0.0,
                value=float(initial.azure_retrieval_usd_per_mb * 1_000_000), format="%.6f",
                key="azure_retrieval_per_tb",
                help="Capacity charge for retrieving one TB from the Archive tier.",
            )
            azure_storage_per_tb_month = st.number_input(
                "Storage (USD/TB/month)", min_value=0.0,
                value=float(initial.azure_storage_usd_per_mb_month * 1_000_000), format="%.6f",
                key="azure_storage_per_tb_month",
                help="Recurring monthly charge to retain one TB in the Archive tier.",
            )

        with st.expander("Tape on-premise assumptions"):
            st.caption("Price reference")
            tape_base_year = st.number_input(
                "Tape price base year", min_value=2000, max_value=2500,
                value=initial.tape_price_base_year, key="tape_base_year",
                help="Calendar year to which the tape media, hardware, and energy prices apply.",
            )
            tape_durability = st.number_input(
                "Tape durability (years)", min_value=1, max_value=1_000,
                value=initial.tape_durability_years, key="tape_durability",
                help="Years between complete tape media replacement writes.",
            )
            st.caption("Base-year prices")
            tape_media_per_tb = st.number_input(
                "Tape media (USD/TB)", min_value=0.0,
                value=float(initial.tape_media_usd_per_tb), format="%.6f",
                key="tape_media_per_tb",
                help="Media purchase cost for one TB, charged on the initial write and every replacement.",
            )
            tape_hardware_per_tb = st.number_input(
                "Hardware (USD/TB)", min_value=0.0,
                value=float(initial.tape_hardware_usd_per_tb), format="%.6f",
                key="tape_hardware_per_tb",
                help="Hardware cost allocated per TB and amortized over the selected durability period.",
            )
            tape_energy_per_tb_year = st.number_input(
                "Energy (USD/TB/year)", min_value=0.0,
                value=float(initial.tape_energy_usd_per_tb_year), format="%.6f",
                key="tape_energy_per_tb_year",
                help="Annual energy cost to retain one TB in the tape system.",
            )
            st.caption("Annual price declines")
            tape_media_decline = st.number_input(
                "Tape media decline (%)", min_value=0.0, max_value=99.99,
                value=float(initial.tape_media_decline_percent), key="tape_media_decline",
                help="Annual reduction applied to tape media purchase prices.",
            )
            tape_hardware_decline = st.number_input(
                "Tape hardware decline (%)", min_value=0.0, max_value=99.99,
                value=float(initial.tape_hardware_decline_percent), key="tape_hardware_decline",
                help="Annual reduction applied to amortized tape hardware costs.",
            )
            tape_energy_decline = st.number_input(
                "Tape energy decline (%)", min_value=0.0, max_value=99.99,
                value=float(initial.tape_energy_decline_percent), key="tape_energy_decline",
                help="Annual reduction applied to tape energy costs.",
            )

        with st.expander("Custom storage assumptions"):
            custom_name = st.text_input(
                "Display name", value=initial.custom_storage_name, key="custom_name",
                help="Name used for the custom technology in charts, tables, and downloads.",
            )
            custom_base_year = st.number_input(
                "Price base year", min_value=2000, max_value=2500, value=initial.custom_cost_base_year,
                key="custom_base_year", help="Year to which all custom prices apply.",
            )
            custom_write_tb = st.number_input(
                "Initial write cost (USD/TB)", min_value=0.0, value=float(initial.custom_write_cost_per_tb),
                key="custom_write_tb", help="Capacity-based charge to write or replace one TB.",
            )
            custom_write_asset = st.number_input(
                "Write request cost (USD/asset)", min_value=0.0,
                value=float(initial.custom_write_cost_per_asset), key="custom_write_asset",
                help="Per-file or per-object charge applied when the archive is written or replaced.",
            )
            custom_storage_tb_year = st.number_input(
                "Annual storage cost (USD/TB)", min_value=0.0,
                value=float(initial.custom_storage_cost_per_tb_year), key="custom_storage_tb_year",
                help="Recurring cost to retain one TB for one year.",
            )
            custom_retrieval_tb = st.number_input(
                "Retrieval cost (USD/TB)", min_value=0.0,
                value=float(initial.custom_retrieval_cost_per_tb), key="custom_retrieval_tb",
                help="Capacity-based charge for each TB retrieved.",
            )
            custom_retrieval_asset = st.number_input(
                "Retrieval request cost (USD/asset)", min_value=0.0,
                value=float(initial.custom_retrieval_cost_per_asset), key="custom_retrieval_asset",
                help="Per-file or per-object charge for the expected assets retrieved each year.",
            )
            custom_decline = st.number_input(
                "Annual price decline (%)", min_value=0.0, max_value=99.99,
                value=float(initial.custom_decline_percent), key="custom_decline",
                help="Annual percentage reduction applied to every custom price.",
            )
            custom_replacement = st.number_input(
                "Replacement interval (years)", min_value=0, max_value=10_000,
                value=initial.custom_replacement_years, key="custom_replacement",
                help="Years between complete rewrites. Use 0 for a service with no replacement writes.",
            )

        submitted = st.form_submit_button("Calculate scenario", type="primary", width="stretch")

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
if not technologies:
    st.error("Select at least one storage technology.")
    st.stop()

archive_multiplier = {"TB": 1, "PB": 1000, "EB": 1_000_000}[archive_unit_input]
asset_multiplier = {"MB": 1, "GB": 1000}[asset_unit_input]
try:
    scenario = Scenario(
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

if submitted:
    params = scenario.to_query_params()
    params.update({"projection_end": str(projection_end), "log_scale": str(log_scale)})
    st.query_params.from_dict(params)

result = _cached_simulation(scenario)
projection = _cached_projection(scenario, int(projection_end))
dna_curve_end = max(
    scenario.dna_cost_base_year,
    min(2500, scenario.start_year + scenario.horizon_years - 1),
)
dna_costs = _cached_dna_costs(scenario, dna_curve_end)
use_present_value = scenario.discount_rate_percent > 0
value_column = "present_value_usd" if use_present_value else "total_cost_usd"

st.title("DNA Storage Cost Explorer")
st.caption("Lifecycle scenarios for archival DNA, cloud, and tape storage")
st.markdown(
    f"""
    <div class="model-strip">
    Model v{result.metadata['model_version']} &nbsp;|&nbsp; {result.metadata['currency']}
    &nbsp;|&nbsp; assumptions reviewed {result.metadata['last_reviewed']}
    &nbsp;|&nbsp; {result.metadata['disclaimer']}
    </div>
    """,
    unsafe_allow_html=True,
)

totals = result.totals.sort_values(value_column)
cheapest = totals.iloc[0]
dna_rows = totals[totals["technology"] == "DNA"]
dna_total = float(dna_rows.iloc[0][value_column]) if not dna_rows.empty else None

metric_columns = st.columns(4)
metric_columns[0].metric("Archive", f"{scenario.archive_size_tb:,.3g} TB")
metric_columns[1].metric("Assets", f"{scenario.number_of_assets:,.0f}")
metric_columns[2].metric(
    f"Lowest: {cheapest['technology']}", _money(float(cheapest[value_column]))
)
if dna_total is None:
    metric_columns[3].metric("Annual retrieval", f"{scenario.annual_retrieval_percent:,.2f}%")
else:
    lowest_cost = float(cheapest[value_column])
    if lowest_cost > 0:
        ratio = dna_total / lowest_cost
        metric_columns[3].metric(f"DNA cost ({ratio:,.2g}x lowest)", _money(dna_total))
    else:
        metric_columns[3].metric("DNA lifecycle cost", _money(dna_total))

overview_tab, outlook_tab, dna_cost_tab, assumptions_tab = st.tabs(
    ["Lifecycle", "Start-year outlook", "DNA unit costs", "Assumptions"]
)

with overview_tab:
    st.caption(
        "Cumulative lifecycle cost for one archive opened in the selected start year. "
        "Each line adds initial writes, replacement writes, annual storage/operation, and expected retrieval."
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
    st.caption(
        "The component chart separates undiscounted write and replacement costs, retrieval costs, "
        "and recurring storage or operating costs over the full retention horizon."
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
    st.caption(
        "The same archive workload and retention horizon are recalculated for every possible storage "
        "start year. A crossover is the first start year for which DNA's lifecycle cost is no greater "
        "than the comparison technology."
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
    st.caption(
        f"Both curves begin with the editable {scenario.dna_cost_base_year} unit costs and apply the "
        f"selected annual decline rates through {dna_curve_end}. They are unit-cost assumptions, "
        "not lifecycle totals."
    )
    chart_columns = st.columns(2)
    with chart_columns[0]:
        synthesis_title = "DNA synthesis cost trajectory"
        st.plotly_chart(
            dna_unit_cost_chart(
                dna_costs,
                "synthesis_cost_usd_per_mb",
                synthesis_title,
                "#b44737",
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
                "#b44737",
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
    st.subheader("Scenario contract")
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
    st.dataframe(contract, hide_index=True, width="stretch")

    st.subheader("Cost scope")
    st.write(
        "Included: DNA synthesis and sequencing, cloud write/retrieval/storage charges, "
        "tape media/hardware/energy assumptions, and the selected custom capacity/request charges. "
        "Excluded: labor, cloud egress, taxes, "
        "retrieval latency, minimum-storage penalties, facilities, and unmodeled migrations."
    )
    st.subheader("Sources")
    for source in assumptions["sources"].values():
        st.markdown(f"- [{source['label']}]({source['url']})")
