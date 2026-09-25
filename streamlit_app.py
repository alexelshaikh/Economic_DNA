from __future__ import annotations

import html
import json
from pathlib import Path
import pandas as pd
import streamlit as st

# Form context for the main-area widgets (see the comment before the cost
# rail): pinned to streamlit==1.63.0, whose form internals these touch.
from streamlit.delta_generator_singletons import get_dg_singleton_instance
from streamlit.elements.lib.form_utils import FormData as _FormData

from economic_dna import (
    PRESET_SCENARIOS,
    Scenario,
    dna_cost_sensitivity,
    find_breakeven_synthesis_cost,
    find_crossover_years,
    find_lifecycle_crossovers,
    load_observed_sequencing_costs,
    simulate_dna_unit_costs,
    simulate_dna_uncertainty_band,
    simulate_scenario,
    simulate_start_years,
    synthesis_historical_trend,
)
from economic_dna.assumptions import load_assumptions
from economic_dna.visualization import (
    breakdown_chart,
    breakeven_chart,
    dna_unit_cost_chart,
    format_display_number,
    lifecycle_chart,
    palette_for,
    projection_chart,
    sensitivity_chart,
)


st.set_page_config(
    page_title="DNA Storage Cost Explorer",
    layout="wide",
    # "auto" keeps the sidebar expanded on desktop and collapsed on phones,
    # where it renders as a full-screen overlay (see the mobile media query).
    initial_sidebar_state="auto",
)

_main_dg = get_dg_singleton_instance().main_dg
# The app temporarily attaches form data to the main dg so header/panel
# controls can join the sidebar form. A callback-triggered rerun can interrupt
# before the normal cleanup point, so each run must begin outside that form.
_main_dg._form_data = None

@st.cache_data(show_spinner=False, max_entries=1)
def _stylesheet(modified_ns: int) -> str:
    return (Path(__file__).parent / "economic_dna" / "ui.css").read_text(encoding="utf-8")


@st.cache_data(show_spinner=False, max_entries=1)
def _form_script(modified_ns: int) -> str:
    return (Path(__file__).parent / "economic_dna" / "form_controls.js").read_text(encoding="utf-8")


st.html(f"<style>{_stylesheet((Path(__file__).parent / 'economic_dna' / 'ui.css').stat().st_mtime_ns)}</style>")

# Theme: an explicit choice (URL parameter or the sidebar toggle) wins; on a
# first visit Python starts from light and the component below redirects to
# dark when the system prefers it. The rendered marker drives the CSS via
# html:has(...), so every rerun reflects the resolved theme immediately.
theme = str(st.query_params.get("theme", ""))
if theme not in ("light", "dark"):
    theme = st.session_state.get("theme", "light")
    if theme not in ("light", "dark"):
        theme = "light"
st.session_state["theme"] = theme

st.markdown(
    f'<div class="theme-marker theme-{theme}" hidden></div>',
    unsafe_allow_html=True,
)


def _apply_system_dark() -> None:
    st.session_state["theme"] = "dark"
    st.query_params["theme"] = "dark"


# Hidden button the theme-sync iframe clicks programmatically when the system
# prefers dark and no explicit theme is set, so the server learns the
# preference (its iframe is sandboxed against navigating the top window).
st.button(
    label=None,
    key="theme_auto_dark",
    on_click=_apply_system_dark,
    icon=":material/dark_mode:",
)

_THEME_SYNC_JS = """
<script>
(() => {
  const params = new URLSearchParams(parent.location.search);
  if (params.get("theme")) return;
  if (!parent.matchMedia("(prefers-color-scheme: dark)").matches) return;
  // Apply the dark classes to the theme marker immediately: the whole CSS
  // dark palette hangs off html:has(.theme-marker.theme-dark), so flipping
  // the class re-themes the page in the same frame instead of waiting for
  // the server rerun (which is what made the tabs flash light for seconds).
  const tryApply = () => {
    const marker = parent.document.querySelector(".theme-marker");
    if (!marker) return false;
    marker.classList.remove("theme-light");
    marker.classList.add("theme-dark");
    return true;
  };
  const tryClick = () => {
    // Still click the hidden button so the SERVER also learns the theme and
    // the rerun re-renders charts and markers in dark consistently.
    const btn = parent.document.querySelector(".st-key-theme_auto_dark button");
    if (!btn) return false;
    btn.click();
    return true;
  };
  tryApply();
  if (tryClick()) return;
  let tries = 0;
  const timer = setInterval(() => {
    tries += 1;
    if (tryClick() || tries > 100) clearInterval(timer);
  }, 150);
})();
</script>
"""

st.iframe(_THEME_SYNC_JS, width=1, height=1)

# Click-away close for the cost panel: a one-shot iframe attaches a document
# listener on the parent page that selects the rail's ✕ radio option whenever
# a click lands outside the rail and the panel. The flag guards against
# duplicate listeners (the iframe is re-created on every rerun).
_COST_CLICK_AWAY_JS = """
<script>
(() => {
  if (parent.document.__costClickAway) return;
  parent.document.__costClickAway = true;
  let downInside = false;
  const inside = (target) => {
    const rail = parent.document.querySelector(".st-key-cost-rail");
    const panel = parent.document.querySelector(".st-key-cost-panel");
    return (rail && rail.contains(target)) || (panel && panel.contains(target));
  };
  // A drag that STARTS inside the rail or panel must never close the panel,
  // even if it is released outside (the click then fires on the common
  // ancestor, which would look like an outside click).
  parent.document.addEventListener("pointerdown", (event) => {
    downInside = inside(event.target);
  }, true);
  parent.document.addEventListener("click", (event) => {
    const closePanel = () => {
      const closeLabel = parent.document.querySelector(
        '.st-key-cost-rail label:has(input[value="0"])'
      );
      // Deferred: a click synthesized inside another event's capture phase
      // is swallowed, so close after the original click has fully completed.
      if (closeLabel) setTimeout(() => closeLabel.click(), 0);
    };
    const targetInside = inside(event.target);
    if (downInside) {
      // Released outside after pressing inside: keep the panel open.
      if (!targetInside) return;
      // Pressed and released on the already-selected model tab: toggle close.
      const label = event.target.closest(".st-key-cost-rail [role='radiogroup'] label");
      if (label) {
        const input = label.querySelector("input");
        if (input && input.checked && input.value !== "0") closePanel();
      }
      return;
    }
    if (targetInside) return;
    closePanel();
  }, true);
})();
</script>
"""
st.iframe(_COST_CLICK_AWAY_JS, width=1, height=1)

# Close the mobile sheet through its native control, preserving the session.
_SIDEBAR_CLOSE_JS = """
<script>
(() => {
  if (parent.innerWidth > 640) return;
  const sidebar = parent.document.querySelector('[data-testid="stSidebar"][aria-expanded="true"]');
  const close = sidebar?.querySelector('[data-testid="stSidebarCollapseButton"] button');
  if (close) close.click();
})();
</script>
"""

def _toggle_theme() -> None:
    new_theme = "dark" if st.session_state["theme"] == "light" else "light"
    st.session_state["theme"] = new_theme
    st.query_params["theme"] = new_theme


_main_dg._form_data = _FormData("scenario_form")
with st.container(key="theme-toggle-anchor"):
    st.form_submit_button(
        label=None,
        icon=":material/dark_mode:" if theme == "light" else ":material/light_mode:",
        key="theme_toggle",
        help="Switch to dark mode" if theme == "light" else "Switch to light mode",
        on_click=_toggle_theme,
    )
_main_dg._form_data = None

# Cost-assumption models for the right-edge rail. Each key is a widget prefix
# (matching the scenario parameter names), and the label is the tab handle.
COST_MODELS = [
    ("dna", "DNA"),
    ("amazon", "Amazon"),
    ("azure", "Azure"),
    ("tape", "Tape"),
    ("custom", "Custom"),
]

# Widget keys per cost model, used by the per-model "Reset to defaults"
# buttons: only the listed keys are restored to the paper baseline.
MODEL_WIDGET_KEYS = {
    "dna": [
        "dna_cost_base_year", "dna_synthesis_cost", "dna_sequencing_cost",
        "synthesis_decline", "sequencing_decline", "dna_durability",
    ],
    "amazon": [
        "amazon_base_year", "amazon_put_per_1000", "amazon_restore_per_1000",
        "amazon_retrieval_per_tb", "amazon_storage_per_tb_month", "amazon_decline",
    ],
    "azure": [
        "azure_base_year", "azure_write_per_1000", "azure_read_per_1000",
        "azure_retrieval_per_tb", "azure_storage_per_tb_month", "azure_decline",
    ],
    "tape": [
        "tape_base_year", "tape_durability", "tape_media_per_tb",
        "tape_hardware_per_tb", "tape_energy_per_tb_year", "tape_media_decline",
        "tape_hardware_decline", "tape_energy_decline",
    ],
    "custom": [
        "custom_name", "custom_base_year", "custom_write_tb", "custom_write_asset",
        "custom_storage_tb_year", "custom_retrieval_tb", "custom_retrieval_asset",
        "custom_decline", "custom_replacement",
    ],
}



SIDEBAR_WIDGET_KEYS = [
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
    "projection_end",
    "log_scale",
]


PRESET_WIDGET_KEYS = [
    "archive_value",
    "archive_unit",
    "asset_value",
    "asset_unit",
    "retrieval",
    "horizon",
]


PRESET_BUTTON_NAMES = [name for name in PRESET_SCENARIOS if name != "Paper baseline"]


NUMBER_INPUT_FORMATS = {
    "retrieval": "%.4g",
    "dna_synthesis_cost": "%.10g",
    "dna_sequencing_cost": "%.10g",
    "amazon_put_per_1000": "%.6f",
    "amazon_restore_per_1000": "%.6f",
    "amazon_retrieval_per_tb": "%.6f",
    "amazon_storage_per_tb_month": "%.6f",
    "azure_write_per_1000": "%.6f",
    "azure_read_per_1000": "%.6f",
    "azure_retrieval_per_tb": "%.6f",
    "azure_storage_per_tb_month": "%.6f",
    "tape_media_per_tb": "%.6f",
    "tape_hardware_per_tb": "%.6f",
    "tape_energy_per_tb_year": "%.6f",
}


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
    try:
        final_year = int(query.get("projection_end", 2350))
    except (TypeError, ValueError, OverflowError):
        final_year = 2350
    return min(2500, max(scenario.start_year, final_year))


def _query_mapping() -> dict[str, str]:
    return {key: str(value) for key, value in st.query_params.items()}


def _initial_scenario() -> Scenario:
    try:
        scenario = Scenario.from_mapping(_query_mapping())
        _cached_simulation(scenario)
        return scenario
    except (TypeError, ValueError, OverflowError):
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


def _snapshot_widgets() -> dict[str, bool | float | int | str]:
    widgets = {key: st.session_state[key] for key in WIDGET_KEYS}
    widgets["projection_end"] = max(widgets["start_year_widget"], widgets["projection_end"])
    return widgets


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


def _archive_size_label(value_tb: float) -> str:
    # Plain ",.3g" silently switches to "1e+03"-style scientific notation for
    # any round value >= 1000 (any archive at or above 1 PB); the preset
    # scenarios made this visible immediately, so it needed a real fix.
    absolute = abs(value_tb)
    for threshold, suffix in ((1e12, "T"), (1e9, "B"), (1e6, "M"), (1e3, "K")):
        if absolute >= threshold:
            return f"{value_tb / threshold:,.3g}{suffix}"
    return f"{value_tb:,.3g}"


def _section_intro(title: str, description: str) -> None:
    st.markdown(
        f'<div class="tab-intro"><h2>{title}</h2><p>{description}</p></div>',
        unsafe_allow_html=True,
    )


def _themed_table(columns: list[str], rows: list[list[str]]) -> None:
    header = "".join(f"<th>{html.escape(column)}</th>" for column in columns)
    body = "".join(
        "<tr>" + "".join(f"<td>{html.escape(str(cell))}</td>" for cell in row) + "</tr>"
        for row in rows
    )
    st.markdown(
        f'<div class="contract-card"><table class="contract-table">'
        f"<thead><tr>{header}</tr></thead><tbody>{body}</tbody></table></div>",
        unsafe_allow_html=True,
    )


def _chart_downloads(
    key: str,
    csv_data: pd.DataFrame,
    filename_base: str,
    chart_key: str,
) -> None:
    st.markdown('<div class="export-label">Export chart</div>', unsafe_allow_html=True)
    columns = st.columns(3)
    columns[0].download_button(
        "CSV",
        lambda: csv_data.to_csv(index=False).encode("utf-8"),
        f"{filename_base}.csv",
        "text/csv",
        key=f"download_{key}_csv",
        icon=":material/download:",
        width="stretch",
        help="Download the data shown in this graph.",
        on_click="ignore",
    )
    columns[1].button(
        "PNG",
        key=f"download_{key}_png",
        icon=":material/download:",
        width="stretch",
        help="Download this graph as a high-resolution PNG image.",
    )
    columns[2].button(
        "SVG",
        key=f"download_{key}_svg",
        icon=":material/download:",
        width="stretch",
        help="Download this graph as an editable vector image.",
    )
    _bind_live_plotly_image_downloads(key, chart_key, filename_base)


def _bind_live_plotly_image_downloads(key: str, chart_key: str, filename_base: str) -> None:
    config_json = json.dumps(
        {
            "chartKey": chart_key,
            "filename": filename_base,
            "buttons": {
                "png": f"download_{key}_png",
                "svg": f"download_{key}_svg",
            },
        }
    )
    st.iframe(
        f"""
        <script>
        (() => {{
          const config = {config_json};
          const findByKey = (key) => parent.document.querySelector(`.st-key-${{key}}`);
          const findGraph = () => {{
            const root = findByKey(config.chartKey);
            return root ? root.querySelector(".js-plotly-plot") : null;
          }};
          const plotlyApi = () => parent.Plotly || parent.document.defaultView.Plotly;
          const bind = () => {{
            const graph = findGraph();
            const Plotly = plotlyApi();
            let bound = 0;
            for (const [format, buttonKey] of Object.entries(config.buttons)) {{
              const buttonRoot = findByKey(buttonKey);
              const button = buttonRoot ? buttonRoot.querySelector("button") : null;
              if (!button || button.dataset.livePlotlyExport === config.chartKey + format) {{
                if (button) bound += 1;
                continue;
              }}
              button.dataset.livePlotlyExport = config.chartKey + format;
              button.addEventListener("click", (event) => {{
                const liveGraph = findGraph();
                const livePlotly = plotlyApi();
                if (!liveGraph || !livePlotly || !livePlotly.downloadImage) return;
                event.preventDefault();
                event.stopPropagation();
                event.stopImmediatePropagation();
                const bounds = liveGraph.getBoundingClientRect();
                const width = Math.max(320, Math.round(bounds.width || liveGraph._fullLayout?.width || 1200));
                const height = Math.max(240, Math.round(bounds.height || liveGraph._fullLayout?.height || 500));
                livePlotly.downloadImage(liveGraph, {{
                  format,
                  filename: config.filename,
                  width,
                  height,
                  scale: format === "png" ? 2 : 1,
                }});
              }}, true);
              bound += 1;
            }}
            return Boolean(graph && Plotly && Plotly.downloadImage && bound === 2);
          }};
          if (bind()) return;
          let tries = 0;
          const timer = setInterval(() => {{
            tries += 1;
            if (bind() || tries > 100) clearInterval(timer);
          }}, 100);
        }})();
        </script>
        """,
        width=1,
        height=1,
    )


def _bind_copy_link_button(key: str) -> None:
    """Copies the browser's current address bar URL (already synced to the
    committed scenario by Calculate) to the clipboard. Runs in the parent
    document, like the other iframe-JS helpers here, so the Clipboard API
    call is made by the top-level page rather than the sandboxed iframe."""
    st.iframe(
        f"""
        <script>
        (() => {{
          const bind = () => {{
            const root = parent.document.querySelector(".st-key-{key}");
            const buttons = root ? Array.from(root.querySelectorAll("button")) : [];
            // Streamlit can render a hidden 0x0 duplicate around tooltip state;
            // bind the button the user can actually click.
            const button = buttons.find((candidate) => {{
              const rect = candidate.getBoundingClientRect();
              const style = parent.getComputedStyle(candidate);
              return rect.width > 0 && rect.height > 0
                && style.display !== "none"
                && style.visibility !== "hidden";
            }}) || null;
            if (!button || button.dataset.copyLinkBound) return Boolean(button);
            button.dataset.copyLinkBound = "1";
            const original = button.innerHTML;
            const flash = (text) => {{
              button.innerHTML = original;
              button.append(Object.assign(parent.document.createElement("span"), {{ textContent: " " + text }}));
              setTimeout(() => {{ button.innerHTML = original; }}, 1600);
            }};
            // Legacy fallback for browsers/contexts that block the async
            // Clipboard API (older browsers, non-secure origins, some
            // automation and enterprise policies) -- still gives feedback
            // either way instead of a click that silently does nothing.
            const legacyCopy = () => {{
              const input = parent.document.createElement("input");
              input.value = parent.location.href;
              input.style.position = "fixed";
              input.style.opacity = "0";
              parent.document.body.appendChild(input);
              input.focus();
              input.select();
              let ok = false;
              try {{ ok = parent.document.execCommand("copy"); }} catch (err) {{ ok = false; }}
              input.remove();
              flash(ok ? "Copied!" : "Copy failed - copy the address bar");
            }};
            button.addEventListener("click", (event) => {{
              event.preventDefault();
              event.stopPropagation();
              event.stopImmediatePropagation();
              if (!parent.navigator.clipboard || !parent.navigator.clipboard.writeText) {{
                legacyCopy();
                return;
              }}
              parent.navigator.clipboard.writeText(parent.location.href).then(
                () => flash("Copied!"),
                () => legacyCopy(),
              );
            }}, true);
            return true;
          }};
          if (bind()) return;
          let tries = 0;
          const timer = setInterval(() => {{
            tries += 1;
            if (bind() || tries > 100) clearInterval(timer);
          }}, 100);
        }})();
        </script>
        """,
        width=1,
        height=1,
    )


def _plot_config(filename: str) -> dict:
    return {
        "displaylogo": False,
        "toImageButtonOptions": {"format": "png", "filename": filename, "scale": 2},
    }


def _render_chart(figure, *, key: str, filename: str) -> None:
    st.plotly_chart(
        figure,
        key=key,
        width="stretch",
        config=_plot_config(filename),
    )


@st.cache_data(show_spinner=False, max_entries=8, ttl=600)
def _cached_simulation(scenario: Scenario):
    return simulate_scenario(scenario)


@st.cache_data(show_spinner=False, max_entries=8, ttl=600)
def _cached_projection(scenario: Scenario, final_year: int) -> pd.DataFrame:
    return simulate_start_years(scenario, final_year)


@st.cache_data(show_spinner=False, max_entries=32)
def _cached_dna_costs(scenario: Scenario, final_year: int) -> pd.DataFrame:
    return simulate_dna_unit_costs(scenario, final_year)


@st.cache_data(show_spinner=False, max_entries=32)
def _cached_uncertainty_band(scenario: Scenario, use_present_value: bool) -> pd.DataFrame:
    return simulate_dna_uncertainty_band(scenario, use_present_value)


@st.cache_data(show_spinner=False, max_entries=32)
def _cached_sensitivity(scenario: Scenario, use_present_value: bool) -> pd.DataFrame:
    return dna_cost_sensitivity(scenario, use_present_value)


@st.cache_data(show_spinner=False, max_entries=32)
def _cached_breakeven(scenario: Scenario, use_present_value: bool) -> pd.DataFrame:
    return find_breakeven_synthesis_cost(scenario, use_present_value)


@st.cache_data(show_spinner=False, max_entries=1)
def _cached_synthesis_history(base_year: int) -> pd.DataFrame:
    return synthesis_historical_trend(2000, base_year)


@st.cache_data(show_spinner=False, max_entries=1)
def _cached_observed_sequencing_costs() -> pd.DataFrame:
    return load_observed_sequencing_costs()


initial = _initial_scenario()
query = _query_mapping()
initial_widgets = _widget_state_from_scenario(initial, query)
baseline_widgets = _widget_state_from_scenario(Scenario())


def _reset_widget_keys(widget_keys: list[str]) -> None:
    for widget_key in widget_keys:
        st.session_state[widget_key] = baseline_widgets[widget_key]


def _apply_preset(preset_name: str) -> None:
    # Presets set workload fields; the form submits other pending edits too.
    preset_widgets = _widget_state_from_scenario(PRESET_SCENARIOS[preset_name], _query_mapping())
    for widget_key in PRESET_WIDGET_KEYS:
        value = preset_widgets[widget_key]
        st.session_state[widget_key] = value


def _preset_button_label(preset_name: str) -> str:
    scenario = PRESET_SCENARIOS[preset_name]
    archive_value, archive_unit, asset_value, asset_unit = _display_units(scenario)
    archive = f"{archive_value:g} {archive_unit}"
    asset = f"{asset_value:g} {asset_unit}/asset"
    retrieval = f"{scenario.annual_retrieval_percent:g}%/yr"
    retention = f"{scenario.horizon_years:g} yrs"
    return f"{preset_name}\n{archive} | {asset} | {retrieval} | {retention}"


def _bind_form_controls(committed: dict) -> None:
    # Keep Python callbacks as a fallback; the browser normally applies these
    # input-only actions to the existing form without submitting it.
    actions = {
        "global_reset": {"values": baseline_widgets},
        "sidebar_reset": {"values": {key: baseline_widgets[key] for key in SIDEBAR_WIDGET_KEYS}},
    }
    for model, keys in MODEL_WIDGET_KEYS.items():
        actions[f"reset_{model}"] = {"values": {key: baseline_widgets[key] for key in keys}}
    for index, name in enumerate(PRESET_SCENARIOS):
        if name in PRESET_BUTTON_NAMES:
            values = _widget_state_from_scenario(PRESET_SCENARIOS[name])
            actions[f"preset_{index}"] = {"values": {key: values[key] for key in PRESET_WIDGET_KEYS}}
    for key in ("dna_synthesis_cost", "dna_sequencing_cost"):
        for suffix, factor in (("div10", 0.1), ("mul10", 10)):
            actions[f"{key}_{suffix}"] = {"scale": key, "factor": factor}
    displayed = {}
    for key, value in committed.items():
        if isinstance(value, float):
            precision = NUMBER_INPUT_FORMATS.get(key, "%.2f")
            displayed[key] = float(precision % value)
        else:
            displayed[key] = value
    config = json.dumps({"committed": committed, "displayed": displayed, "actions": actions})
    script_path = Path(__file__).parent / "economic_dna" / "form_controls.js"
    script = _form_script(script_path.stat().st_mtime_ns)
    # Escape user-entered names before embedding JSON, so they remain data,
    # never executable HTML. Native st.html keeps handlers in the page realm.
    config = config.replace("<", r"\u003c")
    st.html(
        f'<span class="form-controls-marker" hidden></span>'
        f"<script>(() => {{ const config = {config};\n{script}\n }})();</script>",
        unsafe_allow_javascript=True,
    )


for widget_key, default_value in initial_widgets.items():
    st.session_state.setdefault(widget_key, default_value)
st.session_state["projection_end"] = max(
    st.session_state["start_year_widget"], st.session_state["projection_end"]
)
# The graphs only follow the last calculated inputs: the initial load
# counts as the first calculation. Every widget lives inside one form, so
# edits never trigger reruns — the three Calculate buttons are form submits
# and are the only way to commit changes.
st.session_state.setdefault("committed_widgets", dict(initial_widgets))

# All inputs live in one form: Streamlit batches form widgets client-side and
# reruns the script only on a submit, so editing parameters never re-renders
# the buttons or the charts. The form block lives in the sidebar; the panel
# and header Calculate join it via the main dg's form data below.
with st.sidebar:
    with st.form("scenario_form", border=False, enter_to_submit=False):
        with st.expander("Example scenarios", expanded=False):
            with st.container(key="sidebar-presets"):
                for index, preset_name in enumerate(PRESET_SCENARIOS):
                    if preset_name not in PRESET_BUTTON_NAMES:
                        continue
                    st.form_submit_button(
                        _preset_button_label(preset_name),
                        key=f"preset_{index}",
                        on_click=_apply_preset,
                        args=(preset_name,),
                        width="stretch",
                    )
        with st.container(key="sidebar-reset-btn"):
            st.form_submit_button(
                "Reset workload",
                key="sidebar_reset",
                on_click=_reset_widget_keys,
                args=(SIDEBAR_WIDGET_KEYS,),
                width="stretch",
                icon=":material/restart_alt:",
            )
        with st.container(key="sidebar-advanced-toggle"):
            st.checkbox(
                "Show advanced cost assumptions",
                key="show_advanced",
                help="Reveals price base years, durability, and replacement cycles in the cost "
                "assumption panels on the right edge. Unit costs and decline rates -- the inputs "
                "that drive the projections -- are always visible.",
            )
        input_column, action_column = st.columns([6, 1], gap="small")
        with input_column:
            st.markdown(
                """
                <div class="sidebar-kicker">Scenario builder</div>
                <div class="sidebar-title">Model inputs</div>
                """,
                unsafe_allow_html=True,
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
                archive_unit_input = st.radio(
                    "Unit", ["TB", "PB", "EB"], key="archive_unit", horizontal=True,
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
                asset_unit_input = st.radio(
                    "Unit ", ["MB", "GB"], key="asset_unit", horizontal=True,
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
                    step=0.25, key="retrieval", format=NUMBER_INPUT_FORMATS["retrieval"],
                    help="Expected share of the logical archive retrieved each year. 1% means reading 10 TB per year from a 1 PB archive.",
                )
            with finance_col_b:
                discount = st.number_input(
                    "Discount rate (%)", min_value=0.0, max_value=99.0,
                    step=0.25, key="discount",
                    help=(
                        "Annual real discount rate for present-value results. "
                        "Use 0% to count every future dollar at face value. "
                        "Use a higher rate, such as 2-5%, if future costs should count less "
                        "than costs paid near the start year."
                    ),
                )
    
            st.markdown('<div class="sidebar-section">Display</div>', unsafe_allow_html=True)
            chart_col_a, chart_col_b = st.columns([1.4, 1])
            with chart_col_a:
                projection_end = st.number_input(
                    "Outlook end year", min_value=2025, max_value=2500,
                    key="projection_end",
                    help="Final archive start year included in the start-year outlook chart.",
                )
            with chart_col_b:
                st.markdown('<div class="log-scale-label">Log</div>', unsafe_allow_html=True)
                log_scale = st.toggle(
                    label="Log scale",
                    label_visibility="collapsed",
                    key="log_scale",
                    help="Recommended when technologies differ by several orders of magnitude.",
                )
    
            st.markdown('<div class="sidebar-section">Technologies</div>', unsafe_allow_html=True)
            with st.container():
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
    
        with action_column:
            with st.container(key="scenario-action-rail"):
                # The mobile bottom bar; hidden on desktop, where the wide
                # centered header button is the primary Calculate.
                calculate_scenario = st.form_submit_button(
                    "Calculate",
                    type="primary",
                    key="calculate_scenario",
                    width="stretch",
                )
    
# The panel and header Calculate join the sidebar's form. st.form cannot
# wrap both containers (a form is a single block), so the form id is attached
# to the main dg directly — the same mechanism the form block uses on itself.
# Widgets added through `st.foo` calls read this dg's form data, so the panel
# inputs batch with the sidebar's and only a Calculate submit reruns anything.
# Cost-assumption rail and panel: slim vertical model tabs on the right
# edge of the page (a bottom strip on phones). The tabs are pure-HTML radio
# labels — opening, closing, and switching never trigger a script rerun, so
# the buttons and charts stay untouched. All five model blocks stay mounted
# and CSS shows only the checked one, so edits survive closing and switching.
# Wide centered Calculate button in the top header: it stays visible while
# the sidebar or a cost panel is open, so it works after either kind of edit.
# (Hidden on phones, where the sheet bar and the panel bar cover the flows.)
_main_dg._form_data = _FormData("scenario_form")
with st.container(key="global-reset-btn"):
    st.form_submit_button(
        label=None,
        key="global_reset",
        icon=":material/restart_alt:",
        help="Reset all inputs to the paper baseline",
        on_click=_reset_widget_keys,
        args=(WIDGET_KEYS,),
    )

_main_dg._form_data = _FormData("scenario_form")


def _render_model_reset(model_key: str) -> None:
    with st.container(key=f"model-reset-{model_key}"):
        st.form_submit_button(
            "Reset to defaults",
            key=f"reset_{model_key}",
            on_click=_reset_widget_keys,
            args=(MODEL_WIDGET_KEYS[model_key],),
            icon=":material/refresh:",
        )


def _scale_widget_value(widget_key: str, factor: float) -> None:
    # Same before-the-next-rerun timing as _reset_widget_keys: safe to write
    # even though this is a form widget.
    st.session_state[widget_key] = st.session_state[widget_key] * factor


def _render_order_of_magnitude_steppers(widget_key: str) -> None:
    """A cost spanning many orders of magnitude (DNA synthesis cost runs from
    a fraction of a cent to tens of thousands of dollars per MB) is tedious
    to explore by typing digits one at a time. These buttons jump by a full
    decade in either direction."""
    with st.container(key=f"steppers-{widget_key}"):
        step_columns = st.columns(2)
        step_columns[0].form_submit_button(
            "÷10",
            key=f"{widget_key}_div10",
            on_click=_scale_widget_value,
            args=(widget_key, 0.1),
            width="stretch",
            help="Divide this value by 10.",
        )
        step_columns[1].form_submit_button(
            "×10",
            key=f"{widget_key}_mul10",
            on_click=_scale_widget_value,
            args=(widget_key, 10.0),
            width="stretch",
            help="Multiply this value by 10.",
        )


with st.container(key="calculate-anchor"):
    calculate_header = st.form_submit_button(
        "Calculate",
        icon=":material/calculate:",
        key="calculate_header",
        type="primary",
        width="stretch",
    )

_main_dg._form_data = None
with st.container(key="copy-link-anchor"):
    st.button(
        label=None,
        key="copy_scenario_link",
        icon=":material/link:",
        help="Copy a link to the calculated scenario",
    )
_bind_copy_link_button("copy_scenario_link")
_main_dg._form_data = _FormData("scenario_form")

# The tabs are a Streamlit radio: the frontend manages its checked state
# instantly (no rerun \u2014 it is a form widget), the radio group enforces
# exclusivity natively, and the CSS reads the checked input's value to open
# the matching panel. "\u2715" is the closed state. (Raw HTML radios/details do
# not work: Streamlit's page scripts suppress native form-control activation.)
st.session_state.setdefault("cost_model_radio", "\u2715")

with st.container(key="cost-rail"):
    st.radio(
        "Cost model",
        options=["\u2715", "DNA", "Amazon", "Azure", "Tape", "Custom"],
        key="cost_model_radio",
        horizontal=False,
        label_visibility="collapsed",
    )

with st.container(key="cost-panel"):
    st.markdown(
        '<div class="cost-panel-header">'
        '<div class="sidebar-kicker">Cost assumptions</div>'
        '<div class="cost-panel-title cost-title-none">Cost assumptions</div>'
        '<div class="cost-panel-title cost-title-dna">DNA cost assumptions</div>'
        '<div class="cost-panel-title cost-title-amazon">Amazon Deep Archive assumptions</div>'
        '<div class="cost-panel-title cost-title-azure">Azure Blob Archive assumptions</div>'
        '<div class="cost-panel-title cost-title-tape">Tape on-premise assumptions</div>'
        '<div class="cost-panel-title cost-title-custom">Custom storage assumptions</div>'
        '<div class="cost-panel-summary cost-summary-dna"><span>Synthesis</span><span>Sequencing</span><span>Durability</span></div>'
        '<div class="cost-panel-summary cost-summary-amazon"><span>Requests</span><span>Retrieval</span><span>Storage</span></div>'
        '<div class="cost-panel-summary cost-summary-azure"><span>Operations</span><span>Retrieval</span><span>Archive tier</span></div>'
        '<div class="cost-panel-summary cost-summary-tape"><span>Media</span><span>Hardware</span><span>Energy</span></div>'
        '<div class="cost-panel-summary cost-summary-custom"><span>Write</span><span>Store</span><span>Retrieve</span></div>'
        '</div>',
        unsafe_allow_html=True,
    )

    with st.container(key="cost_model_dna"):
        _render_model_reset("dna")
        with st.container(key="advanced-dna_cost_base_year"):
            dna_cost_base_year = st.number_input(
                "DNA cost base year", min_value=2000, max_value=2500,
                key="dna_cost_base_year",
                help="Year to which the editable synthesis and sequencing unit costs apply.",
            )
        dna_synthesis_cost = st.number_input(
            "Synthesis cost (USD/MB)", min_value=0.0,
            format=NUMBER_INPUT_FORMATS["dna_synthesis_cost"], key="dna_synthesis_cost",
            help="Cost in the DNA cost base year to synthesize enough bases for 1 MB of logical data, before redundancy and indexing overhead. "
            "Shown in significant-figure notation (e.g. 1e-07) so very small values stay visible instead of displaying as 0.",
        )
        _render_order_of_magnitude_steppers("dna_synthesis_cost")
        dna_sequencing_cost = st.number_input(
            "Sequencing cost (USD/MB)", min_value=0.0,
            format=NUMBER_INPUT_FORMATS["dna_sequencing_cost"], key="dna_sequencing_cost",
            help="Cost in the DNA cost base year to sequence 1 MB of retrieved logical data. "
            "Shown in significant-figure notation (e.g. 1e-07) so very small values stay visible instead of displaying as 0.",
        )
        _render_order_of_magnitude_steppers("dna_sequencing_cost")
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
        with st.container(key="advanced-dna_durability"):
            dna_durability = st.number_input(
                "DNA durability (years)", min_value=1, max_value=10_000,
                key="dna_durability",
                help="Years before the archive must be synthesized again. No replacement occurs at the exact end of the horizon.",
            )

    with st.container(key="cost_model_amazon"):
        _render_model_reset("amazon")
        with st.container(key="advanced-amazon_reference"):
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
            format=NUMBER_INPUT_FORMATS["amazon_put_per_1000"],
            key="amazon_put_per_1000",
            help="Charge for 1,000 requests when the archive is initially written.",
        )
        amazon_restore_per_1000 = st.number_input(
            "Bulk restore requests (USD/1,000)", min_value=0.0,
            format=NUMBER_INPUT_FORMATS["amazon_restore_per_1000"],
            key="amazon_restore_per_1000",
            help="Charge for 1,000 bulk restore-job requests. Asset size determines the request count.",
        )
        amazon_retrieval_per_tb = st.number_input(
            "Bulk data retrieval (USD/TB)", min_value=0.0,
            format=NUMBER_INPUT_FORMATS["amazon_retrieval_per_tb"],
            key="amazon_retrieval_per_tb",
            help="Capacity charge for retrieving one TB of archived data.",
        )
        amazon_storage_per_tb_month = st.number_input(
            "Storage (USD/TB/month)", min_value=0.0,
            format=NUMBER_INPUT_FORMATS["amazon_storage_per_tb_month"],
            key="amazon_storage_per_tb_month",
            help="Recurring monthly charge to retain one TB in Deep Archive.",
        )

    with st.container(key="cost_model_azure"):
        _render_model_reset("azure")
        with st.container(key="advanced-azure_reference"):
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
            format=NUMBER_INPUT_FORMATS["azure_write_per_1000"],
            key="azure_write_per_1000",
            help="Charge for 1,000 requests when the archive is initially written.",
        )
        azure_read_per_1000 = st.number_input(
            "Read requests (USD/1,000)", min_value=0.0,
            format=NUMBER_INPUT_FORMATS["azure_read_per_1000"],
            key="azure_read_per_1000",
            help="Charge for 1,000 retrieval requests. Asset size determines the request count.",
        )
        azure_retrieval_per_tb = st.number_input(
            "Data retrieval (USD/TB)", min_value=0.0,
            format=NUMBER_INPUT_FORMATS["azure_retrieval_per_tb"],
            key="azure_retrieval_per_tb",
            help="Capacity charge for retrieving one TB from the Archive tier.",
        )
        azure_storage_per_tb_month = st.number_input(
            "Storage (USD/TB/month)", min_value=0.0,
            format=NUMBER_INPUT_FORMATS["azure_storage_per_tb_month"],
            key="azure_storage_per_tb_month",
            help="Recurring monthly charge to retain one TB in the Archive tier.",
        )

    with st.container(key="cost_model_tape"):
        _render_model_reset("tape")
        with st.container(key="advanced-tape_reference"):
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
        st.caption(
            "Tape cartridge and hardware values are added together. If your hardware estimate already "
            "includes cartridges/media, set tape cartridges to 0 to avoid double counting."
        )
        tape_media_per_tb = st.number_input(
            "Tape cartridges (USD/TB per write)", min_value=0.0,
            format=NUMBER_INPUT_FORMATS["tape_media_per_tb"],
            key="tape_media_per_tb",
            help=(
                "Physical tape cartridge/media cost for one TB. The model charges it on the initial "
                "archive write and on each replacement write. Set this to 0 if the hardware value "
                "already includes cartridges."
            ),
        )
        tape_hardware_per_tb = st.number_input(
            "Tape library/drives (USD/TB amortized)", min_value=0.0,
            format=NUMBER_INPUT_FORMATS["tape_hardware_per_tb"],
            key="tape_hardware_per_tb",
            help=(
                "Tape library, drive, robotics, and supporting hardware cost allocated per TB. The "
                "model spreads this over the selected tape durability period as annual maintenance."
            ),
        )
        tape_energy_per_tb_year = st.number_input(
            "Energy (USD/TB/year)", min_value=0.0,
            format=NUMBER_INPUT_FORMATS["tape_energy_per_tb_year"],
            key="tape_energy_per_tb_year",
            help="Annual energy cost to retain one TB in the tape system.",
        )
        st.caption("Annual price declines")
        tape_media_decline = st.number_input(
            "Tape cartridge decline (%)", min_value=0.0, max_value=99.99,
            key="tape_media_decline",
            help="Annual reduction applied to tape cartridge/media purchase prices.",
        )
        tape_hardware_decline = st.number_input(
            "Tape library/drives decline (%)", min_value=0.0, max_value=99.99,
            key="tape_hardware_decline",
            help="Annual reduction applied to amortized tape library, drive, and robotics costs.",
        )
        tape_energy_decline = st.number_input(
            "Tape energy decline (%)", min_value=0.0, max_value=99.99,
            key="tape_energy_decline",
            help="Annual reduction applied to tape energy costs.",
        )

    with st.container(key="cost_model_custom"):
        _render_model_reset("custom")
        custom_name = st.text_input(
            "Display name", key="custom_name",
            help="Name used for the custom technology in charts, tables, and downloads.",
        )
        with st.container(key="advanced-custom_base_year"):
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
        with st.container(key="advanced-custom_replacement"):
            custom_replacement = st.number_input(
                "Replacement interval (years)", min_value=0, max_value=10_000,
                key="custom_replacement",
                help="Years between complete rewrites. Use 0 for a service with no replacement writes.",
            )

    # Phones: the panel carries its own Calculate bar (sticky at the sheet's
    # bottom) so a cost change commits without closing the panel. Hidden on
    # desktop, where the header button covers it.
    calculate_panel = st.form_submit_button(
        "Calculate",
        key="calculate_panel",
        type="primary",
        width="stretch",
    )

_main_dg._form_data = None

submitted = calculate_header or calculate_scenario or calculate_panel

if submitted:
    try:
        submitted_widgets = _snapshot_widgets()
        candidate = _scenario_from_widgets(submitted_widgets)
        _cached_simulation(candidate)
    except (KeyError, ValueError) as error:
        st.error(f"{error} Your last calculated results are shown below.")
    else:
        st.session_state["committed_widgets"] = submitted_widgets
        params = candidate.to_query_params()
        params.update(
            {
                "projection_end": str(int(submitted_widgets["projection_end"])),
                "log_scale": str(bool(submitted_widgets["log_scale"])),
                "theme": theme,
            }
        )
        st.query_params.from_dict(params)
        st.session_state["close_sidebar_on_mobile"] = True

committed_widgets = st.session_state["committed_widgets"]
_bind_form_controls(committed_widgets)

# Close the mobile input sheet after a successful calculation.
if st.session_state.pop("close_sidebar_on_mobile", False):
    st.iframe(_SIDEBAR_CLOSE_JS, width=1, height=1)

# All graphs render from the last calculated inputs; widget edits live inside
# the form and do not rerun the script, so the graphs only change when a
# Calculate submit commits a new snapshot.
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
dna_curve_end = max(
    scenario.dna_cost_base_year,
    min(2500, scenario.start_year + scenario.horizon_years - 1),
)
use_present_value = scenario.discount_rate_percent > 0
value_column = "present_value_usd" if use_present_value else "total_cost_usd"
chart_palette = palette_for(theme)

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
        <span class="model-item"><a href="https://doi.org/10.48550/arXiv.2608.26342" target="_blank" rel="noopener">Paper / About</a></span>
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
st.markdown(
    f"""
    <div class="scenario-bar">
        <span class="scenario-label">Active scenario</span>
        <span><strong>{scenario.start_year}-{period_end}</strong> archive period</span>
        <span><strong>{len(scenario.technologies)}</strong> technologies</span>
        <span>{cost_basis}</span>
        <span>{scale_label}</span>
        <span class="pending-notice" role="status" aria-live="polite"></span>
    </div>
    """,
    unsafe_allow_html=True,
)

metric_columns = st.columns(4)
metric_columns[0].metric("Archive", f"{_archive_size_label(scenario.archive_size_tb)} TB")
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
            "DNA lifecycle cost", _money(dna_total), f"{format_display_number(ratio)}x lowest", delta_color="off"
        )
    else:
        metric_columns[3].metric("DNA lifecycle cost", _money(dna_total))

st.markdown('<div class="workspace-kicker">Analysis workspace</div>', unsafe_allow_html=True)
def _remember_chart_option(key: str) -> None:
    st.session_state[f"saved_{key}"] = st.session_state[key]


@st.fragment
def _render_analysis() -> None:
    overview_tab, outlook_tab, dna_cost_tab, sensitivity_tab, assumptions_tab, about_tab = st.tabs(
        ["Lifecycle", "Start-year outlook", "DNA unit costs", "Sensitivity", "Assumptions", "About"],
        key="analysis_tabs", on_change="rerun",
    )

    if overview_tab.open:
        with overview_tab:
            _section_intro(
                "Lifecycle comparison",
                "Cumulative lifecycle cost for one archive opened in the selected start year. "
                "Each line includes initial and replacement writes, storage or operation, and expected retrieval.",
            )
            show_uncertainty = (
                "DNA" in scenario.technologies
                and st.checkbox(
                    "Show DNA uncertainty band",
                    key="show_uncertainty_band",
                    value=st.session_state.get("saved_show_uncertainty_band", False),
                    on_change=_remember_chart_option,
                    args=("show_uncertainty_band",),
                    help="Shades the P10-P90 range from sampling the synthesis and sequencing "
                    "decline rates +/-30% around your chosen values, everything else held fixed.",
                )
            )
            dna_uncertainty = (
                _cached_uncertainty_band(scenario, use_present_value) if show_uncertainty else None
            )
            lifecycle_crossovers = find_lifecycle_crossovers(result, use_present_value)
            lifecycle_figure = lifecycle_chart(
                result,
                use_present_value,
                log_scale,
                theme=theme,
                dna_uncertainty=dna_uncertainty,
                crossovers=lifecycle_crossovers,
            )
            _render_chart(
                lifecycle_figure,
                key="chart_lifecycle",
                filename="dna-storage-lifecycle",
            )
            _chart_downloads(
                "lifecycle",
                result.yearly,
                "dna-storage-lifecycle",
                "chart_lifecycle",
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
            breakdown_figure = breakdown_chart(result, log_scale, theme=theme)
            _render_chart(
                breakdown_figure,
                key="chart_breakdown",
                filename="dna-storage-cost-components",
            )
            breakdown_columns = ["technology", *[component for component in (
                "write_cost_usd", "read_cost_usd", "maintenance_cost_usd", "total_cost_usd"
            )]]
            _chart_downloads(
                "breakdown",
                result.totals[breakdown_columns],
                "dna-storage-cost-components",
                "chart_breakdown",
            )

    if outlook_tab.open:
        with outlook_tab:
            projection = _cached_projection(scenario, projection_end)
            _section_intro(
                "Start-year sensitivity",
                "The same archive workload and retention horizon are recalculated for every possible storage "
                "start year. A crossover is the first start year for which DNA's lifecycle cost is no greater "
                "than the comparison technology.",
            )
            projection_figure = projection_chart(projection, use_present_value, log_scale, theme=theme)
            _render_chart(
                projection_figure,
                key="chart_projection",
                filename="dna-storage-start-year-outlook",
            )
            _chart_downloads(
                "projection",
                projection,
                "dna-storage-start-year-outlook",
                "chart_projection",
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
                _themed_table(
                    ["Comparison", "First start year"],
                    [[row["Comparison"], row["First start year"]] for row in crossover_rows],
                )
            else:
                st.info("Include DNA and at least one comparison technology to calculate crossover years.")

    if dna_cost_tab.open:
        with dna_cost_tab:
            dna_costs = _cached_dna_costs(scenario, dna_curve_end)
            _section_intro(
                "DNA unit economics",
                f"Both curves begin with the editable {scenario.dna_cost_base_year} unit costs and apply the "
                f"selected annual decline rates through {dna_curve_end}. They are unit-cost assumptions, "
                "not lifecycle totals.",
            )
            show_history = st.checkbox(
                "Show historical context",
                value=st.session_state.get("saved_show_dna_history", True),
                key="show_dna_history",
                on_change=_remember_chart_option,
                args=("show_dna_history",),
                help="Synthesis: the paper's fitted historical trend back to 2000. Sequencing: "
                "NHGRI's measured cost per Mb, 2001-2022.",
            )
            synthesis_history = (
                _cached_synthesis_history(scenario.dna_cost_base_year) if show_history else None
            )
            observed_sequencing = _cached_observed_sequencing_costs() if show_history else None

            chart_columns = st.columns(2)
            with chart_columns[0]:
                synthesis_title = "DNA synthesis cost trajectory"
                synthesis_figure = dna_unit_cost_chart(
                    dna_costs,
                    "synthesis_cost_usd_per_mb",
                    synthesis_title,
                    chart_palette["unit_cost_colors"]["synthesis"],
                    log_scale,
                    theme=theme,
                    history=synthesis_history,
                )
                _render_chart(
                    synthesis_figure,
                    key="chart_dna_synthesis",
                    filename="dna-synthesis-cost-trajectory",
                )
                st.caption(
                    "Modeled cost to synthesize enough DNA bases for 1 MB of logical data. "
                    "Redundancy, indexing, and coding overhead are not added here. The dotted "
                    "trend is the paper's historical fit, not a raw dataset: no public per-MB "
                    "synthesis price series exists the way NHGRI's sequencing table does."
                )
                synthesis_data = dna_costs[["year", "synthesis_cost_usd_per_mb"]]
                _chart_downloads(
                    "dna_synthesis",
                    synthesis_data,
                    "dna-synthesis-cost-trajectory",
                    "chart_dna_synthesis",
                )
            with chart_columns[1]:
                sequencing_title = "DNA sequencing cost trajectory"
                sequencing_figure = dna_unit_cost_chart(
                    dna_costs,
                    "sequencing_cost_usd_per_mb",
                    sequencing_title,
                    chart_palette["unit_cost_colors"]["sequencing"],
                    log_scale,
                    theme=theme,
                    observed=observed_sequencing,
                )
                _render_chart(
                    sequencing_figure,
                    key="chart_dna_sequencing",
                    filename="dna-sequencing-cost-trajectory",
                )
                st.caption(
                    "Modeled cost to sequence and retrieve 1 MB of logical data. "
                    "It is applied to the share of the archive retrieved each year. Observed "
                    "markers are NHGRI's reported nominal-USD cost per Mb at each date, not "
                    "adjusted to this model's constant-USD convention."
                )
                sequencing_data = dna_costs[["year", "sequencing_cost_usd_per_mb"]]
                _chart_downloads(
                    "dna_sequencing",
                    sequencing_data,
                    "dna-sequencing-cost-trajectory",
                    "chart_dna_sequencing",
                )

    if sensitivity_tab.open:
        with sensitivity_tab:
            _section_intro(
                "Synthesis price needed to break even",
                f"Maximum synthesis price in {scenario.dna_cost_base_year} at which DNA matches each "
                "alternative over the archive lifetime. All other assumptions stay fixed.",
            )
            if "DNA" not in scenario.technologies:
                st.info("Select DNA to see its break-even price and cost sensitivity.")
            else:
                breakeven_frame = _cached_breakeven(scenario, use_present_value)
                if breakeven_frame.empty:
                    st.info("Select at least one comparison technology to calculate a break-even price.")
                else:
                    breakeven_figure = breakeven_chart(
                        breakeven_frame, scenario.dna_synthesis_cost_per_mb, theme=theme
                    )
                    _render_chart(
                        breakeven_figure,
                        key="chart_breakeven",
                        filename="dna-synthesis-breakeven",
                    )
                    st.caption(
                        "Not reachable means DNA retrieval alone exceeds the alternative's total cost, "
                        "even with free synthesis. The preservation archive example uses 0.001% annual retrieval."
                    )
                    def _breakeven_cell(breakeven: float) -> str:
                        if pd.isna(breakeven):
                            return "Not reachable"
                        return f"${format_display_number(breakeven)}/MB"

                    def _direction_cell(breakeven: float, reduction_factor: float) -> str:
                        if pd.isna(breakeven):
                            return "—"
                        if pd.isna(reduction_factor):
                            return "Needs free synthesis"
                        if breakeven < scenario.dna_synthesis_cost_per_mb:
                            return f"{format_display_number(reduction_factor)}x lower needed"
                        if breakeven > scenario.dna_synthesis_cost_per_mb:
                            if scenario.dna_synthesis_cost_per_mb == 0:
                                return "Already cheaper (free synthesis)"
                            headroom = breakeven / scenario.dna_synthesis_cost_per_mb
                            return f"Already cheaper ({format_display_number(headroom)}x headroom)"
                        return "At parity today"

                    breakeven_rows = [
                        [
                            row["technology"],
                            _breakeven_cell(row["breakeven_synthesis_cost_usd_per_mb"]),
                            _direction_cell(row["breakeven_synthesis_cost_usd_per_mb"], row["reduction_factor"]),
                        ]
                        for _, row in breakeven_frame.iterrows()
                    ]
                    _themed_table(
                        ["Comparison", "Break-even synthesis cost", "Today vs. break-even"],
                        breakeven_rows,
                    )
                    _chart_downloads(
                        "breakeven",
                        breakeven_frame,
                        "dna-synthesis-breakeven",
                        "chart_breakeven",
                    )

                st.markdown(
                    """
                    <div class="chart-divider">
                        <h3>What drives DNA's own cost</h3>
                        <p>Each bar re-runs the model with one input moved 50% below and 50% above your
                    current value, within valid input limits. Other assumptions stay fixed.</p>
                    </div>
                    """,
                    unsafe_allow_html=True,
                )
                sensitivity_frame = _cached_sensitivity(scenario, use_present_value)
                sensitivity_figure = sensitivity_chart(sensitivity_frame, use_present_value, theme=theme)
                _render_chart(
                    sensitivity_figure,
                    key="chart_sensitivity",
                    filename="dna-cost-sensitivity",
                )
                st.caption(
                    "The dotted line marks DNA's cost at your actual inputs. A bar with no visible width "
                    "means that input has no effect on DNA's cost in this scenario — average asset size, "
                    "for example, never enters DNA's own cost formula. Hover a bar for the exact low/high "
                    "values and resulting costs."
                )
                _chart_downloads(
                    "sensitivity",
                    sensitivity_frame,
                    "dna-cost-sensitivity",
                    "chart_sensitivity",
                )

    if assumptions_tab.open:
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
                _themed_table(
                    list(contract.columns),
                    [list(row) for row in contract.itertuples(index=False)],
                )
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

    if about_tab.open:
        with about_tab:
            _section_intro(
                "About this explorer",
                "An interactive implementation of a DNA storage cost model for comparing long-run archival "
                "economics across DNA, cloud archive, tape, and user-defined storage systems.",
            )
            about_column, paper_column = st.columns([1.25, 1], gap="large")
            with about_column:
                st.subheader("Purpose")
                st.write(
                    "The explorer turns the model assumptions into adjustable controls, then recalculates "
                    "lifecycle cost, start-year sensitivity, and DNA unit-cost trajectories for the selected "
                    "archive workload."
                )
                st.subheader("Model boundary")
                st.write(
                    "The results are scenario estimates, not procurement quotes. They use the visible inputs "
                    "and listed assumptions, and exclude operational details such as labor, taxes, cloud egress, "
                    "retrieval latency, and migration execution risk."
                )
            with paper_column:
                st.subheader("Reference paper")
                st.markdown(
                    "[DNA Storage Cost Model](https://doi.org/10.48550/arXiv.2608.26342)"
                )
                st.write(
                    "Use the linked paper as the source reference for the model framing and baseline assumptions."
                )
                st.subheader("Contact")
                st.markdown("[alex@el-shaikh.com](mailto:alex@el-shaikh.com)")


_render_analysis()
