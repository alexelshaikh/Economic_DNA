from __future__ import annotations

from io import BytesIO

import matplotlib.pyplot as plt
import pandas as pd
import plotly.graph_objects as go

from .simulation import COMPONENTS, SimulationResult


COLORS = {
    "DNA": "#b44737",
    "Amazon Deep Archive": "#176b87",
    "Azure Blob Archive": "#2e7d57",
    "Tape On-premise": "#6d5b8c",
    "Custom storage": "#5c6b2f",
}
COMPONENT_LABELS = {
    "write_cost_usd": "Write & replacement",
    "read_cost_usd": "Retrieval",
    "maintenance_cost_usd": "Storage & operation",
}


def technology_color(technology: str) -> str:
    return COLORS.get(technology, COLORS["Custom storage"])
COMPONENT_COLORS = {
    "write_cost_usd": "#b44737",
    "read_cost_usd": "#176b87",
    "maintenance_cost_usd": "#d19b2a",
}


def lifecycle_chart(result: SimulationResult, use_present_value: bool, log_scale: bool) -> go.Figure:
    value = "cumulative_present_value_usd" if use_present_value else "cumulative_cost_usd"
    figure = go.Figure()
    for technology in result.totals["technology"]:
        data = result.yearly[result.yearly["technology"] == technology]
        figure.add_trace(
            go.Scatter(
                x=data["year"],
                y=data[value],
                name=technology,
                mode="lines",
                line={"color": technology_color(technology), "width": 2.5},
                hovertemplate="%{x}<br>%{y:$,.3s}<extra>%{fullData.name}</extra>",
            )
        )
    figure.update_layout(
        title="Cumulative lifecycle cost",
        xaxis_title="Calendar year",
        yaxis_title="Present value (USD)" if use_present_value else "Cumulative cost (USD)",
        yaxis_type="log" if log_scale else "linear",
        legend={"orientation": "h", "y": 1.12, "x": 0},
        hovermode="x unified",
        margin={"l": 12, "r": 12, "t": 80, "b": 12},
        height=460,
    )
    return style_figure(figure)


def breakdown_chart(result: SimulationResult, log_scale: bool) -> go.Figure:
    totals = result.totals
    figure = go.Figure()
    for component in COMPONENTS:
        figure.add_trace(
            go.Bar(
                x=totals["technology"],
                y=totals[component],
                name=COMPONENT_LABELS[component],
                marker_color=COMPONENT_COLORS[component],
                hovertemplate="%{x}<br>%{y:$,.3s}<extra>%{fullData.name}</extra>",
            )
        )
    figure.update_layout(
        title="Undiscounted cost components",
        yaxis_title="Lifecycle cost (USD)",
        yaxis_type="log" if log_scale else "linear",
        barmode="stack",
        legend={"orientation": "h", "y": 1.12, "x": 0},
        margin={"l": 12, "r": 12, "t": 80, "b": 12},
        height=430,
    )
    return style_figure(figure)


def projection_chart(
    projection: pd.DataFrame, use_present_value: bool, log_scale: bool
) -> go.Figure:
    value = "present_value_usd" if use_present_value else "total_cost_usd"
    figure = go.Figure()
    for technology in projection["technology"].drop_duplicates():
        data = projection[projection["technology"] == technology]
        figure.add_trace(
            go.Scatter(
                x=data["start_year"],
                y=data[value],
                name=technology,
                mode="lines",
                line={"color": technology_color(technology), "width": 2.5},
                hovertemplate="Start %{x}<br>%{y:$,.3s}<extra>%{fullData.name}</extra>",
            )
        )
    figure.update_layout(
        title="Lifecycle cost by storage start year",
        xaxis_title="Storage start year",
        yaxis_title="Present value (USD)" if use_present_value else "Lifecycle cost (USD)",
        yaxis_type="log" if log_scale else "linear",
        legend={"orientation": "h", "y": 1.12, "x": 0},
        hovermode="x unified",
        margin={"l": 12, "r": 12, "t": 80, "b": 12},
        height=500,
    )
    return style_figure(figure)


def style_figure(figure: go.Figure) -> go.Figure:
    figure.update_layout(
        template="plotly_white",
        font={"family": "Inter, Arial, sans-serif", "color": "#202326"},
        title_font={"size": 19},
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="#ffffff",
        colorway=list(COLORS.values()),
    )
    figure.update_xaxes(showgrid=False, linecolor="#d9dcdf")
    figure.update_yaxes(gridcolor="#e8eaec", linecolor="#d9dcdf")
    return figure


def dna_unit_cost_chart(
    costs: pd.DataFrame,
    value_column: str,
    title: str,
    color: str,
    log_scale: bool,
) -> go.Figure:
    figure = go.Figure(
        go.Scatter(
            x=costs["year"],
            y=costs[value_column],
            mode="lines",
            line={"color": color, "width": 2.8},
            fill="tozeroy" if not log_scale else None,
            hovertemplate="%{x}<br>%{y:$,.4g} per MB<extra></extra>",
        )
    )
    use_log = log_scale and bool((costs[value_column] > 0).any())
    figure.update_layout(
        title=title,
        xaxis_title="Calendar year",
        yaxis_title="USD per MB",
        yaxis_type="log" if use_log else "linear",
        margin={"l": 12, "r": 12, "t": 62, "b": 12},
        height=390,
        showlegend=False,
    )
    if use_log:
        figure.update_yaxes(tickformat=".1e")
    return style_figure(figure)


def _figure_exports(fig: plt.Figure) -> dict[str, bytes]:
    exports: dict[str, bytes] = {}
    for file_format in ("png", "svg"):
        output = BytesIO()
        fig.savefig(
            output,
            format=file_format,
            dpi=220,
            bbox_inches="tight",
            facecolor="white",
        )
        exports[file_format] = output.getvalue()
    plt.close(fig)
    return exports


def _finish_export_figure(fig: plt.Figure, ax: plt.Axes, metadata: dict[str, str]) -> None:
    ax.grid(axis="y", color="#e2e4e6", linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.text(
        0.01,
        0.01,
        f"Economic DNA model v{metadata['model_version']} | {metadata['currency']} | "
        f"{metadata['disclaimer']}",
        fontsize=8,
        color="#55595d",
    )
    fig.tight_layout(rect=(0, 0.04, 1, 1))


def lifecycle_exports(
    result: SimulationResult, use_present_value: bool, log_scale: bool
) -> dict[str, bytes]:
    value = "cumulative_present_value_usd" if use_present_value else "cumulative_cost_usd"
    fig, ax = plt.subplots(figsize=(10, 5.6))
    for technology in result.totals["technology"]:
        data = result.yearly[result.yearly["technology"] == technology]
        ax.plot(
            data["year"],
            data[value],
            label=technology,
            color=technology_color(technology),
            linewidth=2.2,
        )
    if log_scale:
        ax.set_yscale("log")
    ax.set_title("Cumulative lifecycle cost", loc="left", fontweight="bold")
    ax.set_xlabel("Calendar year")
    ax.set_ylabel("Present value (USD)" if use_present_value else "Cumulative cost (USD)")
    ax.legend(frameon=False, ncol=2)
    _finish_export_figure(fig, ax, result.metadata)
    return _figure_exports(fig)


def breakdown_exports(result: SimulationResult, log_scale: bool) -> dict[str, bytes]:
    totals = result.totals
    fig, ax = plt.subplots(figsize=(10, 5.4))
    bottom = pd.Series(0.0, index=totals.index)
    for component in COMPONENTS:
        ax.bar(
            totals["technology"],
            totals[component],
            bottom=bottom,
            label=COMPONENT_LABELS[component],
            color=COMPONENT_COLORS[component],
        )
        bottom += totals[component]
    if log_scale:
        ax.set_yscale("log")
    ax.set_title("Undiscounted cost components", loc="left", fontweight="bold")
    ax.set_ylabel("Lifecycle cost (USD)")
    ax.tick_params(axis="x", rotation=12)
    ax.legend(frameon=False, ncol=3)
    _finish_export_figure(fig, ax, result.metadata)
    return _figure_exports(fig)


def projection_exports(
    projection: pd.DataFrame,
    use_present_value: bool,
    log_scale: bool,
    metadata: dict[str, str],
) -> dict[str, bytes]:
    value = "present_value_usd" if use_present_value else "total_cost_usd"
    fig, ax = plt.subplots(figsize=(10, 5.8))
    for technology in projection["technology"].drop_duplicates():
        data = projection[projection["technology"] == technology]
        ax.plot(
            data["start_year"],
            data[value],
            label=technology,
            color=technology_color(technology),
            linewidth=2.2,
        )
    if log_scale:
        ax.set_yscale("log")
    ax.set_title("Lifecycle cost by storage start year", loc="left", fontweight="bold")
    ax.set_xlabel("Storage start year")
    ax.set_ylabel("Present value (USD)" if use_present_value else "Lifecycle cost (USD)")
    ax.legend(frameon=False, ncol=2)
    _finish_export_figure(fig, ax, metadata)
    return _figure_exports(fig)


def dna_unit_cost_exports(
    costs: pd.DataFrame,
    value_column: str,
    title: str,
    color: str,
    log_scale: bool,
    metadata: dict[str, str],
) -> dict[str, bytes]:
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    ax.plot(costs["year"], costs[value_column], color=color, linewidth=2.6)
    if log_scale and bool((costs[value_column] > 0).any()):
        ax.set_yscale("log")
    else:
        ax.fill_between(costs["year"], costs[value_column], color=color, alpha=0.12)
    ax.set_title(title, loc="left", fontweight="bold")
    ax.set_xlabel("Calendar year")
    ax.set_ylabel("USD per MB")
    _finish_export_figure(fig, ax, metadata)
    return _figure_exports(fig)


def lifecycle_export(
    result: SimulationResult,
    use_present_value: bool,
    file_format: str,
    log_scale: bool = True,
) -> bytes:
    """Backward-compatible single-format lifecycle export."""
    return lifecycle_exports(result, use_present_value, log_scale)[file_format]
