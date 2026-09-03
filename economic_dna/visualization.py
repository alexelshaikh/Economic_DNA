from __future__ import annotations

from io import BytesIO

import matplotlib.pyplot as plt
import pandas as pd
import plotly.graph_objects as go
from matplotlib.ticker import StrMethodFormatter

from .simulation import COMPONENTS, SimulationResult


PALETTES = {
    "light": {
        "technology_colors": {
            "DNA": "#bd4b38",
            "Amazon Deep Archive": "#176b87",
            "Azure Blob Archive": "#2e7d57",
            "Tape On-premise": "#6d5b8c",
            "Custom storage": "#5c6b2f",
        },
        "component_colors": {
            "write_cost_usd": "#bd4b38",
            "read_cost_usd": "#176b87",
            "maintenance_cost_usd": "#d19b2a",
        },
        "unit_cost_colors": {"synthesis": "#bd4b38", "sequencing": "#176b87"},
        "figure": {
            "font_color": "#34413d",
            "title_color": "#18211f",
            "plot_bgcolor": "#ffffff",
            "legend_font_color": "#475650",
            "axis_line": "#cbd3cf",
            "tick_color": "#64706c",
            "axis_title_color": "#475650",
            "grid_color": "#e9edeb",
            "hover_bgcolor": "#18211f",
            "hover_font_color": "#ffffff",
        },
        "export": {
            "facecolor": "#ffffff",
            "grid_color": "#e2e4e6",
            "text_color": "#55595d",
            "title_color": "#18211f",
        },
    },
    "dark": {
        "technology_colors": {
            "DNA": "#e07a5f",
            "Amazon Deep Archive": "#4aa3c7",
            "Azure Blob Archive": "#52a878",
            "Tape On-premise": "#9d8bc2",
            "Custom storage": "#a8b45f",
        },
        "component_colors": {
            "write_cost_usd": "#e07a5f",
            "read_cost_usd": "#4aa3c7",
            "maintenance_cost_usd": "#e0b84f",
        },
        "unit_cost_colors": {"synthesis": "#e07a5f", "sequencing": "#4aa3c7"},
        "figure": {
            "font_color": "#b7c4c0",
            "title_color": "#e7edeb",
            "plot_bgcolor": "#171f20",
            "legend_font_color": "#b7c4c0",
            "axis_line": "#38464a",
            "tick_color": "#93a19d",
            "axis_title_color": "#b7c4c0",
            "grid_color": "#232e31",
            "hover_bgcolor": "#0b1112",
            "hover_font_color": "#e7edeb",
        },
        "export": {
            "facecolor": "#171f20",
            "grid_color": "#2b363a",
            "text_color": "#93a19d",
            "title_color": "#e7edeb",
        },
    },
}
DEFAULT_THEME = "light"
# Backward-compatible aliases for the light theme.
COLORS = PALETTES["light"]["technology_colors"]
COMPONENT_COLORS = PALETTES["light"]["component_colors"]
COMPONENT_LABELS = {
    "write_cost_usd": "Write & replacement",
    "read_cost_usd": "Retrieval",
    "maintenance_cost_usd": "Storage & operation",
}


def palette_for(theme: str) -> dict:
    return PALETTES.get(theme, PALETTES[DEFAULT_THEME])


def technology_color(technology: str, theme: str = DEFAULT_THEME) -> str:
    colors = palette_for(theme)["technology_colors"]
    return colors.get(technology, colors["Custom storage"])


def lifecycle_chart(
    result: SimulationResult, use_present_value: bool, log_scale: bool, theme: str = DEFAULT_THEME
) -> go.Figure:
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
                line={"color": technology_color(technology, theme), "width": 2.5},
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
        margin={"l": 12, "r": 12, "t": 88, "b": 12},
        height=460,
    )
    return style_figure(figure, theme)


def breakdown_chart(result: SimulationResult, log_scale: bool, theme: str = DEFAULT_THEME) -> go.Figure:
    component_colors = palette_for(theme)["component_colors"]
    totals = result.totals
    figure = go.Figure()
    for component in COMPONENTS:
        figure.add_trace(
            go.Bar(
                x=totals["technology"],
                y=totals[component],
                name=COMPONENT_LABELS[component],
                marker_color=component_colors[component],
                hovertemplate="%{x}<br>%{y:$,.3s}<extra>%{fullData.name}</extra>",
            )
        )
    figure.update_layout(
        title="Undiscounted cost components",
        yaxis_title="Lifecycle cost (USD)",
        yaxis_type="log" if log_scale else "linear",
        barmode="stack",
        legend={"orientation": "h", "y": 1.12, "x": 0},
        margin={"l": 12, "r": 12, "t": 88, "b": 12},
        height=430,
    )
    return style_figure(figure, theme)


def projection_chart(
    projection: pd.DataFrame, use_present_value: bool, log_scale: bool, theme: str = DEFAULT_THEME
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
                line={"color": technology_color(technology, theme), "width": 2.5},
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
        margin={"l": 12, "r": 12, "t": 88, "b": 12},
        height=500,
    )
    return style_figure(figure, theme)


def style_figure(figure: go.Figure, theme: str = DEFAULT_THEME) -> go.Figure:
    palette = palette_for(theme)
    colors = palette["figure"]
    technology_colors = palette["technology_colors"]
    figure.update_layout(
        template="plotly_white",
        font={
            "family": "Aptos, Segoe UI, Arial, sans-serif",
            "size": 13,
            "color": colors["font_color"],
        },
        title={"x": 0.5, "xanchor": "center"},
        title_font={"size": 21, "color": colors["title_color"]},
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor=colors["plot_bgcolor"],
        colorway=list(technology_colors.values()),
        hoverlabel={
            "bgcolor": colors["hover_bgcolor"],
            "bordercolor": colors["hover_bgcolor"],
            "font": {"color": colors["hover_font_color"], "family": "Aptos, Segoe UI, Arial, sans-serif"},
        },
        legend={"font": {"size": 14, "color": colors["legend_font_color"]}},
    )
    figure.update_xaxes(
        showgrid=False,
        linecolor=colors["axis_line"],
        tickcolor=colors["axis_line"],
        ticks="outside",
        tickfont={"color": colors["tick_color"]},
        title_font={"color": colors["axis_title_color"], "size": 15},
        automargin=True,
        zeroline=False,
    )
    figure.update_yaxes(
        gridcolor=colors["grid_color"],
        linecolor=colors["axis_line"],
        tickcolor=colors["axis_line"],
        tickfont={"color": colors["tick_color"]},
        title_font={"color": colors["axis_title_color"], "size": 15},
        automargin=True,
        zeroline=False,
    )
    return figure


def dna_unit_cost_chart(
    costs: pd.DataFrame,
    value_column: str,
    title: str,
    color: str,
    log_scale: bool,
    theme: str = DEFAULT_THEME,
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
        margin={"l": 12, "r": 12, "t": 70, "b": 12},
        height=390,
        showlegend=False,
    )
    if use_log:
        figure.update_yaxes(tickformat=".1e")
    return style_figure(figure, theme)


def _figure_exports(fig: plt.Figure, palette: dict) -> dict[str, bytes]:
    exports: dict[str, bytes] = {}
    for file_format in ("png", "svg"):
        output = BytesIO()
        fig.savefig(
            output,
            format=file_format,
            dpi=220,
            bbox_inches="tight",
            facecolor=palette["export"]["facecolor"],
        )
        exports[file_format] = output.getvalue()
    plt.close(fig)
    return exports


def _finish_export_figure(
    fig: plt.Figure, ax: plt.Axes, metadata: dict[str, str], palette: dict
) -> None:
    colors = palette["export"]
    ax.grid(axis="y", color=colors["grid_color"], linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(colors["grid_color"])
    ax.tick_params(colors=colors["text_color"])
    ax.xaxis.label.set_color(colors["text_color"])
    ax.yaxis.label.set_color(colors["text_color"])
    ax.title.set_color(colors["title_color"])
    legend = ax.get_legend()
    if legend:
        for text in legend.get_texts():
            text.set_color(colors["text_color"])
    if ax.get_yscale() == "log":
        # Avoid MathText tick labels, which fail in some Streamlit/Matplotlib renderers.
        ax.yaxis.set_major_formatter(StrMethodFormatter("{x:.0e}"))
    fig.text(
        0.01,
        0.01,
        f"Economic DNA model v{metadata['model_version']} | {metadata['currency']} | "
        f"{metadata['disclaimer']}",
        fontsize=8,
        color=colors["text_color"],
    )
    fig.patch.set_facecolor(colors["facecolor"])
    ax.set_facecolor(colors["facecolor"])
    fig.tight_layout(rect=(0, 0.04, 1, 1))


def lifecycle_exports(
    result: SimulationResult, use_present_value: bool, log_scale: bool, theme: str = DEFAULT_THEME
) -> dict[str, bytes]:
    palette = palette_for(theme)
    value = "cumulative_present_value_usd" if use_present_value else "cumulative_cost_usd"
    fig, ax = plt.subplots(figsize=(10, 5.6))
    for technology in result.totals["technology"]:
        data = result.yearly[result.yearly["technology"] == technology]
        ax.plot(
            data["year"],
            data[value],
            label=technology,
            color=technology_color(technology, theme),
            linewidth=2.2,
        )
    if log_scale:
        ax.set_yscale("log")
    ax.set_title("Cumulative lifecycle cost", loc="left", fontweight="bold")
    ax.set_xlabel("Calendar year")
    ax.set_ylabel("Present value (USD)" if use_present_value else "Cumulative cost (USD)")
    ax.legend(frameon=False, ncol=2)
    _finish_export_figure(fig, ax, result.metadata, palette)
    return _figure_exports(fig, palette)


def breakdown_exports(result: SimulationResult, log_scale: bool, theme: str = DEFAULT_THEME) -> dict[str, bytes]:
    palette = palette_for(theme)
    component_colors = palette["component_colors"]
    totals = result.totals
    fig, ax = plt.subplots(figsize=(10, 5.4))
    bottom = pd.Series(0.0, index=totals.index)
    for component in COMPONENTS:
        ax.bar(
            totals["technology"],
            totals[component],
            bottom=bottom,
            label=COMPONENT_LABELS[component],
            color=component_colors[component],
        )
        bottom += totals[component]
    if log_scale:
        ax.set_yscale("log")
    ax.set_title("Undiscounted cost components", loc="left", fontweight="bold")
    ax.set_ylabel("Lifecycle cost (USD)")
    ax.tick_params(axis="x", rotation=12)
    ax.legend(frameon=False, ncol=3)
    _finish_export_figure(fig, ax, result.metadata, palette)
    return _figure_exports(fig, palette)


def projection_exports(
    projection: pd.DataFrame,
    use_present_value: bool,
    log_scale: bool,
    metadata: dict[str, str],
    theme: str = DEFAULT_THEME,
) -> dict[str, bytes]:
    palette = palette_for(theme)
    value = "present_value_usd" if use_present_value else "total_cost_usd"
    fig, ax = plt.subplots(figsize=(10, 5.8))
    for technology in projection["technology"].drop_duplicates():
        data = projection[projection["technology"] == technology]
        ax.plot(
            data["start_year"],
            data[value],
            label=technology,
            color=technology_color(technology, theme),
            linewidth=2.2,
        )
    if log_scale:
        ax.set_yscale("log")
    ax.set_title("Lifecycle cost by storage start year", loc="left", fontweight="bold")
    ax.set_xlabel("Storage start year")
    ax.set_ylabel("Present value (USD)" if use_present_value else "Lifecycle cost (USD)")
    ax.legend(frameon=False, ncol=2)
    _finish_export_figure(fig, ax, metadata, palette)
    return _figure_exports(fig, palette)


def dna_unit_cost_exports(
    costs: pd.DataFrame,
    value_column: str,
    title: str,
    color: str,
    log_scale: bool,
    metadata: dict[str, str],
    theme: str = DEFAULT_THEME,
) -> dict[str, bytes]:
    palette = palette_for(theme)
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    ax.plot(costs["year"], costs[value_column], color=color, linewidth=2.6)
    if log_scale and bool((costs[value_column] > 0).any()):
        ax.set_yscale("log")
    else:
        ax.fill_between(costs["year"], costs[value_column], color=color, alpha=0.12)
    ax.set_title(title, loc="left", fontweight="bold")
    ax.set_xlabel("Calendar year")
    ax.set_ylabel("USD per MB")
    _finish_export_figure(fig, ax, metadata, palette)
    return _figure_exports(fig, palette)


def lifecycle_export(
    result: SimulationResult,
    use_present_value: bool,
    file_format: str,
    log_scale: bool = True,
    theme: str = DEFAULT_THEME,
) -> bytes:
    """Backward-compatible single-format lifecycle export."""
    return lifecycle_exports(result, use_present_value, log_scale, theme)[file_format]
