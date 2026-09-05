from __future__ import annotations

import pandas as pd
import plotly.graph_objects as go

from .simulation import COMPONENTS, SimulationResult


PALETTES = {
    "light": {
        "technology_colors": {
            "DNA": "#bd4b38",
            "Amazon Deep Archive": "#2067b0",
            "Azure Blob Archive": "#2e8f5f",
            "Tape On-premise": "#7558b0",
            "Custom storage": "#75852e",
        },
        "component_colors": {
            "write_cost_usd": "#bd4b38",
            "read_cost_usd": "#2067b0",
            "maintenance_cost_usd": "#d19b2a",
        },
        "unit_cost_colors": {"synthesis": "#bd4b38", "sequencing": "#2067b0"},
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
            "DNA": "#d06a4e",
            "Amazon Deep Archive": "#3a92c9",
            "Azure Blob Archive": "#42ab72",
            "Tape On-premise": "#8f78d6",
            "Custom storage": "#8b9c38",
        },
        "component_colors": {
            "write_cost_usd": "#d06a4e",
            "read_cost_usd": "#3a92c9",
            "maintenance_cost_usd": "#ae821f",
        },
        "unit_cost_colors": {"synthesis": "#d06a4e", "sequencing": "#3a92c9"},
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


def _add_terminal_dot(
    figure: go.Figure,
    x: float,
    y: float,
    color: str,
    surface_color: str,
    name: str,
    legendgroup: str | None = None,
) -> None:
    """Small ring-marked dot on the line's last point — a quiet end-of-line
    accent, hidden from the legend and the unified hover."""
    figure.add_trace(
        go.Scatter(
            x=[x],
            y=[y],
            mode="markers",
            name=name,
            legendgroup=legendgroup,
            showlegend=False,
            hoverinfo="skip",
            marker={"color": color, "size": 8, "line": {"color": surface_color, "width": 2}},
        )
    )


def lifecycle_chart(
    result: SimulationResult, use_present_value: bool, log_scale: bool, theme: str = DEFAULT_THEME
) -> go.Figure:
    value = "cumulative_present_value_usd" if use_present_value else "cumulative_cost_usd"
    surface = palette_for(theme)["figure"]["plot_bgcolor"]
    figure = go.Figure()
    for technology in result.totals["technology"]:
        data = result.yearly[result.yearly["technology"] == technology]
        color = technology_color(technology, theme)
        figure.add_trace(
            go.Scatter(
                x=data["year"],
                y=data[value],
                name=technology,
                legendgroup=technology,
                mode="lines",
                line={"color": color, "width": 3},
                hovertemplate="%{x}<br>%{y:$,.3s}<extra>%{fullData.name}</extra>",
            )
        )
        last = data.iloc[-1]
        _add_terminal_dot(
            figure, last["year"], last[value], color, surface, technology, technology
        )
    figure.update_layout(
        title="Cumulative lifecycle cost",
        xaxis_title="Calendar year",
        yaxis_title="Present value (USD)" if use_present_value else "Cumulative cost (USD)",
        yaxis_type="log" if log_scale else "linear",
        legend={"orientation": "h", "y": 1.12, "x": 0, "groupclick": "togglegroup"},
        hovermode="x unified",
        margin={"l": 12, "r": 12, "t": 88, "b": 12},
        height=460,
    )
    return style_figure(figure, theme)


def breakdown_chart(result: SimulationResult, log_scale: bool, theme: str = DEFAULT_THEME) -> go.Figure:
    palette = palette_for(theme)
    component_colors = palette["component_colors"]
    surface = palette["figure"]["plot_bgcolor"]
    totals = result.totals
    figure = go.Figure()
    for component in COMPONENTS:
        figure.add_trace(
            go.Bar(
                x=totals["technology"],
                y=totals[component],
                name=COMPONENT_LABELS[component],
                marker_color=component_colors[component],
                # 2px surface-colored border: separates stacked segments and
                # adjacent bars the way a surface gap does.
                marker_line={"color": surface, "width": 2},
                hovertemplate="%{x}<br>%{y:$,.3s}<extra>%{fullData.name}</extra>",
            )
        )
    figure.update_layout(
        title="Undiscounted cost components",
        yaxis_title="Lifecycle cost (USD)",
        yaxis_type="log" if log_scale else "linear",
        barmode="stack",
        bargap=0.18,
        legend={"orientation": "h", "y": 1.12, "x": 0, "groupclick": "togglegroup"},
        hovermode="x unified",
        margin={"l": 12, "r": 12, "t": 88, "b": 12},
        height=430,
    )
    return style_figure(figure, theme)


def projection_chart(
    projection: pd.DataFrame, use_present_value: bool, log_scale: bool, theme: str = DEFAULT_THEME
) -> go.Figure:
    value = "present_value_usd" if use_present_value else "total_cost_usd"
    surface = palette_for(theme)["figure"]["plot_bgcolor"]
    figure = go.Figure()
    for technology in projection["technology"].drop_duplicates():
        data = projection[projection["technology"] == technology]
        color = technology_color(technology, theme)
        figure.add_trace(
            go.Scatter(
                x=data["start_year"],
                y=data[value],
                name=technology,
                legendgroup=technology,
                mode="lines",
                line={"color": color, "width": 3},
                hovertemplate="Start %{x}<br>%{y:$,.3s}<extra>%{fullData.name}</extra>",
            )
        )
        last = data.iloc[-1]
        _add_terminal_dot(
            figure, last["start_year"], last[value], color, surface, technology, technology
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
        # Wider hover pick-up radius, so thin lines are easy to hit on touch.
        hoverdistance=20,
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
        # Soft themed crosshair on x-unified hover.
        spikecolor=colors["axis_line"],
        spikethickness=1,
        spikedash="dot",
        spikemode="across",
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
    use_log = log_scale and bool((costs[value_column] > 0).any())
    surface = palette_for(theme)["figure"]["plot_bgcolor"]
    figure = go.Figure(
        go.Scatter(
            x=costs["year"],
            y=costs[value_column],
            mode="lines",
            line={"color": color, "width": 3},
            fill="tozeroy" if not use_log else None,
            fillcolor=color + "26" if not use_log else None,
            hovertemplate="%{x}<br>%{y:$,.4g} per MB<extra></extra>",
        )
    )
    last = costs.iloc[-1]
    _add_terminal_dot(figure, last["year"], last[value_column], color, surface, title)
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


def plotly_image_exports(
    figure: go.Figure,
    *,
    width: int = 1200,
    height: int | None = None,
    scale: int = 2,
) -> dict[str, bytes]:
    """Export the same Plotly figure object that is rendered in the app."""
    export_height = height or int(figure.layout.height or 500)
    return {
        file_format: figure.to_image(
            format=file_format,
            width=width,
            height=export_height,
            scale=scale,
        )
        for file_format in ("png", "svg")
    }


def lifecycle_exports(
    result: SimulationResult, use_present_value: bool, log_scale: bool, theme: str = DEFAULT_THEME
) -> dict[str, bytes]:
    return plotly_image_exports(lifecycle_chart(result, use_present_value, log_scale, theme))


def breakdown_exports(result: SimulationResult, log_scale: bool, theme: str = DEFAULT_THEME) -> dict[str, bytes]:
    return plotly_image_exports(breakdown_chart(result, log_scale, theme))


def projection_exports(
    projection: pd.DataFrame,
    use_present_value: bool,
    log_scale: bool,
    metadata: dict[str, str],
    theme: str = DEFAULT_THEME,
) -> dict[str, bytes]:
    return plotly_image_exports(projection_chart(projection, use_present_value, log_scale, theme))


def dna_unit_cost_exports(
    costs: pd.DataFrame,
    value_column: str,
    title: str,
    color: str,
    log_scale: bool,
    metadata: dict[str, str],
    theme: str = DEFAULT_THEME,
) -> dict[str, bytes]:
    return plotly_image_exports(
        dna_unit_cost_chart(costs, value_column, title, color, log_scale, theme)
    )


def lifecycle_export(
    result: SimulationResult,
    use_present_value: bool,
    file_format: str,
    log_scale: bool = True,
    theme: str = DEFAULT_THEME,
) -> bytes:
    """Backward-compatible single-format lifecycle export."""
    return lifecycle_exports(result, use_present_value, log_scale, theme)[file_format]
