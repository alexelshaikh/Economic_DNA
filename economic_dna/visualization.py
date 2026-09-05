from __future__ import annotations

import math

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


_COMPACT_SUFFIXES = ((1e12, "T"), (1e9, "B"), (1e6, "M"), (1e3, "K"))
_LINEAR_STEP_MULTIPLIERS = (1.0, 2.0, 2.5, 5.0, 10.0)
# Outside this window plain decimals grow unreadably long, so ticks switch
# to power-of-ten notation instead.
_PLAIN_NOTATION_FLOOR = 1e-3
_PLAIN_NOTATION_CEILING = 1e15
_MAX_TICKS = 8


def _compact_number(value: float) -> str:
    """Axis label without a currency symbol — the axis title carries the unit."""
    absolute = abs(value)
    for threshold, suffix in _COMPACT_SUFFIXES:
        if absolute >= threshold:
            return f"{value / threshold:g}{suffix}"
    if value == 0 or absolute >= 1:
        return f"{value:g}"
    decimals = -math.floor(math.log10(absolute)) + 2
    return f"{value:.{decimals}f}".rstrip("0").rstrip(".")


def _power_number(value: float) -> str:
    """`2e-07` reads as 2x10^-7 — one short label instead of a long decimal."""
    if value == 0:
        return "0"
    power = math.floor(math.log10(abs(value)))
    mantissa = value / 10.0**power
    if abs(mantissa - 1.0) < 1e-9:
        return f"10<sup>{power}</sup>"
    return f"{mantissa:g}x10<sup>{power}</sup>"


_SUPERSCRIPT_TRANSLATION = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def _plain_power_number(value: float) -> str:
    """Same style as `_power_number` but with Unicode superscript characters
    instead of an HTML `<sup>` tag, for plain-text contexts (e.g. an
    HTML-escaped table cell) that would otherwise show the literal tag."""
    if value == 0:
        return "0"
    power = math.floor(math.log10(abs(value)))
    mantissa = value / 10.0**power
    exponent = str(power).translate(_SUPERSCRIPT_TRANSLATION)
    if abs(mantissa - 1.0) < 1e-9:
        return f"10{exponent}"
    return f"{mantissa:g}×10{exponent}"


def format_display_number(value: float) -> str:
    """Public one-value formatter for use outside a chart figure (e.g. a
    plain HTML table) that still avoids raw scientific notation ("1.69e+08")
    at extreme magnitudes, without the chart-only HTML tick/hover markup."""
    if value == 0 or _PLAIN_NOTATION_FLOOR <= abs(value) < _PLAIN_NOTATION_CEILING:
        return _compact_number(value)
    return _plain_power_number(value)


def _tick_labels(tickvals: list[float]) -> list[str]:
    """One notation for the whole axis, so labels never mix styles."""
    positive = [value for value in tickvals if value > 0]
    if not positive:
        return [_compact_number(value) for value in tickvals]
    if min(positive) < _PLAIN_NOTATION_FLOOR or max(positive) >= _PLAIN_NOTATION_CEILING:
        return [_power_number(value) for value in tickvals]
    return [_compact_number(value) for value in tickvals]


def _format_pair(low: float, high: float) -> tuple[str, str]:
    """Format two related values (e.g. a sensitivity input's low/high test
    points) with one shared notation, so a hover tooltip never shows one in
    plain decimal and the other in scientific notation."""
    labels = _tick_labels([low, high])
    return labels[0], labels[1]


def _format_money_pair(low: float, high: float) -> tuple[str, str]:
    labels = _format_pair(low, high)
    return f"${labels[0]}", f"${labels[1]}"


def _log_ticks(minimum: float, maximum: float) -> tuple[list[float], list[float]]:
    """Decade ticks over the full data span, thinned so a 10-decade axis stays
    readable. Ticks are laid out from the top down so the highest decade is
    always labelled."""
    start = math.floor(math.log10(minimum))
    end = math.ceil(math.log10(maximum))
    if end <= start:
        end = start + 1
    stride = max(1, math.ceil((end - start) / (_MAX_TICKS - 1)))
    powers = sorted(power for power in range(end, start - 1, -stride))
    return [10.0**power for power in powers], [float(start), float(end)]


def _linear_ticks(
    minimum: float, maximum: float, baseline_zero: bool
) -> tuple[list[float], list[float]]:
    """Ticks on a 1/2/2.5/5 x 10^k step, with the axis bounds snapped onto that
    step so the top and bottom edges are labelled."""
    # Anchor at zero unless the data sits in a narrow band well above it.
    lower = 0.0 if baseline_zero or minimum <= 0.25 * maximum else minimum
    span = maximum - lower
    if span <= 0:
        return [], []
    raw_step = span / (_MAX_TICKS - 1)
    magnitude = 10.0 ** math.floor(math.log10(raw_step))
    step = next(
        (
            multiplier * magnitude
            for multiplier in _LINEAR_STEP_MULTIPLIERS
            if raw_step <= multiplier * magnitude
        ),
        10.0 * magnitude,
    )
    first = math.floor(lower / step)
    last = math.ceil(maximum / step)
    tickvals = [step * index for index in range(first, last + 1)]
    return tickvals, [step * first, step * last]


def _trace_y_values(figure: go.Figure) -> list[float]:
    values: list[float] = []
    bar_totals: dict[object, float] = {}
    all_bars = bool(figure.data) and all(trace.type == "bar" for trace in figure.data)
    for trace in figure.data:
        y_values = getattr(trace, "y", None)
        if y_values is None:
            continue
        x_values = getattr(trace, "x", None)
        for index, raw_value in enumerate(y_values):
            try:
                value = float(raw_value)
            except (TypeError, ValueError):
                continue
            if not math.isfinite(value) or value <= 0:
                continue
            values.append(value)
            if all_bars and x_values is not None:
                key = x_values[index]
                bar_totals[key] = bar_totals.get(key, 0.0) + value
    if bar_totals:
        values.extend(total for total in bar_totals.values() if total > 0)
    return values


def _apply_cost_axis_format(figure: go.Figure) -> None:
    """Label the y axis with evenly spaced, unprefixed ticks and pin the axis
    range to them, so extreme inputs cannot produce Plotly's mixed
    decade/minor-tick labels or an unlabelled edge."""
    values = _trace_y_values(figure)
    if not values:
        return
    minimum, maximum = min(values), max(values)
    if figure.layout.yaxis.type == "log":
        tickvals, axis_range = _log_ticks(minimum, maximum)
    else:
        stacked_bars = bool(figure.data) and all(
            trace.type == "bar" for trace in figure.data
        )
        tickvals, axis_range = _linear_ticks(minimum, maximum, stacked_bars)
    if not tickvals:
        return
    figure.update_yaxes(
        tickmode="array",
        tickvals=tickvals,
        ticktext=_tick_labels(tickvals),
        range=axis_range,
    )


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
    result: SimulationResult,
    use_present_value: bool,
    log_scale: bool,
    theme: str = DEFAULT_THEME,
    *,
    dna_uncertainty: pd.DataFrame | None = None,
    crossovers: dict[str, int | None] | None = None,
) -> go.Figure:
    value = "cumulative_present_value_usd" if use_present_value else "cumulative_cost_usd"
    surface = palette_for(theme)["figure"]["plot_bgcolor"]
    figure = go.Figure()
    if dna_uncertainty is not None and not dna_uncertainty.empty:
        dna_color = technology_color("DNA", theme)
        band_years = list(dna_uncertainty["year"]) + list(dna_uncertainty["year"][::-1])
        band_values = list(dna_uncertainty["p90"]) + list(dna_uncertainty["p10"][::-1])
        figure.add_trace(
            go.Scatter(
                x=band_years,
                y=band_values,
                name="DNA P10-P90 (decline-rate uncertainty)",
                fill="toself",
                fillcolor=dna_color + "26",
                line={"color": "rgba(0,0,0,0)"},
                hoverinfo="skip",
                legendgroup="DNA",
            )
        )
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
    _apply_cost_axis_format(figure)
    if crossovers:
        # Staggered annotation heights (in paper coordinates, so they land
        # just under the legend regardless of the y axis's data range) keep
        # 2-3 simultaneous crossover labels from overlapping.
        for index, (technology, year) in enumerate(
            (technology, year) for technology, year in crossovers.items() if year is not None
        ):
            figure.add_vline(
                x=year,
                line={"color": technology_color(technology, theme), "width": 1.5, "dash": "dash"},
            )
            figure.add_annotation(
                x=year,
                y=1.0 - 0.08 * index,
                yref="paper",
                text=f"DNA ≤ {technology}: {year}",
                showarrow=False,
                font={"color": technology_color(technology, theme), "size": 11},
                bgcolor=surface,
                xanchor="left",
                xshift=6,
            )
    return style_figure(figure, theme)


# Pattern fills supplement (never replace) the component colors, so the
# stacked breakdown stays legible for colorblind viewers and in grayscale
# prints.
COMPONENT_PATTERNS = {
    "write_cost_usd": "",
    "read_cost_usd": "/",
    "maintenance_cost_usd": "x",
}


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
                marker_pattern={
                    "shape": COMPONENT_PATTERNS[component],
                    "fgcolor": surface,
                    "size": 6,
                    "solidity": 0.35,
                },
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
    _apply_cost_axis_format(figure)
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
    _apply_cost_axis_format(figure)
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
    *,
    history: pd.DataFrame | None = None,
    observed: pd.DataFrame | None = None,
) -> go.Figure:
    """`history` is the paper's fitted historical trend (no raw dataset backs
    it); `observed` is real measured data (currently only NHGRI sequencing
    costs). Either shows what the projection is anchored to instead of
    asking a viewer to take the decline rate on faith."""
    use_log = log_scale and bool((costs[value_column] > 0).any())
    surface = palette_for(theme)["figure"]["plot_bgcolor"]
    show_context = history is not None or observed is not None
    figure = go.Figure()
    if history is not None and not history.empty:
        figure.add_trace(
            go.Scatter(
                x=history["year"],
                y=history[value_column],
                name="Historical trend (paper fit)",
                mode="lines",
                line={"color": color, "width": 2, "dash": "dot"},
                opacity=0.6,
                hovertemplate="%{x:.0f}<br>%{y:$,.4g} per MB<extra>Historical trend</extra>",
            )
        )
    if observed is not None and not observed.empty:
        figure.add_trace(
            go.Scatter(
                x=observed["year"],
                y=observed[value_column],
                name="Observed (NHGRI)",
                mode="markers",
                marker={"color": color, "size": 6, "symbol": "circle-open", "line": {"width": 1.5}},
                hovertemplate="%{x:.0f}<br>%{y:$,.4g} per MB<extra>Observed</extra>",
            )
        )
    figure.add_trace(
        go.Scatter(
            x=costs["year"],
            y=costs[value_column],
            name="Your scenario",
            mode="lines",
            line={"color": color, "width": 3},
            fill="tozeroy" if not use_log else None,
            fillcolor=color + "26" if not use_log else None,
            hovertemplate="%{x}<br>%{y:$,.4g} per MB<extra>Your scenario</extra>",
        )
    )
    last = costs.iloc[-1]
    _add_terminal_dot(figure, last["year"], last[value_column], color, surface, title)
    figure.update_layout(
        title=title,
        xaxis_title="Calendar year",
        yaxis_title="USD per MB",
        yaxis_type="log" if use_log else "linear",
        legend={"orientation": "h", "y": 1.14, "x": 0},
        margin={"l": 12, "r": 12, "t": 70, "b": 12},
        height=390,
        showlegend=show_context,
    )
    _apply_cost_axis_format(figure)
    return style_figure(figure, theme)


def sensitivity_chart(
    frame: pd.DataFrame, use_present_value: bool, theme: str = DEFAULT_THEME
) -> go.Figure:
    """Tornado chart: DNA's lifecycle cost swing as each input is perturbed
    +/-50% (or a fixed absolute step when the baseline is 0). Bars are
    sorted so the biggest lever renders at the top."""
    palette = palette_for(theme)
    color = technology_color("DNA", theme)
    figure = go.Figure()
    if frame.empty:
        figure.update_layout(title="What drives DNA's lifecycle cost", height=380)
        return style_figure(figure, theme)

    cost_label = "present value" if use_present_value else "lifecycle cost"
    baseline = float(frame["baseline_cost"].iloc[0])
    ordered = frame.sort_values("swing", ascending=True)
    lows = ordered[["low_cost", "high_cost"]].min(axis=1)
    highs = ordered[["low_cost", "high_cost"]].max(axis=1)

    hover_lines = []
    for row in ordered.itertuples():
        low_value_text, high_value_text = _format_pair(row.low_value, row.high_value)
        low_cost_text, high_cost_text = _format_money_pair(row.low_cost, row.high_cost)
        unit = f" {row.unit}" if row.unit else ""
        swing_text = "$" + _tick_labels([row.swing])[0] if row.swing > 0 else "no effect"
        hover_lines.append(
            f"At {low_value_text}{unit}: {low_cost_text} {cost_label}<br>"
            f"At {high_value_text}{unit}: {high_cost_text} {cost_label}<br>"
            f"Swing: {swing_text}"
        )

    figure.add_trace(
        go.Bar(
            x=highs - lows,
            y=ordered["parameter"],
            base=lows,
            orientation="h",
            marker_color=color,
            hovertext=hover_lines,
            hovertemplate="<b>%{y}</b><br>%{hovertext}<extra></extra>",
        )
    )
    figure.add_vline(
        x=baseline, line={"color": palette["figure"]["axis_line"], "width": 1.5, "dash": "dot"}
    )
    tickvals, axis_range = _linear_ticks(float(lows.min()), float(highs.max()), False)
    if tickvals:
        figure.update_xaxes(
            tickmode="array", tickvals=tickvals, ticktext=_tick_labels(tickvals), range=axis_range
        )
    figure.update_layout(
        title="What drives DNA's lifecycle cost",
        xaxis_title="Present value (USD)" if use_present_value else "Lifecycle cost (USD)",
        margin={"l": 12, "r": 12, "t": 70, "b": 12},
        height=120 + 34 * len(ordered),
        showlegend=False,
    )
    return style_figure(figure, theme)


def breakeven_chart(
    frame: pd.DataFrame, current_synthesis_cost: float, theme: str = DEFAULT_THEME
) -> go.Figure:
    """Current vs. break-even synthesis cost/MB per competing technology, on
    a log axis since the gap commonly spans several orders of magnitude.
    Technologies DNA cannot reach by lowering synthesis cost alone (its other
    costs already exceed them) get a marker instead of a bar."""
    palette = palette_for(theme)
    figure = go.Figure()
    if frame.empty:
        figure.update_layout(title="Synthesis cost needed to break even", height=320)
        return style_figure(figure, theme)

    reachable = frame[frame["breakeven_synthesis_cost_usd_per_mb"].notna()]
    unreachable = frame[frame["breakeven_synthesis_cost_usd_per_mb"].isna()]
    values = [current_synthesis_cost]
    if not reachable.empty:
        breakeven = reachable["breakeven_synthesis_cost_usd_per_mb"]
        lows = breakeven.clip(upper=current_synthesis_cost)
        highs = breakeven.clip(lower=current_synthesis_cost)
        # Zero-cost breakeven cannot sit on a log axis; nudge it to a visible
        # sliver instead of dropping the bar.
        floor = current_synthesis_cost / 1e6 if current_synthesis_cost > 0 else 1e-9
        lows = lows.clip(lower=floor)

        hover_lines = []
        for row in reachable.itertuples():
            breakeven_value = row.breakeven_synthesis_cost_usd_per_mb
            breakeven_text = "$" + _tick_labels([breakeven_value])[0] if breakeven_value > 0 else "$0"
            if breakeven_value <= 0:
                direction = f"DNA would need free synthesis (a floor of $0/MB) to match {row.technology}"
            elif breakeven_value < current_synthesis_cost:
                factor_text = _tick_labels([row.reduction_factor])[0]
                direction = (
                    f"Synthesis would need to fall {factor_text}x from today's "
                    f"price to match {row.technology}"
                )
            elif breakeven_value > current_synthesis_cost:
                headroom = breakeven_value / current_synthesis_cost if current_synthesis_cost > 0 else float("inf")
                headroom_text = _tick_labels([headroom])[0]
                direction = (
                    f"Already cheaper than {row.technology} today -- synthesis could rise "
                    f"{headroom_text}x before losing that edge"
                )
            else:
                direction = f"Exactly at parity with {row.technology} today"
            hover_lines.append(f"Break-even synthesis cost: {breakeven_text}/MB<br>{direction}")

        figure.add_trace(
            go.Bar(
                x=highs - lows,
                y=reachable["technology"],
                base=lows,
                orientation="h",
                marker_color=[technology_color(t, theme) for t in reachable["technology"]],
                hovertext=hover_lines,
                hovertemplate="<b>%{y}</b><br>%{hovertext}<extra></extra>",
            )
        )
        values.extend(lows.tolist())
        values.extend(highs.tolist())
    if not unreachable.empty:
        unreachable_lines = [
            f"Not reachable -- DNA's sequencing and retrieval costs alone already exceed "
            f"{technology}'s lifecycle cost, so no synthesis price (even $0/MB) would make "
            f"DNA cheaper"
            for technology in unreachable["technology"]
        ]
        figure.add_trace(
            go.Scatter(
                x=[current_synthesis_cost] * len(unreachable),
                y=unreachable["technology"],
                mode="markers",
                marker={"symbol": "x-thin", "size": 11, "color": palette["figure"]["tick_color"], "line": {"width": 2}},
                hovertext=unreachable_lines,
                hovertemplate="<b>%{y}</b><br>%{hovertext}<extra></extra>",
                showlegend=False,
            )
        )
    figure.add_vline(
        x=current_synthesis_cost,
        line={"color": palette["figure"]["axis_line"], "width": 1.5, "dash": "dash"},
    )
    figure.add_annotation(
        # Unlike shapes (add_vline above), annotations on a log-typed axis
        # take the position in log10 units, not the raw data value.
        x=math.log10(current_synthesis_cost) if current_synthesis_cost > 0 else 0.0,
        y=1.0,
        yref="paper",
        yanchor="bottom",
        text="Current",
        showarrow=False,
        xanchor="center",
        font={"color": palette["figure"]["tick_color"], "size": 12},
    )
    if min(values) == max(values):
        # Every technology is unreachable: nothing sets a real range, so
        # widen a couple of decades around "current" rather than rendering a
        # zero-width axis with the marker pinned to one edge.
        values = [max(current_synthesis_cost, 1e-9) / 100, max(current_synthesis_cost, 1e-9) * 100]
    tickvals, axis_range = _log_ticks(min(values), max(values))
    if tickvals:
        figure.update_xaxes(
            type="log",
            tickmode="array",
            tickvals=tickvals,
            ticktext=_tick_labels(tickvals),
            range=axis_range,
        )
    figure.update_layout(
        title="Synthesis cost needed to break even",
        xaxis_title="DNA synthesis cost (USD/MB)",
        margin={"l": 12, "r": 12, "t": 92, "b": 12},
        height=160 + 34 * len(frame),
        showlegend=False,
    )
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


def sensitivity_exports(
    frame: pd.DataFrame, use_present_value: bool, theme: str = DEFAULT_THEME
) -> dict[str, bytes]:
    return plotly_image_exports(sensitivity_chart(frame, use_present_value, theme))


def breakeven_exports(
    frame: pd.DataFrame, current_synthesis_cost: float, theme: str = DEFAULT_THEME
) -> dict[str, bytes]:
    return plotly_image_exports(breakeven_chart(frame, current_synthesis_cost, theme))
