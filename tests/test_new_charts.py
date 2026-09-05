import math
import unittest

from economic_dna import (
    Scenario,
    dna_cost_sensitivity,
    find_breakeven_synthesis_cost,
    find_lifecycle_crossovers,
    load_observed_sequencing_costs,
    simulate_dna_unit_costs,
    simulate_dna_uncertainty_band,
    simulate_scenario,
    synthesis_historical_trend,
)
from economic_dna.visualization import (
    COMPONENT_PATTERNS,
    breakdown_chart,
    breakeven_chart,
    dna_unit_cost_chart,
    format_display_number,
    lifecycle_chart,
    sensitivity_chart,
)


class HistoricalOverlayChartTests(unittest.TestCase):
    def test_history_and_observed_traces_show_a_legend_and_widen_the_axis(self):
        scenario = Scenario(horizon_years=5)
        costs = simulate_dna_unit_costs(scenario, 2030)
        history = synthesis_historical_trend(2000, scenario.dna_cost_base_year, points=5)

        plain = dna_unit_cost_chart(costs, "synthesis_cost_usd_per_mb", "Synthesis", "#b44737", True)
        with_history = dna_unit_cost_chart(
            costs, "synthesis_cost_usd_per_mb", "Synthesis", "#b44737", True, history=history
        )

        self.assertFalse(plain.layout.showlegend)
        self.assertTrue(with_history.layout.showlegend)
        self.assertEqual(len(with_history.data), len(plain.data) + 1)
        # The earliest historical year is far before the scenario's own data;
        # the axis must stretch to include it, not clip it.
        historical_trace = next(t for t in with_history.data if t.name == "Historical trend (paper fit)")
        self.assertEqual(min(historical_trace.x), 2000)

    def test_observed_markers_use_the_real_nhgri_data(self):
        scenario = Scenario(horizon_years=5)
        costs = simulate_dna_unit_costs(scenario, 2030)
        observed = load_observed_sequencing_costs()
        figure = dna_unit_cost_chart(
            costs, "sequencing_cost_usd_per_mb", "Sequencing", "#1f77b4", True, observed=observed
        )
        observed_trace = next(t for t in figure.data if t.name == "Observed (NHGRI)")
        self.assertEqual(observed_trace.mode, "markers")
        self.assertEqual(len(observed_trace.x), len(observed))


class LifecycleOverlayChartTests(unittest.TestCase):
    def test_uncertainty_band_adds_a_filled_trace_behind_the_lines(self):
        scenario = Scenario(horizon_years=10)
        result = simulate_scenario(scenario)
        band = simulate_dna_uncertainty_band(scenario, False, n_samples=50, seed=0)

        plain = lifecycle_chart(result, False, True)
        with_band = lifecycle_chart(result, False, True, dna_uncertainty=band)

        self.assertEqual(len(with_band.data), len(plain.data) + 1)
        band_trace = with_band.data[0]
        self.assertEqual(band_trace.fill, "toself")
        self.assertEqual(len(band_trace.x), 2 * len(band))

    def test_reachable_crossover_draws_a_vline_and_annotation(self):
        scenario = Scenario(
            technologies=("DNA", "Tape On-premise"),
            dna_synthesis_cost_per_mb=0.0001,
            annual_retrieval_percent=0.0,
            horizon_years=5000,
        )
        result = simulate_scenario(scenario)
        crossovers = find_lifecycle_crossovers(result, False)
        figure = lifecycle_chart(result, False, True, crossovers=crossovers)

        self.assertEqual(len(figure.layout.shapes), 1)
        self.assertEqual(len(figure.layout.annotations), 1)
        self.assertIn("Tape On-premise", figure.layout.annotations[0].text)

    def test_no_crossover_adds_no_shapes(self):
        scenario = Scenario(technologies=("DNA", "Amazon Deep Archive"))
        result = simulate_scenario(scenario)
        crossovers = find_lifecycle_crossovers(result, False)
        figure = lifecycle_chart(result, False, True, crossovers=crossovers)
        self.assertEqual(len(figure.layout.shapes), 0)
        self.assertEqual(len(figure.layout.annotations), 0)


class ColorblindPatternTests(unittest.TestCase):
    def test_each_component_gets_a_distinct_pattern_alongside_its_color(self):
        scenario = Scenario()
        result = simulate_scenario(scenario)
        figure = breakdown_chart(result, True)
        patterns = {trace.name: trace.marker.pattern.shape for trace in figure.data}
        expected = {
            "Write & replacement": COMPONENT_PATTERNS["write_cost_usd"],
            "Retrieval": COMPONENT_PATTERNS["read_cost_usd"],
            "Storage & operation": COMPONENT_PATTERNS["maintenance_cost_usd"],
        }
        self.assertEqual(patterns, expected)
        self.assertEqual(len(set(COMPONENT_PATTERNS.values())), len(COMPONENT_PATTERNS))


class SensitivityChartTests(unittest.TestCase):
    def test_bars_are_based_at_the_lower_of_low_and_high_cost(self):
        scenario = Scenario()
        frame = dna_cost_sensitivity(scenario, False)
        figure = sensitivity_chart(frame, False)
        bar = figure.data[0]
        for base, width, low, high in zip(bar.base, bar.x, frame.sort_values("swing")["low_cost"], frame.sort_values("swing")["high_cost"]):
            self.assertAlmostEqual(base, min(low, high))
            self.assertAlmostEqual(base + width, max(low, high))

    def test_empty_frame_does_not_raise(self):
        scenario = Scenario(technologies=("Amazon Deep Archive",))
        frame = dna_cost_sensitivity(scenario, False)
        figure = sensitivity_chart(frame, False)
        self.assertEqual(len(figure.data), 0)

    def test_hover_text_names_the_unit_and_the_cost_basis_per_row(self):
        scenario = Scenario()
        frame = dna_cost_sensitivity(scenario, False)
        figure = sensitivity_chart(frame, False)
        hover_by_parameter = dict(zip(frame.sort_values("swing")["parameter"], figure.data[0].hovertext))

        synthesis_hover = hover_by_parameter["Synthesis cost"]
        self.assertIn("$/MB", synthesis_hover)
        self.assertIn("lifecycle cost", synthesis_hover)
        self.assertIn("Swing:", synthesis_hover)
        self.assertNotIn("e+", synthesis_hover)
        self.assertNotIn("e-", synthesis_hover)

        # average_asset_size never enters DNA's own cost formula, so its
        # swing is exactly zero -- the hover text should say so in words
        # rather than showing a $0 swing that reads like a rounding artifact.
        self.assertIn("no effect", hover_by_parameter["Average asset size"])

    def test_hover_text_uses_present_value_wording_when_discounted(self):
        scenario = Scenario()
        frame = dna_cost_sensitivity(scenario, True)
        figure = sensitivity_chart(frame, True)
        self.assertTrue(any("present value" in text for text in figure.data[0].hovertext))


class BreakevenChartTests(unittest.TestCase):
    def test_reachable_rows_get_bars_and_a_log_axis_current_marker(self):
        scenario = Scenario(annual_retrieval_percent=0.0)
        frame = find_breakeven_synthesis_cost(scenario, False)
        figure = breakeven_chart(frame, scenario.dna_synthesis_cost_per_mb)

        self.assertEqual(figure.layout.xaxis.type, "log")
        bar = next(t for t in figure.data if t.type == "bar")
        self.assertEqual(len(bar.x), len(frame))
        annotation = figure.layout.annotations[0]
        self.assertAlmostEqual(annotation.x, math.log10(scenario.dna_synthesis_cost_per_mb))

    def test_unreachable_rows_get_markers_not_bars_and_a_sane_axis(self):
        scenario = Scenario()
        frame = find_breakeven_synthesis_cost(scenario, False)
        figure = breakeven_chart(frame, scenario.dna_synthesis_cost_per_mb)

        self.assertFalse(any(t.type == "bar" for t in figure.data))
        marker_trace = next(t for t in figure.data if t.type == "scatter")
        self.assertEqual(marker_trace.marker.symbol, "x-thin")
        low, high = figure.layout.xaxis.range
        self.assertLess(low, high)

    def test_empty_frame_does_not_raise(self):
        figure = breakeven_chart(find_breakeven_synthesis_cost(Scenario(technologies=("Amazon Deep Archive",)), False), 100.0)
        self.assertEqual(len(figure.data), 0)

    def test_reachable_hover_explains_direction_without_raw_scientific_notation(self):
        scenario = Scenario(annual_retrieval_percent=0.0)
        frame = find_breakeven_synthesis_cost(scenario, False)
        figure = breakeven_chart(frame, scenario.dna_synthesis_cost_per_mb)
        bar = next(t for t in figure.data if t.type == "bar")

        for text in bar.hovertext:
            with self.subTest(text=text):
                self.assertIn("Break-even synthesis cost", text)
                self.assertIn("fall", text)
                self.assertIn("x from today's price", text)
                self.assertNotIn("e+", text)
                self.assertNotIn("e-", text)

    def test_hover_flips_wording_when_dna_is_already_cheaper(self):
        scenario = Scenario(annual_retrieval_percent=0.0, dna_synthesis_cost_per_mb=1e-6)
        frame = find_breakeven_synthesis_cost(scenario, False)
        figure = breakeven_chart(frame, scenario.dna_synthesis_cost_per_mb)
        bar = next(t for t in figure.data if t.type == "bar")

        for text in bar.hovertext:
            with self.subTest(text=text):
                self.assertIn("Already cheaper", text)
                self.assertIn("before losing that edge", text)

    def test_unreachable_hover_explains_why(self):
        scenario = Scenario()
        frame = find_breakeven_synthesis_cost(scenario, False)
        figure = breakeven_chart(frame, scenario.dna_synthesis_cost_per_mb)
        marker_trace = next(t for t in figure.data if t.type == "scatter")
        for text in marker_trace.hovertext:
            with self.subTest(text=text):
                self.assertIn("Not reachable", text)
                self.assertIn("sequencing and retrieval", text)


class FormatDisplayNumberTests(unittest.TestCase):
    def test_never_emits_raw_scientific_notation(self):
        for value in (0, 1, 16573.838543736347, 9.82e-5, 1.69e8, 6.62e8, 1e-11, 1e16):
            with self.subTest(value=value):
                text = format_display_number(value)
                self.assertNotIn("e+", text)
                self.assertNotIn("e-", text)

    def test_uses_unicode_superscripts_not_html_tags_at_extreme_magnitudes(self):
        text = format_display_number(9.82e-5)
        self.assertNotIn("<sup>", text)
        self.assertIn("×10", text)


if __name__ == "__main__":
    unittest.main()
