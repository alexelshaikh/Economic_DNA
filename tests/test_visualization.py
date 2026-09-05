import math
import unittest

from economic_dna import Scenario, simulate_dna_unit_costs, simulate_scenario, simulate_start_years
from economic_dna.visualization import (
    breakdown_chart,
    breakdown_exports,
    dna_unit_cost_chart,
    dna_unit_cost_exports,
    lifecycle_chart,
    lifecycle_exports,
    projection_chart,
    projection_exports,
)


class ChartExportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.scenario = Scenario(horizon_years=5)
        cls.result = simulate_scenario(cls.scenario)
        cls.projection = simulate_start_years(cls.scenario, 2027)
        cls.dna_costs = simulate_dna_unit_costs(cls.scenario, 2030)

    def test_each_graph_has_valid_and_distinct_image_exports(self):
        result = self.result
        exports = {
            "lifecycle": lifecycle_exports(result, False, True),
            "breakdown": breakdown_exports(result, True),
            "projection": projection_exports(
                self.projection, False, True, result.metadata
            ),
            "dna_synthesis": dna_unit_cost_exports(
                self.dna_costs,
                "synthesis_cost_usd_per_mb",
                "DNA synthesis cost trajectory",
                "#b44737",
                True,
                result.metadata,
            ),
            "dna_sequencing": dna_unit_cost_exports(
                self.dna_costs,
                "sequencing_cost_usd_per_mb",
                "DNA sequencing cost trajectory",
                "#176b87",
                True,
                result.metadata,
            ),
        }

        png_files = []
        svg_files = []
        for graph_exports in exports.values():
            self.assertTrue(graph_exports["png"].startswith(b"\x89PNG\r\n\x1a\n"))
            self.assertIn(b"<svg", graph_exports["svg"][:500])
            png_files.append(graph_exports["png"])
            svg_files.append(graph_exports["svg"])

        self.assertEqual(len(set(png_files)), len(exports))
        self.assertEqual(len(set(svg_files)), len(exports))

    def _y_values(self, figure):
        values = []
        for trace in figure.data:
            if getattr(trace, "y", None) is not None:
                values.extend(float(value) for value in trace.y)
        return [value for value in values if value > 0]

    def _all_charts(self, log_scale):
        return {
            "lifecycle": lifecycle_chart(self.result, False, log_scale),
            "breakdown": breakdown_chart(self.result, log_scale),
            "projection": projection_chart(self.projection, False, log_scale),
            "dna_synthesis": dna_unit_cost_chart(
                self.dna_costs, "synthesis_cost_usd_per_mb", "Synthesis", "#b44737", log_scale
            ),
        }

    def test_no_axis_tick_label_carries_a_currency_symbol(self):
        for log_scale in (True, False):
            for name, figure in self._all_charts(log_scale).items():
                with self.subTest(chart=name, log=log_scale):
                    self.assertTrue(figure.layout.yaxis.ticktext)
                    for label in figure.layout.yaxis.ticktext:
                        self.assertNotIn("$", label)

    def test_log_axis_labels_every_gridline_with_a_decade(self):
        figure = lifecycle_chart(self.result, False, True)
        tick_labels = list(figure.layout.yaxis.ticktext)

        self.assertLessEqual(len(tick_labels), 8)
        # Plotly's own log ticks label minor gridlines as bare "2"/"5"; ours
        # only ever carry whole decades.
        self.assertNotIn("2", tick_labels)
        self.assertNotIn("5", tick_labels)
        for value in figure.layout.yaxis.tickvals:
            self.assertAlmostEqual(math.log10(value), round(math.log10(value)))

    def test_axis_range_is_pinned_to_the_top_tick_and_clips_nothing(self):
        for log_scale in (True, False):
            for name, figure in self._all_charts(log_scale).items():
                with self.subTest(chart=name, log=log_scale):
                    low, high = figure.layout.yaxis.range
                    tickvals = list(figure.layout.yaxis.tickvals)
                    if log_scale:
                        low, high = 10**low, 10**high
                    # The top edge is always a labelled gridline.
                    self.assertAlmostEqual(max(tickvals), high, delta=abs(high) * 1e-9)
                    values = self._y_values(figure)
                    self.assertGreaterEqual(min(values), low * (1 - 1e-9))
                    self.assertLessEqual(max(values), high * (1 + 1e-9))

    def test_very_low_synthesis_cost_keeps_the_unit_axis_readable(self):
        scenario = Scenario(dna_synthesis_cost_per_mb=0.000001, horizon_years=60)
        costs = simulate_dna_unit_costs(scenario, 2085)
        figure = dna_unit_cost_chart(
            costs, "synthesis_cost_usd_per_mb", "Synthesis", "#b44737", True
        )
        tick_labels = list(figure.layout.yaxis.ticktext)

        # Previously every minor gridline was labelled in ".1e" form, which
        # produced ~16 labels such as "5.0e-07".
        self.assertLessEqual(len(tick_labels), 8)
        self.assertGreaterEqual(len(tick_labels), 3)
        for label in tick_labels:
            self.assertRegex(label, r"^(\d+(\.\d+)?x)?10<sup>-?\d+</sup>$")

    def test_linear_and_log_axes_share_one_suffix_style(self):
        result = simulate_scenario(Scenario(horizon_years=100))
        linear = set(lifecycle_chart(result, False, False).layout.yaxis.ticktext)
        log = set(lifecycle_chart(result, False, True).layout.yaxis.ticktext)

        # d3's "~s" format used to render billions as "G" on linear axes while
        # the log axis said "B".
        for label in linear | log:
            self.assertNotIn("G", label)
            self.assertNotIn("k", label)


if __name__ == "__main__":
    unittest.main()
