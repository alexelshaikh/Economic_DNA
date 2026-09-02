import unittest

from economic_dna import Scenario, simulate_dna_unit_costs, simulate_scenario, simulate_start_years
from economic_dna.visualization import (
    breakdown_exports,
    dna_unit_cost_exports,
    lifecycle_exports,
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


if __name__ == "__main__":
    unittest.main()
