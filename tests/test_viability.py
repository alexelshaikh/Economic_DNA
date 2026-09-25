from dataclasses import replace
import unittest

import numpy as np

from economic_dna import (
    DNA_SENSITIVITY_PARAMETERS, Scenario, dna_cost_advantage,
    find_breakeven_synthesis_cost, simulate_scenario,
)
from economic_dna.visualization import viability_chart


class ViabilityTests(unittest.TestCase):
    def test_every_driver_matches_full_simulations_of_both_models(self):
        scenario = Scenario(horizon_years=30, discount_rate_percent=2.0)
        for spec in DNA_SENSITIVITY_PARAMETERS:
            for discounted in (False, True):
                with self.subTest(driver=spec.field, discounted=discounted):
                    result = dna_cost_advantage(scenario, spec.field, "Amazon Deep Archive", discounted, focus=False)
                    self.assertLess(len(result.curve), 160)
                    self.assertTrue(np.isfinite(result.curve["gain_usd"]).all())
                    for row in result.curve.iloc[[0, len(result.curve) // 2, -1]].itertuples():
                        value = int(row.parameter_value) if spec.field == "dna_durability_years" else row.parameter_value
                        totals = simulate_scenario(replace(scenario, **{spec.field: value})).totals.set_index("technology_key")
                        column = "present_value_usd" if result.use_present_value else "total_cost_usd"
                        dna, other = totals.loc["DNA", column], totals.loc["Amazon Deep Archive", column]
                        np.testing.assert_allclose([row.dna_cost_usd, row.comparison_cost_usd, row.gain_usd], [dna, other, other - dna], rtol=1e-11, atol=1e-6)

    def test_synthesis_crossings_match_existing_break_even_analysis(self):
        for discounted in (False, True):
            scenario = Scenario(annual_retrieval_percent=0.001, discount_rate_percent=3.0)
            expected = find_breakeven_synthesis_cost(scenario, discounted).set_index("technology")
            for comparison in expected.index:
                result = dna_cost_advantage(scenario, "dna_synthesis_cost_per_mb", comparison, discounted)
                self.assertEqual(len(result.crossings), 1)
                crossing = result.crossings[0]
                self.assertAlmostEqual(crossing.value, expected.loc[comparison, "breakeven_synthesis_cost_usd_per_mb"], places=14)
                at_crossing = result.curve.loc[result.curve["parameter_value"] == crossing.value].iloc[0]
                self.assertAlmostEqual(at_crossing.gain_usd, 0.0, places=8)
                self.assertGreater(result.curve.gain_usd.max(), 0)
                self.assertLess(result.curve.gain_usd.min(), 0)

    def test_baseline_has_no_spurious_price_crossing(self):
        result = dna_cost_advantage(Scenario(), "dna_synthesis_cost_per_mb", "Tape On-premise", False)
        self.assertEqual(result.crossings, ())
        self.assertTrue((result.curve.gain_usd < 0).all())

    def test_sequencing_price_crossing_and_free_price_endpoint(self):
        scenario = Scenario(dna_synthesis_cost_per_mb=0)
        result = dna_cost_advantage(scenario, "dna_sequencing_cost_per_mb", "Tape On-premise", False)
        self.assertEqual(len(result.crossings), 1)
        self.assertGreater(result.crossings[0].value, 0)
        free_tape = replace(scenario, tape_media_usd_per_tb=0, tape_hardware_usd_per_tb=0, tape_energy_usd_per_tb_year=0)
        endpoint = dna_cost_advantage(free_tape, "dna_sequencing_cost_per_mb", "Tape On-premise", False)
        self.assertEqual(endpoint.crossings[0].value, 0)

    def test_zero_retrieval_makes_sequencing_price_irrelevant(self):
        result = dna_cost_advantage(Scenario(annual_retrieval_percent=0), "dna_sequencing_cost_per_mb", "Tape On-premise", False)
        self.assertFalse(result.crossings)
        self.assertEqual(result.curve.gain_usd.nunique(), 1)

    def test_asset_size_changes_request_costs_but_not_dna_costs(self):
        result = dna_cost_advantage(Scenario(), "average_asset_size_mb", "Azure Blob Archive", False, focus=False)
        self.assertEqual(result.curve.dna_cost_usd.nunique(), 1)
        self.assertGreater(result.curve.comparison_cost_usd.nunique(), 1)

    def test_discount_sweep_uses_present_value_even_at_zero_current_discount(self):
        result = dna_cost_advantage(Scenario(), "discount_rate_percent", "Tape On-premise", False)
        self.assertTrue(result.use_present_value)
        self.assertGreater(result.curve.gain_usd.nunique(), 1)

    def test_nonmonotonic_decline_marks_both_crossings(self):
        scenario = Scenario(horizon_years=30, dna_durability_years=1, dna_cost_base_year=2026,
                            dna_synthesis_cost_per_mb=1e-6, dna_sequencing_cost_per_mb=0,
                            technologies=("DNA", "Custom storage"), custom_write_cost_per_tb=10)
        result = dna_cost_advantage(scenario, "synthesis_decline_percent", "Custom storage", False)
        self.assertEqual(len(result.crossings), 2)
        for crossing in result.crossings:
            self.assertFalse(crossing.discrete)
            gain = result.curve.loc[result.curve.parameter_value == crossing.value, "gain_usd"].iloc[0]
            self.assertAlmostEqual(gain, 0, places=8)

    def test_overflowing_sweep_extremes_are_gaps_not_false_crossings(self):
        scenario = Scenario(dna_cost_base_year=2500, synthesis_decline_percent=0, sequencing_decline_percent=0)
        result = dna_cost_advantage(scenario, "synthesis_decline_percent", "Tape On-premise", False)
        self.assertFalse(result.crossings)
        self.assertTrue(result.curve.gain_usd.isna().any())
        self.assertTrue(result.curve.gain_usd.notna().any())
        figure = viability_chart(result)
        self.assertTrue(np.isfinite(figure.layout.yaxis.range).all())
        self.assertFalse(figure.data[1].connectgaps)

    def test_zero_gain_everywhere_has_no_arbitrary_crossing(self):
        scenario = Scenario(technologies=("DNA", "Custom storage"), dna_synthesis_cost_per_mb=0,
                            dna_sequencing_cost_per_mb=0)
        result = dna_cost_advantage(scenario, "annual_retrieval_percent", "Custom storage", False)
        self.assertFalse(result.crossings)
        self.assertTrue((result.curve.gain_usd == 0).all())
        self.assertTrue(np.isfinite(viability_chart(result).layout.yaxis.range).all())

    def test_durability_step_is_not_mislabeled_as_exact_parity(self):
        scenario = Scenario(start_year=2026, horizon_years=30, synthesis_decline_percent=0, dna_synthesis_cost_per_mb=1e-6,
                            dna_sequencing_cost_per_mb=0, tape_price_base_year=2026, tape_media_decline_percent=0,
                            tape_durability_years=1000, tape_energy_usd_per_tb_year=0, tape_hardware_usd_per_tb=0)
        for tape_cost in (6.39, 6.0):
            result = dna_cost_advantage(replace(scenario, tape_media_usd_per_tb=tape_cost), "dna_durability_years", "Tape On-premise", False)
            self.assertEqual(len(result.crossings), 1)
            self.assertEqual(result.crossings[0].value, 5)
            self.assertTrue(result.crossings[0].discrete)
            self.assertTrue(all(float(v).is_integer() for v in result.curve.parameter_value))

    def test_custom_model_and_current_marker_survive_focus_changes(self):
        scenario = Scenario(technologies=("DNA", "Custom storage"), custom_storage_name="My archive", custom_write_cost_per_tb=1000,
                            annual_retrieval_percent=0)
        focused = dna_cost_advantage(scenario, "dna_synthesis_cost_per_mb", "Custom storage", False)
        full = dna_cost_advantage(scenario, "dna_synthesis_cost_per_mb", "Custom storage", False, focus=False)
        self.assertEqual(focused.comparison, "My archive")
        self.assertEqual(focused.crossings, full.crossings)
        self.assertLess(focused.curve.parameter_value.max(), full.curve.parameter_value.max())
        self.assertIn(full.current_value, full.curve.parameter_value.tolist())
        self.assertEqual(full.current_gain, focused.current_gain)

    def test_invalid_driver_and_comparison_are_rejected(self):
        for field, comparison in (("invalid", "Tape On-premise"), ("archive_size_tb", "DNA"), ("archive_size_tb", "Custom storage")):
            with self.assertRaises(ValueError):
                dna_cost_advantage(Scenario(), field, comparison, False)


class ViabilityChartTests(unittest.TestCase):
    def test_signed_axis_zero_line_crossing_and_both_gain_regions(self):
        result = dna_cost_advantage(Scenario(annual_retrieval_percent=0.001), "dna_synthesis_cost_per_mb", "Tape On-premise", False)
        for theme in ("light", "dark"):
            figure = viability_chart(result, theme)
            self.assertLess(figure.layout.yaxis.range[0], 0)
            self.assertGreater(figure.layout.yaxis.range[1], 0)
            self.assertTrue(any(s.y0 == s.y1 == 0 for s in figure.layout.shapes))
            self.assertTrue(any(s.x0 == s.x1 == result.crossings[0].value for s in figure.layout.shapes))
            self.assertIn("Break-even", figure.layout.annotations[0].text)
            self.assertTrue(any(v and v > 0 for v in figure.data[0].y))
            self.assertTrue(any(v and v < 0 for v in figure.data[1].y))

    def test_no_crossing_and_current_input_marker(self):
        result = dna_cost_advantage(Scenario(), "dna_synthesis_cost_per_mb", "Tape On-premise", False)
        figure = viability_chart(result)
        self.assertFalse(figure.layout.annotations)
        marker = next(trace for trace in figure.data if trace.name == "Current inputs")
        self.assertEqual(marker.x[0], result.current_value)
        self.assertEqual(marker.y[0], result.current_gain)


if __name__ == "__main__":
    unittest.main()
