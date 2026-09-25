import unittest
import copy
import importlib.util
import sys
from dataclasses import replace
from unittest.mock import patch

import numpy as np
from streamlit.testing.v1 import AppTest

from economic_dna import (
    PRESET_SCENARIOS, Scenario, dna_cost_sensitivity, find_breakeven_synthesis_cost,
    simulate_dna_uncertainty_band, simulate_scenario, simulate_start_years,
)
from economic_dna.visualization import breakeven_chart, lifecycle_chart
from tests import test_app


class CalculationRegressions(unittest.TestCase):
    def test_yaml_prices_are_used_when_scenario_defaults_load(self):
        import economic_dna.scenario as scenario_module
        from economic_dna.assumptions import load_assumptions
        assumptions = copy.deepcopy(load_assumptions())
        assumptions["dna"]["editable_synthesis_usd_per_mb"] = 123.0
        assumptions["tape"]["media_usd_per_tb"] = 7.0
        spec = importlib.util.spec_from_file_location("economic_dna._yaml_defaults_test", scenario_module.__file__)
        module = importlib.util.module_from_spec(spec)
        with patch.dict(sys.modules, {spec.name: module}), patch("economic_dna.assumptions.load_assumptions", return_value=assumptions):
            spec.loader.exec_module(module)
        self.assertEqual(module.Scenario().dna_synthesis_cost_per_mb, 123.0)
        self.assertEqual(module.Scenario().tape_media_usd_per_tb, 7.0)

    def test_all_free_storage_has_a_visible_zero_cost_axis(self):
        result = simulate_scenario(Scenario(technologies=("Custom storage",)))
        figure = lifecycle_chart(result, False, True)
        self.assertEqual(figure.layout.yaxis.type, "linear")
        self.assertLessEqual(figure.layout.yaxis.range[0], 0)
        self.assertGreater(figure.layout.yaxis.range[1], 0)

    def test_linear_axis_includes_a_free_comparison_technology(self):
        result = simulate_scenario(Scenario(technologies=("DNA", "Custom storage")))
        self.assertEqual(lifecycle_chart(result, False, False).layout.yaxis.range[0], 0)

    def test_projection_matches_independent_simulations_with_all_cost_streams(self):
        scenario = Scenario(
            start_year=2030, horizon_years=81, discount_rate_percent=3.5,
            dna_durability_years=7, tape_durability_years=12,
            tape_hardware_decline_percent=4, tape_energy_decline_percent=9,
            custom_storage_name="My service", custom_write_cost_per_tb=20,
            custom_write_cost_per_asset=0.2, custom_storage_cost_per_tb_year=5,
            custom_retrieval_cost_per_tb=1, custom_retrieval_cost_per_asset=0.04,
            custom_decline_percent=7, custom_replacement_years=8,
            technologies=Scenario().technologies + ("Custom storage",),
        )
        for discount in (0.0, 3.5, 99.0):
            current = replace(scenario, discount_rate_percent=discount)
            projection = simulate_start_years(current, 2050)
            for year in (2030, 2031, 2040, 2050):
                actual = projection[projection.start_year == year].set_index("technology_key")
                expected = simulate_scenario(current.with_start_year(year)).totals.set_index("technology_key")
                columns = ["write_cost_usd", "read_cost_usd", "maintenance_cost_usd", "total_cost_usd", "present_value_usd"]
                np.testing.assert_allclose(actual[columns], expected[columns], rtol=1e-12, atol=1e-12)

    def test_sensitivity_respects_valid_input_limits(self):
        scenarios = (
            Scenario(synthesis_decline_percent=90, sequencing_decline_percent=90, discount_rate_percent=90),
            Scenario(dna_durability_years=1),
            Scenario(archive_size_tb=1e9, annual_retrieval_percent=10000),
            Scenario(archive_size_tb=0.001, average_asset_size_mb=1000),
        )
        fields = {"Synthesis cost": "dna_synthesis_cost_per_mb", "Synthesis decline rate": "synthesis_decline_percent", "DNA durability": "dna_durability_years"}
        for scenario in scenarios:
            with self.subTest(scenario=scenario):
                frame = dna_cost_sensitivity(scenario, True)
                self.assertTrue(np.isfinite(frame.select_dtypes("number")).all().all())
                for _, row in frame.iterrows():
                    if row.parameter not in fields:
                        continue
                    for prefix in ("low", "high"):
                        value = row[f"{prefix}_value"]
                        if row.parameter == "DNA durability":
                            value = int(value)
                        changed = replace(scenario, **{fields[row.parameter]: value})
                        expected = simulate_scenario(changed).totals.set_index("technology_key").loc["DNA", "present_value_usd"]
                        np.testing.assert_allclose(row[f"{prefix}_cost"], expected, rtol=1e-12)

    def test_preservation_preset_has_reachable_prices_for_all_alternatives(self):
        scenario = PRESET_SCENARIOS["1 PB preservation archive, rare retrieval"]
        frame = find_breakeven_synthesis_cost(scenario, False)
        self.assertEqual(len(frame), 3)
        self.assertTrue((frame.breakeven_synthesis_cost_usd_per_mb > 0).all())
        for row in frame.itertuples():
            matched = replace(scenario, dna_synthesis_cost_per_mb=row.breakeven_synthesis_cost_usd_per_mb)
            totals = simulate_scenario(matched).totals.set_index("technology")
            np.testing.assert_allclose(totals.loc["DNA", "total_cost_usd"], totals.loc[row.technology, "total_cost_usd"], rtol=1e-12)

    def test_free_synthesis_chart_has_finite_linear_axis(self):
        scenario = Scenario(dna_synthesis_cost_per_mb=0, annual_retrieval_percent=0)
        figure = breakeven_chart(find_breakeven_synthesis_cost(scenario, False), 0)
        self.assertEqual(figure.layout.xaxis.type, "linear")
        self.assertTrue(np.isfinite(figure.layout.xaxis.range).all())
        self.assertTrue(all(x >= 0 for x in figure.data[0].x))

    def test_tiny_break_even_prices_do_not_produce_negative_bars(self):
        scenario = Scenario(annual_retrieval_percent=0)
        frame = find_breakeven_synthesis_cost(scenario, False)
        figure = breakeven_chart(frame, scenario.dna_synthesis_cost_per_mb)
        np.testing.assert_allclose(np.asarray(figure.data[0].base, dtype=float), frame.breakeven_synthesis_cost_usd_per_mb)

    def test_uncertainty_windows_match_a_full_matrix(self):
        scenario = Scenario(horizon_years=777, dna_durability_years=200, discount_rate_percent=2)
        samples = 31
        rng = np.random.default_rng(0)
        years = np.arange(scenario.start_year, scenario.start_year + scenario.horizon_years)
        synth = np.clip(scenario.synthesis_decline_percent * (1 + rng.uniform(-0.3, 0.3, samples)), 0, 99.999)
        seq = np.clip(scenario.sequencing_decline_percent * (1 + rng.uniform(-0.3, 0.3, samples)), 0, 99.999)
        total = np.power(1 - synth[:, None] / 100, years - scenario.dna_cost_base_year)
        total *= scenario.dna_synthesis_cost_per_mb * scenario.archive_size_mb
        total *= (years - scenario.start_year) % scenario.dna_durability_years == 0
        total += np.power(1 - seq[:, None] / 100, years - scenario.dna_cost_base_year) * scenario.dna_sequencing_cost_per_mb * scenario.archive_size_mb * scenario.annual_retrieval_percent / 100
        total /= np.power(1 + scenario.discount_rate_percent / 100, years - scenario.start_year)
        expected = np.quantile(np.cumsum(total, axis=1), [0.1, 0.5, 0.9], axis=0).T
        actual = simulate_dna_uncertainty_band(scenario, True, n_samples=samples)
        np.testing.assert_allclose(actual[["p10", "p50", "p90"]], expected, rtol=1e-12)

    def test_fractional_durability_is_rejected(self):
        for field in ("horizon_years", "dna_durability_years", "tape_durability_years", "custom_replacement_years"):
            with self.subTest(field=field), self.assertRaises(ValueError):
                Scenario(**{field: 1.5})

    def test_overflowing_decline_assumptions_are_reported(self):
        with self.assertRaisesRegex(ValueError, "numeric range"):
            simulate_scenario(Scenario(dna_cost_base_year=2500, synthesis_decline_percent=99.0))


class AppRegressions(unittest.TestCase):
    def app(self):
        return AppTest.from_file(str(test_app.StreamlitAppTests.APP_PATH), default_timeout=20)

    def test_invalid_display_parameter_falls_back(self):
        app = self.app()
        app.query_params["projection_end"] = "invalid"
        app.run()
        self.assertFalse(app.exception)
        self.assertEqual(app.number_input(key="projection_end").value, 2350)

    def test_start_year_beyond_outlook_end_is_calculable(self):
        app = self.app().run()
        app.number_input(key="start_year_widget").set_value(2400)
        app.button(key="calculate_header").click().run()
        self.assertFalse(app.exception)
        self.assertEqual(app.session_state["committed_widgets"]["projection_end"], 2400)
        test_app.StreamlitAppTests._open_tab(app, "Start-year outlook")
        self.assertFalse(app.exception)
        self.assertEqual(len(app.tabs[1].get("plotly_chart")), 1)

    def test_views_are_prepared_once_and_csv_is_deferred(self):
        import streamlit as st
        st.cache_data.clear()
        with patch("pandas.DataFrame.to_csv", side_effect=AssertionError("CSV was generated before download")):
            app = self.app().run()
        self.assertFalse(app.exception)
        self.assertEqual(len(app.get("plotly_chart")), 8)

        with patch("economic_dna.simulate_start_years", side_effect=AssertionError("Projection cache missed")), patch("economic_dna.dna_cost_sensitivity", side_effect=AssertionError("Sensitivity cache missed")), patch("economic_dna.dna_cost_advantage", side_effect=AssertionError("Viability cache missed")):
            app.run()
        self.assertFalse(app.exception)

    def test_all_analysis_views_are_ready_for_client_navigation(self):
        app = self.app().run()
        for label, count in (("Start-year outlook", 1), ("DNA unit costs", 2), ("Sensitivity", 3), ("Assumptions", 0), ("About", 0), ("Lifecycle", 2)):
            test_app.StreamlitAppTests._open_tab(app, label)
            self.assertFalse(app.exception, label)
            panel = next(tab for tab in app.tabs if tab.label == label)
            self.assertEqual(len(panel.get("plotly_chart")), count, label)
            self.assertEqual(len(panel.get("download_button")), count, label)

    def test_preservation_scenario_with_free_synthesis_renders(self):
        app = self.app().run()
        app.button(key="preset_4").click().run()
        app.number_input(key="dna_synthesis_cost").set_value(0.0)
        app.button(key="calculate_header").click().run()
        test_app.StreamlitAppTests._open_tab(app, "Sensitivity")
        self.assertFalse(app.exception)
        self.assertIn("Already cheaper (free synthesis)", " ".join(m.value for m in app.markdown))

    def test_invalid_submission_keeps_previous_results(self):
        app = self.app().run()
        app.number_input(key="archive_value").set_value(0.001)
        app.number_input(key="asset_value").set_value(100.0)
        app.button(key="calculate_header").click().run()
        self.assertFalse(app.exception)
        self.assertTrue(app.error)
        self.assertEqual(app.metric[0].value, "1 TB")
        self.assertEqual(len(app.get("plotly_chart")), 8)
