import math
import unittest

import numpy as np

import storage_service
from economic_dna import (
    Scenario,
    find_crossover_years,
    simulate_dna_unit_costs,
    simulate_scenario,
    simulate_start_years,
)


class ScenarioTests(unittest.TestCase):
    def test_paper_baseline_maps_to_original_parameters(self):
        scenario = Scenario()
        self.assertEqual(scenario.number_of_assets, 1000)
        self.assertEqual(scenario.annual_assets_retrieved, 10)

    def test_continuous_coefficients_are_exposed_as_annual_percentages(self):
        scenario = Scenario()
        self.assertAlmostEqual(scenario.synthesis_decline_percent / 100, 1 - math.exp(-0.1671))
        self.assertAlmostEqual(scenario.sequencing_decline_percent / 100, 1 - math.exp(-0.4791))

    def test_invalid_inputs_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "archive_size_tb"):
            Scenario(archive_size_tb=0)
        with self.assertRaisesRegex(ValueError, "below 100"):
            Scenario(amazon_decline_percent=100)

    def test_custom_name_is_required_only_when_custom_storage_is_selected(self):
        Scenario(custom_storage_name="", technologies=("DNA",))
        with self.assertRaisesRegex(ValueError, "custom_storage_name"):
            Scenario(custom_storage_name="", technologies=("Custom storage",))

    def test_query_parameter_round_trip(self):
        expected = Scenario(
            archive_size_tb=12.5,
            dna_synthesis_cost_per_mb=123.45,
            amazon_storage_usd_per_mb_month=0.000123,
            azure_retrieval_usd_per_mb=0.000456,
            tape_media_usd_per_tb=42.0,
            custom_storage_name="LTO library",
            custom_storage_cost_per_tb_year=7.5,
            technologies=("DNA", "Tape On-premise", "Custom storage"),
        )
        self.assertEqual(Scenario.from_mapping(expected.to_query_params()), expected)


class SimulationTests(unittest.TestCase):
    def setUp(self):
        self.scenario = Scenario()
        self.result = simulate_scenario(self.scenario)

    def test_horizon_is_number_of_charged_years(self):
        counts = self.result.yearly.groupby("technology")["year"].count()
        self.assertTrue((counts == 100).all())
        self.assertEqual(self.result.yearly["year"].min(), 2025)
        self.assertEqual(self.result.yearly["year"].max(), 2124)

    def test_components_sum_to_total(self):
        yearly = self.result.yearly
        component_total = yearly[
            ["write_cost_usd", "read_cost_usd", "maintenance_cost_usd"]
        ].sum(axis=1)
        np.testing.assert_allclose(component_total, yearly["total_cost_usd"])

    def test_tape_replacements_are_strictly_inside_horizon(self):
        tape = self.result.yearly.query("technology == 'Tape On-premise'")
        replacement_years = tape.loc[tape["write_cost_usd"] > 0, "year"].tolist()
        self.assertEqual(replacement_years, [2025, 2055, 2085, 2115])

    def test_discounting_reduces_present_value(self):
        discounted = simulate_scenario(Scenario(discount_rate_percent=3)).totals
        self.assertTrue((discounted["present_value_usd"] <= discounted["total_cost_usd"]).all())

    def test_numeric_engine_matches_reference_model(self):
        reference_models = {
            "DNA": storage_service.dna(n=1000, k=10, obj_size_mb=1000, start_year=2025, d=100),
            "Amazon Deep Archive": storage_service.amazon(
                n=1000, k=10, obj_size_mb=1000, start_year=2025, d=100
            ),
            "Azure Blob Archive": storage_service.azure(
                n=1000, k=10, obj_size_mb=1000, start_year=2025, d=100
            ),
            "Tape On-premise": storage_service.tape_on_premise_updated(
                n=1000, k=10, obj_size_mb=1000, start_year=2025, d=100
            ),
        }
        totals = self.result.totals.set_index("technology")["total_cost_usd"]
        for technology, model in reference_models.items():
            with self.subTest(technology=technology):
                self.assertAlmostEqual(totals[technology], float(model.eval().doit()), places=5)

    def test_start_year_projection_and_crossover_contract(self):
        projection = simulate_start_years(self.scenario, 2027)
        self.assertEqual(sorted(projection["start_year"].unique().tolist()), [2025, 2026, 2027])
        crossovers = find_crossover_years(projection)
        self.assertEqual(set(crossovers), set(self.scenario.technologies) - {"DNA"})

    def test_paper_baseline_crossover_years(self):
        projection = simulate_start_years(self.scenario, 2350)
        self.assertEqual(
            find_crossover_years(projection),
            {
                "Amazon Deep Archive": 2332,
                "Azure Blob Archive": 2321,
                "Tape On-premise": 2149,
            },
        )

    def test_editable_dna_unit_cost_trajectory(self):
        scenario = Scenario(
            dna_cost_base_year=2026,
            dna_synthesis_cost_per_mb=100,
            dna_sequencing_cost_per_mb=10,
            synthesis_decline_percent=20,
            sequencing_decline_percent=50,
        )
        costs = simulate_dna_unit_costs(scenario, 2027)
        self.assertEqual(costs["synthesis_cost_usd_per_mb"].tolist(), [100, 80])
        self.assertEqual(costs["sequencing_cost_usd_per_mb"].tolist(), [10, 5])

    def test_custom_storage_capacity_and_request_charges(self):
        scenario = Scenario(
            start_year=2026,
            horizon_years=2,
            custom_storage_name="Journalist service",
            custom_cost_base_year=2026,
            custom_write_cost_per_tb=10,
            custom_write_cost_per_asset=0.01,
            custom_storage_cost_per_tb_year=2,
            custom_retrieval_cost_per_tb=5,
            custom_retrieval_cost_per_asset=0.1,
            technologies=("Custom storage",),
        )
        total = simulate_scenario(scenario).totals.iloc[0]
        self.assertEqual(total["technology"], "Journalist service")
        self.assertAlmostEqual(total["write_cost_usd"], 20)
        self.assertAlmostEqual(total["maintenance_cost_usd"], 4)
        self.assertAlmostEqual(total["read_cost_usd"], 2.1)
        self.assertAlmostEqual(total["total_cost_usd"], 26.1)

    def test_amazon_prices_are_editable(self):
        scenario = Scenario(
            horizon_years=1,
            amazon_price_base_year=2025,
            amazon_put_usd_per_request=1,
            amazon_bulk_restore_usd_per_request=2,
            amazon_bulk_retrieval_usd_per_mb=0.001,
            amazon_storage_usd_per_mb_month=0.001,
            amazon_decline_percent=0,
            technologies=("Amazon Deep Archive",),
        )
        total = simulate_scenario(scenario).totals.iloc[0]
        self.assertAlmostEqual(total["write_cost_usd"], 1_000)
        self.assertAlmostEqual(total["read_cost_usd"], 30)
        self.assertAlmostEqual(total["maintenance_cost_usd"], 12_000)

    def test_azure_prices_are_editable(self):
        scenario = Scenario(
            horizon_years=1,
            azure_price_base_year=2025,
            azure_write_usd_per_request=1,
            azure_read_usd_per_request=2,
            azure_retrieval_usd_per_mb=0.001,
            azure_storage_usd_per_mb_month=0.001,
            azure_decline_percent=0,
            technologies=("Azure Blob Archive",),
        )
        total = simulate_scenario(scenario).totals.iloc[0]
        self.assertAlmostEqual(total["write_cost_usd"], 1_000)
        self.assertAlmostEqual(total["read_cost_usd"], 30)
        self.assertAlmostEqual(total["maintenance_cost_usd"], 12_000)

    def test_tape_prices_are_editable(self):
        scenario = Scenario(
            horizon_years=2,
            tape_price_base_year=2025,
            tape_media_usd_per_tb=10,
            tape_hardware_usd_per_tb=20,
            tape_energy_usd_per_tb_year=3,
            tape_media_decline_percent=0,
            tape_hardware_decline_percent=0,
            tape_energy_decline_percent=0,
            tape_durability_years=10,
            technologies=("Tape On-premise",),
        )
        total = simulate_scenario(scenario).totals.iloc[0]
        self.assertAlmostEqual(total["write_cost_usd"], 10)
        self.assertAlmostEqual(total["maintenance_cost_usd"], 10)
        self.assertAlmostEqual(total["total_cost_usd"], 20)


if __name__ == "__main__":
    unittest.main()
