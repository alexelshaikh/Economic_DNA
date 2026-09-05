import math
import unittest
from dataclasses import replace

from economic_dna import (
    PRESET_SCENARIOS,
    Scenario,
    dna_cost_sensitivity,
    find_breakeven_synthesis_cost,
    find_lifecycle_crossovers,
    load_observed_sequencing_costs,
    simulate_dna_uncertainty_band,
    simulate_scenario,
    synthesis_historical_trend,
)


class HistoricalTrendTests(unittest.TestCase):
    def test_observed_sequencing_costs_are_sorted_and_positive(self):
        frame = load_observed_sequencing_costs()
        self.assertGreater(len(frame), 50)
        self.assertTrue((frame["sequencing_cost_usd_per_mb"] > 0).all())
        self.assertTrue(frame["year"].is_monotonic_increasing)
        # NHGRI data starts in 2001 and its cost/Mb has fallen by orders of
        # magnitude since then.
        self.assertLess(frame["year"].min(), 2002)
        self.assertGreater(frame["sequencing_cost_usd_per_mb"].iloc[0] / frame["sequencing_cost_usd_per_mb"].iloc[-1], 100)

    def test_synthesis_historical_trend_matches_the_editable_base_year_cost(self):
        scenario = Scenario()
        trend = synthesis_historical_trend(2000, scenario.dna_cost_base_year, points=2)
        final_cost = trend.iloc[-1]["synthesis_cost_usd_per_mb"]
        self.assertAlmostEqual(final_cost, scenario.dna_synthesis_cost_per_mb, delta=1e-6)

    def test_synthesis_historical_trend_declines_over_time(self):
        trend = synthesis_historical_trend(2000, 2026, points=27)
        self.assertTrue((trend["synthesis_cost_usd_per_mb"].diff().dropna() < 0).all())


class LifecycleCrossoverTests(unittest.TestCase):
    def test_finds_the_year_dna_overtakes_a_cheaper_technology(self):
        scenario = Scenario(technologies=("DNA", "Amazon Deep Archive"))
        result = simulate_scenario(scenario)
        crossovers = find_lifecycle_crossovers(result, use_present_value=False)
        self.assertIn("Amazon Deep Archive", crossovers)
        # DNA's synthesis write cost is front-loaded and enormous relative to
        # Amazon's per-request pricing at the paper baseline, so DNA never
        # catches up within a 100-year horizon.
        self.assertIsNone(crossovers["Amazon Deep Archive"])

    def test_returns_empty_when_reference_technology_is_not_selected(self):
        scenario = Scenario(technologies=("Amazon Deep Archive",))
        result = simulate_scenario(scenario)
        self.assertEqual(find_lifecycle_crossovers(result, use_present_value=False), {})

    def test_finds_a_year_when_a_flat_technology_is_eventually_overtaken(self):
        # DNA's cost within one scenario is essentially flat (one synthesis
        # write, no replacement inside the horizon); Tape's hardware line
        # item has 0% decline by default, so its cumulative cost grows
        # without bound. Given a long enough horizon it must cross a flat
        # line eventually -- this is the one case where a within-lifetime
        # crossover is reachable without contrived inputs.
        scenario = Scenario(
            technologies=("DNA", "Tape On-premise"),
            dna_synthesis_cost_per_mb=0.0001,
            annual_retrieval_percent=0.0,
            horizon_years=5000,
        )
        result = simulate_scenario(scenario)
        crossovers = find_lifecycle_crossovers(result, use_present_value=False)
        year = crossovers["Tape On-premise"]
        self.assertIsNotNone(year)

        yearly = result.yearly.pivot(index="year", columns="technology", values="cumulative_cost_usd")
        self.assertLessEqual(yearly.loc[year, "DNA"], yearly.loc[year, "Tape On-premise"])
        self.assertGreater(yearly.loc[year - 1, "DNA"], yearly.loc[year - 1, "Tape On-premise"])


class BreakevenSynthesisCostTests(unittest.TestCase):
    def test_reachable_breakeven_matches_a_brute_force_bisection(self):
        scenario = Scenario(annual_retrieval_percent=0.0)
        frame = find_breakeven_synthesis_cost(scenario, use_present_value=False)
        self.assertEqual(set(frame["technology"]), {"Amazon Deep Archive", "Azure Blob Archive", "Tape On-premise"})

        value_column = "total_cost_usd"

        def dna_minus_target(x: float, target_technology: str) -> float:
            probe = replace(scenario, dna_synthesis_cost_per_mb=x)
            totals = simulate_scenario(probe).totals.set_index("technology")[value_column]
            return float(totals["DNA"] - totals[target_technology])

        for _, row in frame.iterrows():
            technology = row["technology"]
            breakeven = row["breakeven_synthesis_cost_usd_per_mb"]
            self.assertIsNotNone(breakeven)
            # At the solved price the two totals must be equal; a hair below
            # or above should flip which technology is cheaper.
            self.assertAlmostEqual(dna_minus_target(breakeven, technology), 0.0, delta=1e-6)
            self.assertLess(dna_minus_target(breakeven * 0.999, technology), 0)
            self.assertGreater(dna_minus_target(breakeven * 1.001, technology), 0)

    def test_unreachable_breakeven_is_reported_as_none_not_a_wrong_number(self):
        # High annual retrieval makes DNA's sequencing cost alone (independent
        # of the synthesis price) exceed the paper-baseline competitors.
        scenario = Scenario()
        frame = find_breakeven_synthesis_cost(scenario, use_present_value=False)
        self.assertTrue(frame["breakeven_synthesis_cost_usd_per_mb"].isna().all())
        self.assertTrue(frame["reduction_factor"].isna().all())

    def test_empty_when_dna_is_not_selected(self):
        scenario = Scenario(technologies=("Amazon Deep Archive",))
        frame = find_breakeven_synthesis_cost(scenario, use_present_value=False)
        self.assertTrue(frame.empty)

    def test_reduction_factor_is_current_over_breakeven(self):
        scenario = Scenario(annual_retrieval_percent=0.0)
        frame = find_breakeven_synthesis_cost(scenario, use_present_value=False).set_index("technology")
        breakeven = frame.loc["Amazon Deep Archive", "breakeven_synthesis_cost_usd_per_mb"]
        expected = scenario.dna_synthesis_cost_per_mb / breakeven
        self.assertAlmostEqual(frame.loc["Amazon Deep Archive", "reduction_factor"], expected)


class SensitivityTests(unittest.TestCase):
    def test_tornado_rows_are_sorted_by_descending_swing(self):
        scenario = Scenario()
        frame = dna_cost_sensitivity(scenario, use_present_value=False)
        self.assertGreater(len(frame), 0)
        swings = frame["swing"].tolist()
        self.assertEqual(swings, sorted(swings, reverse=True))
        self.assertTrue((frame["swing"] >= 0).all())

    def test_archive_size_and_synthesis_cost_dominate_at_paper_baseline(self):
        # Archive size scales both the write and read cost linearly (the
        # single biggest lever); synthesis cost scales only the write
        # component, which happens to dominate DNA's total at the paper
        # baseline, so it is a close second. Everything else is noise by
        # comparison.
        scenario = Scenario()
        frame = dna_cost_sensitivity(scenario, use_present_value=False).set_index("parameter")
        self.assertEqual(list(frame.index[:2]), ["Archive size", "Synthesis cost"])
        self.assertGreater(frame.loc["Archive size", "swing"], frame.loc["Synthesis cost", "swing"])
        self.assertGreater(
            frame.loc["Synthesis cost", "swing"] / frame.loc["Sequencing cost", "swing"], 100
        )

    def test_average_asset_size_never_affects_dnas_own_cost(self):
        # DNA's cost formula is driven by archive_size_mb directly, not by
        # asset count, so this parameter is a true zero (not just small).
        scenario = Scenario()
        frame = dna_cost_sensitivity(scenario, use_present_value=False).set_index("parameter")
        self.assertEqual(frame.loc["Average asset size", "swing"], 0.0)

    def test_empty_when_dna_is_not_selected(self):
        scenario = Scenario(technologies=("Amazon Deep Archive",))
        frame = dna_cost_sensitivity(scenario, use_present_value=False)
        self.assertTrue(frame.empty)


class UncertaintyBandTests(unittest.TestCase):
    def test_band_is_ordered_and_brackets_the_deterministic_line(self):
        scenario = Scenario(horizon_years=20)
        band = simulate_dna_uncertainty_band(scenario, use_present_value=False, n_samples=500, seed=1)
        self.assertTrue((band["p10"] <= band["p50"]).all())
        self.assertTrue((band["p50"] <= band["p90"]).all())

        deterministic = simulate_scenario(scenario).yearly
        deterministic = deterministic[deterministic["technology"] == "DNA"].set_index("year")["cumulative_cost_usd"]
        # The deterministic (zero-spread) line uses the scenario's own decline
        # rates, i.e. the sampling distribution's center, so it should sit
        # comfortably inside the P10-P90 band almost everywhere.
        inside = (deterministic >= band.set_index("year")["p10"]) & (deterministic <= band.set_index("year")["p90"])
        self.assertGreater(inside.mean(), 0.9)

    def test_is_reproducible_given_the_same_seed(self):
        scenario = Scenario(horizon_years=10)
        first = simulate_dna_uncertainty_band(scenario, use_present_value=False, n_samples=50, seed=7)
        second = simulate_dna_uncertainty_band(scenario, use_present_value=False, n_samples=50, seed=7)
        self.assertTrue((first["p50"] == second["p50"]).all())

    def test_rejects_too_few_samples(self):
        with self.assertRaises(ValueError):
            simulate_dna_uncertainty_band(Scenario(), use_present_value=False, n_samples=1)


class PresetScenarioTests(unittest.TestCase):
    def test_every_preset_is_a_valid_ready_to_simulate_scenario(self):
        self.assertGreaterEqual(len(PRESET_SCENARIOS), 3)
        for name, scenario in PRESET_SCENARIOS.items():
            with self.subTest(preset=name):
                self.assertIsInstance(scenario, Scenario)
                result = simulate_scenario(scenario)
                self.assertGreater(len(result.yearly), 0)

    def test_paper_baseline_preset_matches_the_default_scenario(self):
        self.assertEqual(PRESET_SCENARIOS["Paper baseline"], Scenario())


if __name__ == "__main__":
    unittest.main()
