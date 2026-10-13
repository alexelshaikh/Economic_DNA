import ast
from contextlib import redirect_stdout
from io import StringIO
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]


class PaperNotebookTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.notebook = json.loads((ROOT / "Figures_Paper_2026.ipynb").read_text(encoding="utf-8"))

    def test_all_code_cells_compile(self):
        for cell in self.notebook["cells"]:
            if cell["cell_type"] == "code":
                ast.parse("".join(cell["source"]))

    def test_inventory_has_two_setups_and_all_40_figures(self):
        tags = [cell["metadata"].get("tags", []) for cell in self.notebook["cells"]]
        self.assertEqual(sum("setup" in value for value in tags), 2)
        self.assertEqual(sum("figure" in value for value in tags), 40)

    def test_combined_parity_has_one_shared_legend(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        from scripts.paper_combined_parity_figure import SOURCE

        cell = next(c for c in self.notebook["cells"] if c.get("id") == "combined-parity-code")
        self.assertEqual("".join(cell["source"]), SOURCE)
        ns = {"__name__": "__main__"}
        with patch("IPython.display.display"), redirect_stdout(StringIO()):
            for setup in self.notebook["cells"]:
                if "setup" in setup.get("metadata", {}).get("tags", []):
                    exec("".join(setup["source"]), ns)
            exec(SOURCE.replace("SAVE = True", "SAVE = False", 1), ns)
        fig = ns["fig"]
        try:
            self.assertEqual(sum(ax.get_legend() is not None for ax in fig.axes), 1)
            self.assertEqual([text.get_text() for text in ns["shared_legend"].get_texts()],
                             ["50%", "55%", "60%", "65%"])
            for ax, frame, column in zip(ns["axes"], (ns["sequencing_data"], ns["retention_data"]),
                                         ("sequencing_decline_percent", "retention_years")):
                for line, rate in zip(ax.lines, ns["PARITY_SYNTHESIS_RATES"]):
                    expected = frame[frame.synthesis_decline_percent == rate].sort_values(column)
                    np.testing.assert_allclose(line.get_ydata(), expected.amazon_decline_percent, equal_nan=True)
                self.assertIsNone(ax.get_legend())
            self.assertGreater(ns["axes"][0].texts[0].get_position()[0], 0)
        finally:
            plt.close(fig)

    def test_s12_long_retention_includes_replacements(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from dataclasses import replace
        import pandas as pd

        ns = {"__name__": "__main__"}
        with patch("IPython.display.display"), redirect_stdout(StringIO()):
            for cell in self.notebook["cells"]:
                if "setup" in cell["metadata"].get("tags", []):
                    exec("".join(cell["source"]), ns)
            baseline = ns["BASE"]
            cell = next(c for c in self.notebook["cells"] if c["cell_type"] == "code"
                        and 'record("S12_' in "".join(c["source"]))
            exec("".join(cell["source"]).replace("SAVE = True", "SAVE = False", 1), ns)
        self.assertEqual(ns["BASE"], baseline)
        self.assertEqual(baseline.horizon_years, 100)
        self.assertEqual(len(ns["fig"].axes), 3)
        self.assertIn("2,000-Year Retention", ns["fig"]._suptitle.get_text())
        audit = ns["s12_replacements"]
        self.assertTrue((audit.retention_years == 2000).all())
        self.assertTrue((audit.replacement_write_usd > 0).all())
        self.assertTrue((audit.replacement_to_initial_write_ratio < 1e-30).all())
        for year, group in audit.groupby("deployment_year"):
            self.assertEqual(group.replacement_count.tolist(), [3, 1, 1])
            self.assertEqual(group.iloc[0].replacement_years,
                             "; ".join(str(year + offset) for offset in (500, 1000, 1500)))
        for key, scenario in ns["s12_scenarios"].items():
            self.assertEqual(scenario.horizon_years, 2000)
            expected = ns["dna_cost_sensitivity"](scenario, use_present_value=False)
            expected = expected[expected.parameter != "Discount rate"].reset_index(drop=True)
            actual = ns["s12"].query("scenario == @key")[expected.columns].reset_index(drop=True)
            pd.testing.assert_frame_equal(expected, actual)
        plt.close("all")
        for cell in self.notebook["cells"]:
            if "figure" in cell["metadata"].get("tags", []):
                self.assertTrue("".join(cell["source"]).startswith("SAVE = True"))

    def test_s10b_parity_matches_full_model(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from dataclasses import replace

        ns = {"__name__": "__main__"}
        with patch("IPython.display.display"), redirect_stdout(StringIO()):
            for cell in self.notebook["cells"]:
                if "setup" in cell["metadata"].get("tags", []):
                    exec("".join(cell["source"]), ns)
            cell = next(c for c in self.notebook["cells"] if c["cell_type"] == "code"
                        and 'record("S10b_' in "".join(c["source"]))
            exec("".join(cell["source"]).replace("SAVE = True", "SAVE = False", 1), ns)
        frame = ns["s10b"]
        self.assertEqual(set(frame.sequencing_decline_percent),
                         {50, 55, 60, 65, ns["BASE"].sequencing_decline_percent})
        for key, group in frame.groupby("scenario"):
            valid = group[group.root_exists]
            self.assertGreater(len(valid), 0)
            for row in valid.iloc[[0, len(valid) // 2, -1]].itertuples():
                scenario = replace(ns["s10b_scenarios"][key],
                    synthesis_decline_percent=row.synthesis_decline_percent,
                    amazon_decline_percent=row.amazon_decline_percent)
                costs = ns["totals"](scenario)
                self.assertAlmostEqual(costs["DNA"] / costs["Amazon Deep Archive"], 1, places=8)
            for row in group.iloc[[0, -1]].itertuples():
                scenario = replace(ns["s10b_scenarios"][key],
                    synthesis_decline_percent=row.synthesis_decline_percent)
                self.assertAlmostEqual(ns["totals"](scenario)["DNA"] / row.dna_cost_usd, 1, places=10)
        ax = ns["fig"].axes[0]
        self.assertEqual(ax.get_xlabel(), "DNA Synthesis Cost Annual Decline (%)")
        self.assertEqual(ns["fig"].legends[0].get_title().get_text(), "DNA Sequencing Cost Annual Decline")
        ns["fig"].canvas.draw()
        plt.close("all")

    def test_s11_style_preserves_results(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import pandas as pd
        from scripts.paper_s10_style import style_s11

        ns = {"__name__": "__main__"}
        with patch("IPython.display.display"), redirect_stdout(StringIO()):
            for cell in self.notebook["cells"]:
                if "setup" in cell["metadata"].get("tags", []):
                    exec("".join(cell["source"]), ns)
            cell = next(c for c in self.notebook["cells"] if c["cell_type"] == "code"
                        and 'record("S11_' in "".join(c["source"]))
            source = "".join(cell["source"])
            self.assertEqual(style_s11(source), source)
            exec(source.replace("SAVE = True", "SAVE = False", 1), ns)
            expected, _ = ns["synthesis_parity_lines"]("retention", ns["s11_x"])
        pd.testing.assert_frame_equal(expected, ns["s11"])
        ax = ns["fig"].axes[0]
        self.assertEqual(ax.get_xscale(), "log")
        self.assertEqual(ax.xaxis.get_major_formatter()(1000, 0), "1,000")
        self.assertEqual(ax.yaxis.get_major_formatter()(10, 0), "10%")
        self.assertEqual(ax.get_xlabel().replace("\n", " "), "Retention (Years)")
        self.assertTrue(all(spine.get_visible() for spine in ax.spines.values()))
        self.assertFalse(any(line.get_visible() for line in ax.get_xgridlines() + ax.get_ygridlines()))
        self.assertEqual(ns["fig"].legends[0].get_title().get_text(), "Synthesis decline")
        self.assertEqual(ns["fig"].legends[0]._ncols, 1)
        self.assertEqual([t.get_text() for t in ns["fig"].legends[0].get_texts()], ["50%", "55%", "60%", "65%"])
        ns["fig"].canvas.draw()
        plt.close("all")

    def test_s10_style_preserves_results(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import pandas as pd

        ns = {"__name__": "__main__"}
        with patch("IPython.display.display"), redirect_stdout(StringIO()):
            for cell in self.notebook["cells"]:
                if "setup" in cell["metadata"].get("tags", []):
                    exec("".join(cell["source"]), ns)
            expected, _ = ns["synthesis_parity_lines"]("sequencing", ns["np"].linspace(0, 80, 321))
            cell = next(c for c in self.notebook["cells"] if c["cell_type"] == "code"
                        and 'record("S10_' in "".join(c["source"]))
            exec("".join(cell["source"]).replace("SAVE = True", "SAVE = False", 1), ns)
        pd.testing.assert_frame_equal(expected, ns["s10"])
        ax = ns["fig"].axes[0]
        self.assertTrue(all(spine.get_visible() for spine in ax.spines.values()))
        self.assertFalse(any(line.get_visible() for line in ax.get_xgridlines() + ax.get_ygridlines()))
        self.assertEqual(ax.xaxis.get_major_formatter()(10, 0), "10%")
        self.assertEqual(ax.yaxis.get_major_formatter()(10, 0), "10%")
        self.assertEqual(ax.xaxis.label.get_fontweight(), "bold")
        self.assertEqual(ns["fig"].legends[0].get_title().get_text(), "Synthesis decline")
        self.assertEqual(ns["fig"].legends[0]._ncols, 1)
        self.assertEqual([t.get_text() for t in ns["fig"].legends[0].get_texts()], ["50%", "55%", "60%", "65%"])
        ns["fig"].canvas.draw()
        plt.close("all")

    def test_s1_s2_s3_style_preserves_results(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import pandas as pd

        original = json.loads((ROOT / "Figures_Sensitivity_Paper.ipynb").read_text(encoding="utf-8"))
        namespace = {"__name__": "__main__"}
        with patch("IPython.display.display"), redirect_stdout(StringIO()):
            for cell in self.notebook["cells"]:
                if "setup" in cell["metadata"].get("tags", []):
                    exec("".join(cell["source"]), namespace)
            for number in (1, 2, 3):
                prefix = f'record("S{number}_'
                old_cell = next(c for c in original["cells"]
                                if c["cell_type"] == "code" and prefix in "".join(c["source"]))
                namespace["SAVE"] = False
                exec("".join(old_cell["source"]), namespace)
                stem = next(key for key in namespace["TABLES"] if key.startswith(f"S{number}_"))
                expected = namespace["TABLES"][stem]
                cell = next(c for c in self.notebook["cells"]
                            if c["cell_type"] == "code" and prefix in "".join(c["source"]))
                exec("".join(cell["source"]).replace("SAVE = True", "SAVE = False", 1), namespace)
                for key, frame in expected.items():
                    actual = namespace["TABLES"][stem][key]
                    if number == 2:
                        actual = actual[actual.scenario != "2100"].reset_index(drop=True)
                    pd.testing.assert_frame_equal(frame, actual)
                fig = namespace["FIGURES"][stem]
                if number == 3:
                    self.assertIn("Break-Even", fig._suptitle.get_text())
                    self.assertIn("Object size: 1,000 MB", fig._suptitle.get_text())
                    self.assertEqual(fig.legends[0].get_title().get_text(), "Annual Synthesis Cost Decline")
                    self.assertEqual(namespace["legend_note"].get_children()[0].get_fontstyle(), "italic")
                    self.assertIn("$y_0$", fig.axes[0].get_xlabel())
                if number == 2:
                    self.assertEqual(len(fig.axes), 3)
                    self.assertEqual(fig._suptitle.get_text(), "Break-even Synthesis Price by Data Retrieval Rate")
                    self.assertIn("$y_0$", fig.axes[0].get_title())
                    self.assertIn("2100", fig.axes[2].get_title())
                    self.assertIn("even free synthesis", namespace["legend_note"].get_text())
                    self.assertEqual(fig.legends[0].get_title().get_text(), "Comparison Storage Model")
                    scenario = namespace["s2_scenarios"]["2100"]
                    reference = namespace["threshold_at_start"](scenario)
                    actual = namespace["s2"]
                    actual = actual[(actual.scenario == "2100") &
                                    (actual.annual_retrieval_percent == scenario.annual_retrieval_percent)]
                    pd.testing.assert_frame_equal(reference.reset_index(drop=True),
                                                  actual.drop(columns="scenario").reset_index(drop=True))
                self.assertEqual(len(fig.legends), 1)
                self.assertTrue(fig.legends[0].get_frame_on())
                for ax in fig.axes:
                    self.assertTrue(all(s.get_visible() for s in ax.spines.values()))
                    self.assertFalse(any(line.get_visible() for line in ax.get_ygridlines()))
                    self.assertEqual(ax.xaxis.label.get_fontweight(), "bold")
                    self.assertEqual([line.get_marker() for line in ax.lines[:3]], ["o", "s", "^"])
        plt.close("all")

    def test_new_retrieval_figures_match_public_model(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        from dataclasses import replace

        namespace = {"__name__": "__main__"}
        with patch("IPython.display.display"), redirect_stdout(StringIO()):
            for cell in self.notebook["cells"]:
                if "setup" in cell["metadata"].get("tags", []):
                    exec("".join(cell["source"]), namespace)
            for prefix in ('record("S2b_', 'record("S2c_'):
                cell = next(c for c in self.notebook["cells"]
                            if c["cell_type"] == "code" and prefix in "".join(c["source"]))
                exec("".join(cell["source"]).replace("SAVE = True", "SAVE = False", 1), namespace)
                self.assertEqual(len(namespace["fig"].axes), 3)
                title = namespace["fig"]._suptitle.get_text()
                self.assertTrue(title.startswith("Break-Even"))
                self.assertIn("Archive: 1 TB", title)
                self.assertIn("Retention: 100 years", title)
                self.assertIn("$y_0$", namespace["fig"].axes[0].get_title())
                note_text = namespace["legend_note"].get_children()[0]
                self.assertEqual(note_text.get_fontstyle(), "italic")
                self.assertLess(note_text.get_fontsize(), namespace["shared_legend"].get_title().get_fontsize())
                for ax in namespace["fig"].axes:
                    self.assertEqual(ax.xaxis.get_major_formatter()(100), "100%")
            thresholds = namespace["s2b"]
            self.assertTrue(thresholds[thresholds.start_year == 2026].synthesis_decline_percent.isna().all())
            for (year, tech), group in thresholds[thresholds.status == "parity"].groupby(["start_year", "technology"]):
                for row in group.iloc[[0, len(group) // 2, -1]].itertuples():
                    scenario = replace(namespace["BASE"], start_year=int(year),
                                       annual_retrieval_percent=row.annual_retrieval_percent,
                                       synthesis_decline_percent=row.synthesis_decline_percent)
                    costs = namespace["totals"](scenario)
                    self.assertAlmostEqual(costs["DNA"] / costs[tech], 1, places=8)
            sizes = namespace["s2c"]
            tape = sizes[sizes.technology == "Tape On-premise"]
            self.assertTrue(tape.empty)
            self.assertNotIn("Tape On-premise", namespace["axes"][0].get_legend_handles_labels()[1])
            for (year, tech), group in sizes[sizes.status == "parity"].groupby(["start_year", "technology"]):
                for row in group.iloc[[0, len(group) // 2, -1]].itertuples():
                    scenario = replace(namespace["BASE"], start_year=int(year),
                        annual_retrieval_percent=row.annual_retrieval_percent,
                        average_asset_size_mb=row.object_size_mb)
                    costs = namespace["totals"](scenario)
                    self.assertAlmostEqual(costs["DNA"] / costs[tech], 1, places=8)
                    smaller = namespace["totals"](replace(scenario, average_asset_size_mb=row.object_size_mb / 2))
                    self.assertLess(smaller["DNA"], smaller[tech])
        plt.close("all")

    def test_s2d_uses_one_2050_price_calibration(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from dataclasses import replace

        namespace = {"__name__": "__main__"}
        with patch("IPython.display.display"), redirect_stdout(StringIO()):
            for cell in self.notebook["cells"]:
                if "setup" in cell["metadata"].get("tags", []):
                    exec("".join(cell["source"]), namespace)
            cell = next(c for c in self.notebook["cells"] if c["cell_type"] == "code"
                        and 'record("S2d_' in "".join(c["source"]))
            exec("".join(cell["source"]).replace("SAVE = True", "SAVE = False", 1), namespace)
        base = namespace["calibrated_base"]
        self.assertIn("2050 calibration: synthesis", namespace["fig"]._suptitle.get_text())
        self.assertIn("Reference object size: 1,000 MB", namespace["fig"]._suptitle.get_text())
        self.assertIn("$y_0$", namespace["fig"].axes[0].get_title())
        costs = namespace["totals"](base)
        self.assertAlmostEqual(costs["DNA"] / costs["Amazon Deep Archive"], 1, places=10)
        self.assertEqual([s.start_year for s in namespace["scenarios"].values()], [2050, 2075, 2100])
        for scenario in namespace["scenarios"].values():
            self.assertEqual(scenario.dna_synthesis_cost_per_mb, base.dna_synthesis_cost_per_mb)
            self.assertEqual(scenario.synthesis_decline_percent, namespace["BASE"].synthesis_decline_percent)
        thresholds = namespace["s2d"]
        point = thresholds[(thresholds.start_year == 2050) &
                           (thresholds.technology == "Amazon Deep Archive") &
                           (thresholds.annual_retrieval_percent == base.annual_retrieval_percent)].iloc[0]
        self.assertAlmostEqual(point.object_size_mb / base.average_asset_size_mb, 1, places=8)
        for (year, tech), group in thresholds.groupby(["start_year", "technology"]):
            for row in group.iloc[[0, len(group) // 2, -1]].itertuples():
                scenario = replace(base, start_year=int(year), annual_retrieval_percent=row.annual_retrieval_percent,
                                   average_asset_size_mb=row.object_size_mb if row.status == "parity" else base.archive_size_mb)
                costs = namespace["totals"](scenario)
                if row.status == "parity":
                    self.assertAlmostEqual(costs["DNA"] / costs[tech], 1, places=8)
                else:
                    self.assertLessEqual(costs["DNA"], costs[tech])
        plt.close("all")

    def test_m1_deployment_years_and_label_order(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        namespace = {"__name__": "__main__"}
        with patch("IPython.display.display"), redirect_stdout(StringIO()):
            for cell in self.notebook["cells"]:
                if "setup" in cell["metadata"].get("tags", []):
                    exec("".join(cell["source"]), namespace)
            cell = next(c for c in self.notebook["cells"]
                        if "figure" in c["metadata"].get("tags", [])
                        and "M1_dna_cost_sensitivity" in "".join(c["source"]))
            exec("".join(cell["source"]).replace("SAVE = True", "SAVE = False", 1), namespace)
        axes = namespace["FIGURES"]["M1_dna_cost_sensitivity"].axes
        self.assertEqual(len(axes), 3)
        self.assertEqual([s.start_year for s in namespace["m1_scenarios"].values()],
                         [2026, 2050, 2100])
        self.assertEqual([t.get_text() for t in axes[0].get_yticklabels()], [
            "Synthesis cost", "Synthesis decline rate", "Sequencing cost",
            "Sequencing decline rate", "DNA durability", "Archive size",
            "Average object size", "Data retrieval rate",
        ])
        for ax, year in zip(axes, [2026, 2050, 2100]):
            self.assertIn(str(year), ax.get_title())
            self.assertEqual(len(ax.collections), 16)
        self.assertEqual(len(namespace["m1"]), 24)
        plt.close("all")

    def test_f2a_reference_style_and_parity_values(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import math

        namespace = {"__name__": "__main__"}
        with patch("IPython.display.display"), redirect_stdout(StringIO()):
            for cell in self.notebook["cells"]:
                if "setup" in cell["metadata"].get("tags", []):
                    exec("".join(cell["source"]), namespace)
            cell = next(c for c in self.notebook["cells"]
                        if "figure" in c["metadata"].get("tags", [])
                        and "F2a_unit_costs" in "".join(c["source"]))
            exec("".join(cell["source"]).replace("SAVE = True", "SAVE = False", 1), namespace)
        ax = namespace["FIGURES"]["F2a_unit_costs"].axes[0]
        lines = [line for line in ax.lines if not line.get_label().startswith("_")]
        self.assertEqual([line.get_marker() for line in lines[:5]], ["o", "^", "*", "d", "x"])
        self.assertEqual(lines[-1].get_color(), namespace["PALETTE"][2])
        for line, percent in zip(lines[:5], namespace["annual_declines"]):
            self.assertAlmostEqual(line.get_ydata()[1] / line.get_ydata()[0], 1 - percent / 100)
        self.assertEqual(namespace["annual_declines"][-2:].tolist(), [70, 80])
        self.assertTrue(math.isnan(namespace["s4_synthesis_parity"](10.0, start_year=2026)))
        self.assertEqual(namespace["annual_declines"][-3:].tolist(), [60, 70, 80])
        self.assertTrue(ax.get_legend().get_frame_on())
        self.assertFalse(any(line.get_visible() for line in ax.get_ygridlines()))
        self.assertTrue(all(spine.get_visible() for spine in ax.spines.values()))
        crossings = namespace["TABLES"]["F2a_unit_costs"]["crossings"]
        self.assertEqual(len(crossings), 4)
        base = namespace["BASE"]
        for row in crossings.itertuples():
            dna = base.dna_synthesis_cost_per_mb * math.exp(-row.coefficient * (row.parity_year - base.dna_cost_base_year))
            amazon = base.amazon_put_usd_per_request / base.average_asset_size_mb * (1 - base.amazon_decline_percent / 100) ** (row.parity_year - base.amazon_price_base_year)
            self.assertAlmostEqual(dna / amazon, 1, places=10)
        plt.close("all")

    def test_f2bc_reference_style_and_costs(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np

        namespace = {"__name__": "__main__"}
        with patch("IPython.display.display"), redirect_stdout(StringIO()):
            for cell in self.notebook["cells"]:
                if "setup" in cell["metadata"].get("tags", []):
                    exec("".join(cell["source"]), namespace)
            for suffix in ("b", "c"):
                name = f"F2{suffix}_unit_costs"
                cell = next(c for c in self.notebook["cells"]
                            if "figure" in c["metadata"].get("tags", [])
                            and name in "".join(c["source"]))
                exec("".join(cell["source"]).replace("SAVE = True", "SAVE = False", 1), namespace)
                ax = namespace["FIGURES"][name].axes[-1] if suffix == "c" else namespace["FIGURES"][name].axes[0]
                lines = [line for line in ax.lines if not line.get_label().startswith("_")]
                self.assertEqual([line.get_marker() for line in lines[:4]], ["o", "*", "d", "x"])
                self.assertEqual(lines[0].get_color(), namespace["PALETTE"][1])
                self.assertTrue(np.isnan(namespace["s8_sequencing_parity"](10, start_year=2026)))
                self.assertEqual(namespace["annual_declines"][-3:].tolist(), [60, 70, 80])
                for line, percent in zip(lines[:4], namespace["annual_declines"]):
                    self.assertAlmostEqual(line.get_ydata()[1] / line.get_ydata()[0], 1 - percent / 100)
                self.assertEqual(lines[-1].get_color(), namespace["PALETTE"][2])
                self.assertTrue(ax.get_legend().get_frame_on())
                self.assertTrue(all(spine.get_visible() for spine in ax.spines.values()))
                self.assertFalse(any(line.get_visible() for line in ax.get_ygridlines()))
                tables = namespace["TABLES"][name]
                self.assertEqual(len(tables["crossings"]), 4)
                if suffix == "c":
                    axes = namespace["FIGURES"][name].axes
                    self.assertEqual(len(axes), 1)
                    self.assertEqual(len(ax.child_axes), 1)
                    zoom = ax.child_axes[0]
                    self.assertGreater(zoom.get_position().x0, ax.get_position().x0)
                    self.assertLess(zoom.get_position().x1, ax.get_position().x1)
                    self.assertGreater(zoom.get_position().y0, ax.get_position().y0)
                    self.assertLess(zoom.get_position().y1, ax.get_position().y1)
                    self.assertEqual(sum(isinstance(p, matplotlib.patches.ConnectionPatch) for p in ax.patches), 2)
                    self.assertEqual(zoom.get_yscale(), "log")
                    self.assertFalse(any("Year" in text.get_text() for text in ax.texts))
                    self.assertEqual(sum("Year" in text.get_text() for text in zoom.texts), 4)
                    self.assertFalse(any("%" in text.get_text() for text in zoom.texts))
                    for row in tables["crossings"].itertuples():
                        self.assertLess(zoom.get_xlim()[0], row.parity_year)
                        self.assertGreater(zoom.get_xlim()[1], row.parity_year)
                        self.assertLess(zoom.get_ylim()[0], row.cost_usd_per_mb)
                        self.assertGreater(zoom.get_ylim()[1], row.cost_usd_per_mb)
                    for original, detail in zip(lines, zoom.lines):
                        np.testing.assert_array_equal(original.get_ydata(), detail.get_ydata())
                        self.assertEqual(detail.get_markevery(), 1)
                        self.assertEqual(detail.get_markersize(), 7)
                base = namespace["BASE"]
                price = base.amazon_bulk_restore_usd_per_request / base.average_asset_size_mb + base.amazon_bulk_retrieval_usd_per_mb
                if suffix == "c":
                    price += base.amazon_storage_usd_per_mb_month * 12 / (base.annual_retrieval_percent / 100)
                years = tables["costs"].year
                np.testing.assert_allclose(tables["costs"].Amazon, price * (1 - base.amazon_decline_percent / 100) ** (years - base.amazon_price_base_year))
                for row in tables["crossings"].itertuples():
                    dna = base.dna_sequencing_cost_per_mb * np.exp(-row.coefficient * (row.parity_year - base.dna_cost_base_year))
                    cloud = price * (1 - base.amazon_decline_percent / 100) ** (row.parity_year - base.amazon_price_base_year)
                    self.assertAlmostEqual(dna / cloud, 1, places=10)
                self.assertEqual(tables["gaps"].year.tolist(), [base.start_year, 2100])
        plt.close("all")

    def test_workload_boundaries_match_full_model(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from dataclasses import replace

        ns = {"__name__": "__main__"}
        with patch("IPython.display.display"), redirect_stdout(StringIO()):
            for cell in self.notebook["cells"]:
                if "setup" in cell["metadata"].get("tags", []):
                    exec("".join(cell["source"]), ns)
            for suffix in ("e", "f", "g"):
                cell = next(c for c in self.notebook["cells"] if c["cell_type"] == "code"
                            and f'record("S4{suffix}_' in "".join(c["source"]))
                source = "".join(cell["source"]).replace("SAVE = True", "SAVE = False", 1)
                source = source.replace("MAX_RETENTION = 10000", "MAX_RETENTION = 120")
                exec(source, ns)
                table = ns["table"]
                expected_years = {2050, 2055, 2060, 2065, 2070, 2075, 2100} if suffix == "f" else {2050, 2075, 2100}
                self.assertEqual(set(table.deployment_year), expected_years)
                if suffix == "f":
                    late = table[table.deployment_year == 2075]
                    self.assertEqual(late[late.status == "parity"].retention_years.tolist(), [1, 2])
                    self.assertTrue((late[late.retention_years >= 3].status == "DNA_competitive_at_every_object_size").all())
                    self.assertIn("$y_0 = 2075$", [text.get_text() for text in ns["ax"].texts])
                if suffix == "e":
                    self.assertEqual(table[table.deployment_year == 2050].minimum_retention_years.unique().tolist(), [100])
                    for year, group in table.groupby("deployment_year"):
                        for row in group.iloc[[0, -1]].itertuples():
                            horizon = int(row.minimum_retention_years)
                            scenario = replace(ns["scenarios"][str(year)], horizon_years=horizon,
                                               archive_size_tb=row.archive_size_tb)
                            costs = ns["totals"](scenario)
                            self.assertLessEqual(costs["DNA"] / costs["Amazon Deep Archive"], 1 + 1e-12)
                            if horizon > 1:
                                before = ns["totals"](replace(scenario, horizon_years=horizon - 1))
                                self.assertGreater(before["DNA"], before["Amazon Deep Archive"])
                else:
                    for year, group in table[table.status == "parity"].groupby("deployment_year"):
                        for row in group.iloc[[0, len(group) // 2, -1]].itertuples():
                            updates = {"average_asset_size_mb": row.object_size_mb}
                            updates.update({"horizon_years": int(row.retention_years)} if suffix == "f"
                                           else {"archive_size_tb": row.archive_size_tb})
                            costs = ns["totals"](replace(ns["scenarios"][str(year)], **updates))
                            self.assertAlmostEqual(costs["DNA"] / costs["Amazon Deep Archive"], 1, places=8)
                exec(source.replace("USE_PARITY_CALIBRATION = True", "USE_PARITY_CALIBRATION = False"), ns)
                self.assertEqual(ns["workload_base"], ns["BASE"])
        plt.close("all")

    def test_s4c_s4d_earliest_deployments_and_parity(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from dataclasses import replace

        ns = {"__name__": "__main__"}
        with patch("IPython.display.display"), redirect_stdout(StringIO()):
            for cell in self.notebook["cells"]:
                if "setup" in cell["metadata"].get("tags", []):
                    exec("".join(cell["source"]), ns)
            for suffix in ("c", "d"):
                cell = next(c for c in self.notebook["cells"] if c["cell_type"] == "code"
                            and f'record("S4{suffix}_' in "".join(c["source"]))
                exec("".join(cell["source"]).replace("SAVE = True", "SAVE = False", 1), ns)
                audit = ns[f"s4{suffix}_audit"]
                self.assertFalse(audit.iloc[:-1].parity_possible.any())
                self.assertTrue(audit.iloc[-1].parity_possible)
                first = ns["first_viable_year"]
                self.assertEqual(first, 2029 if suffix == "c" else 2027)
                self.assertEqual(ns["deployment_years"], [first, first + 5, first + 10, first + 20])
                labels = [text.get_text() for text in ns["ax"].texts]
                self.assertTrue(all(f"$y_0 = {year}$" in labels for year in ns["deployment_years"]
                                    if not (suffix == "d" and year == 2027)))
                self.assertEqual(ns["ax"].get_xlabel().replace("\n", " "), "DNA Synthesis Cost Annual Decline (%)")
                self.assertIn("Cost Annual Decline (%)", ns["ax"].get_ylabel().replace("\n", " "))
                if suffix == "c":
                    self.assertIn("Deployment Year $y_0$", ns["ax"].get_title())
                else:
                    inline_labels = [text for text in ns["ax"].texts
                                     if text.get_text() in [f"$y_0 = {year}$" for year in ns["deployment_years"][1:]]]
                    self.assertTrue(all(0 < text.get_rotation() < 180 for text in inline_labels))
                    self.assertNotIn("$y_0 = 2027$", labels)
                    self.assertFalse(any(isinstance(c, matplotlib.collections.PathCollection)
                                         for c in ns["ax"].collections))
                self.assertEqual(len(ns["fig"].axes), 1)
                self.assertFalse(any(isinstance(c, matplotlib.collections.QuadMesh) for c in ns["ax"].collections))
                self.assertGreater(ns["ax"].get_xlim()[0], 0)
                self.assertLess(ns["SYNTHESIS_RANGE"][1], 100)
                for year, group in ns[f"s4{suffix}_parity"].groupby("deployment_year"):
                    if suffix == "d":
                        self.assertEqual(group.synthesis_decline_percent.max(), ns["SYNTHESIS_RANGE"][1])
                    for row in group.iloc[[0, len(group) // 2, -1]].itertuples():
                        changes = {"start_year": int(year), "synthesis_decline_percent": row.synthesis_decline_percent}
                        if suffix == "c":
                            changes["amazon_decline_percent"] = row.amazon_decline_percent
                        else:
                            changes["sequencing_decline_percent"] = row.sequencing_decline_percent
                        costs = ns["totals"](replace(ns["BASE"], **changes))
                        self.assertAlmostEqual(costs["DNA"] / costs["Amazon Deep Archive"], 1, places=6)
                for year in (first - 1, first):
                    changes = {"start_year": year, "synthesis_decline_percent": ns["SYNTHESIS_RANGE"][1]}
                    changes.update({"amazon_decline_percent": 0.0} if suffix == "c"
                                   else {"sequencing_decline_percent": ns["SEQUENCING_RANGE"][1]})
                    costs = ns["totals"](replace(ns["BASE"], **changes))
                    self.assertEqual(costs["DNA"] <= costs["Amazon Deep Archive"], year == first)
        plt.close("all")

    def test_s4b_joint_prices_match_full_lifecycle_model(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from dataclasses import replace

        ns = {"__name__": "__main__"}
        with patch("IPython.display.display"), redirect_stdout(StringIO()):
            for cell in self.notebook["cells"]:
                if "setup" in cell["metadata"].get("tags", []):
                    exec("".join(cell["source"]), ns)
            cell = next(c for c in self.notebook["cells"] if c["cell_type"] == "code"
                        and 'record("S4b_' in "".join(c["source"]))
            source = "".join(cell["source"]).replace("SAVE = True", "SAVE = False", 1)
            exec(source, ns)
        for row in ns["s4b"].iloc[[0, 1250, 20000, -1]].itertuples():
            probe = replace(ns["s4b_scenario"],
                dna_synthesis_cost_per_mb=row.synthesis_price_anchor_usd_per_mb,
                dna_sequencing_cost_per_mb=row.sequencing_price_anchor_usd_per_mb)
            costs = ns["totals"](probe)
            self.assertAlmostEqual(row.cost_ratio, costs["DNA"] / costs["Amazon Deep Archive"], places=10)
        for row in ns["s4b_parity"].iloc[[0, 100, -1]].itertuples():
            probe = replace(ns["s4b_scenario"],
                dna_synthesis_cost_per_mb=row.synthesis_price_deployment_usd_per_mb / ns["synthesis_factor"],
                dna_sequencing_cost_per_mb=row.sequencing_price_deployment_usd_per_mb / ns["sequencing_factor"])
            costs = ns["totals"](probe)
            self.assertAlmostEqual(costs["DNA"] / costs["Amazon Deep Archive"], 1, places=10)
        threshold = ns["threshold_at_start"](ns["s4b_scenario"]).set_index("technology").loc[
            "Amazon Deep Archive", "threshold_start_usd_per_mb"]
        self.assertAlmostEqual(ns["parity_at_baseline_read"] / threshold, 1, places=10)
        self.assertEqual(len(ns["fig"].axes), 1)
        self.assertIn("Cost parity", [text.get_text() for text in ns["ax"].texts])
        self.assertFalse(any(isinstance(c, matplotlib.collections.QuadMesh) for c in ns["ax"].collections))
        self.assertIn("2050", ns["ax"].get_xlabel())
        with patch("IPython.display.display"), redirect_stdout(StringIO()):
            ns["BASE"] = replace(ns["BASE"], annual_retrieval_percent=0)
            exec(source, ns)
        self.assertEqual(ns["read_coefficient"], 0)
        self.assertTrue((ns["s4b_parity"].synthesis_price_deployment_usd_per_mb == ns["synthesis_intercept"]).all())
        plt.close("all")

    def test_s4_style_preserves_heatmap(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import pandas as pd

        namespace = {"__name__": "__main__"}
        original = json.loads((ROOT / "Figures_Sensitivity_Paper.ipynb").read_text(encoding="utf-8"))
        with patch("IPython.display.display"), redirect_stdout(StringIO()):
            for cell in self.notebook["cells"]:
                if "setup" in cell["metadata"].get("tags", []):
                    exec("".join(cell["source"]), namespace)
            for cells in (original["cells"], self.notebook["cells"]):
                cell = next(c for c in cells if c["cell_type"] == "code"
                            and 'record("S4_joint_decline_sensitivity",' in "".join(c["source"]))
                namespace["SAVE"] = False
                source = "".join(cell["source"]).replace("SAVE = True", "SAVE = False", 1)
                source = source.replace("np.linspace(0, 70, 61)", "np.linspace(0, 70, 3)")
                source = source.replace("np.linspace(0, 20, 41)", "np.linspace(0, 20, 3)")
                exec(source, namespace)
                if cells is original["cells"]:
                    expected = namespace["s4"].copy()
            pd.testing.assert_frame_equal(expected, namespace["s4"])
        ax = namespace["ax"]
        self.assertTrue(all(s.get_visible() for s in ax.spines.values()))
        self.assertEqual(ax.xaxis.label.get_fontweight(), "bold")
        self.assertEqual(ax.xaxis.get_major_formatter()(20), "20%")
        self.assertEqual(namespace["mesh"].get_cmap().name, "RdBu_r")
        self.assertNotIn("Archive:", ax.get_title())
        plt.close("all")

    def test_2026_sensitivity_has_no_parity(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np

        namespace = {"__name__": "__main__"}
        with patch("IPython.display.display"), redirect_stdout(StringIO()):
            for cell in self.notebook["cells"]:
                if "setup" in cell["metadata"].get("tags", []):
                    exec("".join(cell["source"]), namespace)
            for name in ("S4_joint_decline_sensitivity_2026", "S8_sequencing_decline_sensitivity_2026"):
                cell = next(c for c in self.notebook["cells"]
                            if "figure" in c["metadata"].get("tags", []) and name in "".join(c["source"]))
                source = "".join(cell["source"]).replace("SAVE = True", "SAVE = False", 1)
                for original in ("np.linspace(0, 70, 61)", "np.linspace(0, 20, 41)",
                                 "np.linspace(0, 80, 161)", "np.linspace(0, 20, 81)"):
                    source = source.replace(original, original.rsplit(",", 1)[0] + ", 3)")
                exec(source, namespace)
                table = namespace["TABLES"][name]["grid"]
                self.assertTrue((table.dna_cost_usd > table.comparison_cost_usd).all())
                self.assertTrue((table.log10_cost_ratio > 0).all())
                self.assertTrue(all(s["start_year"] == 2026 for s in namespace["SCENARIOS"][name].values()))
                if name.startswith("S8"):
                    self.assertFalse(namespace["TABLES"][name]["parity_boundary"].root_exists.any())
            self.assertTrue(np.isnan(namespace["s4_synthesis_parity"](10, start_year=2026)))
            self.assertTrue(np.isnan(namespace["s8_sequencing_parity"](10, start_year=2026)))
            self.assertAlmostEqual(namespace["s4_synthesis_parity"](10, start_year=2050), 58.925818, places=5)
            self.assertAlmostEqual(namespace["s8_sequencing_parity"](10, start_year=2050), 25.711949, places=5)
        plt.close("all")

    def test_preview_is_write_free_and_can_be_saved_later(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        namespace = {"__name__": "__main__"}
        with tempfile.TemporaryDirectory() as directory, patch.dict(
            os.environ, {"DNA_PAPER_FIGURE_OUTPUT": directory}
        ), patch("IPython.display.display"), redirect_stdout(StringIO()):
            for cell in self.notebook["cells"]:
                if "setup" in cell["metadata"].get("tags", []):
                    exec("".join(cell["source"]), namespace)
            figure = next(c for c in self.notebook["cells"] if "figure" in c["metadata"].get("tags", []))
            source = "".join(figure["source"]).replace("SAVE = True", "SAVE = False", 1)
            exec(source.replace("SHOW_RECORDED_COSTS = True", "SHOW_RECORDED_COSTS = False", 1), namespace)
            ax = namespace["FIGURES"]["F1_dna_trajectories"].axes[0]
            self.assertEqual(["Synthesis Cost Trajectory", "Sequencing Cost Trajectory"], ax.get_legend_handles_labels()[1])
            fig = namespace["FIGURES"]["F1_dna_trajectories"]
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
            self.assertLessEqual(ax.get_legend().get_window_extent(renderer).y1,
                                 ax.get_window_extent(renderer).y1)
            self.assertTrue(all(line.get_linestyle() == "-" for line in ax.lines))
            self.assertTrue(all(spine.get_visible() for spine in ax.spines.values()))
            self.assertFalse(any(line.get_visible() for line in ax.get_ygridlines()))
            self.assertIn(2026, ax.get_xticks())
            self.assertIn(2050, ax.get_xticks())
            gaps = namespace["TABLES"]["F1_dna_trajectories"]["gaps"]
            self.assertEqual(gaps.year.tolist(), [2026, 2050])
            for row, annotation in zip(gaps.itertuples(), ax.texts):
                import math
                self.assertAlmostEqual(row.orders_of_magnitude, math.log10(row.synthesis_to_sequencing_ratio))
                self.assertTrue(annotation.get_text().endswith(rf"$10^{{{round(row.orders_of_magnitude)}}}$"))
                self.assertIsNone(annotation.get_bbox_patch())
            stems = [ax.collections[i].get_segments()[0] for i in (0, 2)]
            for row, segment in zip(gaps.itertuples(), stems):
                self.assertEqual(segment[:, 0].tolist(), [row.year, row.year])
                self.assertAlmostEqual(segment[0, 1], row.sequencing_usd_per_mb)
                self.assertAlmostEqual(segment[1, 1], row.synthesis_usd_per_mb)
            guides = [collection.get_segments()[0] for collection in ax.collections[-2:]]
            for row, segment in zip(gaps.itertuples(), guides):
                self.assertEqual(segment[:, 0].tolist(), [row.year, row.year])
                self.assertEqual(segment[0, 1], ax.get_ylim()[0])
                self.assertAlmostEqual(segment[1, 1], row.sequencing_usd_per_mb)
            self.assertEqual(list(Path(directory).iterdir()), [])
            namespace["save_figure"]("F1_dna_trajectories")
            for extension in ("pdf", "svg", "png", "json"):
                self.assertTrue((Path(directory) / f"F1_dna_trajectories.{extension}").is_file())
            metadata = json.loads((Path(directory) / "F1_dna_trajectories.json").read_text())
            self.assertEqual(metadata["scenarios"]["baseline"]["start_year"], 2026)
            old_costs = namespace["TABLES"]["F1_dna_trajectories"]["costs"].copy()
            old_metrics = namespace["TABLES"]["F1_dna_trajectories"]["fit_metrics"].copy()
            self.assertTrue(namespace["TABLES"]["F1_dna_trajectories"]["observations"].empty)
            exec(source, namespace)
            tables = namespace["TABLES"]["F1_dna_trajectories"]
            observed = tables["observations"]
            namespace["pd"].testing.assert_frame_equal(old_metrics, tables["fit_metrics"])
            for row in tables["fit_metrics"].itertuples():
                points = tables["fit_observations"].query("process == @row.process")
                actual = points.cost_usd_per_logical_mb
                predicted = points.predicted_usd_per_logical_mb
                self.assertEqual(row.n_observations, len(points))
                self.assertAlmostEqual(row.r_squared,
                    1 - ((actual - predicted) ** 2).sum() / ((actual - actual.mean()) ** 2).sum())
                self.assertAlmostEqual(row.rmse_usd_per_logical_mb,
                    float(((actual - predicted) ** 2).mean() ** 0.5))
            self.assertIn("Goodness of fit", namespace["CAPTIONS"]["F1_dna_trajectories"])
            synth = observed[observed.process == "Synthesis"]
            seq = observed[observed.process == "Sequencing"]
            self.assertEqual(len(synth), len(namespace["synthesis_observed"]))
            self.assertEqual(len(seq), len(namespace["load_observed_sequencing_costs"]()))
            self.assertAlmostEqual(synth.iloc[0].cost_usd_per_logical_mb, 0.75 * 4e6)
            self.assertAlmostEqual(seq.iloc[0].cost_usd_per_logical_mb, seq.iloc[0].recorded_cost * 4)
            namespace["pd"].testing.assert_frame_equal(
                old_costs, tables["costs"].query("year >= 2000").reset_index(drop=True))
            ax = namespace["FIGURES"]["F1_dna_trajectories"].axes[0]
            self.assertIn("Synthesis recorded costs", ax.get_legend_handles_labels()[1])
            self.assertIn("Sequencing recorded costs", ax.get_legend_handles_labels()[1])
            self.assertLess(ax.get_xlim()[0], 1995)
            plt.close("all")
