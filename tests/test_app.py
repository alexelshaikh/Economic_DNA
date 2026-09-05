import unittest
from pathlib import Path

from streamlit.delta_generator_singletons import get_dg_singleton_instance
from streamlit.elements.lib.form_utils import FormData
from streamlit.testing.v1 import AppTest

from economic_dna import PRESET_SCENARIOS, Scenario


class StreamlitAppTests(unittest.TestCase):
    APP_PATH = Path(__file__).resolve().parent.parent / "streamlit_app.py"

    @staticmethod
    def _submit_form(app: AppTest) -> None:
        next(button for button in app.button if button.label == "Calculate").click()

    @staticmethod
    def _click_button(app: AppTest, key: str) -> None:
        next(button for button in app.button if button.key == key).click()

    def test_paper_baseline_renders_without_exceptions(self):
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20).run()
        self.assertFalse(app.exception)
        self.assertEqual(app.title[0].value, "DNA Storage Cost Explorer")
        self.assertEqual(
            [tab.label for tab in app.tabs],
            [
                "Lifecycle",
                "Start-year outlook",
                "DNA unit costs",
                "Sensitivity",
                "Assumptions",
                "About",
            ],
        )
        self.assertTrue(
            any(
                "https://doi.org/10.48550/arXiv.2608.26342" in (markdown.value or "")
                for markdown in app.markdown
            )
        )
        self.assertTrue(
            any("Paper / About" in (markdown.value or "") for markdown in app.markdown)
        )
        self.assertTrue(
            any("mailto:alex@el-shaikh.com" in (markdown.value or "") for markdown in app.markdown)
        )
        self.assertEqual(len(app.metric), 4)
        self.assertEqual(len(app.get("plotly_chart")), 7)
        download_buttons = app.get("download_button")
        self.assertEqual(len(download_buttons), 7)
        self.assertEqual(
            {button.key for button in download_buttons},
            {
                f"download_{graph}_csv"
                for graph in (
                    "lifecycle",
                    "breakdown",
                    "projection",
                    "dna_synthesis",
                    "dna_sequencing",
                    "breakeven",
                    "sensitivity",
                )
            },
        )
        image_export_buttons = {
            button.key
            for button in app.button
            if button.key and button.key.startswith("download_")
        }
        self.assertEqual(
            image_export_buttons,
            {
                f"download_{graph}_{file_format}"
                for graph in (
                    "lifecycle",
                    "breakdown",
                    "projection",
                    "dna_synthesis",
                    "dna_sequencing",
                    "breakeven",
                    "sensitivity",
                )
                for file_format in ("png", "svg")
            },
        )

    def test_rerun_starts_outside_synthetic_main_form(self):
        get_dg_singleton_instance().main_dg._form_data = FormData("scenario_form")
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20).run()
        self.assertFalse(app.exception)

    def test_form_submission_recalculates_archive(self):
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20).run()
        app.number_input(key="archive_value").set_value(2.0)
        self._submit_form(app)
        app.run()
        self.assertFalse(app.exception)
        archive_metric = next(metric for metric in app.metric if metric.label == "Archive")
        self.assertEqual(archive_metric.value, "2 TB")

    def test_sidebar_reset_updates_widget_state_before_calculate(self):
        baseline = Scenario()
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20).run()
        self.assertFalse(app.exception)

        app.number_input(key="archive_value").set_value(2.0)
        self._submit_form(app)
        app.run()
        archive_metric = next(metric for metric in app.metric if metric.label == "Archive")
        self.assertEqual(archive_metric.value, "2 TB")

        self._click_button(app, "sidebar_reset")
        app.run()
        self.assertEqual(app.number_input(key="archive_value").value, baseline.archive_size_tb)
        archive_metric = next(metric for metric in app.metric if metric.label == "Archive")
        self.assertEqual(archive_metric.value, "2 TB")

        self._submit_form(app)
        app.run()
        archive_metric = next(metric for metric in app.metric if metric.label == "Archive")
        self.assertEqual(archive_metric.value, "1 TB")

    def test_global_reset_restores_sidebar_and_model_inputs(self):
        baseline = Scenario()
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20).run()
        self.assertFalse(app.exception)

        app.number_input(key="archive_value").set_value(2.0)
        app.number_input(key="dna_synthesis_cost").set_value(123.0)
        self._click_button(app, "global_reset")
        app.run()

        self.assertEqual(app.number_input(key="archive_value").value, baseline.archive_size_tb)
        self.assertEqual(
            app.number_input(key="dna_synthesis_cost").value,
            baseline.dna_synthesis_cost_per_mb,
        )

    def test_custom_storage_can_be_selected_and_named(self):
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20).run()
        app.checkbox(key="tech_custom").check()
        app.text_input(key="custom_name").set_value("Test service")
        self._submit_form(app)
        app.run()
        self.assertFalse(app.exception)
        lowest_metric = next(metric for metric in app.metric if metric.label.startswith("Lowest:"))
        self.assertEqual(lowest_metric.label, "Lowest: Test service")

    def test_builtin_models_have_separate_editable_assumption_panels(self):
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20).run()
        self.assertFalse(app.exception)
        radio = app.radio(key="cost_model_radio")
        self.assertEqual(
            radio.options, ["✕", "DNA", "Amazon", "Azure", "Tape", "Custom"]
        )
        # The panel widgets are always mounted (hidden by CSS when closed).
        self.assertEqual(app.number_input(key="amazon_put_per_1000").value, 0.05)
        self.assertEqual(app.number_input(key="azure_storage_per_tb_month").value, 1.953125)
        self.assertEqual(app.number_input(key="tape_media_per_tb").value, 6.39)
        self.assertEqual(
            app.number_input(key="tape_media_per_tb").label,
            "Tape cartridges (USD/TB per write)",
        )
        self.assertEqual(
            app.number_input(key="tape_hardware_per_tb").label,
            "Tape library/drives (USD/TB amortized)",
        )
        self.assertEqual(
            app.number_input(key="tape_media_decline").label,
            "Tape cartridge decline (%)",
        )
        self.assertEqual(
            app.number_input(key="tape_hardware_decline").label,
            "Tape library/drives decline (%)",
        )
        self.assertTrue(
            any("avoid double counting" in (caption.value or "") for caption in app.caption)
        )

    def test_model_reset_button_restores_model_inputs_only(self):
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20).run()
        self.assertFalse(app.exception)
        baseline = Scenario()

        app.number_input(key="archive_value").set_value(2.0)
        app.number_input(key="dna_synthesis_cost").set_value(123.0)
        self._click_button(app, "reset_dna")
        app.run()

        self.assertEqual(app.number_input(key="archive_value").value, 2.0)
        self.assertEqual(
            app.number_input(key="dna_synthesis_cost").value,
            baseline.dna_synthesis_cost_per_mb,
        )

    def test_sensitivity_tab_shows_breakeven_table_and_tornado_chart(self):
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20).run()
        self.assertFalse(app.exception)
        chart_keys = {chart.key for chart in app.get("plotly_chart")}
        self.assertIn("chart_breakeven", chart_keys)
        self.assertIn("chart_sensitivity", chart_keys)
        table_html = " ".join(markdown.value or "" for markdown in app.markdown)
        self.assertIn("Break-even synthesis cost", table_html)
        self.assertIn("Today vs. break-even", table_html)
        # The paper-baseline default is not reachable by lowering synthesis
        # cost alone (DNA's sequencing/retrieval cost already exceeds every
        # built-in alternative) -- verify that lands in the table, not raw
        # scientific notation or a blank cell.
        self.assertIn("Not reachable", table_html)

    def test_sensitivity_tab_asks_for_dna_when_it_is_not_selected(self):
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20).run()
        app.checkbox(key="tech_dna").uncheck()
        self._submit_form(app)
        app.run()
        self.assertFalse(app.exception)
        chart_keys = {chart.key for chart in app.get("plotly_chart")}
        self.assertNotIn("chart_breakeven", chart_keys)
        self.assertNotIn("chart_sensitivity", chart_keys)
        self.assertTrue(
            any("break-even price" in (info.value or "") for info in app.info)
        )

    def test_uncertainty_band_checkbox_is_only_offered_when_dna_is_selected(self):
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20).run()
        self.assertIn("show_uncertainty_band", {c.key for c in app.checkbox})

        app.checkbox(key="tech_dna").uncheck()
        self._submit_form(app)
        app.run()
        self.assertFalse(app.exception)
        self.assertNotIn("show_uncertainty_band", {c.key for c in app.checkbox})

    def test_preset_button_populates_fields_without_committing_until_calculate(self):
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20).run()
        preset = PRESET_SCENARIOS["1 PB genomics cold archive"]

        self._click_button(app, "preset_1")
        app.run()
        self.assertFalse(app.exception)
        self.assertEqual(app.number_input(key="archive_value").value, 1.0)
        self.assertEqual(app.selectbox(key="archive_unit").value, "PB")
        self.assertEqual(app.number_input(key="asset_value").value, preset.average_asset_size_mb)
        self.assertEqual(app.number_input(key="retrieval").value, preset.annual_retrieval_percent)
        self.assertEqual(app.number_input(key="horizon").value, preset.horizon_years)

        # Fields are populated, but the previous (paper baseline) results are
        # still on screen until Calculate is pressed -- same as Reset.
        archive_metric = next(metric for metric in app.metric if metric.label == "Archive")
        self.assertEqual(archive_metric.value, "1 TB")

        self._submit_form(app)
        app.run()
        self.assertFalse(app.exception)
        archive_metric = next(metric for metric in app.metric if metric.label == "Archive")
        self.assertEqual(archive_metric.value, "1K TB")

    def test_every_preset_button_is_wired_to_a_known_preset(self):
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20).run()
        preset_buttons = [button for button in app.button if button.key and button.key.startswith("preset_")]
        self.assertEqual(len(preset_buttons), len(PRESET_SCENARIOS))
        self.assertEqual({button.label for button in preset_buttons}, set(PRESET_SCENARIOS))

    def test_order_of_magnitude_steppers_scale_the_synthesis_cost_field(self):
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20).run()
        before = app.number_input(key="dna_synthesis_cost").value

        self._click_button(app, "dna_synthesis_cost_mul10")
        app.run()
        self.assertFalse(app.exception)
        self.assertAlmostEqual(app.number_input(key="dna_synthesis_cost").value, before * 10)

        self._click_button(app, "dna_synthesis_cost_div10")
        app.run()
        self.assertFalse(app.exception)
        self.assertAlmostEqual(app.number_input(key="dna_synthesis_cost").value, before)

    def test_dna_cost_fields_use_a_significant_figure_format_that_never_displays_as_zero(self):
        # A fixed-decimal format like "%.6f" silently displays any value
        # below its decimal precision (e.g. 1e-7) as "0.000000" -- the
        # underlying value is untouched (Streamlit's format is display-only),
        # but it looks exactly like data loss to a user. "%.Ng" scales with
        # magnitude instead, so very small values stay visibly non-zero.
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20).run()
        for key in ("dna_synthesis_cost", "dna_sequencing_cost"):
            with self.subTest(key=key):
                widget = app.number_input(key=key)
                self.assertRegex(widget.proto.format, r"^%\.\d+g$")

    def test_setting_a_very_small_synthesis_cost_is_preserved_through_calculate(self):
        # Regression check for the format-only display issue: the committed
        # value must survive Calculate exactly, not just the display text.
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20).run()
        app.number_input(key="dna_synthesis_cost").set_value(1e-7)
        self._submit_form(app)
        app.run()
        self.assertFalse(app.exception)
        self.assertEqual(app.number_input(key="dna_synthesis_cost").value, 1e-7)


    def test_cost_tabs_are_a_form_radio(self):
        # The tabs are a Streamlit radio INSIDE the form: the frontend manages
        # the checked state instantly, so opening/closing never reruns the
        # script (the buttons and charts stay untouched), and the CSS reads
        # the checked input's value to open the matching panel.
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20).run()
        self.assertFalse(app.exception)
        radio = app.radio(key="cost_model_radio")
        self.assertEqual(radio.options, ["✕", "DNA", "Amazon", "Azure", "Tape", "Custom"])
        self.assertEqual(radio.value, "✕")


if __name__ == "__main__":
    unittest.main()
