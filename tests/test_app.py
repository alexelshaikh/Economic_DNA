import unittest
from pathlib import Path

from streamlit.testing.v1 import AppTest

from economic_dna import Scenario


class StreamlitAppTests(unittest.TestCase):
    APP_PATH = Path(__file__).resolve().parent.parent / "streamlit_app.py"

    @staticmethod
    def _submit_form(app: AppTest) -> None:
        next(button for button in app.button if button.label == "Calculate").click()

    def test_paper_baseline_renders_without_exceptions(self):
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20).run()
        self.assertFalse(app.exception)
        self.assertEqual(app.title[0].value, "DNA Storage Cost Explorer")
        self.assertEqual(
            [tab.label for tab in app.tabs],
            ["Lifecycle", "Start-year outlook", "DNA unit costs", "Assumptions"],
        )
        self.assertEqual(len(app.metric), 4)
        self.assertEqual(len(app.get("plotly_chart")), 5)
        download_buttons = app.get("download_button")
        self.assertEqual(len(download_buttons), 15)
        self.assertEqual(
            {button.key for button in download_buttons},
            {
                f"download_{graph}_{file_format}"
                for graph in (
                    "lifecycle",
                    "breakdown",
                    "projection",
                    "dna_synthesis",
                    "dna_sequencing",
                )
                for file_format in ("csv", "png", "svg")
            },
        )

    def test_form_submission_recalculates_archive(self):
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20).run()
        app.number_input(key="archive_value").set_value(2.0)
        self._submit_form(app)
        app.run()
        self.assertFalse(app.exception)
        archive_metric = next(metric for metric in app.metric if metric.label == "Archive")
        self.assertEqual(archive_metric.value, "2 TB")

    def test_reset_to_paper_baseline_updates_form_widgets(self):
        baseline = Scenario()
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20).run()
        app.number_input(key="archive_value").set_value(2.0)
        app.number_input(key="start_year_widget").set_value(2035)
        app.checkbox(key="tech_tape").uncheck()
        self._submit_form(app)
        app.run()

        app.button(key="reset_baseline").click()
        app.run()

        self.assertFalse(app.exception)
        self.assertEqual(app.number_input(key="archive_value").value, baseline.archive_size_tb)
        self.assertEqual(app.selectbox(key="archive_unit").value, "TB")
        self.assertEqual(app.number_input(key="start_year_widget").value, baseline.start_year)
        self.assertTrue(app.checkbox(key="tech_tape").value)

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
        rail_labels = [
            button.label
            for button in app.button
            if button.key and button.key.startswith("cost_tab_")
        ]
        self.assertEqual(rail_labels, ["DNA", "Amazon", "Azure", "Tape", "Custom"])
        # The panel widgets are always mounted (hidden by CSS when closed).
        self.assertEqual(app.number_input(key="amazon_put_per_1000").value, 0.05)
        self.assertEqual(app.number_input(key="azure_storage_per_tb_month").value, 1.953125)
        self.assertEqual(app.number_input(key="tape_media_per_tb").value, 6.39)

    def test_cost_tab_toggles_open_model(self):
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20).run()
        self.assertIsNone(app.session_state["open_cost_model"])
        next(button for button in app.button if button.key == "cost_tab_amazon").click()
        app.run()
        self.assertFalse(app.exception)
        self.assertEqual(app.session_state["open_cost_model"], "amazon")
        # Clicking the same tab again closes the panel.
        next(button for button in app.button if button.key == "cost_tab_amazon").click()
        app.run()
        self.assertIsNone(app.session_state["open_cost_model"])


if __name__ == "__main__":
    unittest.main()
