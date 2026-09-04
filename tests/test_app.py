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

    def test_sidebar_and_global_resets_are_client_side(self):
        baseline = Scenario()
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20).run()

        self.assertFalse(app.exception)
        sidebar_reset = next(
            md
            for md in app.markdown
            if 'class="sidebar-reset-btn"' in (md.value or "")
        )
        self.assertIn("Reset sidebar values to paper baseline", sidebar_reset.value)
        self.assertIn(f"&quot;archive_value&quot;: {baseline.archive_size_tb}", sidebar_reset.value)
        self.assertIn(f"&quot;start_year_widget&quot;: {baseline.start_year}", sidebar_reset.value)
        self.assertIn("&quot;tech_tape&quot;: true", sidebar_reset.value)
        self.assertNotIn("&quot;dna_synthesis_cost&quot;", sidebar_reset.value)

        global_reset = next(
            md
            for md in app.markdown
            if 'class="global-reset-btn"' in (md.value or "")
        )
        self.assertIn("Reset all inputs to paper baseline", global_reset.value)
        self.assertIn(f"&quot;archive_value&quot;: {baseline.archive_size_tb}", global_reset.value)
        self.assertIn("&quot;dna_synthesis_cost&quot;", global_reset.value)

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

    def test_model_reset_buttons_carry_baseline_values(self):
        # The per-model reset is a pure-HTML button whose data-values holds
        # the paper-baseline JSON — the page JS applies it client-side, so no
        # form submit or rerun happens when it is clicked.
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20).run()
        self.assertFalse(app.exception)
        baseline = Scenario()
        dna_btn = next(
            md
            for md in app.markdown
            if 'class="model-reset-btn" data-model="dna"' in (md.value or "")
        )
        self.assertIn(f'"dna_cost_base_year": {baseline.dna_cost_base_year}', dna_btn.value)
        self.assertIn("dna_synthesis_cost", dna_btn.value)
        amazon_btn = next(
            md
            for md in app.markdown
            if 'class="model-reset-btn" data-model="amazon"' in (md.value or "")
        )
        self.assertIn("amazon_storage_per_tb_month", amazon_btn.value)

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
