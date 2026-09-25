import json
import unittest
from pathlib import Path

from streamlit.delta_generator_singletons import get_dg_singleton_instance
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

    @staticmethod
    def _open_tab(app: AppTest, label: str) -> None:
        # All panels now exist up front; actual navigation is browser-only.
        assert label in [tab.label for tab in app.tabs]
        app.run()

    @staticmethod
    def _form_config(app: AppTest) -> dict:
        body = next(element.proto.body for element in app.get("html") if "form-controls-marker" in element.proto.body)
        return json.JSONDecoder().raw_decode(body.split("const config = ", 1)[1])[0]

    def test_pending_state_accounts_for_display_precision(self):
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20).run()
        config = self._form_config(app)
        for widget in app.number_input:
            if widget.key not in config["displayed"]:
                continue
            with self.subTest(key=widget.key):
                self.assertEqual(config["displayed"][widget.key], float(widget.proto.format % widget.value))
        self.assertEqual(config["committed"], app.session_state["committed_widgets"])

    def test_model_widgets_start_with_real_defaults_in_the_browser_payload(self):
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20).run()
        defaults = {widget.key: app.session_state["widget_defaults"][widget.key] for widget in app.number_input}
        for key, value in defaults.items():
            with self.subTest(key=key):
                widget = app.number_input(key=key)
                self.assertEqual(widget.value, value)
                self.assertEqual(widget.proto.default, value)
                self.assertFalse(widget.proto.set_value)
        self.assertEqual(defaults["dna_synthesis_cost"], Scenario().dna_synthesis_cost_per_mb)
        self.assertEqual(defaults["amazon_storage_per_tb_month"], Scenario().amazon_storage_usd_per_mb_month * 1_000_000)

    def test_all_scenario_widgets_declare_their_initial_values(self):
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20).run()
        expected = app.session_state["committed_widgets"]
        seen = set()
        for kind in ("number_input", "radio", "checkbox", "toggle", "text_input"):
            for widget in app.get(kind):
                if widget.key not in expected:
                    continue
                with self.subTest(key=widget.key):
                    value = expected[widget.key]
                    declared = widget.proto.default
                    if kind == "radio":
                        declared = widget.options[declared]
                    self.assertEqual(declared, value)
                    self.assertEqual(widget.value, value)
                    self.assertFalse(widget.proto.set_value)
                    seen.add(widget.key)
        self.assertEqual(seen, set(expected))

    def test_initial_model_defaults_do_not_overwrite_user_edits(self):
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20).run()
        baseline = dict(app.session_state["widget_defaults"])
        app.number_input(key="dna_synthesis_cost").set_value(0.0)
        app.number_input(key="amazon_storage_per_tb_month").set_value(12.5)
        self._submit_form(app)
        app.run()
        self._open_tab(app, "Sensitivity")
        self.assertFalse(app.exception)
        self.assertEqual(app.number_input(key="dna_synthesis_cost").value, 0.0)
        self.assertEqual(app.number_input(key="amazon_storage_per_tb_month").value, 12.5)
        self.assertEqual(app.session_state["widget_defaults"], baseline)

    def test_automatic_theme_sync_preserves_initial_and_pending_inputs(self):
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20).run()
        before = dict(app.session_state["committed_widgets"])
        app.number_input(key="archive_value").set_value(2.0)
        app.radio(key="archive_unit").set_value("PB")
        self._click_button(app, "theme_auto_dark")
        app.run()
        self.assertFalse(app.exception)
        self.assertEqual(app.number_input(key="archive_value").value, 2.0)
        self.assertEqual(app.radio(key="archive_unit").value, "PB")
        self.assertEqual(app.number_input(key="horizon").value, 100)
        self.assertEqual(app.radio(key="asset_unit").value, "GB")
        self.assertTrue(app.checkbox(key="tech_dna").value)
        self.assertEqual(app.session_state["committed_widgets"], before)
        self.assertEqual(app.query_params["theme"], ["dark"])

    def test_small_workload_values_are_not_displayed_as_zero(self):
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20)
        app.query_params["archive_size_tb"] = "0.001"
        app.query_params["average_asset_size_mb"] = "0.001"
        app.run()
        self.assertFalse(app.exception)
        for key in ("archive_value", "asset_value"):
            widget = app.number_input(key=key)
            self.assertEqual(widget.proto.format % widget.value, "0.001")
        self._submit_form(app)
        app.run()
        self.assertFalse(app.exception)
        self.assertEqual(app.session_state["committed_widgets"]["archive_value"], 0.001)

    def test_shared_workload_defaults_survive_reload_and_theme_changes(self):
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20)
        app.query_params.update({
            "archive_size_tb": "3000", "average_asset_size_mb": "500",
            "horizon_years": "50", "annual_retrieval_percent": "0.25",
            "technologies": "Tape On-premise,Custom storage", "log_scale": "False",
            "custom_storage_name": "Shared archive", "theme": "dark",
        })
        app.run()
        self.assertFalse(app.exception)
        self.assertEqual(app.number_input(key="archive_value").proto.default, 3.0)
        self.assertEqual(app.radio(key="archive_unit").proto.default, 1)
        self.assertFalse(app.checkbox(key="tech_dna").proto.default)
        self.assertTrue(app.checkbox(key="tech_custom").proto.default)
        self.assertFalse(app.toggle(key="log_scale").proto.default)
        self.assertEqual(app.text_input(key="custom_name").proto.default, "Shared archive")
        before = dict(app.session_state["committed_widgets"])
        self._click_button(app, "theme_toggle")
        app.run()
        self._submit_form(app)
        app.run()
        self.assertEqual(app.session_state["committed_widgets"], before)
        reloaded = AppTest.from_file(str(self.APP_PATH), default_timeout=20)
        reloaded.query_params.update(app.query_params)
        reloaded.run()
        self.assertFalse(reloaded.exception)
        self.assertEqual(reloaded.session_state["committed_widgets"], before)

    def test_shared_model_prices_are_initial_widget_defaults(self):
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20)
        app.query_params["dna_synthesis_cost_per_mb"] = "0"
        app.query_params["tape_media_usd_per_tb"] = "12.5"
        app.run()
        self.assertFalse(app.exception)
        self.assertEqual(app.number_input(key="dna_synthesis_cost").proto.default, 0.0)
        self.assertEqual(app.number_input(key="tape_media_per_tb").proto.default, 12.5)

    def test_pending_reference_changes_only_after_calculate(self):
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20).run()
        before = self._form_config(app)["committed"]
        self._click_button(app, "preset_4")
        app.run()
        self.assertEqual(self._form_config(app)["committed"], before)
        self._submit_form(app)
        app.run()
        self.assertEqual(self._form_config(app)["committed"]["archive_unit"], "PB")
        self.assertEqual(self._form_config(app)["committed"]["retrieval"], 0.001)

    def test_client_configuration_escapes_user_html(self):
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20).run()
        name = "</script><script>window.untrusted = true</script>"
        app.text_input(key="custom_name").set_value(name)
        self._submit_form(app)
        app.run()
        self.assertEqual(self._form_config(app)["committed"]["custom_name"], name)
        body = next(element.proto.body for element in app.get("html") if "form-controls-marker" in element.proto.body)
        self.assertNotIn(name, body)

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
        self.assertEqual(len(app.metric), 4)
        self.assertEqual(len(app.get("plotly_chart")), 8)
        download_buttons = app.get("download_button")
        self.assertEqual(len(download_buttons), 8)
        self.assertEqual(
            {button.key for button in download_buttons},
            {
                f"download_{graph}_csv"
                for graph in (
                    "lifecycle",
                    "breakdown",
                    "projection", "dna_synthesis", "dna_sequencing", "viability", "breakeven", "sensitivity",
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
                    "projection", "dna_synthesis", "dna_sequencing", "viability", "breakeven", "sensitivity",
                )
                for file_format in ("png", "svg")
            },
        )

    def test_shared_main_generator_never_enters_the_scenario_form(self):
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20).run()
        self.assertFalse(app.exception)
        self.assertIsNone(get_dg_singleton_instance().main_dg._form_data)
        for widget in app.number_input:
            self.assertEqual(widget.proto.form_id, "scenario_form")
        self.assertEqual(app.selectbox(key="viability_driver").proto.form_id, "")
        self.assertEqual(app.text_input(key="viability_min").proto.form_id, "viability_range_form")

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

    def test_outlook_passes_the_same_crossover_years_to_chart_and_table(self):
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20).run()
        self._open_tab(app, "Start-year outlook")
        self.assertFalse(app.exception)
        chart = next(chart for chart in app.get("plotly_chart") if chart.key == "chart_projection")
        layout = json.loads(chart.proto.spec)["layout"]
        self.assertEqual([shape["x0"] for shape in layout["shapes"]], [2332, 2321, 2149])
        table = " ".join(markdown.value or "" for markdown in app.markdown)
        for annotation in layout["annotations"]:
            self.assertIn(str(annotation["x"]), table)

    def test_sensitivity_tab_shows_breakeven_table_and_tornado_chart(self):
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20).run()
        self._open_tab(app, "Sensitivity")
        self.assertFalse(app.exception)
        chart_keys = {chart.key for chart in app.get("plotly_chart")}
        self.assertIn("chart_viability", chart_keys)
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
        self._open_tab(app, "Sensitivity")
        app.checkbox(key="tech_dna").uncheck()
        self._submit_form(app)
        app.run()
        self._open_tab(app, "Sensitivity")
        self.assertFalse(app.exception)
        chart_keys = {chart.key for chart in app.get("plotly_chart")}
        self.assertNotIn("chart_viability", chart_keys)
        self.assertNotIn("chart_breakeven", chart_keys)
        self.assertNotIn("chart_sensitivity", chart_keys)
        self.assertTrue(
            any("break-even price" in (info.value or "") for info in app.info)
        )

    def test_viability_options_persist_without_committing_pending_inputs(self):
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20).run()
        self._open_tab(app, "Sensitivity")
        committed = dict(app.session_state["committed_widgets"])
        app.number_input(key="archive_value").set_value(2.0)
        app.selectbox(key="viability_driver").select("discount_rate_percent")
        self._open_tab(app, "Sensitivity")
        app.selectbox(key="viability_comparison").select("Tape On-premise")
        self._open_tab(app, "Sensitivity")
        app.toggle(key="viability_focus").set_value(False)
        self._open_tab(app, "Sensitivity")
        self.assertFalse(app.exception)
        self.assertEqual(app.session_state["committed_widgets"], committed)
        self.assertEqual(app.number_input(key="archive_value").value, 2.0)
        chart = next(c for c in app.get("plotly_chart") if c.key == "chart_viability")
        layout = json.loads(chart.proto.spec)["layout"]
        self.assertIn("Present-value", layout["yaxis"]["title"]["text"])
        self._open_tab(app, "Lifecycle")
        self._open_tab(app, "Sensitivity")
        self.assertFalse(app.exception)
        self.assertEqual(app.selectbox(key="viability_driver").value, "discount_rate_percent")
        self.assertEqual(app.selectbox(key="viability_comparison").value, "Tape On-premise")
        self.assertFalse(app.toggle(key="viability_focus").value)

    def test_viability_range_apply_validation_reset_and_driver_change(self):
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20).run()
        committed = dict(app.session_state["committed_widgets"])
        defaults = (app.text_input(key="viability_min").value, app.text_input(key="viability_max").value)

        def axis_range():
            chart = next(c for c in app.get("plotly_chart") if c.key == "chart_viability")
            return json.loads(chart.proto.spec)["layout"]["xaxis"]["range"]

        app.number_input(key="archive_value").set_value(2)
        app.text_input(key="viability_min").set_value("-2")
        app.text_input(key="viability_max").set_value("-1")
        app.button(key="viability_apply").click().run()
        self.assertFalse(app.exception)
        self.assertEqual(axis_range(), [-2, -1])
        self.assertEqual(app.session_state["committed_widgets"], committed)
        app.text_input(key="viability_min").set_value("nan")
        app.button(key="viability_apply").click().run()
        self.assertTrue(app.error)
        self.assertEqual(axis_range(), [-2, -1])
        app.button(key="viability_reset").click().run()
        self.assertFalse(app.exception)
        self.assertFalse(app.error)
        self.assertEqual((app.text_input(key="viability_min").value, app.text_input(key="viability_max").value), defaults)
        app.selectbox(key="viability_driver").select("discount_rate_percent").run()
        self.assertIsNone(app.session_state["viability_applied_range"])
        self.assertGreaterEqual(float(app.text_input(key="viability_min").value), 0)

    def test_viability_comparison_falls_back_when_model_is_removed(self):
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20).run()
        self._open_tab(app, "Sensitivity")
        app.selectbox(key="viability_comparison").select("Tape On-premise")
        self._open_tab(app, "Sensitivity")
        app.checkbox(key="tech_tape").uncheck()
        self._submit_form(app)
        self._open_tab(app, "Sensitivity")
        self.assertFalse(app.exception)
        self.assertEqual(app.selectbox(key="viability_comparison").value, "Amazon Deep Archive")
        self.assertIn("chart_viability", {c.key for c in app.get("plotly_chart")})
        for key in ("tech_amazon", "tech_azure"):
            app.checkbox(key=key).uncheck()
        self._submit_form(app)
        self._open_tab(app, "Sensitivity")
        self.assertFalse(app.exception)
        self.assertNotIn("chart_viability", {c.key for c in app.get("plotly_chart")})

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
        app.number_input(key="dna_synthesis_cost").set_value(123.0)
        app.number_input(key="amazon_storage_per_tb_month").set_value(456.0)
        app.number_input(key="azure_retrieval_per_tb").set_value(789.0)
        app.number_input(key="tape_media_per_tb").set_value(12.0)

        self._click_button(app, "preset_1")
        app.run()
        self.assertFalse(app.exception)
        self.assertEqual(app.number_input(key="archive_value").value, 1.0)
        self.assertEqual(app.radio(key="archive_unit").value, "PB")
        self.assertEqual(app.number_input(key="asset_value").value, preset.average_asset_size_mb)
        self.assertEqual(app.number_input(key="retrieval").value, preset.annual_retrieval_percent)
        self.assertEqual(app.number_input(key="horizon").value, preset.horizon_years)
        self.assertEqual(app.number_input(key="dna_synthesis_cost").value, 123.0)
        self.assertEqual(app.number_input(key="amazon_storage_per_tb_month").value, 456.0)
        self.assertEqual(app.number_input(key="azure_retrieval_per_tb").value, 789.0)
        self.assertEqual(app.number_input(key="tape_media_per_tb").value, 12.0)

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
        visible_preset_names = [name for name in PRESET_SCENARIOS if name != "Paper baseline"]
        preset_buttons = [button for button in app.button if button.key and button.key.startswith("preset_")]
        self.assertEqual(len(preset_buttons), len(visible_preset_names))
        self.assertEqual({button.label.splitlines()[0] for button in preset_buttons}, set(visible_preset_names))
        self.assertTrue(all(" | " in button.label for button in preset_buttons))
        self.assertNotIn("Paper baseline", {button.label.splitlines()[0] for button in preset_buttons})

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

    def test_dna_cost_steppers_have_no_obscuring_tooltips(self):
        app = AppTest.from_file(str(self.APP_PATH), default_timeout=20).run()
        for field in ("dna_synthesis_cost", "dna_sequencing_cost"):
            for suffix in ("div10", "mul10"):
                with self.subTest(field=field, suffix=suffix):
                    self.assertFalse(app.button(key=f"{field}_{suffix}").proto.help)

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
