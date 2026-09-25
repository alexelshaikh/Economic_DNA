"""Interaction checks and screenshots for a running local app."""
from __future__ import annotations

import argparse
import re
from pathlib import Path

from PIL import Image, ImageStat
from playwright.sync_api import expect, sync_playwright


def check(url: str, output: Path, channel: str | None) -> None:
    output.mkdir(parents=True, exist_ok=True)
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(channel=channel, headless=True)
        context = browser.new_context(viewport={"width": 1440, "height": 1000}, color_scheme="light")
        page = context.new_page()
        errors = []
        sent_frames = []
        page.on("pageerror", lambda error: errors.append(str(error)))
        page.on("websocket", lambda ws: ws.on("framesent", lambda payload: sent_frames.append(payload)))
        page.goto(url)
        expect(page.locator(".js-plotly-plot")).to_have_count(2)

        def button(key):
            return page.locator(f".st-key-{key} button:visible")

        def model(index):
            return page.locator(f'.st-key-cost-rail [data-testid="stRadioOption"]:has(input[value="{index}"])')

        def tab(label, chart):
            if model(0).is_visible():
                model(0).click()
            page.get_by_role("tab", name=label, exact=True).click()
            expect(page.locator(f".st-key-chart_{chart} .js-plotly-plot")).to_be_visible()

        def shot(name, current=page):
            current.screenshot(path=str(output / f"{name}.png"))

        def pending(changed):
            expect(button("calculate_header")).to_contain_text("Update charts" if changed else "Calculate")
            expect(page.locator(".pending-notice")).to_have_text(
                "Recalculate to update charts." if changed else "Charts up to date"
            )

        def check_initial_model_costs(current):
            for key, value in {
                "dna_synthesis_cost": "16573.83854",
                "dna_sequencing_cost": "0.1347573487",
                "amazon_storage_per_tb_month": "0.966797",
                "azure_storage_per_tb_month": "1.953125",
                "tape_media_per_tb": "6.39",
                "tape_hardware_per_tb": "6.86",
            }.items():
                expect(current.locator(f".st-key-{key} input")).to_have_value(re.compile(re.escape(value) + r"0*$"))

        def check_navigation(current):
            tabs = current.get_by_role("tab")
            expect(tabs).to_have_count(6)
            expect(current.locator('[role="tab"][aria-selected="true"]')).to_have_css("background-color", "rgb(8, 127, 140)")
            assert tabs.evaluate_all("""elements => elements.every(tab => {
                const rect = tab.getBoundingClientRect();
                const list = tab.closest('[role="tablist"]').getBoundingClientRect();
                const text = tab.querySelector('p').getBoundingClientRect();
                return rect.height >= 48 && rect.left >= list.left - 1 && rect.right <= list.right + 1
                    && rect.bottom <= list.bottom + 1 && text.left >= rect.left && text.right <= rect.right;
            })"""), "Analysis tabs are clipped or too small"
            for control in current.locator('.st-key-cost-rail label:visible').all():
                assert control.bounding_box()["height"] >= 44
                expect(control).to_have_css("border-top-width", "1px")

        page.wait_for_function("window.__dnaFormControls !== undefined")
        check_initial_model_costs(page)
        pending(False)
        check_navigation(page)
        initial_frames = len(sent_frames)
        chart = page.locator(".st-key-chart_lifecycle .js-plotly-plot")
        chart.evaluate("e => e.dataset.inputActionCheck = 'original'")
        archive = page.get_by_label("Archive size", exact=True)
        archive.fill("2")
        pending(True)  # Visible before blur or form submission.
        archive.fill("1")
        pending(False)
        archive.press("Tab")
        for key in ("tech_dna", "log_scale"):
            label = page.locator(f".st-key-{key} label:has(input)")
            label.click()
            pending(True)
            label.click()
            pending(False)
        page.get_by_text("Example scenarios", exact=True).click()
        for index, retrieval in ((1, "0.5"), (2, "20"), (3, "2"), (4, "0.001")):
            button(f"preset_{index}").click()
            expect(page.get_by_label("Annual retrieval (%)", exact=True)).to_have_value(retrieval)
            pending(True)
        button("global_reset").click()
        pending(False)
        model(1).click()
        expect(model(1).locator("p")).to_have_css("color", "rgb(255, 255, 255)")
        synthesis = page.get_by_label("Synthesis cost (USD/MB)", exact=True)
        for field in ("dna_synthesis_cost", "dna_sequencing_cost"):
            for suffix in ("div10", "mul10"):
                button(f"{field}_{suffix}").hover()
                page.wait_for_timeout(600)
                expect(page.get_by_role("tooltip")).to_have_count(0)
        synthesis.fill("123")
        button("dna_synthesis_cost_div10").click()
        expect(synthesis).to_have_value("12.3")
        button("dna_synthesis_cost_mul10").evaluate("b => { b.click(); b.click(); b.click(); }")
        expect(synthesis).to_have_value("12300")
        button("reset_dna").click()
        pending(False)
        model(5).click()
        page.get_by_label("Display name", exact=True).fill("My archive")
        pending(True)
        button("reset_custom").click()
        pending(False)
        button("sidebar_reset").click()
        pending(False)
        page.get_by_text("Example scenarios", exact=True).click()
        expect(chart).to_have_attribute("data-input-action-check", "original")
        page.wait_for_timeout(200)
        assert len(sent_frames) == initial_frames, "Input-only actions sent a server request"
        shot("desktop")
        header = page.locator(".st-key-calculate_header")
        header.evaluate("e => e.dataset.stabilityCheck = 'original'")
        page.get_by_label("Archive size", exact=True).fill("2")
        page.get_by_label("Archive size", exact=True).press("Tab")
        model(1).click()
        page.get_by_label("Synthesis cost (USD/MB)", exact=True).fill("123")
        model(2).click()
        model(1).click()
        expect(page.get_by_label("Synthesis cost (USD/MB)", exact=True)).to_have_value("123")
        expect(header).to_have_attribute("data-stability-check", "original")
        shot("cost-panel")
        tab("Sensitivity", "breakeven")
        expect(header).to_have_attribute("data-stability-check", "original")
        expect(page.get_by_label("Archive size", exact=True)).to_have_value(re.compile(r"2(?:\.0+)?"))
        button("calculate_header").click()
        expect(page.locator("[data-testid=stMetricValue]").first).to_have_text("2 TB")
        pending(False)
        expect(page.get_by_role("tab", name="Sensitivity", exact=True)).to_have_attribute("aria-selected", "true")
        frames_after_calculate = len(sent_frames)
        model(1).click()
        button("dna_synthesis_cost_div10").click()
        expect(synthesis).to_have_value("12.3")
        pending(True)
        button("dna_synthesis_cost_mul10").click()
        pending(False)
        assert len(sent_frames) == frames_after_calculate
        page.get_by_label("Archive size", exact=True).fill("3")
        button("theme_toggle").click()
        expect(page.locator(".theme-dark")).to_be_attached()
        expect(page.get_by_role("tab", name="Sensitivity", exact=True)).to_have_attribute("aria-selected", "true")
        expect(page.get_by_label("Archive size", exact=True)).to_have_value(re.compile(r"3(?:\.0+)?"))
        pending(True)
        expander = page.locator('[data-testid="stExpander"] summary')
        expect(expander).to_have_css("background-color", "rgb(27, 30, 35)")
        expect(expander.locator("p")).to_have_css("color", "rgb(233, 237, 242)")
        check_navigation(page)
        page.mouse.move(400, 80)
        shot("dark")
        button("theme_toggle").click()
        expect(page.locator(".theme-light")).to_be_attached()

        page.get_by_text("Example scenarios", exact=True).click()
        button("preset_4").click()
        expect(page.get_by_label("Annual retrieval (%)", exact=True)).to_have_value("0.001")
        model(1).click()
        expect(page.get_by_label("Synthesis cost (USD/MB)", exact=True)).to_have_value("123")
        button("calculate_header").click()
        expect(page.locator("[data-testid=stMetricValue]").first).to_have_text("1K TB")
        pending(False)
        tab("Sensitivity", "breakeven")
        expect(page.locator(".contract-table").first).not_to_contain_text("Not reachable")
        expect(page.locator(".st-key-chart_breakeven .barlayer .point")).to_have_count(3)
        shot("preservation-sensitivity")
        for file_format in ("csv", "png", "svg"):
            with page.expect_download(timeout=15000) as download:
                button(f"download_breakeven_{file_format}").click()
            target = output / download.value.suggested_filename
            download.value.save_as(target)
            assert target.stat().st_size > 100, target
        shared = context.new_page()
        shared.goto(page.url)
        expect(shared.locator("[data-testid=stMetricValue]").first).to_have_text("1K TB")
        expect(shared.get_by_label("Annual retrieval (%)", exact=True)).to_have_value("0.001")
        shared.close()
        tab("DNA unit costs", "dna_synthesis")
        page.get_by_text("Show historical context", exact=True).click()
        expect(page.get_by_label("Show historical context", exact=True)).not_to_be_checked()
        expect(page.get_by_role("tab", name="DNA unit costs", exact=True)).to_have_attribute("aria-selected", "true")
        tab("Sensitivity", "breakeven")
        tab("DNA unit costs", "dna_synthesis")
        expect(page.get_by_label("Show historical context", exact=True)).not_to_be_checked()
        shot("dna-unit-costs")

        page.get_by_label("Average asset size", exact=True).fill("1000000000")
        button("calculate_header").click()
        expect(page.get_by_text(re.compile("no larger than the archive"))).to_be_visible()
        pending(True)
        expect(page.locator("[data-testid=stMetricValue]").first).to_have_text("1K TB")
        page.get_by_label("Average asset size", exact=True).fill("1")
        pending(False)
        button("global_reset").click()
        pending(True)
        expect(page.locator("[data-testid=stMetricValue]").first).to_have_text("1K TB")
        button("calculate_header").click()
        pending(False)
        expect(page.locator("[data-testid=stMetricValue]").first).to_have_text("1 TB")

        dark = browser.new_page(viewport={"width": 1440, "height": 1000}, color_scheme="dark")
        dark.on("pageerror", lambda error: errors.append(str(error)))
        dark.goto(url)
        expect(dark.locator(".theme-dark")).to_be_attached()
        expect(dark.locator(".js-plotly-plot")).to_have_count(2)
        dark.wait_for_function("window.__dnaFormControls !== undefined")
        check_initial_model_costs(dark)
        expect(dark.locator(".pending-notice")).to_have_text("Charts up to date")
        dark.locator(".st-key-calculate_header button:visible").click()
        expect(dark).to_have_url(re.compile(r"dna_synthesis_cost_per_mb=16573"))
        check_initial_model_costs(dark)
        shot("first-visit-dark", dark)
        dark.close()

        for width in (768, 390, 320):
            mobile = browser.new_page(viewport={"width": width, "height": 900}, is_mobile=True, has_touch=True)
            mobile.on("pageerror", lambda error: errors.append(str(error)))
            mobile.goto(url)
            expect(mobile.locator(".js-plotly-plot")).to_have_count(2)
            check_initial_model_costs(mobile)
            mobile.wait_for_timeout(300)
            check_navigation(mobile)
            shot(f"viewport-{width}", mobile)
            assert mobile.locator('[data-testid="stMain"]').evaluate("e => e.scrollWidth <= e.clientWidth + 1")
            if width == 390:
                origin = mobile.evaluate("performance.timeOrigin")
                mobile.locator('[data-testid="stExpandSidebarButton"]').click()
                expect(mobile.locator('[data-testid="stSidebar"]')).to_have_attribute("aria-expanded", "true")
                mobile.wait_for_timeout(300)
                shot("mobile-inputs", mobile)
                mobile.get_by_label("Archive size", exact=True).fill("5")
                expect(mobile.locator('.st-key-calculate_scenario button:visible')).to_contain_text("Update charts")
                mobile.locator('.st-key-calculate_scenario button:visible').click()
                expect(mobile.locator('[data-testid="stSidebar"]')).to_have_attribute("aria-expanded", "false")
                expect(mobile.locator('[data-testid="stSidebar"]')).to_have_css("opacity", "0")
                expect(mobile.locator("[data-testid=stMetricValue]").first).to_have_text("5 TB")
                expect(mobile.locator('.st-key-calculate_header button:visible')).to_contain_text("Calculate")
                assert mobile.evaluate("performance.timeOrigin") == origin, "Calculate reloaded the page"
                mobile.locator('.st-key-chart_lifecycle').scroll_into_view_if_needed()
                mobile.wait_for_timeout(300)
                graph = mobile.locator('.st-key-chart_lifecycle')
                assert graph.evaluate("e => { const legend=e.querySelector('.legend').getBoundingClientRect(); const plot=e.querySelector('.nsewdrag').getBoundingClientRect(); return legend.top >= plot.bottom; }")
                shot("mobile-chart", mobile)
                mobile.locator('.st-key-cost-rail [data-testid="stRadioOption"]:has(input[value="1"])').click()
                expect(mobile.get_by_label("Synthesis cost (USD/MB)", exact=True)).to_be_visible()
                mobile.locator('.st-key-dna_synthesis_cost_div10 button:visible').click()
                expect(mobile.locator('.st-key-calculate_panel button:visible')).to_contain_text("Update charts")
                shot("mobile-cost-panel", mobile)
            mobile.close()
        assert not errors, errors
        browser.close()
    for path in output.glob("*.png"):
        with Image.open(path) as image:
            assert sum(ImageStat.Stat(image.convert("RGB")).var) > 100, f"{path} appears blank"
    print(f"Browser checks passed. Screenshots and downloads: {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--url", default="http://127.0.0.1:8501")
    parser.add_argument("--output", type=Path, default=Path(".tmp/browser-check"))
    parser.add_argument("--channel", default=None, help="Installed browser, for example msedge or chrome")
    args = parser.parse_args()
    check(args.url, args.output, args.channel)
