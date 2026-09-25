"""Verify client-only tab navigation and signed sensitivity range controls."""
from __future__ import annotations

import argparse
from pathlib import Path

from playwright.sync_api import expect, sync_playwright


def check(url: str, output: Path, channel: str | None) -> None:
    output.mkdir(parents=True, exist_ok=True)
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(channel=channel, headless=True)
        # Overlapping sessions must never share synthetic form state.
        visitors = [browser.new_page() for _ in range(3)]
        for index, visitor in enumerate(visitors, start=1):
            visitor.goto(f"{url}?theme=light&archive_size_tb={index}")
        for index, visitor in enumerate(visitors, start=1):
            expect(visitor.locator('.js-plotly-plot')).to_have_count(8, timeout=30000)
            expect(visitor.locator('[data-testid="stException"]')).to_have_count(0)
            expect(visitor.locator('[data-testid="stMetricValue"]').first).to_have_text(f"{index} TB")
            visitor.close()
        for width, theme in ((1440, "light"), (1440, "dark"), (390, "light"), (320, "dark")):
            page = browser.new_page(viewport={"width": width, "height": 1000}, color_scheme=theme)
            sent, errors = [], []
            page.on("pageerror", lambda error: errors.append(str(error)))
            page.on("websocket", lambda ws: ws.on("framesent", lambda payload: sent.append(payload)))
            page.goto(url + "?theme=" + theme)
            expect(page.locator(".js-plotly-plot")).to_have_count(8, timeout=30000)
            expect(page.get_by_text("Contact", exact=True)).to_be_attached()
            page.wait_for_timeout(500)
            if width <= 640:
                page.get_by_role("button", name="Open inputs", exact=True).click()
            page.get_by_label("Archive size", exact=True).fill("2")
            if width <= 640:
                page.get_by_role("button", name="Close inputs", exact=True).click()
            tabs = page.get_by_role("tablist")
            tabs.evaluate("e => e.scrollIntoView({block: 'start'})")
            page.locator('[data-testid="stMain"]').evaluate("e => e.scrollBy(0, -80)")
            page.evaluate("""() => {
                window.navigationSamples = [];
                window.originalTabs = document.querySelector('[role="tablist"]');
                window.originalCharts = [...document.querySelectorAll('.js-plotly-plot')];
                const sample = () => {
                    const tabs = document.querySelector('[role="tablist"]');
                    const panel = document.querySelector('[role="tabpanel"]');
                    let opacity = 1;
                    for (let e = tabs; e; e = e.parentElement) opacity *= Number(getComputedStyle(e).opacity);
                    window.navigationSamples.push({same: tabs === window.originalTabs, opacity,
                        y: tabs.getBoundingClientRect().y, filled: Boolean(panel?.textContent.trim()),
                        retained: window.originalCharts.every(e => e.isConnected)});
                    window.navigationFrame = requestAnimationFrame(sample);
                }; sample();
            }""")
            before = len(sent)
            for _ in range(2):
                for label, charts in (("Start-year outlook", 1), ("DNA unit costs", 2), ("Sensitivity", 3), ("Assumptions", 0), ("About", 0), ("Lifecycle", 2)):
                    page.get_by_role("tab", name=label, exact=True).click()
                    expect(page.locator('[role="tabpanel"] .js-plotly-plot')).to_have_count(charts)
                    page.wait_for_timeout(80)
                    assert len(sent) == before, "Tab navigation sent a server request"
                    expect(page.locator(".pending-notice")).to_have_text("Recalculate to update charts.")
                    expect(page.locator("[data-testid=stMetricValue]").first).to_have_text("1 TB")
            page.get_by_role("tab", name="Lifecycle", exact=True).focus()
            page.keyboard.press("End")
            expect(page.get_by_role("tabpanel")).to_contain_text("About this explorer")
            page.keyboard.press("Home")
            expect(page.get_by_role("tabpanel")).to_contain_text("Lifecycle comparison")
            assert len(sent) == before
            page.evaluate("cancelAnimationFrame(window.navigationFrame)")
            samples = page.evaluate("window.navigationSamples")
            assert all(s["same"] and s["retained"] for s in samples), "Tabs or graphs were remounted"
            assert all(s["opacity"] == 1 and s["filled"] for s in samples), "A blank or faded frame was shown"
            assert max(s["y"] for s in samples) - min(s["y"] for s in samples) < 2, "Tab navigation moved the page"

            page.get_by_role("tab", name="Sensitivity", exact=True).click()
            chart = page.locator(".st-key-chart_viability .js-plotly-plot")
            expect(chart.locator(".annotation-text")).to_contain_text(["Theoretical parity"])
            minimum = page.locator('.st-key-viability_min input')
            maximum = page.locator('.st-key-viability_max input')
            defaults = (minimum.input_value(), maximum.input_value())
            expect(minimum).to_have_css("color", "rgb(233, 237, 242)" if theme == "dark" else "rgb(32, 38, 45)")
            expect(page.locator('.st-key-viability_min label p')).to_have_css("color", "rgb(233, 237, 242)" if theme == "dark" else "rgb(32, 38, 45)")
            expect(page.locator('.st-key-viability_min [data-testid="stTextInputRootElement"]')).to_have_css("background-color", "rgb(30, 33, 38)" if theme == "dark" else "rgb(255, 255, 255)")
            lifecycle = page.locator('.st-key-chart_lifecycle .js-plotly-plot')
            lifecycle.evaluate("e => {window.savedLifecycle = e; window.savedLifecycleData = JSON.stringify(e.data);}")
            initial_range = chart.evaluate("e => e.layout.xaxis.range")
            before = len(sent)
            minimum.fill("-2")
            maximum.fill("-1")
            maximum.blur()
            page.wait_for_timeout(100)
            assert len(sent) == before, "Typing axis limits should not rerun the server"
            assert chart.evaluate("e => e.layout.xaxis.range") == initial_range
            page.locator('.st-key-viability_apply button').click()
            page.wait_for_function("JSON.stringify(document.querySelector('.st-key-chart_viability .js-plotly-plot')?.layout.xaxis.range) === '[-2,-1]'")
            assert page.evaluate("window.savedLifecycle === document.querySelector('.st-key-chart_lifecycle .js-plotly-plot') && window.savedLifecycleData === JSON.stringify(window.savedLifecycle.data)"), "Sensitivity control rebuilt Lifecycle"
            expect(page.locator("[data-testid=stMetricValue]").first).to_have_text("1 TB")
            minimum.fill("5")
            maximum.fill("1")
            page.locator('.st-key-viability_apply button').click()
            expect(page.get_by_text("Range not applied:", exact=False)).to_be_visible()
            assert chart.evaluate("e => e.layout.xaxis.range") == [-2, -1]
            page.locator('.st-key-viability_reset button').click()
            expect(minimum).to_have_value(defaults[0])
            expect(maximum).to_have_value(defaults[1])
            page.wait_for_function("expected => JSON.stringify(document.querySelector('.st-key-chart_viability .js-plotly-plot')?.layout.xaxis.range) === JSON.stringify(expected)", arg=initial_range)
            chart.scroll_into_view_if_needed()
            page.wait_for_timeout(300)
            assert chart.evaluate("e => {const r=e.getBoundingClientRect(); return [...e.querySelectorAll('.annotation-text')].every(a => {const b=a.getBoundingClientRect(); return b.left >= r.left && b.right <= r.right;});}"), "Crossing label exceeds chart width"
            assert page.locator('[data-testid="stMain"]').evaluate("e => e.scrollWidth <= e.clientWidth + 1")
            page.screenshot(path=str(output / f"sensitivity-{width}-{theme}.png"))
            minimum.scroll_into_view_if_needed()
            page.screenshot(path=str(output / f"range-{width}-{theme}.png"))
            assert not errors, errors
            page.close()
        browser.close()
    print(f"Passed: zero-request tabs, retained charts, keyboard, signed ranges, apply/reset, mobile, both themes. Screenshots: {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--url", default="http://127.0.0.1:8501")
    parser.add_argument("--output", type=Path, default=Path(".tmp/browser-navigation"))
    parser.add_argument("--channel", default=None)
    args = parser.parse_args()
    check(args.url, args.output, args.channel)
