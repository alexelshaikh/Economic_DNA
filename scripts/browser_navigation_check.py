"""Check tab transitions frame by frame, including delayed server responses."""
from __future__ import annotations

import argparse
from pathlib import Path

from playwright.sync_api import expect, sync_playwright


def check(url: str, output: Path, channel: str | None) -> None:
    output.mkdir(parents=True, exist_ok=True)
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(channel=channel, headless=True)
        for width, theme in ((1440, "light"), (1440, "dark"), (390, "light"), (320, "dark")):
            page = browser.new_page(viewport={"width": width, "height": 1000}, color_scheme=theme)
            queue = []
            paused = [False]
            errors = []
            sent = []

            def route(socket):
                server = socket.connect_to_server()

                def receive(message):
                    if paused[0]:
                        queue.append((socket, message))
                    else:
                        socket.send(message)

                server.on_message(receive)

            def resume():
                while queue:
                    socket, message = queue.pop(0)
                    socket.send(message)
                    page.wait_for_timeout(8)
                paused[0] = False

            page.route_web_socket("**/_stcore/stream", route)
            page.on("pageerror", lambda error: errors.append(str(error)))
            page.on("websocket", lambda ws: ws.on("framesent", lambda payload: sent.append(payload)))
            page.goto(url + "?theme=" + theme)
            expect(page.locator(".js-plotly-plot")).to_have_count(2)
            page.wait_for_function("window.__dnaAnalysisNavigation !== undefined")
            expect(page.locator(".analysis-ready")).to_have_count(1)
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
                const sample = () => {
                    const root = document.querySelector('.st-key-analysis_tabs');
                    const tabs = root.querySelector('[role="tablist"]');
                    const preview = root.querySelector('.analysis-preview');
                    const panel = root.querySelector('[role="tabpanel"]');
                    const view = preview || panel;
                    let opacity = 1;
                    for (let element = tabs; element; element = element.parentElement) {
                        opacity *= Number(getComputedStyle(element).opacity);
                    }
                    const visible = view && getComputedStyle(view).visibility === 'visible'
                        && getComputedStyle(view).display !== 'none';
                    window.navigationSamples.push({same:tabs===window.originalTabs, opacity,
                        y:tabs.getBoundingClientRect().y, filled:visible && view.textContent.trim().length > 0,
                        previews:root.querySelectorAll('.analysis-preview').length});
                    window.navigationFrame = requestAnimationFrame(sample);
                }; sample();
            }""")

            for label, charts in (("Start-year outlook", 1), ("DNA unit costs", 2), ("Sensitivity", 3), ("Assumptions", 0), ("About", 0), ("Lifecycle", 2)):
                frames_before = len(sent)
                original_bounds = page.locator('[role="tabpanel"] .js-plotly-plot').evaluate_all("nodes => nodes.map(n => {const r=n.getBoundingClientRect(); return [r.x,r.y,r.width,r.height];})")
                paused[0] = True
                page.get_by_role("tab", name=label, exact=True).click()
                expect(page.locator(".analysis-preview")).to_have_count(1)
                expect(page.locator(".analysis-loading-status")).to_have_text(f"Loading {label}...")
                page.wait_for_timeout(650)
                expect(tabs).to_have_css("opacity", "1")
                expect(page.locator(".analysis-preview")).to_have_attribute("inert", "")
                snapshot_bounds = page.locator('.analysis-preview .js-plotly-plot').evaluate_all("nodes => nodes.map(n => {const r=n.getBoundingClientRect(); return [r.x,r.y,r.width,r.height];})")
                assert len(snapshot_bounds) == len(original_bounds)
                assert all(abs(before - after) < 1 for original, snapshot in zip(original_bounds, snapshot_bounds) for before, after in zip(original, snapshot)), f"Retained charts shifted: {label}, {width}: {original_bounds} -> {snapshot_bounds}"
                notice_bounds = page.locator(".analysis-loading-status").bounding_box()
                heading_bounds = page.locator(".analysis-preview .tab-intro h2").first.bounding_box()
                assert notice_bounds["y"] + notice_bounds["height"] <= heading_bounds["y"], "Loading notice overlaps the retained view"
                assert page.locator(".analysis-preview [class*='st-key-'], .analysis-preview iframe").count() == 0
                assert page.locator(".analysis-preview [id]").evaluate_all("nodes => nodes.every(n => n.id.startsWith('analysis-preview-'))")
                assert page.evaluate("""() => [...document.querySelectorAll('.analysis-preview [clip-path]')].every(n => {
                    const id = n.getAttribute('clip-path').match(/url\\(#([^)]*)\\)/)?.[1];
                    return !id || document.getElementById(id)?.closest('.analysis-preview');
                })"""), "Snapshot chart clipping references are broken"
                if label == "Start-year outlook":
                    page.screenshot(path=str(output / f"loading-{width}-{theme}.png"))
                resume()
                expect(page.locator(".analysis-preview")).to_have_count(0)
                expect(page.locator('.st-key-analysis_tabs[aria-busy]')).to_have_count(0)
                expect(page.locator(".js-plotly-plot")).to_have_count(charts)
                assert len(sent) == frames_before + 1, "Switching a tab caused extra server requests"
                expect(page.locator(".pending-notice")).to_have_text("Recalculate to update charts.")
                expect(page.locator("[data-testid=stMetricValue]").first).to_have_text("1 TB")

            # In-flight results for an earlier click must not reveal the wrong view.
            paused[0] = True
            for label in ("Start-year outlook", "Sensitivity", "About"):
                page.get_by_role("tab", name=label, exact=True).click()
            expect(page.locator(".analysis-preview")).to_have_count(1)
            page.wait_for_timeout(150)
            resume()
            expect(page.locator(".analysis-preview")).to_have_count(0)
            expect(page.get_by_role("tab", name="About", exact=True)).to_have_attribute("aria-selected", "true")
            expect(page.get_by_role("tabpanel")).to_contain_text("About this explorer")

            # Arrow-key navigation follows the same guarded transition.
            page.get_by_role("tab", name="About", exact=True).focus()
            paused[0] = True
            page.keyboard.press("Home")
            expect(page.locator(".analysis-preview")).to_have_count(1)
            page.wait_for_timeout(150)
            resume()
            expect(page.locator(".analysis-preview")).to_have_count(0)
            expect(page.get_by_role("tabpanel")).to_contain_text("Lifecycle comparison")
            page.evaluate("cancelAnimationFrame(window.navigationFrame)")
            samples = page.evaluate("window.navigationSamples")
            assert all(s["same"] for s in samples), "The tab bar was remounted"
            assert all(s["opacity"] == 1 for s in samples), "The tab bar faded during navigation"
            assert all(s["filled"] for s in samples), "An empty/hidden content frame was displayed"
            assert max(s["previews"] for s in samples) == 1, "More than one transition snapshot was retained"
            assert max(s["y"] for s in samples) - min(s["y"] for s in samples) < 2, "Tab navigation moved the page"
            page.screenshot(path=str(output / f"complete-{width}-{theme}.png"))
            assert not errors, errors
            page.close()
        browser.close()
    print(f"Navigation checks passed; delayed responses, rapid clicks, keyboard, mobile, and both themes. Screenshots: {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--url", default="http://127.0.0.1:8501")
    parser.add_argument("--output", type=Path, default=Path(".tmp/browser-navigation"))
    parser.add_argument("--channel", default=None, help="Installed browser, for example msedge or chrome")
    args = parser.parse_args()
    check(args.url, args.output, args.channel)
