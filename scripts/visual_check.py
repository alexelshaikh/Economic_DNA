from __future__ import annotations

import argparse
import time
from pathlib import Path

from PIL import Image, ImageStat
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait


def capture(url: str, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    options = webdriver.EdgeOptions()
    options.add_argument("--headless=new")
    options.add_argument("--disable-gpu")
    options.add_argument("--no-first-run")
    driver = webdriver.Edge(options=options)
    try:
        driver.set_window_size(1440, 1100)
        driver.get(url)
        WebDriverWait(driver, 30).until(
            lambda page: len(page.find_elements(By.CSS_SELECTOR, '[data-testid="stMetric"]')) == 4
        )
        WebDriverWait(driver, 30).until(
            lambda page: len(page.find_elements(By.CSS_SELECTOR, '[data-testid="stPlotlyChart"]')) >= 2
        )
        time.sleep(2)
        desktop = output_dir / "desktop.png"
        driver.save_screenshot(str(desktop))

        dna_tab = WebDriverWait(driver, 10).until(
            lambda page: (
                page.find_elements(By.XPATH, "//*[normalize-space(.)='DNA unit costs']")[-1]
                if page.find_elements(By.XPATH, "//*[normalize-space(.)='DNA unit costs']")
                else None
            )
        )
        dna_tab.click()
        WebDriverWait(driver, 20).until(
            lambda page: "DNA synthesis cost trajectory" in page.page_source
            and "DNA sequencing cost trajectory" in page.page_source
        )
        time.sleep(1)
        dna_costs = output_dir / "dna-unit-costs.png"
        driver.save_screenshot(str(dna_costs))

        driver.set_window_size(430, 932)
        time.sleep(2)
        mobile = output_dir / "mobile.png"
        driver.save_screenshot(str(mobile))

        metrics = [item.text.replace("\n", ": ") for item in driver.find_elements(By.CSS_SELECTOR, '[data-testid="stMetric"]')]
        print(f"Rendered metrics: {metrics}")
        for screenshot in (desktop, dna_costs, mobile):
            with Image.open(screenshot) as image:
                variance = sum(ImageStat.Stat(image.convert("RGB")).var)
                if variance < 100:
                    raise RuntimeError(f"{screenshot} appears blank (variance={variance:.2f})")
                print(f"{screenshot}: {image.width}x{image.height}, variance={variance:.2f}")
    finally:
        driver.quit()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--url", default="http://127.0.0.1:8501")
    parser.add_argument("--output", type=Path, default=Path(".tmp") / "visual-check")
    arguments = parser.parse_args()
    capture(arguments.url, arguments.output)
