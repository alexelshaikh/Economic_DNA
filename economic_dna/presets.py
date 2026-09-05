"""One-click starting points, so a first-time visitor does not have to fill
in ~40 fields before seeing a meaningful comparison."""

from __future__ import annotations

from .scenario import Scenario

PRESET_SCENARIOS: dict[str, Scenario] = {
    "Paper baseline": Scenario(),
    "1 PB genomics cold archive": Scenario(
        archive_size_tb=1_000.0,
        average_asset_size_mb=50.0,
        annual_retrieval_percent=0.5,
        horizon_years=50,
    ),
    "10 PB media archive, frequent retrieval": Scenario(
        archive_size_tb=10_000.0,
        average_asset_size_mb=250.0,
        annual_retrieval_percent=20.0,
        horizon_years=25,
    ),
    "100 TB lab archive": Scenario(
        archive_size_tb=100.0,
        average_asset_size_mb=5.0,
        annual_retrieval_percent=2.0,
        horizon_years=15,
    ),
}
