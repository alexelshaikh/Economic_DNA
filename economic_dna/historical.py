"""Real and fitted historical cost trends, shown alongside a scenario's
projected DNA unit costs so a viewer can see what the projection is
anchored to instead of taking the decline rate on faith."""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from .assumptions import load_assumptions

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
SEQUENCING_HISTORY_CSV = DATA_DIR / "Sequencing_Cost_Data_Table_May2022.csv"


@lru_cache(maxsize=1)
def load_observed_sequencing_costs() -> pd.DataFrame:
    """NHGRI's measured cost per Mb of sequenced DNA, nominal USD at each
    reporting date (not adjusted for the model's constant-USD convention)."""
    frame = pd.read_csv(SEQUENCING_HISTORY_CSV, parse_dates=["Date"])
    year = frame["Date"].dt.year + (frame["Date"].dt.dayofyear - 1) / 365.25
    return (
        pd.DataFrame(
            {
                "year": year,
                "sequencing_cost_usd_per_mb": frame["Cost per Mb"].astype(float),
            }
        )
        .sort_values("year")
        .reset_index(drop=True)
    )


def synthesis_historical_trend(start_year: float, end_year: float, points: int = 200) -> pd.DataFrame:
    """The paper's fitted historical synthesis-cost curve (no raw dataset
    backs synthesis cost the way NHGRI's table backs sequencing), evaluated
    with its original fixed decay rate rather than the scenario's editable
    decline rate."""
    dna = load_assumptions()["dna"]
    if end_year <= start_year:
        years = np.array([start_year, end_year], dtype=float)
    else:
        years = np.linspace(start_year, end_year, max(2, points))
    cost = (
        dna["synthesis_baseline_usd_per_base"]
        * np.exp(-dna["synthesis_historical_decay_coefficient"] * (years - dna["synthesis_baseline_year"]))
        * dna["bytes_per_mb"]
        * dna["bases_per_byte"]
    )
    return pd.DataFrame({"year": years, "synthesis_cost_usd_per_mb": cost})
