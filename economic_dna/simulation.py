from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
import pandas as pd

from .assumptions import load_assumptions, model_metadata
from .scenario import Scenario


TECHNOLOGIES = (
    "DNA",
    "Amazon Deep Archive",
    "Azure Blob Archive",
    "Tape On-premise",
    "Custom storage",
)
COMPONENTS = ("write_cost_usd", "read_cost_usd", "maintenance_cost_usd")


@dataclass(frozen=True, slots=True)
class SimulationResult:
    scenario: Scenario
    yearly: pd.DataFrame
    totals: pd.DataFrame
    metadata: dict[str, str]


def _decline_factor(percent: float, year: np.ndarray, anchor_year: float) -> np.ndarray:
    return np.power(1 - percent / 100, year - anchor_year)


def _replacement_mask(years: np.ndarray, start_year: int, durability: int) -> np.ndarray:
    elapsed = years - start_year
    return (elapsed == 0) | ((elapsed > 0) & (elapsed % durability == 0))


def _dna_costs(scenario: Scenario, years: np.ndarray, assumptions: dict) -> tuple[np.ndarray, ...]:
    write_per_mb = scenario.dna_synthesis_cost_per_mb * _decline_factor(
        scenario.synthesis_decline_percent, years, scenario.dna_cost_base_year
    )
    read_per_mb = scenario.dna_sequencing_cost_per_mb * _decline_factor(
        scenario.sequencing_decline_percent, years, scenario.dna_cost_base_year
    )
    replacements = _replacement_mask(years, scenario.start_year, scenario.dna_durability_years)
    write = write_per_mb * scenario.archive_size_mb * replacements
    read = read_per_mb * scenario.archive_size_mb * scenario.annual_retrieval_percent / 100
    return write, read, np.zeros_like(years, dtype=float)


def _amazon_costs(scenario: Scenario, years: np.ndarray, assumptions: dict) -> tuple[np.ndarray, ...]:
    factor = _decline_factor(
        scenario.amazon_decline_percent, years, scenario.amazon_price_base_year
    )
    write = scenario.amazon_put_usd_per_request * scenario.number_of_assets * factor
    read = (
        scenario.amazon_bulk_restore_usd_per_request * scenario.annual_assets_retrieved
        + scenario.amazon_bulk_retrieval_usd_per_mb
        * scenario.archive_size_mb
        * scenario.annual_retrieval_percent
        / 100
    ) * factor
    maintenance = (
        scenario.amazon_storage_usd_per_mb_month * 12 * scenario.archive_size_mb * factor
    )
    write[1:] = 0
    return write, read, maintenance


def _azure_costs(scenario: Scenario, years: np.ndarray, assumptions: dict) -> tuple[np.ndarray, ...]:
    factor = _decline_factor(
        scenario.azure_decline_percent, years, scenario.azure_price_base_year
    )
    write = scenario.azure_write_usd_per_request * scenario.number_of_assets * factor
    read = (
        scenario.azure_read_usd_per_request * scenario.annual_assets_retrieved
        + scenario.azure_retrieval_usd_per_mb
        * scenario.archive_size_mb
        * scenario.annual_retrieval_percent
        / 100
    ) * factor
    maintenance = (
        scenario.azure_storage_usd_per_mb_month * 12 * scenario.archive_size_mb * factor
    )
    write[1:] = 0
    return write, read, maintenance


def _tape_costs(scenario: Scenario, years: np.ndarray, assumptions: dict) -> tuple[np.ndarray, ...]:
    volume_tb = scenario.archive_size_tb
    media_factor = _decline_factor(
        scenario.tape_media_decline_percent, years, scenario.tape_price_base_year
    )
    hardware_factor = _decline_factor(
        scenario.tape_hardware_decline_percent, years, scenario.tape_price_base_year
    )
    energy_factor = _decline_factor(
        scenario.tape_energy_decline_percent, years, scenario.tape_price_base_year
    )
    replacements = _replacement_mask(years, scenario.start_year, scenario.tape_durability_years)
    write = scenario.tape_media_usd_per_tb * volume_tb * media_factor * replacements
    read = np.zeros_like(years, dtype=float)
    maintenance = volume_tb * (
        scenario.tape_energy_usd_per_tb_year * energy_factor
        + scenario.tape_hardware_usd_per_tb
        / scenario.tape_durability_years
        * hardware_factor
    )
    return write, read, maintenance


def _custom_costs(scenario: Scenario, years: np.ndarray, assumptions: dict) -> tuple[np.ndarray, ...]:
    factor = _decline_factor(scenario.custom_decline_percent, years, scenario.custom_cost_base_year)
    if scenario.custom_replacement_years:
        writes = _replacement_mask(years, scenario.start_year, scenario.custom_replacement_years)
    else:
        writes = years == scenario.start_year
    write = (
        scenario.custom_write_cost_per_tb * scenario.archive_size_tb
        + scenario.custom_write_cost_per_asset * scenario.number_of_assets
    ) * factor * writes
    read = (
        scenario.custom_retrieval_cost_per_tb
        * scenario.archive_size_tb
        * scenario.annual_retrieval_percent
        / 100
        + scenario.custom_retrieval_cost_per_asset * scenario.annual_assets_retrieved
    ) * factor
    maintenance = (
        scenario.custom_storage_cost_per_tb_year * scenario.archive_size_tb * factor
    )
    return write, read, maintenance


CALCULATORS: dict[str, Callable[[Scenario, np.ndarray, dict], tuple[np.ndarray, ...]]] = {
    "DNA": _dna_costs,
    "Amazon Deep Archive": _amazon_costs,
    "Azure Blob Archive": _azure_costs,
    "Tape On-premise": _tape_costs,
    "Custom storage": _custom_costs,
}


def _technology_label(scenario: Scenario, technology: str) -> str:
    return scenario.custom_storage_name.strip() if technology == "Custom storage" else technology


def simulate_scenario(scenario: Scenario) -> SimulationResult:
    assumptions = load_assumptions()
    years = np.arange(scenario.start_year, scenario.start_year + scenario.horizon_years, dtype=int)
    discount = np.power(1 + scenario.discount_rate_percent / 100, years - scenario.start_year)
    frames: list[pd.DataFrame] = []

    for technology in scenario.technologies:
        write, read, maintenance = CALCULATORS[technology](scenario, years, assumptions)
        frame = pd.DataFrame(
            {
                "technology_key": technology,
                "technology": _technology_label(scenario, technology),
                "year": years,
                "write_cost_usd": write,
                "read_cost_usd": read,
                "maintenance_cost_usd": maintenance,
            }
        )
        frame["total_cost_usd"] = frame[list(COMPONENTS)].sum(axis=1)
        frame["cumulative_cost_usd"] = frame["total_cost_usd"].cumsum()
        for component in COMPONENTS:
            frame[f"present_value_{component}"] = frame[component] / discount
        frame["present_value_total_cost_usd"] = frame[
            [f"present_value_{component}" for component in COMPONENTS]
        ].sum(axis=1)
        frame["cumulative_present_value_usd"] = frame["present_value_total_cost_usd"].cumsum()
        frames.append(frame)

    yearly = pd.concat(frames, ignore_index=True)
    aggregations = {
        component: (component, "sum") for component in COMPONENTS
    }
    aggregations.update(
        {
            "total_cost_usd": ("total_cost_usd", "sum"),
            "present_value_usd": ("present_value_total_cost_usd", "sum"),
        }
    )
    totals = (
        yearly.groupby(["technology_key", "technology"], sort=False)
        .agg(**aggregations)
        .reset_index()
    )
    totals["total_cost_per_tb_usd"] = totals["total_cost_usd"] / scenario.archive_size_tb
    totals["present_value_per_tb_usd"] = totals["present_value_usd"] / scenario.archive_size_tb
    return SimulationResult(scenario, yearly, totals, model_metadata())


def simulate_start_years(scenario: Scenario, final_start_year: int) -> pd.DataFrame:
    if final_start_year < scenario.start_year:
        raise ValueError("final_start_year must not be earlier than scenario.start_year")
    if final_start_year > 2500:
        raise ValueError("final_start_year must not exceed 2500")
    assumptions = load_assumptions()
    rows: list[dict[str, float | int | str]] = []
    for year in range(scenario.start_year, final_start_year + 1):
        current = scenario.with_start_year(year)
        years = np.arange(year, year + current.horizon_years, dtype=int)
        discount = np.power(1 + current.discount_rate_percent / 100, years - year)
        for technology in current.technologies:
            write, read, maintenance = CALCULATORS[technology](current, years, assumptions)
            total = write + read + maintenance
            rows.append(
                {
                    "start_year": year,
                    "technology_key": technology,
                    "technology": _technology_label(current, technology),
                    "write_cost_usd": float(write.sum()),
                    "read_cost_usd": float(read.sum()),
                    "maintenance_cost_usd": float(maintenance.sum()),
                    "total_cost_usd": float(total.sum()),
                    "present_value_usd": float((total / discount).sum()),
                    "total_cost_per_tb_usd": float(total.sum() / current.archive_size_tb),
                    "present_value_per_tb_usd": float(
                        (total / discount).sum() / current.archive_size_tb
                    ),
                }
            )
    return pd.DataFrame(rows)


def simulate_dna_unit_costs(scenario: Scenario, final_year: int) -> pd.DataFrame:
    if final_year < scenario.dna_cost_base_year:
        raise ValueError("final_year must not be earlier than dna_cost_base_year")
    years = np.arange(scenario.dna_cost_base_year, final_year + 1, dtype=int)
    return pd.DataFrame(
        {
            "year": years,
            "synthesis_cost_usd_per_mb": scenario.dna_synthesis_cost_per_mb
            * _decline_factor(
                scenario.synthesis_decline_percent, years, scenario.dna_cost_base_year
            ),
            "sequencing_cost_usd_per_mb": scenario.dna_sequencing_cost_per_mb
            * _decline_factor(
                scenario.sequencing_decline_percent, years, scenario.dna_cost_base_year
            ),
        }
    )


def find_crossover_years(start_year_totals: pd.DataFrame, reference: str = "DNA") -> dict[str, int | None]:
    pivot = start_year_totals.pivot(index="start_year", columns="technology", values="present_value_usd")
    if reference not in pivot:
        return {}
    crossovers: dict[str, int | None] = {}
    for technology in pivot.columns:
        if technology == reference:
            continue
        matches = pivot.index[pivot[reference] <= pivot[technology]]
        crossovers[technology] = int(matches[0]) if len(matches) else None
    return crossovers


def find_lifecycle_crossovers(
    result: SimulationResult, use_present_value: bool, reference: str = "DNA"
) -> dict[str, int | None]:
    """Within one archive's own timeline (fixed start year), the first
    calendar year in which the reference technology's running total drops to
    or below each comparison technology's. Unlike `find_crossover_years`
    (which varies the start year), this reads `result.yearly` directly."""
    value_column = "cumulative_present_value_usd" if use_present_value else "cumulative_cost_usd"
    pivot = result.yearly.pivot(index="year", columns="technology", values=value_column)
    if reference not in pivot:
        return {}
    crossovers: dict[str, int | None] = {}
    for technology in pivot.columns:
        if technology == reference:
            continue
        matches = pivot.index[pivot[reference] <= pivot[technology]]
        crossovers[technology] = int(matches[0]) if len(matches) else None
    return crossovers


def _dna_synthesis_linear_terms(scenario: Scenario, use_present_value: bool) -> tuple[float, float]:
    """DNA's lifecycle cost is `coefficient * dna_synthesis_cost_per_mb +
    fixed`: only the write/replacement component scales with the synthesis
    price, so `find_breakeven_synthesis_cost` can solve for it directly
    instead of searching."""
    years = np.arange(scenario.start_year, scenario.start_year + scenario.horizon_years, dtype=int)
    decline = _decline_factor(scenario.synthesis_decline_percent, years, scenario.dna_cost_base_year)
    replacements = _replacement_mask(years, scenario.start_year, scenario.dna_durability_years)
    coefficient_terms = decline * replacements * scenario.archive_size_mb

    read_per_mb = scenario.dna_sequencing_cost_per_mb * _decline_factor(
        scenario.sequencing_decline_percent, years, scenario.dna_cost_base_year
    )
    fixed_terms = read_per_mb * scenario.archive_size_mb * scenario.annual_retrieval_percent / 100

    if use_present_value:
        discount = np.power(1 + scenario.discount_rate_percent / 100, years - scenario.start_year)
        coefficient_terms = coefficient_terms / discount
        fixed_terms = fixed_terms / discount

    return float(coefficient_terms.sum()), float(fixed_terms.sum())


BREAKEVEN_COLUMNS = ["technology", "breakeven_synthesis_cost_usd_per_mb", "reduction_factor"]


def find_breakeven_synthesis_cost(scenario: Scenario, use_present_value: bool) -> pd.DataFrame:
    """USD/MB synthesis cost at which DNA's lifecycle cost equals each
    comparison technology's, holding every other input fixed. `None` means
    DNA cannot reach that technology's cost even with free synthesis."""
    if "DNA" not in scenario.technologies:
        return pd.DataFrame(columns=BREAKEVEN_COLUMNS)

    coefficient, fixed_cost = _dna_synthesis_linear_terms(scenario, use_present_value)
    value_column = "present_value_usd" if use_present_value else "total_cost_usd"
    totals = simulate_scenario(scenario).totals.set_index("technology_key")[value_column]

    rows: list[dict[str, float | str | None]] = []
    for technology in scenario.technologies:
        if technology == "DNA":
            continue
        target = float(totals[technology])
        breakeven = (target - fixed_cost) / coefficient if coefficient > 0 else None
        if breakeven is not None and breakeven < 0:
            breakeven = None
        reduction_factor = (
            scenario.dna_synthesis_cost_per_mb / breakeven
            if breakeven is not None and breakeven > 0
            else None
        )
        rows.append(
            {
                "technology": _technology_label(scenario, technology),
                "breakeven_synthesis_cost_usd_per_mb": breakeven,
                "reduction_factor": reduction_factor,
            }
        )
    return pd.DataFrame(rows, columns=BREAKEVEN_COLUMNS)


def simulate_dna_uncertainty_band(
    scenario: Scenario,
    use_present_value: bool,
    n_samples: int = 300,
    relative_spread: float = 0.3,
    seed: int = 0,
) -> pd.DataFrame:
    """P10/P50/P90 band for DNA's cumulative cost, from sampling the
    synthesis and sequencing decline rates within +/- relative_spread of the
    scenario's chosen values. Fully vectorized over (sample, year), no
    per-sample simulation loop."""
    if n_samples < 2:
        raise ValueError("n_samples must be at least 2")
    rng = np.random.default_rng(seed)
    years = np.arange(scenario.start_year, scenario.start_year + scenario.horizon_years, dtype=int)

    def _sampled_decline(base_percent: float) -> np.ndarray:
        sampled = base_percent * (1 + rng.uniform(-relative_spread, relative_spread, n_samples))
        return np.clip(sampled, 0.0, 99.999)[:, None]

    synthesis_factor = _decline_factor(
        _sampled_decline(scenario.synthesis_decline_percent), years[None, :], scenario.dna_cost_base_year
    )
    sequencing_factor = _decline_factor(
        _sampled_decline(scenario.sequencing_decline_percent), years[None, :], scenario.dna_cost_base_year
    )
    replacements = _replacement_mask(years, scenario.start_year, scenario.dna_durability_years)

    write = (
        scenario.dna_synthesis_cost_per_mb
        * synthesis_factor
        * scenario.archive_size_mb
        * replacements[None, :]
    )
    read = (
        scenario.dna_sequencing_cost_per_mb
        * sequencing_factor
        * scenario.archive_size_mb
        * scenario.annual_retrieval_percent
        / 100
    )
    total = write + read
    if use_present_value:
        discount = np.power(1 + scenario.discount_rate_percent / 100, years - scenario.start_year)
        total = total / discount[None, :]

    cumulative = np.cumsum(total, axis=1)
    p10, p50, p90 = np.quantile(cumulative, [0.1, 0.5, 0.9], axis=0)
    return pd.DataFrame({"year": years, "p10": p10, "p50": p50, "p90": p90})
