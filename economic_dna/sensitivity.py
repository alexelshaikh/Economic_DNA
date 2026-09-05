"""Tornado-style sensitivity analysis: how much DNA's lifecycle cost moves
when one input is perturbed while everything else stays fixed."""

from __future__ import annotations

from dataclasses import dataclass, replace

import pandas as pd

from .scenario import Scenario
from .simulation import simulate_scenario


@dataclass(frozen=True, slots=True)
class SensitivityParameter:
    field: str
    label: str
    unit: str = ""
    relative_delta: float = 0.5
    # Used instead of relative_delta when the baseline value is exactly 0
    # (a relative swing around 0 is always 0).
    absolute_fallback: float = 0.0


DNA_SENSITIVITY_PARAMETERS: tuple[SensitivityParameter, ...] = (
    SensitivityParameter("dna_synthesis_cost_per_mb", "Synthesis cost", "$/MB"),
    SensitivityParameter("dna_sequencing_cost_per_mb", "Sequencing cost", "$/MB"),
    SensitivityParameter("synthesis_decline_percent", "Synthesis decline rate", "%/yr"),
    SensitivityParameter("sequencing_decline_percent", "Sequencing decline rate", "%/yr"),
    SensitivityParameter("archive_size_tb", "Archive size", "TB"),
    SensitivityParameter("average_asset_size_mb", "Average asset size", "MB"),
    SensitivityParameter(
        "annual_retrieval_percent", "Annual retrieval", "%/yr", absolute_fallback=1.0
    ),
    SensitivityParameter("discount_rate_percent", "Discount rate", "%/yr", absolute_fallback=5.0),
    SensitivityParameter("dna_durability_years", "DNA durability", "years"),
)

SENSITIVITY_COLUMNS = [
    "parameter",
    "unit",
    "low_value",
    "high_value",
    "low_cost",
    "high_cost",
    "baseline_cost",
    "swing",
]


def _perturbed_value(base: float, sign: int, spec: SensitivityParameter) -> float:
    candidate = base * (1 + sign * spec.relative_delta) if base else sign * spec.absolute_fallback
    return max(candidate, 0.0)


def _dna_cost(scenario: Scenario, value_column: str, fallback: float) -> float:
    try:
        totals = simulate_scenario(scenario).totals.set_index("technology_key")
    except ValueError:
        return fallback
    if "DNA" not in totals.index:
        return fallback
    return float(totals.loc["DNA", value_column])


def dna_cost_sensitivity(
    scenario: Scenario,
    use_present_value: bool,
    parameters: tuple[SensitivityParameter, ...] = DNA_SENSITIVITY_PARAMETERS,
) -> pd.DataFrame:
    """For each parameter, hold everything else fixed and re-simulate at a
    lowered and a raised value. `swing` is the resulting spread in DNA's
    lifecycle cost, the standard tornado-chart sort key."""
    if "DNA" not in scenario.technologies:
        return pd.DataFrame(columns=SENSITIVITY_COLUMNS)

    value_column = "present_value_usd" if use_present_value else "total_cost_usd"
    baseline = _dna_cost(scenario, value_column, fallback=0.0)

    rows = []
    for spec in parameters:
        base_value = getattr(scenario, spec.field)
        low_value = _perturbed_value(base_value, -1, spec)
        high_value = _perturbed_value(base_value, 1, spec)
        low_cost = _dna_cost(replace(scenario, **{spec.field: low_value}), value_column, baseline)
        high_cost = _dna_cost(replace(scenario, **{spec.field: high_value}), value_column, baseline)
        rows.append(
            {
                "parameter": spec.label,
                "unit": spec.unit,
                "low_value": low_value,
                "high_value": high_value,
                "low_cost": low_cost,
                "high_cost": high_cost,
                "baseline_cost": baseline,
                "swing": abs(high_cost - low_cost),
            }
        )
    return pd.DataFrame(rows, columns=SENSITIVITY_COLUMNS).sort_values(
        "swing", ascending=False
    ).reset_index(drop=True)
