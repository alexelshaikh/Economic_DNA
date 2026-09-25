"""One-input sweeps of DNA's lifecycle cost advantage over another model."""

from __future__ import annotations

from dataclasses import dataclass, replace
from functools import lru_cache
import math

import numpy as np
import pandas as pd

from .assumptions import load_assumptions
from .scenario import Scenario
from .sensitivity import DNA_SENSITIVITY_PARAMETERS, SensitivityParameter
from .simulation import CALCULATORS, _dna_synthesis_linear_terms, _technology_label


@dataclass(frozen=True, slots=True)
class ViabilityCrossing:
    value: float
    discrete: bool = False


@dataclass(frozen=True, slots=True)
class ViabilityResult:
    parameter: SensitivityParameter
    comparison: str
    curve: pd.DataFrame
    crossings: tuple[ViabilityCrossing, ...]
    current_value: float
    current_gain: float
    use_present_value: bool


def _bounds(scenario: Scenario, field: str) -> tuple[float, float]:
    value = float(getattr(scenario, field))
    if field == "archive_size_tb":
        return max(scenario.average_asset_size_mb / 1_000_000, value / 100), min(1e9, value * 10)
    if field == "average_asset_size_mb":
        return max(value / 1e6, np.finfo(float).tiny), scenario.archive_size_mb
    if field == "dna_durability_years":
        return 1, min(10_000, max(value * 2, scenario.horizon_years + 1))
    if field == "annual_retrieval_percent":
        return 0, min(10_000, max(100, value * 2))
    if field.endswith("_percent"):
        return 0, max(99.0 if field == "discount_rate_percent" else 99.99, value)
    return 0, max(1.0, value) * 2


def _samples(low: float, high: float, integer: bool) -> np.ndarray:
    if low == high:
        return np.array([low])
    values = np.linspace(low, high, 65)
    positive_low = max(low, high * 1e-9)
    if positive_low > 0:
        values = np.concatenate((values, np.geomspace(positive_low, high, 65)))
    if integer:
        values = np.round(values)
    return np.unique(np.clip(values, low, high))


def dna_cost_advantage(
    scenario: Scenario,
    field: str,
    comparison: str,
    use_present_value: bool,
    *,
    focus: bool = True,
) -> ViabilityResult:
    """Gain is comparison cost minus DNA cost, using the same workload for both.

    Price crossings are solved analytically. Other continuous crossings are
    bracketed by a bounded sample and refined; integer durability transitions
    are marked as steps, not invented points of exact cost equality.
    """
    parameter = next((p for p in DNA_SENSITIVITY_PARAMETERS if p.field == field), None)
    if parameter is None:
        raise ValueError("Unknown DNA cost driver")
    if "DNA" not in scenario.technologies or comparison == "DNA" or comparison not in scenario.technologies:
        raise ValueError("Select DNA and a different active comparison technology")
    use_present_value = use_present_value or field == "discount_rate_percent"
    current = float(getattr(scenario, field))
    integer = field == "dna_durability_years"
    years = np.arange(scenario.start_year, scenario.start_year + scenario.horizon_years)
    elapsed = years - scenario.start_year
    assumptions = load_assumptions()

    @lru_cache(maxsize=512)
    def costs(value: float) -> tuple[float, float, float]:
        try:
            candidate = replace(scenario, **{field: int(value) if integer else float(value)})
            with np.errstate(over="raise", invalid="raise", divide="raise"):
                weights = np.power(1 + candidate.discount_rate_percent / 100, -elapsed) if use_present_value else np.ones(len(years))
                totals = [
                    float(sum(np.dot(stream, weights) for stream in CALCULATORS[technology](candidate, years, assumptions)))
                    for technology in ("DNA", comparison)
                ]
                gain = totals[1] - totals[0]
            if not all(math.isfinite(v) for v in (*totals, gain)):
                raise ValueError("Cost exceeds the numeric range")
            return totals[0], totals[1], gain
        except (ValueError, FloatingPointError, OverflowError):
            # Keep invalid extremes as gaps, never join a line across them.
            return math.nan, math.nan, math.nan

    current_costs = costs(current)
    if not math.isfinite(current_costs[2]):
        raise ValueError("The current comparison exceeds the numeric range")
    low, high = _bounds(scenario, field)
    crossings: list[ViabilityCrossing] = []
    price_field = field in ("dna_synthesis_cost_per_mb", "dna_sequencing_cost_per_mb")
    if price_field:
        if field == "dna_synthesis_cost_per_mb":
            coefficient, fixed = _dna_synthesis_linear_terms(scenario, use_present_value)
        else:
            write_coefficient, coefficient = _dna_synthesis_linear_terms(
                replace(scenario, dna_sequencing_cost_per_mb=1.0), use_present_value
            )
            fixed = write_coefficient * scenario.dna_synthesis_cost_per_mb
        if coefficient > 0:
            root = (current_costs[1] - fixed) / coefficient
            if math.isfinite(root) and root >= 0:
                crossings.append(ViabilityCrossing(root))
                high = max(high, root * 2)
    else:
        values = np.unique(np.append(_samples(low, high, integer), current))
        gains = [costs(float(x))[2] for x in values]
        all_equal = all(gain == 0 for gain in gains)
        if not all_equal:
            for index, (left, right) in enumerate(zip(values[:-1], values[1:])):
                gl, gr = gains[index:index + 2]
                if not (math.isfinite(gl) and math.isfinite(gr)):
                    continue
                if gl == 0 and index == 0:
                    crossings.append(ViabilityCrossing(float(left), integer))
                if gr == 0 and not integer:
                    if gl != 0:
                        crossings.append(ViabilityCrossing(float(right), integer))
                    continue
                if gl == 0 or (gl > 0) == (gr > 0) and gr != 0:
                    continue
                # Refine within this bracket; a geometric midpoint handles
                # thresholds many orders of magnitude below today's input.
                valid_bracket = True
                for _ in range(64):
                    if integer and right - left <= 1:
                        break
                    mid = math.sqrt(left) * math.sqrt(right) if left > 0 and right / left > 4 else left / 2 + right / 2
                    if integer:
                        mid = math.floor(mid)
                    if mid == left or mid == right:
                        break
                    gm = costs(float(mid))[2]
                    if not math.isfinite(gm):
                        valid_bracket = False
                        break
                    if gm == 0 and not integer:
                        left = right = mid
                        break
                    if (gm > 0) == (gl > 0) and gm != 0:
                        left, gl = mid, gm
                    else:
                        right, gr = mid, gm
                if valid_bracket:
                    crossings.append(ViabilityCrossing(float(right if integer else left / 2 + right / 2), integer))

    unique_crossings = []
    for crossing in sorted(crossings, key=lambda c: c.value):
        if not unique_crossings or not math.isclose(crossing.value, unique_crossings[-1].value, rel_tol=1e-8, abs_tol=0):
            unique_crossings.append(crossing)
    crossings = unique_crossings
    if focus and crossings:
        first, last = crossings[0].value, crossings[-1].value
        padding = max((last - first) * 0.5, last * 0.6, 1.0 if integer else high * 1e-12)
        low, high = max(low, first - padding), min(high, last + padding)
        if integer:
            low, high = math.floor(low), math.ceil(high)

    values = list(_samples(low, high, integer))
    if low <= current <= high:
        values.append(current)
    for crossing in crossings:
        values.append(crossing.value)
        if integer and crossing.value > low:
            values.append(crossing.value - 1)
    rows = [
        {"parameter": parameter.label, "unit": parameter.unit, "comparison": _technology_label(scenario, comparison),
         "cost_basis": "present_value" if use_present_value else "undiscounted",
         "parameter_value": float(value), "dna_cost_usd": dna, "comparison_cost_usd": other, "gain_usd": gain}
        for value in sorted(set(values))
        for dna, other, gain in [costs(float(value))]
    ]
    return ViabilityResult(
        parameter, _technology_label(scenario, comparison), pd.DataFrame(rows), tuple(crossings),
        current, current_costs[2], use_present_value,
    )
