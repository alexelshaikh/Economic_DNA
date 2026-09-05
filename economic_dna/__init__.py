"""Public calculation API for the Economic DNA storage model."""

from .historical import load_observed_sequencing_costs, synthesis_historical_trend
from .presets import PRESET_SCENARIOS
from .scenario import Scenario
from .sensitivity import DNA_SENSITIVITY_PARAMETERS, dna_cost_sensitivity
from .simulation import (
    TECHNOLOGIES,
    SimulationResult,
    find_breakeven_synthesis_cost,
    find_crossover_years,
    find_lifecycle_crossovers,
    simulate_scenario,
    simulate_dna_unit_costs,
    simulate_dna_uncertainty_band,
    simulate_start_years,
)

__all__ = [
    "Scenario",
    "SimulationResult",
    "TECHNOLOGIES",
    "PRESET_SCENARIOS",
    "DNA_SENSITIVITY_PARAMETERS",
    "simulate_scenario",
    "simulate_dna_unit_costs",
    "simulate_dna_uncertainty_band",
    "simulate_start_years",
    "find_crossover_years",
    "find_lifecycle_crossovers",
    "find_breakeven_synthesis_cost",
    "dna_cost_sensitivity",
    "load_observed_sequencing_costs",
    "synthesis_historical_trend",
]
