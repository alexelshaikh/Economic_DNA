from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Any

import yaml


ASSUMPTIONS_PATH = Path(__file__).resolve().parent.parent / "assumptions.yaml"


@lru_cache(maxsize=1)
def load_assumptions() -> dict[str, Any]:
    with ASSUMPTIONS_PATH.open("r", encoding="utf-8") as handle:
        assumptions = yaml.safe_load(handle)
    if not isinstance(assumptions, dict):
        raise ValueError("assumptions.yaml must contain a mapping at its root")
    return assumptions


def model_metadata() -> dict[str, str]:
    return dict(load_assumptions()["metadata"])

