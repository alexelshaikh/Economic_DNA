from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import math
from typing import Any, Mapping


DEFAULT_SYNTHESIS_DECLINE = 0.1538850043930884
DEFAULT_SEQUENCING_DECLINE = 0.38065945245876087
DEFAULT_TECHNOLOGIES = ("DNA", "Amazon Deep Archive", "Azure Blob Archive", "Tape On-premise")
AVAILABLE_TECHNOLOGIES = DEFAULT_TECHNOLOGIES + ("Custom storage",)


@dataclass(frozen=True, slots=True)
class Scenario:
    archive_size_tb: float = 1.0
    average_asset_size_mb: float = 1000.0
    annual_retrieval_percent: float = 1.0
    start_year: int = 2025
    horizon_years: int = 100
    discount_rate_percent: float = 0.0
    dna_cost_base_year: int = 2026
    dna_synthesis_cost_per_mb: float = 16573.838543736347
    dna_sequencing_cost_per_mb: float = 0.13475734873958314
    synthesis_decline_percent: float = DEFAULT_SYNTHESIS_DECLINE * 100
    sequencing_decline_percent: float = DEFAULT_SEQUENCING_DECLINE * 100
    amazon_price_base_year: int = 2025
    amazon_put_usd_per_request: float = 0.00005
    amazon_bulk_restore_usd_per_request: float = 0.000025
    amazon_bulk_retrieval_usd_per_mb: float = 0.00000244140625
    amazon_storage_usd_per_mb_month: float = 0.000000966796875
    amazon_decline_percent: float = 10.0
    azure_price_base_year: int = 2025
    azure_write_usd_per_request: float = 0.00001
    azure_read_usd_per_request: float = 0.0005
    azure_retrieval_usd_per_mb: float = 0.00001953125
    azure_storage_usd_per_mb_month: float = 0.000001953125
    azure_decline_percent: float = 10.0
    tape_price_base_year: int = 2025
    tape_media_usd_per_tb: float = 6.39
    tape_hardware_usd_per_tb: float = 6.86
    tape_energy_usd_per_tb_year: float = 0.05
    tape_media_decline_percent: float = 20.0
    tape_hardware_decline_percent: float = 0.0
    tape_energy_decline_percent: float = 15.0
    dna_durability_years: int = 1000
    tape_durability_years: int = 30
    custom_storage_name: str = "Custom storage"
    custom_cost_base_year: int = 2026
    custom_write_cost_per_tb: float = 0.0
    custom_write_cost_per_asset: float = 0.0
    custom_storage_cost_per_tb_year: float = 0.0
    custom_retrieval_cost_per_tb: float = 0.0
    custom_retrieval_cost_per_asset: float = 0.0
    custom_decline_percent: float = 0.0
    custom_replacement_years: int = 0
    technologies: tuple[str, ...] = DEFAULT_TECHNOLOGIES

    def __post_init__(self) -> None:
        errors: list[str] = []
        if not (0 < self.archive_size_tb <= 1_000_000_000):
            errors.append("archive_size_tb must be greater than 0 and at most 1 billion TB")
        if not (0 < self.average_asset_size_mb <= self.archive_size_mb):
            errors.append("average_asset_size_mb must be greater than 0 and no larger than the archive")
        if not (0 <= self.annual_retrieval_percent <= 10_000):
            errors.append("annual_retrieval_percent must be between 0 and 10,000")
        if not (2025 <= self.start_year <= 2500):
            errors.append("start_year must be between 2025 and 2500")
        if not (1 <= self.horizon_years <= 10_000):
            errors.append("horizon_years must be between 1 and 10,000")
        if not (0 <= self.discount_rate_percent < 100):
            errors.append("discount_rate_percent must be at least 0 and below 100")
        if not (2000 <= self.dna_cost_base_year <= 2500):
            errors.append("dna_cost_base_year must be between 2000 and 2500")
        for name, value in {
            "amazon_price_base_year": self.amazon_price_base_year,
            "azure_price_base_year": self.azure_price_base_year,
            "tape_price_base_year": self.tape_price_base_year,
        }.items():
            if not (2000 <= value <= 2500):
                errors.append(f"{name} must be between 2000 and 2500")
        if not (2000 <= self.custom_cost_base_year <= 2500):
            errors.append("custom_cost_base_year must be between 2000 and 2500")

        percentage_fields = {
            "synthesis_decline_percent": self.synthesis_decline_percent,
            "sequencing_decline_percent": self.sequencing_decline_percent,
            "amazon_decline_percent": self.amazon_decline_percent,
            "azure_decline_percent": self.azure_decline_percent,
            "tape_media_decline_percent": self.tape_media_decline_percent,
            "tape_hardware_decline_percent": self.tape_hardware_decline_percent,
            "tape_energy_decline_percent": self.tape_energy_decline_percent,
            "custom_decline_percent": self.custom_decline_percent,
        }
        for name, value in percentage_fields.items():
            if not (0 <= value < 100):
                errors.append(f"{name} must be at least 0 and below 100")

        if self.dna_durability_years < 1 or self.tape_durability_years < 1:
            errors.append("durability values must be positive whole years")
        if not (0 <= self.custom_replacement_years <= 10_000):
            errors.append("custom_replacement_years must be between 0 and 10,000")

        cost_fields = {
            "dna_synthesis_cost_per_mb": self.dna_synthesis_cost_per_mb,
            "dna_sequencing_cost_per_mb": self.dna_sequencing_cost_per_mb,
            "amazon_put_usd_per_request": self.amazon_put_usd_per_request,
            "amazon_bulk_restore_usd_per_request": self.amazon_bulk_restore_usd_per_request,
            "amazon_bulk_retrieval_usd_per_mb": self.amazon_bulk_retrieval_usd_per_mb,
            "amazon_storage_usd_per_mb_month": self.amazon_storage_usd_per_mb_month,
            "azure_write_usd_per_request": self.azure_write_usd_per_request,
            "azure_read_usd_per_request": self.azure_read_usd_per_request,
            "azure_retrieval_usd_per_mb": self.azure_retrieval_usd_per_mb,
            "azure_storage_usd_per_mb_month": self.azure_storage_usd_per_mb_month,
            "tape_media_usd_per_tb": self.tape_media_usd_per_tb,
            "tape_hardware_usd_per_tb": self.tape_hardware_usd_per_tb,
            "tape_energy_usd_per_tb_year": self.tape_energy_usd_per_tb_year,
            "custom_write_cost_per_tb": self.custom_write_cost_per_tb,
            "custom_write_cost_per_asset": self.custom_write_cost_per_asset,
            "custom_storage_cost_per_tb_year": self.custom_storage_cost_per_tb_year,
            "custom_retrieval_cost_per_tb": self.custom_retrieval_cost_per_tb,
            "custom_retrieval_cost_per_asset": self.custom_retrieval_cost_per_asset,
        }
        for name, value in cost_fields.items():
            if not math.isfinite(value) or value < 0:
                errors.append(f"{name} must be a finite value of at least 0")

        custom_name = self.custom_storage_name.strip()
        if not custom_name or len(custom_name) > 60:
            errors.append("custom_storage_name must contain between 1 and 60 characters")
        reserved_names = {name.lower() for name in DEFAULT_TECHNOLOGIES}
        if custom_name.lower() in reserved_names:
            errors.append("custom_storage_name must differ from the built-in technology names")
        if not self.technologies:
            errors.append("at least one technology must be selected")
        unknown = sorted(set(self.technologies) - set(AVAILABLE_TECHNOLOGIES))
        if unknown:
            errors.append(f"unknown technologies: {', '.join(unknown)}")
        if len(set(self.technologies)) != len(self.technologies):
            errors.append("technologies must not contain duplicates")
        if errors:
            raise ValueError("; ".join(errors))

    @property
    def archive_size_mb(self) -> float:
        return self.archive_size_tb * 1_000_000

    @property
    def number_of_assets(self) -> float:
        return self.archive_size_mb / self.average_asset_size_mb

    @property
    def annual_assets_retrieved(self) -> float:
        return self.number_of_assets * self.annual_retrieval_percent / 100

    def with_start_year(self, year: int) -> "Scenario":
        return replace(self, start_year=year)

    def to_query_params(self) -> dict[str, str]:
        values = asdict(self)
        values["technologies"] = ",".join(self.technologies)
        return {name: str(value) for name, value in values.items()}

    @classmethod
    def from_mapping(cls, values: Mapping[str, Any]) -> "Scenario":
        defaults = cls()
        parsed: dict[str, Any] = {}
        int_fields = {
            "start_year",
            "horizon_years",
            "dna_cost_base_year",
            "dna_durability_years",
            "amazon_price_base_year",
            "azure_price_base_year",
            "tape_price_base_year",
            "tape_durability_years",
            "custom_cost_base_year",
            "custom_replacement_years",
        }
        float_fields = {
            "archive_size_tb",
            "average_asset_size_mb",
            "annual_retrieval_percent",
            "discount_rate_percent",
            "dna_synthesis_cost_per_mb",
            "dna_sequencing_cost_per_mb",
            "synthesis_decline_percent",
            "sequencing_decline_percent",
            "amazon_put_usd_per_request",
            "amazon_bulk_restore_usd_per_request",
            "amazon_bulk_retrieval_usd_per_mb",
            "amazon_storage_usd_per_mb_month",
            "amazon_decline_percent",
            "azure_write_usd_per_request",
            "azure_read_usd_per_request",
            "azure_retrieval_usd_per_mb",
            "azure_storage_usd_per_mb_month",
            "azure_decline_percent",
            "tape_media_usd_per_tb",
            "tape_hardware_usd_per_tb",
            "tape_energy_usd_per_tb_year",
            "tape_media_decline_percent",
            "tape_hardware_decline_percent",
            "tape_energy_decline_percent",
            "custom_write_cost_per_tb",
            "custom_write_cost_per_asset",
            "custom_storage_cost_per_tb_year",
            "custom_retrieval_cost_per_tb",
            "custom_retrieval_cost_per_asset",
            "custom_decline_percent",
        }
        for name in int_fields:
            if name in values:
                parsed[name] = int(values[name])
        for name in float_fields:
            if name in values:
                parsed[name] = float(values[name])
        if "technologies" in values:
            raw = values["technologies"]
            parsed["technologies"] = tuple(raw.split(",")) if isinstance(raw, str) else tuple(raw)
        if "custom_storage_name" in values:
            parsed["custom_storage_name"] = str(values["custom_storage_name"])
        return replace(defaults, **parsed)
