from sympy import Symbol, exp
from models import model3
import numpy as np

def _synthesis_cost(year,
                           drop_rate: float = 0.1671,
                           anchor_year: float = 2000):
    # 1) your true historical baseline:
    orig_rate              = 0.1671     # fixed, from the original fit
    baseline_year          = 2000
    baseline_cost_per_base = 0.3193

    # 2) roll that baseline to your chosen anchor, using the ORIGINAL rate
    cost_per_base_at_anchor = (
        baseline_cost_per_base
        * np.exp(-orig_rate * (anchor_year - baseline_year))
    )

    # 3) from that anchor, apply whichever drop_rate you like:
    years_since_anchor = year - anchor_year
    cost = (
        cost_per_base_at_anchor
        * exp(-drop_rate * years_since_anchor)
        * 1_000_000
        * 4
    )
    return cost
def _sequencing_cost(year, drop_rate=0.4791, anchor_year=2000):
    # 1) original “baseline” point:
    baseline_year        = 2011.5
    baseline_cost_per_run = 35.036    # USD per run at baseline_year
    multiplier            = 4         # your original ×4 factor

    # 2) roll that baseline back or forward to your anchor using the *original* rate
    orig_rate = 0.4791
    cost_at_anchor = (
        baseline_cost_per_run
        * np.exp(-orig_rate * (anchor_year - baseline_year))
        * multiplier
    )
    # 3) from the anchor, decay/grow to the target year with the chosen drop_rate
    return cost_at_anchor * exp(-drop_rate * (year - anchor_year))

def dna(n=1000, k=10, obj_size_mb=100, start_year=2025, d=100, device_durability=1000, synthesis_anchor_year=2025, synthesis_drop_rate=0.1671, sequencing_anchor_year=2025, sequencing_drop_rate=0.4791):
    _t = Symbol('t', integer=True)
    dna_model = model3.CostModel("DNA")
    dna_model.cost_write = _synthesis_cost(_t, anchor_year=synthesis_anchor_year, drop_rate=synthesis_drop_rate)
    dna_model.cost_read = _sequencing_cost(_t, anchor_year=sequencing_anchor_year, drop_rate=sequencing_drop_rate)
    dna_model.cost_maint = 0
    dna_model.max_repls = d // device_durability
    dna_model.device_durability = device_durability
    dna_model.n = n
    dna_model.k = k
    dna_model.obj_size = obj_size_mb
    dna_model.start_year = start_year
    dna_model.d = d
    return dna_model

def azure(n=1000, k=10, obj_size_mb=100, start_year=2025, d=100, drop_rate=0.1):
    t = Symbol('t', integer=True)

    # ----------------------------------------------------------------------
    # Base unit prices (US East, PAYG, mid‑2025)
    # https://azure.microsoft.com/en-us/pricing/details/storage/blobs/#pricing
    # Flat Namespace, Region East US 2
    write_per_op       = 0.10  / 10_000        # USD per write request
    read_per_op        = 5     / 10_000        # USD per read request
    retrieval_per_mb   = 0.02  / 1024          # USD per MB retrieved
    storage_per_mb_mon = 0.002 / 1024          # USD per MB‑month stored
    # ----------------------------------------------------------------------

    azure_model = model3.CostModel(f"azure_blob_archive_{drop_rate}")
    # one write request (independent of size)
    azure_model.cost_write = (write_per_op * (1 - drop_rate)**(t - 2025)) / obj_size_mb

    # one read request + obj_size_mb MB of data
    azure_model.cost_read = ((
        read_per_op +
        retrieval_per_mb * obj_size_mb
    ) * (1 - drop_rate)**(t - 2025)) / obj_size_mb

    # storing obj_size_mb for 12 months
    azure_model.cost_maint = ((
        storage_per_mb_mon * obj_size_mb * 12
    ) * (1 - drop_rate)**(t - 2025)) / obj_size_mb

    azure_model.n = n
    azure_model.k = k
    azure_model.obj_size = obj_size_mb
    azure_model.start_year = start_year
    azure_model.d = d
    azure_model.max_repls = 0
    azure_model.device_durability = d

    return azure_model


def amazon(n=1000, k=10, obj_size_mb=100, start_year=2025, d=100, drop_rate=0.1):
    t = Symbol('t', integer=True)

    # ----------------------------------------------------------------------
    # Base unit prices (US East, PAYG, mid‑2025)
    # https://aws.amazon.com/s3/pricing/?nc=sn&loc=4
    # US East (N. Virginia)
    put_per_req = 0.05 / 1_000  # direct PUT to Deep Archive
    bulk_request_per_req = 0.025 / 1_000  # restore-job invocation
    bulk_data_per_mb = 0.0025 / 1024  # data retrieval
    storage_per_mb_mon = 0.00099 / 1024  # storage

    amazon_model = model3.CostModel(f"amazon_deep_archive_{drop_rate}")

    # one direct write request
    amazon_model.cost_write = (put_per_req * (1 - drop_rate) ** (t - 2025)) / obj_size_mb

    # one restore-job request + obj_size_mb MB of data
    amazon_model.cost_read = (
            (bulk_request_per_req / obj_size_mb + bulk_data_per_mb)
            * (1 - drop_rate) ** (t - 2025)
    )

    # storing obj_size_mb for 12 months
    amazon_model.cost_maint = ((
                                      storage_per_mb_mon * 12
                              ) * (1 - drop_rate) ** (t - 2025))
    amazon_model.n = n
    amazon_model.k = k
    amazon_model.obj_size = obj_size_mb
    amazon_model.start_year = start_year
    amazon_model.d = d
    amazon_model.max_repls = 0
    amazon_model.device_durability = d

    return amazon_model


def tape_on_premise(improv_rate_per_year=0.1, write_cost_per_mb=0.00000639, read_cost_per_mb=0, maint_cost_per_mb_per_year=0.00000671, device_durability_years=30, d=100, n=1000, k=10, obj_size_mb=100, start_year=2025, model_name=None):
    _t = Symbol('t', integer=True)

    tape_on_premise = model3.CostModel(f"Tape_on_premise_{improv_rate_per_year}") if model_name is None else model3.CostModel(model_name)

    tape_on_premise.cost_write = write_cost_per_mb * (1 - improv_rate_per_year) ** (_t - 2025)
    tape_on_premise.cost_read = read_cost_per_mb * (1 - improv_rate_per_year) ** (_t - 2025) if read_cost_per_mb > 0 else 0
    tape_on_premise.cost_maint = maint_cost_per_mb_per_year * (1 - improv_rate_per_year) ** (_t - 2025)

    tape_on_premise.n = n
    tape_on_premise.k = k
    tape_on_premise.obj_size = obj_size_mb
    tape_on_premise.start_year = start_year
    tape_on_premise.d = d
    tape_on_premise.device_durability = device_durability_years
    tape_on_premise.max_repls = d // device_durability_years

    return tape_on_premise



def tape_on_premise_updated(
    # --- Base‐Year (2025) Costs ---
    media_cost_per_tb: float           = 6.39,   # $/TB
    hardware_cost_per_tb: float        = 6.86,   # $/TB
    energy_cost_per_tb_per_year: float = 0.05,   # $/TB/yr

    # --- Cost‐Decline Rates ---
    media_decline_per_year: float      = 0.20,   # 20%/yr
    hardware_decline_per_year: float   = 0.00,   # 0%/yr
    energy_decline_per_year: float     = 0.15,   # 15%/yr

    # --- Workloads & Timing ---
    n: int                             = 1000,
    k: int                             = 10,
    obj_size_mb: float                 = 100,
    d: int                             = 100,
    start_year: int                    = 2025,
    device_durability_years: int       = 30,
    model_name: str | None             = None
):
    t = Symbol('t', integer=True)

    # --- Convert to per‐MB and per‐MB/yr ---
    media_per_mb     = media_cost_per_tb / 1e6
    energy_per_mb_yr = energy_cost_per_tb_per_year / 1e6
    # amortize hardware over its lifetime:
    hardware_per_tb_per_year  = hardware_cost_per_tb / device_durability_years
    hardware_per_mb_yr        = hardware_per_tb_per_year / 1e6

    # --- Build model ---
    tape = model3.CostModel(model_name or f"Tape_on_premise_{media_decline_per_year}")

    # --- Cost curves ---
    tape.cost_write = media_per_mb * (1 - media_decline_per_year)**(t - 2025)
    tape.cost_maint = energy_per_mb_yr   * (1 - energy_decline_per_year)**(t - 2025) + hardware_per_mb_yr * (1 - hardware_decline_per_year)**(t - 2025)
    tape.cost_read = 0


    # --- Workload params ---
    tape.n               = n
    tape.k               = k
    tape.obj_size        = obj_size_mb
    tape.start_year      = start_year
    tape.d               = d
    tape.device_durability = device_durability_years
    tape.max_repls       = d // device_durability_years

    return tape

