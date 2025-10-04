import os
from tkinter import getint
from itertools import product
from tqdm import tqdm
import requests
from jupyter_server.utils import ensure_async
from sympy import Symbol, exp
from io import BytesIO
import storage_service
import pandas as pd
import seaborn as sns
import numpy as np
import xlrd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from scipy.optimize import curve_fit
import math
import matplotlib.ticker as mtick
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from matplotlib.patches import Rectangle, ConnectionPatch
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
import matplotlib as mpl
from datetime import timedelta


plt.rcParams['mathtext.fontset'] = 'cm'  # Use Computer Modern for math only
sns.set_theme(
    style    = 'white',   # {'darkgrid', 'whitegrid', 'dark', 'white', 'ticks'}
    context  = 'notebook',        # {'paper', 'notebook', 'talk', 'poster'}
    font_scale = 1.1,         # scale up all font sizes by 10%
)
def save(name):
    os.makedirs("figs", exist_ok=True)
    plt.savefig(f"figs/{name}.svg", transparent=True, bbox_inches="tight")
    plt.savefig(f"figs/{name}.pdf", transparent=True, bbox_inches="tight")
def find_intersection(line1, line2):
    # get raw data
    x1, y1 = line1.get_xdata(), line1.get_ydata()
    x2, y2 = line2.get_xdata(), line2.get_ydata()

    # choose a common x‐grid (union of both sets—for simplicity)
    x_all = np.unique(np.concatenate([x1, x2]))
    # interpolate each curve onto that grid
    y1_i = np.interp(x_all, x1, y1)
    y2_i = np.interp(x_all, x2, y2)

    # compute difference and look for sign changes
    diff = y1_i - y2_i
    sign_changes = np.where(np.diff(np.sign(diff)) != 0)[0]
    intersections = []
    for idx in sign_changes:
        # linear solve between x_all[idx] and x_all[idx+1]
        x0, x1_ = x_all[idx], x_all[idx+1]
        y0d, y1d = diff[idx], diff[idx+1]
        t = -y0d / (y1d - y0d)           # fraction along the segment
        xi = x0 + t * (x1_ - x0)
        yi = np.interp(xi, x_all, y1_i)  # or y2_i, they're (nearly) equal
        intersections.append((xi, yi))
    return intersections
def make_exponential(x1, y1, x2, y2):
    """
    Returns a function f(x)=A*exp(B*x) satisfying
      f(x1)=y1 and f(x2)=y2.
    """
    # Solve for B: e^{B*(x2-x1)} = y2/y1  ==>  B = ln(y2/y1) / (x2-x1)
    B = math.log(y2 / y1) / (x2 - x1)
    # Then A = y1 / exp(B*x1)
    A = y1 * math.exp(-B * x1)
    return lambda x: A * math.exp(B * x)
def obj_size_by_bit_rate(obj_size, bit_rate, copies=1, index_size=30):
    return obj_size * 2 / bit_rate * copies * (1 + index_size / 100)
def replacement_years(start_year, L, d):
    current_year = start_year
    repl_years = []
    while current_year < start_year + d:
        lifetime = L(current_year) if callable(L) else L
        current_year += lifetime
        if current_year < start_year + d:
            repl_years.append(current_year)

    return repl_years
def synthesis_cost_updated(year,
                           drop_rate: float = 0.1671,
                           anchor_year: float = 2000):
    # 1) your true historical baseline:
    orig_rate              = 0.1671     # fixed, from the original fit
    baseline_year          = 2000
    baseline_cost_per_base = 0.3193

    # 2) roll that baseline to chosen anchor using the ORIGINAL rate
    cost_per_base_at_anchor = (
        baseline_cost_per_base
        * np.exp(-orig_rate * (anchor_year - baseline_year))
    )

    # 3) from that anchor, apply whichever drop_rate you like:
    years_since_anchor = year - anchor_year
    cost = (
        cost_per_base_at_anchor
        * np.exp(-drop_rate * years_since_anchor)
        * 1_000_000
        * 4
    )
    return cost

def calc_index_size(n_objects: int,
                    n_primers: int,
                    oligo_size: int,
                    primer_size: int,
                    obj_size_mb: int = None,
                    obj_size_bytes: int = None,
                    base_synthesis_cost: float = 0.004897040777491498):
    # total data‐encoding bases needed (4 bases per byte)
    if obj_size_bytes is None and obj_size_mb is None:
        return {}

    if obj_size_bytes is None:
        obj_size_bytes = obj_size_mb * 1000_000

    total_data_bases = n_objects * obj_size_bytes * 4

    # case 0: no index at all
    payload_size = oligo_size - primer_size
    if payload_size <= 0:
        return {}  # can't even fit one data base

    total_num_oligos = math.ceil(total_data_bases / payload_size)
    num_oligos_to_address = math.ceil(total_num_oligos / n_primers)
    if n_primers >= num_oligos_to_address:
        synthesis_cost_per_oligo = base_synthesis_cost * oligo_size
        total_synthesis_cost = synthesis_cost_per_oligo * total_num_oligos
        return {
            "Index Size": 0,
            "Payload Size": payload_size,
            "Number of Oligos": total_num_oligos,
            "Bytes per Oligo": payload_size / 4,
            "Synthesis Cost Per Oligo": synthesis_cost_per_oligo,
            "Total Synthesis Cost": total_synthesis_cost,
            "Oligo Size": oligo_size,
            "Primer Size": primer_size,
            "Total Data Size (MBs)": (obj_size_bytes * n_objects) / 1e6,
            "Primer Size %": primer_size / oligo_size,
            "Index Size %": 0.0,
            "Payload Size %": payload_size / oligo_size,
            "Primer+Index Sizes %": primer_size / oligo_size,
        }

    # case 1+: try adding index bases
    for index_size in range(1, oligo_size - primer_size):
        payload_size = oligo_size - primer_size - index_size
        if payload_size <= 0:
            break

        total_num_oligos = math.ceil(total_data_bases / payload_size)
        num_oligos_to_address = math.ceil(total_num_oligos / n_primers)
        capacity   = 4**index_size
        if capacity >= num_oligos_to_address:
            synthesis_cost_per_oligo = base_synthesis_cost * oligo_size
            total_synthesis_cost = synthesis_cost_per_oligo * total_num_oligos
            return {
                "Index Size": index_size,
                "Payload Size": payload_size,
                "Number of Oligos": total_num_oligos,
                "Bytes per Oligo": payload_size / 4,
                "Synthesis Cost Per Oligo": synthesis_cost_per_oligo,
                "Total Synthesis Cost": total_synthesis_cost,
                "Oligo Size": oligo_size,
                "Primer Size": primer_size,
                "Total Data Size (MBs)": (obj_size_bytes * n_objects) / 1e6,
                "Primer Size %": primer_size / oligo_size,
                "Index Size %": index_size / oligo_size,
                "Payload Size %": payload_size / oligo_size,
                "Primer+Index Sizes %": (primer_size + index_size) / oligo_size,
            }

    # no workable index size found
    return {}

def get_df_fig2():
    years = list(range(2000, 2500))
    anchor_years = list(range(2000, 2031, 1))
    rows = []
    for anchor_year in anchor_years:
        for year in years:
            for rate in [0.1671, 0.4791, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]:
                rows.append({
                    'Anchor Year': anchor_year,
                    'Trend': rate,
                    'isRealTrend': True if rate == 0.1671 else False,
                    'Year': year,
                    'Cost': 'Synthesis',
                    'Value': synthesis_cost_updated(year, drop_rate=rate, anchor_year=anchor_year),
                })
                rows.append({
                    'Anchor Year': anchor_year,
                    'Trend': rate,
                    'isRealTrend': True if rate == 0.4791 else False,
                    'Year': year,
                    'Cost': 'Sequencing',
                    'Value': sequencing_cost_updated(year, drop_rate=rate, anchor_year=anchor_year),
                })

    return pd.DataFrame(rows)

def sequencing_cost_updated(year, drop_rate=0.4791, anchor_year=2000):
    # 1) original “baseline” point:
    baseline_year         = 2011.5
    baseline_cost_per_run = 35.036    # USD per run at baseline_year
    multiplier            = 4         # 4 DNA abses per byte

    # 2) roll that baseline back or forward to anchor using the original rate
    orig_rate = 0.4791
    cost_at_anchor = (
        baseline_cost_per_run
        * np.exp(-orig_rate * (anchor_year - baseline_year))
        * multiplier
    )

    # 3) from the anchor, decay/grow to the target year with the chosen drop_rate
    return cost_at_anchor * np.exp(-drop_rate * (year - anchor_year))

def get_sequencing_data_available(do_prints=False):
    # URL of the XLS file
    url = "https://www.genome.gov/sites/default/files/media/files/2023-05/Sequencing_Cost_Data_Table_May2022.xls"

    # Create 'data' folder if it doesn't exist
    os.makedirs("data", exist_ok=True)

    # Path for CSV
    csv_path = os.path.join("data", "Sequencing_Cost_Data_Table_May2022.csv")

    if os.path.exists(csv_path):
        if do_prints:
            print("Sequencing data file already exists - was not downloaded again!")
        return pd.read_csv(csv_path, parse_dates=['Date'])
    # Download file into memory
    response = requests.get(url)
    response.raise_for_status()

    # Read the XLS content directly from memory
    df = pd.read_excel(BytesIO(response.content), parse_dates=['Date'])

    # Save it as CSV
    df.to_csv(csv_path, index=False)

    if do_prints:
        print(f"Sequencing data file was was refreshed/downloaded!")
    return df

def get_df_fig3(k=10):
    n = 1000
    d = 100
    start_year = Symbol('start_year', integer=True)
    obj_size = 1000

    dna_model = storage_service.dna(n=n, k=k, start_year=start_year, obj_size_mb=obj_size, d=d, device_durability=1000)
    amazon_model = storage_service.amazon(n=n, k=k, start_year=start_year, obj_size_mb=obj_size, d=d)
    azure_model = storage_service.azure(n=n, k=k, start_year=start_year, obj_size_mb=obj_size, d=d)
    tape_model = storage_service.tape_on_premise_updated(n=n, k=k, start_year=start_year, obj_size_mb=obj_size, d=d, device_durability_years=30, hardware_decline_per_year=0.1)


    years = [2025, 2050, 2075]

    dna_costs   = [float(dna_model.eval().subs(start_year, year).doit()) for year in years]
    amzn_costs  = [float(amazon_model.eval().subs(start_year, year).doit()) for year in years]
    azure_costs = [float(azure_model.eval().subs(start_year, year).doit()) for year in years]
    tape_costs  = [float(tape_model.eval().subs(start_year, year).doit()) for year in years]


    dna_read_costs = [float(dna_model._get_read_sum().subs(start_year, year).doit()) for year in years]
    dna_write_costs = [float(dna_model._get_initial_write_cost().subs(start_year, year).doit()) for year in years]
    dna_run_costs = [float(dna_model._get_maintenance_sum().subs(start_year, year).doit()) for year in years]

    amzn_read_costs = [float(amazon_model._get_read_sum().subs(start_year, year).doit()) for year in years]
    amzn_write_costs = [float(amazon_model._get_initial_write_cost().subs(start_year, year).doit()) for year in years]
    amzn_run_costs = [float(amazon_model._get_maintenance_sum().subs(start_year, year).doit()) for year in years]

    azure_read_costs = [float(azure_model._get_read_sum().subs(start_year, year).doit()) for year in years]
    azure_write_costs = [float(azure_model._get_initial_write_cost().subs(start_year, year).doit()) for year in years]
    azure_run_costs = [float(azure_model._get_maintenance_sum().subs(start_year, year).doit()) for year in years]

    tape_read_costs = [float(tape_model._get_read_sum().subs(start_year, year).doit()) for year in years]
    tape_write_costs = [float(tape_model._get_initial_write_cost().subs(start_year, year).doit()) for year in years]
    tape_run_costs = [float(tape_model._get_maintenance_sum().subs(start_year, year).doit()) for year in years]

    df_wide = pd.DataFrame({
        "Year":                 years,
        "DNA":                  dna_costs,
        "Amazon Deep Archive" : amzn_costs,
        "Azure Blob Archive" :  azure_costs,
        "Tape On-premise":      tape_costs,
    })

    df_long = df_wide.melt(
        id_vars="Year",
        var_name="Storage",
        value_name="Cost"
    )

    def compute_reads(row, percent=False):
        year_index = years.index(int(row['Year']))
        if row['Storage'] == 'DNA':
            return dna_read_costs[year_index] / row['Cost'] if percent else dna_read_costs[year_index]
        elif row['Storage'] == 'Amazon Deep Archive':
            return amzn_read_costs[year_index] / row['Cost'] if percent else amzn_read_costs[year_index]
        elif row['Storage'] == 'Azure Blob Archive':
            return azure_read_costs[year_index] / row['Cost'] if percent else azure_read_costs[year_index]

        return tape_read_costs[year_index] / row['Cost'] if percent else tape_read_costs[year_index]

    def compute_writes(row, percent=False):
        year_index = years.index(row['Year'])
        if row['Storage'] == 'DNA':
            return dna_write_costs[year_index] / row['Cost'] if percent else dna_write_costs[year_index]
        elif row['Storage'] == 'Amazon Deep Archive':
            return amzn_write_costs[year_index] / row['Cost'] if percent else amzn_write_costs[year_index]
        elif row['Storage'] == 'Azure Blob Archive':
            return azure_write_costs[year_index] / row['Cost'] if percent else azure_write_costs[year_index]

        return tape_write_costs[year_index] / row['Cost'] if percent else tape_write_costs[year_index]

    def compute_run(row, percent=False):
        year_index = years.index(row['Year'])
        if row['Storage'] == 'DNA':
            return dna_run_costs[year_index] / row['Cost'] if percent else dna_run_costs[year_index]
        elif row['Storage'] == 'Amazon Deep Archive':
            return amzn_run_costs[year_index] / row['Cost'] if percent else amzn_run_costs[year_index]
        elif row['Storage'] == 'Azure Blob Archive':
            return azure_run_costs[year_index] / row['Cost'] if percent else azure_run_costs[year_index]

        return tape_run_costs[year_index] / row['Cost'] if percent else tape_run_costs[year_index]

    df_long['Cost Read'] = df_long.apply(compute_reads, axis=1)
    df_long['Cost Write'] = df_long.apply(compute_writes, axis=1)
    df_long['Cost Read+Run'] = df_long.apply(lambda row: compute_reads(row, percent=False) + compute_run(row, percent=False), axis=1)
    df_long['Cost Write+Run'] = df_long.apply(lambda row: compute_writes(row, percent=False) + compute_run(row, percent=False), axis=1)

    df_long['Cost Read %'] = df_long.apply(lambda row: compute_reads(row, percent=True), axis=1)
    df_long['Cost Write %'] = df_long.apply(lambda row: compute_writes(row, percent=True), axis=1)

    df_long['Cost Read+Run %'] = df_long.apply(lambda row: compute_reads(row, percent=True) + compute_run(row, percent=True), axis=1)
    df_long['Cost Write+Run %'] = df_long.apply(lambda row: compute_writes(row, percent=True) + compute_run(row, percent=True), axis=1)
    return df_long

def get_df_supp_fig8():
    n = 1_000
    obj_size = 1000
    start_years = [2025, 2050, 2075, 2100, 2125, 2250, 2275, 2300]
    d = 10_000
    t = Symbol('t', integer=True)

    amazon_model = storage_service.amazon(n=n, obj_size_mb=obj_size)
    dna_model = storage_service.dna(n=n, obj_size_mb=obj_size)

    dna_model.cost_maint = 0 * t

    rows = []

    tape_map = {}
    dna_map = {}

    for k in [10, 100, 1000]:
        for start_year in start_years:
            tape_write = 0
            tape_read = 0
            tape_maint = 0

            dna_write = 0
            dna_read = 0
            dna_maint = 0

            amazon_model.start_year = start_year
            dna_model.start_year = start_year

            tape_replacement_years = replacement_years(start_year, amazon_model.device_durability, d)
            dna_replacement_years = replacement_years(start_year, dna_model.device_durability, d)

            for current_year in range(start_year, start_year + d + 1):

                # --- TAPE REPLACEMENT WRITE ---
                if current_year in tape_replacement_years:
                    key = ('Write', current_year, k)
                    tape_write_delta = tape_map.get(key)
                    if tape_write_delta is None:
                        tape_write_delta = float(amazon_model.cost_write.subs(t, current_year).doit()) * n * obj_size
                        tape_map[key] = tape_write_delta
                    tape_write += tape_write_delta
                else:
                    tape_write_delta = 0

                # --- TAPE READ ---
                key = ('Read', current_year, k)
                tape_read_delta = tape_map.get(key)
                if tape_read_delta is None:
                    tape_read_delta = float(amazon_model.cost_read.subs(t, current_year).doit()) * k * obj_size
                    tape_map[key] = tape_read_delta
                tape_read += tape_read_delta

                # --- TAPE MAINTENANCE ---
                key = ('Maint', current_year, k)
                tape_maint_delta = tape_map.get(key)
                if tape_maint_delta is None:
                    tape_maint_delta = float(amazon_model.cost_maint.subs(t, current_year).doit()) * n * obj_size
                    tape_map[key] = tape_maint_delta
                tape_maint += tape_maint_delta

                # --- DNA REPLACEMENT WRITE ---
                if current_year in dna_replacement_years:
                    key = ('Write', current_year)
                    dna_write_delta = dna_map.get(key)
                    if dna_write_delta is None:
                        dna_write_delta = float(dna_model.cost_write.subs(t, current_year).doit()) * n * obj_size
                        dna_map[key] = dna_write_delta
                    dna_write += dna_write_delta
                else:
                    dna_write_delta = 0

                # --- DNA READ ---
                key = ('Read', current_year, k)
                dna_read_delta = dna_map.get(key)
                if dna_read_delta is None:
                    dna_read_delta = float(dna_model.cost_read.subs(t, current_year).doit()) * k * obj_size
                    dna_map[key] = dna_read_delta
                dna_read += dna_read_delta

                # --- DNA MAINTENANCE ---
                key = ('Maint', current_year, k)
                dna_maint_delta = dna_map.get(key)
                if dna_maint_delta is None:
                    dna_maint_delta = float(dna_model.cost_maint.subs(t, current_year).doit()) * n * obj_size
                    dna_map[key] = dna_maint_delta
                dna_maint += dna_maint_delta

                rows.extend([
                    {
                        "Current Year": current_year,
                        "Start Year": start_year,
                        "Cumulative Cost": tape_write,
                        "Delta Cost": tape_write_delta,
                        "Cost Type": "Write",
                        "Storage Tech": "Tape On-premise",
                        "k": k
                    },
                    {
                        "Current Year": current_year,
                        "Start Year": start_year,
                        "Cumulative Cost": tape_read,
                        "Delta Cost": tape_read_delta,
                        "Cost Type": "Read",
                        "Storage Tech": "Tape On-premise",
                        "k": k
                    },
                    {
                        "Current Year": current_year,
                        "Start Year": start_year,
                        "Cumulative Cost": tape_maint,
                        "Delta Cost": tape_maint_delta,
                        "Cost Type": "Maint",
                        "Storage Tech": "Tape On-premise",
                        "k": k
                    },
                    {
                        "Current Year": current_year,
                        "Start Year": start_year,
                        "Cumulative Cost": tape_write + tape_read + tape_maint,
                        "Delta Cost": tape_write_delta + tape_read_delta + tape_maint_delta,
                        "Cost Type": "Total",
                        "Storage Tech": "Tape On-premise",
                        "k": k
                    }
                ])
                rows.extend([
                    {
                        "Current Year": current_year,
                        "Start Year": start_year,
                        "Cumulative Cost": dna_write,
                        "Delta Cost": dna_write_delta,
                        "Cost Type": "Write",
                        "Storage Tech": "DNA",
                        "k": k
                    },
                    {
                        "Current Year": current_year,
                        "Start Year": start_year,
                        "Cumulative Cost": dna_read,
                        "Delta Cost": dna_read_delta,
                        "Cost Type": "Read",
                        "Storage Tech": "DNA",
                        "k": k
                    },
                    {
                        "Current Year": current_year,
                        "Start Year": start_year,
                        "Cumulative Cost": dna_maint,
                        "Delta Cost": dna_maint_delta,
                        "Cost Type": "Maint",
                        "Storage Tech": "DNA",
                        "k": k
                    },
                    {
                        "Current Year": current_year,
                        "Start Year": start_year,
                        "Cumulative Cost": dna_write + dna_read + dna_maint,
                        "Delta Cost": dna_write_delta + dna_read_delta + dna_maint_delta,
                        "Cost Type": "Total",
                        "Storage Tech": "DNA",
                        "k": k
                    }
                ])

    return pd.DataFrame(rows)

def get_df_fig4b():
    n = 1_000
    k = 10
    obj_size = 1000
    start_years = [2025, 2050, 2075, 2100, 2125, 2250, 2275, 2300] + list(range(2300, 2351))
    d = 100
    t = Symbol('t', integer=True)

    amazon_model = storage_service.amazon(n=n, k=k, obj_size_mb=obj_size)
    dna_model = storage_service.dna(n=n, k=k, obj_size_mb=obj_size)

    dna_model.cost_maint = 0 * t


    rows = []

    tape_map = {}
    dna_map = {}

    for start_year in start_years:
        tape_write = 0
        tape_read = 0
        tape_maint = 0

        dna_write = 0
        dna_read = 0
        dna_maint = 0

        amazon_model.start_year = start_year
        dna_model.start_year = start_year

        tape_replacement_years = replacement_years(start_year, amazon_model.device_durability, d)
        dna_replacement_years = replacement_years(start_year, dna_model.device_durability, d)

        for current_year in range(start_year, start_year + d + 1):

            # --- TAPE WRITE AND REPLACEMENT ---
            if current_year == start_year or current_year in tape_replacement_years:
                key = ('Write', current_year)
                tape_write_delta = tape_map.get(key)
                if tape_write_delta is None:
                    tape_write_delta = float(amazon_model.cost_write.subs(t, current_year).doit()) * n * obj_size
                    tape_map[key] = tape_write_delta
                tape_write += tape_write_delta
            else:
                tape_write_delta = 0

            # --- TAPE READ ---
            key = ('Read', current_year)
            tape_read_delta = tape_map.get(key)
            if tape_read_delta is None:
                tape_read_delta = float(amazon_model.cost_read.subs(t, current_year).doit()) * k * obj_size
                tape_map[key] = tape_read_delta
            tape_read += tape_read_delta

            # --- TAPE MAINTENANCE ---
            key = ('Maint', current_year)
            tape_maint_delta = tape_map.get(key)
            if tape_maint_delta is None:
                tape_maint_delta = float(amazon_model.cost_maint.subs(t, current_year).doit()) * n * obj_size
                tape_map[key] = tape_maint_delta
            tape_maint += tape_maint_delta

            # --- DNA WRITE ---
            if current_year == start_year or current_year in dna_replacement_years:
                key = ('Write', current_year)
                dna_write_delta = dna_map.get(key)
                if dna_write_delta is None:
                    dna_write_delta = float(dna_model.cost_write.subs(t, current_year).doit()) * n * obj_size
                    dna_map[key] = dna_write_delta
                dna_write += dna_write_delta
            else:
                dna_write_delta = 0

            # --- DNA READ ---
            key = ('Read', current_year)
            dna_read_delta = dna_map.get(key)
            if dna_read_delta is None:
                dna_read_delta = float(dna_model.cost_read.subs(t, current_year).doit()) * k * obj_size
                dna_map[key] = dna_read_delta
            dna_read += dna_read_delta

            # --- DNA MAINTENANCE ---
            key = ('Maint', current_year)
            dna_maint_delta = dna_map.get(key)
            if dna_maint_delta is None:
                dna_maint_delta = float(dna_model.cost_maint.subs(t, current_year).doit()) * n * obj_size
                dna_map[key] = dna_maint_delta
            dna_maint += dna_maint_delta

            rows.extend([
                {
                    "Current Year": current_year,
                    "Start Year": start_year,
                    "Cumulative Cost": tape_write,
                    "Delta Cost": tape_write_delta,
                    "Cost Type": "Write",
                    "Storage Tech": "Tape On-premise"
                },
                {
                    "Current Year": current_year,
                    "Start Year": start_year,
                    "Cumulative Cost": tape_read,
                    "Delta Cost": tape_read_delta,
                    "Cost Type": "Read",
                    "Storage Tech": "Tape On-premise"
                },
                {
                    "Current Year": current_year,
                    "Start Year": start_year,
                    "Cumulative Cost": tape_maint,
                    "Delta Cost": tape_maint_delta,
                    "Cost Type": "Maint",

                    "Storage Tech": "Tape On-premise"
                },
                {
                    "Current Year": current_year,
                    "Start Year": start_year,
                    "Cumulative Cost": tape_write + tape_read + tape_maint,
                    "Delta Cost": tape_write_delta + tape_read_delta + tape_maint_delta,
                    "Cost Type": "Total",
                    "Storage Tech": "Tape On-premise"
                }
            ])
            rows.extend([
                {
                    "Current Year": current_year,
                    "Start Year": start_year,
                    "Cumulative Cost": dna_write,
                    "Delta Cost": dna_write_delta,
                    "Cost Type": "Write",
                    "Storage Tech": "DNA"
                },
                {
                    "Current Year": current_year,
                    "Start Year": start_year,
                    "Cumulative Cost": dna_read,
                    "Delta Cost": dna_read_delta,
                    "Cost Type": "Read",
                    "Storage Tech": "DNA"
                },
                {
                    "Current Year": current_year,
                    "Start Year": start_year,
                    "Cumulative Cost": dna_maint,
                    "Delta Cost": dna_maint_delta,
                    "Cost Type": "Maint",
                    "Storage Tech": "DNA"
                },
                {
                    "Current Year": current_year,
                    "Start Year": start_year,
                    "Cumulative Cost": dna_write + dna_read + dna_maint,
                    "Delta Cost": dna_write_delta + dna_read_delta + dna_maint_delta,
                    "Cost Type": "Total",
                    "Storage Tech": "DNA"
                }
            ])

    return pd.DataFrame(rows)

def get_df_fig5a():
    n = 1_000
    k = 10
    obj_size = 1000
    d = 100
    start_year = 2025
    bit_rates = [b / 10 for b in range(1, 41, 1)]
    copies = [1]
    index_sizes = [0, 10, 20, 30, 40, 50]
    start_years = [start_year]
    ds = [100, 1000, 10_000]
    rows = []
    for current_year in start_years:
        for br in bit_rates:
            for red in copies:
                for index_size in index_sizes:
                    for d in ds:
                        rows.append({
                            "Start Year": current_year,
                            "Bit Rate": br,
                            "Redundancy": red,
                            "d": d,
                            "Index Size": index_size,
                            "Total Cost": float(
                                storage_service.dna(n=n,
                                                    k=k,
                                                    obj_size_mb=obj_size_by_bit_rate(obj_size, br, red, index_size),
                                                    d=d,
                                                    start_year=start_year
                                                    ).eval().doit()
                            )
                        })
    return pd.DataFrame(rows)

def get_df_fig5bc():
    param_grid = {
        "n_primers":      [1, 2, 5, 10, 20, 30, 50, 100, 1000, 10_000, 100_000],
        "n_objects":      [100, 1000, 10_000, 100_000, 1_000_000, 10_000_000],
        "obj_size_mb":    [1, 100, 1_000, 10_000, 100_000, 1_000_000, 10_000_000],
        "oligo_size":     [100, 150, 200, 250, 300, 400, 500, 1000, 5000, 10_000],
        "primer_size":    [20, 30, 40, 50, 60],
    }

    rows = []

    for (n_primers, n_objects, obj_size_mb, oligo_size, primer_size) in tqdm(product(*param_grid.values()), disable=True):

        infos = calc_index_size(
            n_objects=int(n_objects),
            n_primers=int(n_primers),
            obj_size_mb=int(obj_size_mb),
            oligo_size=int(oligo_size),
            primer_size=int(primer_size),
        )
        if not infos:
            continue  # skip impossible configs

        rows.append({
            "Number of Objects": n_objects,
            "Object Size (MBs)": obj_size_mb,
            "Number of Primers": n_primers,
            "Number of Oligos": infos["Number of Oligos"],
            "Oligo Size": oligo_size,
            "Primer Size": primer_size,
            "Index Size": infos["Index Size"],
            "Payload Size": infos["Payload Size"],
            "Total Data Size (MBs)": infos["Total Data Size (MBs)"],
            "Synthesis Cost Per Oligo": infos["Synthesis Cost Per Oligo"],
            "Total Synthesis Cost": infos["Total Synthesis Cost"],
            "Bytes per Oligo": infos["Bytes per Oligo"],
            "Primer Size %": infos["Primer Size %"],
            "Index Size %": infos["Index Size %"],
            "Payload Size %": infos["Payload Size %"],
            "Primer+Index Sizes %": infos["Primer+Index Sizes %"],
        })

    return pd.DataFrame(rows)