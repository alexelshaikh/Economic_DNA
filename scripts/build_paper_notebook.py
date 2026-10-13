"""Build the consolidated notebook with nbformat; never modify the source notebooks."""
from pathlib import Path
import ast
import textwrap

import nbformat as nbf

ROOT = Path(__file__).resolve().parents[1]
old = nbf.read(ROOT / "Figures_Sensitivity_Paper.ipynb", as_version=4)
cells = []


def markdown(source):
    cells.append(nbf.v4.new_markdown_cell(textwrap.dedent(source).strip()))


def code(source, **metadata):
    cells.append(nbf.v4.new_code_cell(textwrap.dedent(source).strip(), metadata=metadata))


def figure(title, explanation, source):
    markdown(f"## {title}\n\n{explanation}")
    code("SAVE = True  # False previews only; True exports this figure and its data.\n" + textwrap.dedent(source).strip(), tags=["figure"])


markdown("""
    # Paper figures: main article and sensitivity

    Run the **two setup cells** below once. Then run any figure cell independently.
    Each figure cell has a `SAVE` switch and editable settings at the top.
    Exports go to **`figs/paper_2026/`**, leaving the original figures untouched:
    vector PDF/SVG, 300-dpi PNG, CSV data, caption and scenario JSON.
    To save a preview later, call `save_figure("F1_dna_trajectories")` with its figure ID.
    The final optional cell bundles only the figures you have run into a PDF.

    **Contents:** F1; F2a-c unit costs; F3a-c lifecycle components; F4a-b outlook
    and lifecycle; F5a-c encoding/addressing; F6a-b historical fits; M1 tornado;
    S1-S11 sensitivity figures from `Figures_Sensitivity_Paper.ipynb`.

    **Conventions:** the economic figures use the calculator's 2026 baseline and
    decimal MB/TB. A retention of d years charges exactly d years. Annual decline
    percentages are not continuous exponential coefficients. The legacy figure
    topics are retained, but these are recalculations, not pixel-identical reproductions.
    F4a's target-year guide curves are illustrative interpolations, not fitted forecasts.
    F5 retains the legacy addressing/overhead convention; F6 refits historical data
    without changing the economic model's established DNA fits.

    Tape declines still follow the current configuration (media 20%, hardware 0%,
    energy 15%). The suggested 10%/0%/0% alternative has **not** been silently
    adopted. Edit `BASE = replace(BASE, ...)` below when that decision is made.
    No inflation or PPP adjustment is performed. All future trends are scenarios.

    Requirements: repository Python environment plus Jupyter/IPython and nbformat.
    No website server, downloads, browser, or LaTeX installation is required.
""")
markdown("## Setup 1: baseline, style and output")
setup = old.cells[1].source
start = setup.index("# These are explicit")
end = setup.index("FUTURE_START =")
setup = setup[:start] + '''from economic_dna.assumptions import model_metadata
from economic_dna import simulate_dna_unit_costs, find_crossover_years, load_observed_sequencing_costs
from preamble import calc_index_size
from scipy.optimize import curve_fit

# Uses assumptions.yaml, including corrected cloud units and 2026 price anchors.
BASE = Scenario()
# Example alternative, pending the paper's assumption review:
# BASE = replace(BASE, tape_media_decline_percent=10,
#                tape_hardware_decline_percent=0, tape_energy_decline_percent=0)
SAVE = True
''' + setup[end:]
setup = setup.replace('"paper_sensitivity"', '"paper_2026"')
setup = setup.replace('OUT.mkdir(parents=True, exist_ok=True)', '')
setup = setup.replace('"axes.labelweight": "bold"', '"axes.labelweight": "medium"')
setup = setup.replace('"axes.titleweight": "bold"', '"axes.titleweight": "semibold"')
setup = setup.replace('"legend.fontsize": 10, "lines.linewidth": 2,', '''"legend.fontsize": 10, "lines.linewidth": 2,
    "axes.edgecolor": "#aeb8bf", "axes.linewidth": 0.8, "axes.labelpad": 9,
    "text.color": "#202b33", "axes.labelcolor": "#202b33",
    "xtick.color": "#46535d", "ytick.color": "#46535d",
    "legend.handlelength": 2.8, "legend.handletextpad": 0.8,
    "legend.labelspacing": 0.6, "legend.columnspacing": 1.8,
    "legend.edgecolor": "#dce2e6", "legend.facecolor": "white",
    "legend.fancybox": False, "legend.framealpha": 1,''')
setup += '\ndisplay(pd.DataFrame([asdict(BASE)]).T.rename(columns={0: "Baseline value"}))\nprint(model_metadata()["currency"])\n'
code(setup, tags=["setup"])

markdown("## Setup 2: shared plotting and export helpers\n\nRun once. This cell is collapsed to keep the figure cells easy to scan.")
helpers = old.cells[2].source
tree = ast.parse(helpers)
record_node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "record")
lines = helpers.splitlines(keepends=True)
helpers = ''.join(lines[:record_node.lineno - 1]) + ''.join(lines[record_node.end_lineno:])
begin = helpers.index("# Show the baseline components")
end = helpers.index("# S8/S9 configuration")
helpers = helpers[:begin] + helpers[end:]
helpers = helpers.replace("anchored in 2025", "anchored as in BASE")
helpers = helpers.replace('color="0.90", linewidth=0.6', 'color="#e6eaed", linewidth=0.65')
helpers = helpers.replace('return ax.legend(frameon=False, **kwargs)', '''options = {"frameon": True, "borderpad": 0.65}
    options.update(kwargs)
    result = ax.legend(**options)
    result.get_title().set_fontweight("semibold")
    return result''')
helpers += '''

def save_figure(name):
    """Export a previously rendered figure, its exact input scenarios and data."""
    if name not in FIGURES:
        raise KeyError(f"Run the figure cell first: {name}")
    OUT.mkdir(parents=True, exist_ok=True)
    fig = FIGURES[name]
    for extension in ("pdf", "svg", "png"):
        fig.savefig(OUT / f"{name}.{extension}", bbox_inches="tight", facecolor="white")
    for suffix, frame in TABLES[name].items():
        frame.to_csv(OUT / f"{name}_{suffix}.csv", index=False)
    source_names = ["assumptions.yaml", "economic_dna/scenario.py", "economic_dna/simulation.py",
                    "economic_dna/sensitivity.py", "economic_dna/viability.py", "preamble.py",
                    "Figures_Paper_2026.ipynb"]
    metadata = {
        "figure": name, "caption": CAPTIONS[name], "scenarios": SCENARIOS[name],
        "currency": model_metadata()["currency"],
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_sha256": {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in source_names},
        "versions": {p: importlib.metadata.version(p) for p in ("numpy", "pandas", "matplotlib", "scipy", "seaborn")},
    }
    (OUT / f"{name}.json").write_text(json.dumps(metadata, indent=2, allow_nan=False) + "\\n", encoding="utf-8")
    (OUT / f"{name}_caption.txt").write_text(CAPTIONS[name] + "\\n", encoding="utf-8")
    print(f"Saved {name}: PDF, SVG, PNG, CSV and metadata in {OUT}")

def record(name, fig, caption, tables, scenarios):
    """Keep preview and exports tied to this figure, never to plt's current axes."""
    if name in FIGURES:
        plt.close(FIGURES[name])
    FIGURES[name] = fig
    CAPTIONS[name] = caption
    TABLES[name] = {key: value.copy() for key, value in tables.items()}
    SCENARIOS[name] = {key: asdict(value) for key, value in scenarios.items()}
    if SAVE:
        save_figure(name)
    display(fig)
    plt.close(fig)

def mark_parity(ax, years, dna, peer, color="0.4", label="Parity"):
    valid = np.flatnonzero(np.asarray(dna) <= np.asarray(peer))
    if len(valid):
        year = np.asarray(years)[valid[0]]
        ax.axvline(year, color=color, ls=":", lw=1)
        ax.annotate(f"{label}: {year:g}", (year, np.asarray(peer)[valid[0]]),
                    xytext=(5, 12), textcoords="offset points", fontsize=9, color=color)

def addressing_table(scenario, primer_count=10, primer_size=30,
                     oligo_sizes=(100, 150, 300, 500, 1000, 10000)):
    rows = []
    price_per_base = scenario.dna_synthesis_cost_per_mb / 4e6 * (
        1 - scenario.synthesis_decline_percent / 100) ** (scenario.start_year - scenario.dna_cost_base_year)
    for objects in (100, 1000, 10000, 100000, 1000000, 10000000):
        for length in oligo_sizes:
            info = calc_index_size(n_objects=objects, n_primers=primer_count,
                obj_size_mb=scenario.average_asset_size_mb, oligo_size=length,
                primer_size=primer_size, base_synthesis_cost=price_per_base)
            if info:
                rows.append({"Objects": objects, "Oligo bases": str(length), **info})
    return pd.DataFrame(rows)
'''
PARITY_HELPERS = '''
def s4_synthesis_parity(amazon_decline_percent=10.0, start_year=None):
    """Annual synthesis decline at S4 lifecycle parity, not unit-price parity."""
    scenario = replace(BASE, start_year=FUTURE_START if start_year is None else start_year,
                       amazon_decline_percent=amazon_decline_percent,
                       technologies=("DNA", "Amazon Deep Archive"))
    def difference(rate):
        costs = totals(replace(scenario, synthesis_decline_percent=float(rate)))
        return costs["DNA"] - costs["Amazon Deep Archive"]
    rate = parity_root(difference, 0.0, 99.0)
    return rate

def f2_alternative_declines():
    threshold = s4_synthesis_parity(10.0, start_year=2026)
    if not np.isfinite(threshold):
        return [60.0, 70.0, 80.0]
    # Avoid a nearly duplicate illustrative 60% line within two percentage points.
    return [threshold] + ([60.0] if abs(threshold - 60.0) > 2.0 else []) + [70.0, 80.0]

def s8_sequencing_parity(amazon_decline_percent=10.0, start_year=None):
    """Conditional S8(b) lifecycle parity with accelerated synthesis."""
    scenario = replace(BASE, start_year=FUTURE_START if start_year is None else start_year,
                       synthesis_decline_percent=HEATMAP_SYNTHESIS_DECLINE,
                       amazon_decline_percent=amazon_decline_percent,
                       technologies=("DNA", "Amazon Deep Archive"))
    def difference(rate):
        costs = totals(replace(scenario, sequencing_decline_percent=float(rate)))
        return costs["DNA"] - costs["Amazon Deep Archive"]
    rate = parity_root(difference, 0.0, 99.0)
    return rate
'''
helpers += PARITY_HELPERS
code(helpers, tags=["setup"], jupyter={"source_hidden": True})

figure("F1: DNA synthesis and sequencing trajectories", "Historical-fit trajectories extended from the configured price anchor. USD per logical MB, not per nucleotide.", '''
SHOW_RECORDED_COSTS = True  # False hides historical points; fitted curves stay unchanged.
years = np.arange(1995 if SHOW_RECORDED_COSTS else 2000, 2101)
frame = pd.DataFrame({"year": years})
fig, ax = plt.subplots(figsize=(6.4, 4.8), layout="constrained")
line_styles = {"Synthesis": (PALETTE[0], "-"), "Sequencing": (PALETTE[1], "-")}
for name, price, decline in (
    ("Synthesis", BASE.dna_synthesis_cost_per_mb, BASE.synthesis_decline_percent),
    ("Sequencing", BASE.dna_sequencing_cost_per_mb, BASE.sequencing_decline_percent),
):
    frame[name] = price * (1 - decline / 100) ** (years - BASE.dna_cost_base_year)
    color, dash = line_styles[name]
    ax.plot(years, frame[name], label=f"{name} Cost Trajectory", color=color, ls=dash, lw=1.5)
anchor_year = BASE.dna_cost_base_year
gap_years = [anchor_year, 2050]
gaps = []
for year in dict.fromkeys(gap_years):
    synthesis = BASE.dna_synthesis_cost_per_mb * (1 - BASE.synthesis_decline_percent / 100) ** (year - anchor_year)
    sequencing = BASE.dna_sequencing_cost_per_mb * (1 - BASE.sequencing_decline_percent / 100) ** (year - anchor_year)
    orders = np.log10(synthesis / sequencing)
    ax.vlines(year, sequencing, synthesis, color="black", ls=":", lw=1.5)
    ax.hlines([sequencing, synthesis], year - 1, year + 1, color="black", lw=2)
    label = rf"$\times 10^{{{int(round(orders))}}}$"
    ax.annotate(label, (year, np.sqrt(synthesis * sequencing)),
                xytext=(3, 0), textcoords="offset points", va="center", fontsize=11)
    gaps.append({"year": year, "synthesis_usd_per_mb": synthesis,
                 "sequencing_usd_per_mb": sequencing, "synthesis_to_sequencing_ratio": synthesis / sequencing,
                 "orders_of_magnitude": orders})
observations = pd.DataFrame(columns=["year", "process", "recorded_cost", "recorded_unit", "cost_usd_per_logical_mb", "source"])
if SHOW_RECORDED_COSTS:
    from economic_dna.assumptions import load_assumptions
    dna_units = load_assumptions()["dna"]
    synthesis_observed = pd.DataFrame({
        "year": [1995, 2000, 2003, 2006, 2009, 2010, 2011, 2012, 2013, 2014, 2016, 2017, 2023],
        "recorded_cost": [0.75, 0.30, 0.15, 0.10, 0.08, 0.07, 0.07, 0.06, 0.06, 0.06, 0.05, 0.05, 0.09],
    })
    synthesis_observed["process"] = "Synthesis"
    synthesis_observed["recorded_unit"] = "USD/base"
    synthesis_observed["cost_usd_per_logical_mb"] = synthesis_observed.recorded_cost * dna_units["bases_per_byte"] * dna_units["bytes_per_mb"]
    synthesis_observed["source"] = "Figures.ipynb, Figure 6a historical observations"
    sequencing_observed = load_observed_sequencing_costs().rename(columns={"sequencing_cost_usd_per_mb": "recorded_cost"}).copy()
    sequencing_observed["process"] = "Sequencing"
    sequencing_observed["recorded_unit"] = "USD per million sequenced bases"
    sequencing_observed["cost_usd_per_logical_mb"] = sequencing_observed.recorded_cost * dna_units["bases_per_byte"] * dna_units["bytes_per_mb"] / 1e6
    sequencing_observed["source"] = "NHGRI Sequencing_Cost_Data_Table_May2022.csv"
    observations = pd.concat([synthesis_observed, sequencing_observed], ignore_index=True)
    for process, marker in (("Synthesis", "o"), ("Sequencing", "^")):
        points = observations[observations.process == process]
        ax.scatter(points.year, points.cost_usd_per_logical_mb, marker=marker, s=24,
                   facecolors="none", edgecolors=line_styles[process][0], linewidths=0.9,
                   label=f"{process} recorded costs", zorder=4)
ax.set_xticks(sorted({int(years[0]), anchor_year, 2050, 2075, 2100}))
ax.set_yscale("log")
# Freeze the log limits before adding guides so their lower ends meet the axis.
limits = ax.get_ylim()
for gap in gaps:
    ax.vlines(gap["year"], limits[0], gap["sequencing_usd_per_mb"],
              color="0.55", ls="--", lw=1, zorder=0,
              label="_nolegend_")
ax.set_ylim(limits)
ax.set_title("Cost of DNA Synthesis and Sequencing\\nof 1 MB of Data as a Function of Year",
             fontsize=13, fontweight="bold", color="0.15", pad=6)
ax.grid(False, which="both")
ax.tick_params(axis="both", which="major", labelsize=11, colors="0.15", length=0, pad=7)
ax.tick_params(which="minor", length=0)
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_color("0.15")
    spine.set_linewidth(1.1)
ax.set_xlabel("Year", fontsize=12, fontweight="bold", color="0.15", labelpad=4)
ax.set_ylabel("Cost in Dollars ($)", fontsize=12, fontweight="bold", color="0.15", labelpad=4)
ax.legend(loc="upper right", frameon=True, fancybox=True,
          facecolor="white", edgecolor="0.8", framealpha=0.8, fontsize=10.5,
          handlelength=2, handletextpad=0.8, labelspacing=0.4, borderpad=0.4)
record("F1_dna_trajectories", fig, "Historical-fit DNA unit-cost trajectories; extensions are scenarios, not observations. Capped vertical markers show synthesis / sequencing at the price anchor and in 2050, rounded to the nearest power of ten; exact ratios are exported in the gaps table. "
       + ("Historical points are enabled: synthesis observations from Figures.ipynb Figure 6a and NHGRI sequencing data, converted to USD per logical MB using the model's bases-per-byte convention. Observations are not inflation-adjusted and do not replace the established fits." if SHOW_RECORDED_COSTS else "Historical points are disabled."),
       {"costs": frame, "gaps": pd.DataFrame(gaps), "observations": observations}, {"baseline": BASE})
''')

F2A_SOURCE = r'''
SAVE = True  # False previews only; True exports this figure and its data.
years = np.arange(BASE.start_year, 2101)
annual_declines = np.array([BASE.synthesis_decline_percent, BASE.sequencing_decline_percent, *f2_alternative_declines()])
coefficients = -np.log1p(-annual_declines / 100)
markers = ["o", "^", "*", "d", "x", "+"]
colors = [PALETTE[1], PALETTE[3], *([PALETTE[0]] * 4)]
arrow_offsets = [(0, 0), (-41, 42), (7, 19), (-20, -40), (-52, -26), (-75, -60)]
factor = (1 - BASE.amazon_decline_percent / 100) ** (years - BASE.amazon_price_base_year)
cloud = BASE.amazon_put_usd_per_request / BASE.average_asset_size_mb * factor
frame = pd.DataFrame({"year": years, "Amazon": cloud})
fig, ax = plt.subplots(figsize=(10, 6), layout="constrained")
crossings = []
for index, coefficient in enumerate(coefficients):
    decline = 100 * (1 - np.exp(-coefficient))
    values = BASE.dna_synthesis_cost_per_mb * np.exp(-coefficient * (years - BASE.dna_cost_base_year))
    frame[f"DNA: {decline:.2f}% annual decline"] = values
    qualifier = " (baseline fit)" if index == 0 else " (matching sequencing decline)" if index == 1 else ""
    if index == 2 and np.isfinite(s4_synthesis_parity(10.0, start_year=2026)):
        qualifier = " (S4 2026 lifecycle parity, Amazon 10%)"
    label = f"DNA Synthesis: {decline:.2f}% decline/year{qualifier}"
    ax.plot(years, values, color=colors[index], marker=markers[index],
            markevery=5, markersize=5, lw=1.5, ls="-" if index == 0 else "--", label=label)
    # Log interpolation matches the exponential curves between sampled years.
    candidates = np.flatnonzero(values <= cloud)
    if len(candidates) and candidates[0] > 0:
        i = candidates[0]
        before, after = np.log(values[i - 1] / cloud[i - 1]), np.log(values[i] / cloud[i])
        year = years[i - 1] + before / (before - after) * (years[i] - years[i - 1])
        cost = np.exp(np.interp(year, years, np.log(cloud)))
        ax.plot(year, cost, "o", ms=6, mfc="none", mec="black", mew=0.9, zorder=5)
        ax.annotate(rf"Year $\approx${year:.0f}", (year, cost),
                    xytext=arrow_offsets[index], textcoords="offset points", ha="center",
                    va="bottom" if arrow_offsets[index][1] > 0 else "top", fontsize=10,
                    color="0.15", arrowprops={"arrowstyle": "->", "lw": 0.9, "color": "black"})
        crossings.append({"coefficient": coefficient, "annual_decline_percent": decline,
                          "parity_year": year, "cost_usd_per_mb": cost})
ax.plot(years, cloud, color=PALETTE[2], lw=1.5,
        label=f"Amazon Deep Archive Write Cost: {BASE.amazon_decline_percent:g}% decline/year")
gaps = []
baseline_values = frame[f"DNA: {100 * (1 - np.exp(-coefficients[0])):.2f}% annual decline"]
for i in (0, len(years) - 1):
    year, dna_cost, cloud_cost = years[i], baseline_values.iloc[i], cloud[i]
    orders = np.log10(dna_cost / cloud_cost)
    ax.vlines(year, cloud_cost, dna_cost, color="black", ls=":", lw=1.5)
    ax.hlines([cloud_cost, dna_cost], year - 0.6, year + 0.6, color="black", lw=2)
    ax.annotate(rf"$\times 10^{{{round(orders)}}}$", (year, np.sqrt(dna_cost * cloud_cost)),
                xytext=(3 if i == 0 else -3, 0 if i == 0 else 8), textcoords="offset points",
                ha="left" if i == 0 else "right", va="center", fontsize=11)
    gaps.append({"year": int(year), "synthesis_usd_per_mb": dna_cost,
                 "amazon_write_usd_per_mb": cloud_cost, "orders_of_magnitude": orders})
ax.set_yscale("log")
ax.set_xticks(sorted({BASE.start_year, 2050, 2075, 2100}))
ax.set_title("DNA Synthesis and Amazon Deep Archive Write Costs of 1 MB of Data",
             fontsize=13, fontweight="bold", color="0.15", pad=6)
ax.set_xlabel("Year", fontsize=12, fontweight="bold", color="0.15", labelpad=4)
ax.set_ylabel("Cost in Dollars ($)", fontsize=12, fontweight="bold", color="0.15", labelpad=4)
ax.grid(False, which="both")
ax.tick_params(axis="both", which="major", labelsize=11, colors="0.15", length=0, pad=7)
ax.tick_params(which="minor", length=0)
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_color("0.15")
    spine.set_linewidth(1.1)
ax.legend(loc="lower left", frameon=True, fancybox=True,
          facecolor="white", edgecolor="0.8", framealpha=0.8, fontsize=10.5,
          handlelength=2, handletextpad=0.8, labelspacing=0.4, borderpad=0.4)
record("F2a_unit_costs", fig,
       "DNA synthesis and Amazon write costs per logical MB. The original Figure 2a styling is retained, "
       "with the current price baseline and actual annual percentage declines. Circle/arrow markers "
       "show log-interpolated parity years; capped endpoint gaps are rounded to powers of ten.",
       {"costs": frame, "crossings": pd.DataFrame(crossings), "gaps": pd.DataFrame(gaps)}, {"baseline": BASE})
'''

F2BC_PLOT = r'''
sequencing_threshold = s8_sequencing_parity(10.0, start_year=2026)
has_threshold = np.isfinite(sequencing_threshold)
annual_declines = np.array([BASE.sequencing_decline_percent,
                           *([sequencing_threshold] if has_threshold else []), 60, 70, 80])
coefficients = -np.log1p(-annual_declines / 100)
fig, ax = plt.subplots(figsize=(10, 6), layout="constrained")
if include_storage:
    zoom_ax = ax.inset_axes([0.20, 0.34, 0.40, 0.34], zorder=6)
markers = ["o", "^", "*", "d", "x"]
arrow_offsets = [(0, 27.5), (20, 30), (0, 33), (-5, -45), (-55, -25)]
if include_storage:
    arrow_offsets = [(10, 35), (25, 30), (20, 18), (10, -55), (-20, -25)]
if not has_threshold:
    markers.pop(1)
    arrow_offsets.pop(1)
crossings = []
for index, coefficient in enumerate(coefficients):
    decline = 100 * (1 - np.exp(-coefficient))
    values = anchor_price * np.exp(-coefficient * (years - BASE.dna_cost_base_year))
    frame[f"DNA: {decline:.2f}% annual decline"] = values
    ax.plot(years, values, color=PALETTE[1] if index == 0 else PALETTE[3] if index == 1 and has_threshold else PALETTE[0],
            label=f"DNA Sequencing: {decline:.2f}% decline/year" + (" (baseline fit)" if index == 0 else f" (S8 2026 parity; synthesis {HEATMAP_SYNTHESIS_DECLINE:g}%, Amazon 10%)" if index == 1 and has_threshold else ""),
            ls="-" if index == 0 else "--", lw=1.5, marker=markers[index], markevery=5, markersize=5)
    # Locate parity in log space, matching the exponential cost trajectories.
    candidates = np.flatnonzero(values <= cloud)
    if len(candidates) and candidates[0] > 0:
        i = candidates[0]
        before, after = np.log(values[i - 1] / cloud[i - 1]), np.log(values[i] / cloud[i])
        year = years[i - 1] + before / (before - after) * (years[i] - years[i - 1])
        cost = np.exp(np.interp(year, years, np.log(cloud)))
        ax.plot(year, cost, "o", ms=6, mfc="none", mec="black", mew=0.9, zorder=5)
        if not include_storage:
            ax.annotate(rf"Year $\approx${year:.0f}", (year, cost),
                        xytext=arrow_offsets[index], textcoords="offset points", ha="center",
                        va="bottom" if arrow_offsets[index][1] > 0 else "top", fontsize=10,
                        color="0.15", arrowprops={"arrowstyle": "->", "lw": 0.9, "color": "black"})
        crossings.append({"coefficient": coefficient, "annual_decline_percent": decline,
                          "parity_year": year, "cost_usd_per_mb": cost})
cloud_label = "Read + Allocated Storage" if include_storage else "Read"
ax.plot(years, cloud, color=PALETTE[2], lw=1.5,
        label=f"Amazon Deep Archive {cloud_label} Cost: {BASE.amazon_decline_percent:g}% decline/year")
gaps = []
baseline_values = frame[f"DNA: {100 * (1 - np.exp(-coefficients[0])):.2f}% annual decline"]
for i in (0, len(years) - 1):
    year, dna_cost, cloud_cost = years[i], baseline_values.iloc[i], cloud[i]
    orders = abs(np.log10(dna_cost / cloud_cost))
    ax.vlines(year, min(dna_cost, cloud_cost), max(dna_cost, cloud_cost), color="black", ls=":", lw=1.5)
    ax.hlines([cloud_cost, dna_cost], year - 0.6, year + 0.6, color="black", lw=2)
    ax.annotate(rf"$\times 10^{{{round(orders)}}}$", (year, np.sqrt(dna_cost * cloud_cost)),
                xytext=(-9 if include_storage and i == 0 else -3, 0), textcoords="offset points",
                ha="right", va="center", fontsize=11)
    gaps.append({"year": int(year), "sequencing_usd_per_mb": dna_cost,
                 "amazon_usd_per_mb": cloud_cost, "orders_of_magnitude": orders})
ax.set_yscale("log")
ax.margins(y=0.12)
ax.set_xticks(sorted({BASE.start_year, 2050, 2075, 2100}))
title = ("DNA Sequencing and Amazon Deep Archive Read + Storage Costs\nper Retrieved MB"
         if include_storage else "DNA Sequencing and Amazon Deep Archive Read Costs of 1 MB of Data")
ax.set_title(title, fontsize=13, fontweight="bold", color="0.15", pad=6)
ax.set_xlabel("Year", fontsize=12, fontweight="bold", color="0.15", labelpad=4)
ax.set_ylabel("Cost in Dollars ($)", fontsize=12, fontweight="bold", color="0.15", labelpad=4)
ax.grid(False, which="both")
ax.tick_params(axis="both", which="major", labelsize=11, colors="0.15", length=0, pad=7)
ax.tick_params(which="minor", length=0)
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_color("0.15")
    spine.set_linewidth(1.1)
ax.legend(loc="lower left", frameon=True, fancybox=True,
          facecolor="white", edgecolor="0.8", framealpha=0.8, fontsize=10.5,
          handlelength=2, handletextpad=0.8, labelspacing=0.4, borderpad=0.4)
if include_storage:
    zoom_ax.text(0.98, 0.98, "Crossing detail", transform=zoom_ax.transAxes,
                 ha="right", va="top", fontsize=9, fontweight="bold", color="0.15")
    for line in ax.lines:
        if not line.get_label().startswith("_"):
            zoom_ax.plot(line.get_xdata(), line.get_ydata(), color=line.get_color(),
                         ls=line.get_linestyle(), lw=1.5, marker=line.get_marker(),
                         markevery=1, markersize=7)
    ordered = sorted(crossings, key=lambda point: point["parity_year"])
    if ordered:
        left, right = ordered[0]["parity_year"] - 1.5, ordered[-1]["parity_year"] + 2
        lower = min(point["cost_usd_per_mb"] for point in ordered) / 5
        upper = max(point["cost_usd_per_mb"] for point in ordered) * 6
        zoom_ax.set_xlim(left, right)
        zoom_ax.set_ylim(lower, upper)
        ax.add_patch(mpl.patches.Rectangle((left, lower), right - left, upper - lower,
                     fill=False, edgecolor="0.4", lw=1.2))
        for source_x, detail_x in ((left, 0), (right, 1)):
            connector = mpl.patches.ConnectionPatch(
                xyA=(source_x, lower), coordsA=ax.transData,
                xyB=(detail_x, 1), coordsB=zoom_ax.transAxes,
                color="0.55", lw=1.1, clip_on=False, zorder=0.5)
            connector.set_in_layout(False)
            ax.add_artist(connector)
        offsets = [(-2, -52), (64, 4), (60, -24), (2, 18), (0, -25)]
        for index, point in enumerate(ordered):
            year, cost = point["parity_year"], point["cost_usd_per_mb"]
            offset = offsets[index % len(offsets)]
            zoom_ax.plot(year, cost, "o", ms=7, mfc="none", mec="black", mew=1, zorder=5)
            zoom_ax.annotate(rf'Year $\approx${year:.0f}',
                             (year, cost), xytext=offset, textcoords="offset points", ha="center",
                             va="bottom" if offset[1] > 0 else "top", fontsize=8,
                             arrowprops={"arrowstyle": "->", "color": "black", "lw": 0.9})
    else:
        zoom_ax.text(0.5, 0.5, "No crossings in the displayed timeline",
                     transform=zoom_ax.transAxes, ha="center")
    zoom_ax.set_yscale("log")
    zoom_ax.set_xlabel("Year", fontsize=9, fontweight="bold", color="0.15", labelpad=2)
    zoom_ax.tick_params(axis="both", labelsize=8, colors="0.15", pad=2)
    zoom_ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(4, integer=True))
    zoom_ax.grid(False, which="both")
    for spine in zoom_ax.spines.values():
        spine.set_visible(True)
        spine.set_color("0.15")
        spine.set_linewidth(1.1)
record("FIGURE_NAME", fig,
       f"DNA alternative decline scenarios versus Amazon. S8(b) lifecycle parity is checked for deployment 2026 with synthesis decline {HEATMAP_SYNTHESIS_DECLINE:g}%/year and Amazon decline 10%/year. "
       + ("A conditional sequencing threshold is included. " if has_threshold else "No sequencing threshold exists; the illustrative 60/70/80% scenarios are not lifecycle parity rates. ")
       + "Actual annual declines are shown; storage allocation uses annual retrieved volume where included. Circle/arrow markers show log-interpolated unit-cost parity years; capped gaps show the larger-to-smaller cost ratio rounded to powers of ten.",
       {"costs": frame, "crossings": pd.DataFrame(crossings), "gaps": pd.DataFrame(gaps)}, {"baseline": BASE})
'''

for suffix, kind, maintenance in (("a", "Synthesis", False), ("b", "Sequencing", False), ("c", "Sequencing", True)):
    figure(f"F2{suffix}: {kind.lower()} and Amazon unit costs",
           "Alternative continuous coefficients from the original figure are converted to actual annual percentage declines. First sampled parity years are marked when present. " +
           ("Storage is allocated per retrieved MB using the configured retrieval fraction." if maintenance else "Cloud request costs are allocated using the configured object size."), f'''
years = np.arange(BASE.start_year, 2101)
coefficients = {([0.1671, 0.4791, 0.6, 0.7, 0.8] if kind == "Synthesis" else [0.4791, 0.6, 0.7, 0.8])!r}
kind = {kind!r}
include_storage = {maintenance!r}
factor = (1 - BASE.amazon_decline_percent / 100) ** (years - BASE.amazon_price_base_year)
if kind == "Synthesis":
    cloud_price = BASE.amazon_put_usd_per_request / BASE.average_asset_size_mb
    anchor_price = BASE.dna_synthesis_cost_per_mb
else:
    cloud_price = BASE.amazon_bulk_restore_usd_per_request / BASE.average_asset_size_mb + BASE.amazon_bulk_retrieval_usd_per_mb
    anchor_price = BASE.dna_sequencing_cost_per_mb
    if include_storage:
        if BASE.annual_retrieval_percent <= 0:
            raise ValueError("Per-retrieved-MB storage allocation requires positive retrieval.")
        cloud_price += BASE.amazon_storage_usd_per_mb_month * 12 / (BASE.annual_retrieval_percent / 100)
cloud = cloud_price * factor
frame = pd.DataFrame({{"year": years, "Amazon": cloud}})
fig, ax = plt.subplots(figsize=(10, 5.5), layout="constrained")
for index, coefficient in enumerate(coefficients):
    decline = 100 * (1 - np.exp(-coefficient))
    label = f"DNA: {{decline:.2f}}% annual decline"
    values = anchor_price * np.exp(-coefficient * (years - BASE.dna_cost_base_year))
    frame[label] = values
    ax.plot(years, values, color=PALETTE[index], label=label, ls="-" if index == 0 else "--")
    mark_parity(ax, years, values, cloud, PALETTE[index], "Parity")
ax.plot(years, cloud, color="black", label="Amazon" + (" retrieval + allocated storage" if include_storage else " unit cost"))
ax.set_yscale("log")
tidy(ax, "Calendar year", "USD per " + ("written MB" if kind == "Synthesis" else "retrieved MB"), "{kind} and Amazon costs")
legend(ax, fontsize=9)
record("F2{suffix}_unit_costs", fig, "DNA alternative decline scenarios versus Amazon. Actual annual declines are shown; storage allocation uses annual retrieved volume where included.", {{"costs": frame}}, {{"baseline": BASE}})
''')
    if suffix == "a":
        cells[-1].source = textwrap.dedent(F2A_SOURCE).strip()
    else:
        prefix = cells[-1].source.split("fig, ax = plt.subplots", 1)[0]
        prefix = prefix.replace("coefficients = [0.4791, 0.6, 0.7, 0.8]\n", "")
        cells[-1].source = prefix + F2BC_PLOT.strip().replace("FIGURE_NAME", f"F2{suffix}_unit_costs")

F3A_STYLE = r'''
ax.set_title(
    "Total Storage Costs by Deployment Year\n"
    f"({BASE.number_of_assets:,.0f} objects; {BASE.average_asset_size_mb:,.0f} MB/object; "
    f"{BASE.annual_retrieval_percent:g}% retrieval/year; {BASE.horizon_years} years)",
    fontsize=13, fontweight="bold", color="0.15", pad=8)
ax.set_xlabel("")
ax.set_ylabel("Total Costs in Dollars ($)", fontsize=12, fontweight="bold", color="0.15", labelpad=4)
ax.grid(False, which="both")
ax.tick_params(axis="both", which="major", labelsize=11, colors="0.15", length=0, pad=7)
ax.tick_params(which="minor", length=0)
for label in ax.get_xticklabels():
    label.set_fontweight("bold")
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_color("0.15")
    spine.set_linewidth(1.1)
key = ax.legend(title=r"Deployment Year ($y_0$)", loc="upper right",
                frameon=True, fancybox=True, facecolor="white", edgecolor="0.8",
                framealpha=0.8, fontsize=10.5, title_fontsize=11,
                handlelength=2, handletextpad=0.8, labelspacing=0.4, borderpad=0.5)
key.get_title().set_fontweight("semibold")
'''

for suffix, columns, title in (
    ("a", ["total_cost_usd"], "Total lifecycle cost"),
    ("b", ["write_cost_usd", "maintenance_cost_usd"], "Writing and storage cost"),
    ("c", ["read_cost_usd", "maintenance_cost_usd"], "Retrieval and storage cost"),
):
    figure(f"F3{suffix}: {title.lower()} by deployment year", "Every deployment is an independent archive with the same retention and workload.", f'''
start_years = [BASE.start_year, 2050, 2100]
columns = {columns!r}
scenarios = {{str(year): replace(BASE, start_year=year) for year in start_years}}
frames = []
for year, scenario in scenarios.items():
    data = simulate_scenario(scenario).totals.copy()
    data["Start year"] = year
    data["Cost"] = data[columns].sum(axis=1)
    frames.append(data)
frame = pd.concat(frames, ignore_index=True)
fig, ax = plt.subplots(figsize=(10, 5), layout="constrained")
sns.barplot(data=frame, x="technology", y="Cost", hue="Start year", errorbar=None, ax=ax)
ax.set_yscale("log")
tidy(ax, "Storage technology", "Lifecycle cost (USD)", "{title}")
legend(ax, title="Deployment year")
record("F3{suffix}_lifecycle_components", fig, "{title} for independent archives; " + describe(BASE), {{"costs": frame}}, scenarios)
''')
    if suffix == "a":
        cells[-1].source = cells[-1].source.replace(
            'tidy(ax, "Storage technology", "Lifecycle cost (USD)", "Total lifecycle cost")\nlegend(ax, title="Deployment year")',
            F3A_STYLE.strip())
    elif suffix == "b":
        cells[-1].source = cells[-1].source.replace(
            'tidy(ax, "Storage technology", "Lifecycle cost (USD)", "Writing and storage cost")\nlegend(ax, title="Deployment year")',
            F3A_STYLE.strip().replace("Total Storage Costs by Deployment Year", "Writing and Storage Costs by Deployment Year")
            .replace("Total Costs in Dollars ($)", "Costs in Dollars ($)"))
    elif suffix == "c":
        cells[-1].source = cells[-1].source.replace(
            'tidy(ax, "Storage technology", "Lifecycle cost (USD)", "Retrieval and storage cost")\nlegend(ax, title="Deployment year")',
            F3A_STYLE.strip().replace("Total Storage Costs by Deployment Year", "Retrieval and Storage Costs by Deployment Year")
            .replace("Total Costs in Dollars ($)", "Costs in Dollars ($)"))

F4A_SOURCE = r'''
SAVE = True  # False previews only; True exports this figure and its data.
final_year = 2350
target_years = [2035, 2045]
scenario = replace(BASE, technologies=("DNA", "Amazon Deep Archive"))
frame = simulate_start_years(scenario, final_year)
pivot = frame.pivot(index="start_year", columns="technology", values="total_cost_usd")
fig, ax = plt.subplots(figsize=(10, 6), layout="constrained")
ax.plot(pivot.index, pivot.DNA, color=PALETTE[1], lw=1.5, marker="o", markevery=10, ms=4,
        label=f"DNA: synthesis {BASE.synthesis_decline_percent:.2f}%, sequencing {BASE.sequencing_decline_percent:.2f}% decline/year")
ax.plot(pivot.index, pivot["Amazon Deep Archive"], color=PALETTE[2], lw=1.5,
        label=f"Amazon Deep Archive: {BASE.amazon_decline_percent:g}% tariff decline/year")
guides = []
for index, target in enumerate(target_years):
    if not BASE.start_year < target <= final_year:
        raise ValueError("Target year must follow the baseline and lie within the outlook.")
    years = np.arange(BASE.start_year, target + 1)
    values = np.geomspace(pivot.DNA.iloc[0], pivot.loc[target, "Amazon Deep Archive"], len(years))
    guide_decline = 100 * (1 - (values[-1] / values[0]) ** (1 / (target - BASE.start_year)))
    ax.plot(years, values, ls="--", lw=1.7, color=(PALETTE[0], PALETTE[3])[index % 2],
            label=f"Illustrative total-cost guide: {guide_decline:.2f}%/year; parity in {target}")
    ax.plot(target, values[-1], "o", ms=5, mfc="none", mec="black", mew=1, zorder=5)
    guides.append(pd.DataFrame({"year": years, "target": target, "cost_usd": values,
                               "total_cost_decline_percent": guide_decline}))
crossings = []
ratio = np.log(pivot.DNA.to_numpy() / pivot["Amazon Deep Archive"].to_numpy())
candidates = np.flatnonzero(ratio <= 0)
if len(candidates) and candidates[0] > 0:
    i = candidates[0]
    year = pivot.index[i - 1] + ratio[i - 1] / (ratio[i - 1] - ratio[i])
    cost = np.exp(np.interp(year, pivot.index, np.log(pivot["Amazon Deep Archive"])))
    ax.plot(year, cost, "o", ms=6, mfc="none", mec="black", mew=1, zorder=5)
    ax.annotate(rf"Year $\approx${year:.0f}", (year, cost), xytext=(-25, 35),
                textcoords="offset points", ha="center", fontsize=10,
                arrowprops={"arrowstyle": "->", "lw": 0.9, "color": "black"})
    crossings.append({"parity_year": year, "cost_usd": cost})
ax.set_yscale("log")
ax.set_xticks(sorted({BASE.start_year, 2050, *range(2100, final_year + 1, 50)}))
ax.set_title(f"DNA and Amazon Deep Archive: {BASE.horizon_years}-Year Storage Costs\n"
             "by Deployment Year", fontsize=13, fontweight="bold", color="0.15", pad=8)
ax.set_xlabel(r"Deployment Year ($y_0$)", fontsize=12, fontweight="bold", color="0.15", labelpad=4)
ax.set_ylabel("Total Lifecycle Cost in Dollars ($)", fontsize=12, fontweight="bold", color="0.15", labelpad=4)
ax.grid(False, which="both")
ax.tick_params(axis="both", which="major", labelsize=11, colors="0.15", length=0, pad=7)
ax.tick_params(which="minor", length=0)
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_color("0.15")
    spine.set_linewidth(1.1)
ax.legend(loc="upper right", frameon=True, fancybox=True, facecolor="white",
          edgecolor="0.8", framealpha=0.8, fontsize=9.5, handlelength=2,
          handletextpad=0.8, labelspacing=0.5, borderpad=0.5)
record("F4a_start_year_outlook", fig,
       "Independent archive deployments with fixed workload and retention. DNA legend rates apply to synthesis and sequencing unit prices; Amazon's rate applies to all tariffs. "
       "Dashed guides interpolate total lifecycle costs geometrically to Amazon parity in the labeled target year. Their percentages describe reductions in total cost per year of delayed deployment, not synthesis or sequencing decline rates, fitted forecasts, or solved technology trajectories.",
       {"costs": frame, "guides": pd.concat(guides), "crossings": pd.DataFrame(crossings)}, {"baseline": scenario})
'''

figure("F4a: start-year outlook", "Dashed guides show illustrative total-lifecycle-cost reductions per year of delayed deployment, not solved synthesis decline rates.", '''
final_year = 2350
target_years = [2035, 2045]
scenario = replace(BASE, technologies=("DNA", "Amazon Deep Archive"))
frame = simulate_start_years(scenario, final_year)
pivot = frame.pivot(index="start_year", columns="technology", values="total_cost_usd")
fig, ax = plt.subplots(figsize=(10, 5.5), layout="constrained")
for name in scenario.technologies:
    ax.plot(pivot.index, pivot[name], label=name, color=COLORS[name])
guides = []
for target in target_years:
    if not BASE.start_year < target <= final_year:
        raise ValueError("Target year must follow the baseline and lie within the outlook.")
    years = np.arange(BASE.start_year, target + 1)
    values = np.geomspace(pivot.DNA.iloc[0], pivot.loc[target, "Amazon Deep Archive"], len(years))
    ax.plot(years, values, ls="--", label=f"Illustrative guide to {target} parity")
    guides.append(pd.DataFrame({"year": years, "target": target, "cost_usd": values}))
mark_parity(ax, pivot.index, pivot.DNA, pivot["Amazon Deep Archive"])
ax.set_yscale("log")
tidy(ax, "Archive deployment year", "Lifecycle cost (USD)", "Start-year outlook")
legend(ax, fontsize=9)
record("F4a_start_year_outlook", fig, "Independent archive deployments. Dashed guides are geometric interpolations, not forecasts or solved technology trajectories.", {"costs": frame, "guides": pd.concat(guides)}, {"baseline": scenario})
''')
cells[-1].source = F4A_SOURCE.strip()

F4B_STYLE = r'''
ax.set_title("Cumulative Storage Costs Near Deployment-Year Parity\n"
             f"({BASE.number_of_assets:,.0f} objects; {BASE.average_asset_size_mb:,.0f} MB/object; "
             f"{BASE.annual_retrieval_percent:g}% retrieval/year; {BASE.horizon_years} years)",
             fontsize=13, fontweight="bold", color="0.15", pad=8)
ax.set_xlabel("Calendar Year", fontsize=12, fontweight="bold", color="0.15", labelpad=4)
ax.set_ylabel("Cumulative Total Costs in Dollars ($)", fontsize=12, fontweight="bold", color="0.15", labelpad=4)
ax.grid(False, which="both")
ax.tick_params(axis="both", which="major", labelsize=11, colors="0.15", length=0, pad=7)
ax.tick_params(which="minor", length=0)
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_color("0.15")
    spine.set_linewidth(1.1)
technology_handles = [Line2D([], [], color="0.25", marker=marker, lw=1.5,
                             markersize=6, label=name)
                      for name, marker in (("DNA", "^"), ("Amazon Deep Archive", "s"))]
year_handles = [Line2D([], [], color=PALETTE[index], lw=2, label=year)
                for index, year in enumerate(scenarios)]
options = dict(loc="upper left", frameon=True, fancybox=True, facecolor="white",
               edgecolor="0.8", framealpha=0.8, fontsize=10.5, title_fontsize=11,
               handlelength=2, handletextpad=0.8, labelspacing=0.5, borderpad=0.5)
technology_key = ax.legend(handles=technology_handles, title="Storage Technology",
                           bbox_to_anchor=(1.02, 1), **options)
technology_key.get_title().set_fontweight("semibold")
ax.add_artist(technology_key)
year_key = ax.legend(handles=year_handles, title=r"Deployment Year ($y_0$)",
                     bbox_to_anchor=(1.02, 0.65), **options)
year_key.get_title().set_fontweight("semibold")
'''

figure("F4b: lifecycle costs around parity", "Selects three deployments around the first start-year parity; it does not hard-code obsolete crossover years.", '''
scenario = replace(BASE, technologies=("DNA", "Amazon Deep Archive"))
crossing = find_crossover_years(simulate_start_years(scenario, 2500))["Amazon Deep Archive"]
center = crossing if crossing is not None else FUTURE_START
start_years = sorted(set(max(BASE.start_year, min(2500, center + offset)) for offset in (-1, 0, 1)))
scenarios = {str(year): replace(scenario, start_year=year) for year in start_years}
frames = []
fig, ax = plt.subplots(figsize=(10, 5.5), layout="constrained")
for index, (year, current) in enumerate(scenarios.items()):
    data = simulate_scenario(current).yearly.copy()
    data["deployment"] = year
    for name, part in data.groupby("technology", sort=False):
        ax.plot(part.year, part.total_cost_usd.cumsum(), color=PALETTE[index],
                ls="-", lw=1.5, marker="^" if name == "DNA" else "s",
                markersize=6 if name == "DNA" else 4, markevery=5, label=f"{year}: {name}")
    frames.append(data)
ax.set_yscale("log")
tidy(ax, "Calendar year", "Cumulative lifecycle cost (USD)", "Lifecycle costs near start-year parity")
legend(ax, fontsize=9)
record("F4b_lifecycle_near_parity", fig, "Three independent archive deployments near first sampled start-year parity; fallback to the configured future year if none exists.", {"annual_costs": pd.concat(frames)}, scenarios)
''')
cells[-1].source = cells[-1].source.replace(
    'tidy(ax, "Calendar year", "Cumulative lifecycle cost (USD)", "Lifecycle costs near start-year parity")\nlegend(ax, fontsize=9)',
    F4B_STYLE.strip())

F5A_STYLE = r'''
ax.set_title(f"DNA Storage Costs by Encoding Density ({BASE.horizon_years} Years)\n"
             f"({BASE.number_of_assets:,.0f} objects; {BASE.average_asset_size_mb:,.0f} MB/object; "
             f"{BASE.annual_retrieval_percent:g}% retrieval/year)",
             fontsize=13, fontweight="bold", color="0.15", pad=8)
ax.set_xlabel("Encoding Density (bits/base)", fontsize=12, fontweight="bold", color="0.15", labelpad=4)
ax.set_ylabel("Total Cost in Dollars ($)", fontsize=12, fontweight="bold", color="0.15", labelpad=4)
ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5), numticks=12))
ax.yaxis.set_major_formatter(mpl.ticker.LogFormatterSciNotation(
    base=10, labelOnlyBase=False, minor_thresholds=(np.inf, np.inf)))
ax.grid(False, which="both")
ax.tick_params(axis="both", which="major", labelsize=11, colors="0.15", length=0, pad=7)
ax.tick_params(which="minor", length=0)
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_color("0.15")
    spine.set_linewidth(1.1)
key = ax.legend(title="Address Overhead\n(relative to payload)", loc="upper right",
                frameon=True, fancybox=True, facecolor="white", edgecolor="0.8",
                framealpha=0.8, fontsize=10.5, title_fontsize=11,
                handlelength=2, handletextpad=0.8, labelspacing=0.4, borderpad=0.5)
key.get_title().set_fontweight("semibold")
'''

figure("F5a: encoding density and address overhead", "Retains the legacy multiplicative overhead convention. The physical range is limited to at most 2 bits/base; prices are scaled relative to the model's 4 bases/byte convention.", '''
bits_per_base = np.linspace(0.1, 2.0, 40)
overhead_percentages = [0, 10, 20, 30, 40, 50]
copies = 1
baseline_cost = totals(replace(BASE, technologies=("DNA",)))["DNA"]
rows = []
fig, ax = plt.subplots(figsize=(9, 5.5), layout="constrained")
for overhead in overhead_percentages:
    multiplier = 2 / bits_per_base * copies * (1 + overhead / 100)
    costs = baseline_cost * multiplier
    ax.plot(bits_per_base, costs, label=f"{overhead}%", ls="-", lw=1.5,
            marker="o", markersize=4, markevery=2)
    rows.extend({"bits_per_base": b, "overhead_percent": overhead, "copies": copies, "cost_usd": c} for b, c in zip(bits_per_base, costs))
ax.set_yscale("log")
tidy(ax, "Encoding density (bits/base)", "DNA lifecycle cost (USD)", "Encoding density and additive address overhead")
legend(ax, title="Overhead relative to payload")
record("F5a_encoding_overhead", fig, "DNA costs scaled by (2 / bits per base) * copies * (1 + address overhead / 100), following the legacy convention; overhead is relative to payload, not a fraction of total oligo length.", {"costs": pd.DataFrame(rows)}, {"baseline": BASE})
''')
cells[-1].source = cells[-1].source.replace(
    'tidy(ax, "Encoding density (bits/base)", "DNA lifecycle cost (USD)", "Encoding density and additive address overhead")\nlegend(ax, title="Overhead relative to payload")',
    F5A_STYLE.strip())

F5B_STYLE = r'''
ax.set_title("Address Index Length by Number of Data Objects\n"
             f"({primer_count} primers; {primer_bases} bases/primer; {BASE.average_asset_size_mb:,.0f} MB/object)",
             fontsize=13, fontweight="bold", color="0.15", pad=8)
ax.set_xlabel(r"Number of Data Objects ($n$)", fontsize=12, fontweight="bold", color="0.15", labelpad=4)
ax.set_ylabel("Index Length (bases)", fontsize=12, fontweight="bold", color="0.15", labelpad=4)
ax.grid(False, which="both")
ax.tick_params(axis="both", which="major", labelsize=11, colors="0.15", length=0, pad=7)
ax.tick_params(which="minor", length=0)
ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
for label in ax.get_xticklabels():
    label.set_fontweight("bold")
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_color("0.15")
    spine.set_linewidth(1.1)
handles, labels = ax.get_legend_handles_labels()
key = ax.legend(handles, [f"{int(label):,}" for label in labels], title="Oligo Length (bases)",
                loc="upper left", bbox_to_anchor=(1.02, 1), borderaxespad=0,
                frameon=True, fancybox=True, facecolor="white", edgecolor="0.8",
                framealpha=0.8, fontsize=10.5, title_fontsize=11,
                handlelength=2, handletextpad=0.8, labelspacing=0.4, borderpad=0.5)
key.get_title().set_fontweight("semibold")
'''

for suffix, column, ylabel in (("b", "Index Size", "Index length (bases)"), ("c", "Total Synthesis Cost", "Initial synthesis cost (USD)")):
    figure(f"F5{suffix}: {ylabel.lower()}", "Uses the existing oligo-addressing helper on only the plotted combinations. Object count varies, so total archive size varies too; primers/indexes retain the original convention.", f'''
primer_count = 10
primer_bases = 30
oligo_sizes = (100, 150, 300, 500, 1000, 10000)
frame = addressing_table(BASE, primer_count, primer_bases, oligo_sizes)
fig, ax = plt.subplots(figsize=(10, 5.5), layout="constrained")
sns.barplot(data=frame, x="Objects", y={column!r}, hue="Oligo bases", errorbar=None, ax=ax)
{('ax.set_yscale("log")' if suffix == 'c' else '')}
ax.set_xticks(range(6), [r"$10^{{%d}}$" % power for power in range(2, 8)])
tidy(ax, "Number of objects", {ylabel!r}, f"{{primer_count}} primers; {{primer_bases}} primer bases; {{BASE.average_asset_size_mb:g}} MB/object")
legend(ax, title="Oligo length (bases)", fontsize=9)
record("F5{suffix}_addressing", fig, "Legacy oligo-level addressing calculation at the configured deployment-year synthesis price; archive size varies with object count.", {{"addressing": frame}}, {{"baseline": BASE}})
''')
    if suffix == "b":
        source = cells[-1].source
        cells[-1].source = source[:source.index("tidy(ax,")] + F5B_STYLE.strip() + "\n" + source[source.index("record("):]
    elif suffix == "c":
        source = cells[-1].source
        style = F5B_STYLE.strip().replace("Address Index Length by Number of Data Objects", "Initial DNA Synthesis Costs by Number of Data Objects")
        style = style.replace("Index Length (bases)", "Initial Synthesis Cost in Dollars ($)")
        style = style.replace("ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))\n", "")
        cells[-1].source = source[:source.index("tidy(ax,")] + style + "\n" + source[source.index("record("):]

F6_STYLE = r'''
ax.set_title("FIT_TITLE", fontsize=13, fontweight="bold", color="0.15", pad=8)
ax.set_xlabel("Year", fontsize=12, fontweight="bold", color="0.15", labelpad=4)
ax.set_ylabel("FIT_YLABEL", fontsize=12, fontweight="bold", color="0.15", labelpad=4)
ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
ax.grid(False, which="both")
ax.tick_params(axis="both", which="major", labelsize=11, colors="0.15", length=0, pad=7)
ax.tick_params(which="minor", length=0)
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_color("0.15")
    spine.set_linewidth(1.1)
ax.legend(loc="upper right", frameon=True, fancybox=True, facecolor="white",
          edgecolor="0.8", framealpha=0.8, fontsize=10.5,
          handlelength=2, handletextpad=0.8, labelspacing=0.4, borderpad=0.5)
mantissa, exponent = f"{rmse:.2e}".split("e")
metrics = (rf"$R^2 = {r2:.3f}$" + " (original scale)\n"
           + rf"$\mathrm{{RMSE}} = {mantissa}\times10^{{{int(exponent)}}}$" + " FIT_UNIT")
ax.text(0.03, 0.03, metrics, transform=ax.transAxes, ha="left", va="bottom", fontsize=10,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "0.8",
              "linewidth": 0.8, "alpha": 0.92})
'''

figure("F6a: historical synthesis fit", "Refits the historical observations embedded in Figures.ipynb. This diagnostic does not replace the economic model's established synthesis fit.", '''
years = np.array([1995, 2000, 2003, 2006, 2009, 2010, 2011, 2012, 2013, 2014, 2016, 2017, 2023])
costs = np.array([0.75, 0.30, 0.15, 0.10, 0.08, 0.07, 0.07, 0.06, 0.06, 0.06, 0.05, 0.05, 0.09])
def fit_function(year, amplitude, coefficient):
    return amplitude * np.exp(-coefficient * (year - 1995))
parameters, _ = curve_fit(fit_function, years, costs, p0=[1, 0.1])
predicted = fit_function(years, *parameters)
rmse = np.sqrt(np.mean((costs - predicted)**2))
r2 = 1 - np.sum((costs - predicted)**2) / np.sum((costs - costs.mean())**2)
fine_years = np.linspace(years.min(), BASE.start_year, 250)
fig, ax = plt.subplots(figsize=(9, 5), layout="constrained")
ax.scatter(years, costs, label="Synthesis Cost Data", color=PALETTE[0], s=32, zorder=3)
ax.plot(fine_years, fit_function(fine_years, *parameters), label="Exponential Fit", color=PALETTE[1], lw=1.5)
ax.set_yscale("log")
tidy(ax, "Year", "Synthesis price (USD/base)", "Historical synthesis fit")
legend(ax)
caption = f"Legacy synthesis observations; nonlinear least-squares refit in original units. R2={r2:.4f}; RMSE={rmse:.6g} USD/base. This refit is not substituted into the economic model."
record("F6a_synthesis_fit", fig, caption, {"observations": pd.DataFrame({"year": years, "cost_usd_per_base": costs, "fitted": predicted}), "fit": pd.DataFrame({"year": fine_years, "cost_usd_per_base": fit_function(fine_years, *parameters)})}, {"baseline": BASE})
''')
cells[-1].source = cells[-1].source.replace(
    'tidy(ax, "Year", "Synthesis price (USD/base)", "Historical synthesis fit")\nlegend(ax)',
    F6_STYLE.strip().replace("FIT_TITLE", "Fitting DNA Synthesis Cost Over Time")
    .replace("FIT_YLABEL", "Cost per Base in Dollars ($)").replace("FIT_UNIT", "USD/base")
    + '\nax.yaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5), numticks=12))\n'
    + 'ax.yaxis.set_major_formatter(mpl.ticker.LogFormatterSciNotation(base=10, labelOnlyBase=False, minor_thresholds=(np.inf, np.inf)))')

figure("F6b: historical sequencing fit", "Uses the local NHGRI CSV. Units are USD per million sequenced bases (Mb), not USD per logical MB of stored data.", '''
frame = load_observed_sequencing_costs().copy()
years = frame.year.to_numpy()
costs = frame.sequencing_cost_usd_per_mb.to_numpy()
origin = years.min()
slope, intercept = np.polyfit(years - origin, np.log(costs), 1)
fine_years = np.linspace(years.min(), years.max(), 250)
fit_costs = np.exp(intercept + slope * (fine_years - origin))
predicted = np.exp(intercept + slope * (years - origin))
rmse = np.sqrt(np.mean((costs - predicted)**2))
r2 = 1 - np.sum((costs - predicted)**2) / np.sum((costs - costs.mean())**2)
fig, ax = plt.subplots(figsize=(9, 5), layout="constrained")
ax.scatter(years, costs, label="NHGRI Sequencing Cost Data", color=PALETTE[0], s=24, zorder=3)
ax.plot(fine_years, fit_costs, label="Log-linear Exponential Fit", color=PALETTE[1], lw=1.5)
ax.set_yscale("log")
tidy(ax, "Year", "Sequencing cost (USD per million bases)", "Historical sequencing fit")
legend(ax)
caption = f"NHGRI local observations, nominal USD at reporting dates; log-linear refit. Original-scale R2={r2:.4f}; RMSE={rmse:.6g} USD per million bases. Model fit unchanged."
record("F6b_sequencing_fit", fig, caption, {"observations": frame, "fit": pd.DataFrame({"year": fine_years, "cost_usd_per_million_bases": fit_costs})}, {"baseline": BASE})
''')
cells[-1].source = cells[-1].source.replace(
    'tidy(ax, "Year", "Sequencing cost (USD per million bases)", "Historical sequencing fit")\nlegend(ax)',
    F6_STYLE.strip().replace("FIT_TITLE", "Fitting DNA Sequencing Cost Over Time")
    .replace("FIT_YLABEL", "Cost per Million Bases in Dollars ($)").replace("FIT_UNIT", "USD/million bases"))

def style_s1_s2_s3(source):
    if 'record("S1_' in source:
        title = "DNA Cost Advantage by Synthesis Price"
        legend_title = "Comparison Storage Model"
        source = source.replace('figsize=(12.4, 4.9)', 'figsize=(12.8, 5.7)')
        source = source.replace('legend(axes[1], loc="lower left")', '')
        source = source.replace('Cost advantage of DNA (USD)', 'DNA Cost Advantage (USD)')
        source = source.replace('Synthesis price in {BASE.start_year}', 'Synthesis Price in {BASE.start_year}')
        source = source.replace('Retrieval: {scenario.annual_retrieval_percent:g}%/year',
                                'Data retrieval rate: {scenario.annual_retrieval_percent:g}%/year')
        source = source.replace('ax.text(0.03, 0.45, "No parity at a nonnegative synthesis price",',
                                'ax.text(0.03, 0.45, "No parity at a nonnegative synthesis price",')
    elif 'record("S2_' in source:
        title = "Break-Even Synthesis Price by Data Retrieval Rate"
        legend_title = "Comparison Storage Model\nCurves ending before 100% retrieval: at higher retrieval rates, even free synthesis cannot achieve parity."
        source = source.replace('"future": replace(BASE, start_year=FUTURE_START)}',
                                '"future": replace(BASE, start_year=FUTURE_START), "2100": replace(BASE, start_year=2100)}')
        source = source.replace('plt.subplots(1, 2, figsize=(12.4, 4.9)',
                                'plt.subplots(1, 3, figsize=(17, 5.7)')
        source = source.replace('Start years are {BASE.start_year} and {FUTURE_START}',
                                'Deployment years are {BASE.start_year}, {FUTURE_START}, and 2100')
        source = source.replace('lifecycle model; absent values indicate that no nonnegative synthesis price can achieve parity ',
                                'lifecycle model. When a curve ends before the 100% retrieval boundary, higher retrieval rates cannot achieve parity even with free synthesis ')
        source = source.replace("because sequencing alone exceeds the comparator's total. Both axes are logarithmic.",
                                "because sequencing alone exceeds the comparator's total. Curves reaching the right boundary are only limited by the plotted range. Both axes are logarithmic.")
        source = source.replace('fig.legend(handles, labels, loc="outside lower center", ncol=3, frameon=False)', '')
        source = source.replace('Annual retrieval (% of archive/year)', 'Data Retrieval Rate (% of Archive/Year)')
        source = source.replace('Break-even synthesis price (USD/MB)', 'Break-Even Synthesis Price (USD/MB)')
        source = source.replace('Archive starts in {scenario.start_year}', 'Deployment Year: {scenario.start_year}')
        source = source.replace('    ax.set_xscale("log")', '''    ax.set_xscale("log")
    ax.set_xticks([0.00001, 0.001, 0.1, 1, 10, 100])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:.5f}".rstrip("0").rstrip(".") + "%"))
    ax.xaxis.set_minor_formatter(NullFormatter())''')
        source = source.replace('    ax.text(0.03, 0.04, "Gaps indicate no nonnegative price can give parity",\n'
                                '            transform=ax.transAxes, fontsize=9)',
                                '')
        source = source.replace('ax.axvline(BASE.annual_retrieval_percent, color="0.4", ls=":", lw=1)',
                                'ax.axvline(BASE.annual_retrieval_percent, color="0.4", ls=":", lw=1, label="Baseline retrieval rate")')
    else:
        title = "DNA Cost Competitiveness by Deployment Year"
        legend_title = "Annual Synthesis Cost Decline"
        source = source.replace('figsize=(14, 4.8)', 'figsize=(14, 5.7)')
        source = source.replace('legend(axes[0], loc="upper right")', '')
        source = source.replace('Archive start year $y_0$', 'Deployment Year $y_0$')
        source = source.replace('DNA / comparator lifecycle cost', 'DNA / Comparator Lifecycle Cost')
        source = source.replace('ax.axhline(1, color="0.25", lw=1)',
                                'ax.axhline(1, color="0.25", ls=":", lw=1.2, label="Cost parity")')
    styling = '''
for ax in axes:
    ax.grid(False, which="both")
    ax.tick_params(axis="both", which="major", labelsize=10.5, colors="0.15", length=0, pad=7)
    ax.tick_params(which="minor", length=0)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("0.15")
        spine.set_linewidth(1.1)
    ax.set_xlabel(ax.get_xlabel(), fontsize=12, fontweight="bold", color="0.15", labelpad=7)
    ax.set_ylabel(ax.get_ylabel(), fontsize=12, fontweight="bold", color="0.15", labelpad=7)
    panel_title = ax.get_title(loc="left")
    ax.set_title("", loc="left")
    ax.set_title(panel_title, fontsize=12, fontweight="bold", color="0.15", pad=14)
    for line, marker in zip(ax.lines[:3], ("o", "s", "^")):
        line.set_linewidth(2)
        line.set_marker(marker)
        line.set_markersize(4.5)
        line.set_markerfacecolor("none")
        line.set_markeredgewidth(1.1)
        line.set_markevery(max(1, len(line.get_xdata()) // 12))
    for collection in ax.collections:
        collection.set_sizes([55])
        collection.set_facecolor("none")
        collection.set_linewidth(1.6)
handles, labels = axes[0].get_legend_handles_labels()
fig.suptitle(TITLE, fontsize=14, fontweight="bold", color="0.15")
shared_legend = fig.legend(handles, labels, title=LEGEND_TITLE,
           loc="outside lower center", ncol=len(handles), frameon=True,
           fancybox=True, facecolor="white", edgecolor="0.8", framealpha=0.8,
           fontsize=10.5, handlelength=2.8, handletextpad=0.8, borderpad=0.6)
shared_legend.get_title().set_fontweight("bold")
'''.replace('LEGEND_TITLE', repr(legend_title)).replace('TITLE', repr(title))
    start = source.index('record(')
    source = source[:start] + styling + source[start:]
    if 'record("S2_' in source:
        from paper_retrieval_figures import add_parameter_titles
        source = add_parameter_titles(source, "S2")
    elif 'record("S3_' in source:
        from paper_retrieval_figures import refine_s3_titles
        source = refine_s3_titles(source)
    return source


def style_s4(source):
    source = source.replace('figsize=(8.5, 5.6)', 'figsize=(9.5, 6.2)')
    source = source.replace('linewidths=1.5', 'linewidths=1.8')
    source = source.replace('marker="*", s=180', 'marker="*", s=210')
    styling = '''
ax.grid(False, which="both")
ax.set_xlabel(ax.get_xlabel(), fontsize=12, fontweight="bold", color="0.15", labelpad=8)
ax.set_ylabel(ax.get_ylabel(), fontsize=12, fontweight="bold", color="0.15", labelpad=8)
plot_title = ax.get_title(loc="left")
ax.set_title("", loc="left")
ax.set_title(plot_title, fontsize=14, fontweight="bold", color="0.15", pad=16)
ax.tick_params(which="major", labelsize=10.5, colors="0.15", length=0, pad=7)
ax.tick_params(which="minor", length=0)
ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}%"))
ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:g}%"))
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_color("0.15")
    spine.set_linewidth(1.1)
colorbar = mesh.colorbar
colorbar.ax.tick_params(labelsize=10.5, colors="0.15", length=0, pad=7)
colorbar.ax.set_ylabel(colorbar.ax.get_ylabel(), fontsize=12, fontweight="bold", color="0.15", labelpad=10)
colorbar.outline.set_edgecolor("0.15")
colorbar.outline.set_linewidth(1.1)
for text in ax.texts:
    text.set_color("0.15")
    box = text.get_bbox_patch()
    if box is not None:
        box.set_edgecolor("0.8")
        box.set_linewidth(0.8)
        box.set_alpha(0.9)
'''
    start = source.index('s4_threshold =')
    return source[:start] + styling + source[start:]


def style_m1(source):
    source = source.replace('ax.set_xticks([0.1, 0.5, 1, 2, 10])',
        'ax.set_xticks([0.001, 0.1, 1, 10, 1000])')
    source = source.replace('"future": replace(BASE, start_year=FUTURE_START)}',
        '"future": replace(BASE, start_year=FUTURE_START), "2100": replace(BASE, start_year=2100)}')
    source = source.replace('plt.subplots(1, 2, figsize=(12.8, 5.7)',
        'plt.subplots(1, 3, figsize=(16, 5.7)')
    source = source.replace('starting in {BASE.start_year} ', 'starting in {BASE.start_year}, ')
    source = source.replace('and {FUTURE_START}. Workload:', '{FUTURE_START}, and 2100. Workload:')
    source = source.replace('    frame = dna_cost_sensitivity(scenario, use_present_value=False)',
        '    frame = dna_cost_sensitivity(scenario, use_present_value=False)\n'
        '    frame = frame[frame.parameter != "Discount rate"].copy()')
    start = source.index('magnitude = m1.assign(')
    end = source.index('fig, axes =', start)
    source = source[:start] + '''order = [
    "Synthesis cost", "Synthesis decline rate",
    "Sequencing cost", "Sequencing decline rate",
    "DNA durability", "Archive size", "Average asset size", "Annual retrieval",
]
''' + source[end:]
    source = source.replace('p.replace("Average asset size", "Average object size")',
        'p.replace("Average asset size", "Average object size").replace("Annual retrieval", "Data retrieval rate")')
    source = source.replace('Inputs vary by +/-50%, with a 0\u20135% ',
        'Inputs vary by +/-50%. Data retrieval rate is the percentage of archive volume retrieved per year. ')
    source = source.replace('discount-rate fallback; all costs here are undiscounted.',
        'Discount rate is omitted because costs here are undiscounted. Inputs appear in a fixed order across deployment years.')
    source = source.replace('color="0.72", lw=7', 'color="0.72", lw=5')
    source = source.replace('marker="o", s=45', 'marker="o", s=55')
    source = source.replace('marker="s", s=38', 'marker="s", s=48')
    source = source.replace('f"{x:g}"', 'rf"${x:g}\\times$"')
    start = source.index('    tidy(ax, "DNA cost / baseline DNA cost"')
    end = source.index("handles =", start)
    source = source[:start] + '''    ax.set_title(f"({chr(97 + index)}) Deployment Year: {scenario.start_year}",
                 fontsize=12, fontweight="bold", color="0.15", pad=25)
    ax.text(0.5, 1.025, f"Baseline lifecycle cost: ${cost:,.0f}",
            transform=ax.transAxes, fontsize=10, ha="center", color="0.3")
    ax.set_xlabel("DNA Cost Relative to Baseline", fontsize=12, fontweight="bold", color="0.15", labelpad=7)
    ax.grid(False, which="both")
    ax.tick_params(axis="both", which="major", labelsize=10.5, colors="0.15", length=0, pad=7)
    ax.tick_params(which="minor", length=0)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("0.15")
        spine.set_linewidth(1.1)
''' + source[end:]
    source = source.replace('axes[1].legend(handles=handles, loc="lower right", frameon=False)', '''handles.append(Line2D([], [], color="0.25", ls="--", lw=1, label="Baseline cost"))
fig.suptitle("Sensitivity of DNA Lifecycle Costs", fontsize=14, fontweight="bold", color="0.15")
fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=True,
           fancybox=True, facecolor="white", edgecolor="0.8", framealpha=0.8,
           fontsize=10.5, handlelength=2, handletextpad=0.8, borderpad=0.5)''')
    return source

markdown("# Sensitivity figures\n\nThe following figure code is carried over from the sensitivity notebook, using the shared 2026 baseline and export controls above. Each cell is independent of the other figure cells.")
for index in range(4, 27, 2):
    heading = old.cells[index - 1].source.replace("anchored in 2025", "anchored as in BASE").replace("paper's start year", "configured baseline start year")
    markdown(heading)
    code("SAVE = True  # False previews only; True exports this figure and its data.\n" + old.cells[index].source, tags=["figure"])
    if 'record("M1_dna_cost_sensitivity"' in cells[-1].source:
        cells[-1].source = style_m1(cells[-1].source)
    if any(f'record("S{i}_' in cells[-1].source for i in (1, 2, 3)):
        cells[-1].source = style_s1_s2_s3(cells[-1].source)
    if 'record("S4_joint_decline_sensitivity"' in cells[-1].source:
        cells[-1].source = cells[-1].source.replace(
            'record("S4_joint_decline_sensitivity"',
            's4_threshold = s4_synthesis_parity(10.0)\n'
            'print(f"S4 lifecycle parity at Amazon decline 10%/year: synthesis decline {s4_threshold:.6f}%/year")\n'
            'record("S4_joint_decline_sensitivity"',
        ).replace('{"grid": s4}', '{"grid": s4, "parity_at_amazon_10": pd.DataFrame([{"amazon_decline_percent": 10.0, "synthesis_decline_percent": s4_threshold}])}')

def add_2026_sensitivity_variants(target_cells):
    for stem, title in (("S4_joint_decline_sensitivity", "S4 (2026): synthesis and Amazon declines"),
                        ("S8_sequencing_decline_sensitivity", "S8 (2026): sequencing and Amazon declines")):
        index = next(i for i, cell in enumerate(target_cells)
                     if "figure" in cell.metadata.get("tags", []) and f'record("{stem}"' in cell.source)
        source = target_cells[index].source.replace("FUTURE_START", "2026").replace(stem, stem + "_2026")
        if stem.startswith("S4"):
            source = source.replace("s4_synthesis_parity(10.0)", "s4_synthesis_parity(10.0, start_year=2026)")
            source = source.replace('ax.scatter(BASE.synthesis_decline_percent',
                'if z.min() > 0:\n'
                '    ax.text(0.5, 0.83, "No cost parity for a 2026 deployment\\nInitial synthesis is paid at the 2026 anchor price",\n'
                '            transform=ax.transAxes, ha="center", va="center", fontsize=11,\n'
                '            bbox={"facecolor": "white", "alpha": 0.92, "edgecolor": "none"})\n'
                'ax.scatter(BASE.synthesis_decline_percent')
            source = source.replace('print(f"S4 lifecycle parity at Amazon decline 10%/year: synthesis decline {s4_threshold:.6f}%/year")',
                'print(f"Synthesis parity: {s4_threshold:.6f}%/year" if np.isfinite(s4_threshold) else "No synthesis-decline parity for a 2026 deployment at Amazon decline 10%/year.")')
        else:
            source += '\nsequencing_threshold_2026 = s8_sequencing_parity(10.0, start_year=2026)\nprint(f"Sequencing parity: {sequencing_threshold_2026:.6f}%/year" if np.isfinite(sequencing_threshold_2026) else "No sequencing-decline parity for a 2026 deployment at Amazon decline 10%/year.")'
        note = ("Same sweep and price anchors as the 2050 figure, but the archive starts in 2026. "
                "With the baseline 100-year retention, the initial synthesis bill is paid immediately at the 2026 price. "
                "Later declines cannot reduce that bill, and there is no lifecycle break-even rate at Amazon decline 10%/year. "
                "A missing parity contour is intentional, not a plotting error.")
        target_cells[index + 1:index + 1] = [nbf.v4.new_markdown_cell("## " + title + "\n\n" + note),
            nbf.v4.new_code_cell(source, metadata={"tags": ["figure"]})]

add_2026_sensitivity_variants(cells)
for cell in cells:
    if cell.cell_type == "code" and 'record("S4_' in cell.source:
        cell.source = style_s4(cell.source)
from paper_retrieval_figures import add_retrieval_figures
add_retrieval_figures(cells, nbf)
from paper_joint_price_figure import add_joint_price_figure
add_joint_price_figure(cells, nbf)
from paper_deployment_parity_figure import add_deployment_parity_figure
add_deployment_parity_figure(cells, nbf)
from paper_workload_figures import add_workload_figures
add_workload_figures(cells, nbf)
from paper_s10_style import style_s10, style_s11
from paper_f1_diagnostics import add_f1_diagnostics
for cell in cells:
    if cell.cell_type == "code" and 'record("F1_' in cell.source:
        cell.source = add_f1_diagnostics(cell.source)
    if cell.cell_type == "code" and 'record("S10_' in cell.source:
        cell.source = style_s10(cell.source)
    if cell.cell_type == "code" and 'record("S11_' in cell.source:
        cell.source = style_s11(cell.source)
from paper_s10b_figure import add_s10b
from paper_separate_parity_style import compact_parity_legend, simplify_s2_title
for cell in cells:
    if cell.cell_type != "code":
        continue
    if 'record("S10_' in cell.source or 'record("S11_' in cell.source:
        cell.source = compact_parity_legend(cell.source)
    if 'record("S2_retrieval_' in cell.source:
        cell.source = simplify_s2_title(cell.source)
add_s10b(cells, nbf)
from paper_long_retention_figure import add_long_retention_figure
add_long_retention_figure(cells, nbf)

from paper_combined_parity_figure import SOURCE as combined_parity_source
combined_index = next(i for i, cell in enumerate(cells)
                      if cell.cell_type == "code" and 'record("S11_' in cell.source)
cells[combined_index + 1:combined_index + 1] = [
    nbf.v4.new_markdown_cell("## S10 + S11: combined parity figure\n\nTwo panels with one shared legend."),
    nbf.v4.new_code_cell(combined_parity_source, metadata={"tags": ["figure"]}),
]

markdown("## Optional: export all rendered figures\n\nIncludes only figures run in this session, in notebook order. Rerun this cell after changing figures; it does not require running every figure.")
code('''
if not FIGURES:
    raise RuntimeError("Run at least one figure cell first.")
OUT.mkdir(parents=True, exist_ok=True)
notebook = json.loads((ROOT / "Figures_Paper_2026.ipynb").read_text(encoding="utf-8"))
ordered_names = []
for cell in notebook["cells"]:
    source = "".join(cell["source"])
    for name in FIGURES:
        if f'record("{name}"' in source and name not in ordered_names:
            ordered_names.append(name)
with PdfPages(OUT / "Paper_Figures.pdf") as pdf:
    for name in ordered_names:
        save_figure(name)
        pdf.savefig(FIGURES[name], bbox_inches="tight", facecolor="white")
(OUT / "figure_captions.md").write_text(
    "# Figure captions\\n\\n" + "\\n\\n".join(f"## {name}\\n\\n{CAPTIONS[name]}" for name in ordered_names), encoding="utf-8")
display(pd.DataFrame({"Figure": ordered_names, "PDF": [f"{name}.pdf" for name in ordered_names]}))
print(f"Combined PDF: {OUT / 'Paper_Figures.pdf'}")
''', tags=["export"])

from paper_results_export import add_results_export
add_results_export(cells, nbf)

notebook = nbf.v4.new_notebook(cells=cells, metadata={
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python", "version": "3.13"},
})
for index, cell in enumerate(notebook.cells):
    cell.id = f"paper-2026-{index:03d}"
nbf.validate(notebook)
nbf.write(notebook, ROOT / "Figures_Paper_2026.ipynb")
print(f"Created Figures_Paper_2026.ipynb: {len(cells)} cells, 40 figures")
