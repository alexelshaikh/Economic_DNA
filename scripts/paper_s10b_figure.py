"""S10 with synthesis and sequencing decline roles exchanged."""

try:
    from .paper_s10_style import STYLE
except ImportError:
    from paper_s10_style import STYLE


SOURCE = '''SAVE = True  # False previews only; True exports this figure and its data.
SEQUENCING_RATES = (50.0, 55.0, 60.0, 65.0)
s10b_x = np.linspace(0, 80, 321)
s10b_scenarios = {
    f"sequencing_{rate:g}": replace(BASE, start_year=FUTURE_START,
        sequencing_decline_percent=float(rate), technologies=("DNA", "Amazon Deep Archive"))
    for rate in (*SEQUENCING_RATES, float(BASE.sequencing_decline_percent))
}
rows = []
for key, scenario in s10b_scenarios.items():
    # Scale public-model write costs, retaining its workload and replacement years.
    flat = simulate_scenario(replace(scenario, synthesis_decline_percent=0.0,
                                     amazon_decline_percent=0.0)).yearly
    dna = flat[flat.technology == "DNA"]
    amazon = flat[flat.technology == "Amazon Deep Archive"]
    write = dna.write_cost_usd.to_numpy()
    read = dna.read_cost_usd.to_numpy().sum()
    dna_age = dna.year.to_numpy() - scenario.dna_cost_base_year
    cloud_age = amazon.year.to_numpy() - scenario.amazon_price_base_year
    cloud = amazon.total_cost_usd.to_numpy()
    def cloud_total(rate):
        return float(np.dot(cloud, np.power(1 - rate / 100, cloud_age)))
    for synthesis_rate in s10b_x:
        dna_cost = float(np.dot(write, np.power(1 - synthesis_rate / 100, dna_age)) + read)
        def difference(rate):
            return dna_cost - cloud_total(rate)
        cloud_rate = parity_root(difference, *PARITY_AMAZON_RANGE)
        status = ("parity_in_range" if np.isfinite(cloud_rate) else
                  "DNA_more_expensive_throughout_range" if difference(PARITY_AMAZON_RANGE[0]) > 0
                  else "DNA_cheaper_throughout_range")
        rows.append({"scenario": key, "synthesis_decline_percent": float(synthesis_rate),
            "sequencing_decline_percent": scenario.sequencing_decline_percent,
            "retention_years": scenario.horizon_years, "amazon_decline_percent": cloud_rate,
            "dna_cost_usd": dna_cost,
            "comparison_cost_usd": cloud_total(cloud_rate) if np.isfinite(cloud_rate) else np.nan,
            "root_exists": bool(np.isfinite(cloud_rate)), "status": status})
s10b = pd.DataFrame(rows)
fig, ax = plt.subplots(figsize=(10.5, 6.7), layout="constrained")
for index, rate in enumerate(SEQUENCING_RATES):
    part = s10b[s10b.sequencing_decline_percent == rate]
    ax.plot(part.synthesis_decline_percent, part.amazon_decline_percent,
            color=PALETTE[index], ls=("-", "--", "-.", ":")[index], lw=2.2, label=f"{rate:g}%/year")
baseline = s10b[np.isclose(s10b.sequencing_decline_percent, BASE.sequencing_decline_percent)]
if baseline.root_exists.any():
    ax.plot(baseline.synthesis_decline_percent, baseline.amazon_decline_percent,
            color="0.25", ls=(0, (1, 2)), lw=1.8,
            label=f"Baseline: {BASE.sequencing_decline_percent:.2f}%/year")
else:
    ax.text(0.98, 0.98, f"Baseline sequencing ({BASE.sequencing_decline_percent:.2f}%/year):\\n"
            "no parity within this range", transform=ax.transAxes)
ax.scatter(BASE.synthesis_decline_percent, BASE.amazon_decline_percent,
           marker="*", s=160, facecolor="white", edgecolor="black", linewidth=1.2, zorder=5)
ax.annotate("Baseline DNA synthesis\\nand Amazon decline rates",
            (BASE.synthesis_decline_percent, BASE.amazon_decline_percent),
            xytext=(34, 0), textcoords="offset points", ha="left", va="center", zorder=6,
            arrowprops={"arrowstyle": "->", "color": "0.2", "lw": 1.1, "shrinkB": 9},
            bbox={"facecolor": "white", "alpha": 0.92, "edgecolor": "none", "pad": 2.5})
ax.set_xlim(s10b_x.min(), s10b_x.max())
ax.set_ylim(*PARITY_AMAZON_RANGE)
ax.legend()
'''
SOURCE += STYLE.replace("DNA Sequencing Cost Annual Decline (%)", "DNA Synthesis Cost Annual Decline (%)").replace(
    "DNA Synthesis Cost Annual Decline\"", "DNA Sequencing Cost Annual Decline\"").replace(
    "Baseline synthesis (", "Baseline sequencing (").replace("ncol=4", "ncol=3")
SOURCE += '''
record("S10b_synthesis_parity_by_sequencing", fig,
       f"Lifecycle cost parity between DNA and Amazon Deep Archive for deployment in {FUTURE_START}. "
       f"Workload: {BASE.archive_size_tb:g} TB; {BASE.average_asset_size_mb:g} MB/object; "
       f"{BASE.annual_retrieval_percent:g}% retrieval/year; {BASE.horizon_years} years; no discounting. "
       "The x-axis varies annual synthesis cost decline; curves hold annual sequencing cost decline "
       f"at {', '.join(f'{rate:g}%' for rate in SEQUENCING_RATES)}, plus the baseline "
       f"{BASE.sequencing_decline_percent:.2f}%. All DNA price anchors remain fixed. "
       "The y-axis solves Amazon annual decline at equal lifecycle cost. Below a curve DNA is cheaper; "
       "above it Amazon is cheaper. Missing segments have no parity in the Amazon search range "
       f"{PARITY_AMAZON_RANGE[0]:g}-{PARITY_AMAZON_RANGE[1]:g}%. "
       "Curves can nearly coincide where sequencing contributes little to total cost. "
       "The star marks baseline synthesis and Amazon inputs, irrespective of sequencing scenario. "
       "Accelerated sequencing declines are illustrative assumptions, not empirical estimates.",
       {"parity_boundary": s10b}, s10b_scenarios)
display(s10b.groupby("sequencing_decline_percent").root_exists.sum().rename("Solved parity points"))
'''


def add_s10b(cells, nbf):
    index = next(i for i, cell in enumerate(cells)
                 if cell.cell_type == "code" and 'record("S10_' in cell.source)
    cells[index + 1:index + 1] = [
        nbf.v4.new_markdown_cell("## S10b: synthesis break-even curves for several sequencing declines\n\n"
            "S10 with synthesis and sequencing exchanged: synthesis decline varies along the x-axis, "
            "while each curve fixes sequencing decline at 50%, 55%, 60%, or 65%. "
            "The baseline sequencing trajectory is also evaluated. Below each boundary DNA is cheaper; "
            "above it Amazon is cheaper. Nearly overlapping curves indicate weak sensitivity to sequencing "
            "in this range. This cell is independent of S10; run the two setup cells first."),
        nbf.v4.new_code_cell(SOURCE, metadata={"tags": ["figure"]}),
    ]
