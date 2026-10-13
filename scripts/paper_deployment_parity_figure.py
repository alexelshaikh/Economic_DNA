"""Notebook template for deployment-year parity curves on S4's decline axes."""

SOURCE = '''SAVE = True  # False previews only; True exports this figure and its data.
SYNTHESIS_RANGE = (0.0, 70.0)  # Same annual decline ranges as S4, in percent.
AMAZON_RANGE = (0.0, 20.0)
YEAR_OFFSETS = (0, 5, 10, 20)
SEARCH_END = 2500 - max(YEAR_OFFSETS)
SEARCH_START = max(BASE.start_year, BASE.dna_cost_base_year, BASE.amazon_price_base_year)

def s4c_cost_functions(year):
    scenario = replace(BASE, start_year=int(year), technologies=("DNA", "Amazon Deep Archive"))
    reference = simulate_scenario(replace(scenario, synthesis_decline_percent=0,
                                         amazon_decline_percent=0)).yearly
    dna = reference[reference.technology == "DNA"]
    amazon = reference[reference.technology == "Amazon Deep Archive"]
    write = dna.write_cost_usd.to_numpy()
    read = float(dna.read_cost_usd.sum())
    cloud = amazon.total_cost_usd.to_numpy()
    dna_age = dna.year.to_numpy() - scenario.dna_cost_base_year
    cloud_age = amazon.year.to_numpy() - scenario.amazon_price_base_year
    def dna_total(rate):
        return float(np.dot(write, np.power(1 - rate / 100, dna_age)) + read)
    def cloud_total(rate):
        return float(np.dot(cloud, np.power(1 - rate / 100, cloud_age)))
    return scenario, dna_total, cloud_total

# In this nonnegative-age domain, both costs decrease monotonically with their own decline rate.
audit_rows = []
first_viable_year = None
for year in range(SEARCH_START, SEARCH_END + 1):
    scenario, dna_total, cloud_total = s4c_cost_functions(year)
    best_gap = dna_total(SYNTHESIS_RANGE[1]) - cloud_total(AMAZON_RANGE[0])
    worst_gap = dna_total(SYNTHESIS_RANGE[0]) - cloud_total(AMAZON_RANGE[1])
    parity_possible = best_gap <= 0 <= worst_gap
    audit_rows.append({"deployment_year": year, "best_dna_minus_amazon_usd": best_gap,
                       "worst_dna_minus_amazon_usd": worst_gap, "parity_possible": parity_possible})
    if parity_possible:
        first_viable_year = year
        break
if first_viable_year is None:
    raise ValueError("No deployment year has parity in the chosen decline ranges within the search interval.")
deployment_years = [first_viable_year + offset for offset in YEAR_OFFSETS]
s4c_scenarios = {}
fig, ax = plt.subplots(figsize=(10.5, 6.7), layout="constrained")
curve_rows = []
for index, year in enumerate(deployment_years):
    scenario, dna_total, cloud_total = s4c_cost_functions(year)
    s4c_scenarios[str(year)] = scenario
    if index == 0:
        synthesis_grid = np.linspace(*SYNTHESIS_RANGE, 201)
        amazon_grid = np.linspace(*AMAZON_RANGE, 121)
        dna_values = np.array([dna_total(rate) for rate in synthesis_grid])
        amazon_values = np.array([cloud_total(rate) for rate in amazon_grid])
        z = np.log10(dna_values[None, :] / amazon_values[:, None])
        limit = max(1.0, float(np.ceil(np.abs(z).max())))
        mesh = ax.pcolormesh(synthesis_grid, amazon_grid, z, shading="auto", cmap="RdBu_r",
            norm=TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit), rasterized=True)
        grid_x, grid_y = np.meshgrid(synthesis_grid, amazon_grid)
        s4c_grid = pd.DataFrame({"deployment_year": year,
            "synthesis_decline_percent": grid_x.ravel(), "amazon_decline_percent": grid_y.ravel(),
            "log10_cost_ratio": z.ravel()})
    cloud_samples = np.linspace(*AMAZON_RANGE, 201)
    # Include exact intersections with the plot edges, even for a very short first-year curve.
    edge_rates = [parity_root(lambda rate: dna_total(synthesis_rate) - cloud_total(rate), *AMAZON_RANGE)
                  for synthesis_rate in SYNTHESIS_RANGE]
    cloud_samples = np.unique(np.r_[cloud_samples, [rate for rate in edge_rates if np.isfinite(rate)]])
    synthesis_samples = []
    for rate in cloud_samples:
        root = parity_root(lambda synthesis: dna_total(synthesis) - cloud_total(rate), *SYNTHESIS_RANGE)
        # Floating-point error at an exact edge can otherwise discard the endpoint.
        if not np.isfinite(root):
            for edge in SYNTHESIS_RANGE:
                if np.isclose(dna_total(edge), cloud_total(rate), rtol=1e-10, atol=0):
                    root = edge
        synthesis_samples.append(root)
        if np.isfinite(root):
            curve_rows.append({"deployment_year": year, "synthesis_decline_percent": root,
                "amazon_decline_percent": rate, "dna_cost_usd": dna_total(root),
                "amazon_cost_usd": cloud_total(rate)})
    label = f"{year} (first viable)" if index == 0 else str(year)
    line, = ax.plot(synthesis_samples, cloud_samples, color=PALETTE[index], lw=2.1,
        marker=("o", "s", "^", "D")[index % 4], markersize=4.5, markerfacecolor="none",
        markevery=16, label=label, zorder=4)
    import matplotlib.patheffects as path_effects
    line.set_path_effects([path_effects.Stroke(linewidth=3.7, foreground="white"), path_effects.Normal()])
ax.set_xlim(*SYNTHESIS_RANGE)
ax.set_ylim(*AMAZON_RANGE)
ax.set_xlabel("DNA Synthesis Annual Decline (%)", fontsize=12, fontweight="bold", color="0.15")
ax.set_ylabel("Amazon Annual Cost Decline (%)", fontsize=12, fontweight="bold", color="0.15")
ax.set_title("DNA / Amazon Cost Parity by Deployment Year", fontsize=14, fontweight="bold", color="0.15", pad=16)
ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}%"))
ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:g}%"))
ax.tick_params(which="major", labelsize=10.5, colors="0.15", length=0, pad=7)
ax.tick_params(which="minor", length=0)
ax.grid(False, which="both")
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_color("0.15")
    spine.set_linewidth(1.1)
colorbar = fig.colorbar(mesh, ax=ax)
colorbar.set_label(rf"$\\log_{{10}}(C_{{\\mathrm{{DNA}}}}/C_{{\\mathrm{{Amazon}}}})$ in {first_viable_year}",
                   fontsize=12, fontweight="bold", color="0.15", labelpad=10)
colorbar.ax.tick_params(labelsize=10.5, colors="0.15", length=0, pad=7)
colorbar.outline.set_edgecolor("0.15")
colorbar.outline.set_linewidth(1.1)
handles, labels = ax.get_legend_handles_labels()
shared_legend = fig.legend(handles, labels, title=r"Deployment Year $y_0$",
    loc="outside lower center", ncol=len(deployment_years), frameon=True, fancybox=True,
    facecolor="white", edgecolor="0.8", framealpha=0.9, fontsize=10.5, borderpad=0.6)
shared_legend.get_title().set_fontweight("bold")
from matplotlib.offsetbox import TextArea
legend_note = TextArea(f"Heatmap: {first_viable_year} only. Each line marks parity in its labeled year.",
                      textprops={"fontsize": 9, "fontstyle": "italic", "color": "0.3"})
legend_note.set_figure(fig)
shared_legend._legend_box.get_children().insert(0, legend_note)
s4c_parity = pd.DataFrame(curve_rows)
s4c_audit = pd.DataFrame(audit_rows)
record("S4c_deployment_year_parity", fig,
       f"DNA/Amazon lifecycle parity contours for deployments in {', '.join(map(str, deployment_years))}. "
       f"The first viable year is {first_viable_year}, searching integer deployments from {SEARCH_START} "
       f"within synthesis declines {SYNTHESIS_RANGE[0]:g}-{SYNTHESIS_RANGE[1]:g}%/year and Amazon declines "
       f"{AMAZON_RANGE[0]:g}-{AMAZON_RANGE[1]:g}%/year, as in S4. This is not the first viable year under "
       "unrestricted assumptions or at the baseline decline rates. Monotonic corner checks establish "
       "whether any cost-equality combination exists in this rectangle; the search audit is exported. "
       "Curves solve full lifecycle parity, with exact edge intersections included. At fixed Amazon decline, "
       "DNA is cheaper to the right of each year's curve. The heatmap applies ONLY to the earliest "
       "deployment year, not to every contour: negative log ratios favor DNA. All original price anchors, "
       f"sequencing decline, workload, and retention remain unchanged: {describe(BASE)}. These are conditional "
       "scenarios, not forecasts. No default parameters have been added to the figure title.",
       {"parity": s4c_parity, "earliest_year_search": s4c_audit, "reference_grid": s4c_grid}, s4c_scenarios)
print(f"First deployment with parity within S4's decline ranges: {first_viable_year}")
'''


DNA_DECLINES = SOURCE.replace('s4c', 's4d').replace('S4c_', 'S4d_')
DNA_DECLINES = DNA_DECLINES.replace('AMAZON_RANGE', 'SEQUENCING_RANGE').replace(
    'SEQUENCING_RANGE = (0.0, 20.0)', 'SEQUENCING_RANGE = (0.0, 80.0)')
start = DNA_DECLINES.index('def s4d_cost_functions(')
end = DNA_DECLINES.index('# In this nonnegative-age domain', start)
DNA_DECLINES = DNA_DECLINES[:start] + '''def s4d_cost_functions(year):
    scenario = replace(BASE, start_year=int(year), technologies=("DNA", "Amazon Deep Archive"))
    reference = simulate_scenario(replace(scenario, synthesis_decline_percent=0,
                                         sequencing_decline_percent=0)).yearly
    dna = reference[reference.technology == "DNA"]
    amazon = reference[reference.technology == "Amazon Deep Archive"]
    write = dna.write_cost_usd.to_numpy()
    read = dna.read_cost_usd.to_numpy()
    age = dna.year.to_numpy() - scenario.dna_cost_base_year
    amazon_total = float(amazon.total_cost_usd.sum())
    def dna_total(synthesis_rate, sequencing_rate):
        return float(np.dot(write, np.power(1 - synthesis_rate / 100, age))
                     + np.dot(read, np.power(1 - sequencing_rate / 100, age)))
    return scenario, dna_total, amazon_total

''' + DNA_DECLINES[end:]
DNA_DECLINES = DNA_DECLINES.replace(
    '# In this nonnegative-age domain, both costs decrease monotonically with their own decline rate.',
    '# DNA cost decreases with either decline rate; Amazon retains its baseline trajectory.')
DNA_DECLINES = DNA_DECLINES.replace('scenario, dna_total, cloud_total =', 'scenario, dna_total, amazon_total =')
DNA_DECLINES = DNA_DECLINES.replace('dna_total(SYNTHESIS_RANGE[1]) - cloud_total(SEQUENCING_RANGE[0])',
                                    'dna_total(SYNTHESIS_RANGE[1], SEQUENCING_RANGE[1]) - amazon_total')
DNA_DECLINES = DNA_DECLINES.replace('dna_total(SYNTHESIS_RANGE[0]) - cloud_total(SEQUENCING_RANGE[1])',
                                    'dna_total(SYNTHESIS_RANGE[0], SEQUENCING_RANGE[0]) - amazon_total')
DNA_DECLINES = DNA_DECLINES.replace('amazon_grid', 'sequencing_grid')
DNA_DECLINES = DNA_DECLINES.replace(
    'dna_values = np.array([dna_total(rate) for rate in synthesis_grid])\n'
    '        amazon_values = np.array([cloud_total(rate) for rate in sequencing_grid])\n'
    '        z = np.log10(dna_values[None, :] / amazon_values[:, None])',
    'dna_values = np.array([[dna_total(synthesis, sequencing) for synthesis in synthesis_grid]\n'
    '                               for sequencing in sequencing_grid])\n'
    '        z = np.log10(dna_values / amazon_total)')
DNA_DECLINES = DNA_DECLINES.replace('cloud_samples', 'sequencing_samples')
DNA_DECLINES = DNA_DECLINES.replace('dna_total(synthesis_rate) - cloud_total(rate)',
                                    'dna_total(synthesis_rate, rate) - amazon_total')
DNA_DECLINES = DNA_DECLINES.replace('dna_total(synthesis) - cloud_total(rate)',
                                    'dna_total(synthesis, rate) - amazon_total')
DNA_DECLINES = DNA_DECLINES.replace('dna_total(edge), cloud_total(rate)', 'dna_total(edge, rate), amazon_total')
DNA_DECLINES = DNA_DECLINES.replace('"dna_cost_usd": dna_total(root)', '"dna_cost_usd": dna_total(root, rate)')
DNA_DECLINES = DNA_DECLINES.replace('"amazon_cost_usd": cloud_total(rate)', '"amazon_cost_usd": amazon_total')
DNA_DECLINES = DNA_DECLINES.replace('amazon_decline_percent', 'sequencing_decline_percent')
DNA_DECLINES = DNA_DECLINES.replace('Amazon Annual Cost Decline (%)', 'DNA Sequencing Annual Decline (%)')
DNA_DECLINES = DNA_DECLINES.replace('DNA / Amazon Cost Parity by Deployment Year',
                                    'DNA Synthesis / Sequencing Declines at Amazon Cost Parity')
start = DNA_DECLINES.index('record("S4d_')
DNA_DECLINES = DNA_DECLINES[:start] + '''record("S4d_dna_declines_deployment_parity", fig,
       f"DNA synthesis/sequencing decline combinations giving Amazon lifecycle parity for deployments "
       f"in {', '.join(map(str, deployment_years))}. The first viable deployment is {first_viable_year}, "
       f"searching integer years from {SEARCH_START} within synthesis declines "
       f"{SYNTHESIS_RANGE[0]:g}-{SYNTHESIS_RANGE[1]:g}%/year and sequencing declines "
       f"{SEQUENCING_RANGE[0]:g}-{SEQUENCING_RANGE[1]:g}%/year. Amazon's annual decline remains "
       f"{BASE.amazon_decline_percent:g}% throughout. This first year is conditional on these bounds, "
       "not an unrestricted viability date or baseline-rate forecast. Monotonic corner checks establish "
       "whether parity is possible; all earlier searched years are recorded in the audit. "
       "Curves are exact lifecycle parity solutions; higher synthesis or sequencing declines favor DNA. "
       "The heatmap applies ONLY to the earliest deployment year, not all curves. Price anchors, "
       f"archive size, and retention are unchanged: {describe(BASE)}. No default parameters are added to the title.",
       {"parity": s4d_parity, "earliest_year_search": s4d_audit, "reference_grid": s4d_grid}, s4d_scenarios)
print(f"First deployment with parity within the DNA decline ranges: {first_viable_year}")
'''


try:
    from .paper_parity_line_style import line_only_deployments
except ImportError:
    from paper_parity_line_style import line_only_deployments

SOURCE = line_only_deployments(SOURCE, "s4c")
DNA_DECLINES = line_only_deployments(DNA_DECLINES, "s4d")


def add_deployment_parity_figure(cells, nbf):
    index = next(i for i, c in enumerate(cells)
                 if c.cell_type == "code" and 'record("S4b_' in c.source)
    cells[index + 1:index + 1] = [
        nbf.v4.new_markdown_cell("## S4c: cost-parity curves by deployment year\n\nSearches synthesis declines from 0% to just below 100%, keeping sequencing decline at baseline and Amazon price declines within 0-20%. Labels show the earliest deployment and +5, +10, and +20 years. Extreme declines are theoretical scenarios, not forecasts."),
        nbf.v4.new_code_cell(SOURCE, metadata={"tags": ["figure"]}),
        nbf.v4.new_markdown_cell("## S4d: synthesis and sequencing decline parity by deployment year\n\nAmazon price decline stays at baseline. Both DNA price declines range from 0% to just below 100%. Labels show the earliest deployment and +5, +10, and +20 years. Near-100% annual declines are theoretical scenarios, not forecasts."),
        nbf.v4.new_code_cell(DNA_DECLINES, metadata={"tags": ["figure"]}),
    ]
