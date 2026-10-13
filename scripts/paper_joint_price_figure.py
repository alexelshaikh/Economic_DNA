"""Notebook template for the joint synthesis/sequencing price heatmap."""

SOURCE = '''SAVE = True  # False previews only; True exports this figure and its data.
DEPLOYMENT_YEAR = FUTURE_START
SYNTHESIS_MAX_USD_PER_MB = None  # None focuses the plot on the parity boundary.
SEQUENCING_MAX_USD_PER_MB = None
GRID_POINTS = 201
s4b_scenario = replace(BASE, start_year=DEPLOYMENT_YEAR, technologies=("DNA", "Amazon Deep Archive"))
age = s4b_scenario.start_year - s4b_scenario.dna_cost_base_year
synthesis_factor = (1 - s4b_scenario.synthesis_decline_percent / 100) ** age
sequencing_factor = (1 - s4b_scenario.sequencing_decline_percent / 100) ** age
# Unit deployment-year prices give exact lifecycle coefficients from the public model.
unit_scenario = replace(s4b_scenario, dna_synthesis_cost_per_mb=1 / synthesis_factor,
                        dna_sequencing_cost_per_mb=1 / sequencing_factor)
unit = simulate_scenario(unit_scenario)
dna = unit.yearly[unit.yearly.technology == "DNA"]
write_coefficient = float(dna.write_cost_usd.sum())
read_coefficient = float(dna.read_cost_usd.sum())
amazon_total = float(unit.totals.set_index("technology").loc["Amazon Deep Archive", "total_cost_usd"])
if write_coefficient <= 0 or amazon_total <= 0:
    raise ValueError("This heatmap requires positive synthesis exposure and Amazon lifecycle cost.")
synthesis_intercept = amazon_total / write_coefficient
sequencing_intercept = amazon_total / read_coefficient if read_coefficient > 0 else np.nan
baseline_synthesis = s4b_scenario.dna_synthesis_cost_per_mb * synthesis_factor
baseline_sequencing = s4b_scenario.dna_sequencing_cost_per_mb * sequencing_factor
synthesis_max = (1.25 * synthesis_intercept if SYNTHESIS_MAX_USD_PER_MB is None
                 else float(SYNTHESIS_MAX_USD_PER_MB))
sequencing_max = (1.25 * sequencing_intercept if np.isfinite(sequencing_intercept)
                  else max(baseline_sequencing * 1.25, 1e-6)) if SEQUENCING_MAX_USD_PER_MB is None else float(SEQUENCING_MAX_USD_PER_MB)
if synthesis_max <= 0 or sequencing_max <= 0:
    raise ValueError("Both price limits must be positive.")
synthesis_prices = np.linspace(0, synthesis_max, GRID_POINTS)
sequencing_prices = np.linspace(0, sequencing_max, GRID_POINTS)
x, y = np.meshgrid(synthesis_prices, sequencing_prices)
dna_cost = write_coefficient * x + read_coefficient * y
ratio = dna_cost / amazon_total
s4b = pd.DataFrame({"synthesis_price_deployment_usd_per_mb": x.ravel(),
                    "sequencing_price_deployment_usd_per_mb": y.ravel(),
                    "synthesis_price_anchor_usd_per_mb": x.ravel() / synthesis_factor,
                    "sequencing_price_anchor_usd_per_mb": y.ravel() / sequencing_factor,
                    "dna_cost_usd": dna_cost.ravel(), "amazon_cost_usd": amazon_total,
                    "cost_ratio": ratio.ravel()})
fig, ax = plt.subplots(figsize=(9.5, 6.2), layout="constrained")
mesh = ax.pcolormesh(x, y, ratio, shading="auto", cmap="RdBu_r",
                    norm=TwoSlopeNorm(vmin=0, vcenter=1, vmax=max(2, float(ratio.max()))), rasterized=True)
if read_coefficient > 0:
    lower = max(0, (amazon_total - read_coefficient * sequencing_max) / write_coefficient)
    upper = min(synthesis_intercept, synthesis_max)
    parity_x = np.linspace(lower, upper, GRID_POINTS) if lower <= upper else np.array([])
    parity_y = np.maximum(0, (amazon_total - write_coefficient * parity_x) / read_coefficient)
else:
    parity_x = np.full(GRID_POINTS, synthesis_intercept) if synthesis_intercept <= synthesis_max else np.array([])
    parity_y = np.linspace(0, sequencing_max, len(parity_x))
s4b_parity = pd.DataFrame({"synthesis_price_deployment_usd_per_mb": parity_x,
                          "sequencing_price_deployment_usd_per_mb": parity_y})
if len(parity_x):
    ax.plot(parity_x, parity_y, color="0.1", lw=1.8, label="Cost parity")
else:
    ax.text(0.5, 0.5, "No cost parity within the displayed price range",
            transform=ax.transAxes, ha="center", color="0.15")
parity_at_baseline_read = (amazon_total - read_coefficient * baseline_sequencing) / write_coefficient
if 0 <= parity_at_baseline_read <= synthesis_max and 0 <= baseline_sequencing <= sequencing_max:
    ax.scatter(parity_at_baseline_read, baseline_sequencing, marker="*", s=210,
               color="white", edgecolor="black", zorder=5, clip_on=False)
    ax.annotate("Parity at baseline\\nsequencing price", (parity_at_baseline_read, baseline_sequencing),
                xytext=(0.46, 0.18), textcoords="axes fraction", ha="center", fontsize=10,
                arrowprops={"arrowstyle": "->", "color": "0.2", "lw": 1.1, "shrinkB": 10},
                bbox={"facecolor": "white", "alpha": 0.9, "edgecolor": "0.8", "pad": 2.5})
ax.set_xlim(0, synthesis_max)
ax.set_ylim(0, sequencing_max)
ax.set_xlabel(f"DNA Synthesis Price in {DEPLOYMENT_YEAR} (USD/MB)", fontsize=12, fontweight="bold", color="0.15")
ax.set_ylabel(f"DNA Sequencing Price in {DEPLOYMENT_YEAR} (USD/MB)", fontsize=12, fontweight="bold", color="0.15")
ax.set_title(f"DNA / Amazon Cost Parity: Deployment Year $y_0$ = {DEPLOYMENT_YEAR}",
             fontsize=14, fontweight="bold", color="0.15", pad=16)
ax.ticklabel_format(axis="both", style="sci", scilimits=(0, 0), useMathText=True)
ax.tick_params(which="major", labelsize=10.5, colors="0.15", length=0, pad=7)
ax.tick_params(which="minor", length=0)
ax.grid(False, which="both")
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_color("0.15")
    spine.set_linewidth(1.1)
colorbar = fig.colorbar(mesh, ax=ax)
colorbar.set_label("DNA / Amazon Lifecycle Cost", fontsize=12, fontweight="bold", color="0.15", labelpad=10)
colorbar.ax.tick_params(labelsize=10.5, colors="0.15", length=0, pad=7)
colorbar.outline.set_edgecolor("0.15")
colorbar.outline.set_linewidth(1.1)
if len(parity_x):
    shared_legend = ax.legend(loc="upper right", title="Below 1: DNA cheaper; above 1: Amazon cheaper",
        frameon=True, fancybox=True, facecolor="white", edgecolor="0.8", framealpha=0.9,
        fontsize=10.5, borderpad=0.6)
    shared_legend.get_title().set_fontsize(9)
    shared_legend.get_title().set_fontstyle("italic")
s4b_reference = pd.DataFrame([{"deployment_year": DEPLOYMENT_YEAR,
    "baseline_synthesis_deployment_usd_per_mb": baseline_synthesis,
    "baseline_sequencing_deployment_usd_per_mb": baseline_sequencing,
    "synthesis_parity_at_baseline_sequencing_usd_per_mb": parity_at_baseline_read,
    "amazon_lifecycle_usd": amazon_total,
    "write_coefficient": write_coefficient, "read_coefficient": read_coefficient}])
record("S4b_joint_synthesis_sequencing_prices", fig,
       f"Joint synthesis and sequencing price sensitivity for DNA/Amazon lifecycle parity. {describe(s4b_scenario)}. "
       f"Both axes are USD per logical MB at deployment in {DEPLOYMENT_YEAR}, not anchor-year prices. "
       f"Synthesis and sequencing decline remain {s4b_scenario.synthesis_decline_percent:.4g}% and "
       f"{s4b_scenario.sequencing_decline_percent:.4g}% per year. Amazon tariffs and all other inputs remain unchanged. "
       "The heatmap shows the linear lifecycle cost ratio, centered on one (parity); unlike S4 it is not a log ratio. "
       "The exact black boundary includes nonnegative price combinations giving equal full lifecycle costs, "
       "including DNA replacements and retrieval over the retention period. The white star, when in range, "
       "marks parity holding sequencing at its baseline trajectory; it is not the original synthesis baseline. "
       "The default price window focuses on affordable combinations; the fitted baseline synthesis price may "
       "lie outside it and is recorded in the reference CSV. No negative prices are used.",
       {"grid": s4b, "parity": s4b_parity, "reference": s4b_reference}, {"deployment": s4b_scenario})
'''


try:
    from .paper_parity_line_style import line_only_prices
except ImportError:
    from paper_parity_line_style import line_only_prices

SOURCE = line_only_prices(SOURCE)


def add_joint_price_figure(cells, nbf):
    index = max(i for i, c in enumerate(cells)
                if c.cell_type == "code" and 'record("S4_' in c.source)
    cells[index + 1:index + 1] = [
        nbf.v4.new_markdown_cell("## S4b: joint synthesis and sequencing prices at Amazon parity\n\nPrices on both axes are at deployment. Annual decline rates stay fixed. The labeled line marks lifecycle parity."),
        nbf.v4.new_code_cell(SOURCE, metadata={"tags": ["figure"]}),
    ]
