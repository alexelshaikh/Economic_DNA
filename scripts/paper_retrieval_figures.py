"""Independent notebook cells for retrieval-dependent parity figures."""

import ast

PERCENT_AXIS = '''    ax.set_xscale("log")
    ax.set_xlim(retrieval_grid.min(), retrieval_grid.max())
    ax.set_xticks([0.00001, 0.001, 0.1, 1, 10, 100])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:.5f}".rstrip("0").rstrip(".") + "%"))
    ax.xaxis.set_minor_formatter(NullFormatter())
'''

COMMON = '''SAVE = True  # False previews only; True exports this figure and its data.
retrieval_grid = np.unique(np.r_[np.geomspace(1e-5, 100, 101), BASE.annual_retrieval_percent])
scenarios = {str(year): replace(BASE, start_year=year) for year in (BASE.start_year, 2050, 2100)}
fig, axes = plt.subplots(1, 3, figsize=(17, 5.7), layout="constrained")
rows = []
'''

STYLE = '''
for index, ((label, scenario), ax) in enumerate(zip(scenarios.items(), axes)):
''' + PERCENT_AXIS + '''    ax.axvline(BASE.annual_retrieval_percent, color="0.4", ls=":", lw=1,
               label="Baseline retrieval rate")
    ax.set_title(f"({chr(97 + index)}) Deployment Year: {scenario.start_year}",
                 fontsize=12, fontweight="bold", color="0.15", pad=14)
    ax.set_xlabel("Data Retrieval Rate (% of Archive/Year)", fontsize=12, fontweight="bold", color="0.15")
    ax.set_ylabel(YLABEL, fontsize=12, fontweight="bold", color="0.15")
    ax.grid(False, which="both")
    ax.tick_params(which="major", labelsize=10.5, colors="0.15", length=0, pad=7)
    ax.tick_params(which="minor", length=0)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("0.15")
        spine.set_linewidth(1.1)
handles, labels = axes[-1].get_legend_handles_labels()
fig.suptitle(TITLE, fontsize=14, fontweight="bold", color="0.15")
shared_legend = fig.legend(handles, labels, loc="outside lower center", ncol=len(handles),
    title=NOTE, frameon=True, fancybox=True, facecolor="white", edgecolor="0.8",
    framealpha=0.8, fontsize=10.5, handlelength=2.8, borderpad=0.6)
shared_legend.get_title().set_fontweight("bold")
'''

DECLINE = COMMON + '''MAX_DECLINE = 99.999999  # Annual declines must stay below 100%.
for (label, scenario), ax in zip(scenarios.items(), axes):
    # Public-model streams at 0% synthesis decline and 1% retrieval permit exact rescaling.
    reference = simulate_scenario(replace(scenario, synthesis_decline_percent=0,
                                         annual_retrieval_percent=1)).yearly
    dna = reference[reference.technology == "DNA"]
    age = dna.year.to_numpy() - scenario.dna_cost_base_year
    writes = dna.write_cost_usd.to_numpy()
    reads = dna.read_cost_usd.sum()
    for tech, marker in zip(PEERS, ("o", "s", "^")):
        peer = reference[reference.technology == tech]
        peer_fixed = peer.write_cost_usd.sum() + peer.maintenance_cost_usd.sum()
        peer_reads = peer.read_cost_usd.sum()
        rates = []
        for retrieval in retrieval_grid:
            peer_cost = peer_fixed + peer_reads * retrieval
            def difference(rate):
                return np.dot(writes, np.power(1 - rate / 100, age)) + reads * retrieval - peer_cost
            if difference(0) <= 0:
                rate, status = 0.0, "DNA_already_at_or_below_comparator"
            elif difference(MAX_DECLINE) > 0:
                rate, status = np.nan, "no_parity_below_100_percent_decline"
            else:
                rate = parity_root(difference, 0, MAX_DECLINE)
                status = "parity"
            rates.append(rate)
            rows.append({"start_year": scenario.start_year, "technology": tech,
                         "annual_retrieval_percent": retrieval,
                         "synthesis_decline_percent": rate, "status": status})
        ax.plot(retrieval_grid, rates, color=COLORS[tech], ls=STYLES[tech],
                marker=marker, markersize=4.5, markerfacecolor="none", markevery=8, label=tech)
    ax.set_ylim(0, 100)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}%"))
    if not any(np.isfinite(line.get_ydata()).any() for line in ax.lines):
        ax.text(0.36, 0.5, "No feasible synthesis decline rate\\nInitial synthesis alone exceeds\\nthe alternatives' lifecycle costs",
                transform=ax.transAxes, ha="center", va="center", fontsize=10.5, color="0.25")
''' + STYLE.replace('YLABEL', repr('Minimum Annual Synthesis Cost Decline')) \
    .replace('TITLE', repr('Break-Even Synthesis Decline Rate by Data Retrieval Rate')) \
    .replace('NOTE', repr('Missing segments: no synthesis decline below 100% per year achieves parity.')) + '''
s2b = pd.DataFrame(rows)
record("S2b_retrieval_vs_synthesis_decline", fig,
       f"Minimum annual synthesis decline required for DNA lifecycle cost to be no greater than each alternative. "
       f"Deployments: {BASE.start_year}, 2050, 2100. Baseline: {describe(BASE)}. "
       "Only retrieval and synthesis decline vary; sequencing decline and all price anchors remain unchanged. "
       "Positive thresholds solve full lifecycle parity; zero means DNA is already competitive without synthesis decline. "
       "Missing values mean no solution up to 99.999999% annual decline. In the anchor year, future declines "
       "cannot reduce the initial synthesis bill. The retrieval axis is logarithmic with plain percentage labels.",
       {"thresholds": s2b}, scenarios)
'''

OBJECT_SIZE = COMMON + '''
object_comparators = [tech for tech in PEERS if tech != "Tape On-premise"]
for (label, scenario), ax in zip(scenarios.items(), axes):
    reference_scenario = replace(scenario, annual_retrieval_percent=1)
    reference = simulate_scenario(reference_scenario).yearly
    # Removing request fees isolates the volume-dependent cost A in C = A + B/object_size.
    volume_only = simulate_scenario(replace(reference_scenario,
        amazon_put_usd_per_request=0, amazon_bulk_restore_usd_per_request=0,
        azure_write_usd_per_request=0, azure_read_usd_per_request=0)).yearly
    dna = reference[reference.technology == "DNA"]
    dna_cost = dna.write_cost_usd.sum() + dna.read_cost_usd.sum() * retrieval_grid
    for tech, marker in zip(object_comparators, ("o", "s")):
        peer = reference[reference.technology == tech]
        volume = volume_only[volume_only.technology == tech]
        fixed = volume.write_cost_usd.sum() + volume.maintenance_cost_usd.sum()
        volume_cost = fixed + volume.read_cost_usd.sum() * retrieval_grid
        request_write = peer.write_cost_usd.sum() - volume.write_cost_usd.sum()
        request_read = peer.read_cost_usd.sum() - volume.read_cost_usd.sum()
        coefficient = (request_write + request_read * retrieval_grid) * scenario.average_asset_size_mb
        thresholds = []
        for retrieval, dna_total, volume_total, b in zip(retrieval_grid, dna_cost, volume_cost, coefficient):
            gap = dna_total - volume_total
            threshold = np.nan
            if b <= 0:
                status = "no_object_size_dependence"
            elif gap <= 0:
                status = "DNA_competitive_at_every_object_size"
            else:
                candidate = b / gap
                if candidate <= scenario.archive_size_mb:
                    threshold, status = candidate, "parity"
                else:
                    status = "DNA_competitive_at_every_valid_object_size"
            thresholds.append(threshold)
            rows.append({"start_year": scenario.start_year, "technology": tech,
                         "annual_retrieval_percent": retrieval, "object_size_mb": threshold,
                         "object_size_bytes": threshold * 1e6, "status": status})
        ax.plot(retrieval_grid, thresholds, color=COLORS[tech], ls=STYLES[tech],
                marker=marker, markersize=4.5, markerfacecolor="none", markevery=8, label=tech)
    ax.set_yscale("log")
    ax.axhline(1e-6, color="0.4", ls="--", lw=1, label="One byte")
''' + STYLE.replace('YLABEL', repr('Break-Even Average Object Size (MB)')) \
    .replace('TITLE', repr('Average Object Size at Cost Parity')) \
    .replace('NOTE', repr('Smaller objects favor DNA. Thresholds below one byte are theoretical.')) + '''
s2c = pd.DataFrame(rows)
record("S2c_retrieval_vs_object_size", fig,
       f"Average object size giving DNA/comparator lifecycle parity for deployments in {BASE.start_year}, "
       f"2050, and 2100. Baseline: {describe(BASE)}. Archive volume is fixed; the number of objects varies "
       "inversely with average object size. Cloud request charges rise for smaller objects, while DNA cost "
       "is unchanged. Sizes below a plotted threshold favor DNA; larger sizes favor the comparator. "
       "Thresholds below one byte are mathematical extrapolations, not physically attainable objects. "
       "Tape has no per-object charges and therefore no unique object-size threshold. Other absent values "
       "are classified in the CSV. Both axes are logarithmic; decimal MB are used. Storage-service minimum "
       "billable object sizes and extra DNA addressing overhead are not added to the existing model.",
       {"thresholds": s2c}, scenarios)
'''

CALIBRATED_OBJECT_SIZE = OBJECT_SIZE[:OBJECT_SIZE.index('s2c = pd.DataFrame(rows)')]
CALIBRATED_OBJECT_SIZE = CALIBRATED_OBJECT_SIZE.replace(
    'scenarios = {str(year): replace(BASE, start_year=year) for year in (BASE.start_year, 2050, 2100)}',
    '''calibration_scenario = replace(BASE, start_year=2050)
calibration_thresholds = threshold_at_start(calibration_scenario).set_index("technology")
parity_price_2050 = float(calibration_thresholds.loc["Amazon Deep Archive", "threshold_start_usd_per_mb"])
if not np.isfinite(parity_price_2050) or parity_price_2050 < 0:
    raise ValueError("No nonnegative synthesis price gives Amazon parity in 2050 with this baseline.")
calibrated_base = future_at_price(calibration_scenario, parity_price_2050)
scenarios = {str(year): replace(calibrated_base, start_year=year) for year in (2050, 2075, 2100)}
calibration = pd.DataFrame([{
    "deployment_year": 2050, "annual_retrieval_percent": BASE.annual_retrieval_percent,
    "average_object_size_mb": BASE.average_asset_size_mb,
    "synthesis_price_2050_usd_per_mb": parity_price_2050,
    "synthesis_anchor_year": BASE.dna_cost_base_year,
    "synthesis_anchor_usd_per_mb": calibrated_base.dna_synthesis_cost_per_mb,
}])''')
CALIBRATED_OBJECT_SIZE = CALIBRATED_OBJECT_SIZE.replace(
    'ax.axhline(1e-6, color="0.4", ls="--", lw=1, label="One byte")',
    '''ax.axhline(BASE.average_asset_size_mb, color="0.4", ls="--", lw=1, label="Baseline object size")
    if not any(np.isfinite(line.get_ydata()).any() for line in ax.lines[:len(object_comparators)]):
        ax.text(0.35, 0.5, "DNA is competitive with both clouds\\nfor every valid object size\\nthroughout this retrieval range",
                transform=ax.transAxes, ha="center", va="center", fontsize=10, color="0.25")
    if scenario.start_year == 2050:
        ax.scatter(BASE.annual_retrieval_percent, BASE.average_asset_size_mb,
                   color=COLORS["Amazon Deep Archive"], s=65, zorder=5)''')
CALIBRATED_OBJECT_SIZE = CALIBRATED_OBJECT_SIZE.replace(
    "'Average Object Size at Cost Parity'",
    "'Break-Even Object Size with Synthesis Priced for Amazon Parity in 2050'")
CALIBRATED_OBJECT_SIZE = CALIBRATED_OBJECT_SIZE.replace(
    "'Smaller objects favor DNA. Thresholds below one byte are theoretical.'",
    "'Smaller objects favor DNA. Missing curve segments: DNA is competitive at every valid object size.'")
CALIBRATED_OBJECT_SIZE += '''
s2d = pd.DataFrame(rows)
valid_sizes = s2d.object_size_mb.dropna()
for ax in axes:
    ax.set_ylim(min(BASE.average_asset_size_mb, valid_sizes.min()) * 0.5
                if len(valid_sizes) else BASE.average_asset_size_mb / 10,
                BASE.archive_size_mb * 1.1)
record("S2d_retrieval_vs_object_size_calibrated", fig,
       f"Break-even average object size at fixed archive volume for deployments in 2050, 2075, and 2100. "
       f"Synthesis price is calibrated once to full DNA/Amazon lifecycle parity for a 2050 deployment "
       f"at {BASE.annual_retrieval_percent:g}% annual retrieval and {BASE.average_asset_size_mb:g} MB/object: "
       f"USD {parity_price_2050:.10g}/MB in 2050, equivalent to USD "
       f"{calibrated_base.dna_synthesis_cost_per_mb:.10g}/MB at the {BASE.dna_cost_base_year} anchor. "
       f"The original synthesis decline ({BASE.synthesis_decline_percent:g}%/year) and sequencing "
       "trajectory are retained; later deployments are not recalibrated. The orange point marks the "
       "2050 calibration. Below each threshold DNA is cheaper; above it the cloud comparator is cheaper. "
       "Absent thresholds indicate DNA is competitive for all valid object sizes (up to the archive size); "
       "statuses are exported in the CSV. Tape is excluded because it has no per-object costs. "
       "Both axes are logarithmic; no additional object-overhead or minimum-billing assumptions are added.",
       {"thresholds": s2d, "calibration": calibration}, scenarios)
'''


def refine_title_and_legend(source, legend_heading="Comparison Storage Model", legend_note=None):
    for field in ("archive_size_tb", "horizon_years", "dna_durability_years", "average_asset_size_mb"):
        source = source.replace(f'{{BASE.{field}:g}}', f'{{BASE.{field}:,g}}')
    tree = ast.parse(source)
    call = next(node for node in ast.walk(tree) if isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute) and node.func.attr == "legend"
                and any(keyword.arg == "title" for keyword in node.keywords))
    title = next(keyword.value for keyword in call.keywords if keyword.arg == "title")
    note = (ast.literal_eval(title).removeprefix("Comparison Storage Model\n")
            if legend_note is None else legend_note)
    source = source.replace('title=' + ast.get_source_segment(source, title),
                            'title=' + repr(legend_heading), 1)
    anchor = 'shared_legend.get_title().set_fontweight("bold")'
    source = source.replace(anchor, anchor + '''
from matplotlib.offsetbox import TextArea
legend_note = TextArea(NOTE, textprops={"fontsize": 9, "fontstyle": "italic",
                                      "fontweight": "normal", "color": "0.3"})
legend_note.set_figure(fig)
# Pack both text styles inside the legend so constrained layout reserves their space.
shared_legend._legend_box.get_children().insert(0, legend_note)
'''.replace('NOTE', repr(note)), 1)
    return source


def refine_s3_titles(source):
    title = '''fig.suptitle(
    "DNA Cost Ratios and Break-Even by Deployment Year\\n"
    f"Archive: {BASE.archive_size_tb:,g} TB | Retention: {BASE.horizon_years:,g} years | "
    f"DNA durability: {BASE.dna_durability_years:,g} years | Object size: {BASE.average_asset_size_mb:,g} MB\\n"
    f"Synthesis decline: varied (see legend) | Sequencing decline: {BASE.sequencing_decline_percent:.2f}%/year | "
    f"Retrieval: {BASE.annual_retrieval_percent:g}%/year",
    fontsize=14, fontweight="bold", color="0.15")'''
    start = source.index('fig.suptitle(')
    end = source.index('\n', start)
    source = source[:start] + title + source[end:]
    source = source.replace('figsize=(14, 5.7)', 'figsize=(17, 6.7)')
    return refine_title_and_legend(source, "Annual Synthesis Cost Decline",
        "Ratios below 1 favor DNA. Larger circles mark the first deployment year at or below cost parity.")


def add_parameter_titles(source, kind):
    titles = {
        "S2": "Break-Even Synthesis Price by Data Retrieval Rate",
        "S2b": "Break-Even Synthesis Decline Rate by Data Retrieval Rate",
        "S2c": "Break-Even Average Object Size by Data Retrieval Rate",
        "S2d": "Break-Even Average Object Size with Synthesis Priced for Amazon Parity in 2050",
    }
    object_label = ('Object size: {BASE.average_asset_size_mb:g} MB'
                    if kind in ("S2", "S2b") else 'Object size: solved at parity')
    synthesis_label = ('Synthesis decline: solved at parity' if kind == "S2b"
                       else 'Synthesis decline: {BASE.synthesis_decline_percent:.2f}%/year')
    title_code = ('fig.suptitle(\n'
        f'    "{titles[kind]}\\n"\n'
        '    f"Archive: {BASE.archive_size_tb:g} TB | Retention: {BASE.horizon_years:g} years | "\n'
        f'    f"DNA durability: {{BASE.dna_durability_years:g}} years | {object_label}\\n"\n'
        f'    f"{synthesis_label} | Sequencing decline: {{BASE.sequencing_decline_percent:.2f}}%/year | "\n'
        '    f"Retrieval: varied (reference {BASE.annual_retrieval_percent:g}%/year)"\n')
    if kind == "S2d":
        title_code += ('    f"\\n2050 calibration: synthesis ${parity_price_2050:.4g}/MB | "\n'
                       '    f"Reference object size: {BASE.average_asset_size_mb:g} MB"\n')
    title_code += '    , fontsize=14, fontweight="bold", color="0.15")'
    start = source.index('fig.suptitle(')
    end = source.index('\n', start)
    source = source[:start] + title_code + source[end:]
    source = source.replace('Deployment Year:', 'Deployment Year $y_0$:')
    return refine_title_and_legend(source.replace('figsize=(17, 5.7)', 'figsize=(17, 6.7)'))


DECLINE = add_parameter_titles(DECLINE, "S2b")
OBJECT_SIZE = add_parameter_titles(OBJECT_SIZE, "S2c")
CALIBRATED_OBJECT_SIZE = add_parameter_titles(CALIBRATED_OBJECT_SIZE, "S2d")


def add_retrieval_figures(cells, nbf):
    index = next(i for i, c in enumerate(cells)
                 if c.cell_type == "code" and 'record("S2_retrieval_' in c.source)
    cells[index + 1:index + 1] = [
        nbf.v4.new_markdown_cell("## S2b: retrieval vs. synthesis decline needed for parity"),
        nbf.v4.new_code_cell(DECLINE, metadata={"tags": ["figure"]}),
        nbf.v4.new_markdown_cell("## S2c: retrieval vs. break-even average object size\n\nArchive volume stays fixed. Smaller objects increase cloud request charges; tape has no per-object threshold. Sub-byte thresholds are theoretical."),
        nbf.v4.new_code_cell(OBJECT_SIZE, metadata={"tags": ["figure"]}),
        nbf.v4.new_markdown_cell("## S2d: break-even object size with synthesis calibrated to Amazon parity in 2050\n\nOne synthesis price calibration at baseline retrieval and object size is carried forward to 2050, 2075, and 2100 with the original decline rate."),
        nbf.v4.new_code_cell(CALIBRATED_OBJECT_SIZE, metadata={"tags": ["figure"]}),
    ]
