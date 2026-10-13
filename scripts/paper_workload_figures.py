"""Standalone workload-parity notebook cells using the public lifecycle model."""

import json

COMMON = '''SAVE = True  # False previews only; True exports this figure and its data.
USE_PARITY_CALIBRATION = True  # False uses the original fitted synthesis price.
DEPLOYMENT_YEARS = (2050, 2075, 2100)
MAX_RETENTION = 10000
ARCHIVE_RANGE_TB = (0.001, 1000000.0)
workload_base = BASE
if USE_PARITY_CALIBRATION:
    calibration_scenario = replace(BASE, start_year=2050)
    price = threshold_at_start(calibration_scenario).set_index("technology").loc[
        "Amazon Deep Archive", "threshold_start_usd_per_mb"]
    if not np.isfinite(price) or price < 0:
        raise ValueError("No nonnegative synthesis price gives Amazon parity in 2050 with this baseline.")
    workload_base = future_at_price(calibration_scenario, float(price))
scenarios = {str(year): replace(workload_base, start_year=year,
                               technologies=("DNA", "Amazon Deep Archive")) for year in DEPLOYMENT_YEARS}

def workload_costs(scenario, horizon):
    probe = replace(scenario, horizon_years=int(horizon))
    result = simulate_scenario(probe).yearly
    volume = simulate_scenario(replace(probe, amazon_put_usd_per_request=0,
                                      amazon_bulk_restore_usd_per_request=0)).yearly
    dna = result[result.technology == "DNA"].total_cost_usd.to_numpy().cumsum()
    cloud = result[result.technology == "Amazon Deep Archive"].total_cost_usd.to_numpy().cumsum()
    fixed = volume[volume.technology == "Amazon Deep Archive"].total_cost_usd.to_numpy().cumsum()
    # At fixed archive volume, cloud lifecycle cost is A + B / average object size.
    requests = (cloud - fixed) * scenario.average_asset_size_mb
    return dna, fixed, requests

fig, ax = plt.subplots(figsize=(10.5, 6.7), layout="constrained")
rows = []
status_notes = []
from matplotlib.contour import ContourSet

def workload_curve(x, y, year, color):
    points = np.column_stack([x, y])
    valid = np.isfinite(points).all(axis=1) & (points > 0).all(axis=1)
    runs = np.split(np.flatnonzero(valid), np.flatnonzero(np.diff(np.flatnonzero(valid)) > 1) + 1)
    segments = [points[run] for run in runs if len(run) > 1]
    if not segments:
        return
    contour = ContourSet(ax, [0], [segments], colors=[color], linewidths=1.8)
    segment = max(segments, key=len)
    target = segment[len(segment) // 2]
    labels = ax.clabel(contour, fmt={0: f"$y_0 = {year}$"}, inline=True, inline_spacing=5,
                       fontsize=10.5, manual=[tuple(target)])
    for label in labels:
        if label.get_rotation() > 180:
            label.set_rotation((label.get_rotation() + 180) % 360)
'''

STYLE = '''
ax.set_xlabel(XLABEL, fontsize=12, fontweight="bold", color="0.15", labelpad=8)
ax.set_ylabel(YLABEL, fontsize=12, fontweight="bold", color="0.15", labelpad=8)
price_note = "Synthesis calibrated to Amazon parity in 2050" if USE_PARITY_CALIBRATION else "Original fitted synthesis price"
ax.set_title(TITLE + "\\n" + price_note, fontsize=14, fontweight="bold", color="0.15", pad=16)
for axis in (ax.xaxis, ax.yaxis):
    axis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:,.6f}".rstrip("0").rstrip(".")))
    axis.set_minor_formatter(NullFormatter())
ax.grid(False, which="both")
ax.tick_params(which="major", labelsize=10.5, colors="0.15", length=0, pad=7)
ax.tick_params(which="minor", length=0)
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_color("0.15")
    spine.set_linewidth(1.1)
note = NOTE
if status_notes:
    note += "\\n" + "; ".join(status_notes)
fig.supxlabel(note, fontsize=9, fontstyle="italic", color="0.3")
table = pd.DataFrame(rows)
configuration = pd.DataFrame([{
    "uses_2050_price_calibration": USE_PARITY_CALIBRATION,
    "synthesis_anchor_year": workload_base.dna_cost_base_year,
    "synthesis_anchor_usd_per_mb": workload_base.dna_synthesis_cost_per_mb,
    "annual_retrieval_percent": workload_base.annual_retrieval_percent,
    "baseline_archive_tb": workload_base.archive_size_tb,
    "baseline_object_mb": workload_base.average_asset_size_mb,
    "baseline_retention_years": workload_base.horizon_years,
}])
record(STEM, fig,
       CAPTION + f" Price basis: {price_note}. Calibration, when enabled, is performed once at "
       f"{BASE.annual_retrieval_percent:g}% retrieval/year, {BASE.average_asset_size_mb:g} MB/object, "
       f"{BASE.archive_size_tb:g} TB and {BASE.horizon_years:g} years retention. The calibrated anchor price "
       "then follows the original synthesis decline in all deployments. Sequencing, cloud prices, "
       "and their decline rates retain the baseline assumptions. There are no fixed infrastructure costs, "
       "volume discounts, or minimum billable object sizes added to this proportional-cost model. "
       "Sub-byte object thresholds are mathematical extrapolations. Invalid objects larger than their "
       "archive are not plotted. Boundaries and statuses are exported with the configuration.",
       {"boundaries": table, "configuration": configuration}, scenarios)
'''


def finish(body, xlabel, ylabel, title, note, stem, caption):
    style = STYLE
    for key, value in (("XLABEL", xlabel), ("YLABEL", ylabel), ("TITLE", title),
                       ("NOTE", note), ("STEM", stem), ("CAPTION", caption)):
        style = style.replace(key, json.dumps(value))
    return COMMON + body + style


RETENTION_ARCHIVE = finish('''
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(1, MAX_RETENTION)
ax.set_ylim(*ARCHIVE_RANGE_TB)
archive_sizes = np.geomspace(*ARCHIVE_RANGE_TB, 101)
for index, (year, scenario) in enumerate(zip(DEPLOYMENT_YEARS, scenarios.values())):
    dna, fixed, requests = workload_costs(scenario, MAX_RETENTION)
    cloud = fixed + requests / scenario.average_asset_size_mb
    # Whole-year retention is the model's domain. A boundary need not be exact equality.
    competitive = dna <= cloud + 1e-12 * np.maximum(dna, cloud)
    first = int(np.flatnonzero(competitive)[0] + 1) if competitive.any() else None
    if first is None:
        status_notes.append(f"{year}: no competitive retention through {MAX_RETENTION:,} years")
    elif first == 1:
        status_notes.append(f"{year}: competitive from 1 year")
    for archive in archive_sizes:
        valid = archive * 1e6 >= scenario.average_asset_size_mb
        rows.append({"deployment_year": year, "archive_size_tb": archive,
            "minimum_retention_years": first if valid and first is not None else np.nan,
            "status": "object_exceeds_archive" if not valid else
                ("no_competitive_retention" if first is None else "first_competitive_retention")})
    valid_sizes = archive_sizes[archive_sizes * 1e6 >= scenario.average_asset_size_mb]
    if first is not None and first > 1 and len(valid_sizes):
        workload_curve(np.full(len(valid_sizes), first), valid_sizes, year, PALETTE[index])
''', "Retention (Years)", "Archive Size (TB)",
    "Minimum Competitive Retention vs. Archive Size",
    "Archive size does not change the threshold at fixed object size; retention uses whole years.",
    "S4e_retention_vs_archive_size",
    "Minimum whole-year retention with DNA lifecycle cost no greater than Amazon's, varying archive size "
    "while holding average object size fixed. Costs scale proportionally with archive volume, so a vertical "
    "boundary is expected. Each boundary is the first qualifying integer horizon, not an interpolated "
    "fractional-year equality. A relative 1e-12 tolerance handles roundoff at exact calibration parity. ")

RETENTION_OBJECT = finish('''
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(1, MAX_RETENTION)
ax.set_ylim(1e-6, workload_base.archive_size_mb)
retentions = np.arange(1, MAX_RETENTION + 1)
for index, (year, scenario) in enumerate(zip(DEPLOYMENT_YEARS, scenarios.values())):
    dna, fixed, requests = workload_costs(scenario, MAX_RETENTION)
    gap = dna - fixed
    thresholds = np.divide(requests, gap, out=np.full_like(gap, np.nan), where=(gap > 0) & (requests > 0))
    statuses = np.where(gap <= 0, "DNA_competitive_at_every_object_size",
                       np.where(thresholds > scenario.archive_size_mb, "DNA_competitive_at_every_valid_object_size", "parity"))
    invalid = thresholds > scenario.archive_size_mb
    thresholds[invalid] = np.nan
    for retention, threshold, status in zip(retentions, thresholds, statuses):
        rows.append({"deployment_year": year, "retention_years": int(retention),
                     "object_size_mb": threshold, "status": status})
    if not np.isfinite(thresholds).any():
        status_notes.append(f"{year}: DNA competitive for all valid object sizes")
    color = PALETTE[index % len(PALETTE)]
    valid_indices = np.flatnonzero(np.isfinite(thresholds))
    if len(valid_indices):
        last = valid_indices[-1]
        if last < len(retentions) - 1:
            ax.plot(retentions[last], thresholds[last], "o", color=color,
                    markersize=6, markerfacecolor="none", markeredgewidth=1.4)
        if len(valid_indices) <= 3:
            ax.plot(retentions, thresholds, color=color, linewidth=1.8)
            ax.annotate(f"$y_0 = {year}$", (retentions[last], thresholds[last]),
                        xytext=(0, 35) if year == 2075 else (12, 8),
                        textcoords="offset points", color=color,
                        ha="center" if year == 2075 else "left",
                        arrowprops=dict(arrowstyle="-", color=color, lw=0.8) if year == 2075 else None,
                        fontsize=10.5, va="bottom")
        else:
            workload_curve(retentions, thresholds, year, color)
        if year == 2075 and last < len(retentions) - 1 and np.all(gap[last + 1:] <= 0):
            status_notes.append(f"2075: DNA cheaper at every object size from {retentions[last + 1]} years retention")
''', "Retention (Years)", "Break-Even Average Object Size (MB)",
    "Break-Even Average Object Size vs. Retention",
    "Smaller objects favor DNA. Open circles mark the last whole-year parity point of ending curves.",
    "S4f_retention_vs_object_size",
    "Average object size giving DNA/Amazon lifecycle parity at each whole-year retention, at fixed archive "
    "volume. Curves connect integer-horizon results for display; a durability replacement can cause a jump. "
    "Cloud cost is A + B/object_size, solved exactly for object size. Missing thresholds and sub-byte "
    "values remain classified in the CSV. Open circles mark the last whole-year parity point before "
    "a curve ends; curves are not extrapolated into horizons without a valid parity threshold. ")
RETENTION_OBJECT = RETENTION_OBJECT.replace(
    "DEPLOYMENT_YEARS = (2050, 2075, 2100)",
    "DEPLOYMENT_YEARS = (2050, 2055, 2060, 2065, 2070, 2075, 2100)")

ARCHIVE_OBJECT = finish('''
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(*ARCHIVE_RANGE_TB)
ax.set_ylim(1e-6, workload_base.archive_size_mb)
archive_sizes = np.geomspace(*ARCHIVE_RANGE_TB, 151)
ax.fill_between(archive_sizes, archive_sizes * 1e6, ARCHIVE_RANGE_TB[1] * 1e6,
                color="0.94", zorder=0)
ax.plot(archive_sizes, archive_sizes * 1e6, color="0.6", ls=":", lw=1)
for index, (year, scenario) in enumerate(zip(DEPLOYMENT_YEARS, scenarios.values())):
    dna, fixed, requests = workload_costs(scenario, scenario.horizon_years)
    gap = dna[-1] - fixed[-1]
    threshold = requests[-1] / gap if gap > 0 and requests[-1] > 0 else np.nan
    points = []
    for archive in archive_sizes:
        valid = np.isfinite(threshold) and threshold <= archive * 1e6
        points.append(threshold if valid else np.nan)
        rows.append({"deployment_year": year, "archive_size_tb": archive,
            "object_size_mb": threshold if valid else np.nan,
            "status": "parity" if valid else "DNA_competitive_at_every_valid_object_size"})
    if not np.isfinite(points).any():
        status_notes.append(f"{year}: DNA competitive for all valid object sizes")
    workload_curve(archive_sizes, points, year, PALETTE[index])
''', "Archive Size (TB)", "Break-Even Average Object Size (MB)",
    "Break-Even Average Object Size vs. Archive Size",
    "Smaller objects favor DNA. Gray area: object exceeds archive. Archive size does not change the cost ratio.",
    "S4g_archive_vs_object_size",
    "Average object size giving DNA/Amazon lifecycle parity as archive volume varies, at fixed baseline "
    "retention. Volume multiplies every cost component, so it cancels from the parity equation and the "
    "object-size threshold is horizontal. The gray area is physically invalid: an object exceeds its "
    "archive. Where a mathematical threshold exceeds the archive size, every valid object size favors DNA. ")


def add_workload_figures(cells, nbf):
    index = next(i for i, c in enumerate(cells)
                 if c.cell_type == "code" and 'record("S4d_' in c.source)
    added = []
    for title, source in (("S4e: retention vs. archive size", RETENTION_ARCHIVE),
                          ("S4f: retention vs. average object size", RETENTION_OBJECT),
                          ("S4g: archive size vs. average object size", ARCHIVE_OBJECT)):
        added.extend([nbf.v4.new_markdown_cell("## " + title + "\n\nUses S2d's 2050 synthesis-price calibration by default. Set `USE_PARITY_CALIBRATION = False` for the original fit. All other price trajectories stay unchanged."),
                      nbf.v4.new_code_cell(source, metadata={"tags": ["figure"]})])
    cells[index + 1:index + 1] = added
