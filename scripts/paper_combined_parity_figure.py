"""Notebook source for the combined sequencing/retention parity figure."""

SOURCE = r'''SAVE = True  # False previews only; True exports this figure and its data.
FIGURE_SIZE = (10.5, 4.5)  # Width and height in inches, including the shared legend.

sequencing_x = np.linspace(0, 80, 321)
retention_x = np.unique(np.r_[
    np.round(np.geomspace(1, HEATMAP_MAX_RETENTION, 241)).astype(int),
    BASE.horizon_years, BASE.dna_durability_years, BASE.dna_durability_years + 1,
])
retention_x = retention_x[(retention_x >= 1) & (retention_x <= HEATMAP_MAX_RETENTION)]
sequencing_data, sequencing_scenarios = synthesis_parity_lines("sequencing", sequencing_x)
retention_data, retention_scenarios = synthesis_parity_lines("retention", retention_x)

fig = plt.figure(figsize=FIGURE_SIZE, layout="constrained")
grid = fig.add_gridspec(1, 3, width_ratios=(1, 1, 0.48))
axes = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])]
legend_ax = fig.add_subplot(grid[0, 2])
legend_ax.set_axis_off()
for ax, frame, kind, panel in zip(
        axes, (sequencing_data, retention_data), ("sequencing", "retention"), ("a", "b")):
    column = "sequencing_decline_percent" if kind == "sequencing" else "retention_years"
    for index, rate in enumerate(PARITY_SYNTHESIS_RATES):
        part = frame[frame.synthesis_decline_percent == rate].sort_values(column)
        ax.plot(part[column], part.amazon_decline_percent, color=PALETTE[index],
                ls=("-", "--", "-.", ":")[index], lw=2.2, label=f"{rate:g}%")
    baseline = frame[np.isclose(frame.synthesis_decline_percent, BASE.synthesis_decline_percent)]
    if baseline.root_exists.any():
        ax.plot(baseline[column], baseline.amazon_decline_percent, color="0.25",
                ls=(0, (1, 2)), label=f"Baseline: {BASE.synthesis_decline_percent:.2f}%")
    reference_x = BASE.sequencing_decline_percent if kind == "sequencing" else BASE.horizon_years
    ax.scatter(reference_x, BASE.amazon_decline_percent, marker="*", s=140,
               facecolor="white", edgecolor="black", linewidth=1.2, zorder=5)
    annotation = ("Baseline sequencing\nand Amazon declines" if kind == "sequencing"
                  else "Baseline retention\nand Amazon decline")
    ax.annotate(annotation, (reference_x, BASE.amazon_decline_percent),
                xytext=(14, -15) if kind == "sequencing" else (0, -15),
                textcoords="offset points", ha="left" if kind == "sequencing" else "center",
                va="top", fontsize=9, color="0.15", zorder=6,
                arrowprops={"arrowstyle": "->", "color": "black", "lw": 1.1, "shrinkB": 9})
    if kind == "retention":
        ax.set_xscale("log")
        ax.set_xticks([value for value in (1, 10, 100, 1000, 10000) if value <= HEATMAP_MAX_RETENTION])
        ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:,.0f}"))
        ax.xaxis.set_minor_formatter(NullFormatter())
    else:
        ax.set_xticks(np.arange(0, 81, 20))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}%"))
    ax.set_xlim(frame[column].min(), frame[column].max())
    ax.set_ylim(*PARITY_AMAZON_RANGE)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}%"))
    ax.set_xlabel("DNA Sequencing Cost\nAnnual Decline (%)" if kind == "sequencing" else "Retention\n(Years)",
                  fontsize=11, fontweight="bold", color="0.15", labelpad=8)
    ax.set_title(f"({panel}) " + ("Sequencing Decline" if kind == "sequencing" else "Retention"),
                 fontsize=12, fontweight="bold", color="0.15", pad=12)
    ax.tick_params(which="major", labelsize=10, colors="0.15", length=0, pad=7)
    ax.tick_params(which="minor", length=0)
    ax.grid(False, which="both")
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("0.15")
        spine.set_linewidth(1.1)
axes[0].set_ylabel("Amazon Deep Archive Cost\nAnnual Decline (%)", fontsize=11,
                   fontweight="bold", color="0.15", labelpad=8)
handles, labels = [], []
for ax in axes:
    for handle, label in zip(*ax.get_legend_handles_labels()):
        if label not in labels:
            handles.append(handle)
            labels.append(label)
shared_legend = legend_ax.legend(handles, labels, loc="center left",
    title="DNA Synthesis Cost\nAnnual Decline", frameon=False, fontsize=10,
    title_fontsize=10, handlelength=2.6, labelspacing=1.0, borderaxespad=0)
shared_legend.get_title().set_fontweight("bold")
fig.suptitle(f"DNA / Amazon Cost Parity: Deployment Year $y_0 = {FUTURE_START}$",
             fontsize=14, fontweight="bold", color="0.15")
record("S10_S11_combined_parity", fig,
       f"DNA/Amazon lifecycle parity for deployment in {FUTURE_START}. "
       f"(a) Sequencing decline varies at {BASE.horizon_years}-year retention. "
       f"(b) Retention varies at the fitted sequencing decline of {BASE.sequencing_decline_percent:.2f}%/year. "
       "The shared legend identifies annual synthesis declines. Below a curve DNA is cheaper; "
       "above it Amazon is cheaper. Missing segments have no parity within the displayed range. "
       "Stars mark baseline sequencing or retention and baseline Amazon decline, not a synthesis scenario. "
       "The baseline synthesis trajectory is evaluated in both panels. All price anchors are unchanged.",
       {"sequencing_parity": sequencing_data, "retention_parity": retention_data},
       {**{f"sequencing_{key}": value for key, value in sequencing_scenarios.items()},
        **{f"retention_{key}": value for key, value in retention_scenarios.items()}})
'''
