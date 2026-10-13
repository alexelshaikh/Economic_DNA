"""Presentation-only styling for the S10 and S11 parity figures."""

STYLE = '''
ax = fig.axes[0]
fig.set_size_inches(10.5, 6.7)
ax.grid(False, which="both")
ax.set_title("", loc="left")
ax.set_title(f"DNA / Amazon Cost Parity: Deployment Year $y_0 = {FUTURE_START}$",
             fontsize=14, fontweight="bold", color="0.15", pad=16)
ax.set_xlabel("DNA Sequencing Cost Annual Decline (%)",
              fontsize=12, fontweight="bold", color="0.15", labelpad=8)
ax.set_ylabel("Amazon Deep Archive Cost Annual Decline (%)",
              fontsize=12, fontweight="bold", color="0.15", labelpad=8)
ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}%"))
ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:g}%"))
ax.tick_params(which="major", labelsize=10.5, colors="0.15", length=0, pad=7)
ax.tick_params(which="minor", length=0)
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_color("0.15")
    spine.set_linewidth(1.1)
baseline_note = ""
for text in list(ax.texts):
    if text.get_text().startswith("Baseline synthesis ("):
        baseline_note = text.get_text().replace("\\n", " ")
        text.remove()
    else:
        text.set_fontsize(10)
        text.set_color("0.15")
handles, labels = ax.get_legend_handles_labels()
ax.get_legend().remove()
paper_legend = fig.legend(handles, [label.replace("%/year", "%") for label in labels],
    title="DNA Synthesis Cost Annual Decline", loc="outside lower center",
    ncol=4, frameon=True, fancybox=True, facecolor="white", edgecolor="0.8",
    framealpha=0.9, fontsize=10.5, title_fontsize=11, handlelength=3,
    handletextpad=0.8, columnspacing=1.8, borderpad=0.7)
paper_legend.get_title().set_fontweight("bold")
from matplotlib.offsetbox import TextArea
legend_note = "Below a curve: DNA cheaper. Above a curve: Amazon cheaper."
if baseline_note:
    legend_note += "\\n" + baseline_note
note_box = TextArea(legend_note,
    textprops={"fontsize": 9, "fontstyle": "italic", "color": "0.3", "ha": "center"})
note_box.set_figure(fig)
paper_legend._legend_box.get_children().insert(0, note_box)
'''


def style_s10(source):
    marker = 'record("S10_sequencing_parity_by_synthesis"'
    if STYLE.strip() in source:
        return source
    return source.replace(marker, STYLE + "\n" + marker, 1)


def style_s11(source):
    styling = STYLE.replace("DNA Sequencing Cost Annual Decline (%)", "Retention (Years)")
    styling = styling.replace(
        'ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}%"))',
        'ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:,.0f}"))\n'
        'ax.xaxis.set_minor_formatter(NullFormatter())')
    if styling.strip() in source or 'paper_legend = fig.legend(' in source:
        return source
    return source.replace('record("S11_retention_parity_by_synthesis"',
                          styling + '\nrecord("S11_retention_parity_by_synthesis"', 1)
