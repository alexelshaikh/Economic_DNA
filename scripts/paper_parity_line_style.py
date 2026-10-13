"""Convert S4b-d to directly labeled parity curves without heatmaps."""


def line_only_prices(source):
    start = source.index('mesh = ax.pcolormesh(')
    end = source.index('if read_coefficient > 0:', start)
    source = source[:start] + source[end:]
    source = source.replace('    ax.plot(parity_x, parity_y, color="0.1", lw=1.8, label="Cost parity")',
        '''    from matplotlib.contour import ContourSet
    parity_contour = ContourSet(ax, [0], [[np.column_stack([parity_x, parity_y])]],
                               colors=["0.1"], linewidths=1.8)
    midpoint = len(parity_x) // 2
    ax.clabel(parity_contour, fmt={0: "Cost parity"}, inline=True, inline_spacing=5,
              fontsize=10.5, manual=[(parity_x[midpoint], parity_y[midpoint])])''')
    start = source.index('colorbar = fig.colorbar(')
    end = source.index('s4b_reference =', start)
    source = source[:start] + source[end:]
    source = source.replace('The heatmap shows the linear lifecycle cost ratio, centered on one (parity); unlike S4 it is not a log ratio. ',
        'Only the parity boundary is drawn; combinations below it favor DNA and combinations above it favor Amazon. ')
    return source


def line_only_deployments(source, kind):
    source = expand_deployment_search(source, kind)
    start = source.index('        mesh = ax.pcolormesh(')
    end = source.index('        grid_x, grid_y =', start)
    source = source[:start] + source[end:]
    start = source.index('    label = f"{year} (first viable)"')
    end = source.index('ax.set_xlim(', start)
    source = source[:start] + source[end:]
    start = source.index('colorbar = fig.colorbar(')
    end = source.index(f'{kind}_parity =', start)
    source = source[:start] + source[end:]
    start = source.index(f'{kind}_audit =')
    end = source.index('\n', start)
    styling = '''
from matplotlib.contour import ContourSet
visible = PARITY_TABLE.synthesis_decline_percent.dropna()
if len(visible):
    left = max(SYNTHESIS_RANGE[0], np.floor((visible.min() - 2) / 5) * 5)
    ax.set_xlim(left, SYNTHESIS_RANGE[1])
for index, year in enumerate(deployment_years):
    curve = PARITY_TABLE[PARITY_TABLE.deployment_year == year]
    if curve.empty:
        continue
    points = curve[["synthesis_decline_percent", "Y_COLUMN"]].to_numpy()
    color = PALETTE[index]
    parity_contour = ContourSet(ax, [0], [[points]], colors=[color], linewidths=1.8)
    target = points[len(points) // 2]
    # Very short edge segments cannot hold an inline label; keep the year beside the segment.
    relative_height = np.ptp(points[:, 1]) / (ax.get_ylim()[1] - ax.get_ylim()[0])
    if relative_height < 0.08:
        near_top = target[1] > np.mean(ax.get_ylim())
        ax.annotate(str(year), target, xytext=(-35, -16 if near_top else 12), textcoords="offset points",
                    fontsize=10.5, color=color, ha="right", va="top" if near_top else "bottom",
                    arrowprops={"arrowstyle": "-", "color": color, "lw": 0.9})
    else:
        ax.clabel(parity_contour, fmt={0: str(year)}, inline=True, inline_spacing=5,
                  fontsize=10.5, manual=[tuple(target)])
'''.replace('PARITY_TABLE', f'{kind}_parity').replace('Y_COLUMN',
        'amazon_decline_percent' if kind == 's4c' else 'sequencing_decline_percent')
    source = source[:end] + '\n' + styling + source[end:]
    source = source.replace('The heatmap applies ONLY to the earliest deployment year, not all curves. Price anchors, ',
                            'Each directly labeled line refers to its deployment year. Price anchors, ')
    # The reference grid remains available as numerical data, but is no longer drawn.
    source = source.replace('"The heatmap applies ONLY to the earliest "\n'
                            '       "deployment year, not to every contour: negative log ratios favor DNA. All original price anchors, "',
                            '"Each line is labeled with its deployment year. All original price anchors, "')
    source = consistent_deployment_labels(source)
    return refine_s4d_labels(source) if kind == "s4d" else source


def refine_s4d_labels(source):
    source = preserve_s4d_endpoints(source)
    source = omit_2027_callout(source)
    source = source.replace('xytext=(-35, -16 if near_top else 12)',
                            'xytext=(-15, -30 if near_top else 12)')
    source = source.replace('parity_contour = ContourSet(ax, [0], [[points]], colors=[color], linewidths=1.8)',
        'parity_contour = ContourSet(ax, [0], [[points]], colors=[color], linewidths=1.8, zorder=5)')
    source = source.replace('    if relative_height < 0.08:\n', '''    if relative_height < 0.08:
        # The earliest curve lies almost on the 100% spine; mark its true location visibly.
        ax.scatter(*target, s=42, facecolors="white", edgecolors=[color], linewidths=1.5,
                   zorder=6, clip_on=False)
''')
    source = source.replace('        ax.clabel(parity_contour, fmt={0: f"$y_0 = {year}$"}',
                            '        inline_labels = ax.clabel(parity_contour, fmt={0: f"$y_0 = {year}$"}')
    source = source.replace('fontsize=10.5, manual=[tuple(target)])', '''fontsize=10.5, manual=[tuple(target)])
        for label in inline_labels:
            if label.get_rotation() > 180:
                label.set_rotation((label.get_rotation() + 180) % 360)''')
    return source


def omit_2027_callout(source):
    return source.replace(
        '    relative_height = np.ptp(points[:, 1]) / (ax.get_ylim()[1] - ax.get_ylim()[0])',
        '''    if year == 2027:
        continue  # Its near-100% location is indistinguishable from the boundary at this scale.
    relative_height = np.ptp(points[:, 1]) / (ax.get_ylim()[1] - ax.get_ylim()[0])''')


def preserve_s4d_endpoints(source):
    return source.replace(
        '        root = parity_root(lambda synthesis: dna_total(synthesis, rate) - amazon_total, *SYNTHESIS_RANGE)',
        '''        # Preserve the known edge coordinate: re-solving a nearly flat cost function loses it.
        edge_matches = [edge for edge, edge_rate in zip(SYNTHESIS_RANGE, edge_rates)
                        if np.isfinite(edge_rate) and rate == edge_rate]
        root = edge_matches[-1] if edge_matches else parity_root(
            lambda synthesis: dna_total(synthesis, rate) - amazon_total, *SYNTHESIS_RANGE)''')


def consistent_deployment_labels(source):
    source = source.replace('DNA / Amazon Cost Parity by Deployment Year"',
                            'DNA / Amazon Cost Parity by Deployment Year $y_0$"')
    source = source.replace('DNA Synthesis Annual Decline (%)', 'DNA Synthesis Cost Annual Decline (%)')
    source = source.replace('Annual Decline in DNA Sequencing Cost per MB (%)',
                            'DNA Sequencing Cost Annual Decline (%)')
    source = source.replace('Annual Decline in Amazon Deep Archive Prices (%)\\nStorage, Upload and Retrieval',
                            'Amazon Deep Archive Cost Annual Decline (%)\\n(Storage, Upload and Retrieval)')
    source = source.replace('ax.annotate(str(year), target', 'ax.annotate(f"$y_0 = {year}$", target')
    source = source.replace('fmt={0: str(year)}', 'fmt={0: f"$y_0 = {year}$"}')
    return source


def expand_deployment_search(source, kind):
    source = source.replace('SYNTHESIS_RANGE = (0.0, 70.0)  # Same annual decline ranges as S4, in percent.',
        'SYNTHESIS_RANGE = (0.0, float(np.nextafter(100.0, 0.0)))  # All nonnegative declines below 100%.')
    source = source.replace('SEQUENCING_RANGE = (0.0, 80.0)',
        'SEQUENCING_RANGE = (0.0, float(np.nextafter(100.0, 0.0)))')
    source = source.replace('Amazon Annual Cost Decline (%)',
        'Annual Decline in Amazon Deep Archive Prices (%)\\nStorage, Upload and Retrieval')
    source = source.replace('DNA Sequencing Annual Decline (%)',
        'Annual Decline in DNA Sequencing Cost per MB (%)')
    source = source.replace('within synthesis declines {SYNTHESIS_RANGE[0]:g}-{SYNTHESIS_RANGE[1]:g}%/year and Amazon declines ',
        'within synthesis declines from 0% up to (but not including) 100%/year and Amazon declines ')
    source = source.replace('%/year, as in S4. This is not the first viable year under ',
        '%/year. This is not the first viable year under ')
    source = source.replace('f"{SYNTHESIS_RANGE[0]:g}-{SYNTHESIS_RANGE[1]:g}%/year and sequencing declines "\n'
        '       f"{SEQUENCING_RANGE[0]:g}-{SEQUENCING_RANGE[1]:g}%/year. Amazon\'s annual decline remains "',
        '"from 0% up to (but not including) 100%/year for both DNA price components. Amazon\'s annual decline remains "')
    source = source.replace("First deployment with parity within S4's decline ranges:",
        "First deployment with parity allowing synthesis declines below 100%:")
    return source
