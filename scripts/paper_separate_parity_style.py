"""Compact standalone parity plots for the four-panel Results figure."""

import ast


def simplify_s2_title(source):
    tree = ast.parse(source)
    call = next(node for node in ast.walk(tree) if isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute) and node.func.attr == "suptitle")
    lines = source.splitlines(keepends=True)
    lines[call.lineno - 1:call.end_lineno] = [
        'fig.suptitle("Break-even Synthesis Price by Data Retrieval Rate",\n',
        '             fontsize=14, fontweight="bold", color="0.15")\n',
    ]
    return "".join(lines)


def compact_parity_legend(source):
    if 'from matplotlib.offsetbox import TextArea' in source:
        start = source.index('from matplotlib.offsetbox import TextArea')
        end = source.index('\nrecord(', start)
        source = source[:start] + source[end:]
    for size in ('10.5, 6.7', '4.4, 5.2', '4.4, 4.7'):
        source = source.replace(f'fig.set_size_inches({size})', 'fig.set_size_inches(5.8, 4.4)')
    source = source.replace('ncol=4, frameon=True', 'ncol=1, frameon=True')
    source = source.replace('ncol=2, frameon=True', 'ncol=1, frameon=True')
    source = source.replace('title="DNA Synthesis Cost Annual Decline"',
                            'title="Synthesis decline"')
    source = source.replace('title="DNA Synthesis Cost\\nAnnual Decline"',
                            'title="Synthesis decline"')
    source = source.replace('loc="outside lower center"', 'loc="outside right center"')
    source = source.replace('handlelength=3', 'handlelength=2')
    if 'record("S10_' in source:
        source = source.replace('text.set_position((0, -13))', 'text.set_position((14, -15))')
        source = source.replace('text.set_ha("center")', 'text.set_ha("left")')
    return source
