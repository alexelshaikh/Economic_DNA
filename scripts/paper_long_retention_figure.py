"""Supplementary M1 counterpart with a 2,000-year retention period."""

import re


def long_retention_source(m1_source):
    source = re.sub(r"\bBASE\b", "long_retention_base", m1_source)
    source = source.replace("m1", "s12").replace("M1_dna_cost_sensitivity", "S12_dna_cost_sensitivity_long_retention")
    first, body = source.split("\n", 1)
    source = first + "\nlong_retention_base = replace(BASE, horizon_years=2000)\n" + body
    source = source.replace('"Sensitivity of DNA Lifecycle Costs"',
                            '"Sensitivity of DNA Lifecycle Costs: 2,000-Year Retention"')
    audit = '''# Isolate replacement costs, which can be too small to change a float-valued total.
replacement_rows = []
for scenario in s12_scenarios.values():
    for lifetime in (500, 1000, 1500):
        probe = replace(scenario, dna_durability_years=lifetime, technologies=("DNA",))
        annual = simulate_scenario(probe).yearly
        initial_write = float(annual.loc[annual.year == probe.start_year, "write_cost_usd"].sum())
        replacement_write = float(annual.loc[annual.year > probe.start_year, "write_cost_usd"].sum())
        replacement_rows.append({"deployment_year": probe.start_year, "retention_years": probe.horizon_years,
            "dna_lifetime_years": lifetime,
            "replacement_count": (probe.horizon_years - 1) // lifetime,
            "replacement_years": "; ".join(str(y) for y in range(
                probe.start_year + lifetime, probe.start_year + probe.horizon_years, lifetime)),
            "initial_write_usd": initial_write, "replacement_write_usd": replacement_write,
            "replacement_to_initial_write_ratio": replacement_write / initial_write})
s12_replacements = pd.DataFrame(replacement_rows)
print("S12: 2,000-year retention; replacement writes are included, not omitted.")
print(s12_replacements.to_string(index=False, float_format=lambda value: f"{value:.6e}"))
'''
    source = source.replace('record("S12_', audit + '\nrecord("S12_', 1)
    source = source.replace('{"sensitivity": s12}', '{"sensitivity": s12, "replacement_audit": s12_replacements}')
    source = source.replace('Small or overlapping effects are quantified in the accompanying CSV.',
        'Small or overlapping effects are quantified in the accompanying CSV. '
        'Retention is fixed at 2,000 years. DNA lifetimes of 500, 1,000, and 1,500 years '
        'require three, one, and one replacements, respectively; no replacement is charged at '
        'the retention endpoint. Continued exponential synthesis decline makes these future '
        'writes negligible relative to the initial write. This result depends on the absence '
        'of a synthesis-price floor; it does not establish that durability is generally irrelevant.')
    return source


def add_long_retention_figure(cells, nbf):
    original = next(cell.source for cell in cells
                    if cell.cell_type == "code" and 'record("M1_' in cell.source)
    index = next(i for i, cell in enumerate(cells)
                 if cell.cell_type == "code" and 'record("S11_' in cell.source)
    cells[index + 1:index + 1] = [
        nbf.v4.new_markdown_cell("## S12: DNA cost sensitivity with 2,000-year retention\n\n"
            "Supplementary counterpart to M1, retaining its styling and deployment years. "
            "The archive is retained for 2,000 years; baseline DNA lifetime remains 1,000 years "
            "and is varied by +/-50% in its sensitivity row. Replacement writes are included. "
            "The console and exported audit quantify why lifetime still has a negligible cost effect "
            "under the baseline exponential decline without a price floor. Run the two setup cells first; "
            "this cell does not require M1 and does not change BASE."),
        nbf.v4.new_code_cell(long_retention_source(original), metadata={"tags": ["figure"]}),
    ]
