"""Render every figure independently, validating notebook structure and exports."""
from contextlib import redirect_stdout
from io import StringIO
import os
from pathlib import Path
import sys

os.environ["MPLBACKEND"] = "Agg"
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
os.chdir(ROOT)
os.environ["DNA_PAPER_FIGURE_OUTPUT"] = str(ROOT / ".tmp" / "paper_notebook_validation")

import nbformat
import matplotlib.pyplot as plt

notebook = nbformat.read(ROOT / "Figures_Paper_2026.ipynb", as_version=4)
nbformat.validate(notebook)
setups = [c for c in notebook.cells if "setup" in c.metadata.get("tags", [])]
figures = [c for c in notebook.cells if "figure" in c.metadata.get("tags", [])]
assert len(figures) == 40
collected = {name: {} for name in ("FIGURES", "CAPTIONS", "TABLES", "SCENARIOS")}
for index, cell in enumerate(figures, 1):
    namespace = {"__name__": "__main__"}
    with redirect_stdout(StringIO()):
        for setup in setups:
            exec(compile(setup.source, f"setup:{setup.id}", "exec"), namespace)
        namespace["display"] = lambda *args, **kwargs: None
        exec(compile(cell.source, f"figure:{cell.id}", "exec"), namespace)
    assert len(namespace["FIGURES"]) == 1, cell.id
    name, fig = next(iter(namespace["FIGURES"].items()))
    assert fig.axes, name
    for extension in ("pdf", "svg", "png", "json"):
        path = namespace["OUT"] / f"{name}.{extension}"
        assert path.is_file() and path.stat().st_size > 100, path
    assert namespace["TABLES"][name], name
    for registry in collected:
        collected[registry].update(namespace[registry])
    print(f"{index:02d}/{len(figures)} independent render and export: {name}", flush=True)
namespace.update(collected)
with redirect_stdout(StringIO()):
    for cell in notebook.cells:
        if "export" in cell.metadata.get("tags", []):
            exec(compile(cell.source, f"export:{cell.id}", "exec"), namespace)
assert (namespace["OUT"] / "Paper_Figures.pdf").stat().st_size > 1000
# Reuse the original sensitivity notebook's independent numerical comparisons.
namespace["TABLES"] = {f"{name}_{suffix}": frame
                       for name, tables in collected["TABLES"].items()
                       for suffix, frame in tables.items()}
for prefix, name in (
    ("s6", "S6_durability_and_replacements"),
    ("s8", "S8_sequencing_decline_sensitivity"),
    ("s9", "S9_retention_decline_sensitivity"),
    ("s10", "S10_sequencing_parity_by_synthesis"),
    ("s11", "S11_retention_parity_by_synthesis"),
):
    namespace[f"{prefix}_scenarios"] = {
        key: namespace["Scenario"](**{**value, "technologies": tuple(value["technologies"])})
        for key, value in collected["SCENARIOS"][name].items()
    }
    if prefix in ("s10", "s11"):
        namespace[prefix] = collected["TABLES"][name]["parity_boundary"]
original = nbformat.read(ROOT / "Figures_Sensitivity_Paper.ipynb", as_version=4)
checks = next(c.source for c in original.cells if c.cell_type == "code" and c.source.startswith("checks = []"))
exec(compile(checks, "original-sensitivity-numerical-checks", "exec"), namespace)
plt.close("all")
print(f"Validated {len(figures)} independent figure cells, per-figure exports and combined PDF.")
