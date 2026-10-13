"""Self-contained export function embedded in the final paper notebook cell."""

import inspect


def export_results_figures(
    output_dir=r"C:\Users\Alex\papers\An_Economic_Analysis_of_DNA_based_Data_Storage_Systems\figures\results",
    figure_sizes=None,
    default_size=None,
    notebook_path=None,
):
    """Regenerate the Results figures from saved notebook cells; export PDF/EPS.

    Sizes are (width, height) in inches. None retains each figure's own size;
    figure_sizes (e.g. {"F1": (9, 6)}) overrides default_size for individual plots.
    Exports use tight bounding boxes, so final page dimensions include the labels
    and may differ from the requested canvas size. Save notebook edits first.
    This reruns saved setup defaults, not unsaved in-memory parameter changes.
    Existing selected exports are replaced; other files are left untouched.
    """
    import ast
    import json
    import math
    import os
    import sys
    import tempfile
    from pathlib import Path
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import pandas as pd
    from IPython.utils.capture import capture_output

    selected = {
        "F1": "F1_dna_trajectories",
        "F2a": "F2a_unit_costs", "F2b": "F2b_unit_costs", "F2c": "F2c_unit_costs",
        "F3a": "F3a_lifecycle_components", "F3b": "F3b_lifecycle_components",
        "F3c": "F3c_lifecycle_components", "M1": "M1_dna_cost_sensitivity",
        "S2": "S2_retrieval_vs_synthesis_threshold", "S4": "S4_joint_decline_sensitivity",
        "S4c": "S4c_deployment_year_parity", "S4d": "S4d_dna_declines_deployment_parity",
        "S10": "S10_sequencing_parity_by_synthesis", "S11": "S11_retention_parity_by_synthesis",
        "S10_S11": "S10_S11_combined_parity",
    }
    figure_sizes = {} if figure_sizes is None else dict(figure_sizes)
    unknown = set(figure_sizes) - selected.keys()
    if unknown:
        raise ValueError(f"Unknown figure size keys: {sorted(unknown)}")

    def checked_size(size):
        if size is None:
            return None
        try:
            width, height = map(float, size)
        except (TypeError, ValueError) as error:
            raise ValueError("Figure sizes must be (width, height) pairs in inches.") from error
        if not all(math.isfinite(v) and v > 0 for v in (width, height)):
            raise ValueError("Figure dimensions must be finite and positive.")
        return width, height

    default_size = checked_size(default_size)
    figure_sizes = {key: checked_size(value) for key, value in figure_sizes.items()}
    if notebook_path is None:
        root = next((p for p in (Path.cwd(), *Path.cwd().parents)
                     if (p / "Figures_Paper_2026.ipynb").is_file()), None)
        if root is None:
            raise FileNotFoundError("Run from the Economic_DNA repository or pass notebook_path.")
        notebook_path = root / "Figures_Paper_2026.ipynb"
    notebook_path = Path(notebook_path).resolve()
    notebook = json.loads(notebook_path.read_text(encoding="utf-8"))
    setups, sources = [], {}
    for cell in notebook["cells"]:
        if cell["cell_type"] != "code":
            continue
        source = "".join(cell["source"])
        if "setup" in cell.get("metadata", {}).get("tags", []):
            setups.append(source)
        if "figure" not in cell.get("metadata", {}).get("tags", []):
            continue
        tree = ast.parse(source)
        names = [node.args[0].value for node in ast.walk(tree)
                 if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                 and node.func.id == "record" and node.args
                 and isinstance(node.args[0], ast.Constant)]
        for name in names:
            if name in selected.values():
                if name in sources:
                    raise ValueError(f"Duplicate figure cell: {name}")
                # Disable normal PNG/CSV exports; only the requested PDF/EPS files are written.
                for node in ast.walk(tree):
                    if isinstance(node, ast.Assign) and any(
                        isinstance(target, ast.Name) and target.id == "SAVE" for target in node.targets
                    ):
                        node.value = ast.Constant(False)
                sources[name] = compile(ast.fix_missing_locations(tree), str(notebook_path), "exec")
    missing = set(selected.values()) - sources.keys()
    if missing or not setups:
        raise ValueError(f"Missing setup or figure cells: {sorted(missing)}")

    destination = Path(output_dir).expanduser().resolve()
    destination.mkdir(parents=True, exist_ok=True)
    namespace = {"__name__": "__paper_results_export__"}
    original_cwd, original_path = Path.cwd(), list(sys.path)
    existing_figures = set(plt.get_fignums())
    rows = []
    try:
        os.chdir(notebook_path.parent)
        with mpl.rc_context(), plt.ioff(), tempfile.TemporaryDirectory(
            prefix=".results-export-", dir=destination
        ) as staging:
            with capture_output():
                for source in setups:
                    exec(compile(source, str(notebook_path), "exec"), namespace)
            namespace["display"] = lambda *args, **kwargs: None
            for label, stem in selected.items():
                print(f"Regenerating {label}...")
                with capture_output():
                    exec(sources[stem], namespace)
                fig = namespace["FIGURES"][stem]
                size = figure_sizes.get(label, default_size)
                if size is not None:
                    fig.set_size_inches(*size, forward=True)
                fig.savefig(Path(staging) / f"{label}.pdf", bbox_inches="tight", facecolor="white")
                # EPS cannot encode alpha: use opaque artwork, retaining hidden artists.
                translucent = [(artist, artist.get_alpha()) for artist in fig.findobj()
                               if artist.get_alpha() is not None and 0 < artist.get_alpha() < 1]
                try:
                    for artist, _ in translucent:
                        artist.set_alpha(1)
                    fig.savefig(Path(staging) / f"{label}.eps", bbox_inches="tight", facecolor="white")
                finally:
                    for artist, alpha in translucent:
                        artist.set_alpha(alpha)
                width, height = fig.get_size_inches()
                rows.append({"Figure": label, "Width (in)": width, "Height (in)": height,
                             "PDF": str(destination / f"{label}.pdf"),
                             "EPS": str(destination / f"{label}.eps")})
                plt.close(fig)
            # Publish only after every figure has regenerated and saved successfully.
            for label in selected:
                for extension in ("pdf", "eps"):
                    filename = f"{label}.{extension}"
                    os.replace(Path(staging) / filename, destination / filename)
    finally:
        os.chdir(original_cwd)
        sys.path[:] = original_path
        for number in set(plt.get_fignums()) - existing_figures:
            plt.close(number)
    print(f"Exported {len(rows)} figures as PDF and EPS to {destination}")
    return pd.DataFrame(rows)


def add_results_export(cells, nbf):
    cells.extend([
        nbf.v4.new_markdown_cell("## Export Results figures to the paper\n\n"
            "Save this notebook first. Run the final cell to define `export_results_figures()`, then "
            "uncomment one of the example calls. It regenerates the selected figures from the saved "
            "setup and figure cells, without requiring earlier cells to have run. Existing `F1.pdf`, "
            "`F1.eps`, etc. in the paper's `figures/results` folder are replaced. Other files are untouched.\n\n"
            "Sizes are canvas dimensions in inches; tight cropping can change final page dimensions. "
            "Omit sizes to keep each figure's defaults. EPS uses opaque styling where PDF supports transparency."),
        nbf.v4.new_code_cell(inspect.getsource(export_results_figures) + '''

# Keep every figure's default size:
# exported_results = export_results_figures()

# Override selected figures only (width, height in inches):
# exported_results = export_results_figures(figure_sizes={"F1": (9, 6), "M1": (16, 5.7)})

# Or set one canvas size for all figures, with optional individual overrides:
# exported_results = export_results_figures(default_size=(10.5, 6.7), figure_sizes={"M1": (16, 5.7)})
''', metadata={"tags": ["results-export"]}),
    ])
