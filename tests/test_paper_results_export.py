import ast
import json
from pathlib import Path
import tempfile
import unittest

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scripts.paper_results_export import export_results_figures

ROOT = Path(__file__).resolve().parents[1]


class ResultsExportTests(unittest.TestCase):
    def test_invalid_sizes_do_not_write(self):
        for kwargs in ({"figure_sizes": {"typo": (9, 6)}},
                       {"default_size": (0, 6)}, {"default_size": (float("nan"), 6)},
                       {"figure_sizes": {"F1": (9,)}}):
            with self.assertRaises(ValueError):
                export_results_figures(**kwargs)

    def test_notebook_export_regenerates_all_selected_figures(self):
        notebook = json.loads((ROOT / "Figures_Paper_2026.ipynb").read_text(encoding="utf-8"))
        cell = notebook["cells"][-1]
        self.assertIn("results-export", cell["metadata"]["tags"])
        ns = {}
        tree = ast.parse("".join(cell["source"]))
        function = next(node for node in tree.body
                        if isinstance(node, ast.FunctionDef) and node.name == "export_results_figures")
        # Notebook examples may be enabled by the author; never run them during tests.
        exec(compile(ast.Module(body=[function], type_ignores=[]), "<export-test>", "exec"), ns)
        sentinel = plt.figure()
        original_font = matplotlib.rcParams["font.size"]
        try:
            with tempfile.TemporaryDirectory(dir=ROOT) as folder:
                untouched = Path(folder) / "unrelated.txt"
                untouched.write_text("keep", encoding="utf-8")
                (Path(folder) / "F1.pdf").write_bytes(b"stale")
                frame = ns["export_results_figures"](
                    output_dir=folder, notebook_path=ROOT / "Figures_Paper_2026.ipynb",
                    default_size=(10.5, 6.7), figure_sizes={"F1": (9, 6), "M1": (16, 5.7)})
                self.assertEqual(frame.Figure.tolist(), [
                    "F1", "F2a", "F2b", "F2c", "F3a", "F3b", "F3c", "M1",
                    "S2", "S4", "S4c", "S4d", "S10", "S11", "S10_S11"])
                self.assertEqual(len(list(Path(folder).glob("*.pdf"))), 15)
                self.assertEqual(len(list(Path(folder).glob("*.eps"))), 15)
                self.assertEqual(len(list(Path(folder).iterdir())), 31)
                self.assertEqual(untouched.read_text(encoding="utf-8"), "keep")
                for row in frame.to_dict("records"):
                    self.assertTrue(Path(row["PDF"]).read_bytes().startswith(b"%PDF-"))
                    self.assertTrue(Path(row["EPS"]).read_bytes().startswith(b"%!PS-Adobe"))
                    expected = (9, 6) if row["Figure"] == "F1" else (16, 5.7) if row["Figure"] == "M1" else (10.5, 6.7)
                    self.assertEqual((row["Width (in)"], row["Height (in)"]), expected)
            self.assertTrue(plt.fignum_exists(sentinel.number))
            self.assertEqual(matplotlib.rcParams["font.size"], original_font)
        finally:
            plt.close(sentinel)
