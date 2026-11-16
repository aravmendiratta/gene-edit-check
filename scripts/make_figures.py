"""Re-draw the README figures from a saved run. Needs no API key and makes no API calls.

    python scripts/make_figures.py [outputs/report_crc-pd1-exhaustion.json]
"""

from __future__ import annotations

import sys
from pathlib import Path

from gene_edit_check.models import OptimizationRunResult
from gene_edit_check.ui.plots import export_figure_set, export_optimization_figure

report = Path(sys.argv[1] if len(sys.argv) > 1 else "outputs/report_crc-pd1-exhaustion.json")
result = OptimizationRunResult.model_validate_json(report.read_text(encoding="utf-8"))
out = report.parent
for p in export_figure_set(result, out):
    print("wrote", p)
print("wrote", export_optimization_figure(result, out / f"figure_{report.stem.removeprefix('report_')}.png"))
