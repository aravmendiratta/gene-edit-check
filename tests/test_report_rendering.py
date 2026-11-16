"""Smoke tests for the figure and terminal report.

These are rendering-only paths that the unit tests never touched, so a typo in an f-string
reached the CLI instead of the suite. Rendering a real saved run catches that.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from rich.console import Console

from crispr_reflexion.models import OptimizationRunResult
from crispr_reflexion.ui.plots import export_optimization_figure
from crispr_reflexion.ui.terminal_hud import render_optimization_report

REPORT = Path(__file__).resolve().parents[1] / "outputs" / "report_crc-pd1-exhaustion.json"


@pytest.fixture(scope="module")
def saved_run() -> OptimizationRunResult:
    if not REPORT.exists():
        pytest.skip(f"{REPORT.name} not present (run: crispr-reflexion optimize)")
    return OptimizationRunResult.model_validate_json(REPORT.read_text(encoding="utf-8"))


def test_figure_export_renders_every_panel(saved_run: OptimizationRunResult, tmp_path: Path) -> None:
    out = export_optimization_figure(saved_run, output_path=tmp_path / "fig.png")
    assert out.exists() and out.stat().st_size > 50_000  # a real 4-panel figure, not a blank canvas


def test_terminal_report_renders_without_raising(saved_run: OptimizationRunResult) -> None:
    console = Console(file=open(__import__("os").devnull, "w", encoding="utf-8"), width=200)
    render_optimization_report(saved_run, console=console)


def test_saved_run_carries_the_new_analyses(saved_run: OptimizationRunResult) -> None:
    """Guards the report schema the README and the published write-up quote."""
    assert saved_run.ism_prescreen is not None
    assert saved_run.hazard_coupling.get("verdict")
    assert saved_run.hazard_coupling.get("coupling_class")
    t2 = saved_run.iterations[0].report.tier2_hazards
    assert t2.limiting_evidence_tier in {"observed", "imputed_supported", "imputed_exploratory"}
    assert t2.reference_biosample_name


def test_figure_handles_a_run_with_no_coupling_analysis(saved_run: OptimizationRunResult, tmp_path: Path) -> None:
    stripped = saved_run.model_copy(deep=True)
    stripped.hazard_coupling = {}
    stripped.efficacy_null = None
    out = export_optimization_figure(stripped, output_path=tmp_path / "bare.png")
    assert out.exists()
