"""Figures (`ui/plots.py`).

Clean, light, print-friendly charts. The design follows a few fixed rules:

- one blue accent for the thing being discussed, neutral grey for context, orange only where
  something fails or breaks a limit;
- thin marks, hairline solid gridlines, no chart borders, no dual axes;
- direct labels on the few points that matter, never a number on every point;
- pass/fail is carried by shape (tick / cross) as well as colour.

`export_figure_set` writes the individual figures used in the README; `export_optimization_figure`
writes one combined summary figure and is what the CLI calls.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from gene_edit_check.models import OptimizationRunResult

# ---- design tokens (light surface) ---------------------------------------------------------
SURFACE = "#ffffff"
INK = "#111111"        # primary text
INK_2 = "#55534f"      # secondary text
INK_3 = "#8a8985"      # muted text
GRID = "#ebebe8"       # hairlines
NEUTRAL = "#b9b8b3"    # context marks
NEUTRAL_SOFT = "#e4e3df"
BLUE = "#2a78d6"       # accent / the subject
BLUE_SOFT = "#cde2fb"
ORANGE = "#eb6834"     # fail / limit
ORANGE_SOFT = "#fbe1d6"

_FONT = ["Inter", "SF Pro Display", "Segoe UI", "Helvetica Neue", "Arial", "DejaVu Sans"]
_DPI = 200


def _theme() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": _FONT,
            "figure.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
            "savefig.facecolor": SURFACE,
            "text.color": INK,
            "axes.edgecolor": GRID,
            "axes.labelcolor": INK_2,
            "xtick.color": INK_3,
            "ytick.color": INK_3,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
            "axes.grid": False,
            "xtick.major.size": 0,
            "ytick.major.size": 0,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "axes.labelsize": 9.5,
            "legend.frameon": False,
        }
    )


def _clean_axes(ax: Any, grid_axis: str = "y") -> None:
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.grid(True, axis=grid_axis, color=GRID, linewidth=0.8, linestyle="-")
    ax.set_axisbelow(True)


def _title(fig: Any, title: str, subtitle: str = "", x: float = 0.04) -> None:
    """Left-aligned title block. Spacing is set in inches so it holds at any figure height."""
    h = fig.get_figheight()
    top = 1 - 0.22 / h
    fig.text(x, top, title, fontsize=15, fontweight="semibold", color=INK, va="top", ha="left")
    if subtitle:
        fig.text(x, top - 0.36 / h, subtitle, fontsize=9.8, color=INK_2, va="top", ha="left")


def _footer(fig: Any, result: OptimizationRunResult | None, extra: str = "") -> None:
    parts = []
    if result is not None:
        prov = result.provenance
        parts.append(f"{prov.model_version}")
        parts.append(f"seed {prov.random_seed}")
    parts.append("model predictions, not lab measurements")
    if extra:
        parts.append(extra)
    fig.text(0.04, 0.012, "  ·  ".join(parts), fontsize=7.6, color=INK_3, ha="left", va="bottom")


# ---- data helpers --------------------------------------------------------------------------

_CHECKS = [
    ("Lowers the gene\nin the target cells", ("TIER1_DIRECTION", "TIER1_EFFICACY", "TIER1_GOODHART")),
    ("Leaves healthy\ncells alone", ("TIER1_SELECTIVITY",)),
    ("No splicing\ndamage", ("TIER2_SPLICE",)),
    ("Methylation\nstays stable", ("TIER2_METHYLATION",)),
    ("Neighbouring\ngenes untouched", ("TIER1_NEIGHBOR",)),
    ("No unwanted\nbystander edits", ("TIER2_BYSTANDER",)),
]


def _check_row(reasons: list[str]) -> list[bool]:
    """True = passed. One entry per column in `_CHECKS`."""
    return [not any(r.startswith(p) for r in reasons for p in prefixes) for _, prefixes in _CHECKS]


def _candidate_label(step: Any) -> tuple[str, str]:
    cand = step.report.candidate
    v = cand.target_edit
    editor = getattr(cand.editor, "value", str(cand.editor))
    return f"{v.chromosome}:{v.position:,}", f"{v.reference} to {v.alternate}  ·  {editor}"


# ---- panels --------------------------------------------------------------------------------


def _draw_tick(ax: Any, x: float, y: float, color: str, s: float = 0.16) -> None:
    ax.plot([x - s, x - s * 0.25, x + s], [y - s * 0.05, y + s * 0.65, y - s * 0.75],
            color=color, lw=1.8, solid_capstyle="round", solid_joinstyle="round", zorder=4)


def _draw_cross(ax: Any, x: float, y: float, color: str, s: float = 0.15) -> None:
    ax.plot([x - s, x + s], [y - s, y + s], color=color, lw=1.8, solid_capstyle="round", zorder=4)
    ax.plot([x - s, x + s], [y + s, y - s], color=color, lw=1.8, solid_capstyle="round", zorder=4)


def _panel_scorecard(ax: Any, result: OptimizationRunResult) -> None:
    steps = sorted(
        result.iterations,
        key=lambda s: (len(s.report.failure_reasons), abs(s.report.tier2_hazards.cpg_methylation_mean_abs_shift)),
    )
    n, m = len(steps), len(_CHECKS)
    ax.set_xlim(-2.9, m - 0.4)
    ax.set_ylim(n + 0.2, -1.25)
    ax.axis("off")

    fail_counts = [0] * m
    for i, step in enumerate(steps):
        row = _check_row(step.report.failure_reasons)
        if i % 2 == 0:
            ax.add_patch(FancyBboxPatch((-2.95, i - 0.5), m + 2.5, 1.0, boxstyle="round,pad=0,rounding_size=0.12",
                                        fc="#f7f7f5", ec="none", zorder=0))
        pos, sub = _candidate_label(step)
        ax.text(-2.8, i - 0.1, pos, fontsize=9.3, color=INK, va="center", ha="left", fontweight="semibold")
        ax.text(-2.8, i + 0.26, sub, fontsize=7.8, color=INK_3, va="center", ha="left")
        for j, ok in enumerate(row):
            if ok:
                _draw_tick(ax, j, i, NEUTRAL)
            else:
                fail_counts[j] += 1
                ax.scatter([j], [i], s=330, color=ORANGE_SOFT, zorder=2, linewidths=0)
                _draw_cross(ax, j, i, ORANGE)

    for j, (label, _) in enumerate(_CHECKS):
        ax.text(j, -0.8, label, fontsize=8.4, color=INK_2, va="center", ha="center", linespacing=1.25)
    # Column totals: the bar of failures is the result.
    ax.plot([-2.9, m - 0.45], [n - 0.35, n - 0.35], color=GRID, lw=1.0)
    ax.text(-2.8, n + 0.02, "Candidates that fail", fontsize=8.6, color=INK_2, va="center", ha="left")
    worst = max(fail_counts) if fail_counts else 0
    for j, c in enumerate(fail_counts):
        hot = c == worst and c > 0
        ax.text(j, n + 0.02, f"{c} of {n}", fontsize=10.5 if hot else 9, va="center", ha="center",
                color=ORANGE if hot else INK_2, fontweight="semibold" if hot else "normal")


def _panel_coupling(ax: Any, result: OptimizationRunResult, hazard: str = "cpg_methylation") -> None:
    null = result.efficacy_null
    cand_x = np.array([abs(s.report.tier1_selectivity.target_delta_sketch_b) for s in result.iterations])
    if hazard == "cpg_methylation":
        cand_y = np.array([s.report.tier2_hazards.cpg_methylation_mean_abs_shift for s in result.iterations])
        null_y = np.array(null.cpg_methylation_shift_values) if null else np.array([])
        limit, limit_label, ylabel = 0.12, "safety limit 0.12", "Predicted change in CpG methylation"
    else:
        cand_y = np.array([s.report.tier2_hazards.max_splice_junction_delta for s in result.iterations])
        null_y = np.array(null.splice_delta_values) if null else np.array([])
        limit, limit_label, ylabel = 0.15, "safety limit 0.15", "Predicted change in splicing"
    null_x = np.array(null.abs_lfc_values) if null else np.array([])

    _clean_axes(ax, "y")
    if len(null_x) and len(null_x) == len(null_y):
        ax.scatter(null_x, null_y, s=26, color=NEUTRAL, alpha=0.75, linewidths=0, zorder=2,
                   label=f"{len(null_x)} random edits")
    ax.scatter(cand_x, cand_y, s=62, color=BLUE, edgecolors=SURFACE, linewidths=1.6, zorder=3,
               label=f"{len(cand_x)} candidates tested")
    ax.axhline(limit, color=ORANGE, lw=1.0, zorder=1)
    ax.text(ax.get_xlim()[1], limit, f" {limit_label}", color=ORANGE, fontsize=8.2, va="bottom", ha="right")

    if hazard == "cpg_methylation":
        from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

        ax.set_yscale("log")
        ax.yaxis.set_major_locator(FixedLocator([0.01, 0.1, 1.0]))
        ax.yaxis.set_minor_locator(NullLocator())
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        ax.set_title("Methylation", loc="left", fontsize=10.5, fontweight="semibold", color=INK, pad=10)
    else:
        ax.set_title("Splicing", loc="left", fontsize=10.5, fontweight="semibold", color=INK, pad=10)
        ax.set_ylim(bottom=-0.02)
    ax.set_xlabel("How strongly the edit lowers the gene")
    ax.set_ylabel(ylabel)
    ax.legend(loc="upper left", fontsize=8.4, labelcolor=INK_2, handletextpad=0.4, borderaxespad=0.2)

    # Name the single most informative point: the closest miss on this hazard.
    if hazard == "cpg_methylation" and len(cand_y):
        over = [(y, x) for x, y in zip(cand_x, cand_y) if y > limit]
        if over:
            y, x = min(over)
            ax.annotate(f"closest miss\n{y:.3f}", (x, y), xytext=(x + 0.02, y * 0.22),
                        fontsize=8.2, color=INK_2, ha="left", va="center",
                        arrowprops=dict(arrowstyle="-", color=INK_3, lw=0.8, shrinkA=0, shrinkB=4))


def _panel_forest(ax: Any, result: OptimizationRunResult) -> None:
    hc = result.hazard_coupling or {}
    couplings = {c["hazard"]: c for c in hc.get("couplings", [])}
    cpg, spl = couplings.get("cpg_methylation"), couplings.get("splice_junction")
    if not cpg or "ci_low" not in cpg:
        ax.axis("off")
        ax.text(0.5, 0.5, "Run again to compute the confidence intervals", ha="center", va="center",
                color=INK_3, fontsize=9, transform=ax.transAxes)
        return
    rows = [
        ("Methylation\nall edits", cpg["spearman_rho"], cpg["ci_low"], cpg["ci_high"], BLUE, True),
        ("Methylation\nrandom edits only", cpg["spearman_rho_null_only"], cpg["ci_low_null_only"], cpg["ci_high_null_only"], BLUE, False),
    ]
    if spl:
        rows += [
            ("Splicing\nall edits", spl["spearman_rho"], spl["ci_low"], spl["ci_high"], NEUTRAL, True),
            ("Splicing\nrandom edits only", spl["spearman_rho_null_only"], spl["ci_low_null_only"], spl["ci_high_null_only"], NEUTRAL, False),
        ]
    _clean_axes(ax, "x")
    ax.axvline(0, color=INK_3, lw=1.0, zorder=1)
    for i, (label, rho, lo, hi, color, filled) in enumerate(rows):
        ax.plot([lo, hi], [i, i], color=color, lw=2.2, solid_capstyle="round", zorder=2)
        ax.scatter([rho], [i], s=78, color=color if filled else SURFACE, edgecolors=color, linewidths=2.0, zorder=3)
        ax.text(1.0, i, f"{rho:+.2f}", fontsize=9, color=INK, va="center", ha="right", fontweight="semibold",
                transform=ax.get_yaxis_transform())
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows], fontsize=9, color=INK_2, linespacing=1.3)
    ax.set_ylim(len(rows) - 0.4, -0.6)
    ax.set_xlim(-0.35, 1.0)
    ax.set_xlabel("Correlation with how strongly the edit lowers the gene\n(bar = 95% interval)")
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0, pad=6)


# ---- public figures ------------------------------------------------------------------------


def fig_scorecard(result: OptimizationRunResult) -> Any:
    _theme()
    n = len(result.iterations)
    fig = plt.figure(figsize=(9.6, 3.0 + 0.46 * n), dpi=_DPI)
    ax = fig.add_axes([0.03, 0.05, 0.94, 0.80])
    _panel_scorecard(ax, result)
    verdict = "No candidate passed every check." if not result.converged else "A candidate passed every check."
    _title(fig, f"{result.gene_symbol}: {verdict}",
           f"{n} single-letter edits tested against six safety checks. Orange crosses are failures.")
    _footer(fig, result)
    return fig


def fig_coupling(result: OptimizationRunResult) -> Any:
    _theme()
    fig = plt.figure(figsize=(11, 4.9), dpi=_DPI)
    gs = GridSpec(1, 2, figure=fig, left=0.075, right=0.975, top=0.78, bottom=0.17, wspace=0.22)
    _panel_coupling(fig.add_subplot(gs[0]), result, "cpg_methylation")
    _panel_coupling(fig.add_subplot(gs[1]), result, "splice_junction")
    _title(fig, "Edits that work also tend to disturb methylation",
           "Each dot is one edit. Grey dots are random edits tried before any candidate was judged.")
    _footer(fig, result)
    return fig


def fig_forest(result: OptimizationRunResult) -> Any:
    _theme()
    fig = plt.figure(figsize=(9.6, 3.5), dpi=_DPI)
    ax = fig.add_axes([0.20, 0.22, 0.77, 0.52])
    _panel_forest(ax, result)
    _title(fig, "A moderate link, not a precise one",
           "Spearman correlation between how strongly an edit lowers the gene and how much it shifts each hazard.")
    _footer(fig, result)
    return fig


def fig_loci(sweep: list[dict[str, Any]]) -> Any:
    """Per-locus correlation: all edits vs random edits only (dumbbell, one axis)."""
    _theme()
    rows = sorted(sweep, key=lambda r: r.get("meth_rho", 0.0))
    n = len(rows)
    fig = plt.figure(figsize=(9.6, 1.9 + 0.5 * n), dpi=_DPI)
    ax = fig.add_axes([0.16, 0.14, 0.78, 0.68])
    _clean_axes(ax, "x")
    ax.axvline(0, color=INK_3, lw=1.0, zorder=1)
    for i, r in enumerate(rows):
        a, b = r.get("meth_rho_null_only", 0.0), r.get("meth_rho", 0.0)
        ax.plot([a, b], [i, i], color=NEUTRAL_SOFT, lw=3.0, solid_capstyle="round", zorder=1)
        ax.scatter([b], [i], s=80, color=BLUE, edgecolors=SURFACE, linewidths=1.6, zorder=3)
        ax.scatter([a], [i], s=80, color=SURFACE, edgecolors=BLUE, linewidths=2.0, zorder=3)
    ax.set_yticks(range(n))
    ax.set_yticklabels([r["gene"] for r in rows], fontsize=10, color=INK, fontweight="semibold")
    ax.set_ylim(n - 0.4, -0.6)
    ax.set_xlim(-0.2, 0.9)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0, pad=8)
    ax.set_xlabel("Correlation between how strongly an edit lowers the gene and its methylation shift")
    ax.scatter([], [], s=80, color=BLUE, edgecolors=SURFACE, linewidths=1.6, label="all scored edits")
    ax.scatter([], [], s=80, color=SURFACE, edgecolors=BLUE, linewidths=2.0, label="random edits only")
    ax.legend(loc="lower right", fontsize=8.6, labelcolor=INK_2, borderaxespad=0.1)
    _title(fig, "The same pattern shows up at six exhaustion-related genes",
           "Positive at every gene, and consistently weaker among random edits, so treat it as suggestive.")
    _footer(fig, None, "seed 42")
    return fig


def fig_pipeline() -> Any:
    """A four-step picture of how a candidate edit is proposed, measured and judged."""
    _theme()
    fig = plt.figure(figsize=(11, 3.3), dpi=_DPI)
    ax = fig.add_axes([0.02, 0.04, 0.96, 0.68])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 30)
    ax.axis("off")
    steps = [
        ("1", "Propose", "Find places in the gene where\na change could turn it down.", "Patient group A"),
        ("2", "Measure", "Try each change in the model and\nkeep the ones that really work.", "AlphaGenome 2"),
        ("3", "Check on new patients", "Re-test on patients the search\nnever saw, plus healthy cells.", "Patient group B"),
        ("4", "Safety checks", "Splicing, methylation, neighbouring\ngenes, bystander edits, other tissues.", "Fixed limits"),
    ]
    w, gap = 21.4, 3.7
    for i, (num, head, body, tag) in enumerate(steps):
        x0 = i * (w + gap) + 0.5
        ax.add_patch(FancyBboxPatch((x0, 3), w, 24, boxstyle="round,pad=0,rounding_size=1.6",
                                    fc="#f7f7f5", ec=GRID, lw=1.0))
        ax.scatter([x0 + 3.0], [22.6], s=330, color=BLUE, linewidths=0, zorder=3)
        ax.text(x0 + 3.0, 22.55, num, color="white", fontsize=9.5, fontweight="semibold", ha="center", va="center")
        ax.text(x0 + 6.4, 22.6, head, fontsize=10.6, fontweight="semibold", color=INK, ha="left", va="center")
        ax.text(x0 + 2.0, 14.0, body, fontsize=8.5, color=INK_2, ha="left", va="center", linespacing=1.45)
        ax.text(x0 + 2.0, 6.0, tag.upper(), fontsize=7.2, color=BLUE, ha="left", va="center", fontweight="semibold")
        if i < len(steps) - 1:
            ax.add_patch(FancyArrowPatch((x0 + w + 0.5, 15), (x0 + w + gap - 0.5, 15), arrowstyle="-|>",
                                         mutation_scale=9, color=INK_3, lw=1.1))
    _title(fig, "How a candidate edit is judged",
           "The patients used to propose an edit are never the ones used to judge it.")
    return fig


def _save(fig: Any, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=_DPI)
    plt.close(fig)
    return path


def export_figure_set(
    result: OptimizationRunResult,
    out_dir: Path | str = "outputs",
    sweep_json: Path | str | None = None,
) -> list[Path]:
    """Write the README figures. The per-locus chart needs the saved sweep JSON; skipped if absent."""
    out = Path(out_dir)
    paths = [
        _save(fig_pipeline(), out / "fig_pipeline.png"),
        _save(fig_scorecard(result), out / "fig_scorecard.png"),
        _save(fig_coupling(result), out / "fig_coupling.png"),
        _save(fig_forest(result), out / "fig_forest.png"),
    ]
    sweep_path = Path(sweep_json) if sweep_json else out / "sweep_coupling_seed42.json"
    if sweep_path.exists():
        paths.append(_save(fig_loci(json.loads(sweep_path.read_text(encoding="utf-8"))), out / "fig_loci.png"))
    return paths


def export_optimization_figure(
    result: OptimizationRunResult,
    output_path: Path | str = "outputs/gene_edit_check_report.png",
) -> Path:
    """One combined summary figure: scorecard on top, coupling and its uncertainty below."""
    _theme()
    n = len(result.iterations)
    fig = plt.figure(figsize=(11.5, 4.6 + 0.42 * n + 3.9), dpi=_DPI)
    top = 0.935
    gs = GridSpec(
        2, 2, figure=fig, height_ratios=[0.42 * n + 2.2, 3.9], width_ratios=[1.0, 1.0], left=0.075, right=0.975, top=top - 0.04,
        bottom=0.10, hspace=0.26, wspace=0.55,
    )
    ax_sc = fig.add_subplot(gs[0, :])
    _panel_scorecard(ax_sc, result)
    _panel_coupling(fig.add_subplot(gs[1, 0]), result, "cpg_methylation")
    _panel_forest(fig.add_subplot(gs[1, 1]), result)
    verdict = "no candidate passed every check" if not result.converged else "a candidate passed every check"
    _title(fig, f"{result.gene_symbol}: {verdict}",
           f"{result.locus.chromosome}:{result.locus.start:,}-{result.locus.end:,}  ·  goal: "
           f"{result.desired_direction}-regulate  ·  {n} candidates, {result.efficacy_null.n_random_edits if result.efficacy_null else 0} random edits")
    _footer(fig, result)
    return _save(fig, Path(output_path))
