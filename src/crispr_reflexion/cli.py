"""Primary CLI Entrypoint (`crispr_reflexion/cli.py`).

Commands:
- `crispr-reflexion optimize`: donor-held-out coreset -> gradient-propose -> scoped-ISM ->
  two-tier reflexion, on real single-cell data (default) or synthetic demo data.
- `crispr-reflexion ablate-sketch`: AG2-latent coreset vs random subsampling.
- `crispr-reflexion validate-fidelity`: does AlphaGenome 2, conditioned on real cell sketches,
  reproduce the fold-changes observed in those cells' full populations? (LIVE backend only)
- `crispr-reflexion conditions`: which reference conditions the model actually *observed* vs
  imputed, per modality -- this is what decides where Tier 2 hazard calls come from.
- `crispr-reflexion sweep-coupling`: run the same pipeline across several checkpoint loci to
  test whether the efficacy-vs-methylation coupling is locus-specific or general.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import typer
from dotenv import load_dotenv
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from crispr_reflexion.ag_backend.live_client import LiveAlphaGenomeBackend
from crispr_reflexion.ag_backend.protocol import VerifierBackend
from crispr_reflexion.ag_backend.simulator import SurrogateSimulatorBackend
from crispr_reflexion.benchmarks.state_fidelity import run_state_fidelity
from crispr_reflexion.budget import RunBudgetTracker
from crispr_reflexion.cache import ContentAddressedCache
from crispr_reflexion.conditions import refresh_conditions_table, search_conditions
from crispr_reflexion.models import BackendMode
from crispr_reflexion.pipeline import build_all_sketches, load_cohort, run_scenario
from crispr_reflexion.scenarios.registry import (
    COUPLING_SWEEP_GENES,
    IMMUNE_ENSEMBL_PANEL,
    SCENARIOS,
    SWEEP_SCENARIO_IDS,
)
from crispr_reflexion.ui.plots import export_optimization_figure
from crispr_reflexion.ui.terminal_hud import render_optimization_report

load_dotenv()

DEFAULT_SCENARIO = "crc-pd1-exhaustion"

app = typer.Typer(
    name="crispr-reflexion",
    help="Two-tier single-cell selectivity & sequence-hazard verifier and gradient-guided CRISPR base-editing optimizer for AlphaGenome 2.",
    add_completion=False,
)
console = Console()

SIMULATOR_WARNING = (
    "SIMULATOR BACKEND: a SCRIPTED offline stand-in that hard-codes where the 'selective' site is. "
    "Its output is NOT a prediction and must not be presented as a result. Use it for plumbing tests only."
)


def _resolve_backend(
    use_simulator: bool,
    budget: RunBudgetTracker,
    cache_dir: str = ".ag2_cache",
    seed: int = 42,
) -> VerifierBackend:
    api_key = os.environ.get("ALPHAGENOME_API_KEY", "").strip()
    if not use_simulator and api_key:
        cache = ContentAddressedCache(cache_dir=cache_dir)
        return LiveAlphaGenomeBackend(api_key=api_key, budget=budget, cache=cache)
    console.print(Panel(SIMULATOR_WARNING, border_style="bold red", title="WARNING"))
    return SurrogateSimulatorBackend(budget=budget, seed=seed)


def _get_spec(scenario: str):  # type: ignore[no-untyped-def]
    if scenario not in SCENARIOS:
        console.print(f"[bold red]Unknown scenario '{scenario}'. Available: {list(SCENARIOS)}[/bold red]")
        raise typer.Exit(code=1)
    return SCENARIOS[scenario]


@app.command("optimize")
def cmd_optimize(
    scenario: str = typer.Option(DEFAULT_SCENARIO, help="Scenario ID from registry."),
    sketch_size: int = typer.Option(128, min=4, max=256, help="Cells per disjoint sketch (<=256)."),
    pool_size: int = typer.Option(600, min=40, help="Candidate cells per state before sketching."),
    max_iterations: int = typer.Option(6, min=1, max=12, help="Max reflexion iterations."),
    null_n: int = typer.Option(12, min=0, help="Random PAM-valid edits for the efficacy null (0 disables)."),
    seed: int = typer.Option(42, help="Random seed."),
    simulator: bool = typer.Option(False, "--simulator", help="Force the SCRIPTED offline simulator (plumbing only)."),
    output_dir: str = typer.Option("outputs", help="Directory for JSON report and figure."),
) -> None:
    """Run the full pipeline on a scenario and write a JSON report + figure."""
    spec = _get_spec(scenario)
    budget = RunBudgetTracker()
    backend = _resolve_backend(use_simulator=simulator, budget=budget, seed=seed)
    result = run_scenario(
        spec, backend, budget, sketch_size=sketch_size, max_iterations=max_iterations,
        seed=seed, pool_size=pool_size, null_n=null_n,
    )
    render_optimization_report(result, console=console)

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / f"report_{spec.scenario_id}.json"
    fig_path = out_dir / f"figure_{spec.scenario_id}.png"
    json_path.write_text(result.model_dump_json(indent=2), encoding="utf-8")
    export_optimization_figure(result, output_path=fig_path)
    console.print(f"[bold green]Saved JSON report -> {json_path} | Figure -> {fig_path}[/bold green]")


@app.command("ablate-sketch")
def cmd_ablate_sketch(
    scenario: str = typer.Option(DEFAULT_SCENARIO),
    sketch_size: int = typer.Option(128, min=4, max=256),
    pool_size: int = typer.Option(600, min=40),
    seed: int = typer.Option(42),
    simulator: bool = typer.Option(False, "--simulator"),
) -> None:
    """Compare facility-location coresets in AG2 latent space vs uniform random subsampling."""
    spec = _get_spec(scenario)
    budget = RunBudgetTracker()
    backend = _resolve_backend(use_simulator=simulator, budget=budget, seed=seed)
    cohort = load_cohort(spec, pool_size=pool_size, seed=seed)
    _, _, metrics = build_all_sketches(backend, spec, cohort, sketch_size, seed)

    tbl = Table(title=f"AG2 latent coreset ablation [{backend.backend_mode.value}] seed={seed}", show_lines=True)
    for col in ("Cell state", "Pool -> Sketch A/B", "Split", "Donors A/B (overlap)", "MMD² coreset", "MMD² random", "Reduction", "A↔B cosine"):
        tbl.add_column(col)
    for m in metrics:
        tbl.add_row(
            m.cell_state_name, f"{m.total_source_cells} -> 2x{m.sketch_size}", f"by {m.split_by}",
            f"{m.n_donors_sketch_a}/{m.n_donors_sketch_b} ({m.donors_overlap_a_b})",
            f"{m.ag2_latent_mmd_stratified:.5f}", f"{m.ag2_latent_mmd_unstratified_random:.5f}",
            f"{m.mmd_reduction_pct:.1f}%", f"{m.sketch_a_B_cosine_consistency:.4f}",
        )
    console.print(tbl)


@app.command("validate-fidelity")
def cmd_validate_fidelity(
    scenario: str = typer.Option(DEFAULT_SCENARIO),
    sketch_size: int = typer.Option(128, min=4, max=256),
    pool_size: int = typer.Option(600, min=40),
    seed: int = typer.Option(42),
) -> None:
    """Predicted vs observed tumor-vs-normal fold-change across a fixed 24-gene immune panel."""
    spec = _get_spec(scenario)
    if spec.data_source != "real":
        console.print("[bold red]validate-fidelity needs a real-data scenario.[/bold red]")
        raise typer.Exit(code=1)
    budget = RunBudgetTracker(max_gradient_calls=2 * len(IMMUNE_ENSEMBL_PANEL) + 4)
    backend = _resolve_backend(use_simulator=False, budget=budget, seed=seed)
    if backend.backend_mode != BackendMode.LIVE_ALPHAGENOME_V2:
        console.print("[bold red]Fidelity is only meaningful on the LIVE backend (set ALPHAGENOME_API_KEY).[/bold red]")
        raise typer.Exit(code=1)

    cohort = load_cohort(spec, pool_size=pool_size, seed=seed)
    assert cohort.observed_cp10k is not None
    target, healthy, _ = build_all_sketches(backend, spec, cohort, sketch_size, seed)
    h_name = spec.healthy_state_names[0]
    res = run_state_fidelity(
        backend, target.sketch_a, healthy[h_name].sketch_a, cohort.observed_cp10k,
        spec.target_state_name, h_name, list(IMMUNE_ENSEMBL_PANEL),
    )
    tbl = Table(
        title=(
            f"Cell-state fidelity [{backend.backend_mode.value}] | Spearman rho={res.spearman_rho:+.3f} "
            f"(p={res.spearman_p:.4f}) | sign concordance={res.sign_concordance_pct:.0f}% "
            f"of {res.n_concordance_genes} genes with |obs lnFC|>=0.5 | slope pred/obs={res.slope_pred_vs_obs:.2f}"
        ),
        show_lines=False,
    )
    for col in ("Gene", "Observed CP10k (tumor / normal)", "Observed ln FC", "Predicted ln FC"):
        tbl.add_column(col, justify="right")
    for r in sorted(res.records, key=lambda r: -abs(float(r["observed_ln_fc"]))):  # type: ignore[arg-type]
        tbl.add_row(
            str(r["gene"]), f"{float(r['obs_cp10k_target']):.2f} / {float(r['obs_cp10k_healthy']):.2f}",  # type: ignore[arg-type]
            f"{float(r['observed_ln_fc']):+.2f}", f"{float(r['predicted_ln_fc']):+.2f}",  # type: ignore[arg-type]
        )
    console.print(tbl)


@app.command("sweep-coupling")
def cmd_sweep_coupling(
    genes: str = typer.Option(
        "", help="Comma-separated gene symbols (default: all six checkpoint/exhaustion loci)."
    ),
    sketch_size: int = typer.Option(128, min=4, max=256),
    pool_size: int = typer.Option(600, min=40),
    max_iterations: int = typer.Option(5, min=1, max=12),
    null_n: int = typer.Option(10, min=0),
    prescreen_windows: int = typer.Option(8, min=0),
    seed: int = typer.Option(42),
    simulator: bool = typer.Option(False, "--simulator"),
    output_dir: str = typer.Option("outputs", help="Directory for the sweep JSON."),
) -> None:
    """Is the efficacy-vs-methylation coupling specific to PDCD1, or general to this cell state?

    Runs the identical pipeline -- same cells, sketches, thresholds and reference conditions --
    across several checkpoint/exhaustion genes and compares the coupling at each locus. A result
    that only appears at PDCD1 says something about that locus; one that appears everywhere says
    something about base editing in immune regulatory regions generally. Either is informative,
    which is why this is worth running rather than assuming.
    """
    symbols = (
        [g.strip().upper() for g in genes.split(",") if g.strip()]
        if genes
        else [s for _, s in COUPLING_SWEEP_GENES]
    )
    rows: list[dict[str, object]] = []

    for sym in symbols:
        scenario_id = SWEEP_SCENARIO_IDS.get(sym)
        if scenario_id is None:
            console.print(f"[yellow]Skipping unknown gene '{sym}'.[/yellow]")
            continue
        spec = SCENARIOS[scenario_id]
        console.rule(f"[bold cyan]{sym}[/bold cyan] ({spec.gene_id})")
        budget = RunBudgetTracker()
        backend = _resolve_backend(use_simulator=simulator, budget=budget, seed=seed)
        try:
            result = run_scenario(
                spec, backend, budget, sketch_size=sketch_size, max_iterations=max_iterations,
                seed=seed, pool_size=pool_size, null_n=null_n,
                ism_prescreen_windows=prescreen_windows, off_tissue_panel=False,
            )
        except Exception as exc:  # one failing locus must not lose the whole sweep
            console.print(f"[bold red]{sym} failed: {type(exc).__name__}: {exc}[/bold red]")
            rows.append({"gene": sym, "error": f"{type(exc).__name__}: {exc}"})
            continue

        hc = result.hazard_coupling or {}
        meth = next((c for c in hc.get("couplings", []) if c["hazard"] == "cpg_methylation"), None)
        splice = next((c for c in hc.get("couplings", []) if c["hazard"] == "splice_junction"), None)
        best = min(
            (abs(s.report.tier1_selectivity.target_delta_sketch_b) for s in result.iterations),
            default=0.0,
        )
        strongest = max(
            (abs(s.report.tier1_selectivity.target_delta_sketch_b) for s in result.iterations),
            default=0.0,
        )
        row = {
            "gene": sym,
            "scenario_id": scenario_id,
            "converged": result.converged,
            "n_iterations": len(result.iterations),
            "verdict": hc.get("verdict"),
            "coupling_class": hc.get("coupling_class"),
            "meth_minus_splice_rho": hc.get("methylation_minus_splice_rho"),
            "binding_constraint": hc.get("binding_constraint"),
            "n_failing_binding": hc.get("n_candidates_failing_binding_constraint"),
            "meth_rho": meth["spearman_rho"] if meth else None,
            "meth_p": meth["spearman_p"] if meth else None,
            "meth_rho_null_only": meth["spearman_rho_null_only"] if meth else None,
            "meth_p_null_only": meth["spearman_p_null_only"] if meth else None,
            "splice_rho": splice["spearman_rho"] if splice else None,
            "splice_p": splice["spearman_p"] if splice else None,
            "strongest_abs_lfc": round(strongest, 6),
            "weakest_abs_lfc": round(best, 6),
            "gradient_sign_agreement": (
                result.ism_prescreen.gradient_sign_agreement if result.ism_prescreen else None
            ),
            "api_calls": result.provenance.total_api_calls,
        }
        rows.append(row)
        console.print(
            f"  {sym}: verdict={row['verdict']} | methylation rho={row['meth_rho']} "
            f"(p={row['meth_p']}) | splice rho={row['splice_rho']} | "
            f"blocked by {row['binding_constraint']} in {row['n_failing_binding']}/{row['n_iterations']}"
        )

    tbl = Table(
        title="Efficacy-vs-hazard coupling across checkpoint loci (identical cells, sketches and thresholds)",
        show_lines=False,
    )
    for col in (
        "Gene", "Verdict", "Binding check", "Blocked", "CpG rho", "CpG p",
        "CpG rho (random only)", "Splice rho", "Gap", "Coupling class",
    ):
        tbl.add_column(col, justify="right")
    for r in rows:
        if r.get("error"):
            tbl.add_row(str(r["gene"]), "[red]error[/red]", str(r["error"])[:40], "", "", "", "", "", "", "")
            continue
        cls = str(r.get("coupling_class") or "-")
        cls_style = {
            "methylation_specific": "bold magenta",
            "generic_hazard_coupling": "yellow",
            "uncoupled": "dim",
        }.get(cls, "dim")
        rho = r["meth_rho"]
        rho_s = f"[bold]{rho:+.2f}[/bold]" if isinstance(rho, float) and rho > 0.5 else f"{rho:+.2f}" if isinstance(rho, float) else "-"
        tbl.add_row(
            str(r["gene"]),
            str(r["verdict"] or "-").replace("_", " "),
            str(r["binding_constraint"] or "-"),
            f"{r['n_failing_binding']}/{r['n_iterations']}",
            rho_s,
            f"{r['meth_p']:.2g}" if isinstance(r["meth_p"], float) else "-",
            f"{r['meth_rho_null_only']:+.2f}" if isinstance(r["meth_rho_null_only"], float) else "-",
            f"{r['splice_rho']:+.2f}" if isinstance(r["splice_rho"], float) else "-",
            f"{r['meth_minus_splice_rho']:+.2f}" if isinstance(r.get("meth_minus_splice_rho"), float) else "-",
            f"[{cls_style}]{cls.replace('_', ' ')}[/]",
        )
    console.print(tbl)

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"sweep_coupling_seed{seed}.json"
    path.write_text(json.dumps(rows, indent=2), encoding="utf-8")
    console.print(f"[bold green]Saved sweep -> {path}[/bold green]")


@app.command("conditions")
def cmd_conditions(
    query: str = typer.Argument("CD8-positive", help="Substring/regex matched against biosample names."),
    modality: str = typer.Option(
        "methylation", help="Sort by evidence for this modality (methylation, splice_junctions, rna_seq, ...)."
    ),
    limit: int = typer.Option(25, min=1, help="Rows to show."),
    refresh: bool = typer.Option(False, "--refresh", help="Re-download the table from the API."),
) -> None:
    """Explore AlphaGenome 2 reference conditions by how well the model actually knows them.

    `supported_conditions()` reports, per condition and modality, whether the model has
    `observed` training tracks, a validated imputation (`imputed_supported`), or an unvalidated
    one (`imputed_exploratory`). Only 185 of 2,799 human conditions have observed methylation,
    so this is what decides which reference a Tier 2 hazard call should be made under.
    """
    if refresh:
        api_key = os.environ.get("ALPHAGENOME_API_KEY", "").strip()
        if not api_key:
            console.print("[bold red]--refresh needs ALPHAGENOME_API_KEY.[/bold red]")
            raise typer.Exit(code=1)
        path = refresh_conditions_table(api_key)
        console.print(f"[green]Refreshed {path}[/green]")

    hits = search_conditions(query, modality=modality)
    if hits.empty:
        console.print(f"[yellow]No conditions matched '{query}'.[/yellow]")
        raise typer.Exit(code=0)

    style = {"observed": "bold green", "imputed_supported": "yellow", "imputed_exploratory": "dim red"}
    tbl = Table(
        title=f"AlphaGenome 2 reference conditions matching '{query}' ({len(hits)} hits, sorted by {modality} evidence)",
        show_lines=False,
    )
    tbl.add_column("CURIE", style="cyan")
    tbl.add_column("Biosample")
    shown = ("methylation", "splice_junctions", "rna_seq", "chip_tf", "sc_atac")
    for m in shown:
        tbl.add_column(m.replace("_", "\n"), justify="center")

    for _, row in hits.head(limit).iterrows():
        cells = []
        for m in shown:
            tier = str(row.get(m, "unknown"))
            cells.append(f"[{style.get(tier, 'dim')}]{tier.replace('imputed_', 'imp_')}[/]")
        tbl.add_row(str(row["condition"]), str(row["biosample_name"]), *cells)
    console.print(tbl)
    console.print(
        "[dim]observed = real training tracks | imp_supported = validated imputation | "
        "imp_exploratory = unvalidated imputation (weakest evidence)[/dim]"
    )


if __name__ == "__main__":
    app()
