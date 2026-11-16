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
from dotenv import find_dotenv, load_dotenv
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from crispr_reflexion import __version__, paths
from crispr_reflexion.ag_backend.live_client import LiveAlphaGenomeBackend
from crispr_reflexion.ag_backend.protocol import VerifierBackend
from crispr_reflexion.ag_backend.simulator import SurrogateSimulatorBackend
from crispr_reflexion.benchmarks.state_fidelity import run_state_fidelity
from crispr_reflexion.benchmarks.universality import run_multiseed_and_context_invariance
from crispr_reflexion.budget import RunBudgetTracker
from crispr_reflexion.cache import ContentAddressedCache
from crispr_reflexion.conditions import refresh_conditions_table, resolve_reference_pair, search_conditions
from crispr_reflexion.datasets.generic_h5ad import GenericH5adSpec
from crispr_reflexion.models import BackendMode
from crispr_reflexion.pipeline import build_all_sketches, load_cohort, resolve_scenario_references, run_scenario
from crispr_reflexion.scenarios.registry import (
    COUPLING_SWEEP_GENES,
    IMMUNE_ENSEMBL_PANEL,
    SCENARIOS,
    SWEEP_SCENARIO_IDS,
    UNIVERSALITY_CONTEXT_LABELS,
    ScenarioSpec,
)
from crispr_reflexion.ui.plots import export_optimization_figure
from crispr_reflexion.ui.terminal_hud import render_optimization_report

load_dotenv(find_dotenv(usecwd=True))

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
    cache_dir: str | None = None,
    seed: int = 42,
    require_live: bool = False,
) -> VerifierBackend:
    """Pick the live AlphaGenome 2 backend, or the scripted simulator when explicitly allowed.

    `require_live` is set for real-data runs. Without an API key those used to fall back to the
    simulator and print an ordinary-looking report, which a newcomer could mistake for a result.
    Now they stop and say what to do. The simulator is still available with `--simulator`.
    """
    api_key = os.environ.get("ALPHAGENOME_API_KEY", "").strip()
    if not use_simulator and api_key:
        cache = ContentAddressedCache(cache_dir=cache_dir or paths.api_cache_dir())
        return LiveAlphaGenomeBackend(api_key=api_key, budget=budget, cache=cache)
    if not use_simulator and require_live:
        console.print(
            Panel(
                "This run uses real data, so it needs the live AlphaGenome 2 model and no API key was found.\n\n"
                "1. Get a free non-commercial key: https://deepmind.google.com/science/alphagenome\n"
                "2. Set it:  export ALPHAGENOME_API_KEY=...   (PowerShell: $env:ALPHAGENOME_API_KEY='...')\n"
                "   or put ALPHAGENOME_API_KEY=... in a .env file in this folder.\n"
                "3. Check your setup any time with:  crispr-reflexion doctor\n\n"
                "To exercise the plumbing without a key, add --simulator. Its output is scripted and is not a result.",
                title="ALPHAGENOME_API_KEY not set",
                border_style="yellow",
            )
        )
        raise typer.Exit(code=2)
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
    backend = _resolve_backend(use_simulator=simulator, budget=budget, seed=seed, require_live=spec.data_source != "synthetic")
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
    backend = _resolve_backend(use_simulator=simulator, budget=budget, seed=seed, require_live=spec.data_source != "synthetic")
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
    backend = _resolve_backend(use_simulator=False, budget=budget, seed=seed, require_live=True)
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
        backend = _resolve_backend(use_simulator=simulator, budget=budget, seed=seed, require_live=spec.data_source != "synthetic")
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


@app.command("doctor")
def cmd_doctor() -> None:
    """Check your setup: version, API key, where data and cache live, and what is already downloaded."""
    from crispr_reflexion.conditions import load_conditions_table
    from crispr_reflexion.datasets import crc_pelka

    key = os.environ.get("ALPHAGENOME_API_KEY", "").strip()
    try:
        n_cond = len(load_conditions_table())
        cond = f"{n_cond:,} reference conditions shipped with the package"
        cond_ok = True
    except Exception as exc:  # pragma: no cover - means a broken install
        cond, cond_ok = f"missing ({type(exc).__name__})", False
    raw = crc_pelka.RAW_PATH
    cache = paths.api_cache_dir()
    cache_mb = sum(f.stat().st_size for f in cache.rglob("*") if f.is_file()) / 1e6 if cache.exists() else 0.0

    ok, no, warn = "[green]ok[/green]", "[red]missing[/red]", "[yellow]not yet[/yellow]"
    tbl = Table(title=f"crispr-reflexion {__version__}: setup check", show_header=False, show_lines=False)
    tbl.add_column(style="bold")
    tbl.add_column()
    tbl.add_column()
    tbl.add_row("API key", ok if key else no,
                "ALPHAGENOME_API_KEY is set" if key else "needed for real runs; free key: https://deepmind.google.com/science/alphagenome")
    tbl.add_row("Reference conditions", ok if cond_ok else no, cond)
    tbl.add_row("Atlas download", ok if raw.exists() else warn,
                str(raw) if raw.exists() else f"downloads on first real run (about 396 MB) to {raw}")
    tbl.add_row("API response cache", ok if cache.exists() else warn, f"{cache} ({cache_mb:.0f} MB)")
    tbl.add_row("Data / cache home", "", f"{paths.home()}  ({'source checkout' if paths.is_source_checkout() else 'per-user cache'}; override with CRISPR_REFLEXION_HOME)")
    console.print(tbl)

    if not key:
        console.print("\n[bold]Without a key you can still:[/bold] run `crispr-reflexion conditions`, run the test suite, "
                      "and try the plumbing with `--simulator` (scripted output, never a result).")
    else:
        console.print("\n[bold green]Ready.[/bold green] Try: crispr-reflexion optimize")


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


@app.command("run-h5ad")
def cmd_run_h5ad(
    h5ad: str = typer.Option(..., "--h5ad", help="Local path or HTTP(S) URL to any .h5ad file."),
    gene_id: str = typer.Option(..., "--gene-id", help="Target Ensembl gene ID (e.g. ENSG00000188389)."),
    gene_symbol: str = typer.Option(..., "--gene-symbol", help="Target gene symbol (e.g. PDCD1)."),
    target_query: str = typer.Option(
        ..., "--target-query", help="Pandas .obs.query(...) expression (or expr_<GENE> predicate) for target cells."
    ),
    healthy_query: str = typer.Option(
        ..., "--healthy-query", help="Pandas .obs.query(...) expression (or expr_<GENE> predicate) for control cells."
    ),
    target_name: str = typer.Option("Target_State", "--target-name"),
    healthy_name: str = typer.Option("Healthy_State", "--healthy-name"),
    dataset_name: str = typer.Option("Custom_SingleCell_Cohort", "--dataset-name"),
    donor_col: str | None = typer.Option(None, "--donor-col", help="Donor/batch column in .obs (auto-detected if omitted)."),
    ensembl_col: str | None = typer.Option(None, "--ensembl-col", help="Ensembl column in .var (auto-detected if omitted)."),
    cell_type: str = typer.Option("T cell", "--cell-type", help="Cell-type query for evidence-ranked Tier 2 reference."),
    tissue: str | None = typer.Option(None, "--tissue", help="Optional tissue query for anatomical cross-check."),
    desired_direction: str = typer.Option("down", "--desired-direction", help="'down' or 'up'."),
    sketch_size: int = typer.Option(128, min=4, max=256),
    pool_size: int = typer.Option(600, min=40),
    max_iterations: int = typer.Option(6, min=1, max=12),
    null_n: int = typer.Option(12, min=0),
    seed: int = typer.Option(42),
    simulator: bool = typer.Option(False, "--simulator"),
    output_dir: str = typer.Option("outputs"),
) -> None:
    """Dataset-agnostic runner: execute the full Two-Tier CRISPR-Reflexion pipeline on any `.h5ad` file."""
    resolved_ref = resolve_reference_pair(cell_type, tissue)
    console.print(
        f"[cyan]Auto-selected Tier 2 primary reference:[/cyan] {resolved_ref.primary_biosample_name} "
        f"({list(resolved_ref.primary_curies)}, limiting tier: [bold]{resolved_ref.primary_limiting_tier}[/bold])"
    )
    if resolved_ref.tissue_matched_curies:
        console.print(
            f"[cyan]Auto-selected tissue-matched cross-check:[/cyan] {resolved_ref.tissue_matched_biosample_name} "
            f"({list(resolved_ref.tissue_matched_curies)}, limiting tier: {resolved_ref.tissue_matched_limiting_tier})"
        )

    safe_sid = f"custom-{gene_symbol.lower()}-h5ad"
    spec = ScenarioSpec(
        scenario_id=safe_sid,
        title=f"{gene_symbol} ({desired_direction}-regulate) on {dataset_name}",
        gene_id=gene_id,
        gene_symbol=gene_symbol,
        target_state_name=target_name,
        healthy_state_names=(healthy_name,),
        reference_curies=resolved_ref.primary_curies,
        tissue_matched_curies=resolved_ref.tissue_matched_curies,
        data_source="generic",
        desired_direction=desired_direction,
        cell_type_query=cell_type,
        tissue_query=tissue,
        generic_h5ad_spec=GenericH5adSpec(
            h5ad_path_or_url=h5ad,
            dataset_name=dataset_name,
            target_state_name=target_name,
            healthy_state_name=healthy_name,
            target_obs_query=target_query,
            healthy_obs_query=healthy_query,
            donor_col=donor_col,
            ensembl_var_col=ensembl_col,
        ),
    )
    budget = RunBudgetTracker()
    backend = _resolve_backend(use_simulator=simulator, budget=budget, seed=seed, require_live=True)
    result = run_scenario(
        spec,
        backend,
        budget,
        sketch_size=sketch_size,
        max_iterations=max_iterations,
        seed=seed,
        pool_size=pool_size,
        null_n=null_n,
    )
    render_optimization_report(result, console=console)

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / f"report_{safe_sid}.json"
    fig_path = out_dir / f"figure_{safe_sid}.png"
    json_path.write_text(result.model_dump_json(indent=2), encoding="utf-8")
    export_optimization_figure(result, output_path=fig_path)
    console.print(f"[bold green]Saved JSON report -> {json_path} | Figure -> {fig_path}[/bold green]")


@app.command("sweep-universality")
def cmd_sweep_universality(
    seeds: str = typer.Option("42,7,99", help="Comma-separated random seeds for multi-seed stability."),
    sketch_size: int = typer.Option(128, min=4, max=256),
    pool_size: int = typer.Option(600, min=40),
    max_iterations: int = typer.Option(5, min=1, max=12),
    null_n: int = typer.Option(40, min=0),
    prescreen_windows: int = typer.Option(8, min=0),
    simulator: bool = typer.Option(False, "--simulator"),
    output_dir: str = typer.Option("outputs"),
) -> None:
    """Prove invariance across seeds, tumor subtypes, T-cell lineages, an independent breast cancer cohort, and GRCh38 CpG architecture."""
    seed_tuple = tuple(int(s.strip()) for s in seeds.split(",") if s.strip())
    sweep_triples = tuple((gid, sym, SWEEP_SCENARIO_IDS[sym]) for gid, sym in COUPLING_SWEEP_GENES)

    report = run_multiseed_and_context_invariance(
        scenario_specs=SCENARIOS,
        backend_factory=lambda budget, seed: _resolve_backend(
            use_simulator=simulator, budget=budget, seed=seed, require_live=True
        ),
        run_scenario_fn=run_scenario,
        load_cohort_fn=load_cohort,
        build_sketches_fn=build_all_sketches,
        primary_scenario_id=DEFAULT_SCENARIO,
        context_labels=UNIVERSALITY_CONTEXT_LABELS,
        sweep_genes=sweep_triples,
        seeds=seed_tuple,
        sketch_size=sketch_size,
        pool_size=pool_size,
        max_iterations=max_iterations,
        null_n=null_n,
        prescreen_windows=prescreen_windows,
    )

    t_seed = Table(
        title=f"Part A: Multi-Seed Stability on {DEFAULT_SCENARIO} (mean CpG rho = {report.seed_meth_rho_mean} +/- {report.seed_meth_rho_sd})",
        show_lines=False,
    )
    for col in ("Seed", "CpG rho", "CpG p", "Random-null rho", "Splice rho", "Gap", "Grad Sign Agree", "Class", "Verdict"):
        t_seed.add_column(col, justify="right")
    for sr in report.seed_summaries:
        t_seed.add_row(
            str(sr.seed),
            f"{sr.meth_rho:+.3f}" if sr.meth_rho is not None else "-",
            f"{sr.meth_p:.2g}" if sr.meth_p is not None else "-",
            f"{sr.meth_rho_null_only:+.3f}" if sr.meth_rho_null_only is not None else "-",
            f"{sr.splice_rho:+.3f}" if sr.splice_rho is not None else "-",
            f"{sr.meth_minus_splice_rho:+.3f}" if sr.meth_minus_splice_rho is not None else "-",
            f"{sr.gradient_sign_agreement:.0%}" if sr.gradient_sign_agreement is not None else "-",
            sr.coupling_class.replace("_", " "),
            sr.verdict.replace("_", " "),
        )
    console.print(t_seed)

    t_ctx = Table(
        title="Part B: Cross-Cancer, Cross-Subtype & Cross-Lineage Invariance (PDCD1)",
        show_lines=False,
    )
    for col in ("Context", "Dataset", "Cells", "Tier 2 Ref [Tier]", "CpG rho", "CpG p", "Null rho", "Splice rho", "Gap", "Class"):
        t_ctx.add_column(col, justify="right")
    for cr in report.context_summaries:
        t_ctx.add_row(
            cr.context_label,
            cr.dataset_name[:34],
            f"{cr.target_cells_full:,}",
            f"{cr.primary_reference[:24]} [{cr.primary_evidence_tier}]",
            f"{cr.meth_rho:+.3f}" if cr.meth_rho is not None else "-",
            f"{cr.meth_p:.2g}" if cr.meth_p is not None else "-",
            f"{cr.meth_rho_null_only:+.3f}" if cr.meth_rho_null_only is not None else "-",
            f"{cr.splice_rho:+.3f}" if cr.splice_rho is not None else "-",
            f"{cr.meth_minus_splice_rho:+.3f}" if cr.meth_minus_splice_rho is not None else "-",
            cr.coupling_class.replace("_", " "),
        )
    console.print(t_ctx)

    t_arch = Table(
        title="Part C: Physical GRCh38 Sequence CpG Architecture (+/-2.5 kb cis-window) vs AlphaGenome Coupling",
        show_lines=False,
    )
    for col in ("Gene", "GC%", "CpG count", "CpG/100bp", "CpG Obs/Exp", "Local +-50bp CpG vs Meth rho", "Efficacy vs CpG rho", "Class"):
        t_arch.add_column(col, justify="right")
    for ar in report.locus_cpg_architecture:
        t_arch.add_row(
            str(ar["gene_symbol"]),
            f"{100.0 * float(ar['gc_fraction']):.1f}%",
            str(ar["cpg_dinucleotide_count"]),
            f"{float(ar['cpg_per_100bp']):.2f}",
            f"{float(ar['cpg_observed_over_expected']):.3f}",
            f"{ar['local_cpg_vs_methylation_rho']:+.2f}" if ar.get("local_cpg_vs_methylation_rho") is not None else "-",
            f"{ar['efficacy_vs_cpg_rho']:+.2f}" if ar.get("efficacy_vs_cpg_rho") is not None else "-",
            str(ar.get("coupling_class") or "-").replace("_", " "),
        )
    console.print(t_arch)

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "universality_report.json"
    out_path.write_text(report.model_dump_json(indent=2), encoding="utf-8")
    console.print(f"[bold green]Saved universality report -> {out_path}[/bold green]")


if __name__ == "__main__":
    app()
