"""Rich Terminal Telemetry & Verification HUD (`ui/terminal_hud.py`).

Renders structured tables with explicit provenance stamps (`LIVE_ALPHAGENOME_V2` vs
`SURROGATE_SIMULATOR`), Tier 1 vs Tier 2 separation, and Sketch A -> Held-Out Sketch B
transfer metrics.
"""

from __future__ import annotations

from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

from gene_edit_check import __version__
from gene_edit_check.models import BackendMode, OptimizationRunResult


def render_optimization_report(
    result: OptimizationRunResult,
    console: Console | None = None,
) -> None:
    con = console or Console()
    prov = result.provenance
    is_live = prov.backend_mode == BackendMode.LIVE_ALPHAGENOME_V2
    badge_color = "bold #4a90e2" if is_live else "bold #eb6834"

    header_text = Text()
    header_text.append(f"GENE-EDIT-CHECK v{__version__} ", style="bold #4a90e2")
    header_text.append(f"[{prov.backend_mode.value}] ", style=badge_color)
    header_text.append(
        f"| Model: {prov.model_version} | Seed: {prov.random_seed} "
        f"| API Calls: {prov.total_api_calls} (Cache Hits: {prov.cache_hits}) "
        f"| Scoped ISM Variants: {prov.total_ism_variants_scored}\n",
        style="dim",
    )
    header_text.append(prov.research_use_only_notice, style="italic dim")
    con.print(Panel(header_text, title=f"Scenario: {result.scenario_id} ({result.gene_symbol})", border_style="dim"))

    dp = result.data_provenance
    if dp:
        lines = [f"Data: {dp.get('data_source', 'unknown')}"]
        if dp.get("citation"):
            lines.append(str(dp["citation"]))
        for state, info in (dp.get("states") or {}).items():  # type: ignore[union-attr]
            lines.append(
                f"  {state}: {info['n_cells_full_population']:,} cells / {info['n_donors_full_population']} patients "
                f"(pool {info['n_cells_pool']} cells, {info['n_donors_pool']} patients)"
            )
        con.print(Panel(chr(10).join(lines), title="Dataset", border_style="blue"))
    ps = result.ism_prescreen
    if ps:
        con.print(
            f"[bold]ISM pre-screen:[/bold] {ps.n_windows} windows x {ps.window_width_bp} bp -> "
            f"{ps.n_variants_measured} substitutions measured, "
            f"{ps.n_candidates_with_measurement} candidates covered "
            f"([#4a90e2]{ps.n_correct_direction} correct direction[/#4a90e2], "
            f"[#eb6834]{ps.n_wrong_direction_filtered} filtered out[/#eb6834]); "
            f"best measured LFC {ps.best_measured_lfc:+.5f}. "
            f"Gradient sign agreed with measured ISM sign {ps.gradient_sign_agreement:.0%} of the time "
            f"[dim](the heuristic this pre-screen replaces)[/dim]"
        )

    en = result.efficacy_null
    if en:
        con.print(
            f"[bold]Efficacy null:[/bold] {en.n_random_edits} random PAM-valid edits, |LFC| max={en.abs_lfc_max:.5f} "
            f"mean={en.abs_lfc_mean:.5f} -> cutoff {en.cutoff_used:.5f} "
            f"(a pass means beating all of them, empirical p <= {en.empirical_p_floor}); "
            f"goal: {result.desired_direction}-regulate {result.gene_symbol}"
        )

    # Table 1: Batch-Stratified AG2-Latent Coreset Ablation
    sketch_tbl = Table(
        title="Phase 1: Donor-Stratified Coreset Selection in AlphaGenome 2 `encode_cells()` Latent Space (256-D)",
        show_lines=True,
    )
    sketch_tbl.add_column("Cell State", style="bold")
    sketch_tbl.add_column("Pool -> Sketch (A/B)", justify="right")
    sketch_tbl.add_column("A/B split", justify="center")
    sketch_tbl.add_column("Donors A/B (overlap)", justify="right")
    sketch_tbl.add_column("MMD^2 (Stratified Coreset)", justify="right", style="#4a90e2")
    sketch_tbl.add_column("MMD^2 (Random Subsample)", justify="right", style="#eb6834")
    sketch_tbl.add_column("MMD Reduction", justify="right", style="bold #4a90e2")
    sketch_tbl.add_column("Sketch A/B Cosine", justify="right")

    for sm in result.sketch_metrics:
        sketch_tbl.add_row(
            sm.cell_state_name,
            f"{sm.total_source_cells} -> 2x{sm.sketch_size}",
            f"by {sm.split_by}",
            f"{sm.n_donors_sketch_a}/{sm.n_donors_sketch_b} ({sm.donors_overlap_a_b})",
            f"{sm.ag2_latent_mmd_stratified:.6f}",
            f"{sm.ag2_latent_mmd_unstratified_random:.6f}",
            f"-{sm.mmd_reduction_pct:.1f}%",
            f"{sm.sketch_a_B_cosine_consistency:.4f}",
        )
    con.print(sketch_tbl)

    # Table 2: Two-Tier Search Trace
    trace_tbl = Table(
        title="Phase 2: Gradients-Propose -> Scoped-ISM -> Two-Tier Verification Loop",
        show_lines=True,
    )
    trace_tbl.add_column("Iter", justify="center", style="bold")
    trace_tbl.add_column("Candidate (Editor / Pos / Strand)")
    trace_tbl.add_column("Protospacer [20nt] + PAM")
    trace_tbl.add_column("Tier 1: scRNA LFC\n(Sketch A -> Held-Out B)", justify="right")
    trace_tbl.add_column("Tier 1: Healthy\nControl LFC (B)", justify="right")
    trace_tbl.add_column("Tier 1: Worst\nNeighbour Gene", justify="right")
    trace_tbl.add_column("Tier 2: Hazards\n(Splice / CpG / Bystanders)", justify="right")
    trace_tbl.add_column("Status & failure analysis Update")

    for step in result.iterations:
        rep = step.report
        cand = rep.candidate
        t1 = rep.tier1_selectivity
        t2 = rep.tier2_hazards
        h_val = next(iter(t1.healthy_control_deltas.values()), 0.0)
        n_bys = len(cand.bystander_edits_in_window)

        status_badge = "[bold #4a90e2]PASS[/bold #4a90e2]" if rep.all_checks_passed else "[bold #eb6834]FAIL[/bold #eb6834]"
        trace_tbl.add_row(
            str(step.iteration),
            f"{cand.editor.value}\n{cand.target_edit.chromosome}:{cand.target_edit.position} ({cand.strand})\n"
            f"{cand.target_edit.reference}->{cand.target_edit.alternate} (p{cand.target_edit.protospacer_index})",
            f"{cand.protospacer_20nt} [{cand.pam_sequence}]\n"
            f"[dim]{cand.disrupted_motif or 'Regulatory_Peak'}[/dim]",
            f"A: {t1.target_delta_sketch_a:+.5f}\nB: {t1.target_delta_sketch_b:+.5f}\n(SE={t1.target_sketch_bootstrap_se:.5f})",
            f"{h_val:+.5f}",
            (
                f"{t1.worst_neighbor_gene}\n{t1.worst_neighbor_lfc:+.5f}"
                if t1.worst_neighbor_gene
                else "[dim]n/a[/dim]"
            ),
            f"Splice: {t2.max_splice_junction_delta:.4f}\n"
            f"CpG: {t2.cpg_methylation_mean_abs_shift:.4f}\n"
            f"Bystanders: {n_bys}",
            f"{status_badge} {step.diagnosis}\n[#4a90e2]Rule: {step.constraint_added}[/#4a90e2]",
        )
    con.print(trace_tbl)

    # Why the search ended where it did: a weak optimiser, or a locus where efficacy and a
    # hazard are inseparable? These call for opposite next steps, so the run says which it is.
    hc = result.hazard_coupling
    if hc and hc.get("couplings"):
        lines: list[str] = []
        for c in hc["couplings"]:
            star = "[bold]" if c["spearman_p"] < 0.05 and c["spearman_rho"] > 0 else "[dim]"
            ctrl = (
                f" | random-edit control: rho={c['spearman_rho_null_only']:+.2f} "
                f"p={c['spearman_p_null_only']:.3g}"
                if c["n_null_only"]
                else ""
            )
            lines.append(
                f"{star}{c['hazard']}: Spearman rho={c['spearman_rho']:+.2f} "
                f"(p={c['spearman_p']:.3g}, n={c['n_total']}){ctrl}[/]"
            )
        verdict = str(hc.get("verdict", "inconclusive"))
        color = {"target_constrained": "bold #eb6834", "search_limited": "bold #4a90e2"}.get(verdict, "dim")
        lines.append("")
        lines.append(f"Verdict: [{color}]{verdict.upper().replace('_', ' ')}[/] -- {hc.get('verdict_explanation', '')}")
        con.print(
            Panel(
                "\n".join(lines),
                title="Efficacy-Hazard Coupling: is the limit the search, or the target?",
                border_style="dim",
            )
        )

    # Panel 3: Selected Candidate Summary
    if result.selected_candidate is not None:
        sel = result.selected_candidate
        c = sel.candidate
        summary = (
            f"Candidate ID: {c.candidate_id} | Converged: {result.converged}\n"
            f"Genomic Edit: {c.target_edit.chromosome}:{c.target_edit.position} {c.target_edit.reference}->{c.target_edit.alternate} "
            f"| Editor: {c.editor.value} ({c.strand} strand, protospacer pos {c.target_edit.protospacer_index})\n"
            f"Protospacer (5'->3'): {c.protospacer_20nt} | PAM: {c.pam_sequence} | Bystanders in [4..8]: {len(c.bystander_edits_in_window)}\n"
            f"Tier 1 (Single-Cell Conditioned `SC_GENE_EXPRESSION_DENOISED`): "
            f"Sketch A={sel.tier1_selectivity.target_delta_sketch_a:+.5f}, "
            f"Held-Out Sketch B={sel.tier1_selectivity.target_delta_sketch_b:+.5f} "
            f"(Healthy Control={next(iter(sel.tier1_selectivity.healthy_control_deltas.values()), 0.0):+.5f})\n"
            f"Tier 2 (Reference CURIEs {sel.tier2_hazards.reference_curies}): "
            f"SpliceJunctionScorer={sel.tier2_hazards.max_splice_junction_delta:.4f}, "
            f"MethylationScorer(1201)={sel.tier2_hazards.cpg_methylation_mean_abs_shift:.4f}\n"
            f"AlphaGenome Atlas Deep-Link: {c.target_edit.atlas_url}"
        )

        t2s = sel.tier2_hazards
        if t2s.reference_biosample_name:
            tiers = ", ".join(f"{k}={v}" for k, v in t2s.evidence_tiers.items())
            summary += (
                f"\nHazard reference: {t2s.reference_biosample_name} "
                f"[{t2s.limiting_evidence_tier.upper()}] ({tiers})"
            )
        if t2s.tissue_matched_agrees is not None:
            verdict = "AGREES" if t2s.tissue_matched_agrees else "DISAGREES"
            style = "#4a90e2" if t2s.tissue_matched_agrees else "bold #eb6834"
            summary += (
                f"\nTissue-matched cross-check ({t2s.tissue_matched_biosample_name}, weaker evidence): "
                f"[{style}]{verdict}[/{style}] "
                f"(splice={t2s.tissue_matched_splice_delta:.4f}, CpG={t2s.tissue_matched_cpg_shift:.4f})"
            )

        t1s = sel.tier1_selectivity
        if t1s.neighbor_gene_lfcs:
            ranked = sorted(t1s.neighbor_gene_lfcs.items(), key=lambda kv: -abs(kv[1]))[:4]
            neigh = ", ".join(f"{g} {v:+.4f}" for g, v in ranked)
            flag = "PASS" if t1s.neighbor_specificity_passed else "FAIL"
            summary += (
                f"\nNeighbour genes in interval ({len(t1s.neighbor_gene_lfcs)} scored, {flag}): {neigh}"
            )

        con.print(Panel(summary, title="Selected Candidate gRNA Specification (Research-Use-Only)", border_style="#4a90e2"))

        panel = sel.off_tissue_panel
        if panel:
            ot = Table(
                title=(
                    "Off-Tissue Safety Panel (Tier 2 reference conditions, NOT single-cell): "
                    f"bar is {panel.threshold:.5f} = 50% of the on-target effect"
                ),
                show_lines=False,
            )
            ot.add_column("Tissue")
            ot.add_column("|LFC|", justify="right")
            ot.add_column("RNA-seq evidence", justify="center")
            ot.add_column("Verdict", justify="center")
            for tissue, val in sorted(panel.tissues.items(), key=lambda kv: -kv[1]):
                ok = val <= panel.threshold
                ot.add_row(
                    tissue,
                    f"{val:.5f}",
                    panel.evidence_tiers.get(tissue, "?"),
                    "[#4a90e2]ok[/#4a90e2]" if ok else "[bold #eb6834]FLAG[/bold #eb6834]",
                )
            con.print(ot)
