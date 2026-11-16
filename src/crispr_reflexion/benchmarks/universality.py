"""Sequence-Intrinsic CpG Architecture & Multi-Seed / Cross-Context Universality (`benchmarks/universality.py`).

Three tests that separate a physical regulatory mechanism from a single-dataset fluke:
1. Physical GRCh38 CpG Architecture: counts raw `CG` dinucleotides, GC fraction, and
   Gardiner-Garden & Frommer CpG observed/expected ratio directly on the recovered
   reference DNA sequence (zero dependence on AlphaGenome 2's output heads), and tests
   whether local `+-50 bp` CpG density tracks AlphaGenome 2's `MethylationScorer` shift.
2. Multi-Seed Stability: re-runs the pipeline across multiple random seeds (changing both
   the patient-disjoint Sketch A / Sketch B draw and the random PAM-valid null edits).
3. Cross-Cancer, Cross-Subtype, and Cross-Lineage Invariance: evaluates the same
   target locus across tumor subtypes (MMRd vs MMRp CRC), immune lineages (CD8+,
   CD4+ CXCL13+, gamma-delta T cells), and an independent cancer cohort.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from pydantic import BaseModel, Field
from scipy import stats

from crispr_reflexion.ag_backend.protocol import GradientProfile, VerifierBackend
from crispr_reflexion.budget import RunBudgetTracker
from crispr_reflexion.models import OptimizationRunResult


class SequenceCpGArchitecture(BaseModel):
    """Physical sequence composition measured directly on recovered GRCh38 DNA."""

    gene_symbol: str = ""
    window_bp: int
    gc_fraction: float
    cpg_dinucleotide_count: int
    cpg_per_100bp: float
    cpg_observed_over_expected: float
    local_cpg_vs_methylation_rho: float | None = None
    local_cpg_vs_methylation_p: float | None = None


class SeedRunSummary(BaseModel):
    seed: int
    converged: bool
    meth_rho: float | None = None
    meth_p: float | None = None
    meth_rho_null_only: float | None = None
    splice_rho: float | None = None
    meth_minus_splice_rho: float | None = None
    coupling_class: str = "unknown"
    verdict: str = "inconclusive"
    gradient_sign_agreement: float | None = None
    best_abs_lfc: float = 0.0
    mmd_reduction_target_pct: float = 0.0


class ContextInvarianceSummary(BaseModel):
    scenario_id: str
    context_label: str
    dataset_name: str
    gene_symbol: str
    target_state: str
    healthy_state: str
    target_cells_full: int = 0
    target_donors_full: int = 0
    primary_reference: str = ""
    primary_evidence_tier: str = "unknown"
    converged: bool
    meth_rho: float | None = None
    meth_p: float | None = None
    meth_rho_null_only: float | None = None
    splice_rho: float | None = None
    meth_minus_splice_rho: float | None = None
    coupling_class: str = "unknown"
    verdict: str = "inconclusive"
    binding_constraint: str | None = None
    n_failing_binding: int = 0
    n_iterations: int = 0


class UniversalitySweepReport(BaseModel):
    gene_symbol: str
    seeds_tested: list[int]
    seed_summaries: list[SeedRunSummary] = Field(default_factory=list)
    seed_meth_rho_mean: float | None = None
    seed_meth_rho_sd: float | None = None
    seed_null_rho_mean: float | None = None
    context_summaries: list[ContextInvarianceSummary] = Field(default_factory=list)
    locus_cpg_architecture: list[dict[str, Any]] = Field(default_factory=list)


def compute_sequence_cpg_architecture(
    reference_seq: str,
    window_center: int | None = None,
    radius_bp: int = 2500,
    gene_symbol: str = "",
) -> SequenceCpGArchitecture:
    """Measure GC fraction, CpG count/100bp, and Gardiner-Garden CpG O/E on DNA."""
    seq_u = reference_seq.upper()
    if window_center is None:
        window_center = len(seq_u) // 2
    lo = max(0, window_center - radius_bp)
    hi = min(len(seq_u), window_center + radius_bp)
    win = seq_u[lo:hi]

    valid_len = max(1, sum(1 for ch in win if ch in ("A", "C", "G", "T")))
    n_c = win.count("C")
    n_g = win.count("G")
    n_cpg = sum(1 for i in range(len(win) - 1) if win[i] == "C" and win[i + 1] == "G")

    gc_frac = (n_c + n_g) / valid_len
    cpg_100 = 100.0 * n_cpg / valid_len
    expected = (n_c * n_g) / valid_len
    obs_exp = (n_cpg / expected) if expected > 0.0 else 0.0

    return SequenceCpGArchitecture(
        gene_symbol=gene_symbol,
        window_bp=len(win),
        gc_fraction=round(gc_frac, 4),
        cpg_dinucleotide_count=n_cpg,
        cpg_per_100bp=round(cpg_100, 4),
        cpg_observed_over_expected=round(obs_exp, 4),
    )


def local_cpg_count_within_50bp(reference_seq: str, rel_idx: int, half_window: int = 50) -> int:
    """Count physical `CG` dinucleotides within `+-half_window` bp of a relative index."""
    lo = max(0, rel_idx - half_window)
    hi = min(len(reference_seq), rel_idx + half_window + 1)
    sub = reference_seq[lo:hi].upper()
    return sum(1 for i in range(len(sub) - 1) if sub[i] == "C" and sub[i + 1] == "G")


def compute_local_cpg_methylation_correlation(
    grad_profile: GradientProfile,
    result: OptimizationRunResult,
) -> tuple[float | None, float | None]:
    """Correlate physical `+-50 bp` GRCh38 CpG density with AlphaGenome 2 MethylationScorer."""
    seq = grad_profile.reference_sequence
    start0 = grad_profile.interval.start
    cpg_counts: list[int] = []
    meth_shifts: list[float] = []
    for step in result.iterations:
        rel_idx = step.report.candidate.target_edit.position - 1 - start0
        if 0 <= rel_idx < len(seq):
            cpg_counts.append(local_cpg_count_within_50bp(seq, rel_idx))
            meth_shifts.append(step.report.tier2_hazards.cpg_methylation_mean_abs_shift)

    if len(cpg_counts) < 3 or np.ptp(cpg_counts) == 0 or np.ptp(meth_shifts) == 0:
        return None, None
    rho, p = stats.spearmanr(cpg_counts, meth_shifts)
    if not np.isfinite(rho) or not np.isfinite(p):
        return None, None
    return round(float(rho), 4), float(f"{p:.3g}")


def _extract_coupling_pair(result: OptimizationRunResult) -> tuple[float | None, float | None, float | None, float | None]:
    hc = result.hazard_coupling or {}
    couplings = hc.get("couplings", [])
    meth = next((c for c in couplings if c.get("hazard") == "cpg_methylation"), None)
    splice = next((c for c in couplings if c.get("hazard") == "splice_junction"), None)
    return (
        meth.get("spearman_rho") if meth else None,
        meth.get("spearman_p") if meth else None,
        meth.get("spearman_rho_null_only") if meth else None,
        splice.get("spearman_rho") if splice else None,
    )


def summarize_seed_run(result: OptimizationRunResult, seed: int) -> SeedRunSummary:
    meth_rho, meth_p, null_rho, splice_rho = _extract_coupling_pair(result)
    hc = result.hazard_coupling or {}
    best = max(
        (abs(s.report.tier1_selectivity.target_delta_sketch_b) for s in result.iterations),
        default=0.0,
    )
    target_mmd = result.sketch_metrics[0].mmd_reduction_pct if result.sketch_metrics else 0.0
    return SeedRunSummary(
        seed=seed,
        converged=result.converged,
        meth_rho=meth_rho,
        meth_p=meth_p,
        meth_rho_null_only=null_rho,
        splice_rho=splice_rho,
        meth_minus_splice_rho=hc.get("methylation_minus_splice_rho"),
        coupling_class=str(hc.get("coupling_class", "unknown")),
        verdict=str(hc.get("verdict", "inconclusive")),
        gradient_sign_agreement=(
            result.ism_prescreen.gradient_sign_agreement if result.ism_prescreen else None
        ),
        best_abs_lfc=round(best, 6),
        mmd_reduction_target_pct=round(target_mmd, 2),
    )


def summarize_context_run(
    result: OptimizationRunResult,
    context_label: str,
) -> ContextInvarianceSummary:
    meth_rho, meth_p, null_rho, splice_rho = _extract_coupling_pair(result)
    hc = result.hazard_coupling or {}
    dp = result.data_provenance or {}
    states = dp.get("states") or {}
    target_sm = result.sketch_metrics[0] if result.sketch_metrics else None
    healthy_sm = result.sketch_metrics[1] if len(result.sketch_metrics) > 1 else None
    t_name = target_sm.cell_state_name if target_sm else ""
    h_name = healthy_sm.cell_state_name if healthy_sm else ""
    t_info = states.get(t_name, {}) if isinstance(states, dict) else {}

    first_t2 = result.iterations[0].report.tier2_hazards if result.iterations else None
    return ContextInvarianceSummary(
        scenario_id=result.scenario_id,
        context_label=context_label,
        dataset_name=str(dp.get("data_source", result.scenario_id)),
        gene_symbol=result.gene_symbol,
        target_state=t_name,
        healthy_state=h_name,
        target_cells_full=int(t_info.get("n_cells_full_population", 0)),
        target_donors_full=int(t_info.get("n_donors_full_population", 0)),
        primary_reference=first_t2.reference_biosample_name if first_t2 else "",
        primary_evidence_tier=first_t2.limiting_evidence_tier if first_t2 else "unknown",
        converged=result.converged,
        meth_rho=meth_rho,
        meth_p=meth_p,
        meth_rho_null_only=null_rho,
        splice_rho=splice_rho,
        meth_minus_splice_rho=hc.get("methylation_minus_splice_rho"),
        coupling_class=str(hc.get("coupling_class", "unknown")),
        verdict=str(hc.get("verdict", "inconclusive")),
        binding_constraint=hc.get("binding_constraint"),
        n_failing_binding=int(hc.get("n_candidates_failing_binding_constraint", 0)),
        n_iterations=len(result.iterations),
    )


def run_multiseed_and_context_invariance(
    scenario_specs: dict[str, Any],
    backend_factory: Any,
    run_scenario_fn: Any,
    load_cohort_fn: Any,
    build_sketches_fn: Any,
    primary_scenario_id: str = "crc-pd1-exhaustion",
    context_labels: dict[str, str] | None = None,
    sweep_genes: tuple[tuple[str, str, str], ...] = (),
    seeds: tuple[int, ...] = (42, 7, 99),
    sketch_size: int = 128,
    pool_size: int = 600,
    max_iterations: int = 5,
    null_n: int = 10,
    prescreen_windows: int = 6,
) -> UniversalitySweepReport:
    """Execute multi-seed stability, cross-context invariance, and CpG sequence architecture."""
    labels = context_labels or {}
    primary_spec = scenario_specs[primary_scenario_id]
    run_cache: dict[tuple[str, int], OptimizationRunResult] = {}

    def _get_or_run(spec: Any, run_seed: int) -> OptimizationRunResult:
        cache_key = (str(spec.scenario_id), int(run_seed))
        if cache_key in run_cache:
            return run_cache[cache_key]
        budget = RunBudgetTracker()
        backend = backend_factory(budget=budget, seed=run_seed)
        res = run_scenario_fn(
            spec,
            backend,
            budget,
            sketch_size=sketch_size,
            max_iterations=max_iterations,
            seed=run_seed,
            pool_size=pool_size,
            null_n=null_n,
            ism_prescreen_windows=prescreen_windows,
            off_tissue_panel=False,
        )
        run_cache[cache_key] = res
        return res

    seed_summaries: list[SeedRunSummary] = []
    for s in seeds:
        print(f"[sweep-universality] Part A seed={s} ({primary_scenario_id})...", flush=True)
        res = _get_or_run(primary_spec, s)
        summary = summarize_seed_run(res, seed=s)
        seed_summaries.append(summary)
        print(f"  -> seed={s}: CpG rho={summary.meth_rho}, class={summary.coupling_class}", flush=True)

    valid_rhos = [r.meth_rho for r in seed_summaries if r.meth_rho is not None]
    valid_null_rhos = [r.meth_rho_null_only for r in seed_summaries if r.meth_rho_null_only is not None]

    context_summaries: list[ContextInvarianceSummary] = []
    for sid, label in labels.items():
        if sid not in scenario_specs:
            continue
        print(f"[sweep-universality] Part B context={sid} ({label})...", flush=True)
        spec = scenario_specs[sid]
        res = _get_or_run(spec, seeds[0])
        ctx_sum = summarize_context_run(res, context_label=label)
        context_summaries.append(ctx_sum)
        print(f"  -> {sid}: CpG rho={ctx_sum.meth_rho}, class={ctx_sum.coupling_class}", flush=True)

    locus_arch_rows: list[dict[str, Any]] = []
    if sweep_genes:
        cohort = load_cohort_fn(primary_spec, pool_size=pool_size, seed=seeds[0])
        budget_sk = RunBudgetTracker()
        backend_sk = backend_factory(budget=budget_sk, seed=seeds[0])
        target_sk, _, _ = build_sketches_fn(backend_sk, primary_spec, cohort, sketch_size, seeds[0])

        for gid, sym, _sid in sweep_genes:
            print(f"[sweep-universality] Part C locus CpG architecture={sym} ({_sid})...", flush=True)
            b_l = RunBudgetTracker()
            be_l = backend_factory(budget=b_l, seed=seeds[0])
            grad: GradientProfile = be_l.compute_interval_gradients(
                gene_id=gid,
                gene_symbol=sym,
                cell_sketch=target_sk.sketch_a,
                condition_name=target_sk.cell_state_name,
            )
            arch = compute_sequence_cpg_architecture(
                grad.reference_sequence,
                gene_symbol=sym,
            )
            row = arch.model_dump()
            if _sid in scenario_specs:
                r_locus = _get_or_run(scenario_specs[_sid], seeds[0])
                l_rho, l_p = compute_local_cpg_methylation_correlation(grad, r_locus)
                row["local_cpg_vs_methylation_rho"] = l_rho
                row["local_cpg_vs_methylation_p"] = l_p
                hc = r_locus.hazard_coupling or {}
                meth = next((c for c in hc.get("couplings", []) if c.get("hazard") == "cpg_methylation"), None)
                row["efficacy_vs_cpg_rho"] = meth.get("spearman_rho") if meth else None
                row["coupling_class"] = hc.get("coupling_class")
            locus_arch_rows.append(row)

    return UniversalitySweepReport(
        gene_symbol=primary_spec.gene_symbol,
        seeds_tested=list(seeds),
        seed_summaries=seed_summaries,
        seed_meth_rho_mean=round(float(np.mean(valid_rhos)), 4) if valid_rhos else None,
        seed_meth_rho_sd=round(float(np.std(valid_rhos, ddof=1)), 4) if len(valid_rhos) > 1 else 0.0,
        seed_null_rho_mean=round(float(np.mean(valid_null_rhos)), 4) if valid_null_rhos else None,
        context_summaries=context_summaries,
        locus_cpg_architecture=locus_arch_rows,
    )
