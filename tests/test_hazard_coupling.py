"""The coupling analysis must tell a constrained target apart from a weak search."""

from __future__ import annotations

from gene_edit_check.benchmarks.hazard_coupling import analyze_hazard_coupling
from gene_edit_check.models import (
    BaseEditorType,
    CandidateGuideEdit,
    EfficacyNullCalibration,
    GenomicInterval,
    OptimizationRunResult,
    ProvenanceRecord,
    BackendMode,
    ReflexionStepTrace,
    SingleNucleotideVariant,
    Tier1SelectivityResult,
    Tier2HazardResult,
    TwoTierVerificationReport,
)


def _step(i: int, lfc: float, cpg: float, splice: float, reasons: list[str]) -> ReflexionStepTrace:
    v = SingleNucleotideVariant(chromosome="chr2", position=1000 + i, reference="A", alternate="G", protospacer_index=5)
    cand = CandidateGuideEdit(
        candidate_id=f"c{i}", gene_symbol="PDCD1", editor=BaseEditorType.ABE8E, strand="+",
        protospacer_20nt="A" * 20, pam_sequence="AGG", genomic_start=990, genomic_end=1013, target_edit=v,
    )
    rep = TwoTierVerificationReport(
        candidate=cand,
        tier1_selectivity=Tier1SelectivityResult(
            target_condition_name="T", target_delta_sketch_a=lfc, target_delta_sketch_b=lfc,
            target_sketch_bootstrap_se=0.0, target_efficacy_passed=True,
            healthy_selectivity_passed=True, cross_sketch_generalization_passed=True,
        ),
        tier2_hazards=Tier2HazardResult(
            reference_curies=["X"], max_splice_junction_delta=splice,
            cpg_methylation_mean_abs_shift=cpg, splicing_hazard_passed=splice <= 0.15,
            methylation_hazard_passed=cpg <= 0.12, bystander_clean_passed=True,
        ),
        composite_objective_score=0.0, all_checks_passed=not reasons, failure_reasons=reasons,
    )
    return ReflexionStepTrace(iteration=i, report=rep, critic_diagnosis="", episodic_memory_rule_added="")


def _result(steps: list[ReflexionStepTrace], null_lfc: list[float], null_cpg: list[float], null_splice: list[float]) -> OptimizationRunResult:
    return OptimizationRunResult(
        scenario_id="s", gene_symbol="PDCD1",
        locus=GenomicInterval(chromosome="chr2", start=0, end=10_000),
        provenance=ProvenanceRecord(
            backend_mode=BackendMode.SURROGATE_SIMULATOR, model_version="v", random_seed=42,
            total_api_calls=0, cache_hits=0, total_ism_variants_scored=0, gradient_calls_used=0,
        ),
        sketch_metrics=[], iterations=steps, converged=False,
        efficacy_null=EfficacyNullCalibration(
            n_random_edits=len(null_lfc), abs_lfc_values=null_lfc, abs_lfc_max=max(null_lfc),
            abs_lfc_mean=sum(null_lfc) / len(null_lfc), cutoff_used=max(null_lfc),
            empirical_p_floor=1 / (len(null_lfc) + 1),
            cpg_methylation_shift_values=null_cpg, splice_delta_values=null_splice,
        ),
    )


def test_coupled_hazard_yields_target_constrained_verdict() -> None:
    """Methylation rises with efficacy and blocks almost everything -> the locus is the limit."""
    steps = [
        _step(i, -0.05 * i, 0.10 * i, 0.01, ["TIER2_METHYLATION_HAZARD: x"] if i > 1 else [])
        for i in range(1, 7)
    ]
    res = _result(steps, [0.01, 0.02, 0.03, 0.04, 0.05], [0.02, 0.04, 0.06, 0.08, 0.10], [0.01] * 5)
    rep = analyze_hazard_coupling(res)

    assert rep.verdict == "target_constrained"
    assert rep.binding_constraint == "cpg_methylation"
    meth = next(c for c in rep.couplings if c.hazard == "cpg_methylation")
    assert meth.spearman_rho > 0.9 and meth.is_significant
    assert meth.survives_null_control  # holds among randomly chosen edits too
    assert "not an artifact" in meth.interpretation


def test_uncoupled_hazard_yields_search_limited_verdict() -> None:
    """The same failure rate, but methylation is unrelated to efficacy -> keep searching."""
    cpgs = [0.30, 0.01, 0.31, 0.02, 0.32, 0.03]
    steps = [
        _step(i + 1, -0.05 * (i + 1), cpgs[i], 0.01, ["TIER2_METHYLATION_HAZARD: x"])
        for i in range(6)
    ]
    res = _result(steps, [0.01, 0.02, 0.03, 0.04, 0.05], [0.3, 0.01, 0.3, 0.01, 0.3], [0.01] * 5)
    rep = analyze_hazard_coupling(res)

    assert rep.verdict == "search_limited"
    assert "search problem" in rep.verdict_explanation


def test_uncoupled_modality_is_reported_as_independent() -> None:
    """Splicing held flat must not be reported as coupled -- this is the control that rules out
    'the model's output heads simply move together'."""
    steps = [_step(i, -0.05 * i, 0.10 * i, 0.01, []) for i in range(1, 7)]
    res = _result(steps, [0.01, 0.02, 0.03], [0.02, 0.04, 0.06], [0.01, 0.01, 0.01])
    rep = analyze_hazard_coupling(res)
    splice = next(c for c in rep.couplings if c.hazard == "splice_junction")
    assert not splice.is_significant
    assert "independent" in splice.interpretation


def test_mixed_failures_are_inconclusive() -> None:
    steps = [
        _step(1, -0.1, 0.01, 0.01, ["TIER1_SELECTIVITY: x"]),
        _step(2, -0.2, 0.02, 0.01, ["TIER2_SPLICE_HAZARD: x"]),
        _step(3, -0.3, 0.03, 0.01, ["TIER1_DIRECTION: x"]),
        _step(4, -0.4, 0.04, 0.01, ["TIER2_BYSTANDER_HAZARD: x"]),
    ]
    res = _result(steps, [0.01, 0.02, 0.03], [0.02, 0.04, 0.06], [0.01, 0.02, 0.03])
    assert analyze_hazard_coupling(res).verdict == "inconclusive"


def test_handles_a_run_with_no_null_calibration() -> None:
    steps = [_step(i, -0.05 * i, 0.1 * i, 0.01, []) for i in range(1, 5)]
    res = _result(steps, [0.01, 0.02], [0.02, 0.04], [0.01, 0.01])
    res.efficacy_null = None
    rep = analyze_hazard_coupling(res)  # must not raise
    assert rep.couplings and all(c.n_null_only == 0 for c in rep.couplings)


def test_bootstrap_interval_brackets_a_known_correlation_and_is_deterministic() -> None:
    import numpy as np

    from gene_edit_check.benchmarks.hazard_coupling import bootstrap_rho_gap_ci, bootstrap_spearman_ci

    rng = np.random.default_rng(0)
    x = rng.normal(size=60)
    y = x + rng.normal(scale=0.5, size=60)  # strongly related
    lo, hi = bootstrap_spearman_ci(x, y, n_boot=500, seed=3)
    assert 0.0 < lo < hi <= 1.0
    assert (lo, hi) == bootstrap_spearman_ci(x, y, n_boot=500, seed=3)

    noise = rng.normal(size=60)  # unrelated: interval must straddle zero
    lo_n, hi_n = bootstrap_spearman_ci(x, noise, n_boot=500, seed=3)
    assert lo_n < 0.0 < hi_n

    # Gap between a related and an unrelated hazard should be clearly positive.
    g_lo, g_hi = bootstrap_rho_gap_ci(x, y, noise, n_boot=500, seed=3)
    assert g_lo > 0.0
    assert bootstrap_spearman_ci(np.ones(10), np.arange(10.0)) == (0.0, 0.0)  # degenerate input
