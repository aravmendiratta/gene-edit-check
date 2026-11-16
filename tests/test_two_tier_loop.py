"""End-to-end test of the Two-Tier Verifier & Gradient-ISM search loop."""

from __future__ import annotations

from gene_edit_check.ag_backend.simulator import SurrogateSimulatorBackend
from gene_edit_check.search.search_loop import run_gradient_ism_search
from gene_edit_check.budget import RunBudgetTracker
from gene_edit_check.models import BackendMode
from gene_edit_check.scenarios.registry import SCENARIOS, build_til_cohort_anndata
from gene_edit_check.sc_sketch.sketch_variance import build_disjoint_sketches_and_ablation


def test_end_to_end_two_tier_search_loop() -> None:
    spec = SCENARIOS["pd1-exhaustion"]
    budget = RunBudgetTracker()
    backend = SurrogateSimulatorBackend(budget=budget, seed=42)
    cohort = build_til_cohort_anndata(n_cells_per_state=64, seed=42)

    t_adata = cohort[spec.target_state_name]
    t_emb = backend.encode_cells(t_adata, spec.target_state_name)
    t_sketches, t_metrics = build_disjoint_sketches_and_ablation(
        t_adata, t_emb, spec.target_state_name, sketch_size=16, seed=42
    )

    h_adata = cohort["Healthy_Naive_CD8"]
    h_emb = backend.encode_cells(h_adata, "Healthy_Naive_CD8")
    h_sketches, h_metrics = build_disjoint_sketches_and_ablation(
        h_adata, h_emb, "Healthy_Naive_CD8", sketch_size=16, seed=42
    )

    result = run_gradient_ism_search(
        scenario_id=spec.scenario_id,
        gene_id=spec.gene_id,
        gene_symbol=spec.gene_symbol,
        backend=backend,
        budget=budget,
        target_sketches=t_sketches,
        healthy_sketches_map={"Healthy_Naive_CD8": h_sketches},
        sketch_metrics=[t_metrics, h_metrics],
        reference_curies=list(spec.reference_curies),
        max_iterations=5,
        random_seed=42,
    )

    assert result.provenance.backend_mode == BackendMode.SURROGATE_SIMULATOR
    assert len(result.iterations) >= 1
    assert result.selected_candidate is not None
    assert (
        result.selected_candidate.tier1_selectivity.resolution_tier
        == "TIER1_SINGLE_CELL_CONDITIONED"
    )
    assert (
        result.selected_candidate.tier2_hazards.resolution_tier
        == "TIER2_ONTOLOGY_REFERENCE_CONDITION"
    )


def _run_with_regions(regions):  # type: ignore[no-untyped-def]
    spec = SCENARIOS["pd1-exhaustion"]
    budget = RunBudgetTracker()
    backend = SurrogateSimulatorBackend(budget=budget, seed=42)
    cohort = build_til_cohort_anndata(n_cells_per_state=64, seed=42)
    sk, mets = {}, []
    for name in (spec.target_state_name, "Healthy_Naive_CD8"):
        ad = cohort[name]
        s, m = build_disjoint_sketches_and_ablation(ad, backend.encode_cells(ad, name), name, sketch_size=16, seed=42)
        sk[name], _ = s, mets.append(m)
    return run_gradient_ism_search(
        scenario_id=spec.scenario_id, gene_id=spec.gene_id, gene_symbol=spec.gene_symbol, backend=backend,
        budget=budget, target_sketches=sk[spec.target_state_name],
        healthy_sketches_map={"Healthy_Naive_CD8": sk["Healthy_Naive_CD8"]}, sketch_metrics=mets,
        reference_curies=list(spec.reference_curies), max_iterations=2, random_seed=42, null_n=0,
        scan_regions=regions,
    )


def test_scan_region_confines_candidates_and_is_recorded() -> None:
    result = _run_with_regions([(1000, 2200)])
    tss = result.locus.start + (result.locus.end - result.locus.start) // 2
    for it in result.iterations:
        rel = it.report.candidate.target_edit.position - tss
        assert 950 <= rel <= 2250
    assert result.data_provenance["scan_regions_tss_relative_bp"] == [[1000, 2200]]


def test_scan_region_outside_the_interval_is_rejected() -> None:
    import pytest

    with pytest.raises(ValueError, match="outside"):
        _run_with_regions([(10_000_000, 10_000_100)])
