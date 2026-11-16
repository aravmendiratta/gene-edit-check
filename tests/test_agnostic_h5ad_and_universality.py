"""Tests for dataset-agnostic .h5ad loading, automatic reference resolution,
GRCh38 CpG sequence architecture, and multi-seed / cross-context CLI execution."""

from __future__ import annotations

from pathlib import Path

import anndata
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from typer.testing import CliRunner

from crispr_reflexion.ag_backend.protocol import GradientProfile
from crispr_reflexion.benchmarks.universality import (
    compute_local_cpg_methylation_correlation,
    compute_sequence_cpg_architecture,
    local_cpg_count_within_50bp,
)
from crispr_reflexion.cli import app
from crispr_reflexion.conditions import resolve_reference_pair
from crispr_reflexion.datasets.generic_h5ad import GenericH5adSpec, build_generic_cohort
from crispr_reflexion.models import (
    BackendMode,
    BaseEditorType,
    CandidateGuideEdit,
    GenomicInterval,
    OptimizationRunResult,
    ProvenanceRecord,
    ReflexionStepTrace,
    SingleNucleotideVariant,
    Tier1SelectivityResult,
    Tier2HazardResult,
    TwoTierVerificationReport,
)
from crispr_reflexion.pipeline import resolve_scenario_references
from crispr_reflexion.scenarios.registry import SCENARIOS


def test_cross_context_scenarios_registered_and_references_auto_resolved() -> None:
    expected = (
        "crc-pd1-mmrd",
        "crc-pd1-mmrp",
        "crc-pd1-cd4-cxcl13",
        "crc-pd1-gammadelta",
        "tnbc-pd1-exhaustion",
    )
    for sid in expected:
        assert sid in SCENARIOS
        primary, _tm = resolve_scenario_references(SCENARIOS[sid])
        assert len(primary) >= 1

    cd4_prim, cd4_tm = resolve_scenario_references(SCENARIOS["crc-pd1-cd4-cxcl13"])
    assert cd4_prim == ["UBERON:0000178", "CL:0000624"]
    assert cd4_tm == ["UBERON:0001155", "CL:0000624"]


def test_resolve_reference_pair_picks_highest_evidence_across_lineages() -> None:
    cd8 = resolve_reference_pair("CD8-positive", "colon")
    assert cd8.primary_curies == ("UBERON:0000178", "CL:0000625")
    assert cd8.primary_limiting_tier == "observed"
    assert cd8.tissue_matched_curies == ("UBERON:0001155", "CL:0000625")
    assert cd8.tissue_matched_limiting_tier == "imputed_exploratory"

    cd4 = resolve_reference_pair("CD4-positive, alpha-beta T cell", "colon")
    assert cd4.primary_curies == ("UBERON:0000178", "CL:0000624")
    assert cd4.primary_limiting_tier == "observed"
    assert cd4.tissue_matched_curies == ("UBERON:0001155", "CL:0000624")

    b_cell = resolve_reference_pair("B cell")
    assert len(b_cell.primary_curies) >= 1
    assert b_cell.primary_biosample_name

    hep = resolve_reference_pair("hepatocyte")
    assert hep.primary_curies == ("CL:0000182",)
    assert hep.primary_limiting_tier in {"observed", "imputed_supported"}


def _make_symbol_indexed_h5ad(
    path: Path,
    *,
    use_layers: bool = False,
    use_raw: bool = False,
    n_cells: int = 96,
) -> Path:
    rng = np.random.default_rng(123)
    genes = [
        ("PDCD1", "ENSG00000188389.12"),
        ("LAG3", "ENSG00000089692.8"),
        ("IL7R", "ENSG00000168685.15"),
        ("CD8A", "ENSG00000153563.4"),
        ("ACTB", "ENSG00000075624.18"),
    ]
    n_genes = len(genes)
    counts = rng.poisson(3, size=(n_cells, n_genes)).astype(np.int32)
    counts[: n_cells // 2, 0] += 8
    counts[: n_cells // 2, 1] += 6
    counts[n_cells // 2 :, 1] = 0
    counts[n_cells // 2 :, 2] += 7

    obs = pd.DataFrame(
        {
            "patient_id": [f"P{i % 8}" for i in range(n_cells)],
            "cluster": ["Tumor_Exh"] * (n_cells // 2) + ["Normal_Mem"] * (n_cells // 2),
        },
        index=[f"cell_{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(
        {
            "gene_symbol": [g[0] for g in genes],
            "ensembl_id": [g[1] for g in genes],
        },
        index=[g[0] for g in genes],
    )

    if use_layers:
        # Normalized float in .X, raw UMI integers in .layers["counts"]
        x_norm = np.log1p(counts.astype(np.float32))
        ad = anndata.AnnData(X=x_norm, obs=obs, var=var)
        ad.layers["counts"] = sp.csr_matrix(counts)
    elif use_raw:
        x_norm = np.log1p(counts.astype(np.float32))
        ad = anndata.AnnData(X=x_norm, obs=obs, var=var)
        ad.raw = anndata.AnnData(X=sp.csr_matrix(counts), obs=obs.copy(), var=var.copy())
    else:
        ad = anndata.AnnData(X=sp.csr_matrix(counts), obs=obs, var=var)

    ad.write_h5ad(path)
    return path


def test_build_generic_cohort_reads_layers_raw_and_maps_ensembl(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("crispr_reflexion.datasets.generic_h5ad.PROCESSED_DIR", tmp_path / "proc")

    h5_layers = _make_symbol_indexed_h5ad(tmp_path / "layers.h5ad", use_layers=True)
    spec_l = GenericH5adSpec(
        h5ad_path_or_url=h5_layers,
        dataset_name="test-layers",
        target_state_name="Exh",
        healthy_state_name="Mem",
        target_obs_query="cluster == 'Tumor_Exh'",
        healthy_obs_query="cluster == 'Normal_Mem'",
    )
    coh_l = build_generic_cohort(
        spec_l, pool_size=40, seed=42, marker_gene_ids=("ENSG00000188389",), use_cache=False
    )
    assert coh_l.provenance["count_matrix_source"] == "layers[counts]"
    assert coh_l.provenance["donor_column_used"] == "patient_id"
    assert list(coh_l.states["Exh"].var_names)[:2] == ["ENSG00000188389", "ENSG00000089692"]
    assert coh_l.observed_cp10k["Exh"]["ENSG00000188389"] > coh_l.observed_cp10k["Mem"]["ENSG00000188389"]

    h5_raw = _make_symbol_indexed_h5ad(tmp_path / "raw.h5ad", use_raw=True)
    spec_r = GenericH5adSpec(
        h5ad_path_or_url=h5_raw,
        dataset_name="test-raw-expr",
        target_state_name="Exh",
        healthy_state_name="Mem",
        target_obs_query="expr_LAG3 > 0",
        healthy_obs_query="expr_LAG3 == 0 and expr_IL7R > 0",
    )
    coh_r = build_generic_cohort(spec_r, pool_size=40, seed=42, use_cache=False)
    assert coh_r.provenance["count_matrix_source"] == "raw.X"
    assert coh_r.states["Exh"].n_obs > 0 and coh_r.states["Mem"].n_obs > 0


def test_compute_sequence_cpg_architecture_exact_counts() -> None:
    # 20 bp string: "CGCGCGCGCGAAAAAAAAAA" -> C=5, G=5, GC=50%, CG count=5 (positions 0,2,4,6,8)
    seq = "CGCGCGCGCGAAAAAAAAAA"
    arch = compute_sequence_cpg_architecture(seq, window_center=10, radius_bp=50, gene_symbol="TEST")
    assert arch.window_bp == 20
    assert arch.gc_fraction == pytest.approx(0.5)
    assert arch.cpg_dinucleotide_count == 5
    assert arch.cpg_per_100bp == pytest.approx(25.0)
    # O/E = (5 * 20) / (5 * 5) = 4.0
    assert arch.cpg_observed_over_expected == pytest.approx(4.0)
    assert local_cpg_count_within_50bp(seq, rel_idx=4, half_window=50) == 5


def test_cli_run_h5ad_simulator_smoke(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("crispr_reflexion.datasets.generic_h5ad.PROCESSED_DIR", tmp_path / "proc")
    h5_path = _make_symbol_indexed_h5ad(tmp_path / "cli_demo.h5ad", n_cells=120)
    out_dir = tmp_path / "out"

    runner = CliRunner()
    res = runner.invoke(
        app,
        [
            "run-h5ad",
            "--h5ad", str(h5_path),
            "--gene-id", "ENSG00000188389",
            "--gene-symbol", "PDCD1",
            "--target-query", "cluster == 'Tumor_Exh'",
            "--healthy-query", "cluster == 'Normal_Mem'",
            "--cell-type", "CD8-positive",
            "--tissue", "colon",
            "--sketch-size", "12",
            "--pool-size", "50",
            "--max-iterations", "2",
            "--null-n", "4",
            "--simulator",
            "--output-dir", str(out_dir),
        ],
    )
    assert res.exit_code == 0, res.output
    assert (out_dir / "report_custom-pdcd1-h5ad.json").exists()
    assert (out_dir / "figure_custom-pdcd1-h5ad.png").exists()


def test_resolve_reference_pair_handles_no_tissue_and_fallback() -> None:
    res_no_tissue = resolve_reference_pair("CD8-positive")
    assert res_no_tissue.primary_curies == ("UBERON:0000178", "CL:0000625")
    assert res_no_tissue.tissue_matched_curies == ()

    res_fallback = resolve_reference_pair("nonexistent_cell_type_xyz_123")
    assert len(res_fallback.primary_curies) >= 1


def test_build_generic_cohort_rejects_float_normalized_x(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("crispr_reflexion.datasets.generic_h5ad.PROCESSED_DIR", tmp_path / "proc")
    ad = anndata.AnnData(
        X=np.array([[1.23, 0.45], [0.12, 2.34]], dtype=np.float32),
        obs=pd.DataFrame({"donor_id": ["d1", "d2"], "cluster": ["A", "B"]}, index=["c1", "c2"]),
        var=pd.DataFrame(index=["ENSG00000188389", "ENSG00000089692"]),
    )
    bad_path = tmp_path / "float_only.h5ad"
    ad.write_h5ad(bad_path)
    spec = GenericH5adSpec(
        h5ad_path_or_url=bad_path,
        dataset_name="bad-float",
        target_state_name="A",
        healthy_state_name="B",
        target_obs_query="cluster == 'A'",
        healthy_obs_query="cluster == 'B'",
    )
    with pytest.raises(ValueError, match="integer UMI counts"):
        build_generic_cohort(spec, use_cache=False)


def test_local_cpg_methylation_correlation() -> None:
    seq = ("CG" * 20) + ("AT" * 50) + ("AAAA" * 30)
    iv = GenomicInterval(chromosome="chr2", start=1000, end=1000 + len(seq))
    gp = GradientProfile(
        gene_id="ENSG0",
        gene_name="PDCD1",
        condition_name="T",
        interval=iv,
        reference_sequence=seq,
        signed_gradients_1bp=np.zeros(len(seq), dtype=np.float32),
        scalar_score=1.0,
    )
    positions_and_meth = [(1020, 0.80), (1050, 0.35), (1080, 0.10), (1180, 0.01)]
    steps = []
    for idx, (pos1, cpg_shift) in enumerate(positions_and_meth, start=1):
        v = SingleNucleotideVariant(chromosome="chr2", position=pos1, reference="A", alternate="G", protospacer_index=5)
        cand = CandidateGuideEdit(
            candidate_id=f"c{idx}", gene_symbol="PDCD1", editor=BaseEditorType.ABE8E, strand="+",
            protospacer_20nt="A" * 20, pam_sequence="AGG", genomic_start=pos1 - 10, genomic_end=pos1 + 13, target_edit=v,
        )
        rep = TwoTierVerificationReport(
            candidate=cand,
            tier1_selectivity=Tier1SelectivityResult(
                target_condition_name="T", target_delta_sketch_a=-0.1, target_delta_sketch_b=-0.1,
                target_sketch_bootstrap_se=0.0, target_efficacy_passed=True,
                healthy_selectivity_passed=True, cross_sketch_generalization_passed=True,
            ),
            tier2_hazards=Tier2HazardResult(
                reference_curies=["X"], max_splice_junction_delta=0.0,
                cpg_methylation_mean_abs_shift=cpg_shift, splicing_hazard_passed=True,
                methylation_hazard_passed=True, bystander_clean_passed=True,
            ),
            composite_objective_score=0.0, all_checks_passed=True,
        )
        steps.append(ReflexionStepTrace(iteration=idx, report=rep, critic_diagnosis="", episodic_memory_rule_added=""))

    run_res = OptimizationRunResult(
        scenario_id="s", gene_symbol="PDCD1", locus=iv,
        provenance=ProvenanceRecord(
            backend_mode=BackendMode.SURROGATE_SIMULATOR, model_version="v", random_seed=42,
            total_api_calls=0, cache_hits=0, total_ism_variants_scored=0, gradient_calls_used=0,
        ),
        sketch_metrics=[], iterations=steps, converged=True,
    )
    rho, p = compute_local_cpg_methylation_correlation(gp, run_res)
    assert rho is not None and rho > 0.9
    assert p is not None and p < 0.1
