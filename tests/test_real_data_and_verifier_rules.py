"""Tests for the donor-held-out split, verifier direction/null rules, sparse hashing and real data."""

from __future__ import annotations

import anndata
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from crispr_reflexion.adata_utils import counts_fingerprint, drop_unexpressed_genes
from crispr_reflexion.ag_backend.protocol import VariantMultiModalScore
from crispr_reflexion.models import (
    BackendMode,
    BaseEditorType,
    CandidateGuideEdit,
    GenomicInterval,
    SingleNucleotideVariant,
)
from crispr_reflexion.sc_sketch.coreset import DisjointCellSketches
from crispr_reflexion.sc_sketch.sketch_variance import build_disjoint_sketches_and_ablation
from crispr_reflexion.verifier.two_tier import VerificationThresholds, evaluate_candidate_two_tier


def _adata(n_cells: int = 120, n_genes: int = 30, n_donors: int = 10, seed: int = 0) -> anndata.AnnData:
    rng = np.random.default_rng(seed)
    x = sp.csr_matrix(rng.poisson(0.7, size=(n_cells, n_genes)).astype(np.int32))
    obs = pd.DataFrame({"donor_id": [f"d{i % n_donors}" for i in range(n_cells)]}, index=[f"c{i}" for i in range(n_cells)])
    return anndata.AnnData(X=x, obs=obs, var=pd.DataFrame(index=[f"ENSG{i:011d}" for i in range(n_genes)]))


def test_sketch_split_is_held_out_by_donor() -> None:
    a = _adata(n_donors=10)
    emb = np.random.default_rng(1).normal(size=(a.n_obs, 16)).astype(np.float32)
    sk, m = build_disjoint_sketches_and_ablation(a, emb, "S", sketch_size=12, seed=3)
    assert m.split_by == "donor"
    assert m.donors_overlap_a_b == 0
    assert set(sk.sketch_a.obs["donor_id"]).isdisjoint(set(sk.sketch_b.obs["donor_id"]))
    assert set(sk.sketch_a.obs_names).isdisjoint(set(sk.sketch_b.obs_names))


def test_few_donors_falls_back_to_cell_split_and_says_so() -> None:
    a = _adata(n_donors=3)
    emb = np.random.default_rng(1).normal(size=(a.n_obs, 16)).astype(np.float32)
    _, m = build_disjoint_sketches_and_ablation(a, emb, "S", sketch_size=8, seed=3)
    assert m.split_by == "cell"


def test_fingerprint_is_sparse_dense_consistent_and_content_sensitive() -> None:
    a = _adata()
    dense = anndata.AnnData(X=a.X.toarray(), obs=a.obs.copy(), var=a.var.copy())
    assert counts_fingerprint(a) == counts_fingerprint(dense)
    arr = a.X.toarray()
    arr[0, 0] += 1
    b = anndata.AnnData(X=sp.csr_matrix(arr), obs=a.obs.copy(), var=a.var.copy())
    assert counts_fingerprint(a) != counts_fingerprint(b)


def test_drop_unexpressed_genes_preserves_library_size() -> None:
    a = _adata(n_genes=40)
    arr = a.X.toarray()
    arr[:, :10] = 0
    a = anndata.AnnData(X=sp.csr_matrix(arr), obs=a.obs.copy(), var=a.var.copy())
    slim = drop_unexpressed_genes(a)
    assert slim.n_vars <= 30
    assert np.allclose(np.asarray(slim.X.sum(1)), np.asarray(a.X.sum(1)))


class _FakeBackend:
    """Returns fixed LFCs: target sketches -> `target_lfc`, healthy sketches -> `healthy_lfc`."""

    backend_mode = BackendMode.SURROGATE_SIMULATOR
    model_version = "fake"

    def __init__(self, target_lfc: float, healthy_lfc: float) -> None:
        self.t, self.h = target_lfc, healthy_lfc

    def score_variant_multimodal(self, parent_interval, variant, gene_id, cell_sketch, condition_name, reference_curies=None):  # type: ignore[no-untyped-def]
        lfc = self.t if "target" in condition_name.lower() else self.h
        return VariantMultiModalScore(variant, condition_name, lfc, 0.001, 0.001, 0.001)


def _candidate() -> CandidateGuideEdit:
    v = SingleNucleotideVariant(chromosome="chr1", position=1000, reference="A", alternate="G", protospacer_index=6)
    return CandidateGuideEdit(
        candidate_id="c", gene_symbol="G", editor=BaseEditorType.ABE8E, strand="+",
        protospacer_20nt="A" * 20, pam_sequence="NGG", genomic_start=990, genomic_end=1013, target_edit=v,
    )


def _sketches(name: str) -> DisjointCellSketches:
    a = _adata(n_cells=8, n_donors=2)
    emb = np.zeros((8, 4), dtype=np.float32)
    return DisjointCellSketches(name, a, a, emb, emb, emb, True)


def _run(target_lfc: float, healthy_lfc: float, **thr_kwargs):  # type: ignore[no-untyped-def]
    return evaluate_candidate_two_tier(
        backend=_FakeBackend(target_lfc, healthy_lfc),  # type: ignore[arg-type]
        parent_interval=GenomicInterval(chromosome="chr1", start=0, end=2000),
        candidate=_candidate(), gene_id="ENSG0", target_sketches=_sketches("target"),
        healthy_sketches_map={"healthy": _sketches("healthy")}, reference_curies=["X"],
        thresholds=VerificationThresholds(**thr_kwargs),
    )


def test_silencing_goal_accepts_a_selective_decrease() -> None:
    assert _run(-0.20, -0.01, desired_direction="down").all_checks_passed


def test_silencing_goal_rejects_a_large_increase() -> None:
    # Old behaviour used |LFC| and would have PASSED this. Raising PD-1 is the opposite of the goal.
    rep = _run(+0.20, +0.01, desired_direction="down")
    assert not rep.all_checks_passed
    assert any(r.startswith("TIER1_DIRECTION") for r in rep.failure_reasons)


def test_null_cutoff_raises_the_efficacy_bar() -> None:
    assert _run(-0.05, -0.001, desired_direction="down", null_abs_lfc_cutoff=0.0).all_checks_passed
    rep = _run(-0.05, -0.001, desired_direction="down", null_abs_lfc_cutoff=0.08)
    assert not rep.all_checks_passed
    assert any(r.startswith("TIER1_EFFICACY") for r in rep.failure_reasons)


def test_selectivity_ratio_default_is_strict() -> None:
    # healthy shift is 50% of target shift: passed the old 0.55 limit, fails the pre-registered 0.35.
    rep = _run(-0.20, -0.10, desired_direction="down")
    assert not rep.tier1_selectivity.healthy_selectivity_passed


def test_real_cohort_has_the_expected_biology() -> None:
    from crispr_reflexion.datasets import crc_pelka as crc

    if not (crc.PROCESSED_DIR / "pelka_crc_cd8_pool600_seed42.h5ad").exists() and not crc.RAW_PATH.exists():
        pytest.skip("real cohort not downloaded (run: crispr-reflexion ablate-sketch)")
    pdcd1 = "ENSG00000188389"
    coh = crc.build_cohort(pool_size=600, seed=42, marker_gene_ids=(pdcd1,))
    obs = coh.observed_cp10k
    assert obs[crc.TARGET_STATE][pdcd1] > 3 * obs[crc.HEALTHY_STATE][pdcd1]  # PD-1 up in exhausted TILs
    for adata in coh.states.values():
        assert adata.var_names[0].startswith("ENSG")
        x = adata.X
        assert np.allclose(x.data, np.round(x.data))  # raw integer UMIs
        assert adata.obs["donor_id"].nunique() >= 20


class _NeighborBackend:
    """Target moves `target_lfc`; a neighbouring gene moves `neighbor_lfc`."""

    backend_mode = BackendMode.SURROGATE_SIMULATOR
    model_version = "fake"

    def __init__(self, target_lfc: float, neighbor_lfc: float) -> None:
        self.t, self.n = target_lfc, neighbor_lfc

    def score_variant_multimodal(self, parent_interval, variant, gene_id, cell_sketch, condition_name, reference_curies=None):  # type: ignore[no-untyped-def]
        lfc = self.t if cell_sketch is not None and "target" in condition_name.lower() else 0.0
        return VariantMultiModalScore(
            variant, condition_name, lfc, 0.001, 0.001, 0.001,
            neighbor_gene_lfcs={"NEIGHBOR_GENE": self.n},
        )


def _run_neighbor(target_lfc: float, neighbor_lfc: float):  # type: ignore[no-untyped-def]
    return evaluate_candidate_two_tier(
        backend=_NeighborBackend(target_lfc, neighbor_lfc),  # type: ignore[arg-type]
        parent_interval=GenomicInterval(chromosome="chr1", start=0, end=2000),
        candidate=_candidate(), gene_id="ENSG0", target_sketches=_sketches("target"),
        healthy_sketches_map={"healthy": _sketches("healthy")}, reference_curies=["UBERON:0000178", "CL:0000625"],
        thresholds=VerificationThresholds(desired_direction="down"),
    )


def test_neighbor_gene_collateral_fails_the_candidate() -> None:
    # The target drops 0.20 but a neighbouring gene in the same interval drops 0.15 (75%):
    # the edit is hitting a shared regulatory element, not a PDCD1-specific one.
    rep = _run_neighbor(-0.20, -0.15)
    assert not rep.all_checks_passed
    assert not rep.tier1_selectivity.neighbor_specificity_passed
    assert rep.tier1_selectivity.worst_neighbor_gene == "NEIGHBOR_GENE"
    assert any(r.startswith("TIER1_NEIGHBOR_GENE") for r in rep.failure_reasons)


def test_a_gene_specific_edit_passes_the_neighbor_check() -> None:
    rep = _run_neighbor(-0.20, -0.01)
    assert rep.tier1_selectivity.neighbor_specificity_passed
    assert rep.tier1_selectivity.neighbor_gene_lfcs == {"NEIGHBOR_GENE": -0.01}


def test_hazard_report_stamps_the_evidence_tier_of_its_reference() -> None:
    rep = _run_neighbor(-0.20, -0.01)
    t2 = rep.tier2_hazards
    assert "CD8" in t2.reference_biosample_name
    assert t2.evidence_tiers["methylation"] == "observed"
    assert t2.limiting_evidence_tier == "observed"
