"""Tests for evidence-tiered condition selection, neighbour-gene specificity and the
off-tissue safety panel."""

from __future__ import annotations

import pytest

from crispr_reflexion.ag_backend.protocol import VariantMultiModalScore
from crispr_reflexion.conditions import (
    EvidenceTier,
    best_evidence_condition,
    describe_condition,
    evidence_for,
    evidence_profile,
    load_conditions_table,
    rank_by_evidence,
    search_conditions,
    weakest_hazard_tier,
)
from crispr_reflexion.models import BackendMode, BaseEditorType, GenomicInterval
from crispr_reflexion.verifier.off_tissue import run_off_tissue_panel

BLOOD_CD8 = "UBERON:0000178_CL:0000625"
COLON_CD8 = "UBERON:0001155_CL:0000625"


def test_conditions_table_ships_with_the_package() -> None:
    df = load_conditions_table()
    assert len(df) > 2000
    assert {"condition", "biosample_name", "methylation", "splice_junctions"} <= set(df.columns)


def test_evidence_tiers_are_ordered_so_max_picks_the_strongest() -> None:
    assert EvidenceTier.OBSERVED > EvidenceTier.IMPUTED_SUPPORTED > EvidenceTier.IMPUTED_EXPLORATORY
    assert EvidenceTier.parse("observed") is EvidenceTier.OBSERVED
    assert EvidenceTier.parse("nonsense") is EvidenceTier.UNKNOWN
    assert EvidenceTier.OBSERVED.is_trustworthy
    assert not EvidenceTier.IMPUTED_EXPLORATORY.is_trustworthy


def test_blood_cd8_is_the_best_evidenced_cd8_hazard_reference() -> None:
    """The documented reason the project uses blood CD8 rather than the tissue-matched colon
    condition: it is the only CD8 T-cell condition with observed methylation AND splicing."""
    assert evidence_for(BLOOD_CD8, "methylation") is EvidenceTier.OBSERVED
    assert evidence_for(BLOOD_CD8, "splice_junctions") is EvidenceTier.OBSERVED
    assert evidence_for(COLON_CD8, "methylation") is EvidenceTier.IMPUTED_EXPLORATORY
    assert best_evidence_condition([COLON_CD8, BLOOD_CD8]) == BLOOD_CD8
    assert weakest_hazard_tier(BLOOD_CD8) is EvidenceTier.OBSERVED
    assert weakest_hazard_tier(COLON_CD8) is EvidenceTier.IMPUTED_EXPLORATORY


def test_describe_and_profile_are_report_ready() -> None:
    assert "CD8" in describe_condition(BLOOD_CD8)
    assert describe_condition("CL:9999999") == "CL:9999999"  # unknown falls back to the CURIE
    prof = evidence_profile(BLOOD_CD8)
    assert prof["methylation"] == "observed"


def test_rank_by_evidence_orders_candidates() -> None:
    ranked = rank_by_evidence([COLON_CD8, BLOOD_CD8])
    assert ranked[0][0] == BLOOD_CD8
    assert ranked[0][1] > ranked[1][1]


def test_search_sorts_by_evidence_for_the_requested_modality() -> None:
    hits = search_conditions("CD8-positive", modality="methylation")
    assert len(hits) > 10
    assert hits.iloc[0]["methylation"] == "observed"


class _TissueBackend:
    """Reports a large shift in one named tissue and nothing elsewhere."""

    backend_mode = BackendMode.SURROGATE_SIMULATOR
    model_version = "fake"

    def __init__(self, loud_tissue_curie: str, loud_lfc: float) -> None:
        self.loud, self.lfc = loud_tissue_curie, loud_lfc
        self.calls: list[str] = []

    def score_variant_multimodal(self, parent_interval, variant, gene_id, cell_sketch, condition_name, reference_curies=None):  # type: ignore[no-untyped-def]
        curie = "_".join(reference_curies or [])
        self.calls.append(curie)
        lfc = self.lfc if curie == self.loud else 0.001
        return VariantMultiModalScore(variant, condition_name, lfc, 0.0, 0.0, 0.0)


def _cand():  # type: ignore[no-untyped-def]
    from crispr_reflexion.models import CandidateGuideEdit, SingleNucleotideVariant

    v = SingleNucleotideVariant(chromosome="chr2", position=500, reference="A", alternate="G", protospacer_index=5)
    return CandidateGuideEdit(
        candidate_id="c", gene_symbol="PDCD1", editor=BaseEditorType.ABE8E, strand="+",
        protospacer_20nt="A" * 20, pam_sequence="AGG", genomic_start=490, genomic_end=513, target_edit=v,
    )


def _panel(loud_curie: str, loud_lfc: float, on_target: float):  # type: ignore[no-untyped-def]
    return run_off_tissue_panel(
        backend=_TissueBackend(loud_curie, loud_lfc),  # type: ignore[arg-type]
        parent_interval=GenomicInterval(chromosome="chr2", start=0, end=1000),
        candidate=_cand(), gene_id="ENSG00000188389", on_target_abs_lfc=on_target,
    )


def test_off_tissue_panel_passes_a_tissue_silent_edit() -> None:
    p = _panel("CL:0000182", 0.001, on_target=0.20)
    assert p.passed
    assert len(p.tissues) == 6
    assert "hepatocyte" in p.tissues


def test_off_tissue_panel_flags_a_liver_effect() -> None:
    p = _panel("CL:0000182", 0.15, on_target=0.20)  # 75% of the on-target effect, in liver
    assert not p.passed
    assert p.worst_tissue == "hepatocyte"
    assert p.max_abs_lfc == pytest.approx(0.15)


def test_off_tissue_bar_is_relative_to_the_on_target_effect() -> None:
    # The same absolute liver shift is acceptable beside a large on-target effect and not
    # beside a small one.
    assert _panel("CL:0000182", 0.05, on_target=0.40).passed
    assert not _panel("CL:0000182", 0.05, on_target=0.06).passed


def test_off_tissue_panel_records_evidence_tiers() -> None:
    p = _panel("CL:0000182", 0.001, on_target=0.2)
    assert p.evidence_tiers["hepatocyte"] == "observed"
