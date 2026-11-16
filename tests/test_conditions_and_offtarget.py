"""Tests for evidence-tiered condition selection, neighbour-gene specificity and the
off-tissue safety panel."""

from __future__ import annotations

import pytest

from gene_edit_check.ag_backend.protocol import VariantMultiModalScore
from gene_edit_check.conditions import (
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
from gene_edit_check.models import BackendMode, BaseEditorType, GenomicInterval
from gene_edit_check.verifier.off_tissue import run_off_tissue_panel

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
    from gene_edit_check.models import CandidateGuideEdit, SingleNucleotideVariant

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


def test_mcp_tool_functions_work_without_the_mcp_package() -> None:
    """The tool bodies are plain functions, so they are testable and usable without MCP installed."""
    from gene_edit_check.mcp_server.server import list_scenarios_tool, run_optimization_tool, search_conditions_tool

    ids = {s["scenario_id"] for s in list_scenarios_tool()}
    assert "crc-pd1-exhaustion" in ids and "pd1-exhaustion" in ids
    hits = search_conditions_tool("blood CD8", modality="methylation", limit=3)
    assert hits and hits[0]["methylation"] == "observed"
    with pytest.raises(ValueError, match="Unknown scenario"):
        run_optimization_tool(scenario_id="nope")


def test_installed_package_never_writes_into_site_packages(monkeypatch, tmp_path) -> None:  # type: ignore[no-untyped-def]
    """Data and cache go to a per-user directory (or GENE_EDIT_CHECK_HOME), never the install tree."""
    from gene_edit_check import paths

    monkeypatch.setenv("GENE_EDIT_CHECK_HOME", str(tmp_path))
    assert paths.home() == tmp_path
    assert paths.data_dir() == tmp_path / "data"
    assert tmp_path in paths.api_cache_dir().parents


def test_real_data_run_without_a_key_refuses_instead_of_faking(monkeypatch) -> None:  # type: ignore[no-untyped-def]
    import typer
    from gene_edit_check.budget import RunBudgetTracker
    from gene_edit_check.cli import _resolve_backend

    monkeypatch.delenv("ALPHAGENOME_API_KEY", raising=False)
    with pytest.raises(typer.Exit):
        _resolve_backend(use_simulator=False, budget=RunBudgetTracker(), require_live=True)
    # An explicit --simulator and a synthetic scenario are still allowed (with the red warning).
    assert _resolve_backend(use_simulator=True, budget=RunBudgetTracker(), require_live=True).backend_mode.value == "SURROGATE_SIMULATOR"
    assert _resolve_backend(use_simulator=False, budget=RunBudgetTracker(), require_live=False).backend_mode.value == "SURROGATE_SIMULATOR"


def test_doctor_runs_without_a_key(monkeypatch) -> None:  # type: ignore[no-untyped-def]
    from typer.testing import CliRunner
    from gene_edit_check.cli import app

    monkeypatch.delenv("ALPHAGENOME_API_KEY", raising=False)
    result = CliRunner().invoke(app, ["doctor"])
    assert result.exit_code == 0
    assert "setup check" in result.output
