"""Tests for content-addressed SQLite caching & deterministic seed reproducibility."""

from __future__ import annotations

from pathlib import Path
from gene_edit_check.budget import RunBudgetTracker
from gene_edit_check.cache import ContentAddressedCache
from gene_edit_check.ag_backend.simulator import SurrogateSimulatorBackend
from gene_edit_check.scenarios.registry import build_til_cohort_anndata


def test_content_addressed_cache_roundtrip(tmp_path: Path) -> None:
    cache = ContentAddressedCache(cache_dir=tmp_path / "cache", model_version="v2.0-test")
    key1 = cache.make_key("score_variant", "chr2:100-200|150:A>G", "cond_hash_1", "SC_EXPR")
    key2 = cache.make_key("score_variant", "chr2:100-200|150:A>G", "cond_hash_1", "SC_EXPR")
    assert key1 == key2
    assert cache.get(key1) is None

    payload = {"target_gene_lfc": -0.0412, "ok": True}
    cache.put(key1, "score_variant", payload)
    loaded = cache.get(key1)
    assert loaded == payload


def test_deterministic_seed_reproducibility() -> None:
    c1 = build_til_cohort_anndata(n_cells_per_state=40, seed=123)
    c2 = build_til_cohort_anndata(n_cells_per_state=40, seed=123)

    b1 = SurrogateSimulatorBackend(RunBudgetTracker(), seed=123)
    b2 = SurrogateSimulatorBackend(RunBudgetTracker(), seed=123)

    e1 = b1.encode_cells(c1["Exhausted_CD8_TIL"], "Exhausted_CD8_TIL")
    e2 = b2.encode_cells(c2["Exhausted_CD8_TIL"], "Exhausted_CD8_TIL")
    assert (e1 == e2).all()
