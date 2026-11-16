"""Provenance tracking & explicit backend badging (`provenance.py`).

Ensures every CLI output, JSON artifact, and figure clearly records whether results
came from `LIVE_ALPHAGENOME_V2` or `SURROGATE_SIMULATOR`, along with seeds, model version,
and call counts.
"""

from __future__ import annotations

from gene_edit_check.budget import RunBudgetTracker
from gene_edit_check.models import BackendMode, ProvenanceRecord, RESEARCH_USE_ONLY_NOTICE


def build_provenance_record(
    backend_mode: BackendMode,
    model_version: str,
    random_seed: int,
    budget: RunBudgetTracker,
) -> ProvenanceRecord:
    return ProvenanceRecord(
        backend_mode=backend_mode,
        model_version=model_version,
        random_seed=random_seed,
        total_api_calls=budget.total_api_calls,
        cache_hits=budget.cache_hits,
        total_ism_variants_scored=budget.ism_variants_scored,
        gradient_calls_used=budget.gradient_calls_used,
        research_use_only_notice=RESEARCH_USE_ONLY_NOTICE,
    )
