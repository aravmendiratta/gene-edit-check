"""Off-tissue safety panel (`verifier/off_tissue.py`).

Tier 1 asks whether an edit is selective *between cell states* (exhausted tumour T cells vs
normal-colon T cells). That is the scientific question, but it is not the whole safety question.

A base editor is delivered to cells, not to a cell state. Whatever vector carries it will
reach tissues nobody intended. So for a finalist edit we also ask: under reference conditions
for tissues this therapy has no business touching -- liver, brain, heart, lung, kidney -- does
the model predict the locus moves at all?

This is a Tier 2-style check: it uses `ModelCondition.from_curies` reference conditions, never
single-cell sketches, so it must not be described as single-cell resolution. It costs one API
call per tissue, which is why it runs for finalists only.
"""

from __future__ import annotations

from gene_edit_check.ag_backend.protocol import VerifierBackend
from gene_edit_check.conditions import describe_condition, evidence_for
from gene_edit_check.models import (
    CandidateGuideEdit,
    GenomicInterval,
    OffTissueSafetyPanel,
)

# Non-immune reference conditions spanning the organs a systemic delivery vector reaches first.
# Chosen for breadth of lineage (endoderm/ectoderm/mesoderm) rather than for a specific hypothesis.
DEFAULT_OFF_TISSUE_CURIES: tuple[tuple[str, ...], ...] = (
    ("CL:0000182",),       # hepatocyte -- first-pass organ for most systemic vectors
    ("UBERON:0000955",),   # brain
    ("UBERON:0000948",),   # heart
    ("UBERON:0002048",),   # lung
    ("UBERON:0002113",),   # kidney
    ("UBERON:0002106",),   # spleen -- lymphoid, where off-target immune effects would show
)


def run_off_tissue_panel(
    backend: VerifierBackend,
    parent_interval: GenomicInterval,
    candidate: CandidateGuideEdit,
    gene_id: str,
    on_target_abs_lfc: float,
    tissue_curies: tuple[tuple[str, ...], ...] = DEFAULT_OFF_TISSUE_CURIES,
    max_fraction_of_on_target: float = 0.50,
) -> OffTissueSafetyPanel:
    """Score `candidate` under each non-target tissue reference condition.

    The bar is relative, not absolute: an off-tissue shift must stay below
    `max_fraction_of_on_target` of the on-target effect. A tiny absolute shift in liver is
    uninteresting if the on-target effect is large, and alarming if it is not.
    """
    tissues: dict[str, float] = {}
    tiers: dict[str, str] = {}

    for curies in tissue_curies:
        name = describe_condition("_".join(curies))
        score = backend.score_variant_multimodal(
            parent_interval=parent_interval,
            variant=candidate.target_edit,
            gene_id=gene_id,
            cell_sketch=None,
            condition_name=f"OffTissue::{name}",
            reference_curies=list(curies),
        )
        tissues[name] = round(abs(score.target_gene_lfc), 6)
        tiers[name] = evidence_for("_".join(curies), "rna_seq").label

    threshold = max_fraction_of_on_target * abs(on_target_abs_lfc)
    worst_tissue = max(tissues, key=lambda t: tissues[t]) if tissues else ""
    max_abs = tissues.get(worst_tissue, 0.0)

    return OffTissueSafetyPanel(
        tissues=tissues,
        evidence_tiers=tiers,
        max_abs_lfc=max_abs,
        worst_tissue=worst_tissue,
        threshold=round(threshold, 6),
        passed=bool(max_abs <= threshold),
    )
