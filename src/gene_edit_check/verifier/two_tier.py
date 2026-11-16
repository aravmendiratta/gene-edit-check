"""Two-Tier Multi-Cell-State & Sequence-Hazard Verifier (`verifier/two_tier.py`).

Enforces strict separation between:
- Tier 1 (Single-Cell Conditioned Selectivity):
  Evaluates candidate base edits on Proposal `Sketch A` AND Held-Out `Sketch B`
  (disjoint donors/cells) for target efficacy, cross-sketch generalization, and
  healthy cell-state sparing.
- Tier 2 (Sequence-Intrinsic & Reference Ontology Hazards):
  Evaluates `SpliceJunctionScorer()`, `MethylationScorer(width=1201)`, and
  bystander edits in the deaminase window `[4..8]`.
"""

from __future__ import annotations

from typing import Literal

import numpy as np

from gene_edit_check.ag_backend.protocol import VerifierBackend
from gene_edit_check.conditions import describe_condition, evidence_profile, weakest_hazard_tier
from gene_edit_check.models import (
    CandidateGuideEdit,
    GenomicInterval,
    Tier1SelectivityResult,
    Tier2HazardResult,
    TwoTierVerificationReport,
)
from gene_edit_check.sc_sketch.coreset import DisjointCellSketches


class VerificationThresholds:
    def __init__(
        self,
        min_target_abs_lfc: float = 0.0025,
        max_healthy_abs_lfc_ratio: float = 0.35,
        max_splice_junction_delta: float = 0.15,
        max_cpg_methylation_shift: float = 0.12,
        require_bystander_clean: bool = True,
        desired_direction: Literal["down", "up"] = "down",
        null_abs_lfc_cutoff: float = 0.0,
        max_neighbor_lfc_ratio: float = 0.50,
    ) -> None:
        self.desired_direction = desired_direction
        self.null_abs_lfc_cutoff = null_abs_lfc_cutoff
        self.min_target_abs_lfc = min_target_abs_lfc
        self.max_healthy_abs_lfc_ratio = max_healthy_abs_lfc_ratio
        self.max_splice_junction_delta = max_splice_junction_delta
        self.max_cpg_methylation_shift = max_cpg_methylation_shift
        self.require_bystander_clean = require_bystander_clean
        # No neighbouring gene in the window may move more than half as much as the target.
        self.max_neighbor_lfc_ratio = max_neighbor_lfc_ratio


def evaluate_candidate_two_tier(
    backend: VerifierBackend,
    parent_interval: GenomicInterval,
    candidate: CandidateGuideEdit,
    gene_id: str,
    target_sketches: DisjointCellSketches,
    healthy_sketches_map: dict[str, DisjointCellSketches],
    reference_curies: list[str],
    thresholds: VerificationThresholds | None = None,
    tissue_matched_curies: list[str] | None = None,
) -> TwoTierVerificationReport:
    cfg = thresholds or VerificationThresholds()
    failures: list[str] = []

    # 1. Tier 1: Evaluate on Target Sketch A (proposal) and Held-Out Target Sketch B (verification)
    score_target_a = backend.score_variant_multimodal(
        parent_interval=parent_interval,
        variant=candidate.target_edit,
        gene_id=gene_id,
        cell_sketch=target_sketches.sketch_a,
        condition_name=f"{target_sketches.cell_state_name}_SketchA",
    )
    score_target_b = backend.score_variant_multimodal(
        parent_interval=parent_interval,
        variant=candidate.target_edit,
        gene_id=gene_id,
        cell_sketch=target_sketches.sketch_b,
        condition_name=f"{target_sketches.cell_state_name}_SketchB",
    )

    lfc_a = score_target_a.target_gene_lfc
    lfc_b = score_target_b.target_gene_lfc
    bootstrap_se = float(abs(lfc_a - lfc_b) / np.sqrt(2.0))

    cutoff = max(cfg.min_target_abs_lfc, cfg.null_abs_lfc_cutoff)
    direction_ok = (lfc_b < 0) if cfg.desired_direction == "down" else (lfc_b > 0)
    magnitude_ok = abs(lfc_b) >= cutoff
    target_efficacy_ok = direction_ok and magnitude_ok
    if not direction_ok:
        failures.append(
            f"TIER1_DIRECTION: Held-out Sketch B LFC={lfc_b:+.5f} moves the target gene the wrong way "
            f"(goal: {cfg.desired_direction}-regulate)"
        )
    elif not magnitude_ok:
        failures.append(
            f"TIER1_EFFICACY: Held-out Sketch B |LFC|={abs(lfc_b):.5f} < cutoff {cutoff:.5f} "
            f"(null-calibrated; must beat random PAM-valid edits)"
        )

    # Check sign consistency and magnitude retention on Held-Out Sketch B
    same_sign = (lfc_a * lfc_b) > 0 or (abs(lfc_a) < 1e-6 and abs(lfc_b) < 1e-6)
    cross_sketch_ok = bool(same_sign and abs(lfc_b) >= 0.65 * abs(lfc_a))
    if not cross_sketch_ok:
        failures.append(
            f"TIER1_GOODHART_GUARD: Failed Sketch A ({lfc_a:+.5f}) -> Held-Out Sketch B ({lfc_b:+.5f}) transfer"
        )

    # Evaluate on Healthy Control Held-Out Sketches
    healthy_deltas: dict[str, float] = {}
    healthy_selectivity_ok = True
    for h_name, h_sketches in healthy_sketches_map.items():
        h_score = backend.score_variant_multimodal(
            parent_interval=parent_interval,
            variant=candidate.target_edit,
            gene_id=gene_id,
            cell_sketch=h_sketches.sketch_b,
            condition_name=f"{h_name}_SketchB",
        )
        healthy_deltas[h_name] = h_score.target_gene_lfc
        if abs(lfc_b) > 1e-7:
            ratio = abs(h_score.target_gene_lfc) / abs(lfc_b)
        else:
            ratio = 999.0
        if ratio > cfg.max_healthy_abs_lfc_ratio:
            healthy_selectivity_ok = False
            failures.append(
                f"TIER1_SELECTIVITY: Off-target shift in healthy '{h_name}' "
                f"({h_score.target_gene_lfc:+.5f}, ratio={ratio:.2f} > {cfg.max_healthy_abs_lfc_ratio:.2f})"
            )

    # Neighbour-gene specificity. The expression scorer returned every gene in the interval
    # alongside the target, so this check is free: it only requires not discarding those rows.
    neighbor_lfcs = score_target_b.neighbor_gene_lfcs
    worst = score_target_b.worst_neighbor
    neighbor_ok = True
    worst_gene, worst_lfc = (worst if worst else (None, 0.0))
    if worst is not None and abs(lfc_b) > 1e-7:
        n_ratio = abs(worst_lfc) / abs(lfc_b)
        if n_ratio > cfg.max_neighbor_lfc_ratio:
            neighbor_ok = False
            failures.append(
                f"TIER1_NEIGHBOR_GENE: '{worst_gene}' in the same interval shifts {worst_lfc:+.5f} "
                f"(ratio={n_ratio:.2f} > {cfg.max_neighbor_lfc_ratio:.2f} of the target's "
                f"{lfc_b:+.5f}) -- the edit is not specific to {gene_id}"
            )

    tier1 = Tier1SelectivityResult(
        target_condition_name=target_sketches.cell_state_name,
        target_delta_sketch_a=round(lfc_a, 6),
        target_delta_sketch_b=round(lfc_b, 6),
        target_sketch_bootstrap_se=round(bootstrap_se, 6),
        healthy_control_deltas={k: round(v, 6) for k, v in healthy_deltas.items()},
        target_efficacy_passed=target_efficacy_ok,
        healthy_selectivity_passed=healthy_selectivity_ok,
        cross_sketch_generalization_passed=cross_sketch_ok,
        neighbor_gene_lfcs={k: round(v, 6) for k, v in neighbor_lfcs.items()},
        worst_neighbor_gene=worst_gene,
        worst_neighbor_lfc=round(worst_lfc, 6),
        neighbor_specificity_passed=neighbor_ok,
    )

    # 2. Tier 2: Sequence-Intrinsic & Ontology Reference Hazards (Splicing, CpG Methylation, Bystanders)
    ref_hazard_score = backend.score_variant_multimodal(
        parent_interval=parent_interval,
        variant=candidate.target_edit,
        gene_id=gene_id,
        cell_sketch=None,
        condition_name="Reference_Ontology_Check",
        reference_curies=reference_curies,
    )
    max_splice = max(
        score_target_b.max_splice_junction_delta,
        ref_hazard_score.max_splice_junction_delta,
    )
    mean_cpg = max(
        score_target_b.cpg_methylation_mean_abs_shift,
        ref_hazard_score.cpg_methylation_mean_abs_shift,
    )

    splice_ok = max_splice <= cfg.max_splice_junction_delta
    if not splice_ok:
        failures.append(
            f"TIER2_SPLICE_HAZARD: SpliceJunctionScorer max delta={max_splice:.4f} > {cfg.max_splice_junction_delta:.4f}"
        )

    meth_ok = mean_cpg <= cfg.max_cpg_methylation_shift
    if not meth_ok:
        failures.append(
            f"TIER2_METHYLATION_HAZARD: MethylationScorer(1201) shift={mean_cpg:.4f} > {cfg.max_cpg_methylation_shift:.4f}"
        )

    # Stamp how well the model actually knows this reference condition for the two modalities
    # we are trusting. A reader can then discount an `imputed_exploratory` verdict appropriately.
    primary_curie = "_".join(reference_curies)
    evidence = evidence_profile(primary_curie, ("methylation", "splice_junctions"))
    limiting = weakest_hazard_tier(primary_curie)

    # Tissue-matched cross-check: anatomically right, evidentially weaker. Reported side by side
    # so the two never get averaged into one falsely confident number.
    tm_splice: float | None = None
    tm_cpg: float | None = None
    tm_agrees: bool | None = None
    tm_name = ""
    if tissue_matched_curies:
        tm_name = describe_condition("_".join(tissue_matched_curies))
        tm_score = backend.score_variant_multimodal(
            parent_interval=parent_interval,
            variant=candidate.target_edit,
            gene_id=gene_id,
            cell_sketch=None,
            condition_name="Tissue_Matched_Reference_Check",
            reference_curies=tissue_matched_curies,
        )
        tm_splice = round(tm_score.max_splice_junction_delta, 6)
        tm_cpg = round(tm_score.cpg_methylation_mean_abs_shift, 6)
        tm_agrees = bool(
            (tm_score.max_splice_junction_delta <= cfg.max_splice_junction_delta) == splice_ok
            and (tm_score.cpg_methylation_mean_abs_shift <= cfg.max_cpg_methylation_shift) == meth_ok
        )

    bystander_ok = candidate.is_bystander_clean or (not cfg.require_bystander_clean)
    if not bystander_ok:
        n_bys = len(candidate.bystander_edits_in_window)
        failures.append(
            f"TIER2_BYSTANDER_HAZARD: Protospacer contains {n_bys} bystander editable base(s) in window [4..8]"
        )

    tier2 = Tier2HazardResult(
        reference_curies=reference_curies,
        max_splice_junction_delta=round(max_splice, 6),
        cpg_methylation_mean_abs_shift=round(mean_cpg, 6),
        splicing_hazard_passed=splice_ok,
        methylation_hazard_passed=meth_ok,
        bystander_clean_passed=bystander_ok,
        reference_biosample_name=describe_condition(primary_curie),
        evidence_tiers=evidence,
        limiting_evidence_tier=limiting.label,
        tissue_matched_curies=list(tissue_matched_curies or []),
        tissue_matched_biosample_name=tm_name,
        tissue_matched_splice_delta=tm_splice,
        tissue_matched_cpg_shift=tm_cpg,
        tissue_matched_agrees=tm_agrees,
    )

    max_healthy = max((abs(v) for v in healthy_deltas.values()), default=0.0)
    bystander_penalty = 0.015 * len(candidate.bystander_edits_in_window)
    desired_effect = -lfc_b if cfg.desired_direction == "down" else lfc_b
    composite_score = float(
        desired_effect
        - 0.75 * max_healthy
        - 0.40 * abs(worst_lfc)
        - 0.08 * max_splice
        - 0.08 * mean_cpg
        - bystander_penalty
    )
    all_passed = (
        target_efficacy_ok
        and healthy_selectivity_ok
        and neighbor_ok
        and cross_sketch_ok
        and splice_ok
        and meth_ok
        and bystander_ok
    )

    return TwoTierVerificationReport(
        candidate=candidate,
        tier1_selectivity=tier1,
        tier2_hazards=tier2,
        composite_objective_score=round(composite_score, 6),
        all_checks_passed=all_passed,
        failure_reasons=failures,
    )
