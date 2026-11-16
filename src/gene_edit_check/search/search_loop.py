"""Gradient-guided search with scoped in-silico mutagenesis (`search/search_loop.py`).

The search splits the work so the expensive model call is spent only where it is needed:

1. Input gradients on the target-state and healthy-state cell sketches show where in the locus the
   model is sensitive, and recover the reference sequence.
2. Every PAM-valid ABE8e and CBE edit in the scan region (the promoter by default) is listed and
   ranked by saliency.
3. In-silico mutagenesis runs on short windows around the top candidates to measure the direction of
   effect, which the gradient sign alone does not give reliably.
4. Candidates go through the two-tier verifier (held-out Sketch B for expression, reference
   conditions for splice and methylation hazards). When one fails, the failure reason is recorded and
   the remaining candidates are re-ranked with that constraint.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from gene_edit_check.ag_backend.protocol import VerifierBackend
from gene_edit_check.budget import RunBudgetTracker
from gene_edit_check.crispr.constraints import enumerate_valid_base_edits
from gene_edit_check.crispr.motifs import annotate_motifs_in_window
from gene_edit_check.models import (
    CandidateGuideEdit,
    EfficacyNullCalibration,
    GenomicInterval,
    IsmPrescreenSummary,
    OptimizationRunResult,
    SearchStepTrace,
    SketchQualityMetrics,
    TwoTierVerificationReport,
)
from gene_edit_check.provenance import build_provenance_record
from gene_edit_check.verifier.off_tissue import run_off_tissue_panel
from gene_edit_check.sc_sketch.coreset import DisjointCellSketches
from gene_edit_check.verifier.two_tier import (
    VerificationThresholds,
    evaluate_candidate_two_tier,
)


def _variant_key_of(c: CandidateGuideEdit) -> tuple[int, str, str]:
    return (c.target_edit.position, c.target_edit.reference, c.target_edit.alternate)


def _is_desired_direction(lfc: float, desired_direction: str) -> bool:
    return lfc < 0 if desired_direction == "down" else lfc > 0


def run_gradient_ism_search(
    scenario_id: str,
    gene_id: str,
    gene_symbol: str,
    backend: VerifierBackend,
    budget: RunBudgetTracker,
    target_sketches: DisjointCellSketches,
    healthy_sketches_map: dict[str, DisjointCellSketches],
    sketch_metrics: list[SketchQualityMetrics],
    reference_curies: list[str],
    thresholds: VerificationThresholds | None = None,
    max_iterations: int = 5,
    random_seed: int = 42,
    desired_direction: str = "down",
    null_n: int = 12,
    ism_prescreen_windows: int = 10,
    ism_prescreen_width_bp: int = 8,
    tissue_matched_curies: list[str] | None = None,
    run_off_tissue_panel_on_finalist: bool = True,
    data_provenance: dict[str, Any] | None = None,
    scan_regions: list[tuple[int, int]] | None = None,
) -> OptimizationRunResult:
    # Step 1: Gradients Propose on Sketch A (Target vs. primary Healthy control)
    grad_target = backend.compute_interval_gradients(
        gene_id=gene_id,
        gene_symbol=gene_symbol,
        cell_sketch=target_sketches.sketch_a,
        condition_name=target_sketches.cell_state_name,
    )
    primary_healthy_name = next(iter(healthy_sketches_map.keys()))
    primary_healthy_sketches = healthy_sketches_map[primary_healthy_name]
    grad_healthy = backend.compute_interval_gradients(
        gene_id=gene_id,
        gene_symbol=gene_symbol,
        cell_sketch=primary_healthy_sketches.sketch_a,
        condition_name=primary_healthy_name,
    )

    parent_iv = grad_target.interval
    ref_seq = grad_target.reference_sequence
    abs_t = np.abs(grad_target.signed_gradients_1bp)
    abs_h = np.abs(grad_healthy.signed_gradients_1bp)
    signed_t = grad_target.signed_gradients_1bp
    rng = np.random.default_rng(random_seed)

    # Step 2: Scan for PAM-valid ABE8e/CBE candidates. By default that is +/- 2,500 bp around the TSS
    # (the promoter). `scan_regions` takes TSS-relative (start, end) offsets instead, so the search
    # can be pointed at a distal element such as a published enhancer.
    mid = len(ref_seq) // 2
    regions = scan_regions or [(-2500, 2500)]
    raw_candidates = []
    seen_ids: set[str] = set()
    for lo_off, hi_off in regions:
        scan_start_offset = max(0, mid + lo_off)
        scan_end_offset = min(len(ref_seq), mid + hi_off)
        if scan_end_offset - scan_start_offset < 24:
            raise ValueError(f"Scan region {lo_off}:{hi_off} falls outside the {len(ref_seq)}-bp model interval.")
        for cand in enumerate_valid_base_edits(
            chromosome=parent_iv.chromosome,
            window_start_genomic=parent_iv.start + scan_start_offset,
            window_sequence=ref_seq[scan_start_offset:scan_end_offset],
            gene_symbol=gene_symbol,
            pam_rule="NGG",
        ):
            if cand.candidate_id not in seen_ids:
                seen_ids.add(cand.candidate_id)
                raw_candidates.append(cand)
    data_provenance = {**(data_provenance or {}), "scan_regions_tss_relative_bp": [list(r) for r in regions]}

    # Annotate each candidate with Target Saliency, Healthy Saliency, and overlapping TF motifs
    candidate_healthy_grad: dict[str, float] = {}
    candidate_signed_grad: dict[str, float] = {}
    for cand in raw_candidates:
        rel_idx = cand.target_edit.position - 1 - parent_iv.start
        lo = max(0, rel_idx - 2)
        hi = min(len(abs_t), rel_idx + 3)
        t_sal = float(np.mean(abs_t[lo:hi]))
        h_sal = float(np.mean(abs_h[lo:hi]))
        cand.gradient_saliency_score = round(t_sal, 6)
        candidate_healthy_grad[cand.candidate_id] = h_sal
        candidate_signed_grad[cand.candidate_id] = float(np.mean(signed_t[lo:hi]))

        motifs = annotate_motifs_in_window(
            cand.protospacer_20nt, (cand.target_edit.protospacer_index or 5) - 1
        )
        if motifs:
            cand.disrupted_motif = motifs[0]

    # Sort initial pool by raw target gradient saliency: the gradient's job is to say *where* to
    # look, which it does well. Which way expression moves is a separate question, answered next.
    raw_candidates.sort(key=lambda c: c.gradient_saliency_score, reverse=True)
    shortlist = raw_candidates[:60]

    # Step 2.5: ISM PRE-SCREEN. Measure direction instead of guessing it.
    #
    # The gradient sign is a first-order heuristic ('more of this base -> more expression, so
    # replacing it lowers expression'). Relying on it outright wasted whole verification rounds
    # on edits that turned out to *raise* PD-1. Instead we spend a few calls up front:
    # `score_ism_variants` over compact windows returns the true LFC of every substitution in
    # each window, so one call direction-screens `3 * width` variants at once. The verifier then
    # only ever sees candidates whose direction is already measured to be correct.
    measured_lfc: dict[tuple[int, str, str], float] = {}
    prescreen_summary: IsmPrescreenSummary | None = None
    if ism_prescreen_windows > 0:
        centers: list[int] = []
        half = ism_prescreen_width_bp // 2
        for cand in shortlist:
            pos0 = cand.target_edit.position - 1
            if any(abs(pos0 - c) < ism_prescreen_width_bp for c in centers):
                continue  # keep windows non-overlapping so each call covers new sequence
            centers.append(pos0)
            if len(centers) >= ism_prescreen_windows:
                break

        for center in centers:
            try:
                ism_map = backend.score_ism_window(
                    parent_interval=parent_iv,
                    ism_window=GenomicInterval(
                        chromosome=parent_iv.chromosome,
                        start=center - half,
                        end=center - half + ism_prescreen_width_bp,
                    ),
                    gene_id=gene_id,
                    cell_sketch=target_sketches.sketch_a,
                    condition_name=target_sketches.cell_state_name,
                )
            except Exception:
                break  # out of ISM budget: screen with what we have rather than failing the run
            measured_lfc.update(ism_map)

        # Audit the heuristic we just replaced: how often did the gradient sign get it right?
        agree = total = 0
        for (pos, _ref, _alt), lfc in measured_lfc.items():
            rel = pos - 1 - parent_iv.start
            if not (0 <= rel < len(signed_t)) or abs(lfc) < 1e-9:
                continue
            total += 1
            predicted_down = signed_t[rel] > 0
            if predicted_down == (lfc < 0):
                agree += 1

        with_meas = [c for c in shortlist if _variant_key_of(c) in measured_lfc]
        n_right = sum(
            1
            for c in with_meas
            if _is_desired_direction(measured_lfc[_variant_key_of(c)], desired_direction)
        )
        best = min(measured_lfc.values(), default=0.0) if desired_direction == "down" else max(
            measured_lfc.values(), default=0.0
        )
        prescreen_summary = IsmPrescreenSummary(
            n_windows=len(centers),
            window_width_bp=ism_prescreen_width_bp,
            n_variants_measured=len(measured_lfc),
            n_candidates_with_measurement=len(with_meas),
            n_correct_direction=n_right,
            n_wrong_direction_filtered=len(with_meas) - n_right,
            gradient_sign_agreement=round(agree / total, 4) if total else 0.0,
            best_measured_lfc=round(float(best), 6),
        )

    if thresholds is None:
        thresholds = VerificationThresholds()
    thresholds.desired_direction = desired_direction  # type: ignore[assignment]

    # Empirical null for the efficacy test: score random PAM-valid edits at this locus on the same
    # held-out target sketch the candidates are judged on. Thresholds are fixed BEFORE any candidate
    # is verified (no calibration on the candidate under test).
    efficacy_null: EfficacyNullCalibration | None = None
    if null_n > 0 and len(raw_candidates) > null_n:
        null_pool = [raw_candidates[i] for i in rng.choice(len(raw_candidates), size=null_n, replace=False)]
        null_vals: list[float] = []
        null_meth: list[float] = []
        null_splice: list[float] = []
        for nc in null_pool:
            sc = backend.score_variant_multimodal(
                parent_interval=parent_iv,
                variant=nc.target_edit,
                gene_id=gene_id,
                cell_sketch=target_sketches.sketch_b,
                condition_name=f"{target_sketches.cell_state_name}_SketchB_null",
            )
            null_vals.append(abs(sc.target_gene_lfc))
            null_meth.append(sc.cpg_methylation_mean_abs_shift)
            null_splice.append(sc.max_splice_junction_delta)
        thresholds.null_abs_lfc_cutoff = float(max(null_vals))
        efficacy_null = EfficacyNullCalibration(
            n_random_edits=len(null_vals),
            abs_lfc_values=[round(v, 6) for v in null_vals],
            abs_lfc_max=round(float(max(null_vals)), 6),
            abs_lfc_mean=round(float(np.mean(null_vals)), 6),
            cutoff_used=round(max(thresholds.min_target_abs_lfc, thresholds.null_abs_lfc_cutoff), 6),
            empirical_p_floor=round(1.0 / (len(null_vals) + 1), 4),
            cpg_methylation_shift_values=[round(v, 6) for v in null_meth],
            splice_delta_values=[round(v, 6) for v in null_splice],
        )

    # Constraint flags updated by the failure analysis as failures occur
    require_bystander_free = False
    healthy_gradient_penalty_weight = 0.0
    wrong_sign_multiplier = 0.25
    evaluated_ids: set[str] = set()
    # Same genomic edit via a different guide gives an identical Tier 1/2 score; only re-try it if the
    # new guide has fewer bystander edits (the one thing that is guide-specific).
    tried_variant_min_bystanders: dict[tuple[int, str, str], int] = {}

    _variant_key = _variant_key_of
    traces: list[SearchStepTrace] = []
    winning_report: TwoTierVerificationReport | None = None

    for step_idx in range(1, max_iterations + 1):
        # Rank unevaluated candidates using current constraints
        def _rank_key(c: CandidateGuideEdit) -> float:
            if require_bystander_free and not c.is_bystander_clean:
                return -999.0
            h_g = candidate_healthy_grad.get(c.candidate_id, 0.0)
            penalty = healthy_gradient_penalty_weight * h_g
            vk = _variant_key_of(c)

            if vk in measured_lfc:
                # Measured beats inferred. The pre-screen already ran this exact substitution
                # through the model, so rank on the real effect size and put anything moving the
                # wrong way strictly below every unmeasured candidate.
                lfc = measured_lfc[vk]
                if not _is_desired_direction(lfc, desired_direction):
                    return -500.0 + abs(lfc)
                return 1000.0 + abs(lfc) - penalty

            # Unmeasured (the pre-screen covers compact windows, not the whole locus). Fall back
            # to the gradient sign as a weak prior, ranked below anything measured.
            g = candidate_signed_grad.get(c.candidate_id, 0.0)
            sign_ok = (g > 0) if desired_direction == "down" else (g < 0)
            mult = 1.0 if sign_ok else wrong_sign_multiplier
            return mult * c.gradient_saliency_score - penalty

        available = [
            c for c in shortlist
            if c.candidate_id not in evaluated_ids
            and not (
                _variant_key(c) in tried_variant_min_bystanders
                and len(c.bystander_edits_in_window) >= tried_variant_min_bystanders[_variant_key(c)]
            )
        ]
        if not available:
            break
        available.sort(key=_rank_key, reverse=True)
        chosen = available[0]
        evaluated_ids.add(chosen.candidate_id)
        vk = _variant_key(chosen)
        tried_variant_min_bystanders[vk] = min(
            tried_variant_min_bystanders.get(vk, 10**6), len(chosen.bystander_edits_in_window)
        )

        # Step 3: Scoped ISM. Usually already measured by the pre-screen, so this costs nothing;
        # a candidate outside every pre-screen window gets its own compact window here.
        if vk in measured_lfc:
            ism_val = measured_lfc[vk]
        else:
            edit_pos_0 = chosen.target_edit.position - 1
            try:
                ism_map = backend.score_ism_window(
                    parent_interval=parent_iv,
                    ism_window=GenomicInterval(
                        chromosome=parent_iv.chromosome,
                        start=edit_pos_0 - 2,
                        end=edit_pos_0 + 4,
                    ),
                    gene_id=gene_id,
                    cell_sketch=target_sketches.sketch_a,
                    condition_name=target_sketches.cell_state_name,
                )
                measured_lfc.update(ism_map)
                ism_val = ism_map.get(vk, 0.0)
            except Exception:
                ism_val = 0.0  # ISM budget exhausted; the two-tier verifier is still authoritative
        chosen.ism_verification_score = round(ism_val, 6)

        # Step 4: Two-Tier Verification (Held-out Sketch B + Reference Hazards)
        report = evaluate_candidate_two_tier(
            backend=backend,
            parent_interval=parent_iv,
            candidate=chosen,
            gene_id=gene_id,
            target_sketches=target_sketches,
            healthy_sketches_map=healthy_sketches_map,
            reference_curies=reference_curies,
            thresholds=thresholds,
            tissue_matched_curies=tissue_matched_curies,
        )

        if report.all_checks_passed:
            diagnosis = (
                f"CONVERGED: Candidate {chosen.candidate_id} passed all Tier 1 (held-out Sketch B "
                f"LFC={report.tier1_selectivity.target_delta_sketch_b:+.5f}) and Tier 2 hazard checks."
            )
            rule_added = "NONE (All Tier 1 & Tier 2 constraints satisfied)"
            traces.append(
                SearchStepTrace(
                    iteration=step_idx,
                    report=report,
                    diagnosis=diagnosis,
                    constraint_added=rule_added,
                )
            )
            winning_report = report
            break

        # failure analysis analyzes failure reasons and updates search constraints
        reasons_joined = "; ".join(report.failure_reasons)
        new_rules: list[str] = []
        if not report.tier2_hazards.bystander_clean_passed:
            require_bystander_free = True
            new_rules.append("ENFORCE_ZERO_BYSTANDERS_IN_WINDOW_[4..8]")
        if not report.tier1_selectivity.healthy_selectivity_passed:
            healthy_gradient_penalty_weight = max(1.5, healthy_gradient_penalty_weight + 1.25)
            new_rules.append(
                f"INCREASE_DIFFERENTIAL_GRADIENT_PENALTY (lambda={healthy_gradient_penalty_weight:.2f})"
            )
        if not report.tier1_selectivity.neighbor_specificity_passed:
            new_rules.append(
                f"AVOID_SHARED_REGULATORY_ELEMENTS (collateral on "
                f"'{report.tier1_selectivity.worst_neighbor_gene}')"
            )
        if any(r.startswith("TIER1_DIRECTION") for r in report.failure_reasons):
            wrong_sign_multiplier = max(0.02, wrong_sign_multiplier * 0.2)
            new_rules.append("STRONGLY_PREFER_GRADIENT_SIGN_THAT_MOVES_TARGET_IN_DESIRED_DIRECTION")
        elif not report.tier1_selectivity.target_efficacy_passed:
            healthy_gradient_penalty_weight = max(0.5, healthy_gradient_penalty_weight * 0.7)
            new_rules.append("PRIORITIZE_HIGH_ISM_JACOBIAN_PEAKS")

        if not new_rules:
            new_rules.append("EXCLUDE_LOCUS_AND_STEP_TO_NEXT_SALIENT_PROTOSpacer")

        traces.append(
            SearchStepTrace(
                iteration=step_idx,
                report=report,
                diagnosis=f"FAILED [{len(report.failure_reasons)} check(s)]: {reasons_joined}",
                constraint_added=" + ".join(new_rules),
            )
        )

    # If no single step passed every strict threshold, select highest composite bystander-clean report
    if winning_report is None and traces:
        clean_traces = [t.report for t in traces if t.report.tier2_hazards.bystander_clean_passed]
        pool = clean_traces if clean_traces else [t.report for t in traces]
        winning_report = max(pool, key=lambda r: r.composite_objective_score)

    # Off-tissue safety panel on the finalist only. Tier 1 asks whether the edit is selective
    # between cell *states*; this asks whether it does anything in organs the therapy would reach
    # but never intends to touch. One API call per tissue, so finalists only.
    if winning_report is not None and run_off_tissue_panel_on_finalist:
        try:
            winning_report.off_tissue_panel = run_off_tissue_panel(
                backend=backend,
                parent_interval=parent_iv,
                candidate=winning_report.candidate,
                gene_id=gene_id,
                on_target_abs_lfc=abs(winning_report.tier1_selectivity.target_delta_sketch_b),
            )
        except Exception:
            pass  # budget exhausted: the panel is supplementary, not a gate on reporting

    prov = build_provenance_record(
        backend_mode=backend.backend_mode,
        model_version=backend.model_version,
        random_seed=random_seed,
        budget=budget,
    )

    return OptimizationRunResult(
        scenario_id=scenario_id,
        gene_symbol=gene_symbol,
        locus=parent_iv,
        provenance=prov,
        sketch_metrics=sketch_metrics,
        iterations=traces,
        selected_candidate=winning_report,
        converged=bool(winning_report and winning_report.all_checks_passed),
        desired_direction=desired_direction,  # type: ignore[arg-type]
        efficacy_null=efficacy_null,
        ism_prescreen=prescreen_summary,
        data_provenance=data_provenance or {},
    )
