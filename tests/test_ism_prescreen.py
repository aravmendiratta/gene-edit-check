"""The ISM pre-screen must measure edit direction before the verifier spends calls on it."""

from __future__ import annotations

import numpy as np

from gene_edit_check.ag_backend.protocol import GradientProfile, VariantMultiModalScore
from gene_edit_check.search.search_loop import run_gradient_ism_search
from gene_edit_check.budget import RunBudgetTracker
from gene_edit_check.models import BackendMode, GenomicInterval
from gene_edit_check.scenarios.registry import SCENARIOS, build_til_cohort_anndata
from gene_edit_check.sc_sketch.sketch_variance import build_disjoint_sketches_and_ablation

RNG = np.random.default_rng(7)


class _SignFlipBackend:
    """A locus where the gradient sign is *systematically wrong* about direction.

    Every ISM measurement says 'this edit raises the gene' (positive LFC) while the gradient
    claims the opposite. A search that trusts gradients burns its whole iteration budget here;
    one that pre-screens with ISM learns the truth in a few calls.
    """

    backend_mode = BackendMode.SURROGATE_SIMULATOR
    model_version = "fake-signflip"

    def __init__(self, budget: RunBudgetTracker, seq: str, interval: GenomicInterval, ism_lfc: float = +0.30) -> None:
        self._budget, self._seq, self._iv = budget, seq, interval
        self._ism_lfc = ism_lfc
        self.ism_calls = 0
        self.variant_calls = 0

    def encode_cells(self, adata, state_name):  # type: ignore[no-untyped-def]
        return RNG.normal(size=(adata.n_obs, 32)).astype(np.float32)

    def compute_interval_gradients(self, gene_id, gene_symbol, cell_sketch, condition_name, sequence_length=131072):  # type: ignore[no-untyped-def]
        # Positive gradients everywhere => the heuristic predicts every edit lowers the gene.
        grads = np.full(len(self._seq), 0.5, dtype=np.float32)
        return GradientProfile(gene_id, gene_symbol, condition_name, self._iv, self._seq, grads, 4.0)

    def score_ism_window(self, parent_interval, ism_window, gene_id, cell_sketch, condition_name):  # type: ignore[no-untyped-def]
        self.ism_calls += 1
        w = ism_window.width
        self._budget.check_and_record_ism_window(window_bp=w, n_variants=w * 3)
        out = {}
        for pos0 in range(ism_window.start, ism_window.end):
            for ref_b in "ACGT":
                for alt_b in "ACGT":
                    if ref_b != alt_b:
                        out[(pos0 + 1, ref_b, alt_b)] = self._ism_lfc
        return out

    def score_variant_multimodal(self, parent_interval, variant, gene_id, cell_sketch, condition_name, reference_curies=None):  # type: ignore[no-untyped-def]
        self.variant_calls += 1
        self._budget.check_and_record_variant_score(1)
        return VariantMultiModalScore(variant, condition_name, self._ism_lfc, 0.001, 0.001, 0.001)


def _sketches(backend, cohort, name, size=16):  # type: ignore[no-untyped-def]
    adata = cohort[name]
    emb = backend.encode_cells(adata, name)
    return build_disjoint_sketches_and_ablation(adata, emb, name, sketch_size=size, seed=42)


def _run(prescreen_windows: int, ism_lfc: float = +0.30):  # type: ignore[no-untyped-def]
    spec = SCENARIOS["pd1-exhaustion"]
    budget = RunBudgetTracker()
    seq = "".join(RNG.choice(list("ACGT"), size=6000))
    iv = GenomicInterval(chromosome="chr2", start=241_850_000, end=241_856_000)
    backend = _SignFlipBackend(budget, seq, iv, ism_lfc=ism_lfc)
    cohort = build_til_cohort_anndata(n_cells_per_state=64, seed=42)
    t_sk, t_m = _sketches(backend, cohort, spec.target_state_name)
    h_sk, h_m = _sketches(backend, cohort, "Healthy_Naive_CD8")
    result = run_gradient_ism_search(
        scenario_id=spec.scenario_id, gene_id=spec.gene_id, gene_symbol=spec.gene_symbol,
        backend=backend, budget=budget, target_sketches=t_sk,
        healthy_sketches_map={"Healthy_Naive_CD8": h_sk}, sketch_metrics=[t_m, h_m],
        reference_curies=list(spec.reference_curies), max_iterations=4, random_seed=42,
        desired_direction="down", null_n=0, ism_prescreen_windows=prescreen_windows,
        run_off_tissue_panel_on_finalist=False,
    )
    return result, backend


def test_prescreen_measures_direction_and_records_the_audit() -> None:
    result, _ = _run(prescreen_windows=6)
    ps = result.ism_prescreen
    assert ps is not None
    assert ps.n_windows == 6
    assert ps.n_variants_measured > 0
    assert ps.n_candidates_with_measurement > 0
    # Every measurement points the wrong way for a silencing goal, so all are filtered.
    assert ps.n_correct_direction == 0
    assert ps.n_wrong_direction_filtered == ps.n_candidates_with_measurement
    # And the audit catches that the gradient sign was wrong every single time.
    assert ps.gradient_sign_agreement == 0.0


def test_prescreen_can_be_switched_off() -> None:
    result, _ = _run(prescreen_windows=0)
    assert result.ism_prescreen is None


def test_known_wrong_candidates_are_ranked_below_unmeasured_ones() -> None:
    """When the pre-screen proves an edit moves the gene the wrong way, the loop should go
    looking elsewhere rather than verify it -- so it spends ISM on fresh windows instead."""
    result, backend = _run(prescreen_windows=6)
    assert backend.ism_calls == 6 + len(result.iterations)
    assert not result.converged  # every edit at this locus raises the gene


def test_measurements_are_reused_rather_than_re_scored() -> None:
    """With the pre-screen pointing the right way, verification iterations add no ISM calls:
    the measurement the loop needs is already in hand."""
    result, backend = _run(prescreen_windows=6, ism_lfc=-0.30)
    assert backend.ism_calls == 6
    assert result.iterations
    ps = result.ism_prescreen
    assert ps is not None and ps.n_correct_direction == ps.n_candidates_with_measurement
    assert ps.gradient_sign_agreement == 1.0
    assert ps.best_measured_lfc == -0.3
