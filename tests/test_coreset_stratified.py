"""Tests for donor-stratified coreset selection in AG2 latent space & Sketch A/B disjointness."""

from __future__ import annotations

from crispr_reflexion.ag_backend.simulator import SurrogateSimulatorBackend
from crispr_reflexion.budget import RunBudgetTracker
from crispr_reflexion.scenarios.registry import build_til_cohort_anndata
from crispr_reflexion.sc_sketch.sketch_variance import build_disjoint_sketches_and_ablation


def test_disjoint_sketches_and_mmd_improvement() -> None:
    cohort = build_til_cohort_anndata(n_cells_per_state=96, seed=42)
    adata = cohort["Exhausted_CD8_TIL"]
    backend = SurrogateSimulatorBackend(RunBudgetTracker(), seed=42)
    emb = backend.encode_cells(adata, "Exhausted_CD8_TIL")

    sketches, metrics = build_disjoint_sketches_and_ablation(
        adata=adata,
        embeddings=emb,
        cell_state_name="Exhausted_CD8_TIL",
        sketch_size=16,
        seed=42,
    )

    # Verify strict cell disjointness between Proposal Sketch A and Held-Out Sketch B
    cells_a = set(sketches.sketch_a.obs_names)
    cells_b = set(sketches.sketch_b.obs_names)
    assert len(cells_a) == 16
    assert len(cells_b) == 16
    assert cells_a.isdisjoint(cells_b)

    # Verify Facility Location coreset achieves lower MMD^2 than random subsampling
    assert metrics.ag2_latent_mmd_stratified < metrics.ag2_latent_mmd_unstratified_random
    assert metrics.mmd_reduction_pct > 0.0
