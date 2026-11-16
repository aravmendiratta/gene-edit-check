"""Single assembly point for a scenario run (used by the CLI and the MCP server)."""

from __future__ import annotations

from dataclasses import dataclass

import anndata

from crispr_reflexion.ag_backend.protocol import VerifierBackend
from crispr_reflexion.agent.reflexion_loop import run_gradient_ism_reflexion_optimization
from crispr_reflexion.benchmarks.hazard_coupling import analyze_hazard_coupling
from crispr_reflexion.budget import RunBudgetTracker
from crispr_reflexion.datasets import crc_pelka
from crispr_reflexion.models import OptimizationRunResult, SketchQualityMetrics
from crispr_reflexion.sc_sketch.coreset import DisjointCellSketches
from crispr_reflexion.sc_sketch.sketch_variance import build_disjoint_sketches_and_ablation
from crispr_reflexion.scenarios.registry import (
    IMMUNE_ENSEMBL_PANEL,
    ScenarioSpec,
    build_til_cohort_anndata,
)

SYNTHETIC_NOTICE = (
    "SYNTHETIC DATA: cells are generated from hand-set Poisson rates over 24 genes. "
    "For tests/demos only; results say nothing about biology."
)


@dataclass
class Cohort:
    states: dict[str, anndata.AnnData]
    provenance: dict[str, object]
    observed_cp10k: dict[str, dict[str, float]] | None


def load_cohort(spec: ScenarioSpec, pool_size: int, seed: int) -> Cohort:
    if spec.data_source == "real":
        real = crc_pelka.build_cohort(
            pool_size=pool_size,
            seed=seed,
            marker_gene_ids=tuple(g for g, _ in IMMUNE_ENSEMBL_PANEL),
        )
        return Cohort(real.states, real.provenance, real.observed_cp10k)
    states = build_til_cohort_anndata(n_cells_per_state=pool_size, seed=seed)
    return Cohort(states, {"data_source": SYNTHETIC_NOTICE}, None)


def build_all_sketches(
    backend: VerifierBackend,
    spec: ScenarioSpec,
    cohort: Cohort,
    sketch_size: int,
    seed: int,
) -> tuple[DisjointCellSketches, dict[str, DisjointCellSketches], list[SketchQualityMetrics]]:
    def _one(name: str) -> tuple[DisjointCellSketches, SketchQualityMetrics]:
        adata = cohort.states[name]
        emb = backend.encode_cells(adata, state_name=name)
        return build_disjoint_sketches_and_ablation(
            adata=adata, embeddings=emb, cell_state_name=name, sketch_size=sketch_size, seed=seed
        )

    target, t_metrics = _one(spec.target_state_name)
    healthy: dict[str, DisjointCellSketches] = {}
    metrics = [t_metrics]
    for h in spec.healthy_state_names:
        healthy[h], m = _one(h)
        metrics.append(m)
    return target, healthy, metrics


def run_scenario(
    spec: ScenarioSpec,
    backend: VerifierBackend,
    budget: RunBudgetTracker,
    sketch_size: int,
    max_iterations: int,
    seed: int,
    pool_size: int,
    null_n: int = 12,
    ism_prescreen_windows: int = 10,
    off_tissue_panel: bool = True,
) -> OptimizationRunResult:
    cohort = load_cohort(spec, pool_size=pool_size, seed=seed)
    target, healthy, metrics = build_all_sketches(backend, spec, cohort, sketch_size, seed)
    result = run_gradient_ism_reflexion_optimization(
        scenario_id=spec.scenario_id,
        gene_id=spec.gene_id,
        gene_symbol=spec.gene_symbol,
        backend=backend,
        budget=budget,
        target_sketches=target,
        healthy_sketches_map=healthy,
        sketch_metrics=metrics,
        reference_curies=list(spec.reference_curies),
        tissue_matched_curies=list(spec.tissue_matched_curies) or None,
        max_iterations=max_iterations,
        random_seed=seed,
        desired_direction=spec.desired_direction,
        null_n=null_n,
        ism_prescreen_windows=ism_prescreen_windows,
        run_off_tissue_panel_on_finalist=off_tissue_panel,
        data_provenance=cohort.provenance,
    )
    # Why did the search end where it did? Distinguishes a weak optimiser from a locus where
    # efficacy and a hazard are inseparable. Uses variants already scored, so it is free.
    result.hazard_coupling = analyze_hazard_coupling(result).model_dump()
    return result
