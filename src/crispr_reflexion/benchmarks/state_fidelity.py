"""Cell-state fidelity check (`benchmarks/state_fidelity.py`).

Replaces the earlier hand-typed "published ground truth" panel, which had no source.

Question answered with real numbers only: when AlphaGenome 2 is conditioned on a sketch of real
tumor-exhausted CD8 cells vs a sketch of real normal-colon CD8 cells, does its predicted
gene-level fold-change match the fold-change *observed in those cells' full populations*?

- Predicted: `GradientProfile.scalar_score` = natural log of the model's predicted expression for
  the gene (alphagenome `InputGradientsOutput` docs), so the difference between two conditions is a
  predicted natural-log fold-change.
- Observed: ln((mean CP10k in target + c) / (mean CP10k in healthy + c)) over the FULL state
  populations (thousands of cells, dozens of patients), not just the sketch.

This validates the *conditioning* (does the model "see" the cell state), not the base-edit
effect sizes, which have no ground truth here. Only meaningful on the LIVE backend.
"""

from __future__ import annotations

from dataclasses import dataclass

import anndata
import numpy as np
from scipy.stats import spearmanr

from crispr_reflexion.ag_backend.protocol import VerifierBackend


@dataclass
class FidelityResult:
    n_genes: int
    spearman_rho: float
    spearman_p: float
    sign_concordance_pct: float
    n_concordance_genes: int
    slope_pred_vs_obs: float  # <1 means the model compresses the true difference
    records: list[dict[str, object]]


def run_state_fidelity(
    backend: VerifierBackend,
    target_sketch: anndata.AnnData,
    healthy_sketch: anndata.AnnData,
    observed_cp10k: dict[str, dict[str, float]],
    target_state: str,
    healthy_state: str,
    genes: list[tuple[str, str]],
    pseudocount: float = 0.1,
    min_abs_obs_for_concordance: float = 0.5,
) -> FidelityResult:
    records: list[dict[str, object]] = []
    for gene_id, symbol in genes:
        gt = backend.compute_interval_gradients(gene_id, symbol, target_sketch, target_state)
        gh = backend.compute_interval_gradients(gene_id, symbol, healthy_sketch, healthy_state)
        obs_t = observed_cp10k[target_state][gene_id]
        obs_h = observed_cp10k[healthy_state][gene_id]
        records.append(
            {
                "gene": symbol,
                "obs_cp10k_target": obs_t,
                "obs_cp10k_healthy": obs_h,
                "observed_ln_fc": float(np.log((obs_t + pseudocount) / (obs_h + pseudocount))),
                "predicted_ln_fc": float(gt.scalar_score - gh.scalar_score),
            }
        )
    obs = np.array([r["observed_ln_fc"] for r in records], dtype=float)
    pred = np.array([r["predicted_ln_fc"] for r in records], dtype=float)
    rho, p = spearmanr(obs, pred)
    strong = np.abs(obs) >= min_abs_obs_for_concordance
    concord = float(np.mean(np.sign(obs[strong]) == np.sign(pred[strong])) * 100) if strong.any() else float("nan")
    slope = float(np.dot(obs, pred) / max(np.dot(obs, obs), 1e-12))
    return FidelityResult(
        n_genes=len(records),
        spearman_rho=float(rho),
        spearman_p=float(p),
        sign_concordance_pct=concord,
        n_concordance_genes=int(strong.sum()),
        slope_pred_vs_obs=slope,
        records=records,
    )
