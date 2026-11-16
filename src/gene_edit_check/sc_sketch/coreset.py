"""Donor-Stratified Facility Location Coreset Sketcher in AlphaGenome 2 Latent Space (`sc_sketch/coreset.py`).

Pipeline:
1. Takes an `AnnData` object with raw integer UMI counts (`.X`) and Ensembl IDs (`ENSG...`) in `.var_names`.
2. Stratifies candidate pools by `donor_id` / `batch` metadata (since `encode_cells()` does not
   batch-correct).
3. Splits the pool into two disjoint partitions (`Pool A` and `Pool B`, stratified across donors)
   so that proposal (`Sketch A`) and final verification (`Sketch B`) never share cells.
4. Calls `encode_cells()` to obtain AlphaGenome 2's 256-D latent embeddings.
5. Runs Lazy-Greedy Facility Location in AG2 latent space to select the optimal subset of size
   `k <= MAX_CONDITIONING_CELLS (256)` for both `Sketch A` and `Sketch B`.
"""

from __future__ import annotations

from dataclasses import dataclass
import heapq
import anndata
import numpy as np
import pandas as pd


@dataclass
class DisjointCellSketches:
    cell_state_name: str
    sketch_a: anndata.AnnData
    sketch_b: anndata.AnnData
    sketch_a_embeddings: np.ndarray
    sketch_b_embeddings: np.ndarray
    full_pool_embeddings: np.ndarray
    donor_stratified: bool
    split_by: str = "cell"


def compute_rbf_similarity_matrix(embeddings: np.ndarray, gamma: float | None = None) -> np.ndarray:
    """Compute pairwise RBF kernel matrix K(z_i, z_j) = exp(-gamma * ||z_i - z_j||^2)."""
    sq_norms = np.sum(embeddings**2, axis=1, keepdims=True)
    sq_dists = np.maximum(0.0, sq_norms + sq_norms.T - 2.0 * (embeddings @ embeddings.T))
    if gamma is None:
        median_dist = float(np.median(sq_dists))
        gamma = 1.0 / max(median_dist, 1e-6)
    return np.exp(-gamma * sq_dists)


def lazy_greedy_facility_location(
    embeddings: np.ndarray,
    k: int,
    donor_labels: np.ndarray | None = None,
) -> np.ndarray:
    """Select `k` indices maximizing the monotone submodular Facility Location objective
    in AlphaGenome 2's `encode_cells()` embedding space, with optional donor-balance penalty.
    """
    n = embeddings.shape[0]
    if k >= n:
        return np.arange(n, dtype=np.int64)

    sim_matrix = compute_rbf_similarity_matrix(embeddings)
    mean_sim_to_pool = np.mean(sim_matrix, axis=1)  # mu_j = (1/N) sum_i K(i, j)
    current_max_sim = np.zeros(n, dtype=np.float64)
    running_kernel_sum = np.zeros(n, dtype=np.float64)
    selected: list[int] = []
    selected_mask = np.zeros(n, dtype=bool)

    # Track donor quotas if donor stratification is enabled
    donor_counts: dict[str, int] = {}
    donor_cap: dict[str, int] = {}
    if donor_labels is not None:
        unique_donors, counts = np.unique(donor_labels, return_counts=True)
        for d, c in zip(unique_donors, counts):
            donor_counts[str(d)] = 0
            donor_cap[str(d)] = max(1, int(np.ceil(k * (c / n))))

    for step in range(k):
        t_next = step + 1
        # Hybrid Submodular Facility Location + Kernel Herding (MMD^2 minimization) objective:
        # Maximizes mode coverage while minimizing empirical MMD^2 in AG2 latent space
        herding_gain = 2.0 * mean_sim_to_pool - (
            (2.0 * running_kernel_sum + np.diag(sim_matrix)) / float(t_next)
        )
        fl_gain = np.mean(np.maximum(0.0, sim_matrix - current_max_sim[:, None]), axis=0)
        composite_gain = herding_gain + 0.35 * fl_gain

        # Mask already selected or donor-capped candidates
        valid_mask = ~selected_mask
        if donor_labels is not None:
            for idx in range(n):
                if valid_mask[idx] and donor_counts[str(donor_labels[idx])] >= donor_cap[str(donor_labels[idx])]:
                    valid_mask[idx] = False
            if not np.any(valid_mask):
                valid_mask = ~selected_mask

        scores = np.where(valid_mask, composite_gain, -1e18)
        best_idx = int(np.argmax(scores))
        selected.append(best_idx)
        selected_mask[best_idx] = True
        current_max_sim = np.maximum(current_max_sim, sim_matrix[:, best_idx])
        running_kernel_sum += sim_matrix[:, best_idx]
        if donor_labels is not None:
            donor_counts[str(donor_labels[best_idx])] += 1

    return np.sort(np.array(selected[:k], dtype=np.int64))
