"""Ablation & Cross-Sketch Variance Metrics (`sc_sketch/sketch_variance.py`).

Quantifies:
1. Maximum Mean Discrepancy (MMD) in AlphaGenome 2's `encode_cells()` latent space:
   Donor-Stratified Facility Location vs. Unstratified Random Subsampling
   (DeepMind's default `gene_counts_condition.subsample`).
2. Cross-Sketch Cosine Consistency between disjoint `Sketch A` and `Sketch B`.
"""

from __future__ import annotations

import anndata
import numpy as np

from crispr_reflexion.adata_utils import drop_unexpressed_genes
from crispr_reflexion.models import SketchQualityMetrics
from crispr_reflexion.sc_sketch.coreset import (
    DisjointCellSketches,
    compute_rbf_similarity_matrix,
    lazy_greedy_facility_location,
)


def compute_mmd_squared(full_embeddings: np.ndarray, subset_indices: np.ndarray) -> float:
    """Compute unbiased RBF Maximum Mean Discrepancy (MMD^2) between full pool and sketch."""
    sim_all = compute_rbf_similarity_matrix(full_embeddings)
    k_xx = float(np.mean(sim_all))
    k_yy = float(np.mean(sim_all[np.ix_(subset_indices, subset_indices)]))
    k_xy = float(np.mean(sim_all[:, subset_indices]))
    return max(0.0, k_xx + k_yy - 2.0 * k_xy)


def build_disjoint_sketches_and_ablation(
    adata: anndata.AnnData,
    embeddings: np.ndarray,
    cell_state_name: str,
    sketch_size: int = 32,
    donor_col: str = "donor_id",
    seed: int = 42,
    min_donors_for_donor_split: int = 4,
) -> tuple[DisjointCellSketches, SketchQualityMetrics]:
    """Split `adata` into disjoint Proposal Sketch A and Held-Out Verification Sketch B
    in AlphaGenome 2's latent space, and run the stratified vs. unstratified random ablation.
    """
    if sketch_size > 256:
        raise ValueError(f"sketch_size ({sketch_size}) cannot exceed MAX_CONDITIONING_CELLS (256).")

    n_total = adata.shape[0]
    rng = np.random.default_rng(seed)
    donor_labels = (
        adata.obs[donor_col].to_numpy()
        if donor_col in adata.obs.columns
        else None
    )

    # Split into Pool A / Pool B. With enough donors the split is BY DONOR (Pool B patients never
    # appear in Pool A), so Sketch A -> Sketch B transfer tests generalization across people.
    # With few donors we fall back to a within-donor cell split (and report split_by="cell").
    pool_a_idx: list[int] = []
    pool_b_idx: list[int] = []
    split_by = "cell"
    if donor_labels is not None and len(np.unique(donor_labels)) >= min_donors_for_donor_split:
        split_by = "donor"
        donors, sizes = np.unique(donor_labels, return_counts=True)
        order = rng.permutation(len(donors))
        order = order[np.argsort(-sizes[order], kind="stable")]  # big donors first, greedy balance
        load_a = load_b = 0
        for i in order:
            idx = np.where(donor_labels == donors[i])[0]
            if load_a <= load_b:
                pool_a_idx.extend(idx.tolist()); load_a += len(idx)
            else:
                pool_b_idx.extend(idx.tolist()); load_b += len(idx)
    elif donor_labels is not None:
        for d in np.unique(donor_labels):
            d_indices = np.where(donor_labels == d)[0]
            rng.shuffle(d_indices)
            mid = len(d_indices) // 2
            pool_a_idx.extend(d_indices[:mid].tolist())
            pool_b_idx.extend(d_indices[mid:].tolist())
    else:
        perm = rng.permutation(n_total)
        mid = n_total // 2
        pool_a_idx = perm[:mid].tolist()
        pool_b_idx = perm[mid:].tolist()

    pool_a_arr = np.sort(np.array(pool_a_idx, dtype=np.int64))
    pool_b_arr = np.sort(np.array(pool_b_idx, dtype=np.int64))

    k_eff = min(sketch_size, len(pool_a_arr), len(pool_b_arr))

    # Select Sketch A and Sketch B using Facility Location in AG2 latent space
    # (Pool A and Pool B are already donor-stratified above)
    local_a = lazy_greedy_facility_location(
        embeddings[pool_a_arr],
        k=k_eff,
        donor_labels=donor_labels[pool_a_arr] if donor_labels is not None else None,
    )
    local_b = lazy_greedy_facility_location(
        embeddings[pool_b_arr],
        k=k_eff,
        donor_labels=donor_labels[pool_b_arr] if donor_labels is not None else None,
    )

    global_idx_a = pool_a_arr[local_a]
    global_idx_b = pool_b_arr[local_b]

    # Compute MMD of our Stratified Coreset vs. 25 trials of Uniform Random Sampling on Pool A
    pool_a_emb = embeddings[pool_a_arr]
    mmd_strat = compute_mmd_squared(pool_a_emb, local_a)
    random_mmds = []
    for _ in range(25):
        rand_idx = rng.choice(len(pool_a_arr), size=k_eff, replace=False)
        random_mmds.append(compute_mmd_squared(pool_a_emb, rand_idx))
    mmd_rand = float(np.mean(random_mmds))
    reduction_pct = (
        100.0 * (mmd_rand - mmd_strat) / max(mmd_rand, 1e-9)
        if mmd_rand > 0
        else 0.0
    )

    # Cosine consistency between Sketch A centroid and Sketch B centroid in AG2 latent space
    mean_a = np.mean(embeddings[global_idx_a], axis=0)
    mean_b = np.mean(embeddings[global_idx_b], axis=0)
    denom = float(np.linalg.norm(mean_a) * np.linalg.norm(mean_b))
    cosine_sim = float(np.dot(mean_a, mean_b) / max(denom, 1e-9))

    sketches = DisjointCellSketches(
        cell_state_name=cell_state_name,
        sketch_a=drop_unexpressed_genes(adata[global_idx_a]),
        sketch_b=drop_unexpressed_genes(adata[global_idx_b]),
        sketch_a_embeddings=embeddings[global_idx_a],
        sketch_b_embeddings=embeddings[global_idx_b],
        full_pool_embeddings=embeddings,
        donor_stratified=donor_labels is not None,
        split_by=split_by,
    )
    metrics = SketchQualityMetrics(
        cell_state_name=cell_state_name,
        total_source_cells=n_total,
        sketch_size=k_eff,
        donor_stratified=donor_labels is not None,
        split_by=split_by,
        n_donors_sketch_a=int(len(np.unique(donor_labels[global_idx_a]))) if donor_labels is not None else 0,
        n_donors_sketch_b=int(len(np.unique(donor_labels[global_idx_b]))) if donor_labels is not None else 0,
        donors_overlap_a_b=int(len(set(donor_labels[global_idx_a]) & set(donor_labels[global_idx_b]))) if donor_labels is not None else 0,
        ag2_latent_mmd_stratified=round(mmd_strat, 6),
        ag2_latent_mmd_unstratified_random=round(mmd_rand, 6),
        mmd_reduction_pct=round(reduction_pct, 2),
        sketch_a_B_cosine_consistency=round(cosine_sim, 5),
    )
    return sketches, metrics
