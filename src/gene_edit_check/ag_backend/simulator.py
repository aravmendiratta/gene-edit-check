"""Offline Surrogate Backend (`ag_backend/simulator.py`).

Implements the `VerifierBackend` Protocol for deterministic unit tests and zero-key CI environments.
All outputs are explicitly stamped with `BackendMode.SURROGATE_SIMULATOR` in provenance records.
"""

from __future__ import annotations

import hashlib
import anndata
import numpy as np

from gene_edit_check.ag_backend.protocol import (
    GradientProfile,
    VariantMultiModalScore,
    VerifierBackend,
)
from gene_edit_check.adata_utils import counts_fingerprint, dense_counts
from gene_edit_check.budget import RunBudgetTracker
from gene_edit_check.models import BackendMode, GenomicInterval, SingleNucleotideVariant


class SurrogateSimulatorBackend(VerifierBackend):
    def __init__(self, budget: RunBudgetTracker, seed: int = 42) -> None:
        self._budget = budget
        self._seed = seed

    @property
    def backend_mode(self) -> BackendMode:
        return BackendMode.SURROGATE_SIMULATOR

    @property
    def model_version(self) -> str:
        return "surrogate-simulator-v2.0 (OFFLINE)"

    def encode_cells(self, adata: anndata.AnnData, state_name: str) -> np.ndarray:
        self._budget.check_and_record_encode_cells(adata.shape[0])
        x = dense_counts(adata)
        log_x = np.log1p(x)
        rng = np.random.default_rng(self._seed)
        proj = rng.normal(0.0, 1.0 / np.sqrt(max(1, x.shape[1])), size=(x.shape[1], 256)).astype(np.float32)
        emb = log_x @ proj
        norms = np.linalg.norm(emb, axis=1, keepdims=True)
        return emb / np.maximum(norms, 1e-6)

    def compute_interval_gradients(
        self,
        gene_id: str,
        gene_symbol: str,
        cell_sketch: anndata.AnnData,
        condition_name: str,
        sequence_length: int = 131072,
    ) -> GradientProfile:
        self._budget.check_and_record_gradient_call(1)
        start = 241_793_329
        end = start + sequence_length
        iv = GenomicInterval(chromosome="chr2", start=start, end=end)

        # Deterministic synthetic sequence with embedded PAMs and TF motifs around center
        rng = np.random.default_rng(self._seed)
        bases = np.array(["A", "C", "G", "T"])
        seq_arr = rng.choice(bases, size=sequence_length)
        mid = sequence_length // 2

        # Plant canonical motifs + PAMs near mid-200 (shared RUNX1) and mid+160 (exhaustion TOX/NR4A1)
        shared_window = list("GCTATGTGGTTAGGCCAGGCTGG")
        selective_window = list("CGATGCAAATGGAAAGGTCAGGG")
        seq_arr[mid - 200 : mid - 200 + len(shared_window)] = shared_window
        seq_arr[mid + 160 : mid + 160 + len(selective_window)] = selective_window
        ref_seq = "".join(seq_arr.tolist())

        grads = rng.normal(0.0, 0.002, size=sequence_length).astype(np.float32)
        is_target = "exhausted" in condition_name.lower() or "target" in condition_name.lower()
        # Shared motif has strong gradient in both target and healthy
        grads[mid - 200 : mid - 180] = 0.085
        # Selective motif has strong gradient ONLY in target state
        if is_target:
            grads[mid + 160 : mid + 182] = 0.115
        else:
            grads[mid + 160 : mid + 182] = 0.003

        return GradientProfile(
            gene_id=gene_id,
            gene_name=gene_symbol,
            condition_name=condition_name,
            interval=iv,
            reference_sequence=ref_seq,
            signed_gradients_1bp=grads,
            scalar_score=4.85 if is_target else 2.10,
        )

    def score_ism_window(
        self,
        parent_interval: GenomicInterval,
        ism_window: GenomicInterval,
        gene_id: str,
        cell_sketch: anndata.AnnData,
        condition_name: str,
    ) -> dict[tuple[int, str, str], float]:
        w = ism_window.width
        self._budget.check_and_record_ism_window(window_bp=w, n_variants=w * 3)
        out: dict[tuple[int, str, str], float] = {}
        mid_genomic = parent_interval.start + (parent_interval.width // 2)
        is_target = "exhausted" in condition_name.lower() or "target" in condition_name.lower()

        for pos_0 in range(ism_window.start, ism_window.end):
            pos_1 = pos_0 + 1
            for ref_b in ("A", "C", "G", "T"):
                for alt_b in ("A", "C", "G", "T"):
                    if ref_b == alt_b:
                        continue
                    dist_sel = abs(pos_0 - (mid_genomic + 165))
                    dist_shr = abs(pos_0 - (mid_genomic - 195))
                    if dist_sel <= 10:
                        val = -0.042 if is_target else -0.001
                    elif dist_shr <= 10:
                        val = -0.038
                    else:
                        val = -0.002
                    out[(pos_1, ref_b, alt_b)] = val
        return out

    def score_variant_multimodal(
        self,
        parent_interval: GenomicInterval,
        variant: SingleNucleotideVariant,
        gene_id: str,
        cell_sketch: anndata.AnnData | None,
        condition_name: str,
        reference_curies: list[str] | None = None,
    ) -> VariantMultiModalScore:
        self._budget.check_and_record_variant_score(1)
        mid_genomic = parent_interval.start + (parent_interval.width // 2)
        pos_0 = variant.position - 1
        dist_sel = abs(pos_0 - (mid_genomic + 165))
        dist_shr = abs(pos_0 - (mid_genomic - 195))
        is_target = "exhausted" in condition_name.lower() or "target" in condition_name.lower()

        # Add slight deterministic sketch hash perturbation to model Sketch A vs Sketch B
        sketch_jitter = 0.0
        if cell_sketch is not None:
            h_int = int(counts_fingerprint(cell_sketch)[:6], 16)
            sketch_jitter = ((h_int % 100) - 50) * 1e-5

        if dist_sel <= 12:
            lfc = (-0.044 if is_target else -0.0012) + sketch_jitter
            cpg = 0.004
            splice = 0.002
        elif dist_shr <= 12:
            lfc = (-0.039 if is_target else -0.036) + sketch_jitter
            cpg = 0.009
            splice = 0.005
        else:
            lfc = -0.004 + sketch_jitter
            cpg = 0.003
            splice = 0.002

        # Two stand-in neighbours so the neighbour-specificity check is exercised offline. The
        # shared-enhancer site drags a neighbour with it; the selective site does not.
        neighbors = {
            "SIM_NEIGHBOR_PROXIMAL": float(lfc * (0.80 if dist_shr <= 12 else 0.05)),
            "SIM_NEIGHBOR_DISTAL": float(lfc * 0.02),
        }

        return VariantMultiModalScore(
            variant=variant,
            condition_name=condition_name,
            target_gene_lfc=float(lfc),
            cpg_methylation_mean_abs_shift=cpg,
            snmc_seq3_max_abs_shift=cpg * 1.4,
            max_splice_junction_delta=splice,
            neighbor_gene_lfcs=neighbors,
        )
