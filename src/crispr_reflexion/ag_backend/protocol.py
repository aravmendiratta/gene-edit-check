"""VerifierBackend Protocol (`ag_backend/protocol.py`).

Decouples the optimization loop, coreset sketcher, and unit tests from whether
predictions come from the live AlphaGenome 2 gRPC service or the offline surrogate.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Protocol, runtime_checkable
import anndata
import numpy as np

from crispr_reflexion.models import BackendMode, GenomicInterval, SingleNucleotideVariant


@dataclass
class GradientProfile:
    gene_id: str
    gene_name: str
    condition_name: str
    interval: GenomicInterval
    reference_sequence: str
    signed_gradients_1bp: np.ndarray  # shape: (interval.width,)
    scalar_score: float


@dataclass
class VariantMultiModalScore:
    variant: SingleNucleotideVariant
    condition_name: str
    target_gene_lfc: float
    cpg_methylation_mean_abs_shift: float
    snmc_seq3_max_abs_shift: float
    max_splice_junction_delta: float
    # Every *other* gene the scorer returned for this interval (`gene_id -> LFC`), which the
    # expression scorer gives us for free in the same response. A base edit meant to silence one
    # gene must not quietly move its neighbours, and a 100-KB window around PDCD1 contains
    # several. Empty when the backend cannot resolve neighbours.
    neighbor_gene_lfcs: dict[str, float] = field(default_factory=dict)

    @property
    def max_abs_neighbor_lfc(self) -> float:
        return max((abs(v) for v in self.neighbor_gene_lfcs.values()), default=0.0)

    @property
    def worst_neighbor(self) -> tuple[str, float] | None:
        if not self.neighbor_gene_lfcs:
            return None
        gid = max(self.neighbor_gene_lfcs, key=lambda g: abs(self.neighbor_gene_lfcs[g]))
        return gid, self.neighbor_gene_lfcs[gid]


@runtime_checkable
class VerifierBackend(Protocol):
    @property
    def backend_mode(self) -> BackendMode: ...

    @property
    def model_version(self) -> str: ...

    def encode_cells(self, adata: anndata.AnnData, state_name: str) -> np.ndarray:
        """Return (n_cells, 256) float32 embeddings in AlphaGenome 2 latent space."""
        ...

    def compute_interval_gradients(
        self,
        gene_id: str,
        gene_symbol: str,
        cell_sketch: anndata.AnnData,
        condition_name: str,
        sequence_length: int = 131072,
    ) -> GradientProfile:
        """Compute 1-bp input gradients centered on `gene_id` TSS and recover GRCh38 sequence."""
        ...

    def score_ism_window(
        self,
        parent_interval: GenomicInterval,
        ism_window: GenomicInterval,
        gene_id: str,
        cell_sketch: anndata.AnnData,
        condition_name: str,
    ) -> dict[tuple[int, str, str], float]:
        """Run scoped ISM (`score_ism_variants`) on `ism_window` (<= 60 bp) and return
        map of `(1based_pos, ref_base, alt_base) -> target_gene_lfc`.
        """
        ...

    def score_variant_multimodal(
        self,
        parent_interval: GenomicInterval,
        variant: SingleNucleotideVariant,
        gene_id: str,
        cell_sketch: anndata.AnnData | None,
        condition_name: str,
        reference_curies: list[str] | None = None,
    ) -> VariantMultiModalScore:
        """Score a candidate SNV for target gene expression LFC, CpG methylation (`WGBS`/`snmC-seq3`),
        and splice junction shifts (`SpliceJunctionScorer`).
        """
        ...
