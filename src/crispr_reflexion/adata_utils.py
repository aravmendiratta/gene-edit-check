"""Sparse-safe helpers shared by the live client, simulator and sketcher."""

from __future__ import annotations

import hashlib

import anndata
import numpy as np
import scipy.sparse as sp


def counts_fingerprint(adata: anndata.AnnData) -> str:
    """Content hash of the count matrix + gene IDs (dense or sparse, order-independent of storage)."""
    h = hashlib.sha256()
    x = adata.X
    # Canonicalise to CSR so dense and sparse copies of the same counts share a cache key.
    x = sp.csr_matrix(x) if not sp.issparse(x) else x.tocsr()
    x.sum_duplicates()
    x.eliminate_zeros()
    x.sort_indices()
    h.update(np.asarray(x.data, dtype=np.int32).tobytes())
    h.update(np.asarray(x.indices, dtype=np.int32).tobytes())
    h.update(np.asarray(x.indptr, dtype=np.int64).tobytes())
    h.update(str(adata.shape).encode("utf-8"))
    h.update("|".join(map(str, adata.var_names)).encode("utf-8"))
    return h.hexdigest()


def dense_counts(adata: anndata.AnnData, dtype: type = np.float32) -> np.ndarray:
    x = adata.X
    return (x.toarray() if sp.issparse(x) else np.asarray(x)).astype(dtype, copy=False)


def drop_unexpressed_genes(adata: anndata.AnnData) -> anndata.AnnData:
    """Remove genes with zero counts in every cell.

    Lossless for conditioning: library size and every non-zero count are unchanged, it only
    shrinks the request payload (e.g. 38k -> ~10k genes for a 128-cell sketch).
    """
    x = adata.X
    n_cells_expressing = np.asarray((x > 0).sum(axis=0)).ravel()
    return adata[:, n_cells_expressing > 0].copy()
