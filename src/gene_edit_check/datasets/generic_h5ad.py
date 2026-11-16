"""Universal Bring-Your-Own `.h5ad` Cohort Loader (`datasets/generic_h5ad.py`).

Loads any single-cell AnnData file (local path or HTTP/HTTPS URL, such as a
CZ CELLxGENE Discover `.h5ad`) into the donor-balanced `RealCohort` contract
required by AlphaGenome 2 (`encode_cells` + `cell_conditioning`):

1. Resolves raw integer UMI counts from `.raw.X`, `.layers["counts"]`,
   `.layers["raw_counts"]`, or `.X`.
2. Maps gene identifiers to clean Ensembl gene IDs (`ENSG...`, version suffix
   stripped) from `var_names` or `var` metadata columns (`feature_id`,
   `ensembl_id`, `gene_id`, `gene_ids`).
3. Auto-detects donor/patient batch metadata (`donor_id`, `patient_id`, etc.).
4. Evaluates pandas `.obs.query(...)` expressions for the target and healthy
   cell states, computes full-population `observed_cp10k` on marker genes,
   draws donor-balanced pools, and caches processed cohorts in `data/processed/`.
"""

from __future__ import annotations

import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import urlretrieve

import anndata
import numpy as np
import pandas as pd
import scipy.sparse as sp

from gene_edit_check.datasets.crc_pelka import (
    PROCESSED_DIR,
    RAW_PATH,
    RealCohort,
    _donor_balanced_sample,
)

RAW_DIR = RAW_PATH.parent
_LAST_RAW_H5AD: tuple[str, anndata.AnnData] | None = None
DONOR_CANDIDATE_COLS: tuple[str, ...] = (
    "donor_id",
    "patient_id",
    "patient",
    "sample_id",
    "donor",
    "biosample_id",
    "individual",
)
ENSEMBL_VAR_COLS: tuple[str, ...] = (
    "ensembl_id",
    "feature_id",
    "gene_id",
    "gene_ids",
)


@dataclass(frozen=True)
class GenericH5adSpec:
    """Declarative specification for loading any `.h5ad` cohort."""

    h5ad_path_or_url: str | Path
    dataset_name: str
    target_state_name: str
    healthy_state_name: str
    target_obs_query: str
    healthy_obs_query: str
    donor_col: str | None = None
    ensembl_var_col: str | None = None
    citation: str = ""
    raw_filename: str | None = None


def resolve_h5ad_path(h5ad_path_or_url: str | Path, raw_filename: str | None = None) -> Path:
    """Return a local Path, downloading HTTP(S) assets into `data/raw/` if absent."""
    raw_str = str(h5ad_path_or_url)
    if raw_str.startswith(("http://", "https://")):
        RAW_DIR.mkdir(parents=True, exist_ok=True)
        fname = raw_filename or Path(urlparse(raw_str).path).name or "external_cohort.h5ad"
        if not fname.endswith(".h5ad"):
            fname = f"{fname}.h5ad"
        dest = RAW_DIR / fname
        if dest.exists() and dest.stat().st_size > 1_000_000:
            return dest
        urlretrieve(raw_str, dest)
        return dest

    path = Path(raw_str)
    if not path.exists():
        raise FileNotFoundError(f"Input .h5ad file not found: {path}")
    return path


def _detect_donor_series(obs: pd.DataFrame, donor_col: str | None) -> tuple[np.ndarray, str]:
    if donor_col is not None:
        if donor_col not in obs.columns:
            raise KeyError(
                f"Requested donor_col '{donor_col}' not found in .obs. "
                f"Available columns: {list(obs.columns)}"
            )
        return obs[donor_col].astype(str).to_numpy(), donor_col

    for cand in DONOR_CANDIDATE_COLS:
        if cand in obs.columns:
            return obs[cand].astype(str).to_numpy(), cand

    return np.asarray(["single_batch"] * len(obs), dtype=object), "fallback_single_batch"


def _is_integer_matrix(x: sp.spmatrix) -> bool:
    if x.nnz == 0:
        return True
    sample = x.data[: min(500_000, len(x.data))]
    return bool(np.allclose(sample, np.round(sample)) and np.all(sample >= 0))


def _extract_raw_counts(src: anndata.AnnData) -> tuple[sp.csr_matrix, pd.DataFrame, str]:
    """Select integer UMI count matrix and matching `.var` DataFrame."""
    candidates: list[tuple[str, object, pd.DataFrame]] = []
    if src.raw is not None and src.raw.X is not None:
        candidates.append(("raw.X", src.raw.X, src.raw.var))
    for layer_name in ("counts", "raw_counts", "UMI"):
        if layer_name in src.layers:
            candidates.append((f"layers[{layer_name}]", src.layers[layer_name], src.var))
    if src.X is not None:
        candidates.append(("X", src.X, src.var))

    for source_name, mat, var_df in candidates:
        csr = mat.tocsr() if sp.issparse(mat) else sp.csr_matrix(mat)
        if _is_integer_matrix(csr):
            return csr, var_df, source_name

    raise ValueError(
        "Could not locate non-negative integer UMI counts in .raw.X, "
        ".layers['counts'/'raw_counts'], or .X. AlphaGenome 2 encode_cells() requires raw integer counts."
    )


def _resolve_ensembl_ids(
    var_df: pd.DataFrame,
    fallback_var: pd.DataFrame,
    ensembl_var_col: str | None,
) -> tuple[np.ndarray, np.ndarray, dict[str, int]]:
    """Extract clean `ENSG...` IDs (no version suffix) and return (keep_col_mask, ensembl_ids, sym_to_col)."""
    raw_ids: pd.Series | None = None
    if ensembl_var_col is not None:
        if ensembl_var_col in var_df.columns:
            raw_ids = var_df[ensembl_var_col].astype(str)
        elif ensembl_var_col in fallback_var.columns and len(fallback_var) == len(var_df):
            raw_ids = fallback_var[ensembl_var_col].astype(str)
        else:
            raise KeyError(
                f"ensembl_var_col '{ensembl_var_col}' not found in .var columns: {list(var_df.columns)}"
            )
    else:
        idx_str = pd.Series(var_df.index.astype(str), index=var_df.index)
        if idx_str.str.startswith("ENSG").mean() > 0.25:
            raw_ids = idx_str
        else:
            for col in ENSEMBL_VAR_COLS:
                if col in var_df.columns and var_df[col].astype(str).str.startswith("ENSG").mean() > 0.25:
                    raw_ids = var_df[col].astype(str)
                    break
                if (
                    col in fallback_var.columns
                    and len(fallback_var) == len(var_df)
                    and fallback_var[col].astype(str).str.startswith("ENSG").mean() > 0.25
                ):
                    raw_ids = fallback_var[col].astype(str)
                    break

    if raw_ids is None:
        raise ValueError(
            "Could not find Ensembl gene IDs ('ENSG...') in .var_names or .var "
            f"columns {ENSEMBL_VAR_COLS}. Pass --ensembl-col explicitly."
        )

    cleaned = raw_ids.str.split(".").str[0].to_numpy()
    seen: set[str] = set()
    keep_mask = np.zeros(len(cleaned), dtype=bool)
    kept_orig_indices: list[int] = []
    for i, gid in enumerate(cleaned):
        if gid.startswith("ENSG") and gid not in seen:
            seen.add(gid)
            keep_mask[i] = True
            kept_orig_indices.append(i)

    kept_ids = cleaned[keep_mask]
    sym_to_col: dict[str, int] = {gid: idx for idx, gid in enumerate(kept_ids)}
    for sym_col in ("feature_name", "gene_symbols", "gene_symbol", "symbol"):
        src_df = var_df if sym_col in var_df.columns else (
            fallback_var if (sym_col in fallback_var.columns and len(fallback_var) == len(var_df)) else None
        )
        if src_df is not None:
            vals = src_df[sym_col].astype(str).to_numpy()[keep_mask]
            for idx, sym in enumerate(vals):
                if sym and sym != "nan" and sym not in sym_to_col:
                    sym_to_col[sym] = idx
    return keep_mask, kept_ids, sym_to_col


def _query_mask(
    obs: pd.DataFrame,
    raw_x: sp.csr_matrix,
    sym_to_col: dict[str, int],
    query_expr: str,
    state_name: str,
) -> np.ndarray:
    """Evaluate a pandas `.obs.query(...)` expression, optionally injecting `expr_<GENE>` columns."""
    if not query_expr.strip():
        raise ValueError(f"Empty obs query for state '{state_name}'.")
    expr_tokens = set(re.findall(r"\bexpr_([A-Za-z0-9_.-]+)\b", query_expr))
    obs_eval = obs
    if expr_tokens:
        obs_eval = obs.copy()
        for tok in expr_tokens:
            if tok not in sym_to_col:
                raise KeyError(
                    f"Expression token 'expr_{tok}' in query for '{state_name}' did not match any Ensembl ID or gene symbol."
                )
            col_i = sym_to_col[tok]
            obs_eval[f"expr_{tok}"] = np.asarray(raw_x[:, col_i].toarray()).ravel()
    selected_idx = obs_eval.query(query_expr).index
    mask = obs_eval.index.isin(selected_idx)
    if int(mask.sum()) == 0:
        raise ValueError(
            f"obs query for state '{state_name}' matched 0 cells: {query_expr!r}"
        )
    return np.asarray(mask, dtype=bool)


def build_generic_cohort(
    spec: GenericH5adSpec,
    pool_size: int = 600,
    seed: int = 42,
    marker_gene_ids: tuple[str, ...] = (),
    use_cache: bool = True,
) -> RealCohort:
    """Load any `.h5ad` file or URL into a donor-balanced `RealCohort`."""
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    spec_hash = hashlib.sha256(
        "|".join(
            [
                str(spec.h5ad_path_or_url),
                spec.target_state_name,
                spec.healthy_state_name,
                spec.target_obs_query,
                spec.healthy_obs_query,
                str(spec.donor_col or ""),
                str(spec.ensembl_var_col or ""),
            ]
        ).encode("utf-8")
    ).hexdigest()[:12]
    safe_name = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in spec.dataset_name.lower())[:32]
    tag = f"generic_{safe_name}_{spec_hash}_pool{pool_size}_seed{seed}"
    h5_path, js_path = PROCESSED_DIR / f"{tag}.h5ad", PROCESSED_DIR / f"{tag}.json"

    if use_cache and h5_path.exists() and js_path.exists():
        meta = json.loads(js_path.read_text(encoding="utf-8"))
        big = anndata.read_h5ad(h5_path)
        states = {
            s: big[big.obs["cell_state"] == s].copy()
            for s in (spec.target_state_name, spec.healthy_state_name)
        }
        return RealCohort(
            states=states,
            observed_cp10k=meta["observed_cp10k"],
            provenance=meta["provenance"],
        )

    local_path = resolve_h5ad_path(spec.h5ad_path_or_url, raw_filename=spec.raw_filename)
    resolved_key = str(local_path.resolve())
    global _LAST_RAW_H5AD
    if _LAST_RAW_H5AD is not None and _LAST_RAW_H5AD[0] == resolved_key:
        src = _LAST_RAW_H5AD[1]
    else:
        src = anndata.read_h5ad(local_path)
        _LAST_RAW_H5AD = (resolved_key, src)

    raw_x, var_df, count_source = _extract_raw_counts(src)
    keep_cols, gene_ids, sym_to_col = _resolve_ensembl_ids(var_df, src.var, spec.ensembl_var_col)
    if not np.all(keep_cols):
        raw_x = raw_x[:, keep_cols].tocsr()
    umi = np.maximum(1.0, np.asarray(raw_x.sum(axis=1)).ravel())

    donors_vec, resolved_donor_col = _detect_donor_series(src.obs, spec.donor_col)
    masks = {
        spec.target_state_name: _query_mask(
            src.obs, raw_x, sym_to_col, spec.target_obs_query, spec.target_state_name
        ),
        spec.healthy_state_name: _query_mask(
            src.obs, raw_x, sym_to_col, spec.healthy_obs_query, spec.healthy_state_name
        ),
    }

    rng = np.random.default_rng(seed)
    parts: list[anndata.AnnData] = []
    observed: dict[str, dict[str, float]] = {}
    prov_states: dict[str, object] = {}

    gene_to_col = {gid: i for i, gid in enumerate(gene_ids)}
    marker_cols = [gene_to_col[g] for g in marker_gene_ids if g in gene_to_col]
    marker_found_ids = [gene_ids[c] for c in marker_cols]

    for state, m in masks.items():
        idx_all = np.where(m)[0]
        donors_all = donors_vec[idx_all]
        if marker_cols:
            cp10k = raw_x[idx_all][:, marker_cols].multiply(1e4 / umi[idx_all][:, None]).tocsr()
            observed[state] = {
                marker_found_ids[k]: float(cp10k[:, k].mean()) for k in range(len(marker_cols))
            }
        else:
            observed[state] = {}

        pick = idx_all[_donor_balanced_sample(donors_all, pool_size, rng)]
        sub = anndata.AnnData(
            X=raw_x[pick].astype(np.int32),
            obs=pd.DataFrame(
                {
                    "donor_id": donors_vec[pick],
                    "cell_state": state,
                },
                index=src.obs_names[pick].astype(str),
            ),
            var=pd.DataFrame(index=pd.Index(gene_ids.astype(str))),
        )
        parts.append(sub)
        prov_states[state] = {
            "n_cells_full_population": int(len(idx_all)),
            "n_donors_full_population": int(len(np.unique(donors_all))),
            "n_cells_pool": int(sub.n_obs),
            "n_donors_pool": int(sub.obs["donor_id"].nunique()),
        }

    big = anndata.concat(parts, join="inner")
    big.var_names_make_unique()
    with local_path.open("rb") as fh:
        sha = hashlib.sha256(fh.read(8 * 1024 * 1024)).hexdigest()[:16]

    provenance: dict[str, object] = {
        "data_source": spec.dataset_name,
        "citation": spec.citation or f"Local/External .h5ad ({local_path.name})",
        "asset_url": str(spec.h5ad_path_or_url),
        "raw_file_head_sha256_16": sha,
        "count_matrix_source": count_source,
        "donor_column_used": resolved_donor_col,
        "target_definition": spec.target_obs_query,
        "healthy_definition": spec.healthy_obs_query,
        "pool_size_per_state": pool_size,
        "seed": seed,
        "states": prov_states,
    }
    big.write_h5ad(h5_path, compression="gzip")
    js_path.write_text(
        json.dumps({"observed_cp10k": observed, "provenance": provenance}, indent=2),
        encoding="utf-8",
    )
    states = {
        s: big[big.obs["cell_state"] == s].copy()
        for s in (spec.target_state_name, spec.healthy_state_name)
    }
    return RealCohort(states=states, observed_cp10k=observed, provenance=provenance)
