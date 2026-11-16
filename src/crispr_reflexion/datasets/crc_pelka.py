"""Real single-cell cohort: colorectal cancer T cells (Pelka et al., *Cell* 2021).

Source: CZ CELLxGENE Discover, collection `4d82bd8e-8827-410b-8fc7-685b7dd92585`
("Spatially organized multicellular immune hubs in human colorectal cancer"),
dataset "T, NK, ILC cells": 76,965 cells, 62 patients, tumor + matched normal colon,
10x 3' v2/v3, raw UMI counts in `.raw.X`, Ensembl IDs in `.raw.var_names`.

Cell states are the *authors'* annotations (`ClusterFull` / `ClusterMidway`), not ours:
- Target  "Exhausted_CD8_TIL":  tumor CD8 T cells in the three `CD8+ CXCL13+` clusters
  (cTNI14/15/16), the tumor-reactive / exhausted state with high PDCD1, LAG3, HAVCR2.
- Healthy "Normal_Colon_CD8":   CD8 T cells (`ClusterMidway == TCD8`) from adjacent
  normal colon of the same cohort.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from urllib.request import urlretrieve

import anndata
import numpy as np
import pandas as pd
import scipy.sparse as sp

from crispr_reflexion import paths

RAW_PATH = paths.data_dir() / "raw" / "crc_tnk.h5ad"
PROCESSED_DIR = paths.data_dir() / "processed"

COLLECTION_ID = "4d82bd8e-8827-410b-8fc7-685b7dd92585"
DATASET_ID = "7674f8f4-7200-4f08-8412-eb114016a66f"
ASSET_URL = "https://datasets.cellxgene.cziscience.com/8a3ea2bf-70e2-4b60-964d-e7b3f9ad75b3.h5ad"
CITATION = "Pelka K. et al. Spatially organized multicellular immune hubs in human colorectal cancer. Cell 184, 4734 (2021). doi:10.1016/j.cell.2021.08.003"

TARGET_STATE = "Exhausted_CD8_TIL"
HEALTHY_STATE = "Normal_Colon_CD8"
TARGET_CLUSTERS = (
    "cTNI14 (CD8+ CXCL13+)",
    "cTNI15 (CD8+ CXCL13+ HSP+)",
    "cTNI16 (CD8+ CXCL13+ prolif)",
)
TUMOR = "colon adenocarcinoma"
NORMAL = "normal"


@dataclass
class RealCohort:
    states: dict[str, anndata.AnnData]
    observed_cp10k: dict[str, dict[str, float]]  # state -> ENSG -> mean counts per 10k UMI (full population)
    provenance: dict[str, object]


def ensure_raw_download(path: Path = RAW_PATH) -> Path:
    """Download the pinned H5AD (396 MB) if it is not already on disk."""
    if path.exists() and path.stat().st_size > 1e8:
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading the Pelka et al. 2021 atlas (about 396 MB, one time) to {path} ...", flush=True)
    part = path.with_suffix(path.suffix + ".part")
    urlretrieve(ASSET_URL, part)
    part.replace(path)  # only a complete download ever carries the final name
    return path


def _donor_balanced_sample(donors: np.ndarray, n: int, rng: np.random.Generator) -> np.ndarray:
    """Round-robin over donors so no single patient dominates the pool."""
    if n >= len(donors):
        return np.arange(len(donors))
    per_donor = {d: list(rng.permutation(np.where(donors == d)[0])) for d in np.unique(donors)}
    order = list(per_donor)
    rng.shuffle(order)
    picked: list[int] = []
    while len(picked) < n and any(per_donor.values()):
        for d in order:
            if per_donor[d] and len(picked) < n:
                picked.append(per_donor[d].pop())
    return np.sort(np.array(picked, dtype=np.int64))


def _masks(obs: pd.DataFrame) -> dict[str, np.ndarray]:
    return {
        TARGET_STATE: (obs["ClusterFull"].isin(TARGET_CLUSTERS) & (obs["disease"] == TUMOR)).to_numpy(),
        HEALTHY_STATE: ((obs["ClusterMidway"] == "TCD8") & (obs["disease"] == NORMAL)).to_numpy(),
    }


def build_cohort(
    pool_size: int = 600,
    seed: int = 42,
    marker_gene_ids: tuple[str, ...] = (),
    use_cache: bool = True,
) -> RealCohort:
    """Build donor-balanced candidate pools (raw UMI counts) for the target and healthy states."""
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    tag = f"pelka_crc_cd8_pool{pool_size}_seed{seed}"
    h5_path, js_path = PROCESSED_DIR / f"{tag}.h5ad", PROCESSED_DIR / f"{tag}.json"
    if use_cache and h5_path.exists() and js_path.exists():
        meta = json.loads(js_path.read_text(encoding="utf-8"))
        big = anndata.read_h5ad(h5_path)
        states = {s: big[big.obs["cell_state"] == s].copy() for s in (TARGET_STATE, HEALTHY_STATE)}
        return RealCohort(states=states, observed_cp10k=meta["observed_cp10k"], provenance=meta["provenance"])

    src = anndata.read_h5ad(ensure_raw_download())
    raw_x = src.raw.X.tocsr() if sp.issparse(src.raw.X) else sp.csr_matrix(src.raw.X)
    gene_ids = np.asarray(src.raw.var_names)
    umi = np.asarray(raw_x.sum(axis=1)).ravel()
    if not np.allclose(raw_x.data, np.round(raw_x.data)):
        raise ValueError("Expected raw integer UMI counts in .raw.X")

    rng = np.random.default_rng(seed)
    masks = _masks(src.obs)
    parts: list[anndata.AnnData] = []
    observed: dict[str, dict[str, float]] = {}
    prov_states: dict[str, object] = {}
    marker_idx = [int(np.where(gene_ids == g)[0][0]) for g in marker_gene_ids if g in set(gene_ids)]
    for state, m in masks.items():
        idx_all = np.where(m)[0]
        donors_all = src.obs["donor_id"].astype(str).to_numpy()[idx_all]
        # Observed expression from the FULL state population (not just the pool): used for validation.
        cp10k = raw_x[idx_all][:, marker_idx].multiply(1e4 / umi[idx_all][:, None]).tocsr()
        observed[state] = {gene_ids[j]: float(cp10k[:, k].mean()) for k, j in enumerate(marker_idx)}
        pick = idx_all[_donor_balanced_sample(donors_all, pool_size, rng)]
        sub = anndata.AnnData(
            X=raw_x[pick].astype(np.int32),
            obs=pd.DataFrame(
                {
                    "donor_id": src.obs["donor_id"].astype(str).to_numpy()[pick],
                    "cell_state": state,
                    "author_cluster": src.obs["ClusterFull"].astype(str).to_numpy()[pick],
                    "MMRStatus": src.obs["MMRStatusTumor"].astype(str).to_numpy()[pick],
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
    sha = hashlib.sha256(RAW_PATH.read_bytes()[: 8 * 1024 * 1024]).hexdigest()[:16]
    provenance: dict[str, object] = {
        "data_source": "CZ CELLxGENE Discover (Pelka et al. 2021, colorectal cancer T/NK/ILC cells)",
        "citation": CITATION,
        "collection_id": COLLECTION_ID,
        "dataset_id": DATASET_ID,
        "asset_url": ASSET_URL,
        "raw_file_head_sha256_16": sha,
        "target_definition": f"disease == '{TUMOR}' & ClusterFull in {list(TARGET_CLUSTERS)}",
        "healthy_definition": f"disease == '{NORMAL}' & ClusterMidway == 'TCD8'",
        "pool_size_per_state": pool_size,
        "seed": seed,
        "states": prov_states,
    }
    big.write_h5ad(h5_path, compression="gzip")
    js_path.write_text(json.dumps({"observed_cp10k": observed, "provenance": provenance}, indent=2), encoding="utf-8")
    states = {s: big[big.obs["cell_state"] == s].copy() for s in (TARGET_STATE, HEALTHY_STATE)}
    return RealCohort(states=states, observed_cp10k=observed, provenance=provenance)
