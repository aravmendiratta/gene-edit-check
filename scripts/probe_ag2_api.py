#!/usr/bin/env python3
"""Phase 0 Live API Probe for AlphaGenome 2 (`alphagenome==2.0.0`).

Executes minimal, low-quota live calls against `alphagenome.models.v2.dna_client`
to empirically verify:
1. `encode_cells()` output shape and latency.
2. Whether `score_variant(..., cell_conditioning=...)` accepts `MethylationScorer`
   and `SpliceJunctionScorer` alongside `SC_GENE_EXPRESSION_DENOISED`, or requires
   Tier 2 reference conditions (`ModelCondition.from_curies`).
3. `predict_interval_gradients(gene_id="ENSG00000188389", ...)` with `cell_conditioning`.
"""

from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path
from typing import Any

import anndata
import numpy as np
import pandas as pd
from dotenv import load_dotenv

load_dotenv(Path(__file__).resolve().parent.parent / ".env")

from alphagenome.data import genome
from alphagenome.models.v2 import dna_client, dna_output, model_condition, variant_scorers


def make_probe_anndata(n_cells: int = 8) -> anndata.AnnData:
    """Create a minimal valid human AnnData object with real immune ENSG IDs."""
    rng = np.random.default_rng(42)
    ensembl_ids = [
        "ENSG00000188389",  # PDCD1 (PD-1)
        "ENSG00000198846",  # TOX
        "ENSG00000123358",  # NR4A1
        "ENSG00000153563",  # CD8A
        "ENSG00000167286",  # CD3D
        "ENSG00000111537",  # IFNG
        "ENSG00000100453",  # GZMB
        "ENSG00000163599",  # CTLA4
        "ENSG00000089692",  # LAG3
        "ENSG00000135077",  # HAVCR2 (TIM-3)
        "ENSG00000159216",  # RUNX1
        "ENSG00000102908",  # NFATC2
    ]
    counts = rng.integers(1, 25, size=(n_cells, len(ensembl_ids)), dtype=np.int32)
    obs = pd.DataFrame(
        {"donor_id": [f"donor_{i % 2}" for i in range(n_cells)]},
        index=[f"cell_{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=ensembl_ids)
    return anndata.AnnData(X=counts, obs=obs, var=var)


def main() -> int:
    api_key = os.environ.get("ALPHAGENOME_API_KEY", "").strip()
    if not api_key:
        print(json.dumps({"error": "ALPHAGENOME_API_KEY missing in .env"}))
        return 1

    client = dna_client.create(api_key)
    adata = make_probe_anndata(n_cells=8)
    results: dict[str, Any] = {"alphagenome_version": "2.0.0", "probes": {}}

    # Probe 1: encode_cells()
    t0 = time.perf_counter()
    try:
        emb = client.encode_cells(adata, progress_bar=False)
        results["probes"]["1_encode_cells"] = {
            "status": "OK",
            "latency_sec": round(time.perf_counter() - t0, 3),
            "shape": list(emb.shape),
            "dtype": str(emb.dtype),
            "sample_norm": float(np.linalg.norm(emb[0])),
        }
    except Exception as exc:
        results["probes"]["1_encode_cells"] = {"status": "ERROR", "error": repr(exc)}

    # Interval of exact supported length SEQUENCE_LENGTH_100KB (131,072 bp) around PDCD1 on chr2
    # PDCD1 is at chr2:241,849,881-241,858,908 (GRCh38)
    center = 241_854_000
    half = dna_client.SEQUENCE_LENGTH_100KB // 2
    interval_100kb = genome.Interval(chromosome="chr2", start=center - half, end=center + half)
    variant = genome.Variant(
        chromosome="chr2", position=center, reference_bases="C", alternate_bases="T"
    )

    # Probe 2A: score_variant with cell_conditioning + SC_GENE_EXPRESSION_DENOISED only
    t0 = time.perf_counter()
    try:
        sc_expr_scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS["SC_GENE_EXPRESSION_DENOISED"]
        scores_sc = client.score_variant(
            interval=interval_100kb,
            variant=variant,
            variant_scorers=[sc_expr_scorer],
            cell_conditioning=adata,
            cell_condition_name="Exhausted_CD8_Probe",
        )
        results["probes"]["2a_score_variant_sc_expression"] = {
            "status": "OK",
            "latency_sec": round(time.perf_counter() - t0, 3),
            "num_tables": len(scores_sc),
            "obs_shape": list(scores_sc[0].shape) if scores_sc else [],
            "var_cols": list(scores_sc[0].var.columns) if scores_sc else [],
        }
    except Exception as exc:
        results["probes"]["2a_score_variant_sc_expression"] = {
            "status": "ERROR",
            "error": repr(exc),
        }

    # Probe 2B: Test if MethylationScorer / SpliceJunctionScorer work under cell_conditioning (with condition=None)
    t0 = time.perf_counter()
    try:
        hazard_scorers = [
            variant_scorers.MethylationScorer(width=1201),
            variant_scorers.SpliceJunctionScorer(),
        ]
        scores_haz_sc = client.score_variant(
            interval=interval_100kb,
            variant=variant,
            variant_scorers=hazard_scorers,
            cell_conditioning=adata,
            cell_condition_name="Exhausted_CD8_Probe",
        )
        results["probes"]["2b_hazard_scorers_under_cell_conditioning_only"] = {
            "status": "OK",
            "latency_sec": round(time.perf_counter() - t0, 3),
            "table_shapes": [list(t.shape) for t in scores_haz_sc],
        }
    except Exception as exc:
        results["probes"]["2b_hazard_scorers_under_cell_conditioning_only"] = {
            "status": "ERROR",
            "error": repr(exc),
        }

    # Probe 2C: Test Combined condition (CURIE) + cell_conditioning in a single call
    t0 = time.perf_counter()
    try:
        cond_curie = model_condition.from_curies(["UBERON:0000178", "CL:0000625"])  # blood CD8 T-cell
        combined_scorers = [
            variant_scorers.RECOMMENDED_VARIANT_SCORERS["SC_GENE_EXPRESSION_DENOISED"],
            variant_scorers.MethylationScorer(width=1201),
            variant_scorers.SpliceJunctionScorer(),
        ]
        scores_comb = client.score_variant(
            interval=interval_100kb,
            variant=variant,
            condition=cond_curie,
            variant_scorers=combined_scorers,
            cell_conditioning=adata,
            cell_condition_name="Exhausted_CD8_Probe",
        )
        results["probes"]["2c_combined_curie_and_cell_conditioning"] = {
            "status": "OK",
            "latency_sec": round(time.perf_counter() - t0, 3),
            "table_shapes": [list(t.shape) for t in scores_comb],
            "table_var_TrackNames": [list(t.var_names[:5]) for t in scores_comb],
        }
    except Exception as exc:
        results["probes"]["2c_combined_curie_and_cell_conditioning"] = {
            "status": "ERROR",
            "error": repr(exc),
        }

    # Probe 3: predict_interval_gradients with cell_conditioning (100KB)
    t0 = time.perf_counter()
    try:
        grads = client.predict_interval_gradients(
            gene_id="ENSG00000188389",  # PDCD1
            requested_output=dna_output.OutputType.SC_GENE_EXPRESSION_DENOISED,
            sequence_length=dna_client.SEQUENCE_LENGTH_100KB,
            cell_conditioning=adata,
            cell_condition_name="Exhausted_CD8_Probe",
        )
        grad_attrs = [a for a in dir(grads) if not a.startswith("_")]
        results["probes"]["3_predict_interval_gradients"] = {
            "status": "OK",
            "latency_sec": round(time.perf_counter() - t0, 3),
            "attributes": grad_attrs,
        }
    except Exception as exc:
        results["probes"]["3_predict_interval_gradients"] = {
            "status": "ERROR",
            "error": repr(exc),
        }

    out_path = Path(__file__).resolve().parent.parent / "phase0_probe_results.json"
    out_path.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(json.dumps(results, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
