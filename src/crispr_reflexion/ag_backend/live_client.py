"""Live AlphaGenome 2 Backend (`ag_backend/live_client.py`).

Wraps `alphagenome.models.v2.dna_client.DnaClient` with:
- Content-addressed SQLite disk caching (`cache.py`)
- Strict API & ISM window budget enforcement (`budget.py`)
- Exact GRCh38 sequence recovery from `InputGradientsOutput.gradients` one-hot channels
"""

from __future__ import annotations

import anndata
import numpy as np

from alphagenome.data import genome
from alphagenome.models.v2 import dna_client, dna_output, model_condition, variant_scorers

from crispr_reflexion.ag_backend.protocol import (
    GradientProfile,
    VariantMultiModalScore,
    VerifierBackend,
)
from crispr_reflexion.adata_utils import counts_fingerprint
from crispr_reflexion.budget import RunBudgetTracker
from crispr_reflexion.cache import ContentAddressedCache
from crispr_reflexion.models import BackendMode, GenomicInterval, SingleNucleotideVariant

CHANNEL_TO_BASE = {0: "A", 1: "C", 2: "G", 3: "T"}


def _hash_anndata_counts(adata: anndata.AnnData) -> str:
    return counts_fingerprint(adata)


def recover_reference_sequence_from_gradients(gradients_4ch: np.ndarray) -> str:
    """Recover exact GRCh38 DNA sequence from `InputGradientsOutput.gradients` (shape `[L, 4]`).

    In AlphaGenome 2, `gradients` is `grad * one_hot_sequence`, so only the reference base
    channel (`0=A, 1=C, 2=G, 3=T`) is non-zero at each position (N-padded regions are all 0).
    """
    abs_g = np.abs(gradients_4ch)
    max_vals = np.max(abs_g, axis=-1)
    arg_max = np.argmax(abs_g, axis=-1)
    chars = [
        CHANNEL_TO_BASE[int(ch)] if val > 0.0 else "N"
        for ch, val in zip(arg_max, max_vals)
    ]
    return "".join(chars)


class LiveAlphaGenomeBackend(VerifierBackend):
    def __init__(
        self,
        api_key: str,
        budget: RunBudgetTracker,
        cache: ContentAddressedCache,
    ) -> None:
        self._client = dna_client.create(api_key)
        self._budget = budget
        self._cache = cache
        self._model_version = "alphagenome-2.0.0-all_folds"

    @property
    def backend_mode(self) -> BackendMode:
        return BackendMode.LIVE_ALPHAGENOME_V2

    @property
    def model_version(self) -> str:
        return self._model_version

    def encode_cells(self, adata: anndata.AnnData, state_name: str) -> np.ndarray:
        cond_hash = _hash_anndata_counts(adata)
        key = self._cache.make_key("encode_cells", str(adata.shape), cond_hash, "CELL_EMBEDDINGS")
        cached = self._cache.get(key)
        if cached is not None:
            self._budget.record_cache_hit()
            return np.asarray(cached["embeddings"], dtype=np.float32)

        self._budget.check_and_record_encode_cells(adata.shape[0])
        emb = self._client.encode_cells(adata, progress_bar=False)
        self._cache.put(key, "encode_cells", {"embeddings": emb.tolist()})
        return np.asarray(emb, dtype=np.float32)

    def compute_interval_gradients(
        self,
        gene_id: str,
        gene_symbol: str,
        cell_sketch: anndata.AnnData,
        condition_name: str,
        sequence_length: int = dna_client.SEQUENCE_LENGTH_100KB,
    ) -> GradientProfile:
        self._budget.check_conditioning_size(cell_sketch.shape[0])
        cond_hash = _hash_anndata_counts(cell_sketch)
        key = self._cache.make_key(
            "predict_interval_gradients",
            f"{gene_id}:{sequence_length}",
            cond_hash,
            "SC_GENE_EXPRESSION_DENOISED",
        )
        cached = self._cache.get(key)
        if cached is not None:
            self._budget.record_cache_hit()
            return GradientProfile(
                gene_id=cached["gene_id"],
                gene_name=cached["gene_name"],
                condition_name=condition_name,
                interval=GenomicInterval(**cached["interval"]),
                reference_sequence=cached["reference_sequence"],
                signed_gradients_1bp=np.asarray(cached["signed_gradients_1bp"], dtype=np.float32),
                scalar_score=float(cached["scalar_score"]),
            )

        self._budget.check_and_record_gradient_call(1)
        out = self._client.predict_interval_gradients(
            gene_id=gene_id,
            requested_output=dna_output.OutputType.SC_GENE_EXPRESSION_DENOISED,
            sequence_length=sequence_length,
            cell_conditioning=cell_sketch,
            cell_condition_name=condition_name,
        )
        ref_seq = recover_reference_sequence_from_gradients(out.gradients)
        signed_1bp = out.gradients.sum(axis=-1).astype(np.float32)
        iv = GenomicInterval(
            chromosome=out.interval.chromosome,
            start=int(out.interval.start),
            end=int(out.interval.end),
        )
        payload = {
            "gene_id": out.gene_id,
            "gene_name": out.gene_name or gene_symbol,
            "interval": iv.model_dump(),
            "reference_sequence": ref_seq,
            "signed_gradients_1bp": signed_1bp.tolist(),
            "scalar_score": float(out.scalar_score),
        }
        self._cache.put(key, "predict_interval_gradients", payload)
        return GradientProfile(
            gene_id=out.gene_id,
            gene_name=out.gene_name or gene_symbol,
            condition_name=condition_name,
            interval=iv,
            reference_sequence=ref_seq,
            signed_gradients_1bp=signed_1bp,
            scalar_score=float(out.scalar_score),
        )

    def score_ism_window(
        self,
        parent_interval: GenomicInterval,
        ism_window: GenomicInterval,
        gene_id: str,
        cell_sketch: anndata.AnnData,
        condition_name: str,
    ) -> dict[tuple[int, str, str], float]:
        self._budget.check_conditioning_size(cell_sketch.shape[0])
        window_bp = ism_window.width
        n_variants = window_bp * 3
        cond_hash = _hash_anndata_counts(cell_sketch)
        iv_str = f"{parent_interval.chromosome}:{parent_interval.start}-{parent_interval.end}|{ism_window.start}-{ism_window.end}"
        key = self._cache.make_key(
            "score_ism_variants",
            iv_str,
            cond_hash,
            f"SC_GENE_EXPRESSION_DENOISED:{gene_id}",
        )
        cached = self._cache.get(key)
        if cached is not None:
            self._budget.record_cache_hit()
            return {
                (int(item["pos"]), str(item["ref"]), str(item["alt"])): float(item["lfc"])
                for item in cached["ism_scores"]
            }

        self._budget.check_and_record_ism_window(window_bp=window_bp, n_variants=n_variants)
        ag_parent = genome.Interval(
            chromosome=parent_interval.chromosome,
            start=parent_interval.start,
            end=parent_interval.end,
        )
        ag_ism = genome.Interval(
            chromosome=ism_window.chromosome,
            start=ism_window.start,
            end=ism_window.end,
        )
        sc_scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS["SC_GENE_EXPRESSION_DENOISED"]
        ism_tables = self._client.score_ism_variants(
            interval=ag_parent,
            ism_interval=ag_ism,
            variant_scorers=[sc_scorer],
            cell_conditioning=cell_sketch,
            cell_condition_name=condition_name,
            progress_bar=False,
        )
        out_map: dict[tuple[int, str, str], float] = {}
        serialized_items = []
        clean_target = gene_id.split(".")[0]
        for variant_list in ism_tables:
            tbl = variant_list[0]
            var_obj = tbl.uns["variant"]
            # Find row matching target gene_id
            gene_col = tbl.obs["gene_id"].astype(str).str.split(".").str[0]
            match_rows = np.where(gene_col == clean_target)[0]
            row_idx = int(match_rows[0]) if len(match_rows) > 0 else 0
            lfc_val = float(np.asarray(tbl.X)[row_idx, 0])
            k_tuple = (int(var_obj.position), str(var_obj.reference_bases), str(var_obj.alternate_bases))
            out_map[k_tuple] = lfc_val
            serialized_items.append(
                {
                    "pos": int(var_obj.position),
                    "ref": str(var_obj.reference_bases),
                    "alt": str(var_obj.alternate_bases),
                    "lfc": lfc_val,
                }
            )
        self._cache.put(key, "score_ism_variants", {"ism_scores": serialized_items})
        return out_map

    def score_variant_multimodal(
        self,
        parent_interval: GenomicInterval,
        variant: SingleNucleotideVariant,
        gene_id: str,
        cell_sketch: anndata.AnnData | None,
        condition_name: str,
        reference_curies: list[str] | None = None,
    ) -> VariantMultiModalScore:
        if cell_sketch is not None:
            self._budget.check_conditioning_size(cell_sketch.shape[0])
            cond_hash = _hash_anndata_counts(cell_sketch)
            scorers = [
                variant_scorers.RECOMMENDED_VARIANT_SCORERS["SC_GENE_EXPRESSION_DENOISED"],
                variant_scorers.MethylationScorer(width=1201),
                variant_scorers.SpliceJunctionScorer(),
            ]
            req_str = f"SC_EXPR+METH+SPLICE:{gene_id}"
        else:
            curies = reference_curies or ["UBERON:0000178", "CL:0000625"]
            cond_hash = "_".join(curies)
            scorers = [
                variant_scorers.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"],
                variant_scorers.MethylationScorer(width=1201),
                variant_scorers.SpliceJunctionScorer(),
            ]
            req_str = f"RNA_SEQ+METH+SPLICE:{gene_id}"

        v_str = f"{parent_interval.chromosome}:{parent_interval.start}-{parent_interval.end}|{variant.position}:{variant.reference}>{variant.alternate}"
        key = self._cache.make_key("score_variant_multimodal", v_str, cond_hash, req_str)
        cached = self._cache.get(key)
        if cached is not None:
            self._budget.record_cache_hit()
            return VariantMultiModalScore(
                variant=variant,
                condition_name=condition_name,
                target_gene_lfc=float(cached["target_gene_lfc"]),
                cpg_methylation_mean_abs_shift=float(cached["cpg_methylation_mean_abs_shift"]),
                snmc_seq3_max_abs_shift=float(cached["snmc_seq3_max_abs_shift"]),
                max_splice_junction_delta=float(cached["max_splice_junction_delta"]),
                neighbor_gene_lfcs={
                    str(k): float(v) for k, v in (cached.get("neighbor_gene_lfcs") or {}).items()
                },
            )

        self._budget.check_and_record_variant_score(1)
        ag_iv = genome.Interval(
            chromosome=parent_interval.chromosome,
            start=parent_interval.start,
            end=parent_interval.end,
        )
        ag_var = genome.Variant(
            chromosome=variant.chromosome,
            position=variant.position,
            reference_bases=variant.reference,
            alternate_bases=variant.alternate,
        )

        if cell_sketch is not None:
            tables = self._client.score_variant(
                interval=ag_iv,
                variant=ag_var,
                variant_scorers=scorers,
                cell_conditioning=cell_sketch,
                cell_condition_name=condition_name,
            )
        else:
            cond_obj = model_condition.from_curies(reference_curies or ["UBERON:0000178", "CL:0000625"])
            tables = self._client.score_variant(
                interval=ag_iv,
                variant=ag_var,
                condition=cond_obj,
                variant_scorers=scorers,
            )

        expr_tbl, meth_tbl, splice_tbl = tables[0], tables[1], tables[2]
        clean_target = gene_id.split(".")[0]
        neighbor_lfcs: dict[str, float] = {}
        if "gene_id" in expr_tbl.obs.columns and expr_tbl.shape[0] > 0:
            gene_col = expr_tbl.obs["gene_id"].astype(str).str.split(".").str[0]
            matches = np.where(gene_col == clean_target)[0]
            row_idx = int(matches[0]) if len(matches) > 0 else 0
            x = np.asarray(expr_tbl.X)
            target_lfc = float(np.mean(x[row_idx]))
            # The scorer already returned every gene in this interval; keeping the non-target rows
            # costs no extra API call and gives us a neighbour-gene specificity check for free.
            symbols = (
                expr_tbl.obs["gene_name"].astype(str).tolist()
                if "gene_name" in expr_tbl.obs.columns
                else [""] * len(gene_col)
            )
            for i, gid in enumerate(gene_col.tolist()):
                if gid == clean_target:
                    continue
                label = symbols[i] if symbols[i] and symbols[i] != "nan" else gid
                neighbor_lfcs[label] = float(np.mean(x[i]))
        else:
            target_lfc = 0.0

        meth_arr = np.abs(np.asarray(meth_tbl.X)) if meth_tbl.shape[0] > 0 else np.zeros((1, 2))
        cpg_mean = float(np.mean(meth_arr))
        snmc_max = float(np.max(meth_arr[:, -1])) if meth_arr.shape[1] >= 2 else float(np.max(meth_arr))

        splice_arr = np.abs(np.asarray(splice_tbl.X)) if splice_tbl.shape[0] > 0 else np.zeros((1, 1))
        max_splice = float(np.max(splice_arr))

        payload = {
            "target_gene_lfc": target_lfc,
            "cpg_methylation_mean_abs_shift": cpg_mean,
            "snmc_seq3_max_abs_shift": snmc_max,
            "max_splice_junction_delta": max_splice,
            "neighbor_gene_lfcs": neighbor_lfcs,
        }
        self._cache.put(key, "score_variant_multimodal", payload)
        return VariantMultiModalScore(
            variant=variant,
            condition_name=condition_name,
            target_gene_lfc=target_lfc,
            cpg_methylation_mean_abs_shift=cpg_mean,
            snmc_seq3_max_abs_shift=snmc_max,
            max_splice_junction_delta=max_splice,
            neighbor_gene_lfcs=neighbor_lfcs,
        )
