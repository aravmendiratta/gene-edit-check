"""Single-Cell Cohort & Scenario Registry (`scenarios/registry.py`).

Constructs multi-donor human `AnnData` single-cell cohorts with real Ensembl gene IDs
(`ENSG...`) and realistic donor batch effects so `encode_cells()` and donor-stratified
coreset sketching can be evaluated honestly.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING
import anndata
import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from gene_edit_check.datasets.generic_h5ad import GenericH5adSpec


# 24 validated human Ensembl IDs spanning T-cell exhaustion, naive/memory programs, and core TFs
IMMUNE_ENSEMBL_PANEL: tuple[tuple[str, str], ...] = (
    ("ENSG00000188389", "PDCD1"),    # PD-1 (Target)
    ("ENSG00000198846", "TOX"),      # Exhaustion master regulator
    ("ENSG00000123358", "NR4A1"),    # Exhaustion nuclear receptor
    ("ENSG00000089692", "LAG3"),     # Exhaustion checkpoint
    ("ENSG00000135077", "HAVCR2"),   # TIM-3
    ("ENSG00000163599", "CTLA4"),    # Checkpoint
    ("ENSG00000181847", "TIGIT"),    # Checkpoint
    ("ENSG00000138185", "ENTPD1"),   # CD39
    ("ENSG00000081059", "TCF7"),     # TCF-1 (Naive/stemness)
    ("ENSG00000126353", "CCR7"),     # Naive homing
    ("ENSG00000168685", "IL7R"),     # Naive survival
    ("ENSG00000188404", "SELL"),     # CD62L
    ("ENSG00000138795", "LEF1"),     # Naive TF
    ("ENSG00000153563", "CD8A"),     # Lineage
    ("ENSG00000167286", "CD3D"),     # Lineage
    ("ENSG00000159216", "RUNX1"),    # Shared core TF
    ("ENSG00000102908", "NFATC2"),   # TCR signaling TF
    ("ENSG00000111537", "IFNG"),     # Effector cytokine
    ("ENSG00000100453", "GZMB"),     # Cytolytic granule
    ("ENSG00000180644", "PRF1"),     # Perforin
    ("ENSG00000168610", "STAT3"),    # Signaling TF
    ("ENSG00000170345", "FOS"),      # AP-1 subunit
    ("ENSG00000177606", "JUN"),      # AP-1 subunit
    ("ENSG00000075624", "ACTB"),     # Housekeeping
)


@dataclass(frozen=True)
class ScenarioSpec:
    scenario_id: str
    title: str
    gene_id: str
    gene_symbol: str
    target_state_name: str
    healthy_state_names: tuple[str, ...]
    # Primary Tier 2 hazard reference: chosen for *evidence strength*, not anatomy. See
    # `conditions.py` -- if left empty, resolved automatically via `resolve_reference_pair()`.
    reference_curies: tuple[str, ...] = ()
    # Anatomically matched but evidentially weaker condition, used as a concordance check.
    tissue_matched_curies: tuple[str, ...] = ()
    data_source: str = "synthetic"  # "real" -> Pelka CRC; "generic" -> any .h5ad; "synthetic" -> demo/tests
    desired_direction: str = "down"
    generic_h5ad_spec: "GenericH5adSpec | None" = None
    cell_type_query: str = "CD8-positive"
    tissue_query: str | None = None


SCENARIOS: dict[str, ScenarioSpec] = {
    "crc-pd1-exhaustion": ScenarioSpec(
        scenario_id="crc-pd1-exhaustion",
        title=(
            "PDCD1 (PD-1) silencing in exhausted tumor-infiltrating CD8+ T cells vs. normal-colon CD8+ T cells "
            "(REAL DATA: Pelka et al. 2021, colorectal cancer, 62 patients)"
        ),
        gene_id="ENSG00000188389",
        gene_symbol="PDCD1",
        target_state_name="Exhausted_CD8_TIL",
        healthy_state_names=("Normal_Colon_CD8",),
        reference_curies=("UBERON:0000178", "CL:0000625"),
        tissue_matched_curies=("UBERON:0001155", "CL:0000625"),  # colon CD8+ alpha-beta T cell
        data_source="real",
    ),
    "pd1-exhaustion": ScenarioSpec(
        scenario_id="pd1-exhaustion",
        title="[SYNTHETIC DEMO DATA - not biology] PDCD1 in generated exhausted vs naive CD8 cells",
        gene_id="ENSG00000188389",
        gene_symbol="PDCD1",
        target_state_name="Exhausted_CD8_TIL",
        healthy_state_names=("Healthy_Naive_CD8",),
        reference_curies=("UBERON:0000178", "CL:0000625"),
    ),
}


# Additional checkpoint/exhaustion genes on the *same* cells, so the efficacy-vs-methylation
# coupling found at PDCD1 can be tested for locus specificity rather than assumed general.
# Each is a real immunotherapy target whose silencing in exhausted T cells is a plausible
# hypothesis, and each uses the identical cohort, sketches, thresholds and references.
COUPLING_SWEEP_GENES: tuple[tuple[str, str], ...] = (
    ("ENSG00000188389", "PDCD1"),   # PD-1   -- the primary scenario
    ("ENSG00000163599", "CTLA4"),   # CTLA-4 -- ipilimumab target
    ("ENSG00000135077", "HAVCR2"),  # TIM-3
    ("ENSG00000089692", "LAG3"),    # LAG-3  -- relatlimab target
    ("ENSG00000181847", "TIGIT"),   # TIGIT
    ("ENSG00000198846", "TOX"),     # exhaustion master transcription factor
)

# PDCD1 already has a scenario; the rest are registered from the table above.
SWEEP_SCENARIO_IDS: dict[str, str] = {"PDCD1": "crc-pd1-exhaustion"}

for _gid, _sym in COUPLING_SWEEP_GENES:
    _sid = SWEEP_SCENARIO_IDS.get(_sym, f"crc-{_sym.lower()}-exhaustion")
    SWEEP_SCENARIO_IDS[_sym] = _sid
    if _sid in SCENARIOS:
        continue
    SCENARIOS[_sid] = ScenarioSpec(
        scenario_id=_sid,
        title=(
            f"{_sym} silencing in exhausted tumor-infiltrating CD8+ T cells vs. normal-colon CD8+ "
            f"T cells (REAL DATA: Pelka et al. 2021, colorectal cancer, 62 patients)"
        ),
        gene_id=_gid,
        gene_symbol=_sym,
        target_state_name="Exhausted_CD8_TIL",
        healthy_state_names=("Normal_Colon_CD8",),
        reference_curies=("UBERON:0000178", "CL:0000625"),
        tissue_matched_curies=("UBERON:0001155", "CL:0000625"),
        data_source="real",
    )


def _register_agnostic_cross_context_scenarios() -> None:
    """Register cross-subtype, cross-lineage, and independent-cancer PDCD1 scenarios."""
    from gene_edit_check.datasets import crc_pelka
    from gene_edit_check.datasets.generic_h5ad import GenericH5adSpec

    tnbc_url = "https://datasets.cellxgene.cziscience.com/4d036b2c-6135-4515-9c44-a5aa7a2e32e5.h5ad"
    tnbc_citation = (
        "Wu S.Z. et al. Stromal cell diversity associated with immune evasion in human "
        "triple-negative breast cancer. EMBO J 39, e104063 (2020). doi:10.15252/embj.2019104063"
    )
    exh_clusters = "['cTNI14 (CD8+ CXCL13+)', 'cTNI15 (CD8+ CXCL13+ HSP+)', 'cTNI16 (CD8+ CXCL13+ prolif)']"

    entries: list[ScenarioSpec] = [
        ScenarioSpec(
            scenario_id="crc-pd1-mmrd",
            title="PDCD1 in MMRd (MSI-H) colorectal cancer exhausted CD8+ TILs vs normal colon CD8+ (Pelka et al. 2021)",
            gene_id="ENSG00000188389",
            gene_symbol="PDCD1",
            target_state_name="Exhausted_CD8_MMRd",
            healthy_state_names=("Normal_Colon_CD8",),
            data_source="generic",
            cell_type_query="CD8-positive",
            tissue_query="colon",
            generic_h5ad_spec=GenericH5adSpec(
                h5ad_path_or_url=crc_pelka.ASSET_URL,
                raw_filename="crc_tnk.h5ad",
                dataset_name="Pelka 2021 CRC (MMRd / MSI-H subtype)",
                target_state_name="Exhausted_CD8_MMRd",
                healthy_state_name="Normal_Colon_CD8",
                target_obs_query=f"disease == 'colon adenocarcinoma' and MMRStatusTumor == 'MMRd' and ClusterFull in {exh_clusters}",
                healthy_obs_query="disease == 'normal' and ClusterMidway == 'TCD8'",
                donor_col="donor_id",
                citation=crc_pelka.CITATION,
            ),
        ),
        ScenarioSpec(
            scenario_id="crc-pd1-mmrp",
            title="PDCD1 in MMRp (MSS) colorectal cancer exhausted CD8+ TILs vs normal colon CD8+ (Pelka et al. 2021)",
            gene_id="ENSG00000188389",
            gene_symbol="PDCD1",
            target_state_name="Exhausted_CD8_MMRp",
            healthy_state_names=("Normal_Colon_CD8",),
            data_source="generic",
            cell_type_query="CD8-positive",
            tissue_query="colon",
            generic_h5ad_spec=GenericH5adSpec(
                h5ad_path_or_url=crc_pelka.ASSET_URL,
                raw_filename="crc_tnk.h5ad",
                dataset_name="Pelka 2021 CRC (MMRp / MSS subtype)",
                target_state_name="Exhausted_CD8_MMRp",
                healthy_state_name="Normal_Colon_CD8",
                target_obs_query=f"disease == 'colon adenocarcinoma' and MMRStatusTumor == 'MMRp' and ClusterFull in {exh_clusters}",
                healthy_obs_query="disease == 'normal' and ClusterMidway == 'TCD8'",
                donor_col="donor_id",
                citation=crc_pelka.CITATION,
            ),
        ),
        ScenarioSpec(
            scenario_id="crc-pd1-cd4-cxcl13",
            title="PDCD1 in exhausted CD4+ CXCL13+ helper TILs vs normal colon CD4+ T cells (Pelka et al. 2021)",
            gene_id="ENSG00000188389",
            gene_symbol="PDCD1",
            target_state_name="Exhausted_CD4_CXCL13",
            healthy_state_names=("Normal_Colon_CD4",),
            data_source="generic",
            cell_type_query="CD4-positive, alpha-beta T cell",
            tissue_query="colon",
            generic_h5ad_spec=GenericH5adSpec(
                h5ad_path_or_url=crc_pelka.ASSET_URL,
                raw_filename="crc_tnk.h5ad",
                dataset_name="Pelka 2021 CRC (CD4+ CXCL13+ lineage)",
                target_state_name="Exhausted_CD4_CXCL13",
                healthy_state_name="Normal_Colon_CD4",
                target_obs_query="disease == 'colon adenocarcinoma' and ClusterFull == 'cTNI07 (CD4+ CXCL13+)'",
                healthy_obs_query="disease == 'normal' and ClusterMidway == 'TCD4'",
                donor_col="donor_id",
                citation=crc_pelka.CITATION,
            ),
        ),
        ScenarioSpec(
            scenario_id="crc-pd1-gammadelta",
            title="PDCD1 in intratumoral PDCD1+ gamma-delta T cells vs normal colon gamma-delta T cells (Pelka et al. 2021)",
            gene_id="ENSG00000188389",
            gene_symbol="PDCD1",
            target_state_name="Tumor_PDCD1_GammaDelta_T",
            healthy_state_names=("Normal_Colon_GammaDelta_T",),
            data_source="generic",
            cell_type_query="gamma-delta",
            tissue_query="colon",
            generic_h5ad_spec=GenericH5adSpec(
                h5ad_path_or_url=crc_pelka.ASSET_URL,
                raw_filename="crc_tnk.h5ad",
                dataset_name="Pelka 2021 CRC (gamma-delta T lineage)",
                target_state_name="Tumor_PDCD1_GammaDelta_T",
                healthy_state_name="Normal_Colon_GammaDelta_T",
                target_obs_query="disease == 'colon adenocarcinoma' and ClusterFull == 'cTNI18 (gd-like T PDCD1+)'",
                healthy_obs_query="disease == 'normal' and ClusterMidway == 'Tgd'",
                donor_col="donor_id",
                citation=crc_pelka.CITATION,
            ),
        ),
        ScenarioSpec(
            scenario_id="tnbc-pd1-exhaustion",
            title="PDCD1 in triple-negative breast cancer exhausted/LAG3+/HAVCR2+/CXCL13+ CD8+ TILs vs IL7R+/CCR7+ CD8+ (Wu et al. 2020)",
            gene_id="ENSG00000188389",
            gene_symbol="PDCD1",
            target_state_name="Exhausted_CD8_TNBC",
            healthy_state_names=("Memory_Naive_CD8_TNBC",),
            data_source="generic",
            cell_type_query="CD8-positive",
            tissue_query="breast",
            generic_h5ad_spec=GenericH5adSpec(
                h5ad_path_or_url=tnbc_url,
                raw_filename="tnbc_wu2020.h5ad",
                dataset_name="Wu et al. 2020 Triple-Negative Breast Cancer (CELLxGENE)",
                target_state_name="Exhausted_CD8_TNBC",
                healthy_state_name="Memory_Naive_CD8_TNBC",
                target_obs_query=(
                    "celltype_final == 'CD8+ T-cells' and "
                    "(expr_ENSG00000089692 > 0 or expr_ENSG00000135077 > 0 or expr_ENSG00000156234 > 0)"
                ),
                healthy_obs_query=(
                    "celltype_final == 'CD8+ T-cells' and "
                    "expr_ENSG00000089692 == 0 and expr_ENSG00000135077 == 0 and expr_ENSG00000156234 == 0 and "
                    "(expr_ENSG00000168685 > 0 or expr_ENSG00000126353 > 0)"
                ),
                donor_col="donor_id",
                citation=tnbc_citation,
            ),
        ),
    ]
    for s in entries:
        SCENARIOS[s.scenario_id] = s


UNIVERSALITY_CONTEXT_LABELS: dict[str, str] = {
    "crc-pd1-exhaustion": "CRC All CD8+ TILs (Pelka 2021)",
    "crc-pd1-mmrd": "CRC MMRd / MSI-H CD8+ TILs",
    "crc-pd1-mmrp": "CRC MMRp / MSS CD8+ TILs",
    "crc-pd1-cd4-cxcl13": "CRC CD4+ CXCL13+ Helper TILs",
    "crc-pd1-gammadelta": "CRC PDCD1+ gamma-delta TILs",
    "tnbc-pd1-exhaustion": "Breast Cancer (TNBC) CD8+ TILs (Wu 2020)",
}


_register_agnostic_cross_context_scenarios()

def build_til_cohort_anndata(
    n_cells_per_state: int = 120,
    seed: int = 42,
) -> dict[str, anndata.AnnData]:
    """Generate a multi-donor single-cell UMI count cohort with realistic imbalanced
    donor representation (60%, 20%, 12%, 8%) and donor-specific technical scaling.
    """
    rng = np.random.default_rng(seed)
    ensembl_ids = [g[0] for g in IMMUNE_ENSEMBL_PANEL]
    gene_symbols = [g[1] for g in IMMUNE_ENSEMBL_PANEL]
    n_genes = len(ensembl_ids)

    # Mean transcriptional rates for Exhausted CD8 TILs vs. Healthy Naive CD8 T-Cells
    exhausted_rates = np.array(
        [45, 38, 32, 40, 35, 30, 28, 34, 2, 1, 2, 1, 2, 50, 48, 22, 26, 12, 28, 24, 20, 25, 25, 95],
        dtype=np.float64,
    )
    naive_rates = np.array(
        [1, 1, 2, 1, 1, 1, 1, 1, 42, 46, 52, 44, 38, 50, 48, 24, 18, 2, 2, 2, 14, 6, 6, 95],
        dtype=np.float64,
    )

    donor_probs = [0.60, 0.20, 0.12, 0.08]
    donor_names = ["donor_1", "donor_2", "donor_3", "donor_4"]
    donor_tech_multipliers = {
        "donor_1": 1.0,
        "donor_2": 1.45,
        "donor_3": 0.68,
        "donor_4": 1.85,
    }

    out: dict[str, anndata.AnnData] = {}
    for state_name, base_rates in (
        ("Exhausted_CD8_TIL", exhausted_rates),
        ("Healthy_Naive_CD8", naive_rates),
    ):
        donors = rng.choice(donor_names, size=n_cells_per_state, p=donor_probs)
        counts = np.zeros((n_cells_per_state, n_genes), dtype=np.int32)
        for i in range(n_cells_per_state):
            d = str(donors[i])
            mult = donor_tech_multipliers[d]
            # Distinct rare transcriptional sub-state (15% of cells) that random subsampling under-covers
            if i % 7 == 0:
                sub_vec = np.where(base_rates > 20, 3.2, 0.25)
            elif i % 3 == 0:
                sub_vec = np.where(base_rates > 20, 0.55, 1.8)
            else:
                sub_vec = np.ones_like(base_rates)
            lam = np.maximum(1.0, base_rates * mult * sub_vec)
            row = rng.poisson(lam).astype(np.int32)
            # Ensure non-zero counts per gene_counts_condition contract
            row = np.maximum(1, row)
            counts[i] = row

        obs = pd.DataFrame(
            {
                "donor_id": donors,
                "cell_state": [state_name] * n_cells_per_state,
            },
            index=[f"{state_name}_cell_{i}" for i in range(n_cells_per_state)],
        )
        var = pd.DataFrame({"gene_symbol": gene_symbols}, index=ensembl_ids)
        out[state_name] = anndata.AnnData(X=counts, obs=obs, var=var)

    return out
