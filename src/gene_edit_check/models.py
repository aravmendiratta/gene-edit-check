"""Core Pydantic v2 data models for gene-edit-check.

Encodes the architectural separation between:
1. Tier 1 Selectivity (`cell_conditioning` on 256-cell sketches -> `SC_GENE_EXPRESSION_DENOISED`).
2. Tier 2 Sequence-Intrinsic Hazards (`ModelCondition.from_curies` -> `SpliceJunctionScorer`,
   `MethylationScorer(width=1201)`, and chromatin accessibility).
"""

from __future__ import annotations

from enum import Enum
from typing import Any, Literal
from pydantic import BaseModel, Field, computed_field


RESEARCH_USE_ONLY_NOTICE = (
    "RESEARCH USE ONLY: Outputs are computational hypotheses generated via "
    "AlphaGenome 2 and are not validated for clinical, diagnostic, or therapeutic use."
)


class BackendMode(str, Enum):
    LIVE_ALPHAGENOME_V2 = "LIVE_ALPHAGENOME_V2"
    SURROGATE_SIMULATOR = "SURROGATE_SIMULATOR"


class BaseEditorType(str, Enum):
    ABE8E = "ABE8e"  # A->G (sense) or T->C (antisense), protospacer window [4..8]
    CBE = "CBE"      # C->T (sense) or G->A (antisense), protospacer window [4..8]


class GenomicInterval(BaseModel):
    chromosome: str = Field(..., description="Chromosome name, e.g. 'chr2'")
    start: int = Field(..., ge=0, description="0-based start coordinate")
    end: int = Field(..., gt=0, description="0-based end coordinate")

    @property
    def width(self) -> int:
        return self.end - self.start


class SingleNucleotideVariant(BaseModel):
    chromosome: str
    position: int = Field(..., ge=0, description="1-based genomic position")
    reference: Literal["A", "C", "G", "T"]
    alternate: Literal["A", "C", "G", "T"]
    protospacer_index: int | None = Field(
        default=None,
        description="1-based position within the 20-nt protospacer (valid window: 4..8)",
    )

    @computed_field
    @property
    def atlas_url(self) -> str:
        """AlphaGenome Atlas deep-link following `google-deepmind/science-skills` specification."""
        win_start = max(1, self.position - 190)
        win_end = self.position + 190
        return (
            f"https://deepmind.google.com/science/alphagenome/atlas"
            f"?q={self.chromosome}:{self.position}:{self.reference}%3E{self.alternate}"
            f"&m=variant&i={self.chromosome}:{win_start}-{win_end}"
        )


class CandidateGuideEdit(BaseModel):
    candidate_id: str
    gene_symbol: str
    editor: BaseEditorType
    strand: Literal["+", "-"]
    protospacer_20nt: str = Field(..., min_length=20, max_length=20)
    pam_sequence: str = Field(..., min_length=2, max_length=4)
    genomic_start: int
    genomic_end: int
    target_edit: SingleNucleotideVariant
    bystander_edits_in_window: list[SingleNucleotideVariant] = Field(default_factory=list)
    disrupted_motif: str | None = None
    gradient_saliency_score: float = 0.0
    ism_verification_score: float = 0.0

    @property
    def is_bystander_clean(self) -> bool:
        return len(self.bystander_edits_in_window) == 0


class Tier1SelectivityResult(BaseModel):
    """Tier 1: Single-Cell Conditioned Selectivity (`OutputType.SC_GENE_EXPRESSION_DENOISED`).

    Explicitly separates Proposal Sketch A from Held-Out Verification Sketch B
    to prevent Goodhart self-confirmation.
    """

    resolution_tier: Literal["TIER1_SINGLE_CELL_CONDITIONED"] = "TIER1_SINGLE_CELL_CONDITIONED"
    output_modality: Literal["SC_GENE_EXPRESSION_DENOISED"] = "SC_GENE_EXPRESSION_DENOISED"
    target_condition_name: str
    target_delta_sketch_a: float = Field(
        ..., description="Fractional log2FC / expression shift on Proposal Sketch A (256 cells)"
    )
    target_delta_sketch_b: float = Field(
        ..., description="Fractional log2FC / expression shift on Held-Out Verification Sketch B (256 cells)"
    )
    target_sketch_bootstrap_se: float = Field(
        ..., ge=0.0, description="Bootstrap standard error across disjoint cell sketches"
    )
    healthy_control_deltas: dict[str, float] = Field(
        default_factory=dict,
        description="Map of healthy cell state name -> fractional expression shift on held-out sketch",
    )
    target_efficacy_passed: bool
    healthy_selectivity_passed: bool
    cross_sketch_generalization_passed: bool

    # Neighbour-gene specificity, read off the same expression response (no extra API call).
    neighbor_gene_lfcs: dict[str, float] = Field(
        default_factory=dict,
        description="Non-target genes in the scored interval -> LFC on held-out target Sketch B",
    )
    worst_neighbor_gene: str | None = None
    worst_neighbor_lfc: float = 0.0
    neighbor_specificity_passed: bool = True


class Tier2HazardResult(BaseModel):
    """Tier 2: Sequence-Intrinsic Hazards (`ModelCondition.from_curies`).

    Uses ontology-matched reference conditions (`UBERON`/`CL` CURIEs) with
    `SpliceJunctionScorer()` and `MethylationScorer(width=1201)`.
    Never claims single-cell sketch resolution.
    """

    resolution_tier: Literal["TIER2_ONTOLOGY_REFERENCE_CONDITION"] = (
        "TIER2_ONTOLOGY_REFERENCE_CONDITION"
    )
    reference_curies: list[str]
    max_splice_junction_delta: float = Field(
        ..., ge=0.0, description="Max delta from SpliceJunctionScorer()"
    )
    cpg_methylation_mean_abs_shift: float = Field(
        ..., ge=0.0, description="Mean abs shift across CpGs from MethylationScorer(width=1201)"
    )
    splicing_hazard_passed: bool
    methylation_hazard_passed: bool
    bystander_clean_passed: bool

    # How well AlphaGenome 2 actually knows this condition for the modality being trusted,
    # from `supported_conditions()`: observed > imputed_supported > imputed_exploratory.
    # Only 185 of 2,799 human conditions have *observed* methylation, so a verdict taken from an
    # exploratory-imputed condition deserves far less weight -- recorded, never silently dropped.
    reference_biosample_name: str = ""
    evidence_tiers: dict[str, str] = Field(
        default_factory=dict, description="modality -> evidence tier of the primary reference"
    )
    limiting_evidence_tier: str = "unknown"

    # Tissue-matched cross-check under a second, anatomically correct but weaker-evidence
    # condition (e.g. colon CD8 T cell). Disagreement is surfaced, not averaged away.
    tissue_matched_curies: list[str] = Field(default_factory=list)
    tissue_matched_biosample_name: str = ""
    tissue_matched_splice_delta: float | None = None
    tissue_matched_cpg_shift: float | None = None
    tissue_matched_agrees: bool | None = Field(
        default=None,
        description="True when the tissue-matched reference reaches the same pass/fail verdict",
    )


class OffTissueSafetyPanel(BaseModel):
    """Does this edit do anything in tissues it was never meant to touch?

    A base edit is delivered to cells, not to a cell *state*. An edit that silences PDCD1 in
    exhausted T cells but also perturbs the locus in hepatocytes or neurons is a different
    risk profile. Each row is one non-target reference condition scored with the same variant.
    """

    tissues: dict[str, float] = Field(
        default_factory=dict, description="biosample name -> |target-gene LFC| under that condition"
    )
    evidence_tiers: dict[str, str] = Field(default_factory=dict)
    max_abs_lfc: float = 0.0
    worst_tissue: str = ""
    threshold: float = 0.0
    passed: bool = True


class TwoTierVerificationReport(BaseModel):
    candidate: CandidateGuideEdit
    tier1_selectivity: Tier1SelectivityResult
    tier2_hazards: Tier2HazardResult
    composite_objective_score: float
    all_checks_passed: bool
    failure_reasons: list[str] = Field(default_factory=list)
    off_tissue_panel: OffTissueSafetyPanel | None = Field(
        default=None, description="Run only for finalists -- it costs one API call per tissue"
    )

    @property
    def n_checks_failed(self) -> int:
        return len(self.failure_reasons)


class ReflexionStepTrace(BaseModel):
    iteration: int
    report: TwoTierVerificationReport
    critic_diagnosis: str
    episodic_memory_rule_added: str


class ProvenanceRecord(BaseModel):
    backend_mode: BackendMode
    model_version: str
    random_seed: int
    total_api_calls: int
    cache_hits: int
    total_ism_variants_scored: int
    gradient_calls_used: int
    research_use_only_notice: str = RESEARCH_USE_ONLY_NOTICE


class SketchQualityMetrics(BaseModel):
    cell_state_name: str
    total_source_cells: int
    sketch_size: int = 256
    donor_stratified: bool
    split_by: Literal["donor", "cell"] = "cell"
    n_donors_sketch_a: int = 0
    n_donors_sketch_b: int = 0
    donors_overlap_a_b: int = 0
    ag2_latent_mmd_stratified: float
    ag2_latent_mmd_unstratified_random: float
    mmd_reduction_pct: float
    sketch_a_B_cosine_consistency: float


class EfficacyNullCalibration(BaseModel):
    """Empirical null for Tier 1 efficacy: |LFC| of random PAM-valid base edits (same locus, same
    held-out target sketch). A candidate must beat every random edit, so a pass means
    empirical one-sided p <= 1/(n+1) against 'any editable base near the TSS'."""

    n_random_edits: int
    abs_lfc_values: list[float]
    abs_lfc_max: float
    abs_lfc_mean: float
    cutoff_used: float
    empirical_p_floor: float
    # Context only (not used as thresholds): what random edits do on the Tier 2 hazard scorers.
    cpg_methylation_shift_values: list[float] = Field(default_factory=list)
    splice_delta_values: list[float] = Field(default_factory=list)


class IsmPrescreenSummary(BaseModel):
    """What the ISM pre-screen bought us.

    Gradients are a first-order signal: the sign tells you which way expression moves if you
    change the reference base, but not reliably, and a wrong guess costs a full verification
    round (3+ API calls). In-silico mutagenesis over a handful of compact windows measures the
    real effect of every substitution in those windows for one call each, so the verifier only
    ever sees candidates whose direction is already known to be right.
    """

    n_windows: int
    window_width_bp: int
    n_variants_measured: int
    n_candidates_with_measurement: int
    n_correct_direction: int
    n_wrong_direction_filtered: int
    # Did the gradient sign actually predict the measured ISM sign? This is the honest audit of
    # the heuristic the previous version of the search relied on outright.
    gradient_sign_agreement: float = Field(
        default=0.0, description="Fraction of measured variants whose gradient sign matched ISM"
    )
    best_measured_lfc: float = 0.0


class OptimizationRunResult(BaseModel):
    scenario_id: str
    gene_symbol: str
    locus: GenomicInterval
    provenance: ProvenanceRecord
    sketch_metrics: list[SketchQualityMetrics]
    iterations: list[ReflexionStepTrace]
    selected_candidate: TwoTierVerificationReport | None = None
    converged: bool
    desired_direction: Literal["down", "up"] = "down"
    efficacy_null: EfficacyNullCalibration | None = None
    ism_prescreen: IsmPrescreenSummary | None = None
    # Populated after the loop by `benchmarks.hazard_coupling`; kept as a plain dict so the
    # models module stays free of a dependency on the analysis layer.
    hazard_coupling: dict[str, Any] = Field(default_factory=dict)
    data_provenance: dict[str, Any] = Field(default_factory=dict)
