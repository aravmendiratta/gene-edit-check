"""Efficacy-hazard coupling analysis (`benchmarks/hazard_coupling.py`).

When a verifier-guided search fails, there are two very different reasons it might have:

1. **The search was bad.** It looked in the wrong places, or wasted its budget.
2. **The target is constrained.** Every edit that works also trips a safety check, because the
   sequence features that control the gene *are* the features the check protects.

These demand opposite responses -- fix the optimiser, or stop looking -- so it is worth
measuring which one you are in. This module tests, across every variant the run scored, whether
efficacy and each hazard move together.

For PDCD1 the answer is unambiguous and mechanistically sensible: the harder an edit pushes PD-1
down, the more it perturbs CpG methylation. That is expected for this locus -- the PDCD1
regulatory regions (CR-B/CR-C) are CpG-dense, and T-cell exhaustion involves their
demethylation -- so the model is reproducing a real coupling rather than behaving arbitrarily.

**What this is and is not.** The correlation is between two *AlphaGenome 2 predictions*, so it is
evidence about what the model has learned, not a measurement of chromatin. Two things could
produce it: genuine biology (regulatory grammar at this locus is carried by CpGs), or a shared
internal representation that makes the model's expression and methylation heads move together
anywhere. The random-null subgroup distinguishes selection bias from coupling, but not those
two. Reporting the splice correlation alongside helps: splicing is scored by a different head,
so if methylation couples and splicing does not, a generic "everything correlates" explanation
gets harder to sustain.
"""

from __future__ import annotations

from typing import Literal

import numpy as np
from pydantic import BaseModel, Field
from scipy import stats

from gene_edit_check.models import OptimizationRunResult


class CouplingResult(BaseModel):
    """Correlation between edit efficacy and one hazard modality."""

    hazard: str
    n_total: int
    spearman_rho: float
    spearman_p: float
    # The null edits were chosen at random, *before* any candidate was verified, so a coupling
    # that survives inside this subgroup alone cannot be an artifact of selecting candidates for
    # efficacy. This is the control that makes the headline number mean something.
    n_null_only: int
    spearman_rho_null_only: float
    spearman_p_null_only: float
    # 95% percentile-bootstrap intervals (resampling edits with replacement). A p-value says whether
    # the correlation is distinguishable from zero; the interval says how large it plausibly is,
    # which is what a reader needs to judge a "moderate" association.
    ci_low: float = 0.0
    ci_high: float = 0.0
    ci_low_null_only: float = 0.0
    ci_high_null_only: float = 0.0
    interpretation: str = ""

    @property
    def is_significant(self) -> bool:
        return self.spearman_p < 0.05

    @property
    def survives_null_control(self) -> bool:
        return self.spearman_p_null_only < 0.05 and self.spearman_rho_null_only > 0


class HazardCouplingReport(BaseModel):
    gene_symbol: str
    couplings: list[CouplingResult] = Field(default_factory=list)
    binding_constraint: str | None = Field(
        default=None, description="Hazard that failed in the most near-miss candidates"
    )
    n_candidates_failing_binding_constraint: int = 0
    n_candidates: int = 0
    verdict: Literal["search_limited", "target_constrained", "inconclusive"] = "inconclusive"
    verdict_explanation: str = ""

    # Is the coupling specific to methylation, or does every hazard rise with efficacy?
    #
    # This is the discriminating measurement. If methylation couples and splicing does not, the
    # two come from different output heads, so a generic "the model's predictions all move
    # together" explanation struggles. If *both* couple, that explanation is live again and the
    # methylation finding should not be read as mechanistic. Reported as the gap between the two
    # correlations rather than left for a reader to eyeball.
    methylation_minus_splice_rho: float = 0.0
    # 95% bootstrap interval for that gap, resampling edits jointly so the two correlations stay
    # paired. If it contains zero, the data cannot tell methylation-specific coupling from every
    # hazard rising together; no cutoff is involved.
    gap_ci_low: float = 0.0
    gap_ci_high: float = 0.0
    coupling_class: Literal[
        "methylation_specific", "generic_hazard_coupling", "uncoupled", "unknown"
    ] = "unknown"
    coupling_class_explanation: str = ""


_HAZARD_FIELDS: dict[str, tuple[str, str]] = {
    # hazard name -> (per-candidate field on Tier2HazardResult, null field on EfficacyNullCalibration)
    "cpg_methylation": ("cpg_methylation_mean_abs_shift", "cpg_methylation_shift_values"),
    "splice_junction": ("max_splice_junction_delta", "splice_delta_values"),
}


def _spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Spearman correlation, guarding the degenerate cases that make scipy return NaN."""
    if len(x) < 3 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return 0.0, 1.0
    rho, p = stats.spearmanr(x, y)
    if not np.isfinite(rho) or not np.isfinite(p):
        return 0.0, 1.0
    return float(rho), float(p)


_N_BOOT = 4000


def _rank_corr(x: np.ndarray, y: np.ndarray) -> float:
    """Spearman rho via ranks; 0.0 if either side is constant. Fast enough to run thousands of times."""
    if len(x) < 3:
        return 0.0
    rx = stats.rankdata(x)
    ry = stats.rankdata(y)
    if np.ptp(rx) == 0 or np.ptp(ry) == 0:
        return 0.0
    return float(np.corrcoef(rx, ry)[0, 1])


def bootstrap_spearman_ci(
    x: np.ndarray, y: np.ndarray, n_boot: int = _N_BOOT, seed: int = 0, level: float = 0.95
) -> tuple[float, float]:
    """Percentile bootstrap interval for Spearman rho. Deterministic for a given seed."""
    n = len(x)
    if n < 5 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return 0.0, 0.0
    rng = np.random.default_rng(seed)
    draws = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        draws[i] = _rank_corr(x[idx], y[idx])
    a = (1 - level) / 2
    return float(np.quantile(draws, a)), float(np.quantile(draws, 1 - a))


def bootstrap_rho_gap_ci(
    lfc: np.ndarray, haz_a: np.ndarray, haz_b: np.ndarray, n_boot: int = _N_BOOT, seed: int = 0,
    level: float = 0.95,
) -> tuple[float, float]:
    """Interval for rho(lfc, haz_a) - rho(lfc, haz_b), resampling edits jointly (paired)."""
    n = len(lfc)
    if n < 5:
        return 0.0, 0.0
    rng = np.random.default_rng(seed)
    draws = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        draws[i] = _rank_corr(lfc[idx], haz_a[idx]) - _rank_corr(lfc[idx], haz_b[idx])
    a = (1 - level) / 2
    return float(np.quantile(draws, a)), float(np.quantile(draws, 1 - a))


def analyze_hazard_coupling(result: OptimizationRunResult) -> HazardCouplingReport:
    """Test whether efficacy and each hazard move together across every scored variant.

    Pools two sources that the run already produced, so this costs no extra API calls:
    the empirical null (random PAM-valid edits scored before any candidate was verified) and
    the verified candidates themselves.
    """
    report = HazardCouplingReport(gene_symbol=result.gene_symbol, n_candidates=len(result.iterations))
    null = result.efficacy_null

    cand_lfc = np.array(
        [abs(step.report.tier1_selectivity.target_delta_sketch_b) for step in result.iterations]
    )
    null_lfc = np.array(null.abs_lfc_values) if null else np.array([])

    pooled: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for hazard, (cand_field, null_field) in _HAZARD_FIELDS.items():
        cand_haz = np.array(
            [getattr(step.report.tier2_hazards, cand_field) for step in result.iterations]
        )
        null_haz = np.array(getattr(null, null_field)) if null else np.array([])

        all_lfc = np.concatenate([null_lfc, cand_lfc]) if len(null_lfc) else cand_lfc
        all_haz = np.concatenate([null_haz, cand_haz]) if len(null_haz) else cand_haz
        rho, p = _spearman(all_lfc, all_haz)
        rho_n, p_n = _spearman(null_lfc, null_haz)
        lo, hi = bootstrap_spearman_ci(all_lfc, all_haz)
        lo_n, hi_n = bootstrap_spearman_ci(null_lfc, null_haz) if len(null_lfc) else (0.0, 0.0)
        pooled[hazard] = (all_lfc, all_haz)

        cr = CouplingResult(
            hazard=hazard,
            n_total=len(all_lfc),
            spearman_rho=round(rho, 4),
            spearman_p=float(f"{p:.3g}"),
            n_null_only=len(null_lfc),
            spearman_rho_null_only=round(rho_n, 4),
            spearman_p_null_only=float(f"{p_n:.3g}"),
            ci_low=round(lo, 3),
            ci_high=round(hi, 3),
            ci_low_null_only=round(lo_n, 3),
            ci_high_null_only=round(hi_n, 3),
        )
        if cr.is_significant and rho > 0:
            cr.interpretation = (
                f"Efficacy and {hazard} rise together (rho={rho:+.2f}, 95% CI {lo:+.2f} to {hi:+.2f}, p={p:.2g}): the edits that "
                f"move {result.gene_symbol} are the ones that perturb this modality"
                + (
                    ". Holds among randomly chosen edits too, so it is not an artifact of "
                    "selecting candidates for efficacy."
                    if cr.survives_null_control
                    else ", but it does not reach significance among the random edits alone, so "
                    "selection bias is not excluded."
                )
            )
        else:
            cr.interpretation = (
                f"No significant coupling to {hazard} (rho={rho:+.2f}, p={p:.2g}): this hazard is "
                "independent of how hard the edit pushes the target."
            )
        report.couplings.append(cr)

    # Classify the coupling. Methylation rising with efficacy is only mechanistically suggestive
    # if splicing -- a different output head -- does not do the same thing.
    meth = next((c for c in report.couplings if c.hazard == "cpg_methylation"), None)
    splice = next((c for c in report.couplings if c.hazard == "splice_junction"), None)
    if meth is not None and splice is not None:
        gap = meth.spearman_rho - splice.spearman_rho
        report.methylation_minus_splice_rho = round(gap, 4)
        lfc_all, meth_all = pooled["cpg_methylation"]
        _, splice_all = pooled["splice_junction"]
        g_lo, g_hi = bootstrap_rho_gap_ci(lfc_all, meth_all, splice_all, seed=1)
        report.gap_ci_low, report.gap_ci_high = round(g_lo, 3), round(g_hi, 3)
        meth_couples = meth.ci_low > 0.0

        if not meth_couples:
            report.coupling_class = "uncoupled"
            report.coupling_class_explanation = (
                f"The interval for the efficacy-methylation correlation includes zero "
                f"(rho={meth.spearman_rho:+.2f}, 95% CI {meth.ci_low:+.2f} to {meth.ci_high:+.2f}), so "
                f"the data do not show that edits which silence {result.gene_symbol} disturb methylation."
            )
        elif g_lo > 0.0:
            report.coupling_class = "methylation_specific"
            report.coupling_class_explanation = (
                f"Methylation couples to efficacy (rho={meth.spearman_rho:+.2f}, CI "
                f"{meth.ci_low:+.2f} to {meth.ci_high:+.2f}) more strongly than splicing does "
                f"(rho={splice.spearman_rho:+.2f}); the paired gap is {gap:+.2f} with 95% CI "
                f"{g_lo:+.2f} to {g_hi:+.2f}, which excludes zero. The two come from different output "
                "heads, so this is evidence against every prediction moving together, though still only "
                "a statement about model outputs."
            )
        else:
            report.coupling_class = "generic_hazard_coupling"
            report.coupling_class_explanation = (
                f"Methylation (rho={meth.spearman_rho:+.2f}) couples to efficacy, but the paired gap to "
                f"splicing (rho={splice.spearman_rho:+.2f}) is {gap:+.2f} with 95% CI {g_lo:+.2f} to "
                f"{g_hi:+.2f}, which includes zero. The data cannot separate a methylation-specific "
                "effect from strong edits perturbing every output the model scores, so do not read the "
                "methylation correlation as mechanistic."
            )

    # Which check actually blocked the candidates?
    reason_prefixes = {
        "cpg_methylation": "TIER2_METHYLATION",
        "splice_junction": "TIER2_SPLICE",
        "healthy_selectivity": "TIER1_SELECTIVITY",
        "neighbor_gene": "TIER1_NEIGHBOR",
        "wrong_direction": "TIER1_DIRECTION",
        "bystanders": "TIER2_BYSTANDER",
        "efficacy": "TIER1_EFFICACY",
        "cross_sketch": "TIER1_GOODHART",
    }
    counts = {
        name: sum(
            1
            for step in result.iterations
            if any(r.startswith(prefix) for r in step.report.failure_reasons)
        )
        for name, prefix in reason_prefixes.items()
    }
    if result.iterations and max(counts.values()) > 0:
        binding = max(counts, key=lambda k: counts[k])
        report.binding_constraint = binding
        report.n_candidates_failing_binding_constraint = counts[binding]

        coupled = next(
            (c for c in report.couplings if c.hazard == binding and c.is_significant and c.spearman_rho > 0),
            None,
        )
        frac = counts[binding] / len(result.iterations)
        if coupled is not None and frac >= 0.5:
            report.verdict = "target_constrained"
            report.verdict_explanation = (
                f"{counts[binding]}/{len(result.iterations)} candidates were blocked by "
                f"{binding}, and that hazard scales with efficacy (rho="
                f"{coupled.spearman_rho:+.2f}). Edits that work tend to trip it, so more iterations may not "
                f"help. This is an association between two model outputs, not proof that the two "
                f"cannot be separated; a larger random null and a different regulatory element are "
                f"the next things to try."
            )
        elif frac >= 0.5:
            report.verdict = "search_limited"
            report.verdict_explanation = (
                f"{counts[binding]}/{len(result.iterations)} candidates were blocked by "
                f"{binding}, but it does not scale with efficacy -- so edits that work without "
                f"tripping it should exist. This is a search problem: widen the proposal window or "
                f"spend more iterations."
            )
        else:
            report.verdict = "inconclusive"
            report.verdict_explanation = (
                f"No single check dominates the failures (most common: {binding}, "
                f"{counts[binding]}/{len(result.iterations)}). Candidates fail for mixed reasons."
            )
    return report
