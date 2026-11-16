"""Evidence-tiered reference-condition selection (`conditions.py`).

AlphaGenome 2's `supported_conditions()` returns ~2,800 human reference conditions and, for
each of 15 output modalities, *how the model knows about it*:

| tier                  | meaning                                                      |
|-----------------------|--------------------------------------------------------------|
| `observed`            | real training tracks exist for this condition + modality     |
| `imputed_supported`   | no tracks, but imputation is validated for this modality     |
| `imputed_exploratory` | imputed with no validation -- weakest evidence               |

This matters a lot for Tier 2 hazard calls. A methylation or splicing verdict taken from an
`imputed_exploratory` condition is a guess about a guess. Only 185 of 2,799 conditions have
*observed* methylation. So instead of hard-coding CURIEs, we pick the reference condition that
has the strongest evidence for the modality we are about to trust, and we **record the tier we
actually got** in the report, so a reader can discount the number appropriately.

Consequence for this project's PD-1 scenario: `UBERON:0000178_CL:0000625` (blood CD8+ alpha-beta
T cell) is the *only* CD8 T-cell condition with `observed` methylation **and** `observed` splice
junctions **and** `observed` RNA-seq **and** `observed` TF ChIP. The tissue-matched alternative,
`UBERON:0001155_CL:0000625` (colon CD8+ alpha-beta T cell), is `imputed_exploratory` for both
hazard modalities. We therefore keep blood CD8 as the *primary* hazard reference and use colon
CD8 only as a tissue-matched **concordance** check -- a disagreement is reported, never averaged
away.
"""

from __future__ import annotations

import functools
from enum import IntEnum
from pathlib import Path

import pandas as pd

# Shipped with the repo so condition selection is reproducible and works offline.
CONDITIONS_TSV = Path(__file__).parent / "data" / "ag2_supported_conditions.tsv.gz"

HAZARD_MODALITIES: tuple[str, ...] = ("methylation", "splice_junctions", "rna_seq", "chip_tf")


class EvidenceTier(IntEnum):
    """Ordered so that `max()` picks the strongest evidence."""

    UNKNOWN = 0
    IMPUTED_EXPLORATORY = 1
    IMPUTED_SUPPORTED = 2
    OBSERVED = 3

    @classmethod
    def parse(cls, raw: object) -> "EvidenceTier":
        key = str(raw).strip().lower()
        return {
            "observed": cls.OBSERVED,
            "imputed_supported": cls.IMPUTED_SUPPORTED,
            "imputed_exploratory": cls.IMPUTED_EXPLORATORY,
        }.get(key, cls.UNKNOWN)

    @property
    def label(self) -> str:
        return self.name.lower()

    @property
    def is_trustworthy(self) -> bool:
        """True when a hazard verdict from this tier is worth acting on."""
        return self >= EvidenceTier.IMPUTED_SUPPORTED


@functools.lru_cache(maxsize=1)
def load_conditions_table() -> pd.DataFrame:
    """The cached `supported_conditions()` table shipped with the repo.

    Refresh it with `crispr-reflexion conditions --refresh` (needs an API key).
    """
    if not CONDITIONS_TSV.exists():
        raise FileNotFoundError(
            f"{CONDITIONS_TSV} is missing. Run `crispr-reflexion conditions --refresh` "
            "to download it from the AlphaGenome 2 API."
        )
    return pd.read_csv(CONDITIONS_TSV, sep="\t")


def refresh_conditions_table(api_key: str) -> Path:
    """Fetch `supported_conditions()` live and write the shipped TSV."""
    from alphagenome.models.v2 import dna_client

    df = dna_client.create(api_key).supported_conditions()
    CONDITIONS_TSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(CONDITIONS_TSV, sep="\t", index=False, compression="gzip")
    load_conditions_table.cache_clear()
    return CONDITIONS_TSV


def describe_condition(curie: str) -> str:
    """Human-readable biosample name for a condition CURIE, or the CURIE itself."""
    df = load_conditions_table()
    hit = df.loc[df["condition"] == curie, "biosample_name"]
    return str(hit.iloc[0]) if len(hit) else curie


def evidence_for(curie: str, modality: str) -> EvidenceTier:
    """Evidence tier of `modality` for one condition CURIE."""
    df = load_conditions_table()
    if modality not in df.columns:
        return EvidenceTier.UNKNOWN
    hit = df.loc[df["condition"] == curie, modality]
    return EvidenceTier.parse(hit.iloc[0]) if len(hit) else EvidenceTier.UNKNOWN


def evidence_profile(curie: str, modalities: tuple[str, ...] = HAZARD_MODALITIES) -> dict[str, str]:
    """`{modality: tier_label}` for one condition -- what gets stamped into a report."""
    return {m: evidence_for(curie, m).label for m in modalities}


def search_conditions(query: str, modality: str | None = None) -> pd.DataFrame:
    """Case-insensitive substring search over biosample names.

    With `modality`, results are sorted strongest-evidence-first for that modality.
    """
    df = load_conditions_table()
    mask = df["biosample_name"].str.contains(query, case=False, na=False, regex=True)
    out = df.loc[mask].copy()
    if modality and modality in out.columns:
        out["_tier"] = out[modality].map(lambda v: int(EvidenceTier.parse(v)))
        out = out.sort_values("_tier", ascending=False).drop(columns="_tier")
    return out


def rank_by_evidence(
    curies: list[str],
    modalities: tuple[str, ...] = HAZARD_MODALITIES,
) -> list[tuple[str, int, dict[str, str]]]:
    """Rank candidate conditions by total evidence across `modalities`, strongest first.

    Returns `(curie, total_score, per_modality_tiers)`. The score is the sum of tier values,
    so a condition with everything `observed` beats one that is `observed` for a single
    modality and `imputed_exploratory` elsewhere.
    """
    scored = [
        (c, sum(int(evidence_for(c, m)) for m in modalities), evidence_profile(c, modalities))
        for c in curies
    ]
    return sorted(scored, key=lambda t: t[1], reverse=True)


def best_evidence_condition(
    curies: list[str],
    modalities: tuple[str, ...] = HAZARD_MODALITIES,
) -> str:
    """The candidate condition with the strongest total evidence for `modalities`."""
    if not curies:
        raise ValueError("no candidate conditions given")
    return rank_by_evidence(curies, modalities)[0][0]


def weakest_hazard_tier(curie: str) -> EvidenceTier:
    """The limiting evidence tier for hazard calls made under this condition.

    A hazard panel is only as good as its weakest modality, so this is what the report
    headline should quote.
    """
    return min(
        (evidence_for(curie, m) for m in ("methylation", "splice_junctions")),
        default=EvidenceTier.UNKNOWN,
    )
