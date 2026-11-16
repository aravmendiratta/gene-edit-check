"""Regulatory Motif & Sequence Feature Scanner (`crispr/motifs.py`).

Scans candidate protospacer windows for canonical transcription factor binding motifs
(JASPAR consensus patterns) to provide mechanistic attribution alongside gradient
and ISM scores.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class MotifMatch:
    motif_name: str
    family: str
    start_offset: int
    end_offset: int
    matched_sequence: str


# Curated immune, exhaustion, neural, and insulator regulatory consensus motifs
REGULATORY_MOTIFS: tuple[tuple[str, str, tuple[str, ...]], ...] = (
    ("TOX/HMG", "Exhaustion_Chromatin_Remodeler", ("GCAAATG", "CAAATGG", "ATTTGC", "CCATTTG")),
    ("NR4A1/NUR77", "Exhaustion_Nuclear_Receptor", ("AAAGGTCA", "TGACCTTT", "AAGGTCA", "TGACCTT")),
    ("NFAT", "TCR_Activation", ("GGAAA", "TTTCC", "TGGAAA", "TTTCCA")),
    ("AP-1 (FOS/JUN)", "Stress_Activation", ("TGACTCA", "TGAGTCA", "TGACT", "AGTCA")),
    ("RUNX1/3", "Core_T_Cell_Lineage", ("TGTGGTT", "AACCACA", "TGTGGT", "ACCACA")),
    ("CTCF", "Chromatin_Insulator", ("CCCTC", "GAGGG", "CCACTAGGGGGC")),
    ("STAT1/3", "Cytokine_Signaling", ("TTCCNGGAA", "TTCCAGGAA", "TTCCTGGAA")),
    ("PD-1_CR-C_Enhancer", "Exhaustion_Cis_Element", ("CAAAG", "CTTTG", "TGAC", "GTCA")),
)


def annotate_motifs_in_window(
    window_seq: str, edit_offset_in_window: int
) -> list[str]:
    """Return all regulatory motifs overlapping `edit_offset_in_window`."""
    seq = window_seq.upper()
    hits: list[str] = []
    for motif_name, family, patterns in REGULATORY_MOTIFS:
        for pat in patterns:
            idx = seq.find(pat)
            while idx != -1:
                if idx - 2 <= edit_offset_in_window <= idx + len(pat) + 1:
                    label = f"{motif_name} ({family})"
                    if label not in hits:
                        hits.append(label)
                idx = seq.find(pat, idx + 1)
    return hits
