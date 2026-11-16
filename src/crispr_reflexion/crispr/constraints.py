"""Biophysical CRISPR Base-Editing Constraint Scanner (`crispr/constraints.py`).

Enforces physical genome-editing rules for:
- Adenine Base Editors (`ABE8e`): A->G on target strand (or T->C on opposite strand)
  within canonical activity window positions [4..8] (1-indexed from PAM-distal 5' end).
- Cytosine Base Editors (`CBE`): C->T on target strand (or G->A on opposite strand)
  within activity window positions [4..8].
- PAM grammar: canonical SpCas9 (`NGG`) or engineered Cas9-NG (`NG`).
- Bystander edit detection: flags any additional editable substrate base inside [4..8]
  of the same 20-nt protospacer.
"""

from __future__ import annotations

from crispr_reflexion.models import (
    BaseEditorType,
    CandidateGuideEdit,
    SingleNucleotideVariant,
)

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


def reverse_complement(seq: str) -> str:
    return "".join(COMPLEMENT.get(b, "N") for b in reversed(seq.upper()))


def _matches_pam(pam_seq: str, pam_rule: str = "NGG") -> bool:
    if len(pam_seq) < len(pam_rule):
        return False
    for b, r in zip(pam_seq.upper(), pam_rule.upper()):
        if r == "N":
            continue
        if b != r:
            return False
    return True


def enumerate_valid_base_edits(
    chromosome: str,
    window_start_genomic: int,
    window_sequence: str,
    gene_symbol: str,
    pam_rule: str = "NGG",
    activity_window: tuple[int, int] = (4, 8),
) -> list[CandidateGuideEdit]:
    """Scan a genomic sequence window for all PAM-valid ABE8e and CBE base edits.

    Args:
        chromosome: Chromosome identifier (e.g. 'chr2').
        window_start_genomic: 0-based genomic start coordinate of `window_sequence`.
        window_sequence: DNA sequence string (A/C/G/T).
        gene_symbol: Target gene symbol.
        pam_rule: PAM motif ('NGG' or 'NG').
        activity_window: Inclusive 1-based protospacer indices where deaminase is active.
    """
    seq = window_sequence.upper()
    pam_len = len(pam_rule)
    candidates: list[CandidateGuideEdit] = []
    win_lo, win_hi = activity_window

    # 1. Scan forward (+) strand: [20-nt protospacer][PAM]
    for i in range(0, len(seq) - (20 + pam_len) + 1):
        protospacer = seq[i : i + 20]
        pam = seq[i + 20 : i + 20 + pam_len]
        if not _matches_pam(pam, pam_rule):
            continue

        for editor, ref_base, alt_base in (
            (BaseEditorType.ABE8E, "A", "G"),
            (BaseEditorType.CBE, "C", "T"),
        ):
            editable_positions = [
                p_idx
                for p_idx in range(win_lo, win_hi + 1)
                if protospacer[p_idx - 1] == ref_base
            ]
            for target_p_idx in editable_positions:
                genomic_pos_1based = window_start_genomic + i + (target_p_idx - 1) + 1
                target_snv = SingleNucleotideVariant(
                    chromosome=chromosome,
                    position=genomic_pos_1based,
                    reference=ref_base,  # type: ignore[arg-type]
                    alternate=alt_base,  # type: ignore[arg-type]
                    protospacer_index=target_p_idx,
                )
                bystanders = [
                    SingleNucleotideVariant(
                        chromosome=chromosome,
                        position=window_start_genomic + i + (b_idx - 1) + 1,
                        reference=ref_base,  # type: ignore[arg-type]
                        alternate=alt_base,  # type: ignore[arg-type]
                        protospacer_index=b_idx,
                    )
                    for b_idx in editable_positions
                    if b_idx != target_p_idx
                ]
                cand_id = f"{gene_symbol}_{editor.value}_pos{genomic_pos_1based}_plus_p{target_p_idx}"
                candidates.append(
                    CandidateGuideEdit(
                        candidate_id=cand_id,
                        gene_symbol=gene_symbol,
                        editor=editor,
                        strand="+",
                        protospacer_20nt=protospacer,
                        pam_sequence=pam,
                        genomic_start=window_start_genomic + i,
                        genomic_end=window_start_genomic + i + 20 + pam_len,
                        target_edit=target_snv,
                        bystander_edits_in_window=bystanders,
                    )
                )

    # 2. Scan reverse (-) strand: [RC(PAM)][RC(20-nt protospacer)]
    for i in range(0, len(seq) - (20 + pam_len) + 1):
        sub = seq[i : i + 20 + pam_len]
        rc_sub = reverse_complement(sub)
        protospacer = rc_sub[:20]
        pam = rc_sub[20 : 20 + pam_len]
        if not _matches_pam(pam, pam_rule):
            continue

        for editor, minus_ref, minus_alt, plus_ref, plus_alt in (
            (BaseEditorType.ABE8E, "A", "G", "T", "C"),
            (BaseEditorType.CBE, "C", "T", "G", "A"),
        ):
            editable_positions = [
                p_idx
                for p_idx in range(win_lo, win_hi + 1)
                if protospacer[p_idx - 1] == minus_ref
            ]
            for target_p_idx in editable_positions:
                # Map minus-strand protospacer index back to forward genomic coordinate
                offset_from_sub_end = target_p_idx - 1
                plus_index_in_sub = (20 + pam_len - 1) - offset_from_sub_end
                genomic_pos_1based = window_start_genomic + i + plus_index_in_sub + 1

                target_snv = SingleNucleotideVariant(
                    chromosome=chromosome,
                    position=genomic_pos_1based,
                    reference=plus_ref,  # type: ignore[arg-type]
                    alternate=plus_alt,  # type: ignore[arg-type]
                    protospacer_index=target_p_idx,
                )
                bystanders = []
                for b_idx in editable_positions:
                    if b_idx == target_p_idx:
                        continue
                    b_plus_idx = (20 + pam_len - 1) - (b_idx - 1)
                    bystanders.append(
                        SingleNucleotideVariant(
                            chromosome=chromosome,
                            position=window_start_genomic + i + b_plus_idx + 1,
                            reference=plus_ref,  # type: ignore[arg-type]
                            alternate=plus_alt,  # type: ignore[arg-type]
                            protospacer_index=b_idx,
                        )
                    )
                cand_id = f"{gene_symbol}_{editor.value}_pos{genomic_pos_1based}_minus_p{target_p_idx}"
                candidates.append(
                    CandidateGuideEdit(
                        candidate_id=cand_id,
                        gene_symbol=gene_symbol,
                        editor=editor,
                        strand="-",
                        protospacer_20nt=protospacer,
                        pam_sequence=pam,
                        genomic_start=window_start_genomic + i,
                        genomic_end=window_start_genomic + i + 20 + pam_len,
                        target_edit=target_snv,
                        bystander_edits_in_window=bystanders,
                    )
                )

    return candidates
