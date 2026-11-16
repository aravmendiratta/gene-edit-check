"""Tests for PAM scanning, ABE8e/CBE activity windows [4..8], and bystander detection."""

from __future__ import annotations

from gene_edit_check.crispr.constraints import enumerate_valid_base_edits
from gene_edit_check.models import BaseEditorType


def test_abe_clean_vs_bystander_detection() -> None:
    # 20-nt protospacer + NGG PAM:
    # Pos 1..3: GCT, Pos 4..8: ATGTG (only one 'A' at pos 4, one 'T' at pos 5 & 7)
    seq_clean_abe = "GCTATGTGGTTAGGCCAGGC" + "AGG"
    cands = enumerate_valid_base_edits(
        chromosome="chr2",
        window_start_genomic=1000,
        window_sequence=seq_clean_abe,
        gene_symbol="PDCD1",
        pam_rule="NGG",
    )
    abe_plus = [c for c in cands if c.editor == BaseEditorType.ABE8E and c.strand == "+"]
    assert len(abe_plus) == 1
    assert abe_plus[0].target_edit.protospacer_index == 4
    assert abe_plus[0].is_bystander_clean

    # Now place two 'A's inside [4..8] (pos 4 and pos 6):
    seq_bystander_abe = "GCTATATGGTTAGGCCAGGC" + "AGG"
    cands_bys = enumerate_valid_base_edits(
        chromosome="chr2",
        window_start_genomic=1000,
        window_sequence=seq_bystander_abe,
        gene_symbol="PDCD1",
        pam_rule="NGG",
    )
    abe_bys = [c for c in cands_bys if c.editor == BaseEditorType.ABE8E and c.strand == "+"]
    assert len(abe_bys) == 2
    assert not abe_bys[0].is_bystander_clean
    assert len(abe_bys[0].bystander_edits_in_window) == 1
