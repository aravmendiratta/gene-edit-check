"""Tests for API & ISM window budget guardrails (`budget.py`)."""

from __future__ import annotations

import pytest
from crispr_reflexion.budget import BudgetExceededError, RunBudgetTracker


def test_oversized_ism_window_raises() -> None:
    budget = RunBudgetTracker(max_ism_window_bp=60, max_total_ism_variants=240)
    with pytest.raises(BudgetExceededError, match="ISM interval width"):
        budget.check_and_record_ism_window(window_bp=120, n_variants=360)


def test_max_conditioning_cells_cap_raises() -> None:
    budget = RunBudgetTracker()
    budget.check_conditioning_size(256)
    with pytest.raises(BudgetExceededError, match="MAX_CONDITIONING_CELLS is 256"):
        budget.check_conditioning_size(300)


def test_encode_cells_is_not_capped_at_256_but_counts_requests() -> None:
    # encode_cells chunks internally (1000 cells/request); only conditioning is capped at 256.
    budget = RunBudgetTracker()
    budget.check_and_record_encode_cells(1500)
    assert budget.encode_cells_calls_used == 2
    assert budget.encode_cells_total == 1500


def test_encode_cells_total_cell_budget_raises() -> None:
    budget = RunBudgetTracker(max_encode_cells_total=1000)
    budget.check_and_record_encode_cells(800)
    with pytest.raises(BudgetExceededError, match="total-cell budget"):
        budget.check_and_record_encode_cells(300)


def test_gradient_call_cap_raises() -> None:
    budget = RunBudgetTracker(max_gradient_calls=2)
    budget.check_and_record_gradient_call(1)
    budget.check_and_record_gradient_call(1)
    with pytest.raises(BudgetExceededError, match="Gradient call budget exceeded"):
        budget.check_and_record_gradient_call(1)
