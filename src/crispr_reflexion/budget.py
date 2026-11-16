"""Hard API & compute budget guardrails (`budget.py`).

Prevents accidental 1-Mb full-window ISM scans (~3M variants) and enforces strict
per-run limits on gradient calls, scoped ISM window width, and total variant evaluations.
"""

from __future__ import annotations

from dataclasses import dataclass, field


class BudgetExceededError(RuntimeError):
    """Raised when an operation would exceed configured API or variant budgets."""


@dataclass
class RunBudgetTracker:
    max_gradient_calls: int = 8
    max_ism_window_bp: int = 60
    max_total_ism_variants: int = 420
    max_variant_score_calls: int = 90
    max_encode_cells_calls: int = 12
    max_encode_cells_total: int = 4000
    encode_cells_per_request: int = 1000  # alphagenome MAX_ENCODE_CELLS_PER_REQUEST

    gradient_calls_used: int = field(default=0, init=False)
    ism_variants_scored: int = field(default=0, init=False)
    variant_score_calls_used: int = field(default=0, init=False)
    encode_cells_calls_used: int = field(default=0, init=False)
    encode_cells_total: int = field(default=0, init=False)
    cache_hits: int = field(default=0, init=False)

    def check_and_record_gradient_call(self, count: int = 1) -> None:
        if self.gradient_calls_used + count > self.max_gradient_calls:
            raise BudgetExceededError(
                f"Gradient call budget exceeded: requested +{count}, "
                f"used {self.gradient_calls_used}/{self.max_gradient_calls}."
            )
        self.gradient_calls_used += count

    def check_and_record_ism_window(self, window_bp: int, n_variants: int) -> None:
        if window_bp > self.max_ism_window_bp:
            raise BudgetExceededError(
                f"ISM interval width ({window_bp} bp) exceeds hard cap "
                f"({self.max_ism_window_bp} bp). Scope ISM strictly to PAM-valid windows."
            )
        if self.ism_variants_scored + n_variants > self.max_total_ism_variants:
            raise BudgetExceededError(
                f"Total ISM variant budget exceeded: requested +{n_variants}, "
                f"used {self.ism_variants_scored}/{self.max_total_ism_variants}."
            )
        self.ism_variants_scored += n_variants

    def check_and_record_variant_score(self, count: int = 1) -> None:
        if self.variant_score_calls_used + count > self.max_variant_score_calls:
            raise BudgetExceededError(
                f"Variant scoring call budget exceeded: requested +{count}, "
                f"used {self.variant_score_calls_used}/{self.max_variant_score_calls}."
            )
        self.variant_score_calls_used += count

    def check_and_record_encode_cells(self, n_cells: int) -> None:
        """`encode_cells()` chunks internally (<= MAX_ENCODE_CELLS_PER_REQUEST per request), so the
        256-cell cap applies to `cell_conditioning`, not here. We budget API *requests* and total cells."""
        n_requests = max(1, -(-n_cells // self.encode_cells_per_request))
        if self.encode_cells_total + n_cells > self.max_encode_cells_total:
            raise BudgetExceededError(
                f"encode_cells() total-cell budget exceeded: +{n_cells}, "
                f"used {self.encode_cells_total}/{self.max_encode_cells_total}."
            )
        if self.encode_cells_calls_used + n_requests > self.max_encode_cells_calls:
            raise BudgetExceededError(
                f"encode_cells() request budget exceeded ({self.max_encode_cells_calls})."
            )
        self.encode_cells_total += n_cells
        self.encode_cells_calls_used += n_requests

    def check_conditioning_size(self, n_cells: int) -> None:
        if n_cells > 256:
            raise BudgetExceededError(
                f"cell_conditioning received {n_cells} cells; AlphaGenome 2 MAX_CONDITIONING_CELLS is 256."
            )

    def record_cache_hit(self) -> None:
        self.cache_hits += 1

    @property
    def total_api_calls(self) -> int:
        return (
            self.gradient_calls_used
            + self.variant_score_calls_used
            + self.encode_cells_calls_used
        )
