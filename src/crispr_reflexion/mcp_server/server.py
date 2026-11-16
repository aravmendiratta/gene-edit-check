"""MCP Server for Coding Agents (`mcp_server/server.py`)."""

from __future__ import annotations

from dotenv import load_dotenv

from crispr_reflexion.budget import RunBudgetTracker
from crispr_reflexion.cli import DEFAULT_SCENARIO, _resolve_backend
from crispr_reflexion.pipeline import run_scenario
from crispr_reflexion.scenarios.registry import SCENARIOS

load_dotenv()


def run_optimization_tool(
    scenario_id: str = DEFAULT_SCENARIO,
    sketch_size: int = 128,
    pool_size: int = 600,
    max_iterations: int = 6,
    seed: int = 42,
    use_simulator: bool = False,
) -> str:
    """Run the two-tier AlphaGenome 2 CRISPR-Reflexion pipeline and return the JSON report."""
    spec = SCENARIOS.get(scenario_id, SCENARIOS[DEFAULT_SCENARIO])
    budget = RunBudgetTracker()
    backend = _resolve_backend(use_simulator=use_simulator, budget=budget, seed=seed)
    result = run_scenario(
        spec, backend, budget, sketch_size=sketch_size, max_iterations=max_iterations,
        seed=seed, pool_size=pool_size,
    )
    return result.model_dump_json(indent=2)
