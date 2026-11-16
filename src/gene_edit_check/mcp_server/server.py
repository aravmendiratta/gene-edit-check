"""MCP server for coding agents (`mcp_server/server.py`).

Run with `gene-edit-check-mcp` (needs the optional extra: `pip install "gene-edit-check[mcp]"`).
The tool functions below are plain Python so they can be tested and reused without MCP installed.
"""

from __future__ import annotations

import os

from dotenv import find_dotenv, load_dotenv

from gene_edit_check.budget import RunBudgetTracker
from gene_edit_check.cli import DEFAULT_SCENARIO, _resolve_backend
from gene_edit_check.conditions import search_conditions
from gene_edit_check.pipeline import run_scenario
from gene_edit_check.scenarios.registry import SCENARIOS

load_dotenv(find_dotenv(usecwd=True))


def run_optimization_tool(
    scenario_id: str = DEFAULT_SCENARIO,
    sketch_size: int = 128,
    pool_size: int = 600,
    max_iterations: int = 6,
    seed: int = 42,
    null_n: int = 12,
    use_simulator: bool = False,
) -> str:
    """Run the two-tier AlphaGenome 2 pipeline on a scenario and return the JSON report.

    Needs ALPHAGENOME_API_KEY for a real run. `use_simulator=True` is a scripted stand-in for
    plumbing tests only; its output is not a prediction.
    """
    if scenario_id not in SCENARIOS:
        raise ValueError(f"Unknown scenario '{scenario_id}'. Available: {sorted(SCENARIOS)}")
    spec = SCENARIOS[scenario_id]
    if not use_simulator and spec.data_source != "synthetic" and not os.environ.get("ALPHAGENOME_API_KEY", "").strip():
        raise RuntimeError(
            "ALPHAGENOME_API_KEY is not set. Real-data runs need the live model; get a free key at "
            "https://deepmind.google.com/science/alphagenome. use_simulator=True gives scripted output for plumbing tests only."
        )
    budget = RunBudgetTracker()
    backend = _resolve_backend(use_simulator=use_simulator, budget=budget, seed=seed)
    result = run_scenario(
        spec, backend, budget, sketch_size=sketch_size, max_iterations=max_iterations,
        seed=seed, pool_size=pool_size, null_n=null_n,
    )
    return result.model_dump_json(indent=2)


def list_scenarios_tool() -> list[dict[str, str]]:
    """List the registered scenarios (gene, data source, target and healthy cell states)."""
    return [
        {
            "scenario_id": s.scenario_id,
            "gene": s.gene_symbol,
            "data_source": s.data_source,
            "target_state": s.target_state_name,
            "healthy_states": ", ".join(s.healthy_state_names),
        }
        for s in SCENARIOS.values()
    ]


def search_conditions_tool(query: str, modality: str = "methylation", limit: int = 15) -> list[dict[str, str]]:
    """Search AlphaGenome 2 reference conditions, strongest evidence (observed > imputed) first."""
    hits = search_conditions(query, modality=modality).head(limit)
    cols = ["condition", "biosample_name", "methylation", "splice_junctions", "rna_seq"]
    return hits[cols].astype(str).to_dict(orient="records")


def main() -> None:
    try:
        from mcp.server.fastmcp import FastMCP
    except ImportError as exc:  # pragma: no cover - depends on the optional extra
        raise SystemExit('The MCP server needs the optional extra: pip install "gene-edit-check[mcp]"') from exc

    server = FastMCP("gene-edit-check")
    server.tool()(run_optimization_tool)
    server.tool()(list_scenarios_tool)
    server.tool()(search_conditions_tool)
    server.run()


if __name__ == "__main__":
    main()
