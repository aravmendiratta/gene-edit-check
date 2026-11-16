"""Do the model's input gradients point at the PDCD1 regulatory elements found in a wet-lab screen?

The reference is a CRISPRi tiling screen in primary human CD8 T cells (Wang et al., bioRxiv 2026,
doi:10.64898/2026.01.27.702006). It reports PD-1-lowering elements at +14.8 kb and +16.1 kb from the
PDCD1 TSS, plus accessible regions at -4.4, -4.65 and -38 kb. The paper gives positions only relative
to the TSS (hg19), so this script uses the same TSS-relative offsets on the model's interval, which is
centred on the TSS. No genome-build conversion is needed.

    python scripts/validate_enhancers.py            # live; uses the response cache
    python scripts/validate_enhancers.py --width 500 --permutations 20000
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import typer

from gene_edit_check.budget import RunBudgetTracker
from gene_edit_check.cli import _resolve_backend
from gene_edit_check.pipeline import build_all_sketches, load_cohort
from gene_edit_check.scenarios.registry import SCENARIOS

# TSS-relative centres (bp) of elements with experimental support, and what the screen showed.
ELEMENTS = {
    "+14.8 kb CRE": (14_800, "CRISPRi lowers PD-1 (15-25%); Cas9 indels in TIL CD8 lower PD-1"),
    "+16.1 kb CRE": (16_100, "CRISPRi lowers PD-1 (15-25%); reporter +2.6x in CD8 T cells"),
    "-4.4 kb ACR": (-4_400, "CRISPRi sgRNAs lower PD-1"),
    "-4.65 kb ACR": (-4_650, "CRISPRi sgRNAs lower PD-1"),
    "-38 kb ACR": (-38_000, "CRISPRi sgRNAs lower PD-1"),
}


def _window_scores(profile: np.ndarray, width: int) -> np.ndarray:
    """Mean |gradient| in every window of `width` bp (stride 1)."""
    c = np.concatenate([[0.0], np.cumsum(np.abs(profile))])
    return (c[width:] - c[:-width]) / width


def main(
    scenario: str = "crc-pd1-exhaustion",
    width: int = typer.Option(500, help="Window width (bp) around each element centre."),
    permutations: int = typer.Option(20_000, help="Random windows for the null."),
    exclude_promoter_bp: int = typer.Option(2_500, help="Null windows avoid +/- this around the TSS."),
    seed: int = 42,
    output: str = "outputs/enhancer_concordance.json",
) -> None:
    spec = SCENARIOS[scenario]
    budget = RunBudgetTracker()
    backend = _resolve_backend(use_simulator=False, budget=budget, seed=seed, require_live=True)
    cohort = load_cohort(spec, pool_size=600, seed=seed)
    target, healthy, _ = build_all_sketches(backend, spec, cohort, 128, seed)
    h_name = next(iter(healthy))

    gt = backend.compute_interval_gradients(spec.gene_id, spec.gene_symbol, target.sketch_a, target.cell_state_name)
    gh = backend.compute_interval_gradients(spec.gene_id, spec.gene_symbol, healthy[h_name].sketch_a, h_name)
    n = len(gt.signed_gradients_1bp)
    mid = n // 2
    tss_genomic = gt.interval.start + mid
    print(f"interval {gt.interval.chromosome}:{gt.interval.start}-{gt.interval.end} ({n} bp), TSS ~ {tss_genomic}")

    profiles = {
        "target saliency (exhausted TIL)": np.abs(gt.signed_gradients_1bp),
        "differential (target - healthy)": np.abs(gt.signed_gradients_1bp) - np.abs(gh.signed_gradients_1bp),
    }
    rng = np.random.default_rng(seed)
    out: dict = {
        "reference": "Wang et al., bioRxiv 2026, doi:10.64898/2026.01.27.702006 (CRISPRi tiling screen, primary human CD8 T cells)",
        "interval": f"{gt.interval.chromosome}:{gt.interval.start}-{gt.interval.end}",
        "tss_genomic_estimate": int(tss_genomic),
        "window_bp": width,
        "permutations": permutations,
        "profiles": {},
    }
    for pname, prof in profiles.items():
        scores = _window_scores(prof, width)
        centres = np.arange(len(scores)) + width // 2  # offset from interval start
        rel = centres - mid
        allowed = np.abs(rel) > exclude_promoter_bp
        pool = scores[allowed]
        rows = {}
        for name, (off, note) in ELEMENTS.items():
            i = mid + off - width // 2
            if i < 0 or i >= len(scores):
                rows[name] = {"in_interval": False, "note": note}
                continue
            s = float(scores[i])
            pct = float((pool < s).mean())
            p = float((np.sum(pool >= s) + 1) / (len(pool) + 1))
            rows[name] = {"offset_bp": off, "score": s, "percentile_vs_null": pct, "p_one_sided": p, "note": note}
        hits = [r for r in rows.values() if r.get("in_interval", True) and "score" in r]
        ranks = [r["percentile_vs_null"] for r in hits]
        # Joint test: mean percentile of the 5 elements vs 5 random windows, 20k times.
        obs = float(np.mean(ranks)) if ranks else float("nan")
        sp = np.sort(pool)
        draws = rng.choice(len(pool), size=(permutations, max(len(hits), 1)))
        null_mean = (np.searchsorted(sp, pool[draws], side="left") / len(sp)).mean(axis=1)
        p_joint = float((np.sum(null_mean >= obs) + 1) / (len(null_mean) + 1))
        out["profiles"][pname] = {
            "elements": rows,
            "mean_percentile": obs,
            "joint_p_one_sided": p_joint,
            "null_windows": int(len(pool)),
            "top_windows_outside_promoter": [
                {"offset_bp": int(rel[allowed][j]), "score": float(pool[j])}
                for j in np.argsort(pool)[::-1][:5]
            ],
        }
        print(f"\n== {pname} (window {width} bp, {len(pool)} null windows)")
        for k, r in rows.items():
            if "score" in r:
                print(f"  {k:14s} offset {r['offset_bp']:>7}  percentile {r['percentile_vs_null']:.3f}  p={r['p_one_sided']:.4f}")
            else:
                print(f"  {k:14s} outside the model interval")
        print(f"  mean percentile {obs:.3f}, joint p ~ {p_joint:.4f}")

    Path(output).write_text(json.dumps(out, indent=2), encoding="utf-8")
    print("\nwrote", output)


if __name__ == "__main__":
    typer.run(main)
