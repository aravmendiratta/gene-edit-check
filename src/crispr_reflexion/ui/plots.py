"""Publication-Grade 4-Panel Figure Exporter (`ui/plots.py`).

Generates a dark-mode summary plot explicitly stamped with:
- Backend Provenance (`LIVE_ALPHAGENOME_V2` vs `SURROGATE_SIMULATOR`)
- Panel A: Donor-Stratified Coreset MMD Reduction in `encode_cells()` 256-D Latent Space
- Panel B: Tier 1 Single-Cell Selectivity (`Sketch A` vs Held-Out `Sketch B` vs Healthy Control)
- Panel C: Tier 2 Sequence-Intrinsic Hazards (`SpliceJunctionScorer`, `MethylationScorer(1201)`, Bystanders)
- Panel D: Efficacy-vs-hazard coupling (random null + candidates), which says whether the
  search or the locus is the binding constraint
"""

from __future__ import annotations

from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from crispr_reflexion.models import OptimizationRunResult


def export_optimization_figure(
    result: OptimizationRunResult,
    output_path: Path | str = "outputs/crispr_reflexion_report.png",
) -> Path:
    out_file = Path(output_path)
    out_file.parent.mkdir(parents=True, exist_ok=True)

    plt.style.use("dark_background")
    fig, axes = plt.subplots(2, 2, figsize=(14, 9), dpi=160)
    fig.patch.set_facecolor("#0b0f19")
    for ax in axes.flat:
        ax.set_facecolor("#111827")

    prov = result.provenance
    data_src = str(result.data_provenance.get("data_source", "unknown"))
    fig.suptitle(
        f"CRISPR-Reflexion: {result.gene_symbol} ({result.locus.chromosome}:{result.locus.start:,}-{result.locus.end:,})"
        f" | goal: {result.desired_direction}-regulate" + chr(10) +
        f"[{prov.backend_mode.value}] {prov.model_version} | seed {prov.random_seed} | " +
        (
            f"API calls {prov.total_api_calls} (cache hits {prov.cache_hits})"
            + (" - fully cached replay of an earlier live run" if prov.total_api_calls == 0 else "")
        ) + chr(10) +
        f"Data: {data_src[:95]}",
        fontsize=10.5,
        fontweight="bold",
        color="#38bdf8",
        y=0.985,
    )

    # Panel A: Coreset MMD^2 Ablation in AlphaGenome 2 Latent Space
    ax_a = axes[0, 0]
    states = [m.cell_state_name for m in result.sketch_metrics]
    mmd_strat = [m.ag2_latent_mmd_stratified for m in result.sketch_metrics]
    mmd_rand = [m.ag2_latent_mmd_unstratified_random for m in result.sketch_metrics]
    x = np.arange(len(states))
    w = 0.35
    ax_a.bar(x - w / 2, mmd_rand, width=w, label="Unstratified Random (AG2 Default)", color="#f59e0b", alpha=0.85)
    ax_a.bar(x + w / 2, mmd_strat, width=w, label="Donor-Stratified Facility Location", color="#10b981", alpha=0.9)
    ax_a.set_xticks(x)
    ax_a.set_xticklabels(states, fontsize=9)
    ax_a.set_ylabel("RBF MMD² in AG2 256-D Latent Space (Lower = Better)")
    ax_a.set_title("A. Coreset Sketching Ablation (`encode_cells()` Space)", fontsize=10, fontweight="bold")
    ax_a.legend(fontsize=8)
    ax_a.grid(axis="y", alpha=0.2)

    # Panel B: Tier 1 Single-Cell Selectivity across Iterations
    ax_b = axes[0, 1]
    iters = [t.iteration for t in result.iterations]
    lfc_a = [abs(t.report.tier1_selectivity.target_delta_sketch_a) for t in result.iterations]
    lfc_b = [abs(t.report.tier1_selectivity.target_delta_sketch_b) for t in result.iterations]
    lfc_h = [
        abs(next(iter(t.report.tier1_selectivity.healthy_control_deltas.values()), 0.0))
        for t in result.iterations
    ]
    ax_b.plot(iters, lfc_a, "o--", label="Target |LFC| (Proposal Sketch A)", color="#60a5fa", linewidth=1.8)
    ax_b.plot(iters, lfc_b, "s-", label="Target |LFC| (Held-Out Sketch B)", color="#34d399", linewidth=2.4)
    ax_b.plot(iters, lfc_h, "^-", label="Healthy Control |LFC| (Held-Out Sketch B)", color="#f87171", linewidth=2.0)
    ax_b.set_xticks(iters)
    ax_b.set_xlabel("Reflexion Iteration")
    ax_b.set_ylabel("|Δ Expression| (`SC_GENE_EXPRESSION_DENOISED`)")
    ax_b.set_title("B. Tier 1 Single-Cell Selectivity & Sketch A→B Transfer", fontsize=10, fontweight="bold")
    ax_b.legend(fontsize=8)
    ax_b.grid(alpha=0.2)

    # Panel C: Tier 2 Sequence Hazards & Bystander Count
    ax_c = axes[1, 0]
    splice_vals = [t.report.tier2_hazards.max_splice_junction_delta for t in result.iterations]
    cpg_vals = [t.report.tier2_hazards.cpg_methylation_mean_abs_shift for t in result.iterations]
    bys_counts = [len(t.report.candidate.bystander_edits_in_window) for t in result.iterations]
    x_it = np.arange(len(iters))
    ax_c.bar(x_it - 0.2, splice_vals, width=0.2, label="SpliceJunctionScorer Max Δ", color="#a78bfa")
    ax_c.bar(x_it, cpg_vals, width=0.2, label="MethylationScorer(1201) Mean |Δ|", color="#38bdf8")
    ax_c_twin = ax_c.twinx()
    ax_c_twin.plot(x_it + 0.2, bys_counts, "D-", color="#fb7185", label="Bystander Edits in [4..8]", linewidth=2.0)
    ax_c_twin.set_ylabel("Bystander Base Count", color="#fb7185")
    ax_c.set_xticks(x_it)
    ax_c.set_xticklabels([f"Iter {i}" for i in iters])
    ax_c.set_ylim(0, max(max(splice_vals, default=0.1), max(cpg_vals, default=0.1)) * 1.28)
    ax_c.set_ylabel("Hazard Scorer Magnitude")
    ax_c.set_title("C. Tier 2 Sequence-Intrinsic & Reference Ontology Hazards", fontsize=10, fontweight="bold")
    ax_c.legend(loc="upper right", fontsize=8)
    ax_c.grid(axis="y", alpha=0.2)

    # Panel D: Efficacy-vs-hazard coupling -- why the search ended where it did.
    # Random edits (scored before any candidate was verified) and candidates are drawn together,
    # so a reader can see the trend is not an artifact of selecting candidates for efficacy.
    ax_d = axes[1, 1]
    hc = result.hazard_coupling or {}
    null = result.efficacy_null

    cand_lfc = [abs(t.report.tier1_selectivity.target_delta_sketch_b) for t in result.iterations]
    cand_cpg = [t.report.tier2_hazards.cpg_methylation_mean_abs_shift for t in result.iterations]
    cand_spl = [t.report.tier2_hazards.max_splice_junction_delta for t in result.iterations]

    if null:
        ax_d.scatter(
            null.abs_lfc_values, null.cpg_methylation_shift_values,
            s=42, c="#64748b", marker="o", alpha=0.85,
            label=f"CpG, {null.n_random_edits} random edits", zorder=3,
        )
        ax_d.scatter(
            null.abs_lfc_values, null.splice_delta_values,
            s=26, c="#64748b", marker="x", alpha=0.55, label="splice, random edits", zorder=2,
        )
    ax_d.scatter(cand_lfc, cand_cpg, s=90, c="#f472b6", marker="o",
                 edgecolors="white", linewidths=0.6, label="CpG, verified candidates", zorder=5)
    ax_d.scatter(cand_lfc, cand_spl, s=52, c="#38bdf8", marker="x",
                 label="splice, verified candidates", zorder=4)

    ax_d.axhline(0.12, color="#ef4444", ls="--", lw=1.1, alpha=0.85)
    ax_d.text(0.985, 0.12, " CpG limit 0.12", color="#ef4444", fontsize=7,
              va="bottom", ha="right", transform=ax_d.get_yaxis_transform())

    meth = next((c for c in hc.get("couplings", []) if c["hazard"] == "cpg_methylation"), None)
    spl = next((c for c in hc.get("couplings", []) if c["hazard"] == "splice_junction"), None)
    if meth and spl:
        ax_d.set_title(
            "D. Does the hazard scale with efficacy?" + chr(10) +
            f"CpG rho={meth['spearman_rho']:+.2f} (p={meth['spearman_p']:.1g})   "
            f"splice rho={spl['spearman_rho']:+.2f} (p={spl['spearman_p']:.2g})",
            fontsize=9.5, fontweight="bold",
        )
        ax_d.text(
            0.03, 0.97,
            f"pooled n={meth['n_total']}\n"
            f"random-edit-only (n={meth['n_null_only']}): rho={meth['spearman_rho_null_only']:+.2f} "
            f"(p={meth['spearman_p_null_only']:.2g})",
            transform=ax_d.transAxes, fontsize=7.5, color="#fbbf24",
            va="top", ha="left", linespacing=1.4,
        )
    else:
        ax_d.set_title("D. Efficacy vs sequence hazard", fontsize=10, fontweight="bold")

    ax_d.set_xlabel(f"|{result.gene_symbol} LFC| on held-out Sketch B  (efficacy)")
    ax_d.set_ylabel("Hazard magnitude")
    ax_d.set_yscale("symlog", linthresh=0.02)
    ax_d.legend(loc="lower right", fontsize=6.5, framealpha=0.3)
    ax_d.grid(alpha=0.2)

    fig.text(
        0.5,
        0.01,
        prov.research_use_only_notice,
        ha="center",
        fontsize=8,
        color="#9ca3af",
        style="italic",
    )
    plt.tight_layout(rect=(0, 0.03, 1, 0.93))
    fig.savefig(out_file)
    plt.close(fig)
    return out_file
