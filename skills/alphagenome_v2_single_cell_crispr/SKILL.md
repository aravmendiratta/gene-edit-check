---
name: alphagenome-v2-single-cell-crispr
description: >-
  Performs AlphaGenome 2 (`alphagenome.models.v2`) single-cell conditioned
  regulatory variant analysis, donor-stratified 256-cell coreset sketching in
  `encode_cells()` latent space, `predict_interval_gradients` saliency proposal,
  and two-tier CRISPR base-editing (`ABE8e` / `CBE`) verification across
  single-cell expression (`SC_GENE_EXPRESSION_DENOISED`), splicing
  (`SpliceJunctionScorer`), and CpG methylation (`MethylationScorer`). Use when
  designing or auditing state-selective non-coding CRISPR edits, sub-selecting
  single-cell `.h5ad` cohorts for `cell_conditioning` (`MAX_CONDITIONING_CELLS =
  256`), or evaluating regulatory variants across distinct cell states.
metadata:
  author: Arav Mendiratta
  version: "0.2.0"
---

# AlphaGenome 2 Single-Cell Conditioned Variant & CRISPR Base-Editing Verifier (`crispr-reflexion`)

## Prerequisites

1. **`uv` Installation**: Ensure `uv` is installed (`curl -LsSf https://astral.sh/uv/install.sh | sh`).
2. **AlphaGenome Non-Commercial Terms**: AlphaGenome 2 weights and API predictions are subject to the [AlphaGenome Non-Commercial Use Terms](https://deepmind.google.com/science/alphagenome/terms). All outputs carry an explicit `RESEARCH USE ONLY` notice and are not validated for clinical, diagnostic, or therapeutic use.
3. **API Key Setup**:
   Configure `ALPHAGENOME_API_KEY` in `.env` (or `~/.env`). Obtain a free non-commercial key at [https://deepmind.google.com/science/alphagenome](https://deepmind.google.com/science/alphagenome). If no key is present, `crispr-reflexion` can run in offline `--simulator` mode with an explicit `[SURROGATE_SIMULATOR]` provenance stamp.

---

## Core Rules

1. **Always Separate Tier 1 (`cell_conditioning`) and Tier 2 (`ModelCondition.from_curies`)**:
   - `alphagenome.models.v2` raises `ValueError('Cannot specify both condition and cell_conditioning.')` if both are passed simultaneously.
   - **Tier 1 (Single-Cell Selectivity):** Use `cell_conditioning` (`AnnData` sketch $\le 256$ cells) with `OutputType.SC_GENE_EXPRESSION_DENOISED`. Always verify candidates proposed on **Proposal Sketch A** against a **disjoint Held-Out Verification Sketch B** to prevent Goodhart self-confirmation.
   - **Tier 2 (Sequence-Intrinsic Hazards):** Evaluate `SpliceJunctionScorer()` and `MethylationScorer(width=1201)` against matched reference ontology conditions (`model_condition.from_curies`), and encode the resolution tier (`TIER1_SINGLE_CELL_CONDITIONED` vs. `TIER2_ONTOLOGY_REFERENCE_CONDITION`) explicitly in all JSON/Markdown reports.
2. **Donor-Stratified `encode_cells()` Sketching Over Uniform Random Subsampling**:
   - Never use unstratified random sampling (`gene_counts_condition.subsample`) on multi-donor `.h5ad` atlases, as `encode_cells()` does not remove batch effects. Always stratify by donor/batch and select coresets in the 256-D `encode_cells()` latent space (`sc_sketch/coreset.py`).
3. **Gradients Propose, Scoped ISM Verifies**:
   - Never run brute-force ISM across the full 131 kb context window ($\sim 393\text{k}$ variants). Use `predict_interval_gradients` (`1` call) to rank PAM-valid base-editor windows (`ABE8e` / `CBE`, deaminase window `[4..8]`), then run scoped ISM ($\le 60\text{ bp}$) only on top candidate windows, enforced by `CostBudgetTracker` (`budget.py`) and content-addressed SQLite caching (`cache.py`).
4. **Mandatory AlphaGenome Atlas Deep-Linking**:
   - Every reported variant (`chr:pos:ref>alt`) must include its clickable [AlphaGenome Atlas](https://deepmind.google.com/science/alphagenome/atlas) deep-link (`SingleNucleotideVariant.atlas_url`), formatted as `https://deepmind.google.com/science/alphagenome/atlas?q={chr}:{pos}:{ref}%3E{alt}&m=variant&i={chr}:{pos-190}-{pos+190}`.
5. **Companion Skills Integration (`google-deepmind/science-skills`)**:
   - Chain with `alphagenome-atlas-website-links` for multi-track Ref vs. Alt Sashimi/RNA-seq views (`/atlas/track-predictions`), `jaspar_database` for TF PWM verification, and `ensembl_database` for gene/transcript coordinate lookup.

---

## Standard Workflows & CLI Commands

### 1. Run Live Two-Tier CRISPR Reflexion Optimization (`optimize`)
```bash
uv run crispr-reflexion optimize --scenario pd1-exhaustion --sketch-size 16 --max-iterations 4
```
- **Outputs:**
  - Structured JSON report: `outputs/report_pd1-exhaustion.json` (with `atlas_url` and `ProvenanceRecord`)
  - 4-panel dark-mode publication figure: `outputs/figure_pd1-exhaustion.png`

### 2. Run Donor-Stratified `encode_cells()` Coreset Ablation (`ablate-sketch`)
```bash
uv run crispr-reflexion ablate-sketch --sketch-size 16
```
- Quantifies empirical $\text{MMD}^2$ reduction (`-77.9%` in `Exhausted_CD8_TIL`, `-69.1%` in `Healthy_Naive_CD8`) of donor-stratified `encode_cells()` coreset selection vs. naive uniform random subsampling.

### 3. Run Ground-Truth Regulatory Benchmark Validation (`validate-gt`)
```bash
uv run crispr-reflexion validate-gt
```
- Evaluates zero-shot AlphaGenome 2 alignment against published experimental regulatory effect sizes across the *PDCD1* locus (Spearman $\rho = +1.000, p < 0.0001$).

---

## Reference Examples (`docs/examples/`)

- **Positive Case Study (Converged State-Selective Edit):** [`../../docs/examples/example_1_pdcd1_abe8e_converged/report.md`](../../docs/examples/example_1_pdcd1_abe8e_converged/report.md)
- **Negative Case Study (Rejected Hazard & Critic Recovery):** [`../../docs/examples/example_2_pdcd1_cbe_rejected_hazard/report.md`](../../docs/examples/example_2_pdcd1_cbe_rejected_hazard/report.md)

## Citation
See [`../../references/citation.bib`](../../references/citation.bib).
