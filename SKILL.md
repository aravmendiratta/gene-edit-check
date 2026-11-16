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

### 1. Optimize (`optimize`), real data, live API
```bash
uv run crispr-reflexion optimize                 # default scenario crc-pd1-exhaustion, ~1 min
```
Outputs `outputs/report_crc-pd1-exhaustion.json` (provenance, data source, efficacy null, per-iteration Tier 1/Tier 2 results) and a 4-panel figure.
Pre-set rules: direction (silencing must lower the target), efficacy must beat a 12-edit random-edit null, healthy/target ratio <= 0.35, splice <= 0.15, CpG methylation <= 0.12, zero bystanders. **Do not loosen thresholds to make a candidate pass.**

### 2. Validate that conditioning reproduces real biology (`validate-fidelity`)
```bash
uv run crispr-reflexion validate-fidelity
```
Predicted vs observed tumor-vs-normal ln fold-change over a fixed 24-gene panel (live backend only). Last run: Spearman rho = +0.65 (p = 0.0005), 94% sign agreement on 17 changing genes, slope 0.56 (model compresses true differences).

### 3. Coreset ablation (`ablate-sketch`)
```bash
uv run crispr-reflexion ablate-sketch
```
Facility-location sketches in `encode_cells()` space vs random subsampling (about -50% MMD^2 on real data). Sketch A and B are split by patient.

## Reporting rules
- State whether the run converged. Report failures and near-misses, not only passes.
- Always state the data source and backend. `--simulator` and `pd1-exhaustion` (synthetic) are plumbing tests only and must never be presented as results.
- Effect sizes are model predictions with no wet-lab validation; the model compresses true differences (slope ~0.56).

## Reference example
[`docs/examples/crc_pd1_real_run.md`](docs/examples/crc_pd1_real_run.md): real colorectal-cancer data, 8 iterations, not converged; two selective candidates fail only the CpG-methylation hazard.

## Citation
See [`references/citation.bib`](references/citation.bib).
