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
3. **Pick Reference Conditions By Evidence Tier, Not By Anatomy** (`conditions.py`):
   - `supported_conditions()` reports, per condition and output type, whether the model has `observed` tracks, an `imputed_supported` imputation, or an unvalidated `imputed_exploratory` one. Only 185 of 2,799 human conditions have **observed** methylation.
   - Choose the hazard reference with the strongest evidence for the modality being trusted, and **stamp the tier into the report** (`Tier2HazardResult.evidence_tiers`, `limiting_evidence_tier`). For CD8 T cells that is `UBERON:0000178_CL:0000625` (blood CD8: observed methylation, splicing, RNA-seq, TF ChIP); the tissue-matched `UBERON:0001155_CL:0000625` (colon CD8) is `imputed_exploratory` for both hazard modalities.
   - Run the tissue-matched condition as a **concordance check** and report disagreement (`tissue_matched_agrees`). Never average a weak-evidence verdict into a strong-evidence one.
4. **Measure Edit Direction, Do Not Infer It From Gradient Sign**:
   - The gradient sign is a first-order heuristic and, audited against ISM on real data, is correct only ~70% of the time. Run an ISM pre-screen (`ism_prescreen_windows`) over compact windows first: one call measures every substitution in the window, so wrong-direction candidates are filtered before they cost verification calls.
5. **Keep The Neighbour Genes The Scorer Already Returned**:
   - The expression scorer returns every gene in the interval. Keeping the non-target rows costs nothing and catches edits that hit a shared regulatory element (`Tier1SelectivityResult.neighbor_gene_lfcs`). At *PDCD1* this caught three candidates dragging *RTP5* along.
6. **Diagnose A Failed Search, Do Not Just Report It** (`benchmarks/hazard_coupling.py`):
   - Test whether efficacy and each hazard are correlated across all scored variants. A hazard that scales with efficacy means the locus is constrained and a better optimiser will not help; one that does not means the search is the limit. Always report the random-edit-only control (rules out selection bias) and the splicing comparison (rules out generic entanglement of model outputs).
7. **Gradients Propose, Scoped ISM Verifies**:
   - Never run brute-force ISM across the full 131 kb context window ($\sim 393\text{k}$ variants). Use `predict_interval_gradients` (`1` call) to rank PAM-valid base-editor windows (`ABE8e` / `CBE`, deaminase window `[4..8]`), then run scoped ISM ($\le 60\text{ bp}$) only on top candidate windows, enforced by `CostBudgetTracker` (`budget.py`) and content-addressed SQLite caching (`cache.py`).
8. **Mandatory AlphaGenome Atlas Deep-Linking**:
   - Every reported variant (`chr:pos:ref>alt`) must include its clickable [AlphaGenome Atlas](https://deepmind.google.com/science/alphagenome/atlas) deep-link (`SingleNucleotideVariant.atlas_url`), formatted as `https://deepmind.google.com/science/alphagenome/atlas?q={chr}:{pos}:{ref}%3E{alt}&m=variant&i={chr}:{pos-190}-{pos+190}`.
9. **Companion Skills Integration (`google-deepmind/science-skills`)**:
   - Chain with `alphagenome-atlas-website-links` for multi-track Ref vs. Alt Sashimi/RNA-seq views (`/atlas/track-predictions`), `jaspar_database` for TF PWM verification, and `ensembl_database` for gene/transcript coordinate lookup.

---

## Standard Workflows & CLI Commands

### 1. Optimize (`optimize`), real data, live API
```bash
uv run crispr-reflexion optimize                 # default scenario crc-pd1-exhaustion, ~1 min
```
Outputs `outputs/report_crc-pd1-exhaustion.json` (provenance, data source, efficacy null, per-iteration Tier 1/Tier 2 results) and a 4-panel figure.
Pre-set rules: direction (silencing must lower the target), efficacy must beat a 40-edit random-edit null, healthy/target ratio <= 0.35, **no neighbour gene above 0.50 of the target's shift**, splice <= 0.15, CpG methylation <= 0.12, zero bystanders. Finalists also get an **off-tissue safety panel** across six non-immune reference tissues. **Do not loosen thresholds to make a candidate pass.**

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

### 4. Is the coupling locus-specific? (`sweep-coupling`)
```bash
uv run crispr-reflexion sweep-coupling           # PDCD1, CTLA4, HAVCR2, LAG3, TIGIT, TOX
```
Runs the identical pipeline (same cells, sketches, thresholds, references) across six checkpoint loci with a 40-edit random null. Last run: CpG rho between +0.44 and +0.65 at every locus (all p < 0.01), positive everywhere. Do not report a locus as methylation-specific, generic or uncoupled: those labels flip with null size and seed.

### 5. Check robustness across seeds, subtypes, lineages & cancers (`sweep-universality`)
```bash
uv run crispr-reflexion sweep-universality
```
Evaluates multi-seed stability (seeds 42/7/99, mean CpG rho = +0.535 +/- 0.049), CRC MMRd vs MMRp subtypes, CD4+ CXCL13+ helper TILs, PDCD1+ gamma-delta TILs, Triple-Negative Breast Cancer CD8+ TILs (Wu et al. 2020), and physical GRCh38 CpG architecture across all 6 loci (`outputs/universality_report.json`). The random-edit-only control is positive but not individually significant (+0.26 to +0.39), and splicing also correlates among random edits, so the result is a moderate association, not a methylation-specific mechanism.

### 6. Run dataset-agnostic on any `.h5ad` (`run-h5ad`)
```bash
uv run crispr-reflexion run-h5ad --h5ad path/or/url.h5ad --gene-id ENSG... --gene-symbol SYMBOL \
  --target-query "..." --healthy-query "..." --cell-type "CD8-positive"
```
Supports `.raw.X`, `.layers["counts"]`, and `.X`, auto-detects donor columns and Ensembl IDs, supports `expr_<GENE>` count predicates inside queries, and automatically resolves the strongest-evidence Tier 2 reference condition via `resolve_reference_pair()`.

### 7. Inspect reference-condition evidence (`conditions`)
```bash
uv run crispr-reflexion conditions "CD8-positive" --modality methylation
uv run crispr-reflexion conditions --refresh      # re-download the shipped table
```

## Reporting rules
- State whether the run converged. Report failures and near-misses, not only passes.
- Always state the data source and backend. `--simulator` and `pd1-exhaustion` (synthetic) are plumbing tests only and must never be presented as results.
- Effect sizes are model predictions with no wet-lab validation; the model compresses true differences (slope ~0.56).
- **Always quote the evidence tier** behind a Tier 2 hazard verdict. An `imputed_exploratory` number is a guess about a guess and must not be reported as equivalent to an `observed` one.
- A correlation between two model predictions is evidence about the model, not a measurement of chromatin. Use a random null of at least 40 edits; a 12-edit null produced a rho of +0.94 that fell to about +0.5. Report the random-edit-only control and the splicing comparison alongside any coupling claim, and do not label loci with a fixed rho cutoff.

## Reference example
[`docs/examples/crc_pd1_real_run.md`](docs/examples/crc_pd1_real_run.md): real colorectal-cancer data, 8 iterations, not converged. The near-misses fail only the CpG-methylation hazard. Efficacy and CpG shift are positively associated (rho = +0.57, p = 2e-5, n = 48), a moderate and not a near-identical relationship.

## Citation
See [`references/citation.bib`](references/citation.bib).
