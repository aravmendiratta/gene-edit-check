# STATUS (2026-11-15): real-data pipeline built; first live result is "not converged"

- Real dataset wired in: Pelka et al. 2021 colorectal cancer T cells (62 patients) via CELLxGENE. States use the authors' annotations. Synthetic cohort kept for tests/CI only.
- Fabricated "published ground truth" benchmark removed; replaced by `validate-fidelity` (predicted vs observed fold-change, rho +0.65).
- Fixed: Sketch B now held out by patient; efficacy test now direction-aware and null-calibrated (was circular); duplicate edit evaluation; encode_cells wrongly capped at 256; sparse-safe hashing.
- Live run: no candidate passed all checks; two selective ones fail only the methylation hazard. See `docs/examples/crc_pd1_real_run.md`.
- Open: stronger proposer (wrong-direction edits fill iterations 5-8), tissue-matched Tier 2 reference, more seeds/genes, larger null, a second independent dataset.

---

# CRISPR-Reflexion (`crispr-reflexion`) — Revised Engineering & Execution Plan (v2)

> **Research-Use-Only Notice:** `crispr-reflexion` is an in-silico research tool for exploring regulatory sequence grammar and single-cell variant effects using Google DeepMind's AlphaGenome 2 API. Following DeepMind's terms, outputs are candidate computational hypotheses and are **not** validated for clinical, diagnostic, or therapeutic use.

---

## 1. Core Architectural Revisions (Addressing v2 API Physics & Scientific Rigor)

### 1.1 Two-Tier Verifier Architecture (`Selectivity` vs. `Sequence Hazards`)
In `alphagenome.models.v2`, custom single-cell conditioning (`cell_conditioning` via `encode_cells()`, capped at `MAX_CONDITIONING_CELLS = 256`) unlocks **`OutputType.SC_GENE_EXPRESSION_DENOISED`**, whereas specialized variant scorers (`SpliceJunctionScorer()`, `MethylationScorer(width=1201)`, `ContactMapScorer()`) operate under ontology-curated reference conditions (`ModelCondition.from_curies(...)`, e.g., `UBERON:0000178_CL:0001054`).

We encode this distinction directly in the type system (`Tier1SelectivityResult` vs. `Tier2HazardResult`) so no output ever implies single-cell resolution where the underlying scorer uses reference tissue/cell CURIEs:
* **Tier 1 — Single-Cell Selectivity (`cell_conditioning`, `SC_GENE_EXPRESSION_DENOISED`):**
  Evaluates differential target gene expression between the **Target Pathological Cell Sketch** and **Healthy Control Cell Sketches**.
* **Tier 2 — Sequence-Intrinsic Hazards (`ModelCondition.from_curies`, Matched Ontology):**
  Evaluates cryptic splice site creation (`SpliceJunctionScorer()`), local CpG methylation disruption (`MethylationScorer(width=1201)`), and chromatin accessibility shifts on the closest matched `UBERON`/`CL` reference conditions.

### 1.2 "Gradients Propose, Scoped ISM Verifies" (Test-Time Compute Split)
Running `score_ism_variants()` across a full 1-Mb window ($\sim 3 \times 10^6$ variants per condition) is computationally infeasible. Instead, we exploit AlphaGenome 2's native input-gradient endpoint (`predict_interval_gradients` / `predict_interval_gradients_batch` at 100 KB / 1 MB resolution):
1. **Proposal Phase (1 Gradient Call per State):** Call `predict_interval_gradients` centered on the target gene's median TSS for `OutputType.SC_GENE_EXPRESSION_DENOISED` in Target vs. Healthy cell sketches.
2. **Differential Saliency Attribution:** Identify short ($\le 60\text{ bp}$) PAM-accessible windows where the nucleotide gradient magnitude $|\nabla_x \text{Expr}(x \mid C_{\text{target}})|$ is high while $|\nabla_x \text{Expr}(x \mid C_{\text{safe}})|$ is near zero.
3. **Verification Phase (Scoped ISM + Multi-Tier Scorers):** Run `score_ism_variants()` **strictly inside the top candidate 20–40 bp protospacer windows** ($\le 120$ variants total, bounded by `budget.py`), followed by Tier 1 and Tier 2 verification on the surviving base-editor designs.

### 1.3 Batch-Stratified Coreset Selection in AlphaGenome's Own Latent Space
AlphaGenome 2's `encode_cells()` maps raw UMI count matrices into the model's internal cell-state representation, and its docstring explicitly notes that `encode_cells()` **does not perform batch correction**. Running facility location solely in PCA space optimizes a metric the model never sees.
Our revised sketching pipeline (`sc_sketch/coreset.py`):
```text
.h5ad (raw UMI counts, gene symbol → ENSG mapping)
  └── Stratify by donor/batch metadata (prevent single-donor dominance)
      └── Cheap PCA pre-filter to ~2,000–5,000 candidate cells
          └── encode_cells() → AlphaGenome 2 Latent Space Embeddings
              └── Lazy-Greedy Facility Location in AG2 Latent Space
                  └── Emit Disjoint Sketches:
                      ├── Sketch A (256 cells) → Used for Gradient Proposal & Inner Loop
                      └── Sketch B (256 cells, held-out donors/cells) → Used for Final Verification
```
* **Built-In Ablation:** We evaluate sketch stability **with vs. without donor/batch stratification** and measure **Cross-Sketch Generalization ($\text{Sketch } A \to \text{Sketch } B$)** to eliminate Goodhart's Law / self-confirmation bias.

### 1.4 Honest Uncertainty & Ground-Truth Benchmarking
* **Cross-Sketch Stability (`sketch_variance.py`):** Replaces overclaimed "conformal prediction" with honest bootstrap/partition variance across $M$ disjoint 256-cell sketches ($\mu \pm \text{SE}$ and cross-donor transfer correlation).
* **Wet-Lab Ground Truth Anchor (`benchmarks/ground_truth.py`):** Includes a validation benchmark comparing model variant scores against published experimental measurements (**ENCODE CRISPRi-FlowFISH enhancer-gene pairs** / **Gosai et al. *Nature* 2024 MPRA**) reporting Spearman $\rho$ and AUPRC.
* **Strict Provenance & Simulator Labeling (`provenance.py`):** Every CLI run, JSON output, and plot carries an explicit provenance banner (`backend: LIVE_ALPHAGENOME_V2` vs. `backend: SURROGATE_SIMULATOR`, `model_version`, `random_seed`, `api_calls_used`, `cache_hits`). No hardcoded or scripted iteration drama—every trace is generated from a committed random seed.

---

## 2. Revised Repository Layout

```text
crispr-reflexion/
├── pyproject.toml                     # Single primary CLI: `crispr-reflexion`
├── README.md                          # Architecture, two-tier defense, RUO notice, citations
├── PROJECT_PLAN.md                    # This revised plan
├── SKILL.md                           # Agent Skill definition (Claude Code / Cursor / Antigravity)
├── scripts/
│   └── probe_ag2_api.py               # Phase 0 live API probe (tests cell_conditioning + scorers + gradients)
├── src/
│   └── crispr_reflexion/
│       ├── __init__.py
│       ├── cli.py                     # CLI: `probe`, `audit`, `optimize`, `ablate-sketch`, `validate-gt`
│       ├── models.py                  # Pydantic v2 types (explicit Tier1Selectivity vs Tier2Hazard)
│       ├── budget.py                  # Hard caps on API calls, ISM window width (<=60bp), and variants/run
│       ├── cache.py                   # Content-addressed SQLite cache: (sha256(seq), cond_hash, outputs, model_ver)
│       ├── provenance.py              # Tracks seeds, model_version, call counts, LIVE vs SIMULATOR badge
│       ├── ag_backend/                # Named `ag_backend` so it never shadows `import alphagenome`
│       │   ├── __init__.py
│       │   ├── protocol.py            # `VerifierBackend` Protocol (encode_cells, gradients, ism, score_variant)
│       │   ├── live_client.py         # Adapter for `alphagenome.models.v2.dna_client` & `variant_scorers`
│       │   └── simulator.py           # Clearly badged surrogate backend for zero-key offline testing/CI
│       ├── sc_sketch/
│       │   ├── __init__.py
│       │   ├── coreset.py             # Donor-stratified AG2-latent Facility Location + Sketch A / Sketch B split
│       │   └── sketch_variance.py     # Cross-sketch variance & stratified vs. unstratified ablation
│       ├── crispr/
│       │   ├── __init__.py
│       │   ├── constraints.py         # SpCas9/Cas9-NG PAMs, ABE8e/CBE [4..8] window & bystander checks
│       │   └── motifs.py              # JASPAR PWM scanner for mechanistic attribution
│       ├── verifier/
│       │   ├── __init__.py
│       │   └── two_tier.py            # Tier 1 (scRNA selectivity on Sketch A/B) + Tier 2 (CURIE splice/CpG hazards)
│       ├── agent/
│       │   ├── __init__.py
│       │   └── reflexion_loop.py      # Deterministic Gradient-Propose → Scoped-ISM → Reflexion engine (LLM opt-in)
│       ├── benchmarks/
│       │   └── ground_truth.py        # Spearman correlation against measured CRISPRi / MPRA ground truth
│       ├── ui/
│       │   ├── terminal_hud.py        # Rich summary tables & provenance banner
│       │   └── plots.py               # Multi-panel figure with explicit Tier 1 / Tier 2 & backend stamps
│       └── mcp_server/
│           └── server.py              # FastMCP server wrapping the deterministic pipeline
└── tests/
    ├── test_determinism_and_cache.py  # Verifies content-addressed caching & seed reproducibility
    ├── test_budget_guardrails.py      # Verifies BudgetExceededError on oversized ISM intervals (>60bp)
    ├── test_coreset_stratified.py     # Verifies Sketch A/B disjointness & AG2 latent facility location
    ├── test_crispr_constraints.py     # Verifies PAM, ABE8e/CBE windows, and bystander flagging
    └── test_two_tier_loop.py          # Verifies Gradient-Propose -> Scoped-ISM -> Held-out Sketch B check
```

---

## 3. Reordered De-Risked Execution Sequence

### Phase 0: Live API Capability & Quota Probe (`scripts/probe_ag2_api.py`)
- [x] Write `scripts/probe_ag2_api.py` to test against a live `ALPHAGENOME_API_KEY`:
  1. `encode_cells()` latency and output embedding shape.
  2. Whether `score_variant(..., cell_conditioning=...)` accepts `MethylationScorer(width=1201)` and `SpliceJunctionScorer()`, or strictly requires `ModelCondition.from_curies(...)`.
  3. `predict_interval_gradients(..., output_type=OutputType.SC_GENE_EXPRESSION_DENOISED)` output tensor shape.

### Phase 1: Core Vertical Slice (Models, Budget, Cache, Provenance, `VerifierBackend` Protocol)
- [ ] Implement `models.py` (with explicit `Tier1SelectivityResult` and `Tier2HazardResult`), `budget.py`, `cache.py`, `provenance.py`, and `ag_backend/protocol.py` + `simulator.py` + `live_client.py`.

### Phase 2: Batch-Stratified AG2-Latent Coreset Sketcher (`Sketch A` / `Sketch B` Split)
- [ ] Implement `sc_sketch/coreset.py` and `sc_sketch/sketch_variance.py` (supporting both donor-stratified facility location in `encode_cells()` latent space and the stratified vs. unstratified ablation).

### Phase 3: Two-Tier Verifier & Deterministic "Gradients Propose, ISM Verifies" Loop
- [ ] Implement `crispr/constraints.py`, `verifier/two_tier.py`, and `agent/reflexion_loop.py` (deterministic rule engine by default; evaluates winning candidates on held-out `Sketch B`).

### Phase 4: Ground-Truth Validation Benchmark, CLI, Plots, MCP & Tests
- [ ] Implement `benchmarks/ground_truth.py`, `cli.py`, `ui/terminal_hud.py`, `ui/plots.py`, `SKILL.md`, and run `pytest`.
