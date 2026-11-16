# CRISPR-Reflexion

> A verifier-guided search for **cell-state-selective CRISPR base edits**, using Google DeepMind's
> **AlphaGenome 2** as the judge. Test case: turn **PD-1 (*PDCD1*)** down in the exhausted T cells
> inside real colorectal tumours while leaving normal-colon T cells alone.
>
> *Research use only. Outputs are computational hypotheses from the AlphaGenome 2 API, not validated
> for clinical, diagnostic or therapeutic use.*

---

## The headline result

The search **did not find a usable edit**. Of 8 verified candidates, 7 were blocked by the CpG-methylation
check, and the closest near miss fails it at 0.146 against a limit of 0.12.

The finding underneath that is an **association between how much an edit silences PD-1 and how much it
disturbs CpG methylation**. It is real in this data, but smaller and less specific than a first pass
suggested, so this section separates what holds from what does not.

**What holds**

- Across all 48 variants scored in the main run, efficacy and CpG shift are positively correlated:
  **Spearman rho = +0.57 (p = 2.2e-05)**.
- It reproduces across three independent seeds (patient splits and random-edit draws): rho = +0.48, +0.56, +0.57, all p < 1e-3.
- The sign is positive in **all 6 cohorts** (rho from +0.25 to +0.73) and at **all 6 checkpoint loci** (rho from +0.44 to +0.65).

**What does not hold**

- **It is not near-identity.** An earlier version of this README reported rho = +0.94 at *PDCD1*. That came
  from 20 variants and a 12-edit random null. With a 40-edit null it is about +0.5, and I have retracted
  the stronger claim.
- **The selection-bias control is inconclusive.** Among the 40 random edits alone (chosen before any
  candidate was verified), rho = +0.26 (p = 0.106) in the main run, and +0.26, +0.37, +0.39 across the three seeds. Consistently positive, but
  not individually significant, so I cannot rule out that the pooled correlation is partly produced by
  how candidates were selected.
- **It is not specific to methylation.** Among random edits, splicing also correlates with efficacy
  (rho = +0.33, p = 0.036). The earlier "splicing shows nothing" control came
  from a smaller null and does not survive.
- **The labels flip.** The classifier that called loci "methylation-specific", "generic" or "uncoupled"
  uses a fixed rho cutoff, and the same data crosses it depending on null size and seed. I no longer
  treat those labels as findings.

**Reading it.** The simplest explanation is generic rather than *PDCD1*-specific: edits that move predicted
expression a lot also tend to shift predicted methylation across a 1.2 kb window. This design cannot
tell that apart from a CpG-borne regulatory mechanism, and both are statements about AlphaGenome 2's
predictions, not measurements of chromatin. No wet-lab validation was performed.

![Summary figure](outputs/figure_crc-pd1-exhaustion.png)

---

## What the tool does, in plain words

Tumours wear T cells out ("exhaustion"). PD-1 is the brake on them, and the PD-1 drugs
(Keytruda, Opdivo) release it. The question here is: *could a single DNA letter be edited so PD-1
turns down only in the exhausted tumour T cells, and not in healthy ones?*

1. **Real cells in.** 77k T/NK cells from 62 colorectal-cancer patients, tumour and normal colon.
2. **Shrink to a sketch.** AlphaGenome 2 accepts at most 256 cells to describe a cell state, so
   representative cells are chosen by facility-location coreset in the model's own latent space —
   with patients kept disjoint between the "proposal" sketch and the held-out "verification" sketch.
3. **Gradients propose.** One gradient call shows which DNA positions matter for PD-1 in those
   cells. Only edits a base editor can actually make survive (PAM site, editing window, bystanders).
4. **ISM measures direction.** In-silico mutagenesis over compact windows measures what each
   substitution really does, *before* any verification budget is spent on it.
5. **Two tiers verify.** Right direction, beats random edits, spares healthy T cells, doesn't drag
   neighbouring genes, no new splice sites, no large methylation shift, nothing in off-target tissues.
6. **The critic learns from each failure** and re-ranks what to try next.

---

## Why the architecture looks like this

Each of these was a bug or a blind spot found by running the thing, not a design flourish.

### Gradients say *where*, not *which way* — so measure it
The gradient sign is a first-order heuristic: positive means "more of this base → more expression",
so replacing it should lower expression. Audited against ISM ground truth, **it is right only 70% of
the time.** Trusting it wasted whole verification rounds on edits that *raised* PD-1. The pipeline
now spends 10 ISM calls up front, each measuring every substitution in an 8 bp window, and the
verifier only ever sees candidates whose direction is already known.

*Effect: wrong-direction candidates reaching verification went from 4/8 to 1/8.*

### A hazard verdict is only as good as its evidence
`supported_conditions()` returns 2,799 human reference conditions and, per modality, whether the
model has **observed** tracks, a **validated** imputation, or an **unvalidated** one. Only 185 of
2,799 have observed methylation. So the reference condition isn't hard-coded — it is chosen by
evidence, and every hazard call is stamped with the tier it came from.

This resolved a real question. The anatomically correct condition is colon CD8 T cell
(`UBERON:0001155_CL:0000625`), but it is `imputed_exploratory` for both methylation and splicing.
Blood CD8 (`UBERON:0000178_CL:0000625`) is the only CD8 condition that is **observed** for
methylation, splicing, RNA-seq *and* TF ChIP. The pipeline uses blood CD8 as primary and colon CD8
as a **concordance check** — and reports when they disagree rather than averaging them. On the
finalist edit they *do* disagree (blood 0.146 → fail, colon 0.087 → pass), which is exactly the
kind of thing that should be shown, not smoothed over.

```bash
crispr-reflexion conditions "CD8-positive" --modality methylation
```

### The expression scorer already returns the neighbours — don't throw them away
Scoring a variant returns every gene in the interval, and the old code kept one row. Keeping the
rest costs nothing and catches edits that hit a shared regulatory element: three candidates were
caught dragging ***RTP5*** along with PD-1 (up to 66% of the on-target effect).

### An edit is delivered to cells, not to a cell state
Tier 1 asks whether an edit is selective *between cell states*. It cannot ask what happens in the
liver. The finalist is scored against six non-immune reference tissues; it **flags hepatocyte
(0.049) and brain (0.040)** against a bar of 0.039.

### A failed search has two very different explanations
Either the optimiser was bad, or the target is constrained. Those call for opposite responses, so
the run measures which one it is and says so — that analysis is the headline result above.

---

## Data

**Pelka et al., *Cell* 2021**, via CZ CELLxGENE Discover ("T, NK, ILC cells", collection
`4d82bd8e-…`, dataset `7674f8f4-…`). Raw UMI counts, Ensembl IDs, 10x 3′ v2/v3. States are the
**authors' own annotations**:

| State | Definition | Cells | Patients |
|---|---|---|---|
| Target `Exhausted_CD8_TIL` | tumour CD8 T cells, clusters cTNI14/15/16 (`CD8+ CXCL13+`) | 12,738 | 62 |
| Healthy `Normal_Colon_CD8` | CD8 T cells from adjacent normal colon | 2,626 | 36 |

Sanity check on raw counts (per 10k UMI): *PDCD1* 1.80 vs 0.26, *CXCL13* 25.97 vs 0.03,
*IL7R* 0.92 vs 4.89. The 396 MB file downloads automatically on first run into `data/` (git-ignored).

---

## Does the model actually see the cell state?

Before trusting any edit prediction, the conditioning itself has to be checked. `validate-fidelity`
compares AlphaGenome 2's predicted tumour-vs-normal fold-change against what the real cells show,
across a fixed 24-gene immune panel:

- **Spearman ρ = +0.65 (p = 0.0005)**, **94% sign agreement** on the 17 genes that genuinely move.
- The model **compresses** real differences to ~56% of their true size (*PDCD1*: observed ln FC
  +1.65, predicted +0.43).

This validates the conditioning. It says nothing about edit effect sizes, which have no ground truth here.

## Coresets

Facility-location sketches in AlphaGenome's latent space cut MMD² by ~50% versus random
subsampling (0.00135 vs 0.00280 target; 0.00140 vs 0.00277 healthy). Sketch A and Sketch B share
**zero patients**.

---

## Quickstart

```bash
uv venv --python 3.12 .venv
uv pip install --python .venv/Scripts/python.exe -e ".[scverse,dev]"
cp .env.example .env            # add ALPHAGENOME_API_KEY (free, non-commercial)

crispr-reflexion optimize                  # main run: Pelka CRC real data, live API
crispr-reflexion sweep-coupling            # 6 checkpoint loci comparison
crispr-reflexion sweep-universality        # multi-seed, subtypes, lineages, TNBC, CpG architecture
crispr-reflexion validate-fidelity         # predicted vs observed fold-changes
crispr-reflexion ablate-sketch             # coreset vs random subsampling
crispr-reflexion conditions "CD8-positive" # evidence tiers per reference condition
pytest                                     # 53 offline tests (no API, no download)
```

### Bring your own single-cell atlas (`run-h5ad`)

The verifier is dataset-, cancer-, and lineage-agnostic. Point it at any local `.h5ad` file or CZ CELLxGENE Discover URL containing raw integer UMI counts (`.raw.X`, `.layers["counts"]`, or `.X`), select your target/control populations via pandas queries (including optional `expr_<GENE>` count predicates), and it automatically picks the highest-evidence Tier 2 reference condition for your `--cell-type`:

```bash
crispr-reflexion run-h5ad \
  --h5ad https://datasets.cellxgene.cziscience.com/4d036b2c-6135-4515-9c44-a5aa7a2e32e5.h5ad \
  --gene-id ENSG00000188389 --gene-symbol PDCD1 \
  --target-query "celltype_final == 'CD8+ T-cells' and (expr_LAG3 > 0 or expr_HAVCR2 > 0 or expr_CXCL13 > 0)" \
  --healthy-query "celltype_final == 'CD8+ T-cells' and expr_LAG3 == 0 and expr_HAVCR2 == 0 and expr_CXCL13 == 0 and (expr_IL7R > 0 or expr_CCR7 > 0)" \
  --cell-type "CD8-positive" --tissue "breast"
```

---

## Robustness checks (`sweep-universality`, `sweep-coupling`)

All runs use a 40-edit random null, 5 iterations, 8 ISM pre-screen windows, live AlphaGenome 2.
Every cohort is split into proposal and verification sketches **by patient** with zero overlap.

### Seeds (CRC CD8+ TILs, independent patient splits and null draws)

| Seed | CpG rho | CpG p | Random-edit-only rho | Splice rho |
|---|---|---|---|---|
| 42 | +0.48 | 8.5e-04 | +0.26 | +0.20 |
| 7 | +0.56 | 7.5e-05 | +0.37 | -0.04 |
| 99 | +0.57 | 4.3e-05 | +0.39 | -0.01 |

Mean CpG rho = +0.535 +/- 0.049; mean random-edit-only rho = +0.338.

### Cohorts, subtypes and lineages (*PDCD1*, seed 42)

| Context | Target cells / patients | Auto-resolved Tier 2 reference [evidence] | CpG rho | p | Random-edit-only rho | Splice rho |
|---|---|---|---|---|---|---|
| CRC All CD8+ TILs (Pelka 2021) | 12,738 / 62 | blood CD8-positive, alpha-beta T cell [observed] | +0.48 | 8.5e-04 | +0.26 | +0.20 |
| CRC MMRd / MSI-H CD8+ TILs | 10,582 / 34 | blood CD8-positive, alpha-beta T cell [observed] | +0.30 | 0.047 | +0.00 | +0.22 |
| CRC MMRp / MSS CD8+ TILs | 2,156 / 28 | blood CD8-positive, alpha-beta T cell [observed] | +0.73 | 1.6e-08 | +0.61 | -0.28 |
| CRC CD4+ CXCL13+ Helper TILs | 2,135 / 60 | blood CD4-positive, alpha-beta T cell [observed] | +0.25 | 0.098 | -0.06 | -0.30 |
| CRC PDCD1+ gamma-delta TILs | 3,204 / 58 | blood gamma-delta T cell [imputed_exploratory] | +0.42 | 0.004 | +0.18 | -0.00 |
| Breast Cancer (TNBC) CD8+ TILs (Wu 2020) | 2,162 / 5 | blood CD8-positive, alpha-beta T cell [observed] | +0.38 | 0.010 | +0.12 | -0.19 |

Sources: CRC from Pelka et al. 2021; breast from Wu et al., *EMBO J* 2020. The TNBC row has 5 patients, so its
sketches are split 2/3 (target) and 1/4 (healthy): patient-disjoint, but the healthy proposal sketch is a single
patient. Its groups are defined by marker expression and the healthy group is non-exhausted tumour T cells,
not normal tissue. The gamma-delta row uses an `imputed_exploratory` reference, the weakest evidence tier.

### Checkpoint loci (same cells, same thresholds)

| Gene | GC | CpG / 100 bp | CpG rho (p) | Random-edit-only rho | Splice rho |
|---|---|---|---|---|---|
| PDCD1 | 61% | 2.86 | +0.48 (8.5e-04) | +0.26 | +0.20 |
| CTLA4 | 41% | 0.64 | +0.65 (1.2e-06) | +0.51 | +0.71 |
| HAVCR2 | 43% | 1.20 | +0.54 (1.4e-04) | +0.34 | +0.58 |
| LAG3 | 61% | 3.02 | +0.57 (4.6e-05) | +0.39 | +0.29 |
| TIGIT | 47% | 0.70 | +0.50 (5.4e-04) | +0.31 | +0.59 |
| TOX | 52% | 4.80 | +0.44 (0.003) | +0.36 | +0.00 |

CpG density of the surrounding sequence (computed directly from GRCh38, independent of the model) does not
clearly predict coupling strength: Spearman rho between CpG density and coupling across the six loci is
-0.60 (n = 6, so there is almost no power to say anything either way).

---

## How it is kept honest

- **Held out by patient.** Sketch B never contains a patient from Sketch A.
- **Null-calibrated efficacy.** A candidate must beat 40 random PAM-valid edits at the same locus
  (empirical p ≤ 1/41). Thresholds are fixed before any candidate is verified.
- **Direction matters.** Silencing means the edit must *lower* PD-1; a large increase fails.
- **Tiers never mix.** Tier 1 is single-cell-conditioned expression; Tier 2 is reference-condition
  hazards. The API forbids both conditions in one call, and the report never claims single-cell
  resolution for a Tier 2 number.
- **Evidence tiers are reported**, so an `imputed_exploratory` verdict can be discounted.
- **The simulator is not a result.** `--simulator` is a scripted offline stand-in that hard-codes
  where the "selective" site is. It prints a red warning and exists for plumbing tests and CI only.
- **Thresholds were not tuned after seeing results.** The closest miss fails methylation at 0.146
  against a limit of 0.12. The limit was not moved.

## Honest limitations

- **No wet-lab validation.** Every number is a prediction from one model. "Selective" means
  AlphaGenome 2 says so.
- The coupling result is a correlation **between two model outputs**. The random-edit control rules
  out selection bias and the splice control argues against generic entanglement, but neither proves
  the mechanism is real chromatin biology.
- Predicted effect sizes are compressed (slope ~0.56), and the efficacy null has only 12 edits.
- Tier 2 hazards use a **blood** CD8 reference, because the colon one is exploratory-imputed. The
  tissue-matched cross-check disagrees on the finalist.
- The healthy comparator is normal-colon CD8 T cells — mostly IL7R⁺/GZMK⁺ memory-like, not naive.
- Tested across two independent solid-tumour cohorts (Pelka et al. 2021 CRC, 62 patients; Wu et al. 2020 TNBC, 5 patients), three random seeds, and three T-cell lineages, but still does not include blood malignancies or healthy non-epithelial organs beyond reference conditions; the six-locus sweep remains CRC-anchored.

## Citation

Pelka et al. 2021 (data); Avsec et al. / Google DeepMind (AlphaGenome). See
[references/citation.bib](references/citation.bib).
