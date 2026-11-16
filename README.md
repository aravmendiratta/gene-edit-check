# CRISPR-Reflexion

> A verifier-guided search for **cell-state-selective CRISPR base edits**, using Google DeepMind's
> **AlphaGenome 2** as the judge. Test case: turn **PD-1 (*PDCD1*)** down in the exhausted T cells
> inside real colorectal tumours while leaving normal-colon T cells alone.
>
> *Research use only. Outputs are computational hypotheses from the AlphaGenome 2 API, not validated
> for clinical, diagnostic or therapeutic use.*

---

## The headline result

The search **did not find a usable edit** — and the interesting part is *why*.

Across 20 scored variants at the *PDCD1* locus, how much an edit silences PD-1 and how much it
disrupts CpG methylation are almost the same quantity:

**Spearman ρ = +0.94 (p = 8 × 10⁻¹⁰, n = 20)**

Three things make that more than a curiosity:

1. **It isn't selection bias.** It holds among the 12 *randomly chosen* edits alone (ρ = +0.73,
   p = 0.0065), which were scored before any candidate was verified.
2. **It isn't "the model's outputs all move together".** Splicing, a different output head, shows
   no coupling at all (ρ = −0.05, p = 0.84).
3. **It isn't universal.** Running the identical pipeline across six checkpoint loci gives a
   structured answer, not a flat one.

| Gene | CpG ρ | random-only ρ | splice ρ | gap | class | blocked by |
|---|---|---|---|---|---|---|
| **PDCD1** | **+0.92** | +0.75 | −0.30 | **+1.22** | **methylation-specific** | methylation 5/5 |
| **LAG3** | **+0.82** | +0.55 | +0.10 | **+0.72** | **methylation-specific** | methylation 5/5 |
| HAVCR2 | +0.90 | +0.68 | +0.88 | +0.02 | generic hazard coupling | methylation 5/5 |
| CTLA4 | +0.76 | +0.21 | +0.81 | −0.05 | generic hazard coupling | healthy selectivity 5/5 |
| TIGIT | +0.62 | +0.18 | +0.65 | −0.03 | generic hazard coupling | healthy selectivity 5/5 |
| TOX | +0.30 (ns) | −0.02 | +0.00 | +0.30 | uncoupled | healthy selectivity 5/5 |

At *PDCD1* and *LAG3*, the sequence features that carry expression appear to be the CpGs
themselves, so a single base edit that silences the gene without moving methylation may not exist
there. That is a mechanistically sensible result — the *PDCD1* regulatory regions are CpG-dense and
their demethylation is part of how T-cell exhaustion is established — and it is a **falsifiable
claim about where to stop looking**, which is more useful than a candidate nobody checked.

**This is a statement about what AlphaGenome 2 has learned, not a measurement of chromatin.** No
wet-lab validation was performed. See [Honest limitations](#honest-limitations).

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

crispr-reflexion optimize                  # the main run: real data, live API
crispr-reflexion sweep-coupling            # the six-locus comparison above
crispr-reflexion validate-fidelity         # predicted vs observed fold-changes
crispr-reflexion ablate-sketch             # coreset vs random subsampling
crispr-reflexion conditions "CD8-positive" # evidence tiers per reference condition
pytest                                     # 45 offline tests (no API, no download)
```

---

## How it is kept honest

- **Held out by patient.** Sketch B never contains a patient from Sketch A.
- **Null-calibrated efficacy.** A candidate must beat 12 random PAM-valid edits at the same locus
  (empirical p ≤ 1/13). Thresholds are fixed before any candidate is verified.
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
- One dataset, one seed, one cancer type. The six-locus sweep shares all of these.

## Citation

Pelka et al. 2021 (data); Avsec et al. / Google DeepMind (AlphaGenome). See
[references/citation.bib](references/citation.bib).
