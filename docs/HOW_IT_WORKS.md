# How it works

The pipeline, and the reasoning behind each design choice. Back to the [README](../README.md).

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
the run correlates efficacy with each hazard across everything it scored and reports the random-edit
control beside it. See [RESULTS.md](RESULTS.md) for what that analysis does and does not support.

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
