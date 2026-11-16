# Related work, and what is and is not new here

Searched in November 2026. This is a scoping note, not a systematic review; check the primary sources before
citing a claim from here.

## What already exists

| Work | What it does | How gene-edit-check differs |
|---|---|---|
| **Ledidi** (Schreiber lab; [bioRxiv 2025](https://www.biorxiv.org/content/10.1101/2025.04.22.650035v2.full)) | Continuous optimisation of a minimal set of edits to a sequence, steered by one or more pretrained models, validated by STARR-seq in Drosophila cells. | It optimises for a model score and does not report a held-out check against that model's own overfitting. It is not built around single-cell cell states or base-editor constraints. |
| **Synthetic enhancer design** (Taskiran et al., *Nature* 2024; [Yin et al., bioRxiv](https://twistbioscience.com/resources/publication/iterative-deep-learning-design-human-enhancers-exploits-condensed-sequence)) | Designs whole new enhancers specific to a cell type, validated experimentally. | Designs new sequences. This project edits an existing endogenous locus with one base edit. |
| **AlphaGenome variant analysis with base editing** ([RHD study, bioRxiv 2026](https://www.biorxiv.org/content/10.64898/2026.01.21.700828v2)) | Uses AlphaGenome to rank non-coding variants, then base-edits the top ones in K562 cells; high-scoring edits changed expression as predicted. | One cell line, no cell-state selectivity, no hazard screen. Its result is the main reason to think AlphaGenome scores can transfer to a bench. |
| **PD-1 enhancer biology** (Weiss et al., [*Nat Immunol* 2024](https://pubmed.ncbi.nlm.nih.gov/39289557/); Wang et al., [bioRxiv 2026](https://www.biorxiv.org/content/10.64898/2026.01.27.702006.full.pdf)) | Deleting an exhaustion-specific *Pdcd1* enhancer in mouse CD8 T cells gave intermediate PD-1 and better function than a knockout. A CRISPRi screen mapped human *PDCD1* enhancers more than 12 kb from the TSS. | These are the experimental context. They motivate the optional dose ceiling and are the natural place to test the edits proposed here. |
| **Oracle-overfitting in design** ([RL regulatory design](https://arxiv.org/pdf/2503.07981); [diffusion reward fine-tuning](https://arxiv.org/pdf/2410.13643); [Overconfident Oracles, ICML 2024](https://icml.cc/virtual/2024/36750)) | Show that optimising against a model drifts from the truth, and recommend a separate evaluation oracle. | Applies the same principle to AlphaGenome 2 by splitting the *cells*: proposals are made on one set of donors and judged on a different set. |

## What this project adds

Stated narrowly, because no single combination below is a claim that nobody has touched any part of it:

1. **Single-cell-conditioned selectivity for an endogenous base edit.** Edits are scored in a defined tumour
   cell state and in a normal-tissue state using AlphaGenome 2's cell conditioning, and an edit must act in the
   first and not the second.
2. **Patient-disjoint proposal and verification.** The cells used to propose an edit and the cells used to
   judge it never share a patient. This is the held-out-oracle idea applied to cells rather than to a second
   model.
3. **An empirical random-edit null as the efficacy bar,** and bootstrap intervals rather than fixed cutoffs for
   claims about how efficacy relates to hazards. This project retracted its own first headline number after
   the larger null, which is documented in `docs/RESULTS.md`.
4. **One verdict that combines selectivity, splicing, CpG methylation, neighbouring genes, bystander bases and
   off-tissue effects,** using fixed thresholds set before any candidate is scored.
5. **An optional dose ceiling** (`--max-effect`) for when partial reduction is the goal, motivated by the mouse
   result above. It is a design option, not a validated claim.

## What it does not claim

- No wet-lab validation. The RHD study is evidence that AlphaGenome scores can predict base-edit outcomes in
  one setting; it is not evidence about *PDCD1* in T cells.
- The core idea of optimising edits with a sequence model is not new (Ledidi and others).
- The best next step is experimental: test the top *PDCD1* candidates against the enhancer positions in the
  papers above.
