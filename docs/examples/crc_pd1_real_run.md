# Worked example: PD-1 silencing in exhausted colorectal-cancer TILs

Live AlphaGenome 2 (`alphagenome==2.0.0`), real data (Pelka et al. 2021), seed 42, sketch 128,
pool 600, 8 iterations, 12-edit efficacy null, 10 ISM pre-screen windows.

```bash
crispr-reflexion optimize --max-iterations 8 --null-n 12
```

**Outcome: not converged.** No candidate passed every check. The run's own coupling analysis says
why, and that is the result worth reporting.

---

## Setup the run chose

| | |
|---|---|
| Locus | `chr2:241,793,329-241,924,401` (131 kb around the *PDCD1* TSS) |
| Target state | `Exhausted_CD8_TIL` — 12,738 cells, 62 patients |
| Healthy state | `Normal_Colon_CD8` — 2,626 cells, 36 patients |
| Sketch split | by patient; Sketch A and B share **0** patients |
| Hazard reference | `UBERON:0000178_CL:0000625` blood CD8 — **observed** methylation *and* splicing |
| Tissue-matched check | `UBERON:0001155_CL:0000625` colon CD8 — `imputed_exploratory` for both |
| Efficacy cutoff | 0.0604 = max \|LFC\| over 12 random PAM-valid edits (empirical p ≤ 1/13) |

## ISM pre-screen

240 substitutions measured over 10 × 8 bp windows, covering 49 candidates. **19 were discarded for
moving PD-1 the wrong way** before costing any verification call.

The audit that justifies the pre-screen: the gradient sign agreed with the measured ISM sign only
**70.4%** of the time. Relying on it, as the earlier version did, spent half the iteration budget on
edits that *raised* PD-1.

## Iteration trace

| # | Editor | Edit | PD-1 (held-out B) | Healthy | Worst neighbour | Splice | CpG | Failed |
|---|---|---|---|---|---|---|---|---|
| 1 | CBE | chr2:241858884 C>T | −0.3543 | −0.3637 | — | 0.656 | 0.861 | 3 |
| 2 | CBE | chr2:241859997 C>T | −0.1322 | −0.0310 | — | 0.000 | 0.692 | **1 (CpG)** |
| 3 | CBE | chr2:241859900 C>T | −0.1144 | **+0.0005** | FAM240C +0.020 | 0.000 | 0.415 | 2 |
| 4 | ABE8e | chr2:241859845 A>G | −0.0771 | −0.0119 | — | 0.000 | **0.146** | **1 (CpG)** |
| 5 | ABE8e | chr2:241858956 T>C | −0.1378 | −0.2665 | RTP5 +0.107 | 0.000 | 0.710 | 3 |
| 6 | ABE8e | chr2:241858962 T>C | −0.1153 | −0.1536 | RTP5 +0.063 | 0.000 | 0.415 | 3 |
| 7 | CBE | chr2:241858777 G>A | **+0.0405** | +0.0484 | — | 0.161 | 0.059 | 3 |
| 8 | CBE | chr2:241858973 G>A | −0.2900 | −0.5890 | RTP5 +0.193 | 0.000 | 1.338 | 3 |

Seven of eight moved PD-1 in the intended direction (the earlier gradient-only search managed four).

**Iteration 3** is the most selective edit found: PD-1 down 11.4% in exhausted TILs and **+0.05%** in
normal colon T cells — a selectivity ratio of 0.004. It fails on methylation (0.415) and on three
bystander bases in the editing window.

**Iterations 5, 6 and 8** cluster around chr2:241,858,9xx and all drag ***RTP5*** with them, up to
66% of the on-target effect. That is a shared regulatory element, and the critic learned to avoid it
(`AVOID_SHARED_REGULATORY_ELEMENTS`). The previous version of the pipeline discarded those rows and
would have reported these as clean.

## The near miss, in full

`chr2:241,859,845 A>G` · ABE8e · + strand · protospacer `TATTAGTGATGAGATGCGTC` · PAM `AGG` ·
position 5 · 0 bystanders

| Check | Value | Verdict |
|---|---|---|
| Direction | lowers PD-1 | pass |
| Efficacy | −0.0771 vs cutoff 0.0604 | pass |
| Healthy sparing | −0.0119, ratio 0.15 vs limit 0.35 | pass |
| Cross-sketch transfer | A −0.0802 → B −0.0771 | pass |
| Neighbour genes | none above threshold | pass |
| Splice | 0.000 vs limit 0.15 | pass |
| Bystanders | 0 | pass |
| **CpG methylation** | **0.146 vs limit 0.12** | **fail** |

Two things complicate it further:

- **The references disagree.** Blood CD8 (observed evidence) gives 0.146 → fail. Colon CD8
  (tissue-matched, exploratory evidence) gives 0.087 → pass. Reported side by side.
- **It is not tissue-silent.** The off-tissue panel flags hepatocyte (0.049) and brain (0.040)
  against a bar of 0.039 (50% of the on-target effect).

The threshold was not moved. Whether 0.146 is tolerable needs wet-lab evidence, not a judgement from
whoever set the limit.

## Why the search stopped where it did

Across all 20 scored variants (12 random + 8 candidates):

| | Spearman ρ | p | random-edits-only ρ | p |
|---|---|---|---|---|
| **CpG methylation vs efficacy** | **+0.94** | 7.9×10⁻¹⁰ | **+0.73** | 0.0065 |
| Splice vs efficacy | −0.05 | 0.84 | −0.28 | 0.37 |

Verdict: **target constrained**, classification **methylation-specific**.

At *PDCD1*, the sequence changes that silence the gene are the ones that disturb CpG methylation.
More search iterations will not fix that. The random-edits-only column rules out selection bias (the
null was scored before any candidate was chosen), and the splice row rules out the model simply
moving every output at once.

Across six checkpoint loci this pattern is specific rather than universal — methylation-specific at
*PDCD1* and *LAG3*, generic at *CTLA4*/*HAVCR2*/*TIGIT*, absent at *TOX*:

```bash
crispr-reflexion sweep-coupling
```

Across three random seeds (`42, 7, 99`, mean CpG ρ = `+0.906 ± 0.037`), CRC MMRd/MMRp subtypes (`+0.729`, `+0.886`), CD4+ CXCL13+ (`+0.689`) and γδ (`+0.718`) T-cell lineages, and an independent **Triple-Negative Breast Cancer** cohort (Wu et al. 2020, `+0.804`, splice `−0.455`), *PDCD1* classifies identically as `methylation_specific` / `target_constrained`:

```bash
crispr-reflexion sweep-universality
```

## Caveat

Every number above is a prediction from one model, with no wet-lab validation. The coupling is a
correlation between two AlphaGenome 2 outputs; the two controls make it harder to explain away, but
they do not establish that the mechanism is real chromatin biology.
